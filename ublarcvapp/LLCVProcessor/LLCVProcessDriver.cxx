#include "LLCVProcessDriver.h"
#include "larcv/core/Base/LArCVBaseUtilFunc.h"

namespace ublarcvapp {
namespace llcv {

  LLCVProcessDriver::LLCVProcessDriver( std::string name )
    : larcv::ProcessDriver(name)
  {
    
  }

  void LLCVProcessDriver::configure( const std::string config_file ) {

    auto main_cfg = larcv::CreatePSetFromFile(config_file);

    
        LARCV_INFO() << "Configuration driver [" << name() << "] config_file=" << config_file << std::endl;
    if(!main_cfg.contains_pset(name())) {
      LARCV_CRITICAL() << "ProcessDriver configuration (" << name() << ") not found in the config file (dump below)" << std::endl
		       << main_cfg.dump()
		       << std::endl;
      throw larcv::larbys();
    }
    const larcv::PSet& cfg = main_cfg.get<larcv::PSet>(name());
    set_verbosity( (::larcv::msg::Level_t)cfg.get<int>("Verbosity") );
    // we call our configure
    configure(cfg);
    
  }


  void LLCVProcessDriver::configure( const larcv::PSet& pset ) {
    larcv::ProcessDriver::configure( pset );

    LARCV_INFO() << "configure larlite" << std::endl;
    
    // get storage manager for larlite
    larcv::PSet cfg_storage_man = pset.get<larcv::PSet>("storage_manager");

    _do_larlite_config( _io_larlite, cfg_storage_man );
    
  }

  void LLCVProcessDriver::_do_larlite_config( larlite::storage_manager& ioman, larcv::PSet& pset ) {

    LARCV_INFO() << "configure larlite::storage_manager" << std::endl;
    
    // mode
    int iomode = pset.get<int>("IOMode");
    ioman.set_io_mode( (larlite::storage_manager::IOMode_t)iomode );

    std::vector<std::string> input_files = pset.get< std::vector<std::string> >( "InputFiles", std::vector<std::string>() );
    for ( auto const& input_file : input_files )
      ioman.add_in_filename( input_file );
    
    std::string outfilename = pset.get<std::string>( "OutFileName", "" );
    if ( iomode==1 || iomode==2 ) {
      if ( ioman.output_filename().empty()) {
	if ( outfilename.empty()) {
          std::string errmsg = "Larlite file is set to write mode, but does not have an output file name.";
	  LARCV_CRITICAL() << errmsg << std::endl;
	  throw std::runtime_error( errmsg );
	}
	ioman.set_out_filename( outfilename );
	assert(!ioman.output_filename().empty());
      }
    }

    // specified read/write datatypes
    auto readonlyvars  = pset.get<std::vector<std::string> >( "ReadOnlyDataTypes", std::vector<std::string>() );
    auto readonlyname  = pset.get<std::vector<std::string> >( "ReadOnlyProducers", std::vector<std::string>() );
    auto writeonlyvars = pset.get<std::vector<std::string> >( "WriteOnlyDataTypes", std::vector<std::string>() );
    auto writeonlyname = pset.get<std::vector<std::string> >( "WriteOnlyProducers", std::vector<std::string>() );

    if ( readonlyvars.size()!=readonlyname.size() ) {
      std::cout << "ERROR: number of read-only data types and names are not the same." << std::endl;
      assert(false);
    }
    if ( writeonlyvars.size()!=writeonlyname.size() ) {
      std::cout << "ERROR: number of write-only data types and names are not the same." << std::endl;
      assert(false);
    }

    for ( size_t i=0; i<readonlyvars.size(); i++ ) {
      std::string rdvar  = readonlyvars.at(i);
      larlite::data::DataType_t datat = _get_enum_fromstring( rdvar );
      if ( datat!=larlite::data::kUndefined ) ioman.set_data_to_read( datat, readonlyname.at(i) );
    }

    for ( size_t i=0; i<writeonlyvars.size(); i++ ) {
      std::string rdvar  = writeonlyvars.at(i);
      larlite::data::DataType_t datat = _get_enum_fromstring( rdvar );
      if ( datat!=larlite::data::kUndefined ) ioman.set_data_to_write( datat, writeonlyname.at(i) );
    }

    int verbosity = pset.get<int>("Verbosity",2);
    LARCV_INFO() << "set larlite::storage_manager verbosity: " << verbosity << std::endl;
    ioman.set_verbosity( (larlite::msg::Level)verbosity );

  }
  
  larlite::data::DataType_t LLCVProcessDriver::_get_enum_fromstring( std::string name ) {
    for (int i=0; i<larlite::data::kDATA_TYPE_MAX; i++) {
      if ( larlite::data::kDATA_TREE_NAME[i]==name )
        return (larlite::data::DataType_t)i;
    }
    return larlite::data::kUndefined;
  }

  void LLCVProcessDriver::initialize() {
    // initialize parent class
    LARCV_INFO() << "initialize parent class" << std::endl;    
    larcv::ProcessDriver::initialize();

    // initialize io manager
    LARCV_INFO() << "initialize larlite::storage_manager" << std::endl;
    _io_larlite.open();
    
  }
  
}
}
