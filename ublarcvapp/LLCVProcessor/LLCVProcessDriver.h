#ifndef __LLCV_PROCESS_DRIVER_H__
#define __LLCV_PROCESS_DRIVER_H__

/**
 * 
 * \brief class provides means for processing both larcv and larlite files.
 *
 */

#include <string>
#include "larcv/core/Processor/ProcessDriver.h"

// larlite
#include "core/DataFormat/storage_manager.h"

#include "ublarcvapp/LArliteHandler/LArliteManager.h"

namespace ublarcvapp {
namespace llcv {

  class LLCVProcessDriver : public larcv::ProcessDriver {

  public:

    LLCVProcessDriver( std::string name );
    virtual ~LLCVProcessDriver() {};

    /// runs parent configure and also looks for larlite storagemanager and llcv processor modules
    void configure( const std::string config_file ); 
    void configure( const larcv::PSet& pset );
                   
    //void intialize(); ///< call parent initialize, but also setup larlite

    /* void override_larlite_input_files( const std::vector<std::string>& flist ); */
    /* void override_larcv_input_files( const std::vector<std::string>& flist ); */
    
    /* void override_larlite_output_file( const std::string fname ); */
    /* void override_larcv_output_file( const std::string fname ); */

    /* void process_entry( bool autosave_entry ); ///< replaces parent process_entry, includes larlite modules */

    larlite::storage_manager& io_larlite() { return _io_larlite; };

  protected:

    ublarcvapp::LArliteManager _io_larlite;

    void _do_larlite_config( larlite::storage_manager& ioman, larcv::PSet& pset );
    larlite::data::DataType_t _get_enum_fromstring( std::string name );
    
  };
  
}
}

#endif
