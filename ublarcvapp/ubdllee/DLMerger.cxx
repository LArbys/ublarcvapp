#include "DLMerger.h"
#include <fstream>
#include <assert.h>
#include <vector>
#include <string>

#include "nlohmann/json.hpp"

#include "MergeDLInteraction.h"

using json = nlohmann::json;

namespace ublarcvapp {
namespace ubdllee {

  void DLMerger::configure( const std::string process_driver_cfg,
                            const std::string json_filelist,
                            const std::string file_dir ) {

    // config the process driver
    larcv::ProcessDriver::configure(process_driver_cfg);

    // set input file directory
    _file_dir = file_dir;
    
    // read a JSON file
    json j;    
    std::ifstream infile(json_filelist);
    infile >> j;

    for ( auto& jfileset : j["filesets"] ) {
      for ( json::iterator it = jfileset.begin(); it != jfileset.end(); ++it ) {
        if ( it.key()=="runindex" || it.key()=="subrunindex" || it.key()=="sample" ) {
          continue;
        }

        std::string filename = it.key();
        if ( _filelist_map.find(filename)==_filelist_map.end() ) {
          _filelist_map[filename] = std::vector<std::string>();
        }

        _filelist_map[filename].push_back( it.value() );
      }
    }

    // check numbers of files. needs to the same for all.
    int nfiles = -1;
    LARCV_NORMAL() << "=====================================================" << std::endl;
    LARCV_NORMAL() << "File sets loaded (via json): " << std::endl;
    for ( auto it=_filelist_map.begin(); it!=_filelist_map.end(); it++ ) {
      LARCV_NORMAL() << "  " << it->first << ": size=" << it->second.size() << std::endl;
      if ( nfiles<0 )
        nfiles = it->second.size();
      else if ( nfiles!=(int)it->second.size() ) {
        LARCV_CRITICAL() << "does not match number of files." << std::endl;
        assert(false);
      }        
    }
    LARCV_NORMAL() << "=====================================================" << std::endl;    
  }


  void DLMerger::initialize() {
    // these define the reco-portions of an interaction
    LARCV_NORMAL() << "load files into IOManager" << std::endl;
    
    std::vector< std::string > larcv_required_files = { "supera",
                                                        "tagger-larcv",
                                                        "ssnet",
                                                        "vertex-larcv" };
                                                         

    std::vector< std::string > input_files;
    for ( auto const& fname : larcv_required_files ) {
      if ( _filelist_map.find(fname)==_filelist_map.end() ) {
        LARCV_CRITICAL() << " missing filetype=" << fname << " that is required." << std::endl;
        throw std::runtime_error("missing required LArCV filetype");
      }
      for ( auto const& f : _filelist_map[fname] ) {
        if ( _file_dir!="" )
          input_files.push_back(_file_dir + "/" + f);
        else
          input_files.push_back(f);          
      }
    }
    
    // these define analysis quantities for the interactions
    std::vector< std::string > larcv_optional_files = { "michelid-larcv",
                                                        "nueid-larcv",
                                                        "larcvtruth" };
    for ( auto const& fname : larcv_optional_files ) {
      if ( _filelist_map.find(fname)!=_filelist_map.end() ) {
        for ( auto const& f : _filelist_map[fname] ) {
          if ( _file_dir!="" )
            input_files.push_back(_file_dir + "/" + f);
          else
            input_files.push_back(f);          
        }
      }
    }


    std::vector< std::string > input_larlite_files;
    std::vector< std::string > larlite_required_files = { "tracker-larlite",
                                                          "shower-reco" };
    for ( auto const& fname : larlite_required_files ) {
      if ( _filelist_map.find(fname)==_filelist_map.end() ) {
        LARCV_CRITICAL() << " missing filetype=" << fname << " that is required." << std::endl;
        throw std::runtime_error("missing required larlite filetype");
      }
      for ( auto const& f : _filelist_map[fname] ) {
        if ( _file_dir!="" )
          input_larlite_files.push_back(_file_dir + "/" + f);
        else
          input_larlite_files.push_back(f);
      }
    }
    
    LARCV_NORMAL() << "Number of LArCV input files: " << input_files.size() << std::endl;
    LARCV_NORMAL() << "Number of larlite input files: " << input_larlite_files.size() << std::endl;    
    larcv::ProcessDriver::override_input_file( input_files );
    
    larcv::ProcessDriver::initialize();

    // load up the larlite manager in the MergerDLInteraction Process

    MergeDLInteraction* merger = (MergeDLInteraction*)larcv::ProcessDriver::process_ptr( larcv::ProcessDriver::process_id( "MergeDLInteraction" ) );
    merger->setup_larlite_io( input_larlite_files );
    
  }
  
}
}
