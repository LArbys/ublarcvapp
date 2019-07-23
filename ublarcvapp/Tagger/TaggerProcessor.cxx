#include "TaggerProcessor.h"

#include "TaggerCROITypes.h"

namespace ublarcvapp {
  namespace tagger {
    
    TaggerProcessor::TaggerProcessor( std::string cfg_path )
      : larcv::ProcessBase( "TaggerProcessor" ),
      _cfg_path(cfg_path),
      _larlite_io(nullptr)
    {
      
    }

    TaggerProcessor::~TaggerProcessor() {
      if ( _larlite_io ) {
        _larlite_io->close();
        delete _larlite_io;
      }
    }
    
    void TaggerProcessor::configure(const larcv::PSet&) {
      
    }

    void TaggerProcessor::add_larlite_input( std::string llinput ) {
      if ( !_larlite_io ) {
        LARCV_CRITICAL() << "Adding larlite file before creating manager. Must run configure first." << std::endl;
        throw std::runtime_error( "adding larlite file before creating manager" );
      }
      _larlite_io->add_in_filename( llinput );
    }
    
    bool TaggerProcessor::process(larcv::IOManager& io) {

      // sync larlite
      if ( !_larlite_io ) {
        LARCV_CRITICAL() << "Processing entry before creating manager. Must run configure first." << std::endl;
        throw std::runtime_error( "adding larlite file before creating manager" );
      }
      _larlite_io->syncEntry( io );
      
      // prepare inputs
      InputPayload input;

      
      
    }
    
    void TaggerProcessor::finalize() {
      
    }

    
  }
}
