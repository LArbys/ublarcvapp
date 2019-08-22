#include "DLTaggerProcess.h"

#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"
#include "larcv/core/DataFormat/EventClusterMask.h"


namespace ublarcvapp {
namespace dltagger {

  static DLTaggerProcessFactory __global_DLTaggerProcessFactory__;
  
  void DLTaggerProcess::configure( const larcv::PSet& pset ) {

    // larcv inputs
    _input_adc_producer      = "wire";
    _input_chstatus_producer = "wire";
    _input_mask_producer     = "mrcnn_masks";

    // larlite inputs
    _input_opflash_producer  = "SimpleFlashBeam";
    _input_ophit_producer    = "OpHitBeam";
    _larlite_files           = std::vector<std::string>();

    // larcv output
    _output_tagged_image     = "thrumu";
    _output_pixel_clusters   = "thrumupixels";

    // larlite output
    _output_larlite_file     = "output_larlite.root";
    _output_tracks           = "thrumu";
  }

  void DLTaggerProcess::initialize() {

  }

  bool DLTaggerProcess::process( larcv::IOManager& mgr ) {

    // GET INPUT
    larcv::EventImage2D* ev_wire
      = (larcv::EventImage2D*) mgr.get_data(larcv::kProductImage2D,  _input_adc_producer );
    
    larcv::EventChStatus* ev_chstatus
      = (larcv::EventChStatus*)mgr.get_data(larcv::kProductChStatus, _input_chstatus_producer );

    larcv::EventClusterMask* ev_clustermask
      = (larcv::EventClusterMask*)mgr.get_data(larcv::kProductClusterMask, _input_mask_producer );


    // RUN PRECUTS

    // RUN CROI

    // RUN TAGGER    
    m_tagger.runTagger( ev_wire->Image2DArray(),
                        *ev_chstatus,
                        ev_clustermask->as_vector() );


    // OUTPUTS

    // precuts
    
    // croi
    
    // tagger
    larcv::EventImage2D* evout_tagged
      = (larcv::EventImage2D*)mgr.get_data( larcv::kProductImage2D, _output_tagged_image );

    larcv::EventPixel2D* evout_pixel_clusters
      = (larcv::EventPixel2D*)mgr.get_data( larcv::kProductPixel2D, _output_pixel_clusters );

    std::vector<larcv::Image2D> tagged_v;
    m_tagger.transferImages( tagged_v );
    evout_tagged->Emplace( std::move(tagged_v) );
    m_tagger.transferPixelClusters( *evout_pixel_clusters );

    
    
    return true;
  }

  void DLTaggerProcess::finalize() {
  }
  
}
}
