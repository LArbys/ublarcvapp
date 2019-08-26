#include "DLTaggerProcess.h"

#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"
#include "larcv/core/DataFormat/EventClusterMask.h"
#include "larcv/core/DataFormat/EventROI.h"

namespace ublarcvapp {
namespace dltagger {

  static DLTaggerProcessFactory __global_DLTaggerProcessFactory__;
  
  void DLTaggerProcess::configure( const larcv::PSet& pset ) {

    // larcv inputs
    _input_adc_producer      = "wire";
    _input_chstatus_producer = "wire";
    _input_mask_producer     = "mrcnn_masks";

    // larlite inputs
    _input_opflash_producer  = "simpleFlashBeam";
    _input_ophit_producer    = "opHitBeam";
    _larlite_files           = std::vector<std::string>();
    _larlite_files.push_back("testset1/larlite_opreco.root");

    // larcv output
    _output_tagged_image       = "thrumu";
    _output_notcosmic_image    = "notcosmic";    
    _output_cosmic_clusters    = "thrumupixels";
    _output_notcosmic_clusters = "notcosmic";    
    _output_croi               = "croi";
    _output_croi_merged        = "croimerged";

    // larlite output
    _output_larlite_file     = "dltaggerout_larlite.root";
    _output_tracks           = "thrumu";

    // BNB window
    _inbeam_win_start_tick = 215;
    _inbeam_win_end_tick   = 345;
  }

  void DLTaggerProcess::initialize() {
    _larlite_io = new ublarcvapp::LArliteManager( larlite::storage_manager::kBOTH, "DLTaggerLarlite" );
    for ( auto const& input_larlite : _larlite_files ) {
      _larlite_io->add_in_filename( input_larlite );
    }
    _larlite_io->set_out_filename( _output_larlite_file );
    _larlite_io->open();
  }

  bool DLTaggerProcess::process( larcv::IOManager& mgr ) {

    m_tagger.set_verbosity( logger().level() );
    _larlite_io->set_verbosity( logger().level() );
    
    // GET LARCV INPUT
    larcv::EventImage2D* ev_wire
      = (larcv::EventImage2D*) mgr.get_data(larcv::kProductImage2D,  _input_adc_producer );
    
    larcv::EventChStatus* ev_chstatus
      = (larcv::EventChStatus*)mgr.get_data(larcv::kProductChStatus, _input_chstatus_producer );

    larcv::EventClusterMask* ev_clustermask
      = (larcv::EventClusterMask*)mgr.get_data(larcv::kProductClusterMask, _input_mask_producer );

    // GET LARLITE INPUT
    _larlite_io->syncEntry( mgr );
    larlite::event_opflash* ev_opflash
      = (larlite::event_opflash*)_larlite_io->get_data(larlite::data::kOpFlash, _input_opflash_producer );
    std::vector< larlite::opflash > intime_opflash_v;
    for ( auto const& opflash : *ev_opflash ) {
      float t_usec = opflash.Time(); // from start of trigger
      // assuming trigger and beam window readout start the same -- they are not!
      int tick = t_usec/(0.015625); // 64 Mhz clock
      if ( tick>=_inbeam_win_start_tick && tick<=_inbeam_win_end_tick ) {
        intime_opflash_v.push_back( opflash );
      }
      LARCV_DEBUG() << "opflash: t_usec = " << t_usec << " t_tick=" << tick
                    << " window=[" << _inbeam_win_start_tick << "," << _inbeam_win_end_tick << "]"
                    << std::endl;
    }
    LARCV_DEBUG() << "number of intime flashes found: " << intime_opflash_v.size() << std::endl;

    // RUN PRECUTS

    // RUN TAGGER    
    m_tagger.runTagger( ev_wire->Image2DArray(),
                        *ev_chstatus,
                        intime_opflash_v,
                        ev_clustermask->as_vector() );


    // OUTPUTS

    // precuts
    
    // croi
    
    // tagger
    larcv::EventImage2D* evout_tagged
      = (larcv::EventImage2D*)mgr.get_data( larcv::kProductImage2D, _output_tagged_image );

    larcv::EventImage2D* evout_notcosmic
      = (larcv::EventImage2D*)mgr.get_data( larcv::kProductImage2D, _output_notcosmic_image );

    larcv::EventPixel2D* evout_cosmic_clusters
      = (larcv::EventPixel2D*)mgr.get_data( larcv::kProductPixel2D, _output_cosmic_clusters );
    larcv::EventPixel2D* evout_notcosmic_clusters
      = (larcv::EventPixel2D*)mgr.get_data( larcv::kProductPixel2D, _output_notcosmic_clusters );

    larcv::EventROI* ev_croi
      = (larcv::EventROI*)mgr.get_data( larcv::kProductROI, _output_croi );
    larcv::EventROI* ev_croi_merged
      = (larcv::EventROI*)mgr.get_data( larcv::kProductROI, _output_croi_merged );

    std::vector<larcv::Image2D> tagged_v;
    std::vector<larcv::Image2D> notcosmic_v;    
    m_tagger.transferImages( tagged_v, notcosmic_v );
    evout_tagged->Emplace( std::move(tagged_v) );
    evout_notcosmic->Emplace( std::move(notcosmic_v) );    
    
    m_tagger.transferPixelClusters( *evout_cosmic_clusters, *evout_notcosmic_clusters );

    m_tagger.transferCROI( *ev_croi, *ev_croi_merged );
    
    return true;
  }

  void DLTaggerProcess::finalize() {
    _larlite_io->close();
    delete _larlite_io;
  }
  
}
}
