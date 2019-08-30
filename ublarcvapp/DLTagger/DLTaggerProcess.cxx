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
    _input_adc_producer      = pset.get<std::string>("InputADCproducer");
    _input_chstatus_producer = pset.get<std::string>("InputChStatusProducer");
    _input_mask_producer     = pset.get<std::string>("InputMRCNNproducer");

    _has_mcinstance_img      = pset.get<bool>("HasMCInstanceImage",false);
    _input_instance_image    = pset.get<std::string>("InputInstanceProducer","instance");

    // larlite inputs
    _input_opflash_producer  = pset.get<std::string>("InputOpFlashProducer");
    _input_ophit_producer    = pset.get<std::string>("InputOpHitProducer");
    _larlite_files           = pset.get<std::vector<std::string> >("LArLiteInputFiles");

    // larcv output
    _output_tagged_image       = pset.get<std::string>("OutputTaggerImage","thrumu");
    _output_notcosmic_image    = pset.get<std::string>("OutputNotCosmicImage","notcosmic");
    _output_cosmic_clusters    = pset.get<std::string>("OutputCosmicPixelCluster","thrumupixels");
    _output_notcosmic_clusters = pset.get<std::string>("OutputNotCosmicPixelCluster", "notcosmic");    
    _output_croi               = pset.get<std::string>("OutputCROI","croi");
    _output_croi_merged        = pset.get<std::string>("OutputMergedCROI","croimerged");

    // larlite output
    _output_larlite_file     = pset.get<std::string>("OutputLArLiteFile","dltaggerout_larlite.root");
    _output_tracks           = pset.get<std::string>("OutputTracks","thrumu");

    // BNB window
    _inbeam_win_start_tick = pset.get<int>("InBeamWindowStartTick");//215;
    _inbeam_win_end_tick   = pset.get<int>("InBeamWindowEndTick");//345;

    _ana_tree = nullptr;
  }

  void DLTaggerProcess::initialize() {
    _larlite_io = new ublarcvapp::LArliteManager( larlite::storage_manager::kBOTH, "DLTaggerLarlite" );
    for ( auto const& input_larlite : _larlite_files ) {
      _larlite_io->add_in_filename( input_larlite );
    }
    _larlite_io->set_out_filename( _output_larlite_file );
    _larlite_io->open();
    
    setupAnaTree();
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

    larcv::EventImage2D* ev_mcinstance = nullptr;
    if ( _has_mcinstance_img ) {
      ev_mcinstance = (larcv::EventImage2D*)mgr.get_data(larcv::kProductImage2D, _input_instance_image);
    }

    LARCV_NORMAL() << "Processing entry[" << mgr.current_entry() << "] "
                   << "rse=(" << mgr.event_id().run() << "," << mgr.event_id().subrun() << "," << mgr.event_id().event() << ")"
                   << std::endl;

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

    if ( _has_mcinstance_img ) {
      m_tagger.recoTruthMatching( ev_mcinstance->Image2DArray() );
    }
    
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

    fillAnaVars();
    
    return true;
  }

  void DLTaggerProcess::setupAnaTree() {
    if ( !has_ana_file() ) {
      LARCV_WARNING() << "No analysis tree defined" << std::endl;
      return;
    }

    LARCV_NORMAL() << "Setup analysis tree" << std::endl;
    ana_file().cd();
    _ana_tree = new TTree("dltaggervars", "DL Tagger Selection Variables");
    _ana_tree->Branch( "numclusters", &_num_clusters, "numclusters/I" );
    _ana_tree->Branch( "numpixels_plane0", &_numpixels_plane0 );
    _ana_tree->Branch( "numpixels_plane1", &_numpixels_plane1 );
    _ana_tree->Branch( "numpixels_plane2", &_numpixels_plane2 );
    _ana_tree->Branch( "dwall_outermost",  &_dwall_outermost );
    _ana_tree->Branch( "dwall_innermost",  &_dwall_innermost );
    _ana_tree->Branch( "astar_complete",   &_astar_complete );
    _ana_tree->Branch( "dtick_outoftime",  &_dtick_outoftime );
    _ana_tree->Branch( "frac_out_croi_tot", &_frac_in_croi_total );
    _ana_tree->Branch( "frac_out_croi_plane0", &_frac_in_croi_plane0);
    _ana_tree->Branch( "frac_out_croi_plane1", &_frac_in_croi_plane1);
    _ana_tree->Branch( "frac_out_croi_plane2", &_frac_in_croi_plane2);
    _ana_tree->Branch( "nufrac_tot", &_nufrac_total );
    _ana_tree->Branch( "nufrac_plane0", &_nufrac_plane0);
    _ana_tree->Branch( "nufrac_plane1", &_nufrac_plane1);
    _ana_tree->Branch( "nufrac_plane2", &_nufrac_plane2);
    _ana_tree->Branch( "outermost_endpt", &_outermost_endpt_v );
    _ana_tree->Branch( "innermost_endpt", &_innermost_endpt_v );    
  }

  void DLTaggerProcess::clearAnaVars() {
    _num_clusters = 0;
    _astar_complete.clear();
    _dwall_outermost.clear();
    _dwall_innermost.clear();
    _dtick_outoftime.clear();
    _frac_in_croi_plane0.clear();
    _frac_in_croi_plane1.clear();
    _frac_in_croi_plane2.clear();
    _frac_in_croi_total.clear();
    _nufrac_plane0.clear();
    _nufrac_plane1.clear();
    _nufrac_plane2.clear();
    _nufrac_total.clear();
    _numpixels_plane0.clear();
    _numpixels_plane1.clear();
    _numpixels_plane2.clear();
    _outermost_endpt_v.clear();
    _innermost_endpt_v.clear();
  }

  void DLTaggerProcess::fillAnaVars() {
    if ( !has_ana_file() ) {
      LARCV_DEBUG() << "No anatree to fill" << std::endl;
      return;
    }

    LARCV_DEBUG() << "Filling ana tree" << std::endl;
    clearAnaVars();
    auto const& sel_vars_v = m_tagger.getSelectionVars();
    _num_clusters = sel_vars_v.size();
    for ( auto const& vars : sel_vars_v ) {
      _astar_complete.push_back(  vars.astar_complete );
      _dwall_outermost.push_back( vars.dwall_outermost );
      _dwall_innermost.push_back( vars.dwall_innermost );
      _dtick_outoftime.push_back( vars.dtick_outoftime );
      _frac_in_croi_total.push_back( vars.total_frac );
      _frac_in_croi_plane0.push_back( vars.frac_per_plane[0] );
      _frac_in_croi_plane1.push_back( vars.frac_per_plane[1] );
      _frac_in_croi_plane2.push_back( vars.frac_per_plane[2] );
      _nufrac_total.push_back( vars.total_nufrac );
      _nufrac_plane0.push_back( vars.nufrac_per_plane[0] );
      _nufrac_plane1.push_back( vars.nufrac_per_plane[1] );
      _nufrac_plane2.push_back( vars.nufrac_per_plane[2] );      
      _numpixels_plane0.push_back( vars.numpixels[0] );
      _numpixels_plane1.push_back( vars.numpixels[1] );
      _numpixels_plane2.push_back( vars.numpixels[2] );
      _outermost_endpt_v.push_back( vars.outermost_endpt_tyz );
      _innermost_endpt_v.push_back( vars.innermost_endpt_tyz );
    }
    _ana_tree->Fill();
  }

  void DLTaggerProcess::finalize() {
    if ( has_ana_file() && _ana_tree ) {
      _ana_tree->Write();
    }
    _larlite_io->close();
    delete _larlite_io;
  }
  
}
}
