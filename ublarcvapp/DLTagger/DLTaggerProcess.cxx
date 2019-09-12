#include "DLTaggerProcess.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

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
    _input_mcpart_larcvtruth = pset.get<std::string>("InputLArCVMCparticleProducer","segment");

    // larlite inputs
    _input_opflash_producer  = pset.get<std::string>("InputOpFlashProducer");
    _input_ophit_producer    = pset.get<std::string>("InputOpHitProducer");
    _larlite_files           = pset.get<std::vector<std::string> >("LArLiteInputFiles");

    // larcv output
    _output_tagged_image       = pset.get<std::string>("OutputTaggerImage","thrumu");
    _output_notcosmic_image    = pset.get<std::string>("OutputNotCosmicImage","notcosmic");
    _output_cosmic_clusters    = pset.get<std::string>("OutputCosmicPixelCluster","thrumupixels");
    _output_notcosmic_clusters = pset.get<std::string>("OutputNotCosmicPixelCluster", "notcosmic");
    _output_allreco_clusters   = pset.get<std::string>("OutputAllRecoCluster","allreco");        
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


    LARCV_NORMAL() << "Processing entry[" << mgr.current_entry() << "] "
                   << "rse=(" << mgr.event_id().run() << "," << mgr.event_id().subrun() << "," << mgr.event_id().event() << ")"
                   << std::endl;
    
    larcv::EventImage2D* ev_mcinstance = nullptr;
    larcv::EventROI*     ev_mcpartroi  = nullptr;
    if ( _has_mcinstance_img ) {
      LARCV_INFO() << "Has MC information: mcinstance and particle ROI" << std::endl;
      ev_mcinstance = (larcv::EventImage2D*)mgr.get_data(larcv::kProductImage2D, _input_instance_image);
      ev_mcpartroi  = (larcv::EventROI*)mgr.get_data(larcv::kProductROI, _input_mcpart_larcvtruth );
      LARCV_INFO() << "MC Particle ROI:" << std::endl;
      for (auto const& roi : ev_mcpartroi->ROIArray() ) {
        LARCV_INFO() << roi.dump() << std::endl;
      }
    }


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
    larcv::EventPixel2D* evout_allreco_clusters
      = (larcv::EventPixel2D*)mgr.get_data( larcv::kProductPixel2D, _output_allreco_clusters );

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

    // copy contents of pixel clusters to all reco
    for ( size_t p=0; p< 3; p++ ) {
      auto const& pixcluster_v     = evout_cosmic_clusters->Pixel2DClusterArray(p);
      auto const& pixclustermeta_v = evout_cosmic_clusters->ClusterMetaArray(p);
      for ( size_t ic=0; ic<pixcluster_v.size(); ic++ ) 
        evout_allreco_clusters->Append( p, pixcluster_v[ic], pixclustermeta_v[ic] );
      auto const& pixcluster2_v     = evout_notcosmic_clusters->Pixel2DClusterArray(p);
      auto const& pixclustermeta2_v = evout_notcosmic_clusters->ClusterMetaArray(p);
      for ( size_t ic=0; ic<pixcluster2_v.size(); ic++ ) 
        evout_allreco_clusters->Append( p, pixcluster2_v[ic], pixclustermeta2_v[ic] );
    }
    
    m_tagger.transferCROI( *ev_croi, *ev_croi_merged );

    fillAnaVars( *ev_wire, ev_mcinstance, ev_mcpartroi, *evout_tagged, *evout_notcosmic  );
    
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
    _ana_tree->Branch( "frac_wholeimg_cosmictag", &_frac_wholeimg_cosmictag, "frac_wholeimg_cosmictag/F" );
    _ana_tree->Branch( "frac_wholeimg_notcosmictag", &_frac_wholeimg_notcosmictag, "frac_wholeimg_notcosmictag/F" );
    _ana_tree->Branch( "frac_wholeimg_alltags", &_frac_wholeimg_alltags, "frac_wholeimg_alltags/F" );
    _ana_tree->Branch( "frac_nupixel_cosmictag",     &_frac_wholeimg_nu_cosmictag,    "frac_nupixel_cosmictag/F" );
    _ana_tree->Branch( "frac_nupixel_notcosmictag",  &_frac_wholeimg_nu_notcosmictag, "frac_nupixel_notcosmictag/F" );
    _ana_tree->Branch( "frac_vtxpixel_cosmictag",    &_frac_wholeimg_vtx_cosmictag,    "frac_vtxpixel_cosmictag/F" );
    _ana_tree->Branch( "frac_vtxpixel_notcosmictag", &_frac_wholeimg_vtx_notcosmictag, "frac_vtxpixel_notcosmictag/F" );
    
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
    _frac_wholeimg_cosmictag = 0;
    _frac_wholeimg_notcosmictag = 0;
    _frac_wholeimg_alltags = 0;
    _frac_wholeimg_nu_cosmictag = 0;
    _frac_wholeimg_nu_notcosmictag = 0;
    _frac_wholeimg_vtx_cosmictag = 0;
    _frac_wholeimg_vtx_notcosmictag = 0;
  }

  void DLTaggerProcess::fillAnaVars( const larcv::EventImage2D& ev_wholeview,
                                     const larcv::EventImage2D* ev_mcinstance,
                                     const larcv::EventROI* ev_mcparticle,                                     
                                     const larcv::EventImage2D& evout_tagged,
                                     const larcv::EventImage2D& evout_notcosmic ) {

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

    // if have instance image, we can evaluate fraction of cosmic pixels tagged
    if ( _has_mcinstance_img ) {
      // MC ANALYSIS

      // get vertex: find neutrino
      std::vector<double> vertex_xyzt(3,0.0);
      double tstart = 0;
      bool foundvertex = false;
      for ( auto const& roi : ev_mcparticle->ROIArray() ) {
        int pdgcode = abs(roi.PdgCode());
        if ( pdgcode==12 || pdgcode==14 || pdgcode==16 ) {
          vertex_xyzt[0] = roi.X();
          vertex_xyzt[1] = roi.Y();
          vertex_xyzt[2] = roi.Z();
          tstart = roi.T();
          foundvertex = true;
          break;
        }
      }
      // vertex image coordinates
      bool inimage = true;
      if ( vertex_xyzt[0]<0 || vertex_xyzt[0]>larutil::Geometry::GetME()->DetHalfWidth()*2 )
        inimage = false;
      if ( vertex_xyzt[1]<-larutil::Geometry::GetME()->DetHalfHeight() || vertex_xyzt[1]>larutil::Geometry::GetME()->DetHalfHeight() )
        inimage = false;
      if ( vertex_xyzt[2]<0 || vertex_xyzt[2]>larutil::Geometry::GetME()->DetLength() )
        inimage = false;

      float tick = 3200 + vertex_xyzt[0]/larutil::LArProperties::GetME()->DriftVelocity()/0.5;
      int true_row = ev_wholeview.Image2DArray().front().meta().row(tick);
      std::vector<int> true_wire_v(3,0);
      std::vector<int> true_col_v(3,0);
      
      for ( size_t p=0; p<3; p++ ) {
        true_wire_v[p] = larutil::Geometry::GetME()->WireCoordinate( vertex_xyzt, p );
        if ( true_wire_v[p]<0 ) true_wire_v[p] = 0;
        if ( true_wire_v[p]>=larutil::Geometry::GetME()->Nwires(p) )
          true_wire_v[p] = (int)larutil::Geometry::GetME()->Nwires(p)-1;
      }

      
      float total_iscosmic = 0.;
      float total_nupix = 0;
      float total_vtxpix = 0;
      _frac_wholeimg_cosmictag = 0;
      _frac_wholeimg_notcosmictag = 0;
      _frac_wholeimg_alltags = 0;
      _frac_wholeimg_nu_cosmictag = 0;
      _frac_wholeimg_nu_notcosmictag = 0;
      _frac_wholeimg_vtx_cosmictag = 0;
      _frac_wholeimg_vtx_notcosmictag = 0;
      for ( size_t p=0; p<ev_wholeview.Image2DArray().size(); p++ ) {
        auto const& adcimg       = ev_wholeview.Image2DArray().at(p);
        auto const& instanceimg  = ev_mcinstance->Image2DArray().at(p);
        auto const& cosmictagimg = evout_tagged.Image2DArray().at(p);
        auto const& notcosmicimg = evout_notcosmic.Image2DArray().at(p);
        int npixthresh = 0; // pixels above threshold
        int npixcosmic = 0; // pixels above threshold + cosmic (i.e. zero instance id)
        int npixtagged = 0; // pixel above threshold + true cosmic + tagged
        int npixrecoed = 0; // pixel above threshold + true cosmic + (tagged cosmic or not-cosmic)
        int npixnu     = 0;
        int npixnu_tag = 0;
        int npixnu_not = 0;
        int npixvtx    = 0;
        int npixvtx_tag = 0;
        int npixvtx_not = 0;
        for ( size_t c=0; c<adcimg.meta().cols(); c++ ) {
          for ( size_t r=0; r<adcimg.meta().rows(); r++ ) {
            if ( adcimg.pixel(r,c)<10.0 ) continue;
            npixthresh++;
            if ( instanceimg.pixel(r,c)>0 ) {
              // on neutrino pixel
              npixnu++;
              if ( cosmictagimg.pixel(r,c)>0 )
                npixnu_tag++;
              if ( notcosmicimg.pixel(r,c)>0 )
                npixnu_not++;

              if ( inimage ) {
                int drow = abs( true_row - (int)r );
                int dcol = abs( true_col_v[p] - (int)c );
                
                if ( drow<=20 || dcol<=20 ) {
                  npixvtx++;
                  if ( cosmictagimg.pixel(r,c)>0 )
                    npixvtx_tag++;
                  if ( notcosmicimg.pixel(r,c)>0 )
                    npixvtx_not++;
                }
              }
            }
            else {
              // on cosmic overlay pixel
              npixcosmic++;
              if ( cosmictagimg.pixel(r,c)>0 )
                npixtagged++;
              if ( cosmictagimg.pixel(r,c)>0 || notcosmicimg.pixel(r,c)>0 )
                npixrecoed++;
            }
          }
        }
        
        total_iscosmic += (float)npixcosmic;
        total_nupix    += (float)npixnu;
        total_vtxpix   += (float)npixvtx;
        
        _frac_wholeimg_cosmictag += (float)npixtagged;
        _frac_wholeimg_alltags   += (float)npixrecoed;

        _frac_wholeimg_nu_cosmictag += (float)npixnu_tag;
        _frac_wholeimg_nu_notcosmictag += (float)npixnu_not;

        _frac_wholeimg_vtx_cosmictag += (float)npixvtx_tag;
        _frac_wholeimg_vtx_notcosmictag += (float)npixvtx_not;        
      }//end of plane loop

      if ( total_iscosmic>0 ) {
        _frac_wholeimg_cosmictag /= total_iscosmic;
        _frac_wholeimg_alltags   /= total_iscosmic;
      }
      if ( total_nupix>0 ) {
        _frac_wholeimg_nu_cosmictag /= total_nupix;
        _frac_wholeimg_nu_notcosmictag /= total_nupix;
      }
      if ( total_vtxpix>0 ) {
        _frac_wholeimg_vtx_cosmictag /= total_vtxpix;
        _frac_wholeimg_vtx_notcosmictag /= total_vtxpix;
      }
      
      LARCV_INFO() << "fraction of cosmic pixels tagged as cosmic: " << _frac_wholeimg_cosmictag << std::endl;
      LARCV_INFO() << "fraction of cosmic pixels reco'd: " << _frac_wholeimg_alltags << std::endl;
      LARCV_INFO() << "fraction of neutrino pixels tagged as cosmic: " << _frac_wholeimg_nu_cosmictag << std::endl;
      LARCV_INFO() << "fraction of neutrino pixels tagged as not-cosmic: " << _frac_wholeimg_nu_notcosmictag << std::endl;
      LARCV_INFO() << "fraction of vertex pixels tagged as cosmic: " << _frac_wholeimg_vtx_cosmictag << std::endl;
      LARCV_INFO() << "fraction of vertex pixels tagged as not-cosmic: " << _frac_wholeimg_vtx_notcosmictag << std::endl;
    }// end of has mc instance

    
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
