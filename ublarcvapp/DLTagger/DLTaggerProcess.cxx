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
    _output_tracks           = pset.get<std::string>("OutputTracks","dltagger");

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
    std::cout << "sync" << std::endl;
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

    // croi    
    larcv::EventROI* ev_croi
      = (larcv::EventROI*)mgr.get_data( larcv::kProductROI, _output_croi );
    larcv::EventROI* ev_croi_merged
      = (larcv::EventROI*)mgr.get_data( larcv::kProductROI, _output_croi_merged );

    m_tagger.transferCROI( *ev_croi, *ev_croi_merged );
    
    // tagger images
    larcv::EventImage2D* evout_tagged
      = (larcv::EventImage2D*)mgr.get_data( larcv::kProductImage2D, _output_tagged_image );

    larcv::EventImage2D* evout_notcosmic
      = (larcv::EventImage2D*)mgr.get_data( larcv::kProductImage2D, _output_notcosmic_image );

    std::vector<larcv::Image2D> tagged_v;
    std::vector<larcv::Image2D> notcosmic_v;    
    m_tagger.transferImages( tagged_v, notcosmic_v );

    evout_tagged->Emplace( std::move(tagged_v) );
    evout_notcosmic->Emplace( std::move(notcosmic_v) );    

    // pixel clusters
    larcv::EventPixel2D* evout_cosmic_clusters
      = (larcv::EventPixel2D*)mgr.get_data( larcv::kProductPixel2D, _output_cosmic_clusters );
    larcv::EventPixel2D* evout_notcosmic_clusters
      = (larcv::EventPixel2D*)mgr.get_data( larcv::kProductPixel2D, _output_notcosmic_clusters );
    larcv::EventPixel2D* evout_allreco_clusters
      = (larcv::EventPixel2D*)mgr.get_data( larcv::kProductPixel2D, _output_allreco_clusters );
    
    m_tagger.transferPixelClusters( *evout_cosmic_clusters, *evout_notcosmic_clusters );

    // copy contents of pixel clusters to all reco
    for ( size_t p=0; p< 3; p++ ) {
      auto const& pixcluster_v     = evout_cosmic_clusters->Pixel2DClusterArray(p);
      if ( pixcluster_v.size()>0 ) {
        auto const& pixclustermeta_v = evout_cosmic_clusters->ClusterMetaArray(p);
        for ( size_t ic=0; ic<pixcluster_v.size(); ic++ ) 
          evout_allreco_clusters->Append( p, pixcluster_v[ic], pixclustermeta_v[ic] );
      }
      auto const& pixcluster2_v     = evout_notcosmic_clusters->Pixel2DClusterArray(p);
      if ( pixcluster2_v.size()>0 ) {
        auto const& pixclustermeta2_v = evout_notcosmic_clusters->ClusterMetaArray(p);
        for ( size_t ic=0; ic<pixcluster2_v.size(); ic++ ) 
          evout_allreco_clusters->Append( p, pixcluster2_v[ic], pixclustermeta2_v[ic] );
      }
    }

    // larlite tracks
    larlite::event_track* evout_track_cosmic
      = (larlite::event_track*)_larlite_io->get_data( larlite::data::kTrack, _output_tracks+"_cosmic" );
    larlite::event_track* evout_track_notcosmic
      = (larlite::event_track*)_larlite_io->get_data( larlite::data::kTrack, _output_tracks+"_notcosmic" );
    larlite::event_track* evout_track_allreco
      = (larlite::event_track*)_larlite_io->get_data( larlite::data::kTrack, _output_tracks+"_allreco" );
    m_tagger.getCosmicTracks( *evout_track_cosmic );
    m_tagger.getNotCosmicTracks( *evout_track_notcosmic );
    m_tagger.getAllGoodTracks( *evout_track_allreco );
      
    fillAnaVars( *ev_wire, ev_mcinstance, ev_mcpartroi, *evout_tagged, *evout_notcosmic  );

    std::cout << "next" << std::endl;    
    _larlite_io->next_event(); // store data
    
    return true;
  }

  void DLTaggerProcess::setupAnaTree() {
    if ( !has_ana_file() ) {
      LARCV_WARNING() << "No analysis tree defined" << std::endl;
      return;
    }

    LARCV_NORMAL() << "Setup RECO analysis tree" << std::endl;
    ana_file().cd();
    _ana_tree = new TTree("dltaggerana", "DL Tagger Ana Variables");

    // CLUSTER RECO VARIABLES    
    _ana_tree->Branch( "numclusters", &_num_clusters, "numclusters/I" );
    _ana_tree->Branch( "clust_numpixels_plane0", &_numpixels_plane0 );
    _ana_tree->Branch( "clust_numpixels_plane1", &_numpixels_plane1 );
    _ana_tree->Branch( "clust_numpixels_plane2", &_numpixels_plane2 );
    _ana_tree->Branch( "clust_dwall_outermost",  &_dwall_outermost );
    _ana_tree->Branch( "clust_dwall_innermost",  &_dwall_innermost );
    _ana_tree->Branch( "clust_astar_complete",   &_astar_complete );
    _ana_tree->Branch( "clust_dtick_outoftime",  &_dtick_outoftime );
    _ana_tree->Branch( "clust_frac_out_croi_tot", &_frac_in_croi_total );
    _ana_tree->Branch( "clust_frac_out_croi_plane0", &_frac_in_croi_plane0);
    _ana_tree->Branch( "clust_frac_out_croi_plane1", &_frac_in_croi_plane1);
    _ana_tree->Branch( "clust_frac_out_croi_plane2", &_frac_in_croi_plane2);
    _ana_tree->Branch( "clust_outermost_endpt", &_outermost_endpt_v );
    _ana_tree->Branch( "clust_innermost_endpt", &_innermost_endpt_v );

    // CLUSTER MC VARIABLES
    _ana_tree->Branch( "clust_mc_nufrac_tot", &_nufrac_total );
    _ana_tree->Branch( "clust_mc_nufrac_plane0", &_nufrac_plane0);
    _ana_tree->Branch( "clust_mc_nufrac_plane1", &_nufrac_plane1);
    _ana_tree->Branch( "clust_mc_nufrac_plane2", &_nufrac_plane2);

    // RECO EVENT VARIABLES
    _ana_tree->Branch( "plane_num_input_mrcnnmasks",       &_num_input_mrcnnmasks);
    _ana_tree->Branch( "plane_num_used_mrcnnmasks",        &_num_used_mrcnnmasks );        
    _ana_tree->Branch( "plane_frac_wholeimg_cosmictag",    &_frac_wholeimg_cosmictag );
    _ana_tree->Branch( "plane_frac_wholeimg_notcosmictag", &_frac_wholeimg_notcosmictag );
    _ana_tree->Branch( "plane_frac_wholeimg_alltags",      &_frac_wholeimg_alltags );

    // VARIABLES FOR MC EVENT ANALYSIS
    _ana_tree->Branch( "plane_mc_frac_cosmicpixel_cosmictag",    &_frac_wholeimg_cosmic_cosmictag );
    _ana_tree->Branch( "plane_mc_frac_cosmicpixel_notcosmictag", &_frac_wholeimg_cosmic_notcosmictag );    
    _ana_tree->Branch( "plane_mc_frac_nupixel_cosmictag",        &_frac_wholeimg_nu_cosmictag );
    _ana_tree->Branch( "plane_mc_frac_nupixel_notcosmictag",     &_frac_wholeimg_nu_notcosmictag );
    _ana_tree->Branch( "plane_mc_frac_vtxpixel_cosmictag",       &_frac_wholeimg_vtx_cosmictag );
    _ana_tree->Branch( "plane_mc_frac_vtxpixel_notcosmictag",    &_frac_wholeimg_vtx_notcosmictag );
    
  }

  void DLTaggerProcess::clearAnaVars() {

    // PER CLUSTER VARS
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

    // RECO PER PLANE VARS
    _num_input_mrcnnmasks.clear();
    _num_used_mrcnnmasks.clear();
    _frac_wholeimg_cosmictag.clear();
    _frac_wholeimg_notcosmictag.clear();
    _frac_wholeimg_alltags.clear();

    // MC PER PLANE VARS
    _frac_wholeimg_cosmic_cosmictag.clear();
    _frac_wholeimg_cosmic_notcosmictag.clear();
    _frac_wholeimg_nu_cosmictag.clear();
    _frac_wholeimg_nu_notcosmictag.clear();
    _frac_wholeimg_vtx_cosmictag.clear();
    _frac_wholeimg_vtx_notcosmictag.clear();
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

    // CLUSTER LEVEL RECO+MC VARS
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
      //_nufrac_total.push_back( vars.total_nufrac );
      //_nufrac_plane0.push_back( vars.nufrac_per_plane[0] );
      //_nufrac_plane1.push_back( vars.nufrac_per_plane[1] );
      //_nufrac_plane2.push_back( vars.nufrac_per_plane[2] );      
      _numpixels_plane0.push_back( vars.numpixels[0] );
      _numpixels_plane1.push_back( vars.numpixels[1] );
      _numpixels_plane2.push_back( vars.numpixels[2] );
      _outermost_endpt_v.push_back( vars.outermost_endpt_tyz );
      _innermost_endpt_v.push_back( vars.innermost_endpt_tyz );
    }

    size_t nplanes = ev_wholeview.Image2DArray().size();
    
    // event level reco variables
    float total_above_thresh = 0.;    
    _frac_wholeimg_cosmictag.resize(nplanes+1,0);
    _frac_wholeimg_notcosmictag.resize(nplanes+1,0);
    _frac_wholeimg_alltags.resize(nplanes+1,0);
    for ( size_t p=0; p<nplanes; p++ ) {

      auto const& adcimg       = ev_wholeview.Image2DArray().at(p);
      auto const& cosmictagimg = evout_tagged.Image2DArray().at(p);
      auto const& notcosmicimg = evout_notcosmic.Image2DArray().at(p);

      int npixthresh = 0; // above thresh pixels this plane
      
      for ( size_t c=0; c<adcimg.meta().cols(); c++ ) {
        for ( size_t r=0; r<adcimg.meta().rows(); r++ ) {
          if ( adcimg.pixel(r,c)<10.0 ) continue;
          npixthresh++;
          total_above_thresh += 1.0;
          
          if ( cosmictagimg.pixel(r,c)>0 ) {
            _frac_wholeimg_cosmictag[p] += 1.0;
            _frac_wholeimg_cosmictag[nplanes] += 1.0;            
          }
          if ( notcosmicimg.pixel(r,c)>0 ) {
            _frac_wholeimg_notcosmictag[p] += 1.0;
            _frac_wholeimg_notcosmictag[nplanes] += 1.0;
          }
          if ( cosmictagimg.pixel(r,c)>0 || notcosmicimg.pixel(r,c)>0 ) {
            _frac_wholeimg_alltags[p] += 1.0;
            _frac_wholeimg_alltags[nplanes] += 1.0;
          }
        }
      }//end of col loop

      if ( npixthresh>0 )  {
        _frac_wholeimg_cosmictag[p]    /= (float)npixthresh;
        _frac_wholeimg_notcosmictag[p] /= (float)npixthresh;
        _frac_wholeimg_alltags[p]      /= (float)npixthresh;
      }
    }//end of plane loop

    if ( total_above_thresh>0 ) {
      _frac_wholeimg_cosmictag[nplanes]    /= total_above_thresh;
      _frac_wholeimg_notcosmictag[nplanes] /= total_above_thresh;
      _frac_wholeimg_alltags[nplanes]      /= total_above_thresh;
    }

    // NUM MASK VARIABLES
    m_tagger.numMaskInfo( _num_input_mrcnnmasks, _num_used_mrcnnmasks );

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


      float total_cosmicpix = 0;
      float total_nupix = 0;
      float total_vtxpix = 0;
      _frac_wholeimg_nu_cosmictag.resize(nplanes+1,0);
      _frac_wholeimg_nu_notcosmictag.resize(nplanes+1,0);
      _frac_wholeimg_cosmic_cosmictag.resize(nplanes+1,0);
      _frac_wholeimg_cosmic_notcosmictag.resize(nplanes+1,0);
      
      if ( foundvertex && inimage ) {
        _frac_wholeimg_vtx_cosmictag.resize(nplanes+1,0);
        _frac_wholeimg_vtx_notcosmictag.resize(nplanes+1,0);
      }
      else {
        _frac_wholeimg_vtx_cosmictag.resize(nplanes+1,-1);
        _frac_wholeimg_vtx_notcosmictag.resize(nplanes+1,-1);
      }
      
      for ( size_t p=0; p<ev_wholeview.Image2DArray().size(); p++ ) {
        auto const& adcimg       = ev_wholeview.Image2DArray().at(p);
        auto const& instanceimg  = ev_mcinstance->Image2DArray().at(p);
        auto const& cosmictagimg = evout_tagged.Image2DArray().at(p);
        auto const& notcosmicimg = evout_notcosmic.Image2DArray().at(p);
        int npixthresh = 0; // pixels above threshold
        int npixcosmic = 0; // pixels above threshold + cosmic (i.e. zero instance id)
        int npixtagged = 0; // pixel above threshold + true cosmic + tagged
        int npixnottagged = 0; // pixel above threshold + true cosmic + (tagged cosmic or not-cosmic)
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
              if ( notcosmicimg.pixel(r,c)>0)
                npixnottagged++;
            }
          } // row loop
        }// col loop
        
        total_cosmicpix += (float)npixcosmic;
        total_nupix     += (float)npixnu;
        total_vtxpix    += (float)npixvtx;

        if ( npixnu>0 ) {
          _frac_wholeimg_nu_cosmictag[p]    = (float)npixnu_tag/(float)npixnu;
          _frac_wholeimg_nu_notcosmictag[p] = (float)npixnu_not/(float)npixnu;
          _frac_wholeimg_nu_cosmictag[nplanes]    += (float)npixnu_tag;
          _frac_wholeimg_nu_notcosmictag[nplanes] += (float)npixnu_not;          
        }

        if ( npixcosmic>0 ) {
          _frac_wholeimg_cosmic_cosmictag[p]    = (float)npixtagged/(float)npixcosmic;
          _frac_wholeimg_cosmic_notcosmictag[p] = (float)npixnottagged/(float)npixcosmic;
          _frac_wholeimg_cosmic_cosmictag[nplanes]    += (float)npixtagged;
          _frac_wholeimg_cosmic_notcosmictag[nplanes] += (float)npixnottagged;                    
        }

        if ( npixvtx ) {
          _frac_wholeimg_vtx_cosmictag[p]    = (float)npixvtx_tag/(float)npixvtx;
          _frac_wholeimg_vtx_notcosmictag[p] = (float)npixvtx_not/(float)npixvtx;
          _frac_wholeimg_vtx_cosmictag[nplanes]    += (float)npixvtx_tag;
          _frac_wholeimg_vtx_notcosmictag[nplanes] += (float)npixvtx_not;
        }
        
      }//end of plane loop
      
      if ( total_cosmicpix>0 ) {
        _frac_wholeimg_cosmic_cosmictag[nplanes]    /= total_cosmicpix;
        _frac_wholeimg_cosmic_notcosmictag[nplanes] /= total_cosmicpix;        
      }
      if ( total_nupix>0 ) {
        _frac_wholeimg_nu_cosmictag[nplanes]    /= total_nupix;
        _frac_wholeimg_nu_notcosmictag[nplanes] /= total_nupix;        
      }
      if ( total_vtxpix>0 ) {
        _frac_wholeimg_vtx_cosmictag[nplanes]    /= total_vtxpix;
        _frac_wholeimg_vtx_notcosmictag[nplanes] /= total_vtxpix;        
      }

      LARCV_INFO() << "Has MC info: total-nu-pixel=" << total_nupix << "  vtx-in-image=" << inimage << " total-vtx=pixel=" << total_vtxpix << std::endl;
      
    }//end of has mc
      
    LARCV_INFO() << "fraction of pixels tagged as cosmic: "
                 << " plane[0]=" << _frac_wholeimg_cosmictag[0]
                 << " plane[1]=" << _frac_wholeimg_cosmictag[1]
                 << " plane[2]=" << _frac_wholeimg_cosmictag[2]
                 << " total=" << _frac_wholeimg_cosmictag[3]
                 << std::endl;

    LARCV_INFO() << "fraction of pixels reco'd: "
                 << " plane[0]=" << _frac_wholeimg_alltags[0]
                 << " plane[1]=" << _frac_wholeimg_alltags[1]
                 << " plane[2]=" << _frac_wholeimg_alltags[2]
                 << " total=" << _frac_wholeimg_alltags[3]
                 << std::endl;

    if ( _has_mcinstance_img ) {

      LARCV_INFO() << "fraction of TRUE COSMIC tagged as [cosmic]: "
                   << " plane[0]=" << _frac_wholeimg_cosmic_cosmictag[0]
                   << " plane[1]=" << _frac_wholeimg_cosmic_cosmictag[1]
                   << " plane[2]=" << _frac_wholeimg_cosmic_cosmictag[2]
                   << " total=" << _frac_wholeimg_cosmic_cosmictag[3]
                   << std::endl;
      LARCV_INFO() << "fraction of TRUE COSMIC reco'd: "
                   << " plane[0]=" << _frac_wholeimg_cosmic_cosmictag[0]+_frac_wholeimg_cosmic_notcosmictag[0]
                   << " plane[1]=" << _frac_wholeimg_cosmic_cosmictag[1]+_frac_wholeimg_cosmic_notcosmictag[1]
                   << " plane[2]=" << _frac_wholeimg_cosmic_cosmictag[2]+_frac_wholeimg_cosmic_notcosmictag[2]
                   << " total=" << _frac_wholeimg_cosmic_cosmictag[3]+_frac_wholeimg_cosmic_notcosmictag[3]
                   << std::endl;

      LARCV_INFO() << "fraction of TRUE NU tagged as [cosmic]: "
                   << " plane[0]=" << _frac_wholeimg_nu_cosmictag[0]
                   << " plane[1]=" << _frac_wholeimg_nu_cosmictag[1]
                   << " plane[2]=" << _frac_wholeimg_nu_cosmictag[2]
                   << " total=" << _frac_wholeimg_nu_cosmictag[3]
                   << std::endl;
      LARCV_INFO() << "fraction of TRUE NU reco'd: "
                   << " plane[0]=" << _frac_wholeimg_nu_cosmictag[0]+_frac_wholeimg_nu_notcosmictag[0]
                   << " plane[1]=" << _frac_wholeimg_nu_cosmictag[1]+_frac_wholeimg_nu_notcosmictag[1]
                   << " plane[2]=" << _frac_wholeimg_nu_cosmictag[2]+_frac_wholeimg_nu_notcosmictag[2]
                   << " total=" << _frac_wholeimg_nu_cosmictag[3]+_frac_wholeimg_nu_notcosmictag[3]
                   << std::endl;

      LARCV_INFO() << "fraction of TRUE VTX tagged as [cosmic]: "
                   << " plane[0]=" << _frac_wholeimg_vtx_cosmictag[0]
                   << " plane[1]=" << _frac_wholeimg_vtx_cosmictag[1]
                   << " plane[2]=" << _frac_wholeimg_vtx_cosmictag[2]
                   << " total=" << _frac_wholeimg_vtx_cosmictag[3]
                   << std::endl;
      LARCV_INFO() << "fraction of TRUE VTX reco'd: "
                   << " plane[0]=" << _frac_wholeimg_vtx_cosmictag[0]+_frac_wholeimg_vtx_notcosmictag[0]
                   << " plane[1]=" << _frac_wholeimg_vtx_cosmictag[1]+_frac_wholeimg_vtx_notcosmictag[1]
                   << " plane[2]=" << _frac_wholeimg_vtx_cosmictag[2]+_frac_wholeimg_vtx_notcosmictag[2]
                   << " total=" << _frac_wholeimg_vtx_cosmictag[3]+_frac_wholeimg_vtx_notcosmictag[3]
                   << std::endl;
      
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
