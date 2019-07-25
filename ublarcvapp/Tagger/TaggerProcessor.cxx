#include "TaggerProcessor.h"

#include "TaggerCROITypes.h"

#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"

#include "ublarcvapp/UBImageMod/EmptyChannelAlgo.h"
#include "BMTCV.h"
#include "BoundaryMuonTaggerAlgo.h"
#include "FlashMuonTaggerAlgo.h"
#include "CACAEndPtFilter.h"

#include "TH2D.h"
#include "TStyle.h"
#include "TLine.h"
#include "larcv/core/ROOTUtil/ROOTUtils.h"

namespace ublarcvapp {
  namespace tagger {
    
    TaggerProcessor::TaggerProcessor( std::string cfg_path )
      : larcv::ProcessBase( "TaggerProcessor" ),
      _cfg_path(cfg_path),
      _larlite_io(nullptr),
      _RunPMTPrecuts(false),
      _ApplyPMTPrecuts(false),
      _FilterBoundaryPoints(true)
    {

      m_time_tracker.resize(kNumStages,0);
      
      // setup larlite file manager
      _larlite_io = new ublarcvapp::LArliteManager( larlite::storage_manager::kREAD, "tagger_larlite_input" );
      
      if (_cfg_path!="")
        configure(_cfg_path);


    }

    TaggerProcessor::~TaggerProcessor() {
      if ( _larlite_io ) {
        delete _larlite_io;
      }
    }
    
    void TaggerProcessor::configure(const larcv::PSet& pset) {
      set_verbosity( (larcv::msg::Level_t) pset.get<int>("Verbosity") );
      _RunPMTPrecuts   = pset.get<bool>("RunPreCuts");   // run the precuts
      _ApplyPMTPrecuts = pset.get<bool>("ApplyPreCuts"); // if True: if precuts fail, we do not run rest of algos (e.g. produce CROIs, thrumu mask)
    }

    /**
     * configure class using a config file
     *
     * affect is to modify the class member, m_config,
     * which holds the configuration parameters.
     *
     * @param[in] cfg_path path to configuration file
     *
     */
    void TaggerProcessor::configure(const std::string cfg_path) {
      m_config = TaggerCROIAlgoConfig::makeConfigFromFile( cfg_path );
      configure( m_config.main_cfg );
    }

    void TaggerProcessor::initialize() {
      _larlite_io->open();
    }
    
    void TaggerProcessor::add_larlite_input( std::string llinput ) {
      if ( !_larlite_io ) {
        LARCV_CRITICAL() << "Adding larlite file before creating manager. Must run configure first." << std::endl;
        throw std::runtime_error( "adding larlite file before creating manager" );
      }
      _larlite_io->add_in_filename( llinput );
    }
    
    bool TaggerProcessor::process(larcv::IOManager& io) {

      // we have to get a piece of data first, before the run, subrun, event
      // in the iomanager gets updated. we need it to sync with the larlite file

      auto event_imgs    = (larcv::EventImage2D*)io.get_data( larcv::kProductImage2D, m_config.larcv_image_producer );
      
      // sync larlite
      if ( !_larlite_io ) {
        LARCV_CRITICAL() << "Processing entry before creating manager. Must run configure first." << std::endl;
        throw std::runtime_error( "adding larlite file before creating manager" );
      }
      //_larlite_io->set_verbosity( larcv::msg::kDEBUG );
      _larlite_io->set_verbosity( larcv::msg::kINFO );
      _larlite_io->syncEntry( io );

      // prepare inputs
      InputPayload input;
      loadInput( io, *_larlite_io, input );

      // precuts
      runPMTPrecuts( _larlite_io );

      // ThruMu
      ThruMuPayload thrumu;
      
      // boudary end point finder
      runBoundaryTagger( input, thrumu );

      // visualize
      saveBoundaryEndpointImage( input, thrumu );
      
    }
    
    void TaggerProcessor::finalize() {
      if (_larlite_io)
        _larlite_io->close();
    }

    /**
     * load data into InputPayload class
     *
     */
    void TaggerProcessor::loadInput( larcv::IOManager& io_larcv,
                                     ublarcvapp::LArliteManager& io_larlite,
                                     InputPayload& input ) {

      // timing variables
      std::clock_t timer;
      
      // get images (from larcv)
      larcv::EventImage2D* event_imgs = NULL;
      try {
        event_imgs    = (larcv::EventImage2D*)io_larcv.get_data( larcv::kProductImage2D, m_config.larcv_image_producer );
        if ( event_imgs->Image2DArray().size()==0) {
          if ( !m_config.skip_empty_events )
            throw std::runtime_error("Number of images=0. LArbys.");
          else {
            std::cout << "Skipping Empty Events." << std::endl;
            throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
          }
        }
        else if ( event_imgs->Image2DArray().size()!=3 ) {
          throw std::runtime_error("Number of Images!=3. Weird.");
        }
      }
      catch (const std::exception &exc) {
        std::cerr << "Error retrieving TPC Images: " << exc.what() << std::endl;
        throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
      }

      input.img_v.clear();
      for ( auto const& img : event_imgs->Image2DArray() )
        input.img_v.push_back( img );

      // Get run, subrun, event, entry
      // quirk of larcv iomanager that this is not filled until first data product loaded
      input.run    = io_larcv.event_id().run();
      input.subrun = io_larcv.event_id().subrun();
      input.event  = io_larcv.event_id().event();      
      input.entry  = io_larcv.current_entry();
      

      // ------------------------------------------------------------------------------------------//
      // MODIFY IMAGES

      // algos
      ublarcvapp::EmptyChannelAlgo emptyalgo;
            
      try {
        // MCC8 modification 
        if ( m_config.DeJebWires ) {
          // raise the adc values for pixels in a certain region
          for ( auto &img : input.img_v ) {
            const larcv::ImageMeta& meta = img.meta();
            for ( int col=0; col<(int)meta.cols(); col++) {
              for (int row=0; row<(int)meta.rows(); row++) {
                float val = img.pixel( row, col );
                if (meta.plane()==0 ) {
                  if ( (col>=2016 && col<=2111) || (col>=2176 && 2212) ) {
                    val *= m_config.jebwiresfactor;
                  }
                }
                img.set_pixel(row,col,val);
              }
            }
          }
        }
      }
      catch (const std::exception& e ) {
        std::cerr << "Error Dejebbing wires: " << e.what() << std::endl;
        throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
      }

      // larlitecv::UnipolarHackAlgo unihackalgo;
      // try {
      //   if ( m_config.apply_unipolar_hack ) {
      //     // find signs that a unipolar pulse on the collection plane was deconvolved poorly
      //     std::vector<int> applyhack(3,0);
      //     applyhack[1] = 1;
      //     std::vector<float> hackthresh(3,-10.0);
      //     input.img_v = unihackalgo.hackUnipolarArtifact( input.img_v, applyhack, hackthresh );
      //   }
      //   if ( input.img_v.size()!=3 )
      //     throw std::runtime_error("Number of Images incorrect.");
      // }
      // catch ( const std::exception& e ) {
      //   std::cerr << "Error making unipolar artifact: " << e.what() << std::endl;
      //   throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
      // }

      // ------------------------------------------------------------------------------------------//
      // LABEL EMPTY CHANNELS
      // larlitecv::EmptyChannelAlgo emptyalgo;
      // std::vector< larcv::Image2D > emptyimgs;
      // try {
      //   for ( auto const &img : event_imgs->Image2DArray() ) {
      //     int p = img.meta().plane();
      //     larcv::Image2D emptyimg = emptyalgo.labelEmptyChannels( m_config.emptych_thresh.at(p), img );
      //     emptyimgs.emplace_back( emptyimg );
      //   }
      // }
      // catch ( const std::exception& e ) {
      //   std::cerr << "Error making empty channels: " << e.what() << std::endl;
      //   throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");      
      // }
      
      // ------------------------------------------------------------------------------------------//
      // LABEL BAD CHANNELS

      input.badch_v.clear();
      try {
        if ( m_config.chstatus_datatype=="LARCV" ) {
          larcv::EventChStatus* ev_status = (larcv::EventChStatus*)io_larcv.get_data( larcv::kProductChStatus, m_config.larcv_chstatus_producer );
          input.badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
          ev_status->clear(); // clear, we copied the info
        }
        else if ( m_config.chstatus_datatype=="LARLITE" ) {
          larlite::event_chstatus* ev_status = (larlite::event_chstatus*)io_larlite.get_data( larlite::data::kChStatus, "chstatus" );
          input.badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
          ev_status->clear(); // clear, we copied the info
        }        
        else if ( m_config.chstatus_datatype=="NONE" ) {
          for ( auto const& img : input.img_v ) {
            larcv::Image2D badch( img.meta() );
            badch.paint(0.0);
            input.badch_v.emplace_back( std::move(badch) );
          }
        }
        else {
          throw std::runtime_error("ERROR: ChStatusDataType must be either LARCV or LARLITE");
        }
      
        if ( input.badch_v.size()!=3 ) {
          throw std::runtime_error("Number of Bad Channels not correct.");
        }
      }
      catch ( const std::exception& e ) {
        std::cerr << "Error making badch image: " << e.what() << std::endl;
        throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
      }

      // ------------------------------------------------------------------------------------------//
      // MAKE GAP CHANNEL IMAGE

      int maxgap = 200; // this seems bad to have this parameter here and not in a config file
      std::vector< larcv::Image2D> gapchimgs_v;
      try {
        gapchimgs_v = emptyalgo.findMissingBadChs( input.img_v, input.badch_v, 5, maxgap );
        // combine with badchs
        for ( size_t p=0; p<gapchimgs_v.size(); p++ ) {
          larcv::Image2D& gapchimg = gapchimgs_v[p];
          gapchimg += input.badch_v[p];
        }

        // Set BadCh in input data
        for ( size_t p=0; p<gapchimgs_v.size(); p++ ) {
          input.gapch_v.emplace_back( std::move(gapchimgs_v.at(p)) );
        }
        if ( input.gapch_v.size()!=3 )
          throw std::runtime_error("Number of Gap Channels not correct.");
      }
      catch ( const std::exception& e ) {
        std::cerr << "Error making gapch images: " << e.what() << std::endl;
        throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");      
      }

      // -------------------------------------------------------------------------------------------//
      // COLLECT FLASH INFORMATION
    
      try {
        for ( auto &flashproducer : m_config.opflash_producers ) {
          larlite::event_opflash* opdata = (larlite::event_opflash*)io_larlite.get_data(larlite::data::kOpFlash, flashproducer );
          LARCV_INFO() << "search for flash hits from " << flashproducer << ": " << opdata->size() << " flashes" << std::endl;
          input.opflashes_v.push_back( opdata );
          for ( auto const& opflash : *opdata ) {
            LARCV_DEBUG() << "  flashtime: usec=" << opflash.Time()
                          << " tick=" << opflash.Time()/0.015625
                          << " pe=" << opflash.TotalPE()
                          << std::endl;
          }
        }
      }
      catch ( const std::exception& e ) {
        std::cerr << "Error retrieving flash information: " << e.what() << std::endl;
        throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
      }
      
      // -------------------------------------------------------------------------------------------//
      // LOAD MC TRACK INFORMATION: USED FOR PERFORMANCE METRICS
      if ( m_config.load_mctrack ) {
        try {
          input.p_ev_mctrack = (larlite::event_mctrack*)io_larlite.get_data(larlite::data::kMCTrack, m_config.mctrack_producer);
        }
        catch (const std::exception& e ) {
          std::cerr << "Error retrieving MC track information upon request: " << e.what() << std::endl;
          throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
        }
      }

      // -------------------------------------------------------------------------------------------//
      // LOAD EVENT TRIGGER
      try {
        input.p_ev_trigger = (larlite::trigger*)io_larlite.get_data(larlite::data::kTrigger, m_config.trigger_producer);
      }
      catch (const std::exception& e ) {
        std::cerr << "Error retrieving MC track information upon request: " << e.what() << std::endl;
        input.p_ev_trigger = NULL;
        throw std::runtime_error("[TaggerCROIAlgo::loadInput] ERROR");
      }

      // ----------------------------------
      // RECORD TIMING FOR THIS PROCESS
      m_time_tracker[kInputTotal] += double(std::clock()-timer)/(double)CLOCKS_PER_SEC;
      
    }
    
    /**
     *
     * run precuts
     *
     * the pmt precut algorithm requires larlite::ophit objects
     * 
     * @param[in] llio pointer to larlite storage manager with event data
     */
    void TaggerProcessor::runPMTPrecuts( larlite::storage_manager* llio) {

      LARCV_INFO() << "== Run PMT Precuts  ===============================" << std::endl;

      // timing variables
      std::clock_t timer;
      
      // conversion into larlite fcllite::PSet
      fcllite::PSet tmp( "tmp",m_config.precut_cfg.data_string() );
      fcllite::PSet precutcfg( tmp.get<fcllite::PSet>("LEEPreCut") );
      larlite::LEEPreCut  precutalgo;
      precutalgo.configure( precutcfg ); 
      precutalgo.initialize();

      if ( _RunPMTPrecuts ) {
        precutalgo.analyze( llio );
        // save data in larlite user data structure
        //larlite::event_user* ev_precutresults = (larlite::event_user*)dataco_out.get_larlite_data( larlite::data::kUserInfo, "precutresults" );
        //larlite::user_info precut_results;
        LARCV_INFO() << "==== PRECUT RESULTS ===============" << std::endl;
        LARCV_INFO() << " PASS: "     << precutalgo.passes() << std::endl;
        LARCV_INFO() << " vetope:   " << precutalgo.vetoPE() << std::endl;
        LARCV_INFO() << " beampe:   " << precutalgo.beamPE() << std::endl;
        LARCV_INFO() << " maxfrac:  " << precutalgo.maxFrac() << std::endl;
        LARCV_INFO() << " beamTick: " << precutalgo.beamFirstTick() << std::endl;
        LARCV_INFO() << " vetoTick: " << precutalgo.vetoFirstTick() << std::endl;
        LARCV_INFO() << "===================================" << std::endl;
	
        // precut_results.store( "pass",    precutalgo.passes() );
        // precut_results.store( "vetoPE",  precutalgo.vetoPE() );
        // precut_results.store( "beamPE",  precutalgo.beamPE() );
        // precut_results.store( "maxFrac", precutalgo.maxFrac() );
        // precut_results.store( "beamFirstTick", precutalgo.beamFirstTick() );
        // precut_results.store( "vetoFirstTick", precutalgo.vetoFirstTick() );      
        // ev_precutresults->emplace_back( std::move(precut_results) );
      }
      else {
        // create dummy values
        // larlite::event_user* ev_precutresults = (larlite::event_user*)dataco_out.get_larlite_data( larlite::data::kUserInfo, "precutresults" );
        // larlite::user_info precut_results;      
        // precut_results.store( "pass",     1   );
        // precut_results.store( "vetoPE",  -1.0 );
        // precut_results.store( "beamPE",  -1.0 );
        // precut_results.store( "maxFrac", -1.0 );
        // precut_results.store( "beamFirstTick", -1 );
        // precut_results.store( "vetoFirstTick", -1 );
        // ev_precutresults->emplace_back( std::move(precut_results) );
      }

      m_time_tracker[kPMTPrecuts] += double(std::clock()-timer)/(double)CLOCKS_PER_SEC;
      
      return;
    }// end of runpmtprecuts

    // void TaggerProcessor::runBoundaryPointFinder() {

    // }

    void TaggerProcessor::runBoundaryTagger( const InputPayload& input, ThruMuPayload& output ) {

      LARCV_INFO() << "== Run Boundary Tagger ===============================" << std::endl;
      

      // configure different stages of the Thrumu Tagger
      std::clock_t timer_total;
      std::clock_t timer;

      BMTCV m_bmtcv_algo;
      
      // (0) make contours
      timer = std::clock();
      m_bmtcv_algo.analyzeImages( input.img_v, input.badch_v, 10.0, 2 );
      m_time_tracker[kThruMuContour] += double(std::clock()-timer)/(double)CLOCKS_PER_SEC;

      // (1) side tagger
      timer = std::clock();
      BoundaryMuonTaggerAlgo sidetagger;
      sidetagger.configure( m_config.sidetagger_cfg );
      //sidetagger.printConfiguration();

      // (2) flash tagger
      FlashMuonTaggerAlgo anode_flash_tagger(   FlashMuonTaggerAlgo::kSearchAnode );
      FlashMuonTaggerAlgo cathode_flash_tagger( FlashMuonTaggerAlgo::kSearchCathode );
      FlashMuonTaggerAlgo imgends_flash_tagger( FlashMuonTaggerAlgo::kSearchOutOfImage );
    
      anode_flash_tagger.configure(   m_config.flashtagger_cfg );
      cathode_flash_tagger.configure( m_config.flashtagger_cfg );
      imgends_flash_tagger.configure( m_config.flashtagger_cfg );
      
      // loading time
      m_time_tracker[kThruMuConfig] = ( std::clock()-timer )/(double)CLOCKS_PER_SEC;

      // run side tagger
      timer = std::clock();
      sidetagger.searchforboundarypixels3D( input.img_v, input.badch_v, output.side_spacepoint_v, output.boundarypixel_image_v, output.realspacehit_image_v );
      int nsides[4] = {0};
      for ( auto const& sp : output.side_spacepoint_v ) {
        nsides[ sp.at(0).type ]++;
      }

      LARCV_INFO() << " Side Tagger End Points: " << output.side_spacepoint_v.size() << std::endl;
      LARCV_INFO() << "   Top: "        << nsides[0] << std::endl;
      LARCV_INFO() << "   Bottom: "     << nsides[1] << std::endl;
      LARCV_INFO() << "   Upstream: "   << nsides[2] << std::endl;
      LARCV_INFO() << "   Downstream: " << nsides[3] << std::endl;

      m_time_tracker[kThruMuBMT] += (std::clock()-timer)/(double)CLOCKS_PER_SEC;

      // run flash tagger
      timer = std::clock();
      anode_flash_tagger.flashMatchTrackEnds(   input.opflashes_v, input.img_v, input.badch_v, output.anode_spacepoint_v );
      cathode_flash_tagger.flashMatchTrackEnds( input.opflashes_v, input.img_v, input.badch_v, output.cathode_spacepoint_v  );
      imgends_flash_tagger.findImageTrackEnds(  input.img_v, input.badch_v, output.imgends_spacepoint_v  );

      int totalflashes = (int)output.anode_spacepoint_v.size() + (int)output.cathode_spacepoint_v.size() + (int)output.imgends_spacepoint_v.size();
      LARCV_INFO() << " Flash Tagger End Points: " << totalflashes << std::endl;
      LARCV_INFO() << "  Anode: "      << output.anode_spacepoint_v.size() << std::endl;
      LARCV_INFO() << "  Cathode: "    << output.cathode_spacepoint_v.size() << std::endl;
      LARCV_INFO() << "  Image Ends: " << output.imgends_spacepoint_v.size() << std::endl;

      LARCV_INFO() << "Anode spacepoint flash indices: " << std::endl;
      for (int i=0; i<(int)output.anode_spacepoint_v.size(); i++) {
        LARCV_INFO() << "    [" << i << "] flashidx="
                     << "(" << output.anode_spacepoint_v.at(i).getFlashIndex().ivec << ","
                     << output.anode_spacepoint_v.at(i).getFlashIndex().idx << ")"	
                     << std::endl;
      }
      LARCV_INFO() << "Cathode spacepoint flash indices: " << std::endl;
      for (int i=0; i<(int)output.cathode_spacepoint_v.size(); i++) {
        LARCV_INFO() << "    [" << i << "] flashidx="
                     << "(" << output.cathode_spacepoint_v.at(i).getFlashIndex().ivec << ","
                     << output.cathode_spacepoint_v.at(i).getFlashIndex().idx << ")"	
                     << std::endl;
      }
      
      m_time_tracker[kThruMuFlash] += (std::clock()-timer)/(double)CLOCKS_PER_SEC;

      // run end point filters

      //larlitecv::RadialEndpointFilter radialfilter;  // remove end points that cannot form a 3d segment nearby [deprecated]
      //larlitecv::PushBoundarySpacePoint endptpusher; // remove endpt [deprecated]
      //larlitecv::EndPointFilter endptfilter; // removes duplicates [deprecated]

      if ( _FilterBoundaryPoints ) {
        timer = std::clock();
        CACAEndPtFilter cacaalgo;
        cacaalgo.setVerbosity(0);

        // we collect pointers to all the end points (make a copy for now)
        std::vector< BoundarySpacePoint > all_endpoints;
        
        // gather endpoints from space points
        for (int isp=0; isp<(int)output.side_spacepoint_v.size(); isp++) {
          const BoundarySpacePoint* pts = &(output.side_spacepoint_v.at( isp ));
          all_endpoints.push_back( *pts );
        }
        for (int isp=0; isp<(int)output.anode_spacepoint_v.size(); isp++) {
          const BoundarySpacePoint* pts = &(output.anode_spacepoint_v.at(isp));
          all_endpoints.push_back( *pts );
        }
        for (int isp=0; isp<(int)output.cathode_spacepoint_v.size(); isp++) {
          const BoundarySpacePoint* pts = &(output.cathode_spacepoint_v.at(isp));
          all_endpoints.push_back( *pts );
        }
        for (int isp=0; isp<(int)output.imgends_spacepoint_v.size(); isp++) {
          const BoundarySpacePoint* pts = &(output.imgends_spacepoint_v.at(isp));
          all_endpoints.push_back( *pts );
        }
        if ( m_config.verbosity>0 )
          std::cout << "number of endpoints pre-filters: " << all_endpoints.size() << std::endl;

        std::vector< const std::vector<BoundarySpacePoint>* > sp_v;
        sp_v.push_back( &all_endpoints );
        std::vector< std::vector<int> > caca_results;    
        cacaalgo.evaluateEndPoints( sp_v, input.opflashes_v, input.img_v, input.badch_v, m_bmtcv_algo.m_plane_atomicmeta_v, 150.0, caca_results );
      
        // prepare the boundary points that pass
        std::vector<BoundarySpacePoint> cacapassing_moved_v = cacaalgo.regenerateFitleredBoundaryPoints( input.img_v );
      
        // clean up
        all_endpoints.clear();
        sp_v.clear();
        m_time_tracker[kThruMuFilter] += double(std::clock()-timer)/(double)CLOCKS_PER_SEC;
    
        
        // collect output
        // remove the filtered end points
        for ( size_t idx=0; idx<cacapassing_moved_v.size(); idx++ ) {
          BoundarySpacePoint& sp = cacapassing_moved_v[idx];
          
          if (sp.type()<=kDownstream ) {
            output.side_filtered_v.emplace_back( std::move(sp) );
          }
          else if (sp.type()==kAnode) {
            output.anode_filtered_v.emplace_back( std::move(sp) );
          }
          else if (sp.type()==kCathode) {
            output.cathode_filtered_v.emplace_back( std::move(sp) );
          }
          else if (sp.type()==kImageEnd) {
            output.imgends_filtered_v.emplace_back( std::move(sp) );
          }
          else {
            std::stringstream ss;
            ss << __FILE__ << ":" << __LINE__ << " unrecognized boundary type" << std::endl;
            throw std::runtime_error(ss.str());
          }
        }
      }
      else {
        // just copy over, if not filtering
        output.side_filtered_v    = output.side_spacepoint_v;
        output.anode_filtered_v   = output.anode_spacepoint_v;
        output.cathode_filtered_v = output.cathode_spacepoint_v;
        output.imgends_filtered_v = output.imgends_spacepoint_v;
      }

      LARCV_INFO() << " Filtered Side Tagger End Points -----" << std::endl;      
      LARCV_INFO() << "   Side: "        << output.side_filtered_v.size() << std::endl;
      LARCV_INFO() << "   Anode: "       << output.anode_filtered_v.size() << std::endl;
      LARCV_INFO() << "   Cathode: "     << output.cathode_filtered_v.size() << std::endl;
      LARCV_INFO() << "   ImageEnds: "   << output.imgends_filtered_v.size() << std::endl;
      
      
      LARCV_INFO() << "== End of Boundary Tagger ===============================" << std::endl;
      
      return;
    }//end of runBoundaryTagger

    /*
    void TaggerProcessor::runThrumu( InputPayload& input, ThruMuPayload& thrumu ) {

      LARCV_INFO() << "== Run ThruMu Tracker ===============================" << std::endl;
      
      // configure different stages of the Thrumu Tagger
      std::clock_t timer = std::clock();

      // thrumu tracker
      ThruMuTracker thrumu_tracker( m_config.thrumu_tracker_cfg );

      m_time_tracker[kThruMuConfig] = ( std::clock()-timer )/(double)CLOCKS_PER_SEC;

      // RUN THE THRUMU ALGOS

      // run side tagger
      timer = std::clock();

      // we collect pointers to all the end points
      std::vector< const BoundarySpacePoint* > all_endpoints;

      // Collect the indices for the flash that determines each of the endpoints and
      // the index for the producer that makes each of the endpoints as well.
      // This is a precaution to ensure that we are collecting this information properly.
      std::vector< int > all_endpoints_flash_idx_v;
      all_endpoints_flash_idx_v.clear();

      std::vector< int > all_endpoints_boundary_type_idx_v;
      all_endpoints_boundary_type_idx_v.clear();

      // When filling the 'all_endpoints_flash_idx_v' and 'all_endpoints_boundary_type_idx_v',
      // the correct flash indices have already been matched to the correct endpoint.
      // The important thing is that these two vectors are filled at the same point as their corresponding endpoint in 'all_endpoints'.

      // gather endpoints from space points
      for (int isp=0; isp<(int)output.side_filtered_v.size(); isp++) {
        const larlitecv::BoundarySpacePoint* pts = &(output.side_filtered_v.at( isp ));
        all_endpoints.push_back( pts );
      }
      for (int isp=0; isp<(int)output.anode_filtered_v.size(); isp++) {
        const larlitecv::BoundarySpacePoint* pts = &(output.anode_filtered_v.at(isp));
        all_endpoints.push_back( pts );
      }
      for (int isp=0; isp<(int)output.cathode_filtered_v.size(); isp++) {
        const larlitecv::BoundarySpacePoint* pts = &(output.cathode_filtered_v.at(isp));
        all_endpoints.push_back( pts );
      }
      for (int isp=0; isp<(int)output.imgends_filtered_v.size(); isp++) {
        const larlitecv::BoundarySpacePoint* pts = &(output.imgends_filtered_v.at(isp));
        all_endpoints.push_back( pts );
      }
      LARCV_INFO() << "number of endpoints to search for thrumu: " << all_endpoints.size() << std::endl;

      // make track clusters
      std::vector<int> used_endpoints( all_endpoints.size(), 0 );
      if ( m_config.run_thrumu_tracker ) {
        timer = std::clock();
        output.trackcluster3d_v.clear();
        output.tagged_v.clear();
        
        thrumu_tracker.makeTrackClusters3D( m_config.croi_selection_cfg.genflashmatch_cfg, input.img_v, input.gapch_v,
                                            all_endpoints, output.trackcluster3d_v, 
                                            output.tagged_v, used_endpoints, input.opflashes_v ); 
        
        m_time_tracker[kThruMuTracker]  +=  (std::clock()-timer)/(double)CLOCKS_PER_SEC;
        if ( m_config.verbosity>0 )
          std::cout << "thrumu tracker search " << all_endpoints.size() << " end points in " << m_time_tracker[kThruMuTracker] << " sec" << std::endl;
      }
      else {
        if ( m_config.verbosity>0 )
          std::cout << "config tells us to skip thrumu track." << std::endl;
      }
    
      // collect unused endpoints
      output.used_spacepoint_v.clear();
      output.unused_spacepoint_v.clear();
      for ( size_t isp=0; isp<all_endpoints.size(); isp++ ) {
        if ( used_endpoints.at(isp)==1 )
          output.used_spacepoint_v.push_back( *(all_endpoints.at(isp)) );
        else {
          const BoundarySpacePoint& sp = *(all_endpoints.at(isp));
          if ( m_config.verbosity>1 )
            std::cout << "unused spacepoint for StopMu: (" << sp.pos()[0] << "," << sp.pos()[1] << "," << sp.pos()[2] << ")" << std::endl;
          output.unused_spacepoint_v.push_back( *(all_endpoints.at(isp)) );
        }
      }
      
      // copy track and pixels into separate containers.
      output.track_v.clear();
      output.pixelcluster_v.clear();
      for ( auto const& bmtrack : output.trackcluster3d_v ) {
        output.track_v.push_back( bmtrack.makeTrack() );
        std::vector< larcv::Pixel2DCluster > cluster_v;
        for ( auto const& track2d : bmtrack.plane_pixels )
          cluster_v.push_back( track2d );
        output.pixelcluster_v.emplace_back( std::move(cluster_v) );
      }
      
      LARCV_INOF() << "== End of ThruMu Tracker ===============================" << std::endl;

    }
    */
    
    void TaggerProcessor::saveBoundaryEndpointImage( InputPayload& input, ThruMuPayload& thrumu ) {
      gStyle->SetOptStat(0);
      for ( size_t p=0; p<3; p++ ) {

        TCanvas c("c","c",1200,600);

        char histname[20];
        sprintf( histname, "hadc_plane%d", (int)p );
        TH2D hadc = larcv::as_th2d( input.img_v.at(p), histname );
        hadc.Draw("colz");
        hadc.GetZaxis()->SetRangeUser(10,100);

        std::vector<TGraph*> g_v;
        std::vector< const std::vector<BoundarySpacePoint>* > sp_vv
          = { &thrumu.side_filtered_v,
              &thrumu.anode_filtered_v,
              &thrumu.cathode_filtered_v,
              &thrumu.imgends_filtered_v };
        int marker_v[4] = { 24, 25, 26, 27 };

        size_t pt_t = 0;
        for ( auto const& p_sp_v : sp_vv ) {
          
          TGraph* g = new TGraph( p_sp_v->size() );
          g->SetMarkerStyle( marker_v[pt_t] );
          for ( size_t ipt=0; ipt<p_sp_v->size(); ipt++ ) {
            auto const& sp = p_sp_v->at(ipt);
            int tick = sp.tick( input.img_v.at(p).meta() );
            std::vector<int> wires = sp.wires( input.img_v.at(p).meta() );
            g->SetPoint( ipt, wires[p], tick );
          }
          g_v.push_back( g );
          pt_t++;
        }

        // draw lines for flashes
        std::vector< TLine* > line_cathode_v;
        std::vector< TLine* > line_anode_v;        
        for ( auto const& p_opdata : input.opflashes_v ) {

          for ( auto const& opflash : *p_opdata ) {
            int tick = 3200 + opflash.Time()*2;
            TLine* lanode  = new TLine( 0, tick, 3455, tick );
            line_anode_v.push_back( lanode );
            TLine* lcathode  = new TLine( 0, tick+255.0/0.114*2, 3455, tick+255.0/0.114*2 );
            line_cathode_v.push_back( lcathode );
          }

        }

        for ( auto& p_line : line_anode_v ) {
          p_line->SetLineColor(kMagenta);
          p_line->Draw("L");
        }
        for ( auto& p_line : line_cathode_v ) {
          p_line->SetLineColor(kCyan);
          p_line->Draw("L");
        }
        for ( auto& p_g : g_v ) p_g->Draw("P");   
        
        c.Draw();
        c.Update();
        
        char canvname[100];
        sprintf( canvname, "bmt_plane%d.pdf", (int)p );
        c.SaveAs( canvname );
        c.Close();

        for ( auto& p_g : g_v )
          delete p_g;
      }
      
    }//end of saveboundary image

    
    
  }//end of tagger namespace
}// end of ublarcvapp namespace
