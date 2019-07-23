#include "TaggerProcessor.h"

#include "TaggerCROITypes.h"

#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"


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

      // algos
      ///larlitecv::EmptyChannelAlgo emptyalgo;
      
      // prepare inputs
      InputPayload input;
      input.clear();

      loadInput( io, *_larlite_io, input );
      
    }
    
    void TaggerProcessor::finalize() {
      
    }

    /**
     * load data into InputPayload class
     *
     */
    void TaggerProcessor::loadInput( larcv::IOManager& io_larcv,
                                     ublarcvapp::LArliteManager& io_larlite,
                                     InputPayload& input ) {

      // Get run, subrun, event, entry
      input.run    = io_larcv.event_id().run();
      input.subrun = io_larcv.event_id().subrun();
      input.event  = io_larcv.event_id().event();      
      input.entry  = io_larcv.current_entry();

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

      // ------------------------------------------------------------------------------------------//
      // MODIFY IMAGES
      
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

      // larlitecv::EmptyChannelAlgo emptyalgo;
      
      input.badch_v.clear();
      try {
        if ( m_config.chstatus_datatype=="LARCV" ) {
          larcv::EventChStatus* ev_status = (larcv::EventChStatus*)io_larcv.get_data( larcv::kProductChStatus, m_config.larcv_chstatus_producer );
          //input.badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
          ev_status->clear(); // clear, we copied the info
        }
        // else if ( m_config.chstatus_datatype=="LARLITE" ) {
        //   larlite::event_chstatus* ev_status = (larlite::event_chstatus*)dataco_input.get_larlite_data( larlite::data::kChStatus, "chstatus" );
        //   input.badch_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
        //   ev_status->clear(); // clear, we copied the info
        // }        
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
        //gapchimgs_v = emptyalgo.findMissingBadChs( input.img_v, input.badch_v, 5, maxgap );
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
          std::cout << "search for flash hits from " << flashproducer << ": " << opdata->size() << " flashes" << std::endl;
          input.opflashes_v.push_back( opdata );
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
      
    }
    
    
  }
}
