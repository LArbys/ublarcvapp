#include "ThruMuTracker.h"

#include "RadialEndpointFilter.h"
#include "AStarNodes2BMTrackCluster3D.h"
#include "ThruMuFoxExtender.h"
#include "PushBoundarySpacePoint.h"

// larcv 
#include "ublarcvapp/UBWireTool/UBWireTool.h"
#include "ublarcvapp/Reco3D/AStar3DAlgo.h"


namespace ublarcvapp {
namespace tagger {

  ThruMuTracker::ThruMuTracker( const ThruMuTrackerConfig& config )
    : m_config(config)
  {}

  // Add additional arguments: 
  //                           
  void ThruMuTracker::makeTrackClusters3D( const std::vector<larcv::Image2D>& img_v,  
					   const std::vector<larcv::Image2D>& badchimg_v,
                                           const std::vector< const BoundarySpacePoint* >& spacepts,
					   std::vector< BMTrackCluster3D >& trackclusters,
					   std::vector< larcv::Image2D >& tagged_v,
                                           std::vector<int>& used_endpoints_indices,
                                           const std::vector< larlite::event_opflash* >& opflashsets ) {

    // This method takes in the list of boundaryspacepoints and pairs them up in order to try to find through-going muons
    // input:
    //   flash_match_config': configuration for 'GeneralFlashMatchAlgo' instance within the 'flashMatchTracks' function.    
    //   img_v: vector images, one for each plane. assumed to be in U,V,Y order
    //   badchimg_v: vector images, one for each plane. assumed to be in U,V,Y order    
    //   spacepts: all possible boundary crossing points
    // output:
    //   trackclusters: thrumu tracks. collection of space points and pixel2dclusters
    //   tagged_v: pixels tagged as thrumu.  the value indicates the pass when the muon was tagged.
    //   used_endpoints_indices: indicates which endpoints in spacepts was used. 0=not used. 1=used.

    // Declare a set, 'impossible_match_endpoints', meant to save the indices of endpoints that together constituted an impossible match.
    // This removes from consideration endpoints that were found to fail the flashmatching stage of reconstruction.                                                                                   
    std::set< std::vector< BoundarySpacePoint > > impossible_match_endpoint_v;

    // Declare a new vector, 'track_endpoint_indices', based off the indices in 'spacepts', for the endpoints that are used to make tracks.
    std::vector < std::vector< BoundarySpacePoint >  > track_endpoint_v;
    track_endpoint_v.clear();

    // Declare a vector that contains the 'BoundaryFlashIndex' information for the flash that the endpoints correspond to.
    std::vector< std::vector< BoundaryFlashIndex > > track_endpoint_flash_v;
    track_endpoint_flash_v.clear();

    // Declare a vector for the indices of the boundary of the detector that this flash corresponds to.
    std::vector< std::vector< BoundaryEnd_t >  > track_endpoint_boundary_type_idx_v;
    track_endpoint_boundary_type_idx_v.clear();

    // Declare a set of the 'BoundaryFlashIndex' information of the flashes that have already been well-matched to a track.
    std::set< BoundaryFlashIndex > already_matched_flash_v;

    // Declare a vector for the well-matched tracks 'well_matched_tracks_idx_v'.
    std::vector < int > well_matched_tracks_idx_v;
    well_matched_tracks_idx_v.clear();

    const int nendpts = (int)spacepts.size();
    used_endpoints_indices.resize( nendpts, 0 );

    // pair up containers
    int ntotsearched = 0;
    int npossible = 0;

    // compress image for astar. store in member class. we'll use these again in runAStarTagger
    m_img_compressed_v.clear();
    m_badch_compressed_v.clear();
    int downsampling_factor = m_config.downsampling_factor;
    for (size_t p=0; p<img_v.size(); p++) {
      larcv::Image2D img_compressed( img_v[p] );
      larcv::Image2D badch_compressed( badchimg_v[p] );
      img_compressed.compress( img_v[p].meta().rows()/downsampling_factor,
                                img_v[p].meta().cols()/downsampling_factor,
                                (larcv::Image2D::CompressionModes_t)m_config.compression_mode );
      badch_compressed.compress( img_v[p].meta().rows()/downsampling_factor,
                                img_v[p].meta().cols()/downsampling_factor,
                                (larcv::Image2D::CompressionModes_t)m_config.compression_mode );
      m_img_compressed_v.emplace_back( std::move(img_compressed) );
      m_badch_compressed_v.emplace_back( std::move(badch_compressed) );
    }

    // tagged image
    if ( tagged_v.size()==0 )  {
      for (size_t p=0; p<img_v.size(); p++) {
        larcv::Image2D tagged( img_v[p].meta() );
        tagged.paint(0);
        tagged_v.emplace_back( std::move(tagged) );
      }
    }


    // poor-man's profiling
    const clock_t begin_time = clock();

    std::vector<int> tracks_per_pass;
    for (int ipass=0; ipass<m_config.num_passes; ipass++) {
      
      // Clear 'well_matched_tracks_idx_v' at the start of the pass.
      well_matched_tracks_idx_v.clear();
      
      const ThruMuTrackerConfig::ThruMuPassConfig& passcfg = m_config.pass_configs.at(ipass);
      runPass( ipass, passcfg, spacepts, img_v, badchimg_v, tagged_v, used_endpoints_indices, trackclusters,
	       track_endpoint_flash_v, track_endpoint_boundary_type_idx_v, track_endpoint_v );

      bool anode_and_cathode_only = true;
      
      int tracks_in_pass = trackclusters.size();
      for ( int i = int( tracks_per_pass.size() - 1 ); i > -1; --i ) {
        tracks_in_pass -= tracks_per_pass.at( i );
      }
      tracks_per_pass.push_back( tracks_in_pass );

      // Depending if this is greater than the first pass, then set 'anode_and_cathode_only' to 'false'.
      if ( ipass > 0 ) 
	anode_and_cathode_only = false;

      // Use the flag from the config file to determine if you want to flashmatch the tracks that survive this pass of the thrumu tracker.
      // if (m_config.thrumu_flashmatch == true ) {

      //   bool single_pass = true;

      //   // Call a function that uses all of the flashmatching infrastructure developed in 'GeneralFlashMatchAlgo'.
      //   flashMatchTracks( flash_match_config, spacepts, opflashsets, trackclusters, impossible_match_endpoint_v, 
      //   		  already_matched_flash_v, well_matched_tracks_idx_v, tracks_in_pass, track_endpoint_flash_v, 
      //   		  track_endpoint_boundary_type_idx_v, track_endpoint_v, anode_and_cathode_only );	
      //   sortOutBadTracks( trackclusters, well_matched_tracks_idx_v, tracks_per_pass, tracks_per_pass.at( tracks_per_pass.size() - 1 ), single_pass ); 

      // }
      
    }

    // tag the pixels
    for (int itrack=0; itrack<(int)trackclusters.size(); itrack++ ) {
      BMTrackCluster3D& track3d = trackclusters[itrack];
      track3d.markImageWithTrack( img_v, badchimg_v, m_config.pixel_threshold, m_config.tag_neighborhood, tagged_v, 0.3, itrack+1 );
    }

    if ( m_config.verbosity>0 ) {
      for (int ipass=0; ipass<m_config.num_passes; ipass++) {
        std::cout << "Number of tracks found in ThruMu pass #" << ipass << ": " << tracks_per_pass.at(ipass) << std::endl;
      }
    }
    
  }

  // Update this function to pass the index of the flash that is matched to the endpoints of the track along with the indices themselves.
  // New arguments: 
  //   track_endpoint_flash_idx_v: This vector of ints corresponds to the index of the flash that is matched to the track endpoints at this point in the vector.
  //   track_endpoint_boundary_type_idx_v: This vector of ints contains the producer of the flash that corresponds to the endpoints of the track in the 'trackclusters' vector.
  void ThruMuTracker::runPass( const int passid, const ThruMuTrackerConfig::ThruMuPassConfig& passcfg,
                               const std::vector< const BoundarySpacePoint* >& spacepts,
			       const std::vector<larcv::Image2D>& img_v,
                               const std::vector<larcv::Image2D>& badchimg_v,
                               std::vector<larcv::Image2D>& tagged_v,
			       std::vector<int>& used_endpoints_indices,
                               std::vector<BMTrackCluster3D>& trackclusters,
			       std::vector< std::vector< BoundaryFlashIndex > >& track_endpoint_flash_v,
			       std::vector< std::vector< BoundaryEnd_t > >& track_endpoint_boundary_type_idx_v,
                               std::vector< std::vector< BoundarySpacePoint > >& track_endpoint_v) {

    // Make a pass to try and connect some end points.
    // A pass consists of running the linear3d tagger and the astar3d tagger.
    // Cuts on the quality of each determine which, if any, of the taggers are run and produce the track.
    // input:
    //   passid: id for the current pass
    //   passcfg: configuration for the pass
    //   spacepts: boundary end points
    //   img_v: original images
    //   tagged_v: tagged images
    // input/output:
    //   used_endpoints_indices: marks which endpoints have been used
    //   trackclusters: holds thrumu tracks made

    if ( m_config.verbosity>0 ) {
      std::cout << "Run Pass #" << passid << ": "
		            << " radialfilter=" << passcfg.run_radial_filter
		            << " linear=" << passcfg.run_linear_tagger
		            << " astar=" << passcfg.run_astar_tagger
		            << std::endl;
    }

    const int nendpts = (int)spacepts.size();
    std::vector< BMTrackCluster3D > pass_track_candidates;
    std::vector< std::vector<int> > pass_end_indices;

    // Declare the vector for the matched track indices up here.
    std::vector< BoundaryFlashIndex > single_track_endpoint_flash_v;
    single_track_endpoint_flash_v.clear();

    std::vector< BoundaryEnd_t > single_track_boundary_type_idx_v;
    single_track_boundary_type_idx_v.clear();

    std::vector< BoundarySpacePoint > single_track_endpoint_v;
    single_track_endpoint_v.clear();

    RadialEndpointFilter radialfilter;

    int num_tracked = 0;

    for (int i=0; i<nendpts; i++) {
      if ( used_endpoints_indices.at(i)==1 ) continue;
      const BoundarySpacePoint& pts_a = *(spacepts[i]);

      bool use_a = true;
      if ( passcfg.run_radial_filter ) {
        int num_segs_a = 0;
        bool within_line_a =  radialfilter.isWithinStraightSegment( pts_a.pos(), img_v, badchimg_v, passcfg.radial_cfg, num_segs_a );
        bool pass_num_segs =  passcfg.radial_cfg.min_segments<=num_segs_a && passcfg.radial_cfg.max_segments>=num_segs_a;
        if ( m_config.verbosity> 1 ) {
          std::cout << "Endpoint #" << i << ": "
                                    << " " << pts_a.printImageCoords( img_v.front().meta() ) << " "
                                    << " within_line=" << within_line_a << " num_segs=" << num_segs_a << " pass_num_segs=" << pass_num_segs
                                    << " run=" << (!within_line_a && pass_num_segs)
                                    << std::endl;
        }

        if ( !within_line_a &&  pass_num_segs )
          use_a = true;
        else
          use_a = false;
      }

      if ( !use_a ) {
        if ( m_config.verbosity>1 ) {
          std::cout << "[ Pass " << passid << " (" << i << ", X). Skip #" << i << " ]" << std::endl;
        }
        continue;
      }

      for (int j=i+1; j<nendpts; j++) {
        if ( used_endpoints_indices.at(j)==1) continue;
        const BoundarySpacePoint& pts_b = *(spacepts[j]);

        if ( pts_a.type()==pts_b.type() ) {
          if ( m_config.verbosity>1 )
            std::cout << "[ Pass " << passid << ":  endpoints (" << i << "," << j << ") skip same type ]" << std::endl;
          continue; // don't connect same type
        }


        float a2b = 0.;
        for (int v=0; v<3; v++) {
          a2b += (pts_a.pos()[v]-pts_b.pos()[v])*(pts_a.pos()[v]-pts_b.pos()[v]);
        }
        a2b = sqrt(a2b);

        if ( a2b<passcfg.min_point_separation ) {
          if ( m_config.verbosity>1 )
            std::cout << "[ Pass " << passid << ":  endpoints (" << i << "," << j << ") below min separation. ]" << std::endl;
          continue;
        }

        bool use_b = true;
        if ( passcfg.run_radial_filter ) {
          int num_segs_b = 0;
          bool within_line_b = radialfilter.isWithinStraightSegment( pts_b.pos(), img_v, badchimg_v, passcfg.radial_cfg, num_segs_b );
          bool pass_num_segs = passcfg.radial_cfg.min_segments<=num_segs_b && passcfg.radial_cfg.max_segments>=num_segs_b;
          if ( !within_line_b && passcfg.radial_cfg.min_segments<=num_segs_b && passcfg.radial_cfg.max_segments>=num_segs_b )
            use_b = true;
          else
            use_b = false;
        }
        if (!use_b) {
          if ( m_config.verbosity>1 ) {
            std::cout << "[ Pass " << passid << " (" << i << "," << j << "). Skip #" << j << " ]" << std::endl;
          }
          continue;
        }

        if ( m_config.verbosity>1 ) {
          std::cout << "[ Pass " << passid << ": path-finding for endpoints (" << i << "," << j << ") "
                    << "of type (" << pts_a.type() << ") -> (" << pts_b.type() << ") ]" << std::endl;
          std::cout << "  start: (" << pts_a.pos()[0] << "," << pts_a.pos()[1] << "," << pts_a.pos()[2] << ") "
                    << "  end: (" << pts_b.pos()[0] << "," << pts_b.pos()[1] << "," << pts_b.pos()[2] << ") " << std::endl;
          std::cout << "  start: " << pts_a.printImageCoords( img_v.front().meta() ) << " end: " << pts_b.printImageCoords( img_v.front().meta() ) << std::endl;
        }

        // empty candidate
        BMTrackCluster3D track3d;
        LinearTaggerInfo linear_result;
        AStarTaggerInfo astar_result;
        FoxTrotExtenderInfo extender_result;
        bool tracked = false;

        // first run the linear charge tagger, if so chosen
        if ( passcfg.run_linear_tagger ) {
          if ( m_config.verbosity>1 ) {
            std::cout << "  Running linear3dchargetagger." << std::endl;
          }
          BMTrackCluster3D linear_track = runLinearChargeTagger( passcfg, pts_a, pts_b, img_v, badchimg_v, linear_result );
          if ( m_config.verbosity>1 ) {
            std::cout << "  good fraction: " << linear_result.goodfrac << std::endl;
            std::cout << "  majority planes w/ charge fraction: " << linear_result.majfrac << std::endl;
	    if ( linear_track.path3d.size()>0 ) {
	      const std::vector<double>& retstart = linear_track.path3d.front();
	      const std::vector<double>& retend = linear_track.path3d.back();
	      std::cout << "  " << retstart.size() << " " << retend.size() << std::endl;
	      std::cout << "  returned start=(" << retstart[0] << "," << retstart[1] << "," << retstart[2] << ")" << std::endl;
	      std::cout << "  returned end=(" << retend[0] << "," << retend[1] << "," << retend[2] << ")" << std::endl;
	    }
          }
          if ( linear_result.isgood ) {
            if ( m_config.verbosity>1 ) std::cout << "  Linear Result is good. length=" << linear_track.path3d.size() << std::endl;
            std::swap(track3d,linear_track);
          }
          else {
            if ( m_config.verbosity>1 ) std::cout << "  Result is bad." << std::endl;
          }
          tracked = true;
        }//end of if run_linear

        // next run astar tagger
        if ( passcfg.run_astar_tagger
             && (linear_result.goodfrac>passcfg.astar3d_min_goodfrac || linear_result.majfrac>passcfg.astar3d_min_majfrac ) ) {
          if ( m_config.verbosity>1 )
            std::cout << "  Running astar tagger." << std::endl;
          BMTrackCluster3D astar_track = runAStarTagger( passcfg, pts_a, pts_b, img_v, badchimg_v, astar_result );
          tracked = true;
          if ( astar_result.isgood ) {
            if ( m_config.verbosity>1 )  {
              std::vector<double>& astar_end = astar_track.path3d.back();
              std::cout << "  AStar end point: (" << astar_end[0] << "," << astar_end[1] << "," << astar_end[2] << ")" << std::endl;
              std::cout << "  AStar Result is good." << std::endl;
              if ( astar_result.goal_reached )
                std::cout << "  AStar reached goal." << std::endl;
              else
                std::cout << "  AStar did not reach goal." << std::endl;
            }
            std::swap(track3d,astar_track);
          }
          else {
            if ( m_config.verbosity>1 ) std::cout << "  AStar Result is bad." << std::endl;
          }
        }

        if ( tracked )
          num_tracked++;

        // run bezier fitter (doesn't exist yet. write this?)

        // could add criteria to filter here
        bool accept_track = false;
        if ( (passcfg.run_linear_tagger && !passcfg.run_astar_tagger) && linear_result.isgood ) {
          // if we only run the linear tagger, then we save
          accept_track = true;
        }
        else if ( (passcfg.run_linear_tagger && passcfg.run_astar_tagger) && astar_result.isgood && astar_result.goal_reached ) {
          // when we run both, we use the linear tagger as a prefilter for astar, which can be expensive
          accept_track = true;
        }

        if ( accept_track ) {
          if ( m_config.verbosity>1 ) {
            std::cout << "  track found. "
                      << " length: " << track3d.path3d.size()
                      << " empty=" << track3d.isempty()
                      << " linear-good=" << linear_result.isgood
                      << " astar-good=" << astar_result.isgood
                      << " astar-reached=" << astar_result.goal_reached
                      << std::endl;
          }
          // extend the good track
          runFoxTrotExtender( passcfg, track3d.path3d, img_v, badchimg_v, tagged_v, extender_result);
        }

        if ( accept_track && !track3d.isempty() ) {
          std::vector<int> indices(2);
          indices[0] = i;
          indices[1] = j;

	  // Declare a vector for the flash index and producer index at the ith and jth points in the 'flash_idx_v' and 'boundary_type_idx_v' vectors.
	  // 'i' necessarily must be less than j because 'j' begins to iterate at 'i+1'.
	  single_track_endpoint_flash_v.resize(2);
	  single_track_endpoint_flash_v[0]          = spacepts.at( i )->getFlashIndex();
	  single_track_endpoint_flash_v[1]          = spacepts.at( j )->getFlashIndex();
	  track_endpoint_flash_v.emplace_back( std::move(single_track_endpoint_flash_v) );

	  single_track_boundary_type_idx_v.resize(2);
	  single_track_boundary_type_idx_v[0] = spacepts.at( i )->type();
	  single_track_boundary_type_idx_v[1] = spacepts.at( j )->type();
	  track_endpoint_boundary_type_idx_v.emplace_back( std::move(single_track_boundary_type_idx_v) );

	  // Save the information for the boundary points in this vector.
	  single_track_endpoint_v.resize(2);
	  single_track_endpoint_v[0] = *spacepts.at( i );
	  single_track_endpoint_v[1] = *spacepts.at( j );
	  track_endpoint_v.emplace_back( std::move(single_track_endpoint_v) );
	  
          if ( m_config.verbosity>1 ) {
            std::cout << "#### Storing track. size=" << track3d.path3d.size() << ". indices (" << indices[0] << "," << indices[1] << ") ####" << std::endl;
          }
          pass_end_indices.emplace_back( std::move(indices) );
          pass_track_candidates.emplace_back( std::move(track3d) );
	  
        }

      }// second loop over end points
    }// first loop over end points

    if ( m_config.verbosity>0 ) {
      std::cout << "Pass #" << passid << ": "
                << " number of pairs tracked: " << num_tracked
                << " number of tracks passed: " << pass_track_candidates.size()
                << " number of indice pairs: " << pass_end_indices.size()
                << std::endl;
    }


    // post-process tracks (nothing yet)
    for ( auto &track : pass_track_candidates ) {
      trackclusters.emplace_back( std::move(track) );
    }
    for ( auto&indices : pass_end_indices ) {
      used_endpoints_indices.at(indices[0]) = 1;
      used_endpoints_indices.at(indices[1]) = 1;
      
    }

    return;
  }


  BMTrackCluster3D ThruMuTracker::runLinearChargeTagger( const ThruMuTrackerConfig::ThruMuPassConfig& pass_cfg,
                                                         const BoundarySpacePoint& pts_a, const BoundarySpacePoint& pts_b,
                                                         const std::vector<larcv::Image2D>& img_v,
                                                         const std::vector<larcv::Image2D>& badchimg_v,
                                                         ThruMuTracker::LinearTaggerInfo& result_info ) {

    // linear 3D track
    Linear3DChargeTagger linetrackalgo( pass_cfg.linear3d_cfg ); // finds charge along a line in 3D space

    // get pixel start and end points
    std::vector<int> a_cols(img_v.size(),0);
    std::vector<int> b_cols(img_v.size(),0);
    for (int p=0; p<(int)img_v.size(); p++) {
      a_cols[p] = pts_a.at(p).col;
      b_cols[p] = pts_b.at(p).col;
    }

    PointInfoList straight_track = linetrackalgo.findpath( img_v, badchimg_v, pts_a.at(0).row, pts_b.at(0).row, a_cols, b_cols );
    result_info.numpts   = straight_track.size();
    result_info.goodfrac = straight_track.fractionGood();
    result_info.majfrac  = straight_track.fractionHasChargeOnMajorityOfPlanes();
    if ( result_info.numpts   > pass_cfg.linear3d_min_tracksize
          && (result_info.goodfrac > pass_cfg.linear3d_min_goodfraction
	        || result_info.majfrac  > pass_cfg.linear3d_min_majoritychargefraction) ) {
      result_info.isgood = true;
    }
    else {
      result_info.isgood = false;
    }

    if ( !result_info.isgood ) {
      BMTrackCluster3D emptytrack;
      return emptytrack;
    }

    std::vector< std::vector<double> > path3d;
    for ( auto const& ptinfo : straight_track ) {
      std::vector<double> xyz(3);

      // Declare a float value for the same vector to use that.
      std::vector<float> xyz_float(3,0.0);
      
      for (int i=0; i<3; i++) {
      	xyz[i] = ptinfo.xyz[i];
        xyz_float[i] = ptinfo.xyz[i];
      }

      // Ensure that the 3 coordinates give a value inside the image before appending this value to the track.
      // Convert the 3D position to a 2D pixel image.
      std::vector<int> imgcoords = ublarcvapp::UBWireTool::getProjectedImagePixel( xyz_float, img_v.front().meta(), img_v.size() );

      // Declare an instance of type 'PushBoundarySpacePoint'.
      PushBoundarySpacePoint point_push;

      bool pixel_in_image = point_push.isPixelWithinImage(img_v, imgcoords);

      if (pixel_in_image)
	path3d.emplace_back( std::move(xyz) );

    }

    BMTrackCluster3D track3d( pts_a, pts_b, path3d );
    
    return track3d;
  }

  BMTrackCluster3D ThruMuTracker::runAStarTagger( const ThruMuTrackerConfig::ThruMuPassConfig& pass_cfg,
                                                  const BoundarySpacePoint& pts_a, const BoundarySpacePoint& pts_b,
                                                  const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimg_v,
                                                  ThruMuTracker::AStarTaggerInfo& result_info ) {

    // collect meta/translate start/goal tick/wires to row/col in compressed image
    std::vector< const larcv::ImageMeta* > meta_compressed_v;
    std::vector<int> start_cols(m_img_compressed_v.size(),0);
    std::vector<int> start_rows(m_img_compressed_v.size(),0);
    std::vector<int> goal_cols(m_img_compressed_v.size(),0);
    std::vector<int> goal_rows(m_img_compressed_v.size(),0);
    std::vector< const larcv::Pixel2D* > start_pix;
    std::vector< const larcv::Pixel2D* > goal_pix;


    for (size_t p=0; p<m_img_compressed_v.size(); p++) {

      // get start/end point informatio in compressed image
      const larcv::ImageMeta* ptr_meta = &(m_img_compressed_v[p].meta()); // compressed image meta
      const larcv::ImageMeta& meta = img_v[p].meta(); // original image meta

      const BoundaryEndPt& start_endpt = pts_a[p];
      start_rows[p] =  ptr_meta->row( meta.pos_y( start_endpt.row ) );
      start_cols[p] =  ptr_meta->col( meta.pos_x( start_endpt.col ) );
      start_pix.push_back( new larcv::Pixel2D( start_endpt.getcol(), start_endpt.getrow() ) );

      const BoundaryEndPt& goal_endpt  = pts_b[p];
      goal_rows[p]  =  ptr_meta->row( meta.pos_y( goal_endpt.row ) );
      goal_cols[p]  =  ptr_meta->col( meta.pos_x( goal_endpt.col ) );
      goal_pix.push_back( new larcv::Pixel2D( goal_endpt.getcol(), goal_endpt.getrow() ) );

      meta_compressed_v.push_back( ptr_meta );
    }

    reco3d::AStar3DAlgo algo( pass_cfg.astar3d_cfg );
    std::vector<reco3d::AStar3DNode> path;
    result_info.goal_reached = false;
    result_info.isgood = true;
    int goalhit = 0;
    try {
      path = algo.findpath( m_img_compressed_v, m_badch_compressed_v, m_badch_compressed_v, // tagged_compressed_v
			    start_rows.front(), goal_rows.front(), start_cols, goal_cols, goalhit );

      if ( goalhit==1 ) {
        result_info.goal_reached = true;
      }
    }
    catch (const std::exception& e) {
      std::cout << "*** [ Exception running astar3dalgo::findpath: " << e.what() << "] ***" << std::endl;
      result_info.isgood = false;
      result_info.goal_reached = false;
      result_info.nbad_nodes = -1;
      result_info.total_nodes = -1;
      std::vector< std::vector<double> > empty_path3d;
      BMTrackCluster3D track3d( pts_a, pts_b, empty_path3d );
      for (int p=0; p<3; p++ ) {
        delete start_pix.at(p);
        delete goal_pix.at(p);
      }
      return track3d; // return empty track
    }

    result_info.nbad_nodes = 0;
    result_info.total_nodes = 0;
    for ( auto& node : path ) {
      if ( node.badchnode )
        result_info.nbad_nodes+=1.0;
      result_info.total_nodes+=1.0;
    }

    // if majority are bad ch nodes, reject this track
    result_info.frac_bad = float(result_info.nbad_nodes)/float(result_info.total_nodes);
    if ( result_info.frac_bad>0.5 || result_info.total_nodes<=3)
      result_info.isgood = false;

    BMTrackCluster3D track3d = AStarNodes2BMTrackCluster3D( path, img_v, pts_a, pts_b, 0.3 );

    for (int p=0; p<3; p++ ) {
      delete start_pix.at(p);
      delete goal_pix.at(p);
    }

    return track3d;
  }

  void ThruMuTracker::runFoxTrotExtender( const ThruMuTrackerConfig::ThruMuPassConfig& pass_cfg,
                                          std::vector<std::vector<double> >& track,
                                          const std::vector<larcv::Image2D>& img_v,
                                          const std::vector<larcv::Image2D>& badch_v,
                                          const std::vector<larcv::Image2D>& tagged_v,
                                          ThruMuTracker::FoxTrotExtenderInfo& result_info ) {
    // try to extend the track.
    // make sure we don't go backwards

    if ( !pass_cfg.run_foxtrot_extender || track.size()<2 ) {
      result_info.isgood = false;
      return;
    }

    ThruMuFoxExtender extender_algo( pass_cfg.foxextend_cfg );

    // need a forward and backward extension
    result_info.isgood  = extender_algo.extendTrack( track, img_v, badch_v, tagged_v );

    return;
  }

}
}
