#include "MRCNNMatch.h"

#include <algorithm>

#include "ublarcvapp/UBImageMod/EmptyChannelAlgo.h"
#include "CropMaskCombo.h"

namespace ublarcvapp {
namespace dltagger {

  /**
   * match mask R-CNN masks across planes.
   *
   * the output of the algorithm is stored in the m_combo_XXXX_v
   * 
   * @param[in] clustermask_vv vector of masks for each plane. output of Mask-RCNN.
   * @param[in] wholeview_v whole images for each plane.
   * @param[in] ev_chstatus channel status info for event.
   *
   */
  void MRCNNMatch::matchMasksAcrossPlanes( const std::vector<std::vector<larcv::ClusterMask>>& clustermask_vv,
                                           const std::vector<larcv::Image2D>& wholeview_v,
                                           const larcv::EventChStatus& ev_chstatus,
                                           bool use_gap_ch ) {

    // clear container variables
    clear();
    
    // make badch image
    ublarcvapp::EmptyChannelAlgo badchalgo;
    auto const& wholeimage_meta = wholeview_v.front().meta();
    
    std::vector<larcv::Image2D> badch_v =
      badchalgo.makeBadChImage( 4, wholeview_v.size(),
                                wholeimage_meta.min_y(),
                                wholeimage_meta.rows()*wholeimage_meta.pixel_height(),
                                wholeimage_meta.cols()*wholeimage_meta.pixel_width(),
                                wholeimage_meta.pixel_height(),
                                wholeimage_meta.pixel_width(),                                
                                ev_chstatus );

    std::vector<larcv::Image2D> gapch_v = badchalgo.findMissingBadChs( wholeview_v, badch_v, 1.0, 200 );
    

    // compile key match criteria for each mask
    LARCV_DEBUG() << "compile MaskMatchData" << std::endl;    
    m_matchdata_vv.clear();
    for ( auto const& clustermask_v : clustermask_vv ) {
      std::vector<MaskMatchData> data_v;
      int idx=0;
      for ( auto const& mask : clustermask_v ) {
        MaskMatchData data( mask.meta.plane(), idx, mask );
        data_v.emplace_back( std::move(data) );
        idx++;
      }
      std::sort( data_v.begin(), data_v.end() );
      LARCV_INFO() << "plane[" << clustermask_v.front().meta.plane() << "] number of masks: " << data_v.size() << std::endl;
      if ( logger().debug() ) {
        // print list of masks and the data
        for ( auto& data : data_v )
          LARCV_DEBUG() << "  " << data << std::endl;
      }
      m_matchdata_vv.emplace_back( std::move(data_v) );
    }

    if ( use_gap_ch ) {
      // add the bad ch
      for ( size_t p=0; p<3; p++ ) {
        auto const& gapch = gapch_v[p];
        auto& badch = badch_v[p];
        for ( size_t c=0; c<gapch.meta().cols(); c++ ) {
          if ( gapch.pixel(0,c)>0 && badch.pixel(0,c)<1 )
            badch.paint_col( c, 50 );
        }
      }
    }

    // make 3-plane matches
    run3PlanePass( clustermask_vv, wholeview_v, badch_v, m_matchdata_vv );
    int npassed = 0;
    for ( auto& isgood : m_pass ) {
      if ( isgood==1 ) npassed++;
    }
    
    
    LARCV_INFO() << "------------------------------------------------" << std::endl;
    LARCV_INFO() << "AFTER 3 PLANE MATCHES" << std::endl;
    LARCV_INFO() << "number of reco'd tracks: " << npassed << std::endl;
    LARCV_INFO() << "Masked used in matches" << std::endl;
    std::vector<int> used_v( 3, 0 );
    for ( size_t p=0; p<3; p++ ) {
      for ( auto const& matchdata : m_matchdata_vv[p] ) {
        if ( matchdata.used ) used_v[p]++;
      }
      LARCV_INFO() << "  Plane[" << p<< "]: " << used_v[p] << " of " << m_matchdata_vv[p].size() << std::endl;
    }
    LARCV_INFO() << "------------------------------------------------" << std::endl;
    
    // make 2-plane matches
    run2PlanePass( clustermask_vv, wholeview_v, badch_v, m_matchdata_vv );
    npassed = 0;
    for ( auto& isgood : m_pass ) {
      if ( isgood==1 ) npassed++;
    }    
    LARCV_INFO() << "------------------------------------------------" << std::endl;
    LARCV_INFO() << "AFTER 2 PLANE MATCHES" << std::endl;
    LARCV_INFO() << "number of reco'd tracks: " << npassed << std::endl;
    LARCV_INFO() << "Masked used in matches" << std::endl;

    for ( size_t p=0; p<3; p++ ) {
      used_v[p] = 0;      
      for ( auto const& matchdata : m_matchdata_vv[p] ) {
        if ( matchdata.used ) used_v[p]++;
      }
      LARCV_INFO() << "  Plane[" << p<< "]: " << used_v[p] << " of " << m_matchdata_vv[p].size() << std::endl;
    }
    LARCV_INFO() << "------------------------------------------------" << std::endl;


  }

  /**
   * match masks across all three planes.
   *
   * serves as first pass
   *
   * @param[in] clustermask_vv  input cluster masks for each plane.
   * @param[in] wholeview_v     whole image for each plane.
   * @param[in] badch_v         whole image indicating location of bad channels. Made by ublarcvapp/UBImageMod/EmptyChannelAlgo
   * @param[inout] matchdata_vv info on mask bounds. also contains flag if mask is claimed, which this method can change.
   */
  void MRCNNMatch::run3PlanePass( const std::vector<std::vector<larcv::ClusterMask>>& clustermask_vv,
                                  const std::vector<larcv::Image2D>& wholeview_v,
                                  const std::vector<larcv::Image2D>& badch_v,
                                  std::vector< std::vector<MaskMatchData> >& matchdata_vv ) {

    /// initial screen for overlapping times
    // start from Y plane, match U, then V
    LARCV_DEBUG() << "make 3 plane matches" << std::endl;

    std::vector<MaskCombo> combo_v;
    
    for ( size_t iy=0; iy<matchdata_vv.at(2).size(); iy++ ) {
      auto const& data_yplane = matchdata_vv.at(2).at(iy);
      auto const& mask_yplane = clustermask_vv.at(2).at( data_yplane.index );

      for ( size_t iu=0; iu<matchdata_vv.at(0).size(); iu++ ) {

        auto const& data_uplane = matchdata_vv.at(0).at(iu);
        auto const& mask_uplane = clustermask_vv.at(0).at(data_uplane.index);


        if ( data_uplane.tick_min > data_yplane.tick_max ) {
          break; // since we sorted, we can cut off
        }

        // create a combo filled with Y-plane to start
        MaskCombo combo;
        combo.addMask( mask_yplane, data_yplane );
        combo.maskdata_indices[2] = iy;

        if ( combo.iscompatible( data_uplane ) ) {
          // add u-plane
          combo.addMask( mask_uplane, data_uplane );
          combo.maskdata_indices[0] = iu;
          //std::cout << "y-u match: " << combo << std::endl;
          // store
          combo_v.emplace_back( std::move(combo) );
        }
      }//end of u-plane loop

    }
    std::sort( combo_v.begin(), combo_v.end() );
    
    LARCV_DEBUG() << "combos after Y-U match: " << combo_v.size() << std::endl;    
    //for ( auto const& combo : combo_v )
    //  std::cout << "  " << combo << std::endl;

    // pair with v-plane
    std::vector< MaskCombo > combo_3plane_v;
    
    for ( auto const& combo : combo_v ) {

      for ( size_t iv=0; iv<matchdata_vv.at(1).size(); iv++ ) {
        auto const& data_vplane = matchdata_vv.at(1).at(iv);
        auto const& mask_vplane = clustermask_vv.at(1).at(data_vplane.index);

        if ( combo.union_tick[1] < data_vplane.tick_min ) break;

        if ( combo.iscompatible( data_vplane ) ) {

          MaskCombo combo_yuv( combo );

          combo_yuv.addMask( mask_vplane, data_vplane );
          combo_yuv.maskdata_indices[1] = iv;          

          LARCV_DEBUG() << "combo Y-U-V: " << combo_yuv << std::endl;
          if ( combo_yuv.iou()>0.20 )
            combo_3plane_v.emplace_back( std::move(combo_yuv) );
          
        }
      }// loop over v-plane cluster/mask

    }// end of combo y-u loop

    std::sort( combo_3plane_v.begin(), combo_3plane_v.end() );
    //std::cout << "combos after Y-U-V match: " << m_combo_3plane_v.size() << std::endl;
    //for ( auto const& combo : m_combo_3plane_v )
    //std::cout << "  " << combo << std::endl;

    // make crops, mask charge, contours on charge, contours on mask
    for ( size_t icombo=0; icombo<combo_3plane_v.size(); icombo++ ) {
      auto& combo = combo_3plane_v[icombo];
      LARCV_INFO() << "----------------------------------" << std::endl;
      LARCV_INFO() << "Evaluate 3-Plane Combo[" << icombo << "] "
                   << "maskindices="
                   << "(" << combo.maskdata_indices[0] << "," << combo.maskdata_indices[1] << "," << combo.maskdata_indices[2] << ")"
                   << std::endl;
      LARCV_INFO() << "----------------------------------" << std::endl;      
      

      // greedy: if already used successfully, we do not reuse
      bool isused = false;
      for ( size_t p=0; p<combo.maskdata_indices.size(); p++ ) {
        if ( matchdata_vv[p][ combo.maskdata_indices[p] ].used ) isused = true;
      }
      if ( isused )
        continue;
      
      // we step through a number of data products

      // make crop around mask. make image of both charge and mask pixels
      CropMaskCombo     cropmaker( combo, wholeview_v, badch_v );
      // extract contours, PCA of mask pixels. hold off extending mask.
      FeaturesMaskCombo features( cropmaker, false );
      // attempt to define 3D endpoints
      Gen3DEndpoints    endpoints( features );

      bool good_tri_area = true;
      for ( auto const& triarea_score : endpoints.endpt_tri_v )
        if ( triarea_score>_config.triarea_maximum ) good_tri_area = false;

      bool isochronous = false;
      if ( cropmaker.crops_v.front().meta().height()<400.0 ) {
        isochronous = true;
        LARCV_INFO() << "  combo evaluated as isochronous" << std::endl;
      }
      
      if ( good_tri_area ) {
        features.extendMaskWithPCAregion(10.0);
        LARCV_INFO() << "  combo has consistency end points built from mask PCA" << std::endl;
      }
      
      // scan for 3D points for later AStar graph
      GenGraphPoints    graphpoints( features, endpoints, larcv::msg::kNORMAL );
      LARCV_INFO() << "  max good pixel gap from graph reco: " << graphpoints.m_maxgapdist << std::endl;

      // run astar only if the 3d points are fairly consistent
      //bool runastar = good_tri_area;
      bool runastar = false;
      
      AStarMaskCombo    astar( endpoints, badch_v, runastar,
                               _config.astar_max_downsample_factor,
                               _config.astar_store_score_image,
                               _config.astar_verbosity );

      if ( good_tri_area || isochronous ) {
        
        if (astar.astar_completed==1 || graphpoints.m_maxgapdist < 50.0 ) {
          m_pass.push_back(1);
          LARCV_INFO() << "  combo passes reco." << std::endl;
          // we mark the mask data in the combo as used
          for ( size_t p=0; p<combo.maskdata_indices.size(); p++ )
            matchdata_vv[p][ combo.maskdata_indices[p] ].used = true;
          
        }
        else {
          LARCV_INFO() << "  combo fails reco." << std::endl;
          m_pass.push_back(0);
        }

        // to kep allignment, if pass fills, we fill objects
        m_combo_3plane_v.emplace_back( std::move(combo) );
        m_combo_crops_v.emplace_back( std::move(cropmaker) );
        m_combo_features_v.emplace_back( std::move(features) );
        m_combo_graphpts_v.emplace_back( std::move(graphpoints) );
        m_combo_endpt3d_v.emplace_back( std::move(endpoints) );
        m_combo_astar_v.emplace_back( std::move(astar) );

      }
      else {
        LARCV_INFO() << "  combo is not recod." << std::endl;
      }

      // for debug
      //if ( m_combo_3plane_v.size()>=3 )
      //break;
    }

    LARCV_INFO() << "Number of 3-plane matches: " << m_combo_3plane_v.size() << std::endl;
    
  }

  /**
   * match across 2 planes only.
   *
   * intended for a second-pass
   *
   * @param[in] clustermask_vv vector of cluster masks for each plane.
   * @param[in] wholeview_v    whole ADC image for each plane.
   * @param[in] badch_v        whole view image where bad channels are marked.
   * @param[inout] matchdata_v data for each cluster masks used for matching. an element is flagged after being used.
   *
   */
  void MRCNNMatch::run2PlanePass( const std::vector<std::vector<larcv::ClusterMask>>& clustermask_vv,
                                  const std::vector<larcv::Image2D>& wholeview_v,
                                  const std::vector<larcv::Image2D>& badch_v,
                                  std::vector< std::vector<MaskMatchData> >& matchdata_vv ) {

    // assemble 2-plane combos
    std::vector<MaskCombo> combo_v;
    
    for ( size_t iy=0; iy<matchdata_vv.at(2).size(); iy++ ) {
      auto const& data_yplane = matchdata_vv.at(2).at(iy);
      auto const& mask_yplane = clustermask_vv.at(2).at( data_yplane.index );

      if ( data_yplane.used )
        continue; // don't reuse
      
      std::vector<int> indices(3);
      indices[2] = iy;

      for ( size_t iu=0; iu<matchdata_vv.at(0).size(); iu++ ) {

        auto const& data_uplane = matchdata_vv.at(0).at(iu);
        auto const& mask_uplane = clustermask_vv.at(0).at(data_uplane.index);

        if ( data_uplane.used )
          continue; // don't reuse

        if ( data_uplane.tick_min > data_yplane.tick_max ) {
          break; // since we sorted, we can cut off
        }

        // create a combo filled with Y-plane to start
        MaskCombo combo;
        combo.addMask( mask_yplane, data_yplane );
        combo.maskdata_indices[0] = iu;
        combo.maskdata_indices[2] = iy;

        if ( combo.iscompatible( data_uplane ) ) {
          // add u-plane
          combo.addMask( mask_uplane, data_uplane );
          //std::cout << "y-u match: " << combo << std::endl;
          // store
          combo_v.emplace_back( std::move(combo) );
        }
      }//end of u-plane loop

      // match with v-plane
      for ( size_t iv=0; iv<matchdata_vv.at(1).size(); iv++ ) {

        auto const& data_vplane = matchdata_vv.at(1).at(iv);
        auto const& mask_vplane = clustermask_vv.at(1).at(data_vplane.index);

        if ( data_vplane.used )
          continue; // don't reuse

        if ( data_vplane.tick_min > data_yplane.tick_max ) {
          break; // since we sorted, we can cut off
        }

        // create a combo filled with Y-plane to start
        MaskCombo combo;
        combo.addMask( mask_yplane, data_yplane );
        combo.maskdata_indices[1] = iv;
        combo.maskdata_indices[2] = iy;

        if ( combo.iscompatible( data_vplane ) ) {
          // add u-plane
          combo.addMask( mask_vplane, data_vplane );
          //std::cout << "y-u match: " << combo << std::endl;
          // store
          combo_v.emplace_back( std::move(combo) );
        }
      }//end of v-plane loop
      
    }
    std::sort( combo_v.begin(), combo_v.end() );
    LARCV_DEBUG() << "NUM OF 2-PLANE COMBOS (2nd pass): " << combo_v.size() << std::endl;
    if ( logger().debug() ) {
      for ( auto const& combo : combo_v ) {
        std::cout << combo << std::endl;
      }
    }
    
    // Analyze 2-plane combos
    int npassed = 0;
    int icombo = -1;
    for ( auto const& combo : combo_v ) {
      icombo++;
      LARCV_INFO() << "----------------------------------" << std::endl;
      LARCV_INFO() << "Evaluate 2-Plane Combo[" << icombo << "]" << std::endl;
      LARCV_INFO() << "----------------------------------" << std::endl;      

      // greedy: if already used successfully, we do not reuse
      bool isused = false;
      for ( size_t p=0; p<combo.maskdata_indices.size(); p++ ) {
        if ( combo.maskdata_indices[p]==-1 ) continue;
        if ( matchdata_vv[p][ combo.maskdata_indices[p] ].used ) isused = true;
      }
      if ( isused ) {
        LARCV_INFO() << "masks in combo are used." << std::endl;        
        continue;
      }

      // make crop around mask. make image of both charge and mask pixels
      CropMaskCombo     cropmaker( combo, wholeview_v, badch_v );
      // extract contours, PCA of mask pixels
      FeaturesMaskCombo features( cropmaker, true );
      // attempt to define 3D endpoints
      Gen3DEndpoints    endpoints( features );

      bool isochronous = false;
      if ( cropmaker.crops_v.front().meta().height()<400.0 ) {
        isochronous = true;
        LARCV_INFO() << "  combo evaluated as isochronous" << std::endl;
      }
      
      bool run = true;
      if (endpoints.endpt_tpc_v[0]==0 || endpoints.endpt_tpc_v[1]==0) {
        run = false;
        LARCV_INFO() << "  skipping astar for this combo (no end points)" << std::endl;
      }
      else {
        LARCV_INFO() << "  running astar for this combo" << std::endl;
      }

      features.extendMaskWithPCAregion(10.0);

      LARCV_INFO() << "  run graph-based track reco" << std::endl;
      GenGraphPoints    graphpoints( features, endpoints, larcv::msg::kNORMAL );
      LARCV_INFO() << "  max good pixel gap from graph reco: " << graphpoints.m_maxgapdist << std::endl;
      run = false;
      
      AStarMaskCombo    astar( endpoints, badch_v, run,
                               _config.astar_max_downsample_factor,
                               _config.astar_store_score_image,
                               _config.astar_verbosity );
      //astar.set_verbosity((larcv::msg::Level_t)0);

      bool didpass = false;
      if (astar.astar_completed==1 || graphpoints.m_maxgapdist<50.0 ) {
        m_pass.push_back(1);
        npassed++;
        didpass = true;
        
        // we mark the mask data in the combo as used
        for ( size_t p=0; p<combo.maskdata_indices.size(); p++ ) {
          if ( combo.maskdata_indices[p]!=-1 )
            matchdata_vv[p][ combo.maskdata_indices[p] ].used = true;
        }

        // for 2 plane matches, we won't have features, which we need later for pixel tagging. we do that now
        m_combo_3plane_v.emplace_back( std::move(combo) );
        m_combo_crops_v.emplace_back( std::move(cropmaker) );
        m_combo_features_v.emplace_back( std::move(features) );
        m_combo_graphpts_v.emplace_back( std::move(graphpoints) );        
        m_combo_endpt3d_v.emplace_back( std::move(endpoints) );
        m_combo_astar_v.emplace_back( std::move(astar) );
      }
      
      if ( didpass==1 ) 
        LARCV_INFO() << " combo passes reco." << std::endl;
      else
        LARCV_INFO() << " combo fails reco" << std::endl;
      
      
    }//end of combo loop
    
    LARCV_INFO() << "Number of 2-plane matches: " << npassed << std::endl;
    
  }
  
  /**
   * clear output containers
   *
   */
  void MRCNNMatch::clear() {
    m_matchdata_vv.clear();
    m_combo_3plane_v.clear();
    m_combo_crops_v.clear();
    m_combo_features_v.clear();
    m_combo_endpt3d_v.clear();
    m_combo_graphpts_v.clear();
    m_combo_astar_v.clear();
    m_pass.clear();
  }
  
}
}
