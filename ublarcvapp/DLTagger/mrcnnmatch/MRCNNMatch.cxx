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
                                           const larcv::EventChStatus& ev_chstatus ) {

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
    

    // compile key match criteria for each mask
    std::vector< std::vector<MaskMatchData> > matchdata_vv;
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
      matchdata_vv.emplace_back( std::move(data_v) );
    }

    // make 3-plane matches
    run3PlanePass( clustermask_vv, wholeview_v, badch_v, matchdata_vv );
    // make 2-plane matches
    run2PlanePass( clustermask_vv, wholeview_v, badch_v, matchdata_vv );
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
    std::vector<int> pass;
    for ( size_t icombo=0; icombo<combo_3plane_v.size(); icombo++ ) {

      auto& combo = combo_3plane_v[icombo];

      // greedy: if already used successfully, we do not reuse
      bool isused = false;
      for ( size_t p=0; p<combo.maskdata_indices.size(); p++ ) {
        if ( matchdata_vv[p][ combo.maskdata_indices[p] ].used ) isused = true;
      }
      if ( isused )
        continue;
      
      // we step through a number of data products

      // make crop around mask. make image of both charge and mask pixels
      CropMaskCombo     cropmaker( combo, wholeview_v );
      // extract contours, PCA of mask pixels
      FeaturesMaskCombo features( cropmaker );
      // attempt to define 3D endpoints
      Gen3DEndpoints    endpoints( features );

      // run astar only if the 3d points are fairly consistent
      bool runastar = true;
      for ( auto const& triarea_score : endpoints.endpt_tri_v )
        if ( triarea_score>_config.triarea_maximum ) runastar = false;
      
      AStarMaskCombo    astar( endpoints, badch_v, runastar,
                               _config.astar_max_downsample_factor,
                               _config.astar_store_score_image,
                               _config.astar_verbosity );

      if ( runastar ) {
        if (astar.astar_completed==1 ) {
          pass.push_back(1);

          // we mark the mask data in the combo as used
          for ( size_t p=0; p<combo.maskdata_indices.size(); p++ )
            matchdata_vv[p][ combo.maskdata_indices[p] ].used = true;
          
        }
        else {
          pass.push_back(0);
        }

        if ( pass.back()==1 || !_config.filter_astar_failures ) {
          m_combo_3plane_v.emplace_back( std::move(combo) );
          m_combo_crops_v.emplace_back( std::move(cropmaker) );
          m_combo_features_v.emplace_back( std::move(features) );
          m_combo_endpt3d_v.emplace_back( std::move(endpoints) );
          m_combo_astar_v.emplace_back( std::move(astar) );
        }
      }

      // for debug
      //if ( icombo==0 )
      //  break;
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
    std::vector<int> pass;
    int npassed = 0;
    for ( auto const& combo : combo_v ) {

      // greedy: if already used successfully, we do not reuse
      bool isused = false;
      for ( size_t p=0; p<combo.maskdata_indices.size(); p++ ) {
        if ( combo.maskdata_indices[p]==-1 ) continue;
        if ( matchdata_vv[p][ combo.maskdata_indices[p] ].used ) isused = true;
      }
      if ( isused ) {
        //std::cout << "combo is used." << std::endl;        
        continue;
      }

      // make crop around mask. make image of both charge and mask pixels
      CropMaskCombo     cropmaker( combo, wholeview_v );
      // extract contours, PCA of mask pixels
      FeaturesMaskCombo features( cropmaker );
      // attempt to define 3D endpoints
      Gen3DEndpoints    endpoints( features );

      bool run = true;
      if (endpoints.endpt_tpc_v[0]==0 || endpoints.endpt_tpc_v[1]==0 )
        run = false;
      
      // astar
      AStarMaskCombo    astar( endpoints, badch_v, run,
                               _config.astar_max_downsample_factor,
                               _config.astar_store_score_image,
                               _config.astar_verbosity );
                               
      if ( run ) {
        if (astar.astar_completed==1 ) {
          pass.push_back(1);
          npassed++;
          
          // we mark the mask data in the combo as used
          for ( size_t p=0; p<combo.maskdata_indices.size(); p++ ) {
            if ( combo.maskdata_indices[p]!=-1 )
              matchdata_vv[p][ combo.maskdata_indices[p] ].used = true;
          }
        }
        else {
          pass.push_back(0);
        }

        if ( pass.back()==1 || !_config.filter_astar_failures ) {
          m_combo_3plane_v.emplace_back( std::move(combo) );
          m_combo_crops_v.emplace_back( std::move(cropmaker) );
          m_combo_features_v.emplace_back( std::move(features) );
          m_combo_endpt3d_v.emplace_back( std::move(endpoints) );
          m_combo_astar_v.emplace_back( std::move(astar) );
          //std::cout << "STORE 2-PLANE COMBO" << std::endl;
        }
      }

    }//end of combo loop

    LARCV_INFO() << "Number of 2-plane matches: " << npassed << std::endl;

  }

}
}
