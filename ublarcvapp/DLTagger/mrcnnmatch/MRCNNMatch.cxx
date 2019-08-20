#include "MRCNNMatch.h"

#include <algorithm>

#include "ublarcvapp/UBImageMod/EmptyChannelAlgo.h"
#include "CropMaskCombo.h"

namespace ublarcvapp {
namespace dltagger {

  /**
   * match mask R-CNN masks across planes
   * 
   * @param[in] clustermask_vv vector of masks for each plane
   * @param[inout] match_indices each match stores as a vector of mask indices for each plane. -1 if no match for that plane
   *
   */
  void MRCNNMatch::matchMasksAcrossPlanes( const std::vector<std::vector<larcv::ClusterMask>>& clustermask_vv,
                                           const std::vector<larcv::Image2D>& wholeview_v,
                                           const larcv::EventChStatus& ev_chstatus,
                                           std::vector< std::vector<int> >& match_indices ) {

    // make badch image
    ublarcvapp::EmptyChannelAlgo badchalgo;
    std::vector<larcv::Image2D> badch_v =
      badchalgo.makeBadChImage( 4, 3, 2400, 1008*6, 3456, 6, 1, ev_chstatus );
    
    
    // first we compile key match criteria for each mask

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
      std::cout << "PLANE " << clustermask_v.front().meta.plane() << " MASKS" << std::endl;
      for ( auto& data : data_v )
        std::cout << "  " << data << std::endl;
      matchdata_vv.emplace_back( std::move(data_v) );
    }

    /// initial screen for overlapping times
    // start from Y plane, match U, then V

    std::vector<MaskCombo> combo_v;
    
    for ( size_t iy=0; iy<matchdata_vv.at(2).size(); iy++ ) {
      auto const& data_yplane = matchdata_vv.at(2).at(iy);
      auto const& mask_yplane = clustermask_vv.at(2).at( data_yplane.index );
      std::vector<int> indices(3);
      indices[2] = iy;

      for ( size_t iu=0; iu<matchdata_vv.at(0).size(); iu++ ) {

        auto const& data_uplane = matchdata_vv.at(0).at(iu);
        auto const& mask_uplane = clustermask_vv.at(0).at(data_uplane.index);


        if ( data_uplane.tick_min > data_yplane.tick_max ) {
          break; // since we sorted, we can cut off
        }

        // create a combo filled with Y-plane to start
        MaskCombo combo;
        combo.addMask( mask_yplane, data_yplane );
        

        if ( combo.iscompatible( data_uplane ) ) {
          // add u-plane
          combo.addMask( mask_uplane, data_uplane );
          //std::cout << "y-u match: " << combo << std::endl;
          // store
          combo_v.emplace_back( std::move(combo) );
        }
      }//end of u-plane loop

    }
    std::sort( combo_v.begin(), combo_v.end() );
    std::cout << "combos after Y-U match: " << combo_v.size() << std::endl;    
    //for ( auto const& combo : combo_v )
    //  std::cout << "  " << combo << std::endl;

    // pair with v-plane
    m_combo_3plane_v.clear();
    
    for ( auto const& combo : combo_v ) {

      for ( size_t iv=0; iv<matchdata_vv.at(1).size(); iv++ ) {
        auto const& data_vplane = matchdata_vv.at(1).at(iv);
        auto const& mask_vplane = clustermask_vv.at(1).at(data_vplane.index);

        if ( combo.union_tick[1] < data_vplane.tick_min ) break;

        if ( combo.iscompatible( data_vplane ) ) {

          MaskCombo combo_yuv( combo );

          combo_yuv.addMask( mask_vplane, data_vplane );

          std::cout << "combo Y-U-V: " << combo_yuv << std::endl;
          if ( combo_yuv.iou()>0.5 )
            m_combo_3plane_v.emplace_back( std::move(combo_yuv) );
          
        }
      }// loop over v-plane cluster/mask

    }// end of combo y-u loop

    std::sort( m_combo_3plane_v.begin(), m_combo_3plane_v.end() );
    //std::cout << "combos after Y-U-V match: " << m_combo_3plane_v.size() << std::endl;
    //for ( auto const& combo : m_combo_3plane_v )
    //std::cout << "  " << combo << std::endl;

    // make crops, mask charge, contours on charge, contours on mask
    for ( size_t icombo=0; icombo<m_combo_3plane_v.size(); icombo++ ) {
      CropMaskCombo     cropmaker( m_combo_3plane_v.at(icombo), wholeview_v );
      FeaturesMaskCombo features( cropmaker );
      Gen3DEndpoints    endpoints( features );

      // run astar only if the 3d points are fairly consistent
      bool runastar = true;
      for ( auto const& triarea_score : endpoints.endpt_tri_v )
        if ( triarea_score>200 ) runastar = false;
      
      AStarMaskCombo    astar( endpoints, badch_v, runastar );
      
      m_combo_crops_v.emplace_back( std::move(cropmaker) );
      m_combo_features_v.emplace_back( std::move(features) );
      m_combo_endpt3d_v.emplace_back( std::move(endpoints) );
      m_combo_astar_v.emplace_back( std::move(astar) );

      // for debug
      //if ( icombo==0 )
      //  break;
    }

    // initial round of selection

    // PASS 2 using only 2 plane matches and remaining masks
        
  }

}
}
