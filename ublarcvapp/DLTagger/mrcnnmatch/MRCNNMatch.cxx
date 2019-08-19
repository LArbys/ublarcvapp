#include "MRCNNMatch.h"

#include <algorithm>

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
                                           std::vector< std::vector<int> >& match_indices ) {

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
    std::cout << "combos after Y-U match: " << combo_v.size() << std::endl;
    std::sort( combo_v.begin(), combo_v.end() );
    for ( auto const& combo : combo_v )
      std::cout << "  " << combo << std::endl;
  }

}
}
                                          
    

