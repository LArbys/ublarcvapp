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
    
    
  }

}
}
                                          
    

