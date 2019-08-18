#include "MRCNNMatch.h"

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

    
    
    
    
  }

}
}
                                          
    

