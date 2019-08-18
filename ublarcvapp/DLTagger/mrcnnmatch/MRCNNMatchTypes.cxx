#include "MRCNNMatchTypes.h"

namespace ublarcvapp {
namespace dltagger {

  MaskMatchData::MaskMatchData( int _plane, int _index, const larcv::ClusterMask& mask )
    : plane(_plane),
      index(_index)
  {

    // we loop over the points in the mask
    std::cout << "meta: " << mask.meta.dump() << std::endl;
    int xoffset = mask.box.min_x();
    int yoffset = mask.box.min_y();
    // for ( size_t ipt=0; ipt<mask.points_v.size(); ipt++ ) {
      
    // }

  }
  
}
}
    
