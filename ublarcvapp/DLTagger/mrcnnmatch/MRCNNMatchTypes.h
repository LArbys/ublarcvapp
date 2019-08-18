#ifndef __MRCNN_MATCH_TYPES_H__
#define __MRCNN_MATCH_TYPES_H__

#include "larcv/core/DataFormat/ClusterMask.h"

namespace ublarcvapp {
namespace dltagger {

  class MaskMatchData {
  public:

    MaskMatchData( int plane, int index, const larcv::ClusterMask& mask );
    virtual ~MaskMatchData() {};

    int plane;
    int index;
    float tick_min;
    float tick_max;
    float detz_min;
    float detz_max;
      
  };
}
}

#endif
