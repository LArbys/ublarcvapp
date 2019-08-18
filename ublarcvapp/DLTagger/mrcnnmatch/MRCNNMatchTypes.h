#ifndef __MRCNN_MATCH_TYPES_H__
#define __MRCNN_MATCH_TYPES_H__

#include <iostream>
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
    float wire_min;
    float wire_max;
    float detz_min;
    float detz_max;

    // comparison operator needed for sorting
    bool operator < (const MaskMatchData& b) const {
      if ( tick_min<b.tick_min ) return true;
      else if ( tick_min==b.tick_min ) {
        if ( detz_min<b.detz_min ) return true;
      }
      return false;
    };

    // stdout streamer
    friend std::ostream& operator<<(std::ostream &os,const MaskMatchData& m);
    
  };


  
}
}

#endif
