#ifndef __UBLARCVAPP_RECO3D_TRACK_REVERSER_H__
#define __UBLARCVAPP_RECO3D_TRACK_REVERSER_H__

#include "larcv/core/Base/larcv_base.h"
#include "larlite/DataFormat/track.h"

namespace ublarcvapp {
namespace reco3d {

  class TrackReverser : public larcv::larcv_base {
  public:

    TrackReverser()
      : larcv::larcv_base("TrackReverser")
      {};
    virtual ~TrackReverser() {};
    
    static larlite::track reverseTrack( const larlite::track& original,
                                        bool flip_direction=false );
    
  };
  
}
}

#endif

