#ifndef __UBLARCVAPP_UBIMAGEMOD_TRACK_IMAGE_MASK_H__
#define __UBLARCVAPP_UBIMAGEMOD_TRACK_IMAGE_MASK_H__

#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "DataFormat/track.h"

namespace ublarcvapp {
namespace ubimagemod {

  class TrackImageMask : public larcv::larcv_base {

  public:
    
    TrackImageMask()
      : larcv::larcv_base( "TrackImageMask" ),
      show_timing(false)
      {};

    int maskTrack( const larlite::track& track,
                   const larcv::Image2D& adc,
                   larcv::Image2D& mask,
                   const float thresh,                                    
                   const int dcol,
                   const int drow,
                   const float maxstepsize );
    
    bool show_timing;
  };
  
}
}


#endif
