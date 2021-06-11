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
      min_col(0),
      min_row(0),
      max_col(0),
      max_row(0),
      show_timing(false)
      {};

    void clear();
    
    int maskTrack( const larlite::track& track,
                   const larcv::Image2D& adc,
                   larcv::Image2D& mask,
                   const float thresh,                                    
                   const int dcol,
                   const int drow,
                   const float maxstepsize );

    int labelTrackPath( const larlite::track& track,
                        const larcv::Image2D& img,                         
                        larcv::Image2D& smin,
                        larcv::Image2D& smax,
                        const float maxstepsize=0.1 );

    int makePixelList( const larlite::track& track,
                       const larcv::Image2D& img,
                       const float maxstepsize=0.1,                       
                       bool fill_map=true );
    

    // first we make a list of pixels covered by the track
    // we also get the min and max col bounds
    struct Pix_t {
      int col;
      int row;
      float smin;
      float smax;
      Pix_t()
        : col(0),row(0),smin(0),smax(0)
      {};
      Pix_t( int c, int r, float ssmin, float ssmax )
        : col(c), row(r), smin(ssmin), smax(ssmax)
      {};
    };

    // output data members for use
    std::map< std::pair<int,int>, Pix_t > pixel_map; ///< map of pixel (col,row) to Pix_t struct
    std::vector< std::vector<int> > pixel_v; ///< list of (col,row) pixels that follows path of tracko
    int min_col;
    int min_row;
    int max_col;
    int max_row;    

    bool show_timing; ///< set to true to print output
  };
  
}
}


#endif
