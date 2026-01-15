#ifndef __UBLARCVAPP_UBIMAGEMOD_TRACK_IMAGE_MASK_H__
#define __UBLARCVAPP_UBIMAGEMOD_TRACK_IMAGE_MASK_H__

#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "larlite/DataFormat/track.h"

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
    ~TrackImageMask()
    {
      clear();
    }

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
                        const float threshold,
                        const float maxstepsize=0.1 );

    int makePixelList( const larlite::track& track,
                       const larcv::Image2D& img,
                       const float threshold,
                       const float maxstepsize=0.1,                       
                       bool fill_map=true );


    float sumOverPixelList( const std::vector< std::vector<int> >& pixlist_v,
                            const larcv::Image2D& img,
                            const int dcol,
                            const int drow,
                            const float thresh ) const;
    
    bool getExtremaPixelValues( const std::vector< std::vector<int> >& pixlist_v,
                                const larcv::Image2D& img,
                                const int dcol,
                                const int drow,
                                const float thresh,
                                float& pixval_min,
                                float& pixval_max ) const;

    float getMaximumChargeGap( std::vector< std::vector<float> >& gap_points ) const;
    int getMaskedPixels() const { return nmaskedpixels; };
    int getPathPixels() const { return npathpixels; };
    
    
    // first we make a list of pixels covered by the track
    // we also get the min and max col bounds
    struct Pix_t {
      int col;
      int row;
      float smin;
      float smax;
      float pixval;
      Pix_t()
        : col(0),row(0),smin(0),smax(0),pixval(0)
      {};
      Pix_t( int c, int r, float ssmin, float ssmax,float pv )
        : col(c), row(r), smin(ssmin), smax(ssmax),pixval(pv)
      {};
    };

    // output data members for use
    std::map< std::pair<int,int>, Pix_t > pixel_map; ///< map of pixel (col,row) to Pix_t struct
    std::vector< std::vector<int> > pixel_v; ///< list of (col,row) pixels that follows path of track
    std::vector< float > pixel_kernel_pixval_v; ///< list of charge values within the kernel for each pixel
    int min_col;
    int min_row;
    int max_col;
    int max_row;
    int kernel_dcol;
    int kernel_drow;
    int kernel_N;
    int nmaskedpixels;
    int npathpixels;

    bool show_timing; ///< set to true to print output
  };
  
}
}


#endif
