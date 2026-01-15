#ifndef __UBLARCAPP_UBIMAGEMOD_POINT_IMAGE_PROJECTION_H__
#define __UBLARCAPP_UBIMAGEMOD_POINT_IMAGE_PROJECTION_H__

/**
 * 
 * @brief Class with functions for projecting a point into a larcv image
 *
 * Image represents the output of a lartpc wireplane.  
 * We want to be able to take a 3D point inside the TPC and find the pixel that location
 * projects to.  Then we can ask about the content of the pixels in the neighborhood of that pixel.
 *
 *
 */

#include <vector>

#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/Image2D.h"


namespace ublarcvapp {
  namespace ubimagemod {
    
    class PointImageProjection : public larcv::larcv_base
    {

    public:

      PointImageProjection()
	: larcv::larcv_base("PointImageProjection")
	{};

      ~PointImageProjection() {};


      float getPixelSumAroundProjPoint( const std::vector<float>& xyz, const larcv::Image2D& img, int pixel_kernel_radius, float pix_threshold ) const;
	
      
    };
    
  }
}

#endif
