#ifndef __PATH_2_PIXELS_H__
#define __PATH_2_PIXELS_H__

#include <vector>

// larlite
#include "DataFormat/track.h"

// larcv
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/Pixel2DCluster.h"

namespace ublarcvapp {
namespace tagger {

  // w/ badch
  std::vector<larcv::Pixel2DCluster> getTrackPixelsFromImages( const std::vector< std::vector<double> >& path3d,
							       const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchimgs,
							       const std::vector<float>& thresholds, const std::vector<int>& neighborhood_size,
							       const float stepsize );
  // w/o badch
  std::vector<larcv::Pixel2DCluster> getTrackPixelsFromImagesNoBadCh( const std::vector< std::vector<double> >& path3d,
								      const std::vector<larcv::Image2D>& imgs, 
								      const std::vector<float>& thresholds, const std::vector<int>& neighborhood_size,
								      const float stepsize );

  // w/ badch
  std::vector<larcv::Pixel2DCluster> getTrackPixelsFromImages( const larlite::track& lltrack,
							       const std::vector<larcv::Image2D>& imgs, const std::vector<larcv::Image2D>& badchimgs,
							       const std::vector<float>& thresholds, const std::vector<int>& neighborhood_size,
							       const float stepsize );
  

}
}

#endif
