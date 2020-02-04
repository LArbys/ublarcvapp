#ifndef ASTARUTILS_H
#define ASTARUTILS_H

// larcv
#include "larcv/core/DataFormat/Image2D.h"

// ROOT
#include "TVector3.h"

// opencv
#include <opencv2/core.hpp>

namespace ublarcvapp {
namespace reco3d {

  void
  ProjectTo3D(const larcv::ImageMeta& meta,
	      double parent_x,
	      double parent_y,
	      double parent_z,
	      double parent_t,
	      uint plane,
	      double& xpixel, double& ypixel);

#ifdef LARCV_OPENCV
  std::vector<std::vector<cv::Point_<int> > > TrackToPixels(const std::vector<TVector3>& xyz_v,
							      const std::vector<larcv::ImageMeta>& meta_v);
#endif

}
}

#endif
/** @} */ // end of doxygen group
