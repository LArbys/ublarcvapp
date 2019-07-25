#ifndef __PIXEL_2_SPACE_POINT_H__
#define __PIXEL_2_SPACE_POINT_H__

#include <vector>

// larcv
#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/DataFormat/Pixel2D.h"

// tagger
#include "BoundaryMuonTaggerTypes.h"
#include "BoundarySpacePoint.h"

namespace ublarcvapp {
namespace tagger {

  BoundarySpacePoint Pixel2SpacePoint( const std::vector<larcv::Pixel2D>& pixels,
                                       const BoundaryEnd_t endtype,
                                       const larcv::ImageMeta& meta );

}
}

#endif
