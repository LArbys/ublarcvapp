#ifndef __crossingPointsAnaMethods_h__
#define __crossingPointsAnaMethods_h__

#include <vector>

// larlite
#include "DataFormat/mctrack.h"
#include "DataFormat/mctrajectory.h"

class TTree;

namespace larcv {
  class ImageMeta;
}

namespace larutil {
  class SpaceChargeMicroBooNE;
}

namespace ublarcvapp {

  float getTick( const std::vector<float>& step, const float trig_time=4050.0,
                 const larutil::SpaceChargeMicroBooNE* psce=NULL );
  
  float getTick( const larlite::mcstep& step, const float trig_time=4050.0,
                 const larutil::SpaceChargeMicroBooNE* psce=NULL );

  int doesTrackCrossImageBoundary( const larlite::mctrack& track, const larcv::ImageMeta& meta,
                                   const float trig_time, const larutil::SpaceChargeMicroBooNE* psce );

  std::vector<int> getImageBoundaryCrossingPoint( const larlite::mctrack& track, std::vector<float>& crossingpt, const larcv::ImageMeta& meta,
						  const float boundary_tick_buffer, const float trig_time,
                                                  const larutil::SpaceChargeMicroBooNE* psce );

  std::vector<float> getFirstStepPosInsideImage( const larlite::mctrack& track, const larcv::ImageMeta& meta, const float trig_time,
						 const bool startAtstart, const float max_step_size, const float fv_border,
                                                 const larutil::SpaceChargeMicroBooNE* psce );
  
}

#endif
