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
namespace mctools {

  class CrossingPointsAnaMethods {

  public:

    CrossingPointsAnaMethods(){};
    virtual ~CrossingPointsAnaMethods(){};    
    
    
    static float getTick( const std::vector<float>& step, const float trig_time=4050.0,
                          const larutil::SpaceChargeMicroBooNE* psce=nullptr );
  
    static float getTick( const larlite::mcstep& step, const float trig_time=4050.0,
                          const larutil::SpaceChargeMicroBooNE* psce=nullptr );

    static int doesTrackCrossImageBoundary( const larlite::mctrack& track, const larcv::ImageMeta& meta,
                                            const float trig_time, const larutil::SpaceChargeMicroBooNE* psce=nullptr );

    static std::vector<int> getFirstStepPosInsideImage( const larlite::mctrack& track, const larcv::ImageMeta& meta, const float trig_time,
                                                        const bool startAtstart, const float max_step_size, const float fv_border,
                                                        const larutil::SpaceChargeMicroBooNE* psce );

  };
  
}
}

#endif
