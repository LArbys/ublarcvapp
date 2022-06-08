#ifndef __crossingPointsAnaMethods_h__
#define __crossingPointsAnaMethods_h__

#include <vector>

// larlite
#include "larlite/DataFormat/mctrack.h"
#include "larlite/DataFormat/mctrajectory.h"

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
    
    
    static float getTick( const std::vector<float>& step,
			  const int tpcid, const int cryoid, const float trig_time=4050.0,
                          const larutil::SpaceChargeMicroBooNE* psce=nullptr );
  
    static float getTick( const larlite::mcstep& step,
			  const int tpcid, const int cryoid, const float trig_time=4050.0,
                          const larutil::SpaceChargeMicroBooNE* psce=nullptr );

    static int doesTrackCrossImageBoundary( const larlite::mctrack& track, const larcv::ImageMeta& meta,
                                            const float trig_time, const larutil::SpaceChargeMicroBooNE* psce=nullptr );

    static std::vector<int> getFirstStepPosInsideImage( const larlite::mctrack& track,
                                                        const larcv::ImageMeta& meta,
                                                        const float trig_time,
                                                        const bool startAtstart,
                                                        const float max_step_size,
                                                        const float fv_border,
                                                        std::vector<float>& endpt,
                                                        const larutil::SpaceChargeMicroBooNE* psce,
                                                        bool verbose );

  };
  
}
}

#endif
