#ifndef __UBLARCVAPP_MCTOOLS_TRUTHTRACK_SCE_H__
#define __UBLARCVAPP_MCTOOLS_TRUTHTRACK_SCE_H__

#include <vector>
#include "LArUtil/SpaceChargeMicroBooNE.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/track.h"
#include "larcv/core/Base/larcv_base.h"

namespace ublarcvapp {
namespace mctools {

  class TruthTrackSCE : public larcv::larcv_base {

  public:

    TruthTrackSCE();
    TruthTrackSCE(larutil::SpaceChargeMicroBooNE* psce);
    virtual ~TruthTrackSCE();

    larlite::track applySCE( const larlite::mctrack& mct );

    void dist2track( const std::vector<float>& testpt,
                     const larlite::track& track,
                     float& min_r,
                     int& min_step );
    
    float pointLineDistance( const std::vector<float>& linept1,
                             const std::vector<float>& linept2,
                             const std::vector<float>& pt );
    
  protected:

    bool _kown_sce;
    const larutil::SpaceChargeMicroBooNE* _p_sce;
    bool _kdebug;

  };
  
}
}


#endif
