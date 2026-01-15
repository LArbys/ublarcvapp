#ifndef __UBLARCVAPP_MCTOOLS_TRUTHTRACK_SCE_H__
#define __UBLARCVAPP_MCTOOLS_TRUTHTRACK_SCE_H__

#include <vector>
#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"
#include "larlite/DataFormat/mctrack.h"
#include "larlite/DataFormat/track.h"
#include "larcv/core/Base/larcv_base.h"

namespace ublarcvapp {
namespace mctools {

  /**
   * @ingroup MCTools
   * @class TruthTrackSCE
   *
   * Provides a function to take the true track trajectory (stored as a larlite::mctrack object)
   * and applies the space-charge effect. This will provide a 3D 
   * trajectory that matches what the wire-planes in a LArTPC will see.
   *
   */
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

    bool _kown_sce; ///< if true, class assumes it owns memory pointed to by _p_sce and will destroy it in destructor
    const larutil::SpaceChargeMicroBooNE* _p_sce; ///< pointer to class that calculates charge displacement due to space-charge effect
    bool _kdebug; ///< if true, dump verbose output

  };
  
}
}


#endif
