#ifndef __UBLARCVAPP_MCTOOLS_TRUTHSHOWERTRUNKSCE_H__
#define __UBLARCVAPP_MCTOOLS_TRUTHSHOWERTRUNKSCE_H__

#include <vector>
#include "LArUtil/SpaceChargeMicroBooNE.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/track.h"
#include "larcv/core/Base/larcv_base.h"

namespace ublarcvapp {
namespace mctools {

  /**
   * @ingroup MCTools
   * @class TruthShowerTrunkSCE
   *
   * Provides a function to take a larlite::mcshower instance
   * and provides a shower trunk that takes into the space-charge effect when relevant.
   * This will provide a 3D line segment that matches what the wire-planes in a LArTPC will see.
   *
   */
  class TruthShowerTrunkSCE : public larcv::larcv_base {

  public:

    TruthShowerTrunkSCE();
    TruthShowerTrunkSCE(larutil::SpaceChargeMicroBooNE* psce);
    virtual ~TruthShowerTrunkSCE();

    larlite::track applySCE( const larlite::mcshower& mcs );
    
  protected:

    bool _kown_sce; ///< if true, class assumes it owns memory pointed to by _p_sce and will destroy it in destructor
    const larutil::SpaceChargeMicroBooNE* _p_sce; ///< pointer to class that calculates charge displacement due to space-charge effect
    bool _kdebug; ///< if true, dump verbose output

  };
  
  
}
}
  

#endif
