#ifndef __UBLARCVAPP_MCTOOLS_CONVERTMCINFOFORLARCV2_H__
#define __UBLARCVAPP_MCTOOLS_CONVERTMCINFOFORLARCV2_H__

#include <vector>
#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/Particle.h"
#include "larlite/DataFormat/mctrack.h"
#include "larlite/DataFormat/mcshower.h"
#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"

namespace ublarcvapp {
namespace mctools {

  /**
   * @class ConvertMCInfoForLArCV2
   * @brief Convert larlite particle info into larcv::Particle objects
   *
   * This is mostly targeted for lartpc_mlreco3d use.
   * we serve up the info needed to be able to run
   *  lartpc_mlreco3d/mlreco/iotools/parsers/particles.py:parse_particle_points
   *
   */
  class ConvertMCInfoForLArCV2 : public larcv::larcv_base {
    
  public:
    ConvertMCInfoForLArCV2()
      : larcv::larcv_base("ConvertMCInfoForLArCV2"),
	kApplySCE(false),
	kApplyT0drift(false),
	_psce(nullptr)
    {};
    virtual ~ConvertMCInfoForLArCV2();

    std::vector<larcv::Particle> convert( const std::vector<larlite::mctrack>& mctrack,
					  const std::vector<larlite::mcshower>& mcshower  );

    bool kApplySCE; ///< if true, we apply Space charge corrections
    bool kApplyT0drift; ///< if true, we use the time of the event releative to the trigger to shift the position of the event.
    void doApplySCE( bool doit ) { kApplySCE=doit; };
    void doApplyT0drift( bool doit ) { kApplyT0drift=doit; };
    std::vector<double> preparePosition( double x, double y, double z, double t );

    larutil::SpaceChargeMicroBooNE* _psce;
  };

}
}

#endif
