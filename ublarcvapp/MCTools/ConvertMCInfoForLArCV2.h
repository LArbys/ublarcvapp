#ifndef __UBLARCVAPP_MCTOOLS_CONVERTMCINFOFORLARCV2_H__
#define __UBLARCVAPP_MCTOOLS_CONVERTMCINFOFORLARCV2_H__

#include <vector>
#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/Particle.h"
#include "larlite/DataFormat/mctrack.h"
#include "larlite/DataFormat/mcshower.h"

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
    : larcv::larcv_base("ConvertMCInfoForLArCV2")
    {};
    ~ConvertMCInfoForLArCV2() {};

    std::vector<larcv::Particle> convert( const std::vector<larlite::mctrack>& mctrack,
					  const std::vector<larlite::mcshower>& mcshower  );
    
  };

}
}

#endif
