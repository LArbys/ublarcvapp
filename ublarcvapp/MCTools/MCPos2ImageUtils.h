#ifndef __UBLARCVAPP_MCTOOLS_MCPOS2IMAGEUTILS_H__
#define __UBLARCVAPP_MCTOOLS_MCPOS2IMAGEUTILS_H__

#include <vector>

#include "larcv/core/Base/larcv_base.h"

#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"
#include "larlite/DataFormat/mctrack.h"

namespace ublarcvapp {
namespace mctools {

  /**
   * The purpose of this singleton class is to standardize common conversions between
   *   the detector geometry and image locations.
   * This requires detector-specific information that includes the geometry, drift direction, 
   *   electronics clock, field, and space-charge effect.
   * This utility class specializes in the conversions that involve MC values and objects.
   *
   * This is an attempt to standardize tools for later use across SBN.
   *
   */

  class MCPos2ImageUtils : public larcv::larcv_base {
    
  public:

    // default constructor
    MCPos2ImageUtils()
      : larcv::larcv_base("MCPos2ImageUtils"),
      psce(nullptr)
      {
	psce = new larutil::SpaceChargeMicroBooNE;
      };


    // default destructor
    virtual ~MCPos2ImageUtils() {
      delete psce;
    };

    static MCPos2ImageUtils* Get() {
      if ( !_p_singleton ) {
	_p_singleton = new MCPos2ImageUtils;
      }
      return _p_singleton;
    };
    

    std::vector< std::vector<float> >
    getRecoSpacepoints( const larlite::mctrack& track,
			bool apply_t0_shift=true,
			bool apply_sce=true );

    std::vector<float>
    truepos_to_recopos( const float x_cm, const float y_cm, const float z_cm,
			const float t_ns,
			bool apply_t0_shift=true,
			bool apply_sce=true );

    std::vector<float>
    get_sce_shifted_pos( const float x, const float y, const float z );


    larutil::SpaceChargeMicroBooNE* psce;

  private:

    static MCPos2ImageUtils* _p_singleton;
    
  };
  
}
}

#endif
