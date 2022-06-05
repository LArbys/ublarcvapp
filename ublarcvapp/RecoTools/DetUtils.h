#ifndef __UBLARCVAPP_RECOUTILS_DETUTILS_H__
#define __UBLARCVAPP_RECOUTILS_DETUTILS_H__

#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "larlite/LArUtil/Geometry.h"

namespace ublarcvapp {
namespace recotools {

  class DetUtils : public larcv::larcv_base {
    
  public:
    DetUtils() {};
    virtual ~DetUtils() {};
      
    static std::vector< const larcv::Image2D* > getTPCImages( const std::vector<larcv::Image2D>& adc_v,
							      const int tpcid, const int cryoid );
    
  };

}
}

#endif
