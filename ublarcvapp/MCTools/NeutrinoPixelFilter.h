#ifndef __UBLARCVAPP_MCTOOLS_NEUTRINOPIXELFILTER_H__
#define __UBLARCVAPP_MCTOOLS_NEUTRINOPIXELFILTER_H__

#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/IOManager.h"

namespace ublarcvapp {
namespace mctools {

  /**
   * @class NeutrinoPixelFilter
   * @ingroup MCTools
   * @brief Class isolates neutrino pixels in an image using instance+ancestor images
   */  
  class NeutrinoPixelFilter : public larcv::larcv_base {
  public:

    NeutrinoPixelFilter()
      {};

    ~NeutrinoPixelFilter()
      {};

    void process( larcv::IOManager& iolarcv );
    
  };

}
}

#endif
