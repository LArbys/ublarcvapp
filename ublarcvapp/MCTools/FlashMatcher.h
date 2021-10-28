#ifndef __FLASH_MATCHER_H__
#define __FLASH_MATCHER_H__

/*
 * class for flash matching with mc truth
 */

#include <vector>

// larlite/core
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"

namespace ublarcvapp {
namespace mctools {

  class FlashMatcher {
  public:

    FlashMatcher() {};
    virtual ~FlashMatcher() {};

    std::string   matchTrackAndFlash( larlite::storage_manager& ioll );

  };

}
}

#endif
