#ifndef __FLASH_MATCHER_H__
#define __FLASH_MATCHER_H__

/*
 * class for flash matching with mc truth
 */

// larlite/core
#include "larlite/DataFormat/storage_manager.h"

namespace ublarcvapp {
namespace mctools {

  class FlashMatcher {
  public:

    FlashMatcher() {};
    virtual ~FlashMatcher() {};

    static int numTracks( larlite::storage_manager& ioll );
    static std::tuple<double, std::string, Bool_t> grabTickFromMCTrack( larlite::storage_manager& ioll, int i );
    static std::vector<double> grabTickFromOpflash( larlite::storage_manager& opio, std::string producer );
    static double matchTicks( double mctrack_tick, std::vector<double> flash_ticks, Bool_t isCosmic );

  };

}
}

#endif
