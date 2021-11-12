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

    FlashMatcher() {

        isCosmic = 0;

    };
    virtual ~FlashMatcher() {};

    int numTracks( larlite::storage_manager& ioll );
    int numShowers( larlite::storage_manager& ioll );
    double grabTickFromMCTrack( larlite::storage_manager& ioll, int i );
    double grabTickFromMCShower( larlite::storage_manager& ioll, int i );
    std::vector<double> grabTickFromOpflash( larlite::storage_manager& opio );
    double matchTicks( double mctrack_tick, std::vector<double> flash_ticks );

    Bool_t isCosmic;
    std::string producer;

  };

}
}

#endif
