#include "FlashMatcher.h"

#include <algorithm>

#include "larlite/LArUtil/Geometry.h"
#include "larlite/DataFormat/mctrack.h"
#include "larlite/DataFormat/opflash.h"

#include "crossingPointsAnaMethods.h"

namespace ublarcvapp {
namespace mctools {

  /*
   * grab time coordinate -> tick from mctrack mcstep
   *
   * @param[in] ioll The larlite storage manager that contains mctruth class
   * @return tick
   */

  std::tuple<double, std::string, Bool_t> FlashMatcher::grabTickFromMCTrack( larlite::storage_manager& ioll )
  {
    larlite::event_mctrack* ev_mctrack
      = (larlite::event_mctrack*)ioll.get_data(larlite::data::kMCTrack,"mcreco");

    auto const& mctrack = ev_mctrack->at(0);

    std::cout << "Origin: " << mctrack.Origin() << std::endl;

    std::string producer;
    Bool_t isCosmic;
    if ( mctrack.Origin() == 1 ) {
      producer = "simpleFlashBeam";
    } else {
      producer = "simpleFlashCosmic";
      isCosmic = 1;
    }

    const larlite::mcstep& start = mctrack.Start();

    //std::cout << "Time: " << start.T() << std::endl;

    larutil::SpaceChargeMicroBooNE* _sce = nullptr;
    double tick = CrossingPointsAnaMethods::getTick(start, 4050.0, _sce);

    return std::make_tuple( tick, producer, isCosmic );

  }

  std::vector<double> FlashMatcher::grabTickFromOpflash( larlite::storage_manager& opio, std::string producer ) {

    //std::cout << "If this is a neutrino, this would print 1: " << isNeutrino << std::endl;

    std::vector<double> tickContainer = {};

    larlite::event_opflash* ev_opflash
      = (larlite::event_opflash*)opio.get_data(larlite::data::kOpFlash, producer);

    for (auto const& opflash : *ev_opflash) {
      double time = opflash.Time();
      //std::cout << time << std::endl;
      double tick = time/0.5 + 3200.0;
      tickContainer.push_back(tick);
    }

  std::sort( tickContainer.begin(), tickContainer.end() );

  return tickContainer;
  }

  double FlashMatcher::matchTicks( double mctrackTick, std::vector<double> flashTicks, Bool_t isCosmic ) {

    // for cosmic tracks, will return -999.999 if there is no opflash found within the threshold

    double threshold;
    if (isCosmic == 1) {
        threshold = 400.0; // 1 tick = 0.5 us
    } else {
      threshold = std::numeric_limits<double>::infinity();
    }

    auto match = std::lower_bound( flashTicks.begin(), flashTicks.end(), mctrackTick );

    double a = *(match - 1);
    double b = *(match);

    if (match == flashTicks.begin() && fabs(b-mctrackTick) <= threshold) {
      return flashTicks[0];
    }

    if (fabs(mctrackTick - a) < fabs(mctrackTick - b) && fabs(mctrackTick - a) <= threshold) {
      return flashTicks [ match - flashTicks.begin() - 1 ];
    }

    if ( fabs(mctrackTick - b) <= threshold  )
      return flashTicks[ match - flashTicks.begin() ];

    return -999.999;

  }


}
}
