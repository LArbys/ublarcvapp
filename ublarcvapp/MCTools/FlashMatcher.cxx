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

  double FlashMatcher::grabTickFromMCTrack( larlite::storage_manager& ioll )
  {
    larlite::event_mctrack* ev_mctrack
      = (larlite::event_mctrack*)ioll.get_data(larlite::data::kMCTrack,"mcreco");

    auto const& mctrack = ev_mctrack->at(0);
    const larlite::mcstep& start = mctrack.Start();

    std::cout << "Time: " << start.T() << std::endl;

    larutil::SpaceChargeMicroBooNE* _sce = nullptr;
    double tick = CrossingPointsAnaMethods::getTick(start, 4050.0, _sce);

    return tick;

  }

  std::vector<double> FlashMatcher::grabTickFromOpflash( larlite::storage_manager& opio , std::string producer) {

    std::vector<double> tickContainer = {};

    larlite::event_opflash* ev_opflash
      = (larlite::event_opflash*)opio.get_data(larlite::data::kOpFlash, producer);

    for (auto const& opflash : *ev_opflash) {
      double time = opflash.Time();
      std::cout << time << std::endl;
      double tick = time/0.5 + 3200.0;
      tickContainer.push_back(tick);
    }

  std::sort( tickContainer.begin(), tickContainer.end() );

  return tickContainer;
  }

double FlashMatcher::matchTicks( double mctrackTick, std::vector<double> flashTicks ) {

  auto match = std::lower_bound( flashTicks.begin(), flashTicks.end(), mctrackTick );

  if (match == flashTicks.begin()) {
    return flashTicks[0];
  }

  double a = *(match - 1);
  double b = *(match);

  if (fabs(mctrackTick - a) < fabs(mctrackTick - b)) {
        return flashTicks [ match - flashTicks.begin() - 1 ];
    }

  return flashTicks[ match - flashTicks.begin() ];

}

}
}
