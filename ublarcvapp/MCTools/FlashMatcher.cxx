#include "FlashMatcher.h"

#include "larlite/DataFormat/mctrack.h"
#include "larlite/DataFormat/opflash.h"

#include "crossingPointsAnaMethods.h"

namespace ublarcvapp {
namespace mctools {

  /*
   * grab time coordinate from mctrack mcstep -> convert to tick
   *
   * @param[in] ioll The larlite storage manager that contains mctruth class
   * @return tuple with tick, producer string, and cosmic flag
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

    larutil::SpaceChargeMicroBooNE* _sce = nullptr;
    double tick = CrossingPointsAnaMethods::getTick(start, 4050.0, _sce);

    return std::make_tuple( tick, producer, isCosmic );

  }

  /*
   * grab time coordinate from opflash -> convert to tick
   *
   * @param[in] opio The larlite storage manager that contains opflash class
   * @param[in] producer String labeling the producer, e.g. "simpleFlashBeam"
   * @return Vector containing all opflash times for the track in ticks
   */
  std::vector<double> FlashMatcher::grabTickFromOpflash( larlite::storage_manager& opio, std::string producer ) {

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

  /*
   * matches mctrack tick to opflash tick
   * for cosmics, will return -999.999 if there is no opflash found within the threshold
   * for beam tracks, will find the closest matching opflash
   *
   * @param[in] mctrackTick Time in ticks for the mctrack to be matched to
   * @param[in] flashTicks Vector of potential opflash matches in ticks
   * @param[in] isCosmic Cosmic flag to determine threshold
   *
   * @return Value of the closest matching opflash tick
   */
  double FlashMatcher::matchTicks( double mctrackTick, std::vector<double> flashTicks, Bool_t isCosmic ) {

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
