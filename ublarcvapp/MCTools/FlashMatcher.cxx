#include "FlashMatcher.h"

#include "larlite/LArUtil/Geometry.h"
#include "larlite/DataFormat/mctruth.h"
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

  Float_t FlashMatcher::grabTickFromMCTrack( larlite::storage_manager& ioll )
  {
    larlite::event_mctruth* ev_mctruth
      = (larlite::event_mctruth*)ioll.get_data(larlite::data::kMCTruth,"generator");

    auto const& mctruth = ev_mctruth->at(0);
    const larlite::mcstep& start = mctruth.GetNeutrino().Nu().Trajectory().front();

    std::cout << "Time: " << start.T() << std::endl;

    larutil::SpaceChargeMicroBooNE* _sce = nullptr;
    Float_t tick = CrossingPointsAnaMethods::getTick(start, 4050.0, _sce);

    return tick;

  }

  std::vector<double> FlashMatcher::grabTickFromOpflash( larlite::storage_manager& opio )
  {
  std::vector<double> tickContainer = {};

  larlite::event_opflash* ev_opflash_cosmic
    = (larlite::event_opflash*)opio.get_data(larlite::data::kOpFlash,"simpleFlashCosmic");

  for (auto const& opflash : *ev_opflash_cosmic) {
  double time = opflash.Time();
  std::cout << time << std::endl;
  double tick = time*1.0e-3 + 3200.0;
  tickContainer.push_back(tick);
  }

  //larutil::SpaceChargeMicroBooNE* _sce = nullptr;
  //tick = CrossingPointsAnaMethods::getTick(start, 4050.0, _sce);

  return tickContainer;
  }

Float_t FlashMatcher::matchTicksFromTrackAndFlash( larlite::storage_manager& opio )
{


return 0.0;

}

}
}
