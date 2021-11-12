#include "FlashMatcher.h"

#include "larlite/DataFormat/mcshower.h"
#include "larlite/DataFormat/mctrack.h"
#include "larlite/DataFormat/opflash.h"

#include "crossingPointsAnaMethods.h"

namespace ublarcvapp {
namespace mctools {


  int FlashMatcher::numTracks( larlite::storage_manager& ioll ) {
    larlite::event_mctrack* ev_mctrack
      = (larlite::event_mctrack*)ioll.get_data(larlite::data::kMCTrack,"mcreco");

    int numTracks = ev_mctrack->size();

    return numTracks;
  }

  int FlashMatcher::numShowers( larlite::storage_manager& ioll ) {

    larlite::event_mcshower* ev_mcshower
      = (larlite::event_mcshower*)ioll.get_data(larlite::data::kMCShower,"mcreco");

    int numShowers = ev_mcshower->size();

    return numShowers;

  }
  /*
   * grab time coordinate from mctrack mcstep -> convert to tick
   *
   * @param[in] ioll The larlite storage manager that contains mctruth class
   * @return tuple with tick, producer string, and cosmic flag
   */
  double FlashMatcher::grabTickFromMCTrack( larlite::storage_manager& ioll, int i ) {
    larlite::event_mctrack* ev_mctrack
      = (larlite::event_mctrack*)ioll.get_data(larlite::data::kMCTrack,"mcreco");

    std::cout << "Number of tracks in event: " << ev_mctrack->size() << std::endl;

    auto const& mctrack = ev_mctrack->at(i);

    std::cout << "TrackID is: " << mctrack.TrackID() << std::endl;
    std::cout << "Mother TrackID is: " << mctrack.MotherTrackID() << std::endl;
    std::cout << "PDG is: " << mctrack.PdgCode() << std::endl;

    std::cout << "Origin: " << mctrack.Origin() << std::endl;

    if ( mctrack.Origin() == 1 ) {
      producer = "simpleFlashBeam";
      isCosmic = 0;
    } else {
      producer = "simpleFlashCosmic";
      isCosmic = 1;
    }

    const larlite::mcstep& start = mctrack.Start();

    larutil::SpaceChargeMicroBooNE* _sce = nullptr;
    double tick = CrossingPointsAnaMethods::getTick(start, 4050.0, _sce);

    // check for primary muons (cosmic)
    if ( isCosmic == 1 && mctrack.TrackID() != mctrack.MotherTrackID() && mctrack.PdgCode() != 13 ) {
      tick = -999.997;
    }

    return tick;

  }

  double FlashMatcher::grabTickFromMCShower( larlite::storage_manager& ioll, int i ) {

    larlite::event_mcshower* ev_mcshower
      = (larlite::event_mcshower*)ioll.get_data(larlite::data::kMCShower,"mcreco");

    std::cout << "Number of showers in event: " << ev_mcshower->size() << std::endl;

    auto const& mcshower = ev_mcshower->at(i);

    std::cout << "TrackID is: " << mcshower.TrackID() << std::endl;
    std::cout << "Mother TrackID is: " << mcshower.MotherTrackID() << std::endl;
    std::cout << "PDG is: " << mcshower.PdgCode() << std::endl;

    std::cout << "Origin: " << mcshower.Origin() << std::endl;

    if ( mcshower.Origin() == 1 ) {
      producer = "simpleFlashBeam";
      isCosmic = 0;
    } else {
      producer = "simpleFlashCosmic";
      isCosmic = 1;
    }

    const larlite::mcstep& start = mcshower.Start();

    larutil::SpaceChargeMicroBooNE* _sce = nullptr;
    double tick = CrossingPointsAnaMethods::getTick(start, 4050.0, _sce);

    // check for primary muons (cosmic)
    if ( isCosmic == 1 && mcshower.TrackID() != mcshower.MotherTrackID() ) {
      tick = -999.997;
    }

    return tick;


  }

  /*
   * grab time coordinate from opflash -> convert to tick
   *
   * @param[in] opio The larlite storage manager that contains opflash class
   * @param[in] producer String labeling the producer, e.g. "simpleFlashBeam"
   * @return Vector containing all opflash times for the track in ticks
   */
  std::vector<double> FlashMatcher::grabTickFromOpflash( larlite::storage_manager& opio ) {

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
  double FlashMatcher::matchTicks( double mctrackTick, std::vector<double> flashTicks ) {

    double threshold;
    if (isCosmic == 1) {
        threshold = 2.0; // 1 tick = 0.5 us
    } else {
      threshold = std::numeric_limits<double>::infinity();
    }

    if (flashTicks.empty() == 1)
      return 999.999;

    auto match = std::lower_bound( flashTicks.begin(), flashTicks.end(), mctrackTick );
    double b = *(match);

    if (match == flashTicks.begin() && fabs(b - mctrackTick) <= threshold) {
      return flashTicks[0];
    }

    double a = *(match - 1);

    if (fabs(mctrackTick - a) < fabs(mctrackTick - b) && fabs(mctrackTick - a) <= threshold) {
      return flashTicks [ match - flashTicks.begin() - 1 ];
    }

    if ( fabs(mctrackTick - b) <= threshold  )
      return flashTicks[ match - flashTicks.begin() ];

    return -999.999;

  }


}
}
