#include "FlashMatcher.h"

#include "larlite/DataFormat/mcshower.h"
#include "larlite/DataFormat/mctrack.h"
#include "larlite/DataFormat/opflash.h"
#include "larlite/LArUtil/LArProperties.h"

#include "crossingPointsAnaMethods.h"

namespace ublarcvapp {
namespace mctools {


  void FlashMatcher::initialize()
  {
    _fm_tree = new TTree("fmtree","Flashmatched Tree");
    bindAnaVariables(_fm_tree);

    std::cout << "FLASHMATCHER: Made tree and bound vars" << std::endl;

    //_voxelTree = new TTree("voxtree","Voxeldata Tree");
    //_voxelTree = mgr.CloneTree()

  }

  void FlashMatcher::bindAnaVariables( TTree* fm_tree )
  {

    // event indexing
    fm_tree->Branch("flashmatcher_run",     &_run,     "flashmatcher_run/I");
    fm_tree->Branch("flashmatcher_subrun",  &_subrun,  "flashmatcher_subrun/I");
    fm_tree->Branch("flashmatcher_event",   &_event,   "flashmatcher_event/I");
    fm_tree->Branch("flashmatcher_ancestorid"  , &_ancestorID,   "flashmatcher_ancestorid/I");
    //fm_tree->Branch("flashmatcher_trackid"  , &_trackID,   "flashmatcher_trackid/I");
    fm_tree->Branch("flashmatcher_clustertick"  , &_clusterTick,   "flashmatcher_clustertick/D");
    fm_tree->Branch("flashmatcher_flashtick"  , &_flashTick,   "flashmatcher_flashtick/D");
    fm_tree->Branch("flashmatcher_origin"  , &_origin,   "flashmatcher_origin/I");

  }

  bool FlashMatcher::process(larlite::storage_manager& mgr)
  {

    Clear();

    _run    = mgr.run_id();
    _subrun = mgr.subrun_id();
    _event  = mgr.event_id();
    _ancestorID = ancestorID;
    //_trackID = trackID;
    _clusterTick = clusterTick;
    _flashTick = flashTick;
    _origin = origin;

    if ( _fm_tree )
      _fm_tree->Fill();

    return true;

  }

  void FlashMatcher::finalize()
  {
    if ( _fm_tree )
      _fm_tree->Write();
  }

  void FlashMatcher::Clear()
  {
    _run    = 0;
    _subrun = 0;
    _event  = 0;
    _ancestorID = 0;
    //_trackID = 0;
    _clusterTick = 0;
    _flashTick = 0;
    _origin = 0;

  }

  int FlashMatcher::numTracks( larlite::storage_manager& ioll ) {
    larlite::event_mctrack* ev_mctrack
      = (larlite::event_mctrack*)ioll.get_data(larlite::data::kMCTrack,"mcreco");

    int numTracks = ev_mctrack->size();

    return numTracks;
  }

/*
  int FlashMatcher::trackAncestorID( larlite::storage_manager& ioll, int i ) {

    larlite::event_mctrack* ev_mctrack
      = (larlite::event_mctrack*)ioll.get_data(larlite::data::kMCTrack,"mcreco");

    auto const& mctrack = ev_mctrack->at(i);

    auto ancestorID = mctrack.AncestorTrackID();

    return ancestorID;

  }
  */

  int FlashMatcher::trackAncestorID() {

    if (!ancestorID)
      return -9999;

    return ancestorID;

  }

  int FlashMatcher::getTrackID() {

    if (!trackID)
      return -9999;

    return trackID;

  }

  int FlashMatcher::trackOrigin() {

    if (!origin)
      return -9999;

    return origin;

  }

/*
  int FlashMatcher::trackOrigin( larlite::storage_manager& ioll, int i ) {


  }
  */

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

    ///std::cout << "Number of tracks in event: " << ev_mctrack->size() << std::endl;

    auto const& mctrack = ev_mctrack->at(i);

    trackID = mctrack.TrackID();
    ancestorID = mctrack.AncestorTrackID();

    ///std::cout << "TrackID is: " << mctrack.TrackID() << std::endl;
    ///std::cout << "Mother TrackID is: " << mctrack.MotherTrackID() << std::endl;
    ///std::cout << "PDG is: " << mctrack.PdgCode() << std::endl;

    ///std::cout << "Origin: " << mctrack.Origin() << std::endl;
    origin = mctrack.Origin();

    if ( mctrack.Origin() == 1 ) {
      producer = "simpleFlashBeam";
      isCosmic = 0;
    } else {
      producer = "simpleFlashCosmic";
      isCosmic = 1;
    }

    ///std::cout << "mctrack.size() = " << mctrack.size() << std::endl;

    if ( mctrack.size() == 0 ) {
      ///std::cout << "Empty vector of mcsteps (probably a cosmic that didn't cross the TPC)" << std::endl;
      return -9999.998;
    }

    //const larlite::mcstep& start = mctrack.Start();
    const larlite::mcstep& start = mctrack.at(0);

    ///std::cout << "First mcstep X positiom is " << start.X() << " while the track started at: " << mctrack.Start().X() << std::endl;

    larutil::SpaceChargeMicroBooNE* _sce = nullptr;

    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    double xPos = start.X();
    double xPos2 = mctrack.Start().X();

    double tick = CrossingPointsAnaMethods::getTick(start, 4050.0, _sce);
    tick = tick - xPos / cm_per_tick;

    double tick2 = CrossingPointsAnaMethods::getTick(mctrack.Start(), 4050.0, _sce);
    tick2 = tick2 - xPos2 / cm_per_tick;

    ///std::cout << "Tick for first mcstep in TPC: " << tick << " and tick for starting up in the sky: " << tick2 << std::endl;

    // check for primaries
    if ( mctrack.TrackID() != mctrack.MotherTrackID() ) {
      return -9999.997;
    }

    ///std::cout << "Ancestor ID: " << mctrack.AncestorTrackID() << std::endl;

    return tick;

  }

  double FlashMatcher::grabTickFromMCShower( larlite::storage_manager& ioll, int i ) {

    larlite::event_mcshower* ev_mcshower
      = (larlite::event_mcshower*)ioll.get_data(larlite::data::kMCShower,"mcreco");

    ///std::cout << "Number of showers in event: " << ev_mcshower->size() << std::endl;

    auto const& mcshower = ev_mcshower->at(i);

    ///std::cout << "TrackID is: " << mcshower.TrackID() << std::endl;
    ///std::cout << "Mother TrackID is: " << mcshower.MotherTrackID() << std::endl;
    ///std::cout << "PDG is: " << mcshower.PdgCode() << std::endl;

    ///std::cout << "Origin: " << mcshower.Origin() << std::endl;
    origin = mcshower.Origin();

    if ( mcshower.Origin() == 1 ) {
      producer = "simpleFlashBeam";
      isCosmic = 0;
    } else {
      producer = "simpleFlashCosmic";
      isCosmic = 1;
    }

    const larlite::mcstep& start = mcshower.Start();

    larutil::SpaceChargeMicroBooNE* _sce = nullptr;

    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    double xPos = start.X();

    double tick = CrossingPointsAnaMethods::getTick(start, 4050.0, _sce);
    tick = tick - xPos / cm_per_tick;

    trackID = mcshower.TrackID();
    ancestorID = mcshower.AncestorTrackID();

    // check for primaries
    if ( mcshower.TrackID() != mcshower.MotherTrackID() ) {
      return -9999.997;
    }

    ///std::cout << "Ancestor ID: " << mcshower.AncestorTrackID() << std::endl;

    return tick;


  }

  /*
   * grab time coordinate from opflash -> convert to tick
   *
   * @param[in] opio The larlite storage manager that contains opflash class
   * @param[in] producer String labeling the producer, e.g. "simpleFlashBeam"
   * @return Vector containing all opflash times for the track in ticks
   */
 std::vector< std::pair<double, int> > FlashMatcher::grabTickFromOpflash( larlite::storage_manager& opio ) {

    std::vector< std::pair<double, int> > tickContainer;

    larlite::event_opflash* ev_opflash
      = (larlite::event_opflash*)opio.get_data(larlite::data::kOpFlash, producer);

    int counter = 0;

    for (auto const& opflash : *ev_opflash) {
      double time = opflash.Time();
      //std::cout << time << std::endl;
      double tick = time/0.5 + 3200.0;
      tickContainer.push_back(std::make_pair(tick,counter));
      counter++;
    }

    /*
    std::cout << "Before sort: " << std::endl;
    for (int i = 0; i < tickContainer.size(); i++) {
      std::cout << tickContainer[i].first << "\t"
               << tickContainer[i].second << std::endl;
      }
    */
  //std::cout << "Tick container: " << tickContainer << std::endl;

  std::sort( tickContainer.begin(), tickContainer.end() );

  /*
  std::cout << "After sort: " << std::endl;
    for (int i = 0; i < tickContainer.size(); i++) {
        std::cout << tickContainer[i].first << "\t"
             << tickContainer[i].second << std::endl;
    }
  */

  //std::cout << "Tick container: " << tickContainer << std::endl;

  return tickContainer;
  }

  /*
   * matches mctrack tick to opflash tick
   * for cosmics, will return -9999.999 if there is no opflash found within the threshold
   * for beam tracks, will find the closest matching opflash
   *
   * @param[in] mctrackTick Time in ticks for the mctrack to be matched to
   * @param[in] flashTicks Vector of potential opflash matches in ticks
   * @param[in] isCosmic Cosmic flag to determine threshold
   *
   * @return Value of the closest matching opflash tick
   */
std::pair<double, int> FlashMatcher::matchTicks( double mctrackTick, std::vector< std::pair<double, int> > flashTicks ) {

    double threshold;
    if (isCosmic == 1) {
        threshold = 2.0; // 1 tick = 0.5 us
    } else {
      threshold = std::numeric_limits<double>::infinity();
    }

    if (flashTicks.empty() == 1) {
      clusterTick = mctrackTick;
      flashTick = 9999.999;
      return std::make_pair(9999.999,9999);
    }

    auto match = std::lower_bound( flashTicks.begin(), flashTicks.end(), std::make_pair( mctrackTick, std::numeric_limits<int>::min()) );
    double b = (*(match)).first;

    //std::cout << "match [iterator] is: " << match << std::endl;
    //std::cout << "*match [dereferenced iterator, shoudl be a pair] is: " << (*(match)) << std::endl;
    ///std::cout << "(*(match)).first [actual value of the tick] is: " << b << std::endl;

    if (match == flashTicks.begin() && fabs(b - mctrackTick) <= threshold) {
      clusterTick = mctrackTick;
      flashTick = flashTicks[0].first;
      return flashTicks[0];
    }

    double a = (*(match - 1)).first;

    if (fabs(mctrackTick - a) < fabs(mctrackTick - b) && fabs(mctrackTick - a) <= threshold) {
      clusterTick = mctrackTick;
      flashTick = flashTicks [ match - flashTicks.begin() - 1 ].first;
      return flashTicks [ match - flashTicks.begin() - 1 ];
    }

    if ( fabs(mctrackTick - b) <= threshold  ) {
      clusterTick = mctrackTick;
      flashTick = flashTicks[ match - flashTicks.begin() ].first;
      return flashTicks[ match - flashTicks.begin() ];
    }

    return std::make_pair(-9999.999,9999);

  }


}
}
