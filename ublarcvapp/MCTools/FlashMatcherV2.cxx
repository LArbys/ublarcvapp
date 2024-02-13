#include "FlashMatcherV2.h"

#include <algorithm>

#include "larlite/DataFormat/mcshower.h"
#include "larlite/DataFormat/mctrack.h"
#include "larlite/DataFormat/opflash.h"
#include "larlite/LArUtil/LArProperties.h"

// ublarcvapp/mctools
#include "crossingPointsAnaMethods.h"

namespace ublarcvapp {
namespace mctools {



  bool FlashMatcherV2::process(larlite::storage_manager& mgr)
  {

    // first get list of times
    buildRecoFlashPool( mgr );

    // match truth tracks to flashes
    matchTracksAndFlashes( mgr );

    // list results
    if ( _verbose_level>=1 )
      printMatches();


    return true;

  }


  int FlashMatcherV2::numTracks( larlite::storage_manager& ioll ) {
    larlite::event_mctrack* ev_mctrack
      = (larlite::event_mctrack*)ioll.get_data(larlite::data::kMCTrack,"mcreco");

    int numTracks = ev_mctrack->size();

    return numTracks;
  }

  int FlashMatcherV2::numShowers( larlite::storage_manager& ioll ) {

    larlite::event_mcshower* ev_mcshower
      = (larlite::event_mcshower*)ioll.get_data(larlite::data::kMCShower,"mcreco");

    int numShowers = ev_mcshower->size();

    return numShowers;

  }

  // /*
  //  * grab time coordinate from mctrack mcstep -> convert to tick
  //  *
  //  * @param[in] ioll The larlite storage manager that contains mctruth class
  //  * @return tuple with tick, producer string, and cosmic flag
  //  */
  // double FlashMatcherV2::grabTickFromMCTrack( larlite::storage_manager& ioll, int i ) {
  //   larlite::event_mctrack* ev_mctrack
  //     = (larlite::event_mctrack*)ioll.get_data(larlite::data::kMCTrack,"mcreco");

  //   std::cout << "Number of tracks in event: " << ev_mctrack->size() << std::endl;

  //   auto const& mctrack = ev_mctrack->at(i);

  //   std::cout << "TrackID is: " << mctrack.TrackID() << std::endl;
  //   std::cout << "Mother TrackID is: " << mctrack.MotherTrackID() << std::endl;
  //   std::cout << "PDG is: " << mctrack.PdgCode() << std::endl;

  //   std::cout << "Origin: " << mctrack.Origin() << std::endl;
  //   origin = mctrack.Origin();

  //   if ( mctrack.Origin() == 1 ) {
  //     producer = "simpleFlashBeam";
  //     isCosmic = 0;
  //   } else {
  //     producer = "simpleFlashCosmic";
  //     isCosmic = 1;
  //   }

  //   const larlite::mcstep& start = mctrack.Start();

  //   larutil::SpaceChargeMicroBooNE* _sce = nullptr;

  //   const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
  //   double xPos = start.X();

  //   double tick = CrossingPointsAnaMethods::getTick(start, 4050.0, _sce);
  //   tick = tick - xPos / cm_per_tick;

  //   // check for primaries
  //   if ( mctrack.TrackID() != mctrack.MotherTrackID() ) {
  //     return -999.997;
  //   }

  //   std::cout << "Ancestor ID: " << mctrack.AncestorTrackID() << std::endl;

  //   trackID = mctrack.TrackID();
  //   ancestorID = mctrack.AncestorTrackID();

  //   return tick;

  // }

  // double FlashMatcherV2::grabTickFromMCShower( larlite::storage_manager& ioll, int i ) {

  //   larlite::event_mcshower* ev_mcshower
  //     = (larlite::event_mcshower*)ioll.get_data(larlite::data::kMCShower,"mcreco");

  //   std::cout << "Number of showers in event: " << ev_mcshower->size() << std::endl;

  //   auto const& mcshower = ev_mcshower->at(i);

  //   std::cout << "TrackID is: " << mcshower.TrackID() << std::endl;
  //   std::cout << "Mother TrackID is: " << mcshower.MotherTrackID() << std::endl;
  //   std::cout << "PDG is: " << mcshower.PdgCode() << std::endl;

  //   std::cout << "Origin: " << mcshower.Origin() << std::endl;
  //   origin = mcshower.Origin();

  //   if ( mcshower.Origin() == 1 ) {
  //     producer = "simpleFlashBeam";
  //     isCosmic = 0;
  //   } else {
  //     producer = "simpleFlashCosmic";
  //     isCosmic = 1;
  //   }

  //   const larlite::mcstep& start = mcshower.Start();

  //   larutil::SpaceChargeMicroBooNE* _sce = nullptr;

  //   const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
  //   double xPos = start.X();

  //   double tick = CrossingPointsAnaMethods::getTick(start, 4050.0, _sce);
  //   tick = tick - xPos / cm_per_tick;

  //   // check for primaries
  //   if ( mcshower.TrackID() != mcshower.MotherTrackID() ) {
  //     return -999.997;
  //   }

  //   std::cout << "Ancestor ID: " << mcshower.AncestorTrackID() << std::endl;

  //   return tick;


  // }


  /**
   * @brief inspect reco flash containers and store RecoFlash_t objects that store flash times
   *
   * When this is run, the member container \ref FlashMatcherV2.recoflash_v will be populated.
   * This is a vector of RecoFlash_t objects that stores a list of trackIDs to each associated flash.
   *
   *
   *
   * params
   * ------
   * @param[in] ioll larlite::storage_manager where an entry has already been loaded.
   *
   */
  void FlashMatcherV2::buildRecoFlashPool( larlite::storage_manager& ioll ) {

    // clear flash pools
    recoflash_v.clear();
    
    std::vector<std::string> producer_v = {"simpleFlashBeam","simpleFlashCosmic"} ;
    
    for (int iproducer=0; iproducer<(int)producer_v.size(); iproducer++) {
      std::string producer = producer_v.at(iproducer);

      larlite::event_opflash* ev_opreco
	= (larlite::event_opflash*)ioll.get_data(larlite::data::kOpFlash, producer );
      if ( _verbose_level>=2 )
	std::cout << "Reco flashes in producer[" << producer << "]: " << ev_opreco->size() << std::endl;

      for (int iflash=0; iflash<(int)ev_opreco->size(); iflash++) {

	auto const& dataflash = ev_opreco->at(iflash);
	
	RecoFlash_t flash;
	flash.producerid = iproducer;
	flash.index = iflash;
	flash.used = 0;
	flash.time_us = dataflash.Time();
	flash.tick    = dataflash.Time()/0.5 + 3200; // assuming optical time origin set at TPC clock tick 3200
	recoflash_v.emplace_back( std::move(flash) );
      };
    }

    // sort all the optical flashes by time
    std::sort( recoflash_v.begin(), recoflash_v.end() );
    if ( _verbose_level>=1 )
      std::cout << "Number of reco flashes in the (truth) matching pool: " << recoflash_v.size() << std::endl;
    
  };

  /**
   * @brief match flashes and to MC tracks
   *
   * 
   */
  void FlashMatcherV2::matchTracksAndFlashes( larlite::storage_manager& ioll )
  {

    // use the mc particle graph tool
    MCPixelPGraph mcpg;
    mcpg.buildgraphonly( ioll );

    associateTruthTrackIDs2recoFlashes( ioll, mcpg );
    buildNullFlashesToTrackIDs( ioll, mcpg );

    std::sort( recoflash_v.begin(), recoflash_v.end() );
    
  }

  /**
   * @brief associate trackIDs from the mctrack and mcshower containers to the reco flashes
   *
   */
  void FlashMatcherV2::associateTruthTrackIDs2recoFlashes( larlite::storage_manager& ioll,
							 ublarcvapp::mctools::MCPixelPGraph& mcpg )
  {
    
    if ( _verbose_level>=2 ) 
      mcpg.printGraph(nullptr,false);

    // get the primary list
    auto node_v = mcpg.getPrimaryParticles();

    // loop over the primary nodes: these are the recorded list of primary particles
    // defined as trackid=ancestorid
    for ( auto const& node : node_v ) {
      if ( _verbose_level>=2 )
	std::cout << "primary node: trackid=" << node->tid << " ancestorid=" << node->aid << " origin=" << node->origin << std::endl;

      std::vector<float> txyz = { node->start[3] , node->start[0], node->start[1], node->start[2] };

      float tpctick_nodrift = CrossingPointsAnaMethods::getTrueTick( txyz, 4050.0, NULL );
      if ( _verbose_level>=2 )      
	std::cout << "  tpctick_nodrift=" << tpctick_nodrift << std::endl;

      int closest_flashidx = -1;
      float min_dtick = 1.0e9;
      for ( int idx=0; idx<(int)recoflash_v.size(); idx++ ) {
	float dtick = fabs(recoflash_v[idx].tick-tpctick_nodrift);
	if (dtick<min_dtick ) {
	  min_dtick = dtick;
	  closest_flashidx = idx;
	}
      }
      
      if ( _verbose_level>=2 )      
	std::cout << "  closest flashidx=" << closest_flashidx << "  dtick=" << min_dtick << std::endl;
      
      if ( min_dtick < _dtick_threshold && closest_flashidx>=0 ) {
	// we've matched this track to a flash

	// get the reco flash that best matched to the primary
	auto& flash = recoflash_v[closest_flashidx];

	// assign this track ID to the flash
	if ( node->origin==2 ) {
	  // cosmic origin track
	  if ( flash.ancestorid<0 ) {
	    flash.ancestorid = node->aid;
	  }
	  else if (flash.ancestorid>=0 && flash.ancestorid!=node->aid) {
	    std::cout << "  WARNING: flash already matched to node with another ancestorid! old=" << flash.ancestorid << std::endl;
	  }
	}
	else if (node->origin==1 ) {
	  // neutrino origin track, set the flash as neutrino origin
	  flash.ancestorid = 0;
	}
	flash.trackid_v.insert( node->tid );
	
	// get all the descendent particle records matched to this primary
	auto node_et_daughters = mcpg.getNodeAndDescendentsFromTrackID( node->tid );
	int nadded = 0;
	for ( auto& dnode : node_et_daughters ) {
	  flash.trackid_v.insert( dnode->tid );
	  nadded++;
	}
	if ( _verbose_level>=2 )
	  std::cout << "  inserted " << nadded << " trackids to flash. Flash len(track_id)=" << flash.trackid_v.size() << std::endl;
	
      }//end of if tick difference is below max threshold
      
    }//end of loop over vector of primary node pointers

    return;
  }

  /**
   * @brief find trackIDs that did not have a reco flash
   * 
   */
  void FlashMatcherV2::buildNullFlashesToTrackIDs( larlite::storage_manager& ioll,
						 ublarcvapp::mctools::MCPixelPGraph& mcpg )
  {

    // we need a list of ancestor ids where we've already matched
    std::set< int > matched_ancestor_ids;
    for ( auto const& recoflash : recoflash_v ) {
      if ( recoflash.ancestorid>0 ) {
	matched_ancestor_ids.insert( recoflash.ancestorid );
      }
      else if ( recoflash.ancestorid==0 ) {
	// for neutrino cluster, the trackids are their own ancestorids
	for (auto const& trackid : recoflash.trackid_v ) {
	  matched_ancestor_ids.insert( trackid );
	}
      }
    }

    // container for null flashes
    std::vector< RecoFlash_t > null_v; 
    
    // now look for ancestor ids not in the above list
    auto node_v = mcpg.getPrimaryParticles();
    for ( auto const& node : node_v ) {
      
      auto it_aid = matched_ancestor_ids.find( node->aid );
      if ( it_aid!=matched_ancestor_ids.end() )
	continue; // we've already matched this ancestor id
      
      // not in the list, make a null flash and associate
      std::vector<float> txyz = { node->start[3] , node->start[0], node->start[1], node->start[2] };
      float tpctick_nodrift = CrossingPointsAnaMethods::getTrueTick( txyz, 4050.0, NULL );
      
      RecoFlash_t nullflash;
      nullflash.producerid = -1;
      nullflash.index = int(null_v.size());
      nullflash.used = 0;
      nullflash.ancestorid = node->aid;
      nullflash.tick = tpctick_nodrift;
      nullflash.time_us = (tpctick_nodrift-3200.0)*(0.5); // (tick diff from trigger tick)*(0.5 usec per tick)
      nullflash.trackid_v.insert( node->tid );
      
      // put in daughters
      auto node_et_daughters = mcpg.getNodeAndDescendentsFromTrackID( node->tid );      
      int nadded = 0;
      for ( auto& dnode : node_et_daughters ) {
	nullflash.trackid_v.insert( dnode->tid );
	nadded++;
      }

      if ( _verbose_level>=2 ) {
	std::cout << "  created null flash for ancestorid=" << node->aid << std::endl;
	std::cout << "    inserted " << nadded << " trackids to flash. Flash len(track_id)=" << nullflash.trackid_v.size() << std::endl;
      }

      null_v.emplace_back( std::move(nullflash) );
    }

    if ( _verbose_level>=1 ) {
      std::cout << "created " << null_v.size() << " null flashes" << std::endl;
    }

    for (int ii=0; ii<(int)null_v.size(); ii++) {
      recoflash_v.push_back( null_v.at(ii) );
    }
      
  }

  void FlashMatcherV2::printMatches() {
    
    std::cout << "=================================" << std::endl;
    std::cout << " FlashMatcherV2::printMatches" << std::endl;
    
    for (int iflash=0; iflash<(int)recoflash_v.size(); iflash++) {
      auto const& flash = recoflash_v.at(iflash);
      std::stringstream flashinfo;
      flashinfo << " flash[" << iflash << "]"
		<< " producer[" << flash.producerid << "]"
		<< " time_us=" << flash.time_us
		<< " tick=" << flash.tick
		<< " aid=" << flash.ancestorid
		<< " matched={ ";
      for ( auto const& tid : flash.trackid_v ) {
	flashinfo << tid << " ";
      }
      flashinfo << "}";
      std::cout << flashinfo.str() << std::endl;
    }
    std::cout << "==================================" << std::endl;
  }


}
}
