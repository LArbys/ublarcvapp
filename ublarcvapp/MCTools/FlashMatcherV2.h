#ifndef __FLASH_MATCHER_V2_H__
#define __FLASH_MATCHER_V2_H__

#include <string>
#include <set>
#include "TTree.h"

// larlite/core
#include "larlite/DataFormat/storage_manager.h"

// ublarcvapp/mctools
#include "ublarcvapp/MCTools/MCPixelPGraph.h"

namespace ublarcvapp {
namespace mctools {

  /**
   * @class RecoFlash_t
   *
   * class for storing mc tracks associated to a flash
   *
   * can represent null flashes as well for charge depositing tracks
   * that produce zero light in the PMTs. in this case the producerid
   * will be 'null'.
   *
   */
  class RecoFlash_t {
  public:
    
    // default constructor
    RecoFlash_t()
      : producerid(-1),
	index(-1),
	used(0),
	ancestorid(-1),
	time_us(0.0)
    {
      trackid_v.clear();
    };
    ~RecoFlash_t() {};
    
    int producerid;
    int index;
    int used;
    int ancestorid;  ///< truth-matched ancestor ID of the particle cascade that made the flash
    float time_us;   ///< relative to optical system trigger
    float tick;      ///< TPC clock tick
    
    float operator<( const RecoFlash_t& rhs ) const {
      if ( time_us < rhs.time_us ) return true;
      return false;
    };
    
    std::set<int> trackid_v; // geant track ids
    
  };


  /**
   * @class FlashMatcherV2
   *
   * class for matching reco optical flashes to true track information.
   *
   *
   */  
  class FlashMatcherV2 {
  public:

    FlashMatcherV2()
      : _dtick_threshold(1.0),
	_verbose_level(0)
    {
      recoflash_v.clear();
    };
    virtual ~FlashMatcherV2() {};

    bool process(larlite::storage_manager& mgr);
    void printMatches();

    int numTracks( larlite::storage_manager& ioll );
    int numShowers( larlite::storage_manager& ioll );    

    void buildRecoFlashPool( larlite::storage_manager& ioll );
    void matchTracksAndFlashes( larlite::storage_manager& ioll );
    void associateTruthTrackIDs2recoFlashes( larlite::storage_manager& ioll,
					     ublarcvapp::mctools::MCPixelPGraph& mgpg );
    void buildNullFlashesToTrackIDs( larlite::storage_manager& ioll,
				     ublarcvapp::mctools::MCPixelPGraph& mgpg );
    
    void setVerboseLevel( int level ) { _verbose_level=level; };

  protected:

    float _dtick_threshold; ///< absolute time difference between flash and mc track object to match
    int _verbose_level; ///< control output to stdout. [0] quiet (default), [1] info, [2] debug

  public:


    // container with matches
    std::vector<RecoFlash_t> recoflash_v;
    
  };

}
}

#endif
