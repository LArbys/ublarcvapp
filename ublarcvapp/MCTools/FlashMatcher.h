#ifndef __FLASH_MATCHER_H__
#define __FLASH_MATCHER_H__

/*
 * class for flash matching with mc truth
 */

// larlite/core
#include "larlite/DataFormat/storage_manager.h"

#include "TTree.h"

namespace ublarcvapp {
namespace mctools {

  class FlashMatcher {
  public:

    FlashMatcher() {

        isCosmic = 0;
        if (_fm_tree)
          delete _fm_tree;
        _fm_tree = nullptr;

    };
    virtual ~FlashMatcher() {};

    void initialize();
    void bindAnaVariables( TTree* );
    bool process(larlite::storage_manager& mgr);
    void finalize();
    void Clear();

    int numTracks( larlite::storage_manager& ioll );
    //int trackAncestorID( larlite::storage_manager& ioll, int i );

    // Getters for meta information
    int trackAncestorID();
    int getTrackID();
    int trackOrigin();

    int numShowers( larlite::storage_manager& ioll );
    double grabTickFromMCTrack( larlite::storage_manager& ioll, int i );
    double grabTickFromMCShower( larlite::storage_manager& ioll, int i );
    std::vector<double> grabTickFromOpflash( larlite::storage_manager& opio );
    double matchTicks( double mctrack_tick, std::vector<double> flash_ticks );

    Bool_t isCosmic;
    std::string producer;

  protected:

    TTree* _fm_tree;
    //TTree* _voxelTree;

  public:

    // Variables in TTree
    int _run;
    int _subrun;
    int _event;
    int _ancestorID;
    //int _trackID;
    double _clusterTick;
    double _flashTick;
    int _origin;

    // Vars not in tree but in loop
    int ancestorID;
    int trackID;
    double clusterTick;
    double flashTick;
    int origin;


  };

}
}

#endif
