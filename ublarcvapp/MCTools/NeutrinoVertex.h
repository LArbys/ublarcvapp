#ifndef __NEUTRINO_VERTEX_H__
#define __NEUTRINO_VERTEX_H__

/*
 * class providing some static functions to quickly get neutrino vertex
 * either in detector 3D coordinates or in the image
 *
 */

#include <vector>

// larlite/core
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"

namespace ublarcvapp {
namespace mctools {

  class NeutrinoVertex {
  public:
    
    NeutrinoVertex() {};
    virtual ~NeutrinoVertex() {};

    static std::vector<int>   getImageCoords( larlite::storage_manager& ioll );
    static std::vector<float> getPos3D( larlite::storage_manager& ioll );
    static std::vector<float> getPos3DwSCE( larlite::storage_manager& ioll, larutil::SpaceChargeMicroBooNE* psce=0 );
    
  };
  
}
}

#endif
