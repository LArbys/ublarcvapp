#ifndef __ASTAR_MASK_COMBO_H__
#define __ASTAR_MASK_COMBO_H__

#include <vector>

// larlite
#include "DataFormat/track.h"

// larcv
#include "larcv/core/DataFormat/EventChStatus.h"
#include "larcv/core/DataFormat/Image2D.h"

// ublarcvapp
#include "ublarcvapp/Reco3D/AStar3DAlgo.h"

// dltagger
#include "Gen3DEndpoints.h"

namespace ublarcvapp {
namespace dltagger {

  class AStarMaskCombo {
  public:
    AStarMaskCombo()
      : pEndpoint(nullptr),
      astar_completed(0)
    {};
    virtual ~AStarMaskCombo() {};

    AStarMaskCombo( const Gen3DEndpoints& endpt, const larcv::EventChStatus& evchstatus, bool run=true );
    AStarMaskCombo( const Gen3DEndpoints& endpt, const std::vector<larcv::Image2D>& whole_badch_v, bool run=true );
    
    std::vector<ublarcvapp::reco3d::AStar3DNode> astar_path;
    int astar_completed;
    std::vector< larcv::Image2D > badch_crop_v; //< badch image
    std::vector< larcv::Image2D > astar_crop_v; //< adc image with some mods using pcaline and masks

    const Gen3DEndpoints *pEndpoint;
    
  protected:
    
    void _prep_badch_crop( const std::vector<larcv::Image2D>& whole_badch_v );
    void _run_astar();
    
  };
  
}
}
    

#endif
