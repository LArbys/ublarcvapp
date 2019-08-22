#ifndef __ASTAR_MASK_COMBO_H__
#define __ASTAR_MASK_COMBO_H__

#include <vector>

// larlite
#include "DataFormat/track.h"

// larcv
#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/EventChStatus.h"
#include "larcv/core/DataFormat/Image2D.h"

// ublarcvapp
#include "ublarcvapp/Reco3D/AStar3DAlgo.h"

// dltagger
#include "Gen3DEndpoints.h"

namespace ublarcvapp {
namespace dltagger {

  class AStarMaskCombo : public larcv::larcv_base {
  public:
    AStarMaskCombo()
      : larcv::larcv_base("AStarMaskCombo"),
      pEndpoint(nullptr),
      astar_completed(0)
    {};
    virtual ~AStarMaskCombo() {};

    AStarMaskCombo( const Gen3DEndpoints& endpt,
                    const larcv::EventChStatus& evchstatus,
                    bool run=true,
                    int max_downsample_factor=16,
                    int store_score_image=0,
                    int astar_verbosity=0);

    AStarMaskCombo( const Gen3DEndpoints& endpt,
                    const std::vector<larcv::Image2D>& whole_badch_v,
                    bool run=true,
                    int max_downsample_factor=16,
                    int store_score_image=0,
                    int astar_verbosity=0);
    
    std::vector<ublarcvapp::reco3d::AStar3DNode> astar_path;
    int astar_completed;
    std::vector< larcv::Image2D > badch_crop_v; //< badch image
    std::vector< larcv::Image2D > astar_crop_v; //< adc image with some mods using pcaline and masks
    std::vector< larcv::Image2D > score_crop_v; //< adc image showing score

    const Gen3DEndpoints *pEndpoint;
    
  protected:

    void _prep_crop_image( std::vector<larcv::Image2D>& input_v );
    void _prep_badch_crop( const std::vector<larcv::Image2D>& whole_badch_v,
                           const std::vector<larcv::Image2D>& input_crop_v );
    void _run_astar( const std::vector<larcv::Image2D>& input_v );

    // parameters
    int _max_downsample_factor; ///< maximum downsample factor to use when running AStar
    int _store_score_image;     ///< if 1, have AStar3DAlgo make and store score image for debug
    int _astar_verbosity;       ///< verbosity of the AStar3DAlgo
    
    
  };
  
}
}
    

#endif
