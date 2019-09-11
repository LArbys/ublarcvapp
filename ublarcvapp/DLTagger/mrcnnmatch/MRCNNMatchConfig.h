#ifndef __MRCNNMatchConfig_H__
#define __MRCNNMatchConfig_H__

namespace ublarcvapp {
namespace dltagger {

  class MRCNNMatchConfig {

  public:

    MRCNNMatchConfig();
    virtual ~MRCNNMatchConfig() {};

    void setDefaults();

    float triarea_maximum; //< maximum score that measures 2D (y,z) detector position consistency of 3 image points, one on each plane.
    bool filter_astar_failures; //< if true, do not save matches that fail the astar track reco

    // sub-algo parameters
    int astar_max_downsample_factor;
    int astar_store_score_image;
    int astar_verbosity;
    
  };
  
}
}

    


#endif
