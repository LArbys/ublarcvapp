#include "MRCNNMatchConfig.h"

namespace ublarcvapp {
namespace dltagger {

  MRCNNMatchConfig::MRCNNMatchConfig() {
    setDefaults();
  }


  void MRCNNMatchConfig::setDefaults() {
    
    triarea_maximum = 200;
    filter_astar_failures = false;

    astar_max_downsample_factor = 16;
    astar_store_score_image = 0;
    astar_verbosity = 1;
    
  }

  
}
}
