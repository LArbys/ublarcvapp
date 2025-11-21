#ifndef __UBLARCVAPP_MCTOOLS_EVENT_MCPIXEL_LABELS_H__
#define __UBLARCVAPP_MCTOOLS_EVENT_MCPIXEL_LABELS_H__

#include <vector>
#include <map>
#include <array>

#include "MCPixelLabels.h"

namespace ublarcvapp {
namespace mctools {

class EventMCPixelLabels {

public:

  EventMCPixelLabels(){};
  ~EventMCPixelLabels(){};

  void clear();

  std::vector<MCPixelLabels> _triplets_v; //< container of triplet info
  std::map< std::array<int,4>, unsigned long > _imgcoord_to_tripindex; //< (u,v,y,row) to position in _triplets_v

  //MCPixelLabels& get( int u, int v, int y, )

};

}
}

#endif