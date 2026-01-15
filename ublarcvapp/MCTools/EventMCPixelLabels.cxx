#include "EventMCPixelLabels.h"

namespace ublarcvapp {
namespace mctools{

    void EventMCPixelLabels::clear()
    {
        _triplets_v.clear();
        _imgcoord_to_tripindex.clear();
    }

}
}