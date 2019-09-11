#ifndef __BOUNDARYMUONTAGGERTYPES__
#define __BOUNDARYMUONTAGGERTYPES__

#include <string>

namespace ublarcvapp {
namespace tagger {

  typedef enum { kUndefined=-1, kTop, kBottom, kUpstream, kDownstream, kAnode, kCathode, kImageEnd, kNumEndTypes } BoundaryEnd_t;

  std::string BoundaryEndNames( const BoundaryEnd_t );

}
}

#endif
