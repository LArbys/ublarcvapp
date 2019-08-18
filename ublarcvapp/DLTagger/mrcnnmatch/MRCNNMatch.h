#ifndef __MRCNN_MATCH_H__
#define __MRCNN_MATCH_H__

#include <vector>

#include "larcv/core/DataFormat/ClusterMask.h"

#include "MRCNNMatch.h"

namespace ublarcvapp {
namespace dltagger {

  class MRCNNMatch {
  public:
    MRCNNMatch() {};
    virtual ~MRCNNMatch() {};

    void matchMasksAcrossPlanes( const std::vector<std::vector<larcv::ClusterMask>>& clustermask_vv,
                                 std::vector< std::vector<int> >& match_indices );

  };
}
}

#endif
