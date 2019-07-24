#ifndef __ublarcvapp_DBScan_h__
#define __ublarcvapp_DBScan_h__

#include "DBScanTypes.h"

namespace ublarcvapp {

  class DBScan {

  public:
    DBScan() {};
    virtual ~DBScan() {};
    

    static std::vector< dbscan::Cluster_t > makeCluster( const float maxdist, const float minhits, const int maxkdneighbors,
                                                         const std::vector<std::vector<float> >& clust );
    
    
  };

}

#endif
