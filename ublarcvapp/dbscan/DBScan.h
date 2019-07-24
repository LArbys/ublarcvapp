#ifndef __ublarcvapp_DBScan_h__
#define __ublarcvapp_DBScan_h__

#include "DBScanTypes.h"

namespace ublarcvapp {
namespace dbscan {

  class DBScan {

  public:
    DBScan() {};
    virtual ~DBScan() {};
    
    template< typename T >
      static std::vector< dbCluster > makeCluster( const float maxdist, const float minhits, const int maxkdneighbors,
                                                   const std::vector<std::vector<T> >& pointsxyz );
    static std::vector< dbCluster > makeCluster3f( const float maxdist, const float minhits, const int maxkdneighbors,
                                                   const std::vector<std::vector<float> >& pointsxyz );
    static std::vector< dbCluster > makeCluster3d( const float maxdist, const float minhits, const int maxkdneighbors,
                                                   const std::vector<std::vector<double> >& pointsxyz );
  };

}
}

#endif
