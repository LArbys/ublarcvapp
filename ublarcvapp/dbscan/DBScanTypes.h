#ifndef __UBLARCVAPP_DBSCANTYPES_H__
#define __UBLARCVAPP_DBSCANTYPES_H__

#include <vector>

namespace ublarcvapp {
namespace dbscan {
  
  typedef std::vector< std::vector<double> > dbPoints; // list of (x,y,z,....) points
  typedef std::vector<int> dbCluster; // list of indices to dbPoints
  typedef std::vector< std::vector<int> >  dbClusters; // list of list of indices to provided dbPoints
  
}
}



#endif
