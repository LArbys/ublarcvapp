#ifndef __ASTAR_NODES_2_BMTRACKCLUSTER3D_H__
#define __ASTAR_NODES_2_BMTRACKCLUSTER3D_H__

#include <vector>

// larcv
#include "larcv/core/DataFormat/Image2D.h"

// ublarcvapp
#include "ublarcvapp/Reco3D/AStar3DAlgo.h"
#include "BMTrackCluster3D.h"



namespace ublarcvapp {
namespace tagger {
  
  BMTrackCluster3D AStarNodes2BMTrackCluster3D( const std::vector<reco3d::AStar3DNode>& path,
                                                const std::vector<larcv::Image2D>& img_v,
						const BoundarySpacePoint& start_pt,
                                                const BoundarySpacePoint& end_pt,
                                                const float link_step_size );
  
}  
}

#endif
