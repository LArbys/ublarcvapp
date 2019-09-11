#ifndef __ASTAR_3D_ALGO_H__
#define __ASTAR_3D_ALGO_H__

/** 

AStar algorithm assuming 3D!! grid points. 

Uses Image2D to hold image.

 **/

#include <iostream>
#include <queue>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <array>

// larcv
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/Base/PSet.h"

#include "AStar3DAlgoConfig.h"
#include "AStar3DTypes.h"
#include "Lattice.h"

namespace ublarcvapp {
namespace reco3d {

  // ALGO 
  class AStar3DAlgo {

    AStar3DAlgo() { verbose=2; };
  public:

    AStar3DAlgo( AStar3DAlgoConfig config ) { 
      _config = config; 
      setVerbose( _config.verbosity );
    };
    virtual ~AStar3DAlgo() {};
    
    void setVerbose( int v ) { verbose = v; };

    std::vector<AStar3DNode> findpath( const std::vector<larcv::Image2D>& img_v,
                                       const std::vector<larcv::Image2D>& badch_v,
                                       const std::vector<larcv::Image2D>& tagged_v,
                                       const int start_row, const int goal_row,
                                       const std::vector<int>& start_cols,
                                       const std::vector<int>& goal_cols,
                                       int& goal_reached );

    std::vector<AStar3DNode> downsampleAndFindPath( const int downsampling_factor, const std::vector<larcv::Image2D>& img_v,
						    const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v,
						    const int start_row, const int goal_row,
						    const std::vector<int>& start_cols, const std::vector<int>& goal_cols, int& goal_reached,
						    int compression_mode=-1 );
    
    const std::vector<larcv::Image2D>& getScoreImages() { return m_visualizedimgs; }
    
    // supporting 
    std::vector<AStar3DNode> makeRecoPath( AStar3DNode* start, AStar3DNode* goal, bool& path_completed );  

    void evaluateNeighborNodes( AStar3DNode* current, const AStar3DNode* start, const AStar3DNode* goal,
      const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v,
      AStar3DNodePtrList& openset, AStar3DNodePtrList& closedset, Lattice& lattice );

    bool evaluteLatticePoint( const A3DPixPos_t& latticept, AStar3DNode* current, const AStar3DNode* start, const AStar3DNode* goal,
      const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v,
      AStar3DNodePtrList& openset, AStar3DNodePtrList& closedset, Lattice& lattice );

    float distanceFromCentralLine( const std::vector<float>& start_tyz, const std::vector<float>& end_tyz, const std::vector<float>& testpt_tyz );


  protected:

    AStar3DAlgoConfig _config;
    int verbose;

    std::vector<larcv::Image2D> visualizeScores( std::string score_name, const std::vector<larcv::Image2D>& orig_img_v, Lattice& lattice );    
    
    // void evaluateBadChNeighbors( AStar3DNode* current, const AStar3DNode* start, const AStar3DNode* goal,
   //    AStar3DNodePtrList& openset, AStar3DNodePtrList& closedset, 
   //    const int neighborhood_size, const int min_c, const int min_r, const int win_c, const int win_r, 
   //    const larcv::Image2D& img, const larcv::ImageMeta& meta, const bool use_bad_chs, std::map< A3DPixPos_t, AStar3DNode* >& position_lookup);

    std::vector< larcv::Image2D > m_visualizedimgs;
  };


}  
}

#endif
