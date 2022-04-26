#ifndef __MCPIXELPGRAPH_H__
#define __MCPIXELPGRAPH_H__

#include <vector>

// larcv
#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/Image2D.h"

// larlite
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/DataFormat/mcshower.h"
#include "larlite/DataFormat/mctrack.h"
#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"

/**
 * Determine particle graph. Collect pixels for each particle.
 * In particular, handle the showers.
 *
 */

namespace ublarcvapp {
namespace mctools {

  class MCPixelPGraph : public larcv::larcv_base {
  public:

    MCPixelPGraph()
      : larcv::larcv_base("MCPixelPGraph"),
      adc_tree("wire")
      {};
    virtual ~MCPixelPGraph() {};

    void buildgraph( larcv::IOManager& iolcv, larlite::storage_manager& ioll );
    
    void buildgraph( const std::vector<larcv::Image2D>& adc_v,
                     const std::vector<larcv::Image2D>& segment_v,
                     const std::vector<larcv::Image2D>& instance_v,
                     const std::vector<larcv::Image2D>& ancestor_v,
                     const larlite::event_mcshower& shower_v,
                     const larlite::event_mctrack&  track_v,
                     const larlite::event_mctruth&  mctruth_v );

    void buildgraphonly( larlite::storage_manager& ioll );    
    void buildgraphonly( const larlite::event_mcshower& shower_v,
                         const larlite::event_mctrack&  track_v,
                         const larlite::event_mctruth&  mctruth_v );
    


    struct Node_t {
      int nodeidx;    // position in node_v
      int type;       // track=0, shower=1
      int vidx;       // position in mcshower or mctrack vector
      int tid;        // geant4 track-ID
      int aid;        // ancestor geant4 trackid
      int mtid;       // mother geant4 trackid
      int pid;        // particle ID
      Node_t* mother; // pointer to Mother Node_t
      int mid;        // mother nodeidx
      float E_MeV;    // energy
      std::string process; // creating process
      std::vector<int>      daughter_idx_v; // daughter node indices in node_v
      std::vector<Node_t*>  daughter_v;     // pointer to daughters 
      std::vector< std::vector<int> > pix_vv; // pixels in each plane. pixels stored in (tick,wire) coordinates
      std::vector<float> start;   //< (x,y,z,t) before sce
      std::vector<float> imgpos4; //< (x,y,z,tick) after sce
      std::vector< std::vector<float> > plane_bbox_twHW_vv; /// bounding box for pixels in each plane
      int origin;

      Node_t()
      : nodeidx(-1),
        type(-1),
        vidx(-1),        
        tid(-1),
        aid(-1),
	mtid(-1),
        pid(-1),
        mother(nullptr),
        mid(-1),
        E_MeV(-1.0),
	process("null"),
        start({0,0,0}),
        imgpos4({0,0,0,0}),	
        origin(-1)
      {};
        
      Node_t(int _nodeidx, int _type, int _tid, int _vidx,
	     int _pid,
	     Node_t* _mother=nullptr,
	     int _mid=-1,
	     float _energy=-1.0,
	     std::string proc="null")
      : nodeidx(_nodeidx),
        type(_type),
        vidx(_vidx),        
        tid(_tid),
	aid(-1),
	mtid(-1),
        pid(_pid),
        mother(_mother),
        mid(_mid),
        E_MeV(_energy),
	process(proc),
        start({0,0,0,0}),
        imgpos4({0,0,0,0}),
        origin(-1)
      {};

      bool operator<( const Node_t& rhs ) const {
        if ( tid < rhs.tid ) return true;
        return false;
      };
    };

    // list of nodes
    std::vector< Node_t > node_v;
    std::vector< std::vector<int> > _unassigned_pixels_vv;

    // number of planes
    size_t _nplanes; // set when _scanPixelData is run

    // search methods
    Node_t* findTrackID( int trackid );

    // print info methods
    void printAllNodeInfo();
    void printNodeInfo( const Node_t& node );
    std::string strNodeInfo( const Node_t& node );
    void printGraph( Node_t* start_node=nullptr, bool visible_only=true );

    // get pixels
    std::vector< std::vector<int> > getPixelsFromParticleAndDaughters( int trackid );

    // graph traversal
    std::vector<Node_t*> getNodeAndDescendentsFromTrackID( const int& trackid );
    void recursiveGetNodeAndDescendents( Node_t* node, std::vector<Node_t*>& nodelist );

    // get primary list
    std::vector<Node_t*> getPrimaryParticles( bool exclude_neutrons=true );
    std::vector<Node_t*> getNeutrinoPrimaryParticles( bool exclude_neutrons=true );    
    std::vector<Node_t*> getNeutrinoParticles( bool exclude_neutrons=true );
    
  protected:
    
    void _recursivePrintGraph( Node_t* node, int& depth, bool visible_only=true );
    void _scanPixelData( const std::vector<larcv::Image2D>& adc_v,
                         const std::vector<larcv::Image2D>& segment_v,
                         const std::vector<larcv::Image2D>& instance_v,
                         const std::vector<larcv::Image2D>& ancestor_v,
                         const std::vector<float> threshold_v );
    void _get_imgpos( std::vector<float>& realpos4,
                      std::vector<float>& imgpos4,
                      larutil::SpaceChargeMicroBooNE& sce,
		      bool apply_sce=true );
    
  public:
    
    std::map<int,int> _shower_daughter2mother;
    void _fill_shower_daughter2mother_map( const std::vector<larlite::mcshower>& mcsh_v );
      
  public:

    // configuration parameters
    // name of trees
    std::string adc_tree; ///< default is wire
    void set_adc_treename( std::string name ) { adc_tree = name; };
    
  };
  
}
}

#endif
