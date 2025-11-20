#ifndef __UBLARCVAPP_MCTOOLS_MCPARTICLE_GRAPH_H__
#define __UBLARCVAPP_MCTOOLS_MCPARTICLE_GRAPH_H__

/**
 * @ingroup ublarcvapp::mctools
 * @class MCParticleGraph
 * @brief Store simulated particle information in an event in graph form
 * 
 * We build this graph from larlite::mcshower and larlite::mctrack information.
 * These data products come from the larsorft MCReco producer.
 */

#include <vector>
#include <map>

#include "larcv/core/Base/larcv_base.h"
#include "MCPGNode.h"

// forward declarations
namespace larlite {
  class storage_manager;
  class event_mcshower;
  class event_mctrack;
  class event_mctruth;
  class mcshower;
  class mctrack;
}

namespace ublarcvapp {
namespace mctools {

  class MCParticleGraph : public larcv::larcv_base {

  public:

    MCParticleGraph()
    : larcv::larcv_base("MCParticleGraph"),
      _eventRootNode(nullptr),
      _cluster_neutrino_particles(false),
    	_kNuVertexDistCutoff_cm(0.3)
    {};

    ~MCParticleGraph()
    {};

    void clear();
    void cluster_nu_particles( bool doit=true ) { _cluster_neutrino_particles=doit; };
    void buildgraph( larlite::storage_manager& ioll );
    void buildgraph( const larlite::event_mcshower& shower_v,
                     const larlite::event_mctrack&  track_v,
                     const larlite::event_mctruth&  mctruth_v );

    MCPGNode* findTrackID( long trackid );
    std::vector<MCPGNode*> getNodeAndDescendentsFromTrackID( const int& trackid );
    std::vector<MCPGNode*> getPrimaryParticles( bool exclude_neutrons );
    std::vector<MCPGNode*> getNeutrinoPrimaryParticles( bool exclude_neutrons );
    std::vector<MCPGNode*> getNeutrinoParticles( bool exclude_neutrons );

    void printAllNodeInfo();
    std::string strNodeInfo( const MCPGNode& node );
    void printNodeInfo( const MCPGNode& node );
    void printGraph( MCPGNode* rootnode, bool visible_only );
    

  protected:
    void recursiveGetNodeAndDescendents( MCPGNode* node, std::vector<MCPGNode*>& nodelist );
    void _recursivePrintGraph( MCPGNode* node, int& depth, bool visible_only );
    const larlite::mctrack&  _retrieve_mctrackobject(  const MCPGNode* node, const larlite::event_mctrack&  ev_track_v );
    const larlite::mcshower& _retrieve_mcshowerobject( const MCPGNode* node, const larlite::event_mcshower& ev_shower_v );    
    const larlite::mctrack&  _retrieve_mctrackobject(  const MCPGNode* node, larlite::storage_manager& ioll, std::string producername="mcreco" );
    const larlite::mcshower& _retrieve_mcshowerobject( const MCPGNode* node, larlite::storage_manager& ioll, std::string producername="mcreco" );


  public:

    std::vector< MCPGNode > node_v; //< collection of nodes
    MCPGNode* _eventRootNode;
    bool _cluster_neutrino_particles; //< if true, nu primary particles clustered together
    float _kNuVertexDistCutoff_cm;
    std::vector< std::vector<float> > _nu_vertices_v;
    std::map<long,long> _tid_to_node_v; ///< map from trackid to position in node_v

    std::map<int,int> _shower_daughter2mother;

  protected:
    void _fill_shower_daughter2mother_map( const std::vector<larlite::mcshower>& mcsh_v );
    int _define_neutrino_interaction_nodes( larlite::storage_manager& ioll );
    int _define_neutrino_interaction_nodes( const larlite::event_mctrack& ev_mctrack,
							 const larlite::event_mcshower& ev_mcshower );




  };

}
}

#endif