#ifndef __MCPIXELPGRAPH_H__
#define __MCPIXELPGRAPH_H__

#include <vector>

// larcv
#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/EventImage2D.h"

// larlite
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/DataFormat/mcshower.h"
#include "larlite/DataFormat/mctrack.h"
#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"

// ROOT
#include "TH2D.h"
#include "TText.h"
#include "TMarker.h"
#include "TCanvas.h"

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
	node_v(),
	_unassigned_pixels_vv(),
	_eventRootNode(nullptr),
	_nplanes(3),
	photon_start_edep_radius_cm(10.0),
	photon_start_pixval_threshold(30.0),
	photon_start_min_cluster_size(10),
	ev_mctrack(nullptr),
	ev_mcshower(nullptr),
	ev_mctruth(nullptr),
	ev_adc(nullptr),
	ev_seg(nullptr),
	ev_ins(nullptr),
	ev_anc(nullptr),
	_kNuVertexDistCutoff_cm(0.3),
	_nu_vertices_v(),
	_cluster_neutrino_particles(false),
	cvis(nullptr),
	_shower_daughter2mother(),
	adc_tree("wire")
      {};
    virtual ~MCPixelPGraph() {};

    void buildgraph( larcv::IOManager& iolcv, larlite::storage_manager& ioll );
    
    void buildgraph( const std::vector<larcv::Image2D>& adc_v,
                     const std::vector<larcv::Image2D>& segment_v,
                     const std::vector<larcv::Image2D>& instance_v,
                     const std::vector<larcv::Image2D>& ancestor_v,
		     const std::vector<larcv::Image2D>& larflow_v,
                     const larlite::event_mcshower& shower_v,
                     const larlite::event_mctrack&  track_v,
                     const larlite::event_mctruth&  mctruth_v );

    void buildgraphonly( larlite::storage_manager& ioll );    
    void buildgraphonly( const larlite::event_mcshower& shower_v,
                         const larlite::event_mctrack&  track_v,
                         const larlite::event_mctruth&  mctruth_v );

    void set_cluster_neutrino_particles( bool doit ) { _cluster_neutrino_particles=doit; };
    


    struct Node_t {
      int nodeidx;    // position in node_v
      int type;       // track=0, shower=1, nu-vertex=2, genie_fs=3
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
      std::vector< std::vector<float> > pixval_vv; // pixel values in each plane, aligned with pix_vv 
      std::vector<float> pixsum_v; // sum of pixel values in each plane
      std::vector<float> start;   //< (x,y,z,t) before sce, true start of particle
      std::vector<float> first_edep_pos; //< (x,y,z,t) before sce, first step that leaves edep in cryostat
      std::vector<float> first_tpc_pos;  //< (x,y,z,t) before sce, first step inside the TPC, visible in the image
      std::vector<float> first_img_pos;  //< (x,y,z,t) before sce, first step inside the TPC, visible in the image      
      std::vector<float> imgpos4; //< (x,y,z,tick) after sce // the image position corresponding to first_tpc_pos
      std::vector<float> imgpos4_edep; //< (x,y,z,tick) after sce // the image position corresponding to the first_edep_pos
      std::vector<float> imgpos4_start; //< (x,y,z,tick) after sce // the image position corresponding to the start pos after SCE
      std::vector< std::vector<float> > plane_bbox_twHW_vv; /// bounding box for pixels in each plane
    
      int origin; // 1=neutrino, 2=cosmic, 0=unassigned, -1=unassigned

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
        start({0,0,0,0}),
        first_edep_pos({0,0,0,0}),	
        first_tpc_pos({0,0,0,0}),
        first_img_pos({0,0,0,0}),		
        imgpos4({0,0,0,0}),
        imgpos4_edep({0,0,0,0}),
        imgpos4_start({0,0,0,0}),
        origin(-1)
      {
	      daughter_v.clear();
      };
        
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
        first_edep_pos({0,0,0,0}),	
        first_tpc_pos({0,0,0,0}),
        first_img_pos({0,0,0,0}),	
        imgpos4({0,0,0,0}),
        imgpos4_edep({0,0,0,0}),
        imgpos4_start({0,0,0,0}),		
        origin(-1)
      {};

      bool operator<( const Node_t& rhs ) const {
        if ( tid < rhs.tid ) return true;
        return false;
      };

      bool isTrackObject() const {
	      if ( type==0 )
	        return true;
	      return false;
      };

      bool isShowerObject() const {
	      if (type==1)
	        return true;
	      return false;
      };

      bool isNuVertexObject() const {
	      if (type==2)
	        return true;
	      return false;
      };

      bool isGenieFinalStateObject() const {
	      if (type==3)
	        return true;
	      return false;
      };            

    };

    // list of nodes
    std::vector< Node_t > node_v;
    std::vector< std::vector<int> > _unassigned_pixels_vv;
    Node_t* _eventRootNode;

    // number of planes
    size_t _nplanes; // set when _scanPixelData is run

    // search methods
    Node_t* findTrackID( int trackid );
    int getAncestorID( int trackid );
    int getShowerMotherID( int trackid );
    int getParticleID( int trackid );

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

    // get positions
    std::vector<float> getParticleEDepPos( const int& trackid );

    // make visualization of nu particles
    // std::vector< TH2D > makeTH2D( std::string hist_stem_name );

    std::vector< std::vector<float> >
    makeNode3DpointsFromLArFlowTruth( Node_t& node,
				      const std::vector<larcv::Image2D>& adc_v,						   
				      const std::vector<larcv::Image2D>& larflow_v );

    // clear the state
    void clear();

  public:

    // variables and functions for making better quantities to
    // characterize true photons
    typedef std::vector< std::vector<float> > pointList; ///< just a redefinition for convenience
    typedef std::set< std::pair<int,int> > PixelSet_t;   ///< just a redefinition for convenience
    typedef std::vector< larcv::Image2D > ImageSet_t;    ///< just a redefinition for convenience
    
  protected:
    float photon_start_edep_radius_cm;
    float photon_start_pixval_threshold;
    int   photon_start_min_cluster_size;
    std::vector< pointList > _true_photon_v; /// list of 3d points representing the shower trunk
    std::vector< std::vector<larcv::Image2D> > _true_photon_plane_trunkimg_vv;
    std::vector< std::vector< PixelSet_t > >   _true_photon_plane_pixset_vv;
    std::vector< std::vector< float > >        _true_photon_plane_pixsum_vv;
    std::map< int, int > _nodeidx_to_photonlist_index_v; /// map from node.nodeidx to index in _true_photon_v
    pointList _empty_photon_pointlist_v;
    std::vector<float> _empty_pixsum_v;
    std::vector< MCPixelPGraph::PixelSet_t > _empty_pixelset_v;
    std::vector< larcv::Image2D > _empty_imagemask_v;
    
    std::vector<float> fixingPhotonStartPoints( Node_t& node,
						const std::vector<larcv::Image2D>& instance_v,
						const std::vector<larcv::Image2D>& ancestor_v,
						const std::vector<larcv::Image2D>& adc_v,
						const std::vector<larcv::Image2D>& larflow_v );
    std::vector<float> getPlanePixelSumsFromPointList( const std::vector< std::vector<float> >& pointlist,
						       const std::vector<larcv::Image2D>& adc_v );
    std::vector< MCPixelPGraph::PixelSet_t > _getPlanePixelSetsAndPixelSums( const MCPixelPGraph::pointList& pt_v,
									     const std::vector<larcv::Image2D>& adc_v,
									     std::vector<float>& pixsum_v );

    //* Functions to get info about the true photon trunk energy deposits and pixels */
  public:
    const pointList&   getTruePhotonTrunk3DPoints( Node_t& node );
    const pointList&   getTruePhotonTrunk3DPoints( int trackid );    
    std::vector<float> getTruePhotonTrunkPlanePixelSums( int trackid );
    const std::vector< MCPixelPGraph::PixelSet_t >& getTruePhotonTrunkPlanePixelSets( int trackid );
    const std::vector< larcv::Image2D >& getTruePhotonTrunkPlaneImage2DMasks( int trackid );

  protected:

    larlite::event_mctrack*  ev_mctrack;
    larlite::event_mcshower* ev_mcshower;
    larlite::event_mctruth*  ev_mctruth;
    larcv::EventImage2D*    ev_adc;
    larcv::EventImage2D*    ev_seg;
    larcv::EventImage2D*    ev_ins;
    larcv::EventImage2D*    ev_anc;
    larcv::EventImage2D*    ev_larflow;
    
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

    void _adoptNeutrinoOrphans( const larlite::event_mctruth* ev_mctruth );    

    virtual const larlite::mctrack&  _retrieve_mctrackobject( const Node_t* node,  const larlite::event_mctrack&  ev_track_v );
    virtual const larlite::mcshower& _retrieve_mcshowerobject( const Node_t* node, const larlite::event_mcshower& ev_shower_v );    
    virtual const larlite::mctrack&  _retrieve_mctrackobject( const Node_t* node,  larlite::storage_manager& ioll, std::string producername="mcreco" );
    virtual const larlite::mcshower& _retrieve_mcshowerobject( const Node_t* node, larlite::storage_manager& ioll, std::string producername="mcreco" );

    float _kNuVertexDistCutoff_cm;
    std::vector< std::vector<float> > _nu_vertices_v;

    //std::map< int, int > _map_trackid_to_nu_ancestor_v;
    bool _cluster_neutrino_particles;
    int _define_neutrino_interaction_nodes( larlite::storage_manager& ioll );
    int _define_neutrino_interaction_nodes( const larlite::event_mctrack& ev_track_v, const larlite::event_mcshower& ev_shower_v );

    // visualization
    TCanvas* cvis;
    std::vector< TMarker* > vis_markers_v;
    std::vector< TText* >   vis_label_v;
    
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
