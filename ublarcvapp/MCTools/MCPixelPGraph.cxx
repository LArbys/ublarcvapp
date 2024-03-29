#include "MCPixelPGraph.h"

#include <sstream>

// larcv
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/DataFormatTypes.h"

// larlite
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctruth.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"

#include "crossingPointsAnaMethods.h"

namespace ublarcvapp {
namespace mctools {

  void MCPixelPGraph::buildgraph( larcv::IOManager& iolcv,
                                  larlite::storage_manager& ioll ) {
    
    larcv::EventImage2D* ev_adc = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, adc_tree );
    larcv::EventImage2D* ev_seg = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, "segment" );
    larcv::EventImage2D* ev_ins = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, "instance" );
    larcv::EventImage2D* ev_anc = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, "ancestor" );

    if ( ev_adc->Image2DArray().size()==0 ) {
      throw std::runtime_error("No ADC images!");
    }
    if ( ev_seg->Image2DArray().size()==0 ) {
      throw std::runtime_error("No segment images!");
    }
    if ( ev_ins->Image2DArray().size()==0 ) {
      throw std::runtime_error("No instance images!");
    }
    if ( ev_anc->Image2DArray().size()==0 ) {
      throw std::runtime_error("No ancestor images!");
    }
    
    
    larlite::event_mctrack*  ev_mctrack  = (larlite::event_mctrack*) ioll.get_data( larlite::data::kMCTrack,  "mcreco" );
    larlite::event_mcshower* ev_mcshower = (larlite::event_mcshower*)ioll.get_data( larlite::data::kMCShower, "mcreco" );
    larlite::event_mctruth*  ev_mctruth  = (larlite::event_mctruth*) ioll.get_data( larlite::data::kMCTruth,  "generator" );

    buildgraph( ev_adc->Image2DArray(),
                ev_seg->Image2DArray(),
                ev_ins->Image2DArray(),
                ev_anc->Image2DArray(),
                *ev_mcshower, *ev_mctrack, *ev_mctruth );
  }

  /**
   * @brief build only the particle graph (no pixel scan)
   *
   */
  void MCPixelPGraph::buildgraphonly( larlite::storage_manager& ioll )
  {
    larlite::event_mctrack*  ev_mctrack  = (larlite::event_mctrack*) ioll.get_data( larlite::data::kMCTrack,  "mcreco" );
    larlite::event_mcshower* ev_mcshower = (larlite::event_mcshower*)ioll.get_data( larlite::data::kMCShower, "mcreco" );
    larlite::event_mctruth*  ev_mctruth  = (larlite::event_mctruth*) ioll.get_data( larlite::data::kMCTruth,  "generator" );

    buildgraphonly( *ev_mcshower, *ev_mctrack, *ev_mctruth );
  }

  /**
   * @brief build the truth graph and assign pixels
   *
   */
  void MCPixelPGraph::buildgraph( const std::vector<larcv::Image2D>& adc_v,
                                  const std::vector<larcv::Image2D>& segment_v,
                                  const std::vector<larcv::Image2D>& instance_v,
                                  const std::vector<larcv::Image2D>& ancestor_v,
                                  const larlite::event_mcshower& shower_v,
                                  const larlite::event_mctrack&  track_v,
                                  const larlite::event_mctruth&  mctruth_v )
  {

    buildgraphonly( shower_v, track_v, mctruth_v );

    // fill the daugher to mother shower ID map
    _fill_shower_daughter2mother_map( shower_v );    
    
    std::vector<float> threshold_v(adc_v.size(),10.0);
    _scanPixelData( adc_v, segment_v, instance_v, ancestor_v, threshold_v );
  }
  
  /**
   * @brief build the particle graph (no pixel assignments)
   *
   */
  void MCPixelPGraph::buildgraphonly( const larlite::event_mcshower& shower_v,
                                      const larlite::event_mctrack&  track_v,
                                      const larlite::event_mctruth&  mctruth_v ) {

    // how do we build this graph?
    // we want to order N    

    // (0) create root node
    // (1) loop through track and shower, creating node objects
    // (2) loop through node objects connecting daughters to mothers
    // (3) (optional) get depth of each node by doing breath-first traversal
    // (4) sort vector pointers by depth (necessary?)

    node_v.clear();
    node_v.reserve( shower_v.size()+track_v.size() );

    // Create ROOT node
    Node_t neutrino ( node_v.size(), -1, 0, 0, -1 );

    // if there is a neutrino, then we add the start position
    if ( mctruth_v.size()>0 ) {
      const larlite::mctruth& mct = mctruth_v.front();
      neutrino.E_MeV = mct.GetNeutrino().Nu().Trajectory().front().E()*1000.0;
      neutrino.start.resize(4);
      neutrino.start[0] = mct.GetNeutrino().Nu().Trajectory().front().X();
      neutrino.start[1] = mct.GetNeutrino().Nu().Trajectory().front().Y();
      neutrino.start[2] = mct.GetNeutrino().Nu().Trajectory().front().Z();
      neutrino.start[3] = mct.GetNeutrino().Nu().Trajectory().front().T();
    }
    
    node_v.emplace_back( std::move(neutrino) );

    // load spacechargemicroboone
    larutil::SpaceChargeMicroBooNE sce;

    for (int vidx=0; vidx<(int)track_v.size(); vidx++ ) {
      const larlite::mctrack& mct = track_v[vidx];
      LARCV_DEBUG() << "track[" << vidx << "] origin=" << mct.Origin()
                    << " tid=" << mct.TrackID()
                    << " mid=" << mct.MotherTrackID()
                    << " aid=" << mct.AncestorTrackID()
                    << " pid=" << mct.PdgCode()
                    << std::endl;

      // toss out neutrons? (sigh...)
      
      //if ( mct.Origin()==1 ) {
      // neutrino origin

      Node_t tracknode( node_v.size(), 0, mct.TrackID(), vidx, mct.PdgCode() );
      tracknode.E_MeV = mct.Start().E();
      tracknode.aid  = mct.AncestorTrackID();
      tracknode.mtid = mct.MotherTrackID();      
      if ( mct.PdgCode()==2212 ) tracknode.E_MeV -= 938.0;
      else if ( mct.PdgCode()==2112 ) tracknode.E_MeV -= 940.0;
      else if ( abs(mct.PdgCode())==13 )   tracknode.E_MeV -= 105.;
      else if ( abs(mct.PdgCode())==211 )  tracknode.E_MeV -= 135.;
      tracknode.origin = mct.Origin();

      // real position, time
      tracknode.start.resize(4);
      tracknode.start[0] = mct.Start().X();
      tracknode.start[1] = mct.Start().Y();
      tracknode.start[2] = mct.Start().Z();
      tracknode.start[3] = mct.Start().T();
      _get_imgpos( tracknode.start, tracknode.imgpos4, sce );
      
      
      node_v.emplace_back( std::move(tracknode) );
    }

    for (int vidx=0; vidx<(int)shower_v.size(); vidx++ ) {
      const larlite::mcshower& mcsh = shower_v[vidx];

      LARCV_DEBUG() << "shower[" << vidx << "] origin=" << mcsh.Origin()
                    << " tid=" << mcsh.TrackID()
                    << " mid=" << mcsh.MotherTrackID()
                    << " aid=" << mcsh.AncestorTrackID()
                    << " pid=" << mcsh.PdgCode()
                    << std::endl;
      
      //if ( mcsh.Origin()==1 ) {
      // neutrino origin

      Node_t showernode( node_v.size(), 1, mcsh.TrackID(), vidx, mcsh.PdgCode() );
      showernode.E_MeV = mcsh.Start().E();
      showernode.origin = mcsh.Origin();      
      showernode.start.resize(4);
      showernode.aid  = mcsh.AncestorTrackID();
      showernode.mtid = mcsh.MotherTrackID();      
      // showernode.start[0] = mcsh.DetProfile().X();
      // showernode.start[1] = mcsh.DetProfile().Y();
      // showernode.start[2] = mcsh.DetProfile().Z();
      // showernode.start[3] = mcsh.DetProfile().T(); 
      showernode.start[0] = mcsh.Start().X();
      showernode.start[1] = mcsh.Start().Y();
      showernode.start[2] = mcsh.Start().Z();
      showernode.start[3] = mcsh.Start().T();
      _get_imgpos( showernode.start, showernode.imgpos4, sce );
      //showernode.imgpos4 = showernode.start;
      //showernode.imgpos4[3] = 3200 + mcsh.DetProfile().X()/larutil::LArProperties::GetME()->DriftVelocity()/0.5;
      
      node_v.emplace_back( std::move(showernode) );
    }

    // sort Node_t object by geant track ID, relabel node IDs
    std::sort( node_v.begin(), node_v.end() );
    for ( size_t nid=0; nid<node_v.size(); nid++ ) {
      node_v[nid].nodeidx = nid;
    }

    // now connect mothers and daughters. we use the mother and ancestor ID to do this.
    for ( auto& node : node_v ) {
      if ( node.tid<0 ) continue; // the root node

      // find the mother node
      Node_t* mothernode = nullptr;
      
      if ( node.type==0 ) {
        // track nodes
        const larlite::mctrack& track = track_v[ node.vidx ];
        if ( track.TrackID()==track.MotherTrackID() ) {
          // primary
          mothernode = &node_v[0];
        }
        else {
          // secondary
          mothernode = findTrackID( track.MotherTrackID() );          
          if( mothernode==nullptr || mothernode->tid!=(int)track.MotherTrackID() ) {
            // try ancestor ID
            mothernode = findTrackID( track.AncestorTrackID() );
            if ( mothernode && mothernode->tid!=(int)track.AncestorTrackID() )
              mothernode = nullptr;
          }
        }
      }
      else if (node.type==1) {
        //shower nodes
        const larlite::mcshower& shower = shower_v[node.vidx];
        if (shower.TrackID()==shower.MotherTrackID() ) {
          //primary
          mothernode = &node_v[0];
        }
        else {
          //secondary
          mothernode = findTrackID( shower.MotherTrackID() );
          if( mothernode==nullptr || mothernode->tid!=(int)shower.MotherTrackID() ) {
            // try ancestor ID
            mothernode = findTrackID( shower.AncestorTrackID() );
            if ( mothernode && mothernode->tid!=(int)shower.AncestorTrackID() )
              mothernode = nullptr;
          }
        }
      }

      if (mothernode) {
        // found mother, connect
        //std::cout << "found mother: " << strNodeInfo(*mothernode) << std::endl;
        node.mother = mothernode;
        node.mid    = mothernode->nodeidx;
        mothernode->daughter_v.push_back( &node );
        mothernode->daughter_idx_v.push_back( node.nodeidx );
      }
      
    }//end of node loop

    
    //printAllNodeInfo();
    //printGraph();
  }

  /**
   * locate Node_t in node_v using trackid (from geant4)
   *
   * @return The node if found, nullptr if not found
   *
   */
  MCPixelPGraph::Node_t* MCPixelPGraph::findTrackID( int trackid ) {
    Node_t dummy;
    dummy.tid = trackid;
    auto it = std::lower_bound( node_v.begin(), node_v.end(), dummy );
    if ( it==node_v.end() || it->tid!=dummy.tid ) {
      // no node, check the daughter ID map
      auto it_showerdaughter = _shower_daughter2mother.find( trackid );
      if ( it_showerdaughter!=_shower_daughter2mother.end() ) {
        // found an id
        LARCV_DEBUG() << "  found map to mother: " << it_showerdaughter->second << std::endl;
        dummy.tid = it_showerdaughter->second;
      }
      else {
        // still nope
        return nullptr;
      }
      // try again
      it = std::lower_bound( node_v.begin(), node_v.end(), dummy );
      if ( it!=node_v.end() )
        LARCV_DEBUG() << "  mother id maps to existing node" << std::endl;
    }
    
    if ( it==node_v.end() || it->tid!=dummy.tid ) { 
      return nullptr;
    }
    //std::cout << "find trackid=" << trackid << ": " << strNodeInfo( *it ) << std::endl;    
    return &*(it+0);
  }

  /**
   * print node info for all nodes stored in node_v
   *
   */
  void MCPixelPGraph::printAllNodeInfo()  {
    for ( auto const& node : node_v ) {
      printNodeInfo(node);
    }
  }

  /**
   * create string with info from a given Node_t
   *
   * @param[in] node Note_t object to make info for.
   *
   */
  std::string MCPixelPGraph::strNodeInfo( const Node_t& node ) {
    std::stringstream ss;
    //ss << "node[" << node.nodeidx << "," << &node << "] "
    ss << "node[" << node.nodeidx << "] "
       << " (type,vidx)=(" << node.type << "," << node.vidx << ") "
       << " origin=" << node.origin
       << " tid=" << node.tid
       << " mid=" << node.mtid      
       << " aid=" << node.aid
       << " pdg=" << node.pid
       << " KE=" << node.E_MeV << " MeV"
      //<< " start=(" << node.start[0] << "," << node.start[1] << "," << node.start[2] << "," << node.start[3] << ")"
       << " imgpos=(" << node.imgpos4[0] << "," << node.imgpos4[1] << "," << node.imgpos4[2] << "," << node.imgpos4[3] << ")"
      //<< " (mid,mother)=(" << node.mid << "," << node.mother << ") "
      //<< " (mid,mother)=(" << node.mid << ") "
       << " ndaughters=" << node.daughter_v.size()
       << " npixs=(";
    for ( size_t i=0; i<node.pix_vv.size(); i++ ) {
      ss << node.pix_vv[i].size()/2;
      if ( i+1<node.pix_vv.size() ) ss << ", ";
    }
    ss << ")";

    return ss.str();
  }

  /**
   * print Node_t info to standard out
   *
   * @param[in] node Node_t object to print info for.
   *
   */
  void MCPixelPGraph::printNodeInfo( const Node_t& node ) {
    std::cout << strNodeInfo(node) << std::endl;
  }

  /**
   * dump graph to standard out
   *
   */
  void MCPixelPGraph::printGraph( Node_t* rootnode, bool visible_only ) {
    //  here we go!
    std::cout << "=======[ MCPixelPGraph::printGraph ]==============" << std::endl;
    int depth = 0;
    if (rootnode==nullptr )
      rootnode = &node_v.front();
    _recursivePrintGraph( rootnode, depth, visible_only );
  }

  /*
   * internal recursive function that prints node info
   *
   */
  void MCPixelPGraph::_recursivePrintGraph( Node_t* node, int& depth, bool visible_only ) {
    if ( depth<0 ) return; // we're done (error?)
    
    // depth first printing of nodes   
    std::string info = strNodeInfo( *node );
    std::string branch = "";
    for ( int i=0; i<depth; i++ )
      branch += " |";
    if ( depth>0 ) 
      branch += "-- ";
    if ( visible_only ) {
      int nvis = 0;
      for ( auto const& pix_v : node->pix_vv )
        nvis += pix_v.size();
      if ( nvis>0 )
        std::cout << branch << info << std::endl;
    }
    else {
      std::cout << branch << info << std::endl;
    }

    // we loop through our daughters
    ++depth;    
    for ( auto& daughter : node->daughter_v ) {
      _recursivePrintGraph( daughter, depth, visible_only );
    }
    --depth;
    return;
  }

  /**
   * internal method to scan the adc and truth images and fill pixel locations in the Node_t objects
   *
   */
  void MCPixelPGraph::_scanPixelData( const std::vector<larcv::Image2D>& adc_v,
                                      const std::vector<larcv::Image2D>& segment_v,
                                      const std::vector<larcv::Image2D>& instance_v,
                                      const std::vector<larcv::Image2D>& ancestor_v,
                                      const std::vector<float> threshold_v ) {

    _nplanes = adc_v.size();

    // need to check that images have same meta (they should!)
    
    // loop through nodes and setup pixel arrays    
    for (auto& node : node_v ) {
      node.pix_vv.resize(_nplanes);
      for ( size_t p=0; p<_nplanes; p++ ) {
        node.pix_vv[p].clear();
      }
    }

    // track how efficient we were in assigning owner to pixel
    std::vector<int> nabove_thresh(_nplanes,0); // all pixels (some will be cosmic with no truth of course)
    std::vector<int> nabove_thresh_withlabel(_nplanes,0); // pixels with labels
    std::vector<int> nassigned(_nplanes,0);               // assignable to node in node_v
    _unassigned_pixels_vv.clear();
    _unassigned_pixels_vv.resize(_nplanes);

    std::set<int> shower_ancestor_ids;
    
    // loop through images, store pixels into graph nodes
    for ( size_t p=0; p<_nplanes; p++ ) {

      _unassigned_pixels_vv[p].clear();

      auto const& meta = adc_v[p].meta();
      const float threshold = threshold_v[p];
      
      for (size_t r=0; r<meta.rows(); r++) {
        int tick = (int)meta.pos_y(r);
        for (size_t c=0; c<meta.cols(); c++ ) {
          int wire = (int)meta.pos_x(c);

          float adc = 0.0;
          try {
            adc = adc_v[p].pixel(r,c,__FILE__,__LINE__);
          }
          catch (...) {
            std::cerr << __FILE__ << ":L" << __LINE__ << " error getting ADC pixel (" << r << "," << c << ")" << std::endl;
            continue;
          }
          if ( adc<threshold )
            continue;

          nabove_thresh[p]++;
          
          // above threshold, now lets find instance or ancestor
          int tid = 0;
          try {            
            tid = instance_v[p].pixel(r,c,__FILE__,__LINE__);
          }
          catch (...) {
            std::cerr << __FILE__ << ":L" << __LINE__ << " error getting instance pixel (" << r << "," << c << ")" << std::endl;
            continue;            
          }
          int aid = 0;
          try {
            aid = ancestor_v[p].pixel(r,c,__FILE__,__LINE__);
          }
          catch (...) {
            std::cerr << __FILE__ << ":L" << __LINE__ << " error getting ancestor pixel (" << r << "," << c << ")" << std::endl;
            continue;                        
          }

          int seg = 0; 
          try {
            seg = segment_v[p].pixel(r,c,__FILE__,__LINE__);
          }
          catch (...) {
            std::cerr << __FILE__ << ":L" << __LINE__ << " error getting segment pixel (" << r << "," << c << ")" << std::endl;
            continue;                                    
          }

          if ( tid<0 && (seg==(int)larcv::kROIEminus || seg==(int)larcv::kROIGamma) )
            tid *= -1;

          if ( tid>0 || aid>0 )
            nabove_thresh_withlabel[p]++;

          if ( seg==(int)larcv::kROIEminus || seg==(int)larcv::kROIGamma ) {
            shower_ancestor_ids.insert( aid );
          }

          Node_t* node = nullptr;

          if ( tid>0 ) {
            // first we use the instance ID          
            node = findTrackID( tid );
            if ( node==nullptr && adc>10.0 )
              LARCV_DEBUG() << "  no node for charge-pixel tid=" << tid << std::endl;
            // if ( node && node->tid!=tid )
            //   node = nullptr; // reset (what's this?)
          }
          
          // use ancestor if we could not find the node
          if ( !node && aid>0 ) {
            node = findTrackID( aid );
            if ( node && node->tid!=aid )
              node = nullptr;
          }

          if ( node ) {

            // if ( node->tid!=tid && node->aid!=aid ) {
            //   std::cout << "pixel assigned without matching tid or aid exactly: "
            //             << " pix-tid=" << tid << " pix-aid=" << aid
            //             << " node-tid=" << node->tid << " node-mid=" << node->mtid << " node-aid=" << node->aid
            //             << std::endl;
            // }
            
            nassigned[p]++;
            node->pix_vv[p].push_back( tick );
            node->pix_vv[p].push_back( wire );
          }
          else {
            _unassigned_pixels_vv[p].push_back( tick );
            _unassigned_pixels_vv[p].push_back( wire );
          }
          
        }//end of loop over columns
      }//end of loop over rows
      
    }//end of loop over planes

    // no make bounding boxes
    for ( auto& node : node_v ) {
      node.plane_bbox_twHW_vv.clear();
      node.plane_bbox_twHW_vv.resize(_nplanes);      
      for (size_t p=0; p<_nplanes; p++ ) {
        node.plane_bbox_twHW_vv[p].resize( 4, 0 );

        int npix = node.pix_vv[p].size()/2;
        if ( npix==0 )
          continue;
        
        float minx = 1e9;
        float maxx = 0;
        float miny = 1e9;
        float maxy = 0.;

        for ( int ipix=0; ipix<npix; ipix++) {
          float wire = node.pix_vv[p][2*ipix+1];
          float tick = node.pix_vv[p][2*ipix];
          if ( minx>wire ) minx = wire;
          if ( maxx<wire ) maxx = wire;
          if ( miny>tick ) miny = tick;
          if ( maxy<tick ) maxy = tick;
        }
        node.plane_bbox_twHW_vv[p][0] = 0.5*(maxy+miny); // middle tick
        node.plane_bbox_twHW_vv[p][1] = 0.5*(maxx+minx); // middle wire
        node.plane_bbox_twHW_vv[p][2] = 0.5*fabs(maxy-miny); // half-height
        node.plane_bbox_twHW_vv[p][3] = 0.5*fabs(maxx-minx); // half-width
      }// end of loop over planes 
    }//end of loop over nodes

    for (size_t p=0; p<_nplanes; p++ ) {
      std::stringstream msg;      
      msg << " plane[" << p << "]"
          << " num above threshold=" << nabove_thresh[p]
          << " and with label=" << nabove_thresh_withlabel[p]
          << " num assigned=" << nassigned[p]
          << " num unassigned=" << _unassigned_pixels_vv[p].size()/2;
      if ( nabove_thresh_withlabel[p]>0 )
        msg << " fraction=" << float(nassigned[p])/float(nabove_thresh_withlabel[p]);
      LARCV_INFO() << msg.str() << std::endl;
    }
    std::stringstream msg;
    msg << "  ancestor list from all shower pixels: [";
    for ( auto& aid : shower_ancestor_ids )
      msg << aid << " ";
    msg << "]";
    LARCV_INFO() << msg.str() << std::endl;
    
  }

  /**
   * get pixels associated with node and its descendents
   * 
   */
  std::vector< std::vector<int> > MCPixelPGraph::getPixelsFromParticleAndDaughters( int trackid ) {
    std::vector< std::vector<int> > pixels_vv(_nplanes);

    std::vector<MCPixelPGraph::Node_t*> nodelist = getNodeAndDescendentsFromTrackID( trackid );
    for ( auto const& pnode : nodelist ) {
      for ( size_t p=0; p<3; p++ ) {
        if ( pnode->pix_vv[p].size()>0 ) {
          pixels_vv[p].insert( pixels_vv[p].end(), pnode->pix_vv[p].begin(), pnode->pix_vv[p].end() );
        }
      }
    }

    return pixels_vv;    
  }

  /**
   * get list of Nodes_t that are decendents of the given trackID
   *
   *
   */
  std::vector<MCPixelPGraph::Node_t*> MCPixelPGraph::getNodeAndDescendentsFromTrackID( const int& trackid ) {
    std::vector<MCPixelPGraph::Node_t*> nodelist;

    Node_t* rootnode = findTrackID( trackid );
    if ( rootnode==nullptr )
      return nodelist;

    nodelist.push_back( rootnode );
    recursiveGetNodeAndDescendents( rootnode, nodelist );
    return nodelist;
  }

  /**
   * recursively get list of Nodes_t that are descendents of the given Node_t*
   *
   * follows depth-first traversal
   */
  void MCPixelPGraph::recursiveGetNodeAndDescendents( Node_t* node, std::vector<Node_t*>& nodelist ) {
    if ( node==nullptr ) return;
    for ( auto& pdaughter : node->daughter_v ) {
      nodelist.push_back( pdaughter );
      recursiveGetNodeAndDescendents( pdaughter, nodelist );
    }
    return;
  }

  /**
   * get list of primary particles
   *
   * by default, neutrons are excluded
   *
   */
  std::vector<MCPixelPGraph::Node_t*> MCPixelPGraph::getPrimaryParticles( bool exclude_neutrons ) {
    std::vector<Node_t*> nodelist;
    Node_t* rootnode = &node_v[0];
    for ( auto& node : node_v ) {
      if ( node.mother==rootnode ) {
        // primary
        if ( !exclude_neutrons || node.pid!=2112 ) {
          nodelist.push_back( &node );
        }
      }
    }
    return nodelist;      
  }

  /**
   * get list of neutrino-only primary particles
   *
   * by default, neutrons are excluded
   *
   */
  std::vector<MCPixelPGraph::Node_t*> MCPixelPGraph::getNeutrinoPrimaryParticles( bool exclude_neutrons ) {
    std::vector<Node_t*> nodelist;
    Node_t* rootnode = &node_v[0];
    for ( auto& node : node_v ) {
      if ( node.mother==rootnode && node.origin==1 ) {
        // primary
        if ( !exclude_neutrons || node.pid!=2112 ) {
          nodelist.push_back( &node );
        }
      }
    }
    return nodelist;      
  }

  /**
   * get list of neutrino-only primary particles
   *
   * by default, neutrons are excluded
   *
   */
  std::vector<MCPixelPGraph::Node_t*> MCPixelPGraph::getNeutrinoParticles( bool exclude_neutrons ) {
    std::vector<Node_t*> nodelist;
    //Node_t* rootnode = &node_v[0];
    for ( auto& node : node_v ) {
      if ( node.origin==1 ) {
        // neutrino particle
        if ( !exclude_neutrons || node.pid!=2112 ) {
          nodelist.push_back( &node );
        }
      }
    }
    return nodelist;      
  }
  
  /**
   * convert real position+time and calculate apparent position
   */
  void MCPixelPGraph::_get_imgpos( std::vector<float>& realpos4,
                                   std::vector<float>& imgpos4,
                                   larutil::SpaceChargeMicroBooNE& sce )
  {

    imgpos4.resize(4,0);
    
    // apparent position according to image
    std::vector<double> dpos(3,0);
    for (int i=0; i<3; i++) {
      dpos[i]   = realpos4[i];
    }
    
    std::vector<double> offset = sce.GetPosOffsets( dpos[0], dpos[1], dpos[2] );
    std::vector<float>  txyz(4,0);    
    dpos[0] = dpos[0] - offset[0] + 0.7;
    dpos[1] = dpos[1] + offset[1];
    dpos[2] = dpos[2] + offset[2];
    for (int i=0; i<3; i++) {
      txyz[1+i] = realpos4[i];
    }
    txyz[0] = realpos4[3];
    float tick = CrossingPointsAnaMethods::getTick( txyz, 4050.0, &sce );
    for (int i=0; i<3; i++) {
      imgpos4[i] = dpos[i];
    }
    imgpos4[3] = tick;

    // now make x an apparent x
    imgpos4[0] = (tick-3200)*0.5*larutil::LArProperties::GetME()->DriftVelocity();
    
  }

  void MCPixelPGraph::_fill_shower_daughter2mother_map( const std::vector<larlite::mcshower>& mcsh_v )
  {
    LARCV_DEBUG() << "daughter2mother fill" << std::endl;
    _shower_daughter2mother.clear();

    
    for (auto const& mcsh : mcsh_v ) {      
      int showerid = mcsh.TrackID();
      std::vector<unsigned int> dlist = mcsh.DaughterTrackID();
      std::sort( dlist.begin(), dlist.end() );
      for (auto const& daughterid : dlist ) {
        _shower_daughter2mother[daughterid]= showerid;
        LARCV_DEBUG() << "  " << daughterid << " -> " << showerid << std::endl;
      }
    }
    LARCV_INFO() << "Num entries in daughter2mother map: " << _shower_daughter2mother.size() << std::endl;
  }
  
}
}
