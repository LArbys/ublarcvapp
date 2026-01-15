#include "MCParticleGraph.h"

#include <sstream>
#include <set>

// larlite
#include "larlite/DataFormat/mctrack.h"
#include "larlite/DataFormat/mcshower.h"
#include "larlite/DataFormat/mctruth.h"
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/LArUtil/LArProperties.h"
#include "larlite/LArUtil/Geometry.h"
#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"

// ublarcvapp
#include "ublarcvapp/dbscan/DBScan.h"

#include "crossingPointsAnaMethods.h"
#include "MCPos2ImageUtils.h"

namespace ublarcvapp {
namespace mctools {

  /**
   * @brief clear all data containers
   */
  void MCParticleGraph::clear() 
  {
      node_v.clear();
      _eventRootNode = nullptr;
      _shower_daughter2mother.clear();
      _tid_to_node_v.clear();
  }

  /**
   * @brief build only the particle graph (no pixel scan)
   *
   */
  void MCParticleGraph::buildgraph( larlite::storage_manager& ioll )
  {
    
    auto ev_mctrack  = (larlite::event_mctrack*) ioll.get_data( larlite::data::kMCTrack,  "mcreco" );
    auto ev_mcshower = (larlite::event_mcshower*)ioll.get_data( larlite::data::kMCShower, "mcreco" );
    auto ev_mctruth  = (larlite::event_mctruth*) ioll.get_data( larlite::data::kMCTruth,  "generator" );

    buildgraph( *ev_mcshower, *ev_mctrack, *ev_mctruth );

  }
  
  /**
   * @brief build the particle graph (no pixel assignments)
   *
   */
  void MCParticleGraph::buildgraph( const larlite::event_mcshower& shower_v,
                                    const larlite::event_mctrack&  track_v,
                                    const larlite::event_mctruth&  mctruth_v ) {

    // how do we build this graph?
    // we want to order N    

    // (0) create root node
    // (1) loop through track and shower, creating node objects
    // (2) loop through node objects connecting daughters to mothers
    // (3) (optional) get depth of each node by doing breath-first traversal
    // (4) sort vector pointers by depth (necessary?)

    // dump mctruth info
    // LARCV_DEBUG() << "MCTruth Dump" << std::endl;
    // int imctruth=0;
    // for ( auto& mctruth : mctruth_v ) {
    //   LARCV_DEBUG() << "MCTRUTH[" << imctruth << "] --------------" << std::endl;
    //   int imcpart = 0;
    //   for ( auto& part : mctruth.GetParticles() ) {
    // 	LARCV_DEBUG() << " mcpart[" << imcpart << "] -------" << std::endl;
    // 	LARCV_DEBUG() << "   status=" << part.StatusCode() << std::endl;
    // 	LARCV_DEBUG() << "   trackid=" << part.TrackId() << std::endl;
    // 	LARCV_DEBUG() << "   pdg=" << part.PdgCode() << std::endl;
    // 	LARCV_DEBUG() << "   motherid=" << part.Mother() << std::endl;
    // 	LARCV_DEBUG() << "   process=" << part.Process() << " endprocess=" << part.EndProcess() << std::endl;
    // 	LARCV_DEBUG() << "   num daughters=" << part.Daughters().size() << std::endl;
    //   }
    // }

    clear();
    
    // fill the daugher to mother shower ID map
    _fill_shower_daughter2mother_map( shower_v ); 

    node_v.clear();
    node_v.reserve( shower_v.size()+track_v.size()+100 );

    // Create ROOT node
    MCPGNode neutrino ( node_v.size(), -1, 0, 0, -1 );

    // define TPC widths
    float tpc_x = larutil::Geometry::GetME()->DetHalfWidth()*2.0;
    float tpc_y = larutil::Geometry::GetME()->DetHalfHeight();
    float tpc_z = larutil::Geometry::GetME()->DetLength();

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

    _eventRootNode = &node_v[0];

    // keep tabs on all track IDs found
    std::set<int> tid_list;

    // load spacechargemicroboone: forward effect
    //larutil::SpaceChargeMicroBooNE sce;

    // struct for primary particles from nu interactions
    struct NuPart_t {
      int geantid;
      int pdg;
      float E_MeV;
      std::vector<float> pos;
    };

    std::vector< NuPart_t > nu_part_v;

    // collect data from mcreco tracks
    int mc_product_type = 0; // larlite::mctrack
    for (int vidx=0; vidx<(int)track_v.size(); vidx++ ) {
      const larlite::mctrack& mct = track_v[vidx];
      LARCV_DEBUG() << "track[" << vidx << "] origin=" << mct.Origin()
                    << " tid=" << mct.TrackID()
                    << " mid=" << mct.MotherTrackID()
                    << " aid=" << mct.AncestorTrackID()
                    << " pid=" << mct.PdgCode()
                    << std::endl;

      MCPGNode tracknode( node_v.size(), mc_product_type, mct.TrackID(), vidx, mct.PdgCode() );
      tracknode.E_MeV = mct.Start().E();
      tracknode.aid  = mct.AncestorTrackID();
      tracknode.mtid = mct.MotherTrackID();
      tracknode.process = mct.Process();
      tracknode.mother_process   = mct.MotherProcess();
      tracknode.ancestor_process = mct.AncestorProcess();

      // convert energy into relativistic KE
      if ( mct.PdgCode()==2212 ) tracknode.E_MeV -= 938.0;
      else if ( mct.PdgCode()==2112 ) tracknode.E_MeV -= 940.0;
      else if ( abs(mct.PdgCode())==13 )   tracknode.E_MeV -= 105.;
      else if ( abs(mct.PdgCode())==211 )  tracknode.E_MeV -= 135.;
      tracknode.origin = mct.Origin();

      // exception for process=neutronInelastic
      if ( mct.PdgCode()==2112 && mct.Process()=="neutronInelastic" ) {
          // these are actually protons made by a neutrino inelastic process
          // just the way the MC reco producer saves these.
          tracknode.E_MeV = mct.Start().E()-938.0;
      }

      tid_list.insert( tracknode.tid );
      
      // set creation position, time
      tracknode.start.resize(4);
      tracknode.start[0] = mct.Start().X();
      tracknode.start[1] = mct.Start().Y();
      tracknode.start[2] = mct.Start().Z();
      tracknode.start[3] = mct.Start().T();
      try {
        tracknode.mom4[0]  = mct.Start().E();
        tracknode.mom4[1]  = mct.Start().Px();
        tracknode.mom4[2]  = mct.Start().Py();
        tracknode.mom4[3]  = mct.Start().Pz();
      }
      catch (...) {
        tracknode.mom4 = std::vector<float>{0,0,0,0};
      }

      if ( tracknode.origin==1 ) {
	      // store nu particle
	      NuPart_t nuparticle;
	      nuparticle.geantid = tracknode.tid;
	      nuparticle.pdg = tracknode.pid;
	      nuparticle.E_MeV = mct.Start().E();
	      nuparticle.pos = tracknode.start;
	      nu_part_v.push_back( nuparticle );
      }
      
      // find the first position
      //  - with cryo edep
      //  - inside the TPC
      tracknode.first_edep_pos = std::vector<float>(4,0);
      tracknode.first_tpc_pos = std::vector<float>(4,0);
	
      if ( mct.size()>=1 ) {
        tracknode.first_edep_pos[0] = mct[0].X();
        tracknode.first_edep_pos[1] = mct[0].Y();
        tracknode.first_edep_pos[2] = mct[0].Z();
        tracknode.first_edep_pos[3] = mct[0].T();
      
        // first TPC point
        std::vector<float> recopos(4,0);
        std::vector<float> edep_imgpos(4,0);		
        std::vector<float> xyz(4,0);
        bool intpc = false;
        bool inimage = false;
        bool first_step = false;
        for (auto const& step : mct ) {
          xyz[0] = (float)step.X();
          xyz[1] = (float)step.Y();
          xyz[2] = (float)step.Z();
          xyz[3] = (float)step.T();

          if ( !first_step ) {
            first_step = true;
            tracknode.first_edep_pos = xyz;
          }
          
          if ( ( xyz[0]>0.0 && xyz[0]<tpc_x )
              && ( fabs(xyz[1])<tpc_y )
              && ( xyz[2]>0 && xyz[2]<tpc_z ) ) {

            if ( intpc==false ) {
              intpc = true;
              tracknode.first_tpc_pos = xyz;
            }

          }
          if ( intpc ) {
            // no need to keep searching
            break;
          }
        }//end of mcstep loop
      }//end of if more than 0 mcsteps

      node_v.emplace_back( std::move(tracknode) );
    }

    // collect from mcshowers
    mc_product_type = 1;
    for (int vidx=0; vidx<(int)shower_v.size(); vidx++ ) {
      const larlite::mcshower& mcsh = shower_v[vidx];

      LARCV_DEBUG() << "shower[" << vidx << "] origin=" << mcsh.Origin()
                    << " tid=" << mcsh.TrackID()
                    << " mid=" << mcsh.MotherTrackID()
                    << " aid=" << mcsh.AncestorTrackID()
                    << " pid=" << mcsh.PdgCode()
                    << " start=(" << (float)mcsh.Start().X() << "," 
                    << (float)mcsh.Start().Y() << "," 
                    << (float)mcsh.Start().Z() << ","
                    << " t=" << (float)mcsh.Start().T() << ")"
                    << std::endl;

      MCPGNode showernode( node_v.size(), mc_product_type, mcsh.TrackID(), vidx, mcsh.PdgCode() );
      showernode.E_MeV = mcsh.Start().E();
      showernode.process = mcsh.Process();
      showernode.origin = mcsh.Origin();      
      showernode.start = std::vector<float>{ (float)mcsh.Start().X(), (float)mcsh.Start().Y(), (float)mcsh.Start().Z(), (float)mcsh.Start().T() };
      showernode.aid  = mcsh.AncestorTrackID();
      showernode.mtid = mcsh.MotherTrackID();

      showernode.first_edep_pos = std::vector<float>(4,0);
      showernode.first_tpc_pos  = std::vector<float>(4,0);

      try {
        showernode.mom4[0]  = mcsh.Start().E();
        showernode.mom4[1]  = mcsh.Start().Px();
        showernode.mom4[2]  = mcsh.Start().Py();
        showernode.mom4[3]  = mcsh.Start().Pz();
      }
      catch (...) {
        showernode.mom4 = std::vector<float>{0,0,0,0};
      }

      bool x_isinf = false;
      if ( mcsh.DetProfile().X()>1.0e100 || mcsh.DetProfile().X()<-1.0e100 )
	      x_isinf = true;
      if ( !x_isinf && !std::isnan( mcsh.DetProfile().X() ) ) {
        //std::cout << "inf test: " << mcsh.DetProfile().X() << " " << std::isinf( mcsh.DetProfile().X() ) << std::endl;
        std::vector<float> detprofile = { (float)mcsh.DetProfile().X(), (float)mcsh.DetProfile().Y(), (float)mcsh.DetProfile().Z(), (float)mcsh.DetProfile().T() };
        showernode.first_edep_pos = detprofile;
        showernode.first_tpc_pos  = detprofile;
      }
      if ( showernode.origin==1 ) {
	      // store nu particle
	      NuPart_t nuparticle;
	      nuparticle.geantid = showernode.tid;
	      nuparticle.pdg = showernode.pid;
	      nuparticle.E_MeV = mcsh.Start().E();
	      nuparticle.pos = std::vector<float>{ (float)mcsh.Start().X(), 
                                             (float)mcsh.Start().Y(), 
                                             (float)mcsh.Start().Z(), 
                                             (float)mcsh.Start().T() };
	      nu_part_v.push_back( nuparticle );
      }
      
      tid_list.insert( showernode.tid );      
      node_v.emplace_back( std::move(showernode) );
    }
    
    // find the geant4 trackid offset for the neutrinos
    long smallest_nu_tid = -1;
    for ( auto& node : node_v ) {
      if ( node.origin==1 ) {
	      if (smallest_nu_tid<0 || node.tid < smallest_nu_tid ) {
	        smallest_nu_tid = node.tid;
	      }
      }
    }
    LARCV_DEBUG() << "Smallest neutrino Geant4 Track ID: " << smallest_nu_tid << std::endl;
    
    // collect from mctruth (data form GENIE generator)
    mc_product_type = 3;

    // we need to match the earliest track or shower object to
    // a genie final state particle to get offset
    int ifs_1 = 0;
    int matched_fs = -1;
    int matched_geant_id = -1;
    for ( auto& mctruth : mctruth_v ) {
      for ( auto& part : mctruth.GetParticles() ) {
	      if ( part.StatusCode()==1 ) {
	        ifs_1 += 1;
	        // does ginal state particel match any of the geant4 particles?
	        for (auto& nupart : nu_part_v) {
	          // match pdg
	          if ( nupart.pdg!=part.PdgCode() )
	            continue;

	          // match energy
	          float dE_MeV = std::fabs(nupart.E_MeV-part.Momentum(0)[3]*1000.0);
	          LARCV_DEBUG() << "mctruth part: " << dE_MeV << std::endl;
	          if ( dE_MeV > 10.0 )
	            continue;

	          matched_fs = ifs_1;
	          matched_geant_id = nupart.geantid;
	          LARCV_DEBUG() << "have first match. "
                          << "genie finalstate id=" << matched_fs 
                          << " geantid="  << matched_geant_id 
                          << std::endl;
	            break;
	        }
	  
	      }//end of if status code
	      if (matched_fs>=0)
	        break;
      }//end of mcpart loop
      if (matched_fs>=0)
	      break;
    }//end of mctruth loop

    int geantid_offset = smallest_nu_tid  - (matched_fs-1);
    LARCV_INFO() << "geantid offset: " << geantid_offset << " smallest_nu_tid=" << smallest_nu_tid << std::endl;
    LARCV_INFO() << " matched_fs=" << matched_fs << " matched_geant4=" << matched_geant_id << std::endl;
	
    int imctruth=0;
    int ifs = 0;
    for ( auto& mctruth : mctruth_v ) {
      int imcpart = 0;
      for ( auto& part : mctruth.GetParticles() ) {
        LARCV_DEBUG() << " mcpart[" << imctruth << "," << imcpart << "] -------" << std::endl;
        LARCV_DEBUG() << "   status=" << part.StatusCode() << std::endl;
        LARCV_DEBUG() << "   trackid=" << part.TrackId() << std::endl;
        LARCV_DEBUG() << "   pdg=" << part.PdgCode() << std::endl;
        LARCV_DEBUG() << "   motherid=" << part.Mother() << std::endl;
        LARCV_DEBUG() << "   process=" << part.Process() << " endprocess=" << part.EndProcess() << std::endl;
        LARCV_DEBUG() << "   num daughters=" << part.Daughters().size() << std::endl;
	
        if ( part.StatusCode()==1 ) {
          int geant_trackid = geantid_offset + ifs;
          LARCV_DEBUG() << "  Stable Final State. Implied Geant4 ID = " << geant_trackid << std::endl;
          ifs++;
          
          auto it_tid = tid_list.find(geant_trackid);
          if (it_tid==tid_list.end() ) {
            int nodeidx = node_v.size();
            int type = 3; // genie-final-state
            int tid  = geant_trackid;
            int vidx = imcpart;
            int pid = part.PdgCode();
            MCPGNode* mother = nullptr;
            int mid = -1;
            float energy = part.Momentum(0)[3]*1e3 - part.Mass()*1.0e3; // in GeV, convert to MeV
            std::vector<float> start { (float)part.Position(0)[0],
                                      (float)part.Position(0)[1],
                                      (float)part.Position(0)[2],
                                      (float)part.Position(0)[3]};
            
            //LARCV_DEBUG() << "Creating Mother node from GENIE final states: tid=" << tid << " type=" << 3 << std::endl;
            MCPGNode fsnode( nodeidx, mc_product_type, tid, vidx, pid, mother, mid, energy, "primary" );
            fsnode.start = start;
            // fsnode.first_edep_pos = start;
            // fsnode.first_tpc_pos  = start;
            fsnode.mtid = tid;
            fsnode.aid  = tid;
            fsnode.origin = 1; // neutrino origin (from genie)

            try {
              for (int i=0; i<3; i++)
                fsnode.mom4[i+1]  = part.Momentum(0)[i];
              fsnode.mom4[0]  = part.Momentum(0)[3];
            }
            catch (...) {
              fsnode.mom4 = std::vector<float>{0,0,0,0};
            }

            LARCV_DEBUG() << "Add Genie Final State Particle to initial List: tid=" << tid << " pdg=" << pid << std::endl;	    
            node_v.emplace_back( std::move(fsnode) );
          }
          else {
            LARCV_DEBUG() << "Genie Final State Partial [tid=" << geant_trackid << "] Already in MCReco List" << std::endl;
          }
        }
        imcpart++;	
      }
      imctruth++;
    }
        
    // try to connect primary neutrino origin nodes to mctruth info
    //_adoptNeutrinoOrphans( &mctruth_v ); // whats this?

    // now we connect the graph

    // sort MCPGNode object by geant track ID, relabel node IDs
    std::sort( node_v.begin(), node_v.end() );
    for ( size_t nid=0; nid<node_v.size(); nid++ ) {
      node_v[nid].nodeidx = nid;
    }

    // now connect mothers and daughters. we use the mother and ancestor ID to do this.
    for ( auto& node : node_v ) {
      if ( node.tid<0 ) continue; // the root node

      // find the mother node
      MCPGNode* mothernode = nullptr;
      
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
      else if (node.type==3) {
	      // genie fs nodes
	      // connet to ROOT node as they are primary by definition
	      mothernode = &(node_v[0]);
      }      
      else {
	      continue;
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

    if ( _cluster_neutrino_particles ) {
      int nnu = _define_neutrino_interaction_nodes( track_v, shower_v );
      LARCV_INFO() << "Rearranged graph to include neutrino vertex nodes. Number of Nu Interactions: " << nnu << std::endl;
    }
    
    // map track id to position in node_v container
    for (size_t i=0; i<node_v.size(); i++){
        auto& node = node_v.at(i);
        _tid_to_node_v[ node.tid ] = (long)i;
    }
    //printAllNodeInfo();
    //printGraph();
  }

  /**
   * @brief wrapper function to retrieve Node's corresponding mctrack object from larlite container
   */
  const larlite::mctrack& 
  MCParticleGraph::_retrieve_mctrackobject( const MCPGNode* node,
								            const larlite::event_mctrack& ev_track_v )
  {
    try {
      const larlite::mctrack& x = ev_track_v.at(node->vidx);
      // to do: check that truth info matches?
      return x;
    }
    catch ( std::exception& ex ) {
      std::stringstream err;
      err << "Error accessing event_mctrack container at index=" << node->vidx << " for node[idx]=nodeidx" << std::endl;
      err << "Node info: " << std::endl;
      err << strNodeInfo( *node ) << std::endl;
      err << "Exception: " << ex.what() << std::endl;
      throw std::runtime_error( err.str() );
    }
  }

  /**
   * @brief wrapper function to retrieve Node's corresponding mcshower object from larlite container
   */
  const larlite::mcshower& 
  MCParticleGraph::_retrieve_mcshowerobject( const MCPGNode* node,
								             const larlite::event_mcshower& ev_shower_v )
  {
    try {
      const larlite::mcshower& x = ev_shower_v.at(node->vidx);
      // to do: check that truth info matches?
      return x;
    }
    catch ( std::exception& ex ) {
      std::stringstream err;
      err << "Error accessing event_mcshower container at index=" << node->vidx << " for node[idx]=nodeidx" << std::endl;
      err << "Node info: " << std::endl;
      err << strNodeInfo( *node ) << std::endl;
      err << "Exception: " << ex.what() << std::endl;
      throw std::runtime_error( err.str() );
    }
  }

  /**
   * @brief convenience function to retrieve Node's corresponding mctrack object from larlite top-level io interface
   */
  const larlite::mctrack&  MCParticleGraph::_retrieve_mctrackobject( const MCPGNode* node,
								   larlite::storage_manager& ioll,
								   std::string producername )
  {
    const larlite::event_mctrack* ev_mctrack
      = (larlite::event_mctrack*)ioll.get_data(larlite::data::kMCTrack, producername );

    return _retrieve_mctrackobject( node, *ev_mctrack );
  }

  /**
   * @brief convenience function to retrieve Node's corresponding mcshower object from larlite top-level io interface
   */
  const larlite::mcshower&  
  MCParticleGraph::_retrieve_mcshowerobject( const MCPGNode* node,
								     larlite::storage_manager& ioll,
								     std::string producername )
  {
    const larlite::event_mcshower* ev_mcshower
      = (larlite::event_mcshower*)ioll.get_data(larlite::data::kMCShower, producername );

    return _retrieve_mcshowerobject( node, *ev_mcshower );
  }
  
  
  /**
   * @brief locate MCPGNode in node_v using trackid (from geant4)
   * 
   * we first search our vector of MCPGNode objects, which are sorted by track id.
   *
   * we then search the keys of the _shower_daughter2mother map, where the keys
   *  are track id of particles made by the initial shower.
   *
   * @return The node if found, nullptr if not found
   *
   */
  MCPGNode* MCParticleGraph::findTrackID( long trackid ) {

    auto it = _tid_to_node_v.find( (long)trackid );
    if ( it==_tid_to_node_v.end())
      return nullptr;

    auto& node = node_v.at(it->second);
    return &node;

    // MCPGNode dummy;
    // dummy.tid = trackid;
    // auto it = std::lower_bound( node_v.begin(), node_v.end(), dummy );
    // if ( it==node_v.end() || it->tid!=dummy.tid ) {
    //   // no node, check the daughter ID map
    //   auto it_showerdaughter = _shower_daughter2mother.find( trackid );
    //   if ( it_showerdaughter!=_shower_daughter2mother.end() ) {
    //     // found an id
    //     //LARCV_DEBUG() << "  found map to mother: " << it_showerdaughter->second << std::endl;
    //     dummy.tid = it_showerdaughter->second;
    //   }
    //   else {
    //     // still nope
    //     return nullptr;
    //   }
    //   // with the mother shower's trackid, try to find the node again
    //   it = std::lower_bound( node_v.begin(), node_v.end(), dummy );
    //   //if ( it!=node_v.end() )
    //   //  LARCV_DEBUG() << "  mother id maps to existing node" << std::endl;
    // }
    
    // if ( it==node_v.end() || it->tid!=dummy.tid ) { 
    //   return nullptr;
    // }
    // //std::cout << "find trackid=" << trackid << ": " << strNodeInfo( *it ) << std::endl;    
    // return &*(it+0);
  }

  /**
   * print node info for all nodes stored in node_v
   *
   */
  void MCParticleGraph::printAllNodeInfo()  {
    for ( auto const& node : node_v ) {
      printNodeInfo(node);
    }
  }

  /**
   * create string with info from a given MCPGNode
   *
   * @param[in] node Note_t object to make info for.
   *
   */
  std::string MCParticleGraph::strNodeInfo( const MCPGNode& node ) 
  {

    int hasmother = ( node.mother ) ? 1 : 0;
    
    std::stringstream ss;
    //ss << "node[" << node.nodeidx << "," << &node << "] "
    ss << "node[" << node.nodeidx << "] "
       << " (type,vidx)=(" << node.type << "," << node.vidx << ") "
       << " origin=" << node.origin
       << " pdg=" << node.pid
       << " p=" << node.process;
    if ( node.mother_process!="" )
       ss << "/" << node.mother_process;
    if ( node.ancestor_process!="" )
       ss << "/" << node.ancestor_process;
    ss << " tid=" << node.tid
       << " mtid=" << node.mtid      
       << " aid=" << node.aid
       << " KE=" << node.E_MeV << " MeV"
      //<< " xyzt=
      //<< " tusec=" << node.start[3]*1.0e-3 << " us"
      //<< " (mid,mother)=(" << node.mid << "," << node.mother << ") "
      //<< " (mid,mother)=(" << node.mid << ") "
       << " hasmother=" << hasmother
       << " ndaughters=" << node.daughter_v.size()
       << std::endl;
    // additional positional info
    if ( node.start.size()>=4 )
      ss << "    start(x,y,z,t)=(" << node.start[0] << "," << node.start[1] << "," << node.start[2] << "," << node.start[3]*1.0e-3 << " us) " << std::endl;
    
    if ( node.first_edep_pos.size()>=4 )
      ss << "    edep-cryo(x,y,z,t)=(" << node.first_edep_pos[0] << ","
          << node.first_edep_pos[1] << ","
          << node.first_edep_pos[2] << ","
          << node.first_edep_pos[3]*1.0e-3
          << " us)" << std::endl;

    if ( node.first_tpc_pos.size()>=4 )
      ss << "    edep-tpc(x,y,z,t)=(" << node.first_tpc_pos[0] << ","
          << node.first_tpc_pos[1] << ","
          << node.first_tpc_pos[2] << ","
          << node.first_tpc_pos[3]*1.0e-3
          << " us)" << std::endl;

    if ( node.mom4.size()>=4 )
      ss << "    mom4(E,px,py,pz)=(" << node.mom4[0] << ","
          << node.mom4[1] << ","
          << node.mom4[2] << ","
          << node.mom4[3] << ") "
          << " MeV" << std::endl;
    

    //ss << std::endl;
    return ss.str();
  }

  /**
   * print MCPGNode info to standard out
   *
   * @param[in] node MCPGNode object to print info for.
   *
   */
  void MCParticleGraph::printNodeInfo( const MCPGNode& node ) {
    std::cout << strNodeInfo(node) << std::endl;
  }

  /**
   * dump graph to standard out
   *
   */
  void MCParticleGraph::printGraph( MCPGNode* rootnode, bool visible_only ) {
    //  here we go!
    std::cout << "=======[ MCParticleGraph::printGraph ]==============" << std::endl;
    int depth = 0;
    if (rootnode==nullptr )
      rootnode = &node_v.front();
    _recursivePrintGraph( rootnode, depth, visible_only );
  }

  /*
   * internal recursive function that prints node info
   *
   */
  void MCParticleGraph::_recursivePrintGraph( MCPGNode* node, int& depth, bool visible_only ) {
    if ( depth<0 ) return; // we're done (error?)
    
    // depth first printing of nodes   
    std::string info = strNodeInfo( *node );
    std::string branch  = "";
    std::string branch2 = "";
    for ( int i=0; i<depth; i++ ) {
      branch +=  " |";
      branch2 += " |";
    }
    if ( depth>0 ) {
      branch +=  "-- ";
      branch2 += "   ";
    }

    std::stringstream ss(info);
    std::string to;
    bool firstline = true;
    while(std::getline(ss,to,'\n'))
    {
        if (firstline) 
            std::cout << branch << to << std::endl;
        else
            std::cout << branch2 << to << std::endl;
        firstline = false;
    }
      
    //std::cout << branch << info;// << std::endl;

    // we loop through our daughters
    ++depth;    
    for ( auto& daughter : node->daughter_v ) {
      _recursivePrintGraph( daughter, depth, visible_only );
    }
    --depth;
    return;
  }

  /**
   * get list of Nodes_t that are decendents of the given trackID
   *
   *
   */
  std::vector<MCPGNode*> MCParticleGraph::getNodeAndDescendentsFromTrackID( const int& trackid ) {
    std::vector<MCPGNode*> nodelist;

    MCPGNode* rootnode = findTrackID( trackid );
    if ( rootnode==nullptr )
      return nodelist;

    nodelist.push_back( rootnode );
    recursiveGetNodeAndDescendents( rootnode, nodelist );
    return nodelist;
  }

  /**
   * recursively get list of Nodes_t that are descendents of the given MCPGNode*
   *
   * follows depth-first traversal
   */
  void MCParticleGraph::recursiveGetNodeAndDescendents( MCPGNode* node, std::vector<MCPGNode*>& nodelist ) {
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
  std::vector<MCPGNode*> MCParticleGraph::getPrimaryParticles( bool exclude_neutrons ) {
    std::vector<MCPGNode*> nodelist;
    MCPGNode* rootnode = &node_v[0];
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
  std::vector<MCPGNode*> MCParticleGraph::getNeutrinoPrimaryParticles( bool exclude_neutrons ) {
    std::vector<MCPGNode*> nodelist;
    MCPGNode* rootnode = &node_v[0];
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
  std::vector<MCPGNode*> MCParticleGraph::getNeutrinoParticles( bool exclude_neutrons ) {
    std::vector<MCPGNode*> nodelist;
    //MCPGNode* rootnode = &node_v[0];
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
  
//   /**
//    * convert real position+time and calculate apparent position
//    */
//   void MCParticleGraph::_get_imgpos( std::vector<float>& realpos4,
//                                    std::vector<float>& imgpos4,
//                                    larutil::SpaceChargeMicroBooNE& sce,
// 				   bool apply_sce )
//   {

//     imgpos4.resize(4,0);
    
//     // apparent position according to image
//     std::vector<double> dpos(3,0);
//     for (int i=0; i<3; i++) {
//       dpos[i]   = realpos4[i];
//     }

//     std::vector<float>  txyz(4,0);
//     if ( apply_sce ) {
//       std::vector<double> offset = sce.GetPosOffsets( dpos[0], dpos[1], dpos[2] );        
//       dpos[0] = dpos[0] - offset[0] + 0.7;
//       dpos[1] = dpos[1] + offset[1];
//       dpos[2] = dpos[2] + offset[2];
//     }
    
//     for (int i=0; i<3; i++) {
//       txyz[1+i] = realpos4[i];
//     }
//     txyz[0] = realpos4[3];
//     float tick = 0.;
//     if (apply_sce)
//       tick = CrossingPointsAnaMethods::getTick( txyz, 4050.0, &sce );
//     else
//       tick = CrossingPointsAnaMethods::getTick( txyz, 4050.0, NULL );
    
//     for (int i=0; i<3; i++) {
//       imgpos4[i] = dpos[i];
//     }
//     imgpos4[3] = tick;

//     // now make x an apparent x
//     imgpos4[0] = (tick-3200)*0.5*larutil::LArProperties::GetME()->DriftVelocity();
    
//   }

  /**
   * @brief makes a map from daughter track IDs to their mother (shower) track ID
   *
   * The info comes from the MCShowerObjecs themselves. They store daughter IDs
   *  in a vector<unsigned int> accessed by mcshower::DaughterTrackID().
   *
   * This map will let us map the track ids stored in each pixel of the
   *   the instance ID map to its original shower.
   * 
   */
  void MCParticleGraph::_fill_shower_daughter2mother_map( const std::vector<larlite::mcshower>& mcsh_v )
  {
    LARCV_DEBUG() << "daughter2mother fill" << std::endl;
    _shower_daughter2mother.clear();

    
    for (auto const& mcsh : mcsh_v ) {      
      int showerid = mcsh.TrackID();
      std::vector<unsigned int> dlist = mcsh.DaughterTrackID();
      std::sort( dlist.begin(), dlist.end() );
      for (auto const& daughterid : dlist ) {
        _shower_daughter2mother[daughterid]= showerid;
        //LARCV_DEBUG() << "  " << daughterid << " -> " << showerid << std::endl;
      }
    }
    LARCV_INFO() << "Num entries in daughter2mother map: " << _shower_daughter2mother.size() << std::endl;
  }

  int MCParticleGraph::_define_neutrino_interaction_nodes( larlite::storage_manager& ioll )
  {
    larlite::event_mctrack* ev_mctrack
      = (larlite::event_mctrack*)ioll.get_data(larlite::data::kMCTrack,"mcreco");
    larlite::event_mcshower* ev_mcshower
      = (larlite::event_mcshower*)ioll.get_data(larlite::data::kMCShower,"mcreco");
    
    return _define_neutrino_interaction_nodes( *ev_mctrack, *ev_mcshower );
  }
  
  int MCParticleGraph::_define_neutrino_interaction_nodes( const larlite::event_mctrack& ev_mctrack,
							 const larlite::event_mcshower& ev_mcshower )
  {
    // first we collect nodes with neutrino origin
    bool exclude_neutrons = false;
    std::vector< MCPGNode* > _nu_primary_v = getNeutrinoPrimaryParticles( exclude_neutrons );
    LARCV_DEBUG() << "Number of nu primaries returned: " << _nu_primary_v.size() << std::endl;

    // if we have access to MCTruth info, we should use it to define vertex locations.
    // otherwise we use some distance threshold. 0.3 mm, the pitch length?
    std::map< int, std::set<int> > _nu_collected_primaries_v;
    _nu_vertices_v.clear();

    
    for ( auto& pnode : _nu_primary_v ) {
      // we test the vertex for every primary
      LARCV_DEBUG() << "Considering Nu Primary node[" << pnode->nodeidx << "] tid=" << pnode->tid << std::endl;

      std::vector<float> start = {0,0,0,0};

      if ( pnode->isTrackObject() ) {
	auto const& track = _retrieve_mctrackobject( pnode, ev_mctrack );
	for (int i=0; i<4; i++) {
	  start[i] = track.Start().Position()[i];
	}
      }
      else if (pnode->isShowerObject()) {
	auto const& shower = _retrieve_mcshowerobject( pnode, ev_mcshower );
	for (int i=0; i<4; i++) {
	  start[i] = shower.Start().Position()[i]; // creation point in geant4, for photon, not the same as conversion point where visible EM cascade begins
	}
      }
      else if (pnode->isGenieFinalStateObject()) {
	start = pnode->start;
      }
      else {
        LARCV_DEBUG() << "Node Primary neither a shower nor track object: " << pnode->nodeidx << std::endl;
	LARCV_DEBUG() << "node: " << strNodeInfo( *pnode ) << std::endl;
	continue;
      }

      bool found_vertex_match = false;
      int idx_matched_vertex = -1;
      float closest_match = 1e9;
      std::vector<float> pos = start;
      
      LARCV_DEBUG() << "Nu primary start: " << pos[0] << " " << pos[1] << " " << pos[2] << " " << pos[3] << std::endl;      
      
      for ( int idx_vertex=0; idx_vertex<(int)_nu_vertices_v.size(); idx_vertex++ ) {
	
	auto& vertex = _nu_vertices_v.at(idx_vertex);
	
	// distance to existing vertex
	float dist = 0.;
	for (int i=0; i<3; i++) {
	  dist += ( pos[i]-vertex[i] )*( pos[i]-vertex[i] );
	}
	dist = sqrt(dist);

	// update closest match
	if ( dist < closest_match ) {
	  idx_matched_vertex = idx_vertex;
	  closest_match = dist;

	  // qualifies as found?
	  if ( dist < _kNuVertexDistCutoff_cm ) {
	    found_vertex_match = true;
	    idx_matched_vertex = idx_vertex;
	  }
	}

      }//end of loop over established nu vertices
      LARCV_DEBUG() << "result of vertex search: closest=" << closest_match << " idx_matched=" << idx_matched_vertex << " found_match=" << found_vertex_match << std::endl;

      if ( !found_vertex_match ) {
	// new neutrino vertex defined using position
	LARCV_DEBUG() << "New vertex defined." << std::endl;
	std::vector<float> new_vertex = { (float)pos[0], (float)pos[1], (float)pos[2], (float)pos[3] };
	_nu_vertices_v.push_back( new_vertex );

	int nu_vertex_id = (int)_nu_vertices_v.size()-1; // using position in _nu_vertices as an id number (should use a struct I know)
	std::set<int> nu_primary_set;
	nu_primary_set.insert( pnode->nodeidx ); // add node index
	_nu_collected_primaries_v[ nu_vertex_id ] = nu_primary_set;
      }
      else {
	// found match
	LARCV_DEBUG() << "Matched primary to existing vertex. IDX=" << idx_matched_vertex << std::endl;
	auto it=_nu_collected_primaries_v.find( idx_matched_vertex );
	it->second.insert( pnode->nodeidx );
      }
	
    }//end of loop over neutrino primaries

    LARCV_DEBUG() << "Number of vertices defined: " << _nu_vertices_v.size() << std::endl;

    if ( _nu_vertices_v.size()==0 )
      return 0;

    // now that we have nu vertices, we need to define a new neutrino ancestor ID, then relabel ancestor IDs for daughters
    long max_geant4_trackid = -1;
    if ( _nu_vertices_v.size()>0 ) {
      for (auto pnode : node_v) {
	if ( pnode.tid>max_geant4_trackid ) {
	  max_geant4_trackid = pnode.tid;
	}
      }
    }
    LARCV_DEBUG() << "Starting with max trackid=" << max_geant4_trackid << " to assign nu vertex nodes" << std::endl;
    std::set<int> nu_attached_v; // gather list of indices that have been attached to neutrinos
    std::vector< MCPGNode* > nu_pnode_v;

    for (int inuvtx=0; inuvtx<(int)_nu_vertices_v.size(); inuvtx++) {
      
      LARCV_DEBUG() << "Build nu vertex node and assign daughters. [NU VTX IDX=" << inuvtx << "]" << std::endl;
      
      // get the next node index
      int nu_node_idx = (int)node_v.size();
      // need an acceptable fake trackID
      long fake_trackid = max_geant4_trackid+1;
      max_geant4_trackid++;
      int type_id = 2;
      int pid = -1;
      int mtid = fake_trackid;
      float energy = 0.0;
      std::string proc = "nuvertex";
      MCPGNode nu_node( nu_node_idx, type_id, fake_trackid, inuvtx,
		      pid, _eventRootNode,
		      mtid, energy, proc );

      nu_node.origin = 1;
      nu_node.aid = fake_trackid;
      nu_node.mtid = -1;

      // insert into node vector
      node_v.emplace_back( std::move(nu_node) );

      // get pointer
      MCPGNode* pnu_node = &(node_v.at(nu_node_idx));
      nu_pnode_v.push_back( pnu_node );

      // now we need to
      // (1) change the nu primaries associated to this interaction
      //     to list their mother node to this node representing the neutrino interaction
      // (2) relabel the ancestor ID of all primary neutrinos to this ID
      auto it_prim = _nu_collected_primaries_v.find( inuvtx );
      if ( it_prim==_nu_collected_primaries_v.end() )
	continue;
      
      for ( auto& nodeidx : it_prim->second ) {

	LARCV_DEBUG() << "Add primary, node=" << nodeidx << ", to nu vertex[" << inuvtx << "]" << std::endl;

	LARCV_DEBUG() << "  reassign mother to nu vertex pnode=" << pnu_node << std::endl;
	MCPGNode* pnode_nuprim = &(node_v.at(nodeidx));
	pnode_nuprim->mother = pnu_node;

	// collect daughters
	std::vector<MCPGNode*> prim_daughters = getNodeAndDescendentsFromTrackID( pnode_nuprim->tid );
	LARCV_DEBUG() << "Collect descendents of nu  primary node[ " << nodeidx << "] ndaughters=" << (int)prim_daughters.size()-1 << std::endl;	

	// reset the ancestor id of all of these nodes to the new fake track ID for nu interaction
	// note: the function above returns the node of the starting track id
	for ( auto& pdnode : prim_daughters ) {
	  pdnode->aid = (int)fake_trackid;
	}
	
	// add this primary to the daughter list of the nu node
	pnu_node->daughter_idx_v.push_back( pnode_nuprim->nodeidx );
	pnu_node->daughter_v.push_back( pnode_nuprim );
	LARCV_DEBUG() << "Add to NuVertex node list of daughters: now " << pnu_node->daughter_v.size() << std::endl;
	nu_attached_v.insert( pnode_nuprim->nodeidx );
	
	// add its E
	pnu_node->E_MeV += pnode_nuprim->E_MeV;
	LARCV_DEBUG() << "Add to NuVertex energy: now " << pnu_node->E_MeV << " MeV" << std::endl;
      }//end of loop over primary node

      nu_attached_v.insert( nu_node_idx );

    }//end of loop over newly creatd neutrino vertices

    // now we have to redefine the root node's connections
    _eventRootNode = &(node_v[0]);
    std::vector<MCPGNode*> all_prim_v = getPrimaryParticles(exclude_neutrons);

    _eventRootNode->daughter_idx_v.clear();
    _eventRootNode->daughter_v.clear();
    
    
    for (auto& pnode : all_prim_v ) {
      auto it_attached = nu_attached_v.find( pnode->nodeidx );
      if ( it_attached==nu_attached_v.end() ) {
	// not attached to neutrino vertex, so add to root node
	_eventRootNode->daughter_idx_v.push_back( pnode->nodeidx );
	_eventRootNode->daughter_v.push_back( pnode );
      }
    }

    // attach the nu vertex nodes
    for (auto& pnode : nu_pnode_v ) {
      _eventRootNode->daughter_idx_v.push_back( pnode->nodeidx );
      _eventRootNode->daughter_v.push_back( pnode );
    }
    
    return _nu_vertices_v.size();
  }

  long MCParticleGraph::getShowerMotherID( long trackid )
  {
    auto it = _shower_daughter2mother.find( trackid );
    if ( it!=_shower_daughter2mother.end() ) {
      return it->second;
    }
    return -1;
  }

  long MCParticleGraph::getParticleID( long trackid )
  {
    auto pnode  = findTrackID( trackid );
    if ( pnode==nullptr )
      return -1;

    return pnode->pid;
  }

  long MCParticleGraph::getAncestorID( long trackid )
  {
    auto pnode  = findTrackID( trackid );
    if ( pnode==nullptr )
      return -1;

    return pnode->aid;
  }

  
}
}
