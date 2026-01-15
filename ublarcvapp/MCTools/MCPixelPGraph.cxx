#include "MCPixelPGraph.h"

#include <sstream>
#include <set>

// larcv
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/DataFormatTypes.h"
//#include "larcv/core/ROOTUtil/ROOTUtils.h"

// larlite
#include "larlite/DataFormat/mctrack.h"
#include "larlite/DataFormat/mcshower.h"
#include "larlite/DataFormat/mctruth.h"
#include "larlite/LArUtil/LArProperties.h"
#include "larlite/LArUtil/Geometry.h"
#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"

// ublarcvapp
#include "ublarcvapp/dbscan/DBScan.h"

#include "crossingPointsAnaMethods.h"
#include "MCPos2ImageUtils.h"

namespace ublarcvapp {
namespace mctools {

  void MCPixelPGraph::buildgraph( larcv::IOManager& iolcv,
                                  larlite::storage_manager& ioll ) {

    ev_adc = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, adc_tree );
    ev_seg = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, "segment" );
    ev_ins = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, "instance" );
    ev_anc = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, "ancestor" );
    ev_larflow = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, "larflow" );

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
    
    
    ev_mctrack  = (larlite::event_mctrack*) ioll.get_data( larlite::data::kMCTrack,  "mcreco" );
    ev_mcshower = (larlite::event_mcshower*)ioll.get_data( larlite::data::kMCShower, "mcreco" );
    ev_mctruth  = (larlite::event_mctruth*) ioll.get_data( larlite::data::kMCTruth,  "generator" );

    buildgraph( ev_adc->Image2DArray(),
                ev_seg->Image2DArray(),
                ev_ins->Image2DArray(),
                ev_anc->Image2DArray(),
		ev_larflow->Image2DArray(),
                *ev_mcshower, *ev_mctrack, *ev_mctruth );

  }

  /**
   * @brief build only the particle graph (no pixel scan)
   *
   */
  void MCPixelPGraph::buildgraphonly( larlite::storage_manager& ioll )
  {
    
    ev_mctrack  = (larlite::event_mctrack*) ioll.get_data( larlite::data::kMCTrack,  "mcreco" );
    ev_mcshower = (larlite::event_mcshower*)ioll.get_data( larlite::data::kMCShower, "mcreco" );
    ev_mctruth  = (larlite::event_mctruth*) ioll.get_data( larlite::data::kMCTruth,  "generator" );

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
				                          const std::vector<larcv::Image2D>& larflow_v,
                                  const larlite::event_mcshower& shower_v,
                                  const larlite::event_mctrack&  track_v,
                                  const larlite::event_mctruth&  mctruth_v )
  {

    clear();
    
    buildgraphonly( shower_v, track_v, mctruth_v );   
    
    std::vector<float> threshold_v(adc_v.size(),10.0);
    _scanPixelData( adc_v, segment_v, instance_v, ancestor_v, threshold_v );
    
    // next we fix photons
    for ( auto& node : node_v ) {
      if ( node.pid==22 || abs(node.pid)==11 ) {
	      std::vector<float > startpt_info
	        = fixingPhotonStartPoints( node,
				     instance_v, ancestor_v,
				     adc_v, larflow_v );
        // update edep position of the node for photons only
        if ( node.pid==22 ) {
          for (int v=0; v<4; v++) {
            node.first_edep_pos[v] = startpt_info[v];    // set the location of the starting edep position (detectort coordinates)
            node.imgpos4_edep[v]   = startpt_info[4+v];  // set the location ion the image
          }
        }
        // if ( node.origin==2 && abs(node.pid)==11 ) {
        //   // the cosmics timing are messed up (at least corsika samples are?)
        //   for (int v=0; v<4; v++) {
        //     node.start[v] = startpt_info[v];    // set the location of the starting edep position (detectort coordinates)
        //     node.first_edep_pos[v] = startpt_info[v];    // set the location of the starting edep position (detectort coordinates)
        //     node.imgpos4_edep[v]   = startpt_info[4+v];  // set the location ion the image
        //   }
        // }
      }
    }
    
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
    Node_t neutrino ( node_v.size(), -1, 0, 0, -1 );


    std::set<int> tid_list;
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

    // load spacechargemicroboone
    larutil::SpaceChargeMicroBooNE sce;

    struct NuPart_t {

      int geantid;
      int pdg;
      float E_MeV;
      std::vector<float> pos;
    };

    std::vector< NuPart_t > nu_part_v;

    // collect from mcreco tracks
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
      tracknode.process = mct.Process();
      if ( mct.PdgCode()==2212 ) tracknode.E_MeV -= 938.0;
      else if ( mct.PdgCode()==2112 ) tracknode.E_MeV -= 940.0;
      else if ( abs(mct.PdgCode())==13 )   tracknode.E_MeV -= 105.;
      else if ( abs(mct.PdgCode())==211 )  tracknode.E_MeV -= 135.;
      tracknode.origin = mct.Origin();

      tid_list.insert( tracknode.tid );
      
      // real position, time
      tracknode.start.resize(4);
      tracknode.start[0] = mct.Start().X();
      tracknode.start[1] = mct.Start().Y();
      tracknode.start[2] = mct.Start().Z();
      tracknode.start[3] = mct.Start().T();

      
      if ( tracknode.origin==1 ) {
	      // store nu particle
	      NuPart_t nuparticle;
	      nuparticle.geantid = tracknode.tid;
	      nuparticle.pdg = tracknode.pid;
	      nuparticle.E_MeV = mct.Start().E();
	      nuparticle.pos = tracknode.start;
	      nu_part_v.push_back( nuparticle );
      }
      


      // set start position.
      // try to move it into the tpc first
      tracknode.first_edep_pos = std::vector<float>(4,0);
      tracknode.first_tpc_pos = std::vector<float>(4,0);
      tracknode.first_img_pos = std::vector<float>(4,0);
      tracknode.imgpos4 = std::vector<float>(4,0);
      tracknode.imgpos4_edep = std::vector<float>(4,0);            
      tracknode.imgpos4_start = ublarcvapp::mctools::MCPos2ImageUtils::Get()->truepos_to_imagepos( tracknode.start[0],
												   tracknode.start[1],
												   tracknode.start[2],
												   tracknode.start[3],
												   true );
	

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

            recopos = ublarcvapp::mctools::MCPos2ImageUtils::Get()->truepos_to_imagepos( xyz[0], xyz[1], xyz[2], xyz[3], true );
            if ( !inimage && recopos[3]>2400.0 && recopos[3]<2400+1008*6 ) {
              // in the image
              inimage = true;
              tracknode.first_img_pos = xyz;
              tracknode.imgpos4 = recopos;
              tracknode.imgpos4_edep = recopos;
            }
          }
          if ( intpc && inimage ) {
            // no need to keep searching
            break;
          }
        }//end of mcstep loop
      }//end of if more than 0 mcsteps

      node_v.emplace_back( std::move(tracknode) );
    }

    // collect from mcshowers
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
      
      //if ( mcsh.Origin()==1 ) {
      // neutrino origin

      Node_t showernode( node_v.size(), 1, mcsh.TrackID(), vidx, mcsh.PdgCode() );
      showernode.E_MeV = mcsh.Start().E();
      showernode.process = mcsh.Process();
      showernode.origin = mcsh.Origin();      
      showernode.start = std::vector<float>{ (float)mcsh.Start().X(), (float)mcsh.Start().Y(), (float)mcsh.Start().Z(), (float)mcsh.Start().T() };
      showernode.aid  = mcsh.AncestorTrackID();
      showernode.mtid = mcsh.MotherTrackID();

      showernode.first_edep_pos = std::vector<float>(4,0);
      showernode.first_tpc_pos  = std::vector<float>(4,0);
      showernode.first_img_pos  = std::vector<float>(4,0);
      showernode.imgpos4 = std::vector<float>(4,0);
      showernode.imgpos4_edep = std::vector<float>(4,0);
      showernode.imgpos4_start = std::vector<float>(4,0);

      bool x_isinf = false;
      if ( mcsh.DetProfile().X()>1.0e100 || mcsh.DetProfile().X()<-1.0e100 )
	      x_isinf = true;
      if ( !x_isinf && !std::isnan( mcsh.DetProfile().X() ) ) {
        //std::cout << "inf test: " << mcsh.DetProfile().X() << " " << std::isinf( mcsh.DetProfile().X() ) << std::endl;
        std::vector<float> detprofile = { (float)mcsh.DetProfile().X(), (float)mcsh.DetProfile().Y(), (float)mcsh.DetProfile().Z(), (float)mcsh.DetProfile().T() };
        showernode.first_edep_pos = detprofile;
        showernode.first_tpc_pos  = detprofile;
        showernode.first_img_pos  = detprofile;
        //_get_imgpos( detprofile, showernode.imgpos4, sce, false );
        showernode.imgpos4 = ublarcvapp::mctools::MCPos2ImageUtils::Get()->truepos_to_imagepos( showernode.first_edep_pos[0],
                              showernode.first_edep_pos[1],
                              showernode.first_edep_pos[2],
                              showernode.first_edep_pos[3],
                              true );
        if ( showernode.imgpos4.size()==0 )
          showernode.imgpos4.resize(4,0);

        if ( abs(showernode.pid)==11 ) {
          showernode.imgpos4_edep = ublarcvapp::mctools::MCPos2ImageUtils::Get()->truepos_to_imagepos( showernode.start[0],
                                    showernode.start[1],
                                    showernode.start[2],
                                    showernode.start[3],
                                    true );
          if ( showernode.imgpos4_edep.size()==0 )
            showernode.imgpos4_edep.resize(4,0);
          
          showernode.first_edep_pos = showernode.start;
        }
        else {
          // photons
          showernode.imgpos4_edep = ublarcvapp::mctools::MCPos2ImageUtils::Get()->to_imagepos( showernode.first_edep_pos[0],
                                  showernode.first_edep_pos[1],
                                  showernode.first_edep_pos[2],
                                  showernode.first_edep_pos[3] );
          if ( showernode.imgpos4_edep.size()==0 ) {
            showernode.imgpos4_edep.resize(4,0);
          }
          else {
            showernode.imgpos4_edep[3] += 12.0; // hacky fix for photons
          }
        }
        
        showernode.imgpos4_start = ublarcvapp::mctools::MCPos2ImageUtils::Get()->truepos_to_imagepos( showernode.start[0],
                                    showernode.start[1],
                                    showernode.start[2],
                                    showernode.start[3],
                                    true );
        if ( showernode.imgpos4_start.size()==0 ) {
          showernode.imgpos4_start.resize(4,0);
        }
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
    
    // collect from mctruth
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
            Node_t* mother = nullptr;
            int mid = -1;
            float energy = part.Momentum(0)[3]*1e3 - part.Mass()*1.0e3; // in GeV, convert to MeV
            std::vector<float> start { (float)part.Position(0)[0],
                                      (float)part.Position(0)[1],
                                      (float)part.Position(0)[2],
                                      (float)part.Position(0)[3]};
            std::vector<float> imgpos4(4,0);
            _get_imgpos( start, imgpos4, sce, true );
            
            //LARCV_DEBUG() << "Creating Mother node from GENIE final states: tid=" << tid << " type=" << 3 << std::endl;
            Node_t fsnode( nodeidx, type, tid, vidx, pid, mother, mid, energy, "primary" );
            fsnode.start = start;
            fsnode.imgpos4 = imgpos4;
            fsnode.mtid = tid;
            fsnode.aid  = tid;
            fsnode.origin = 1; // neutrino origin (from genie)
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
    //_adoptNeutrinoOrphans( &mctruth_v );

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
    
    //printAllNodeInfo();
    //printGraph();
  }

  /**
   * @brief wrapper function to retrieve Node's corresponding mctrack object from larlite container
   */
  const larlite::mctrack& 
  MCPixelPGraph::_retrieve_mctrackobject( const Node_t* node,
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
  const larlite::mcshower&  MCPixelPGraph::_retrieve_mcshowerobject( const Node_t* node,
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
  const larlite::mctrack&  MCPixelPGraph::_retrieve_mctrackobject( const Node_t* node,
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
  const larlite::mcshower&  MCPixelPGraph::_retrieve_mcshowerobject( const Node_t* node,
								     larlite::storage_manager& ioll,
								     std::string producername )
  {
    const larlite::event_mcshower* ev_mcshower
      = (larlite::event_mcshower*)ioll.get_data(larlite::data::kMCShower, producername );

    return _retrieve_mcshowerobject( node, *ev_mcshower );
  }
  
  
  /**
   * @brief locate Node_t in node_v using trackid (from geant4)
   * 
   * we first search our vector of Node_t objects, which are sorted by track id.
   *
   * we then search the keys of the _shower_daughter2mother map, where the keys
   *  are track id of particles made by the initial shower.
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
        //LARCV_DEBUG() << "  found map to mother: " << it_showerdaughter->second << std::endl;
        dummy.tid = it_showerdaughter->second;
      }
      else {
        // still nope
        return nullptr;
      }
      // with the mother shower's trackid, try to find the node again
      it = std::lower_bound( node_v.begin(), node_v.end(), dummy );
      //if ( it!=node_v.end() )
      //  LARCV_DEBUG() << "  mother id maps to existing node" << std::endl;
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
  std::string MCPixelPGraph::strNodeInfo( const Node_t& node ) 
  {

    int hasmother = ( node.mother ) ? 1 : 0;
    
    std::stringstream ss;
    //ss << "node[" << node.nodeidx << "," << &node << "] "
    ss << "node[" << node.nodeidx << "] "
       << " (type,vidx)=(" << node.type << "," << node.vidx << ") "
       << " origin=" << node.origin
       << " p=" << node.process
       << " tid=" << node.tid
       << " mtid=" << node.mtid      
       << " aid=" << node.aid
       << " pdg=" << node.pid
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
      ss << "    edep-tpcpos(x,y,z,t)=(" << node.first_edep_pos[0] << ","
          << node.first_edep_pos[1] << ","
          << node.first_edep_pos[2] << ","
          << node.first_edep_pos[3]*1.0e-3
          << " us)" << std::endl;
    
    if ( node.imgpos4.size()>=4 )
      ss << "    tpc-imgpos4(u,v,y,tick)=(" << node.imgpos4[0] << "," 
         << node.imgpos4[1] << "," 
         << node.imgpos4[2] << "," 
         << node.imgpos4[3] << " tick) " << std::endl;
    
    if (node.imgpos4_edep.size()>=4)
      ss << "    edep-imgpos4(u,v,y,tick)=(" << node.imgpos4_edep[0] << "," 
         << node.imgpos4_edep[1] << "," 
         << node.imgpos4_edep[2] << "," 
         << node.imgpos4_edep[3] << " tick) " 
         << std::endl;
    
    ss << "    npixs=(";
    for ( size_t i=0; i<node.pix_vv.size(); i++ ) {
      ss << node.pix_vv[i].size()/2;
      if ( i+1<node.pix_vv.size() ) ss << ", ";      
    }
    ss << ")";
    ss << std::endl;

    ss << "    pixelsum=(";
    for ( size_t i=0; i<node.pixsum_v.size(); i++ ) {
      ss << node.pixsum_v[i] << ", ";      
    }
    ss << ")";
    ss << std::endl;
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
   * @brief scan the adc and truth images and associate them with the different particle Node_t
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
      node.pixval_vv.resize(_nplanes);
      node.pixsum_v.resize(_nplanes,0.0);
      for ( size_t p=0; p<_nplanes; p++ ) {
        node.pix_vv[p].clear();
        node.pixval_vv[p].clear();
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
      int num_neg_shower_ids = 0;

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
          // if ( adc<threshold )
          //   continue;

          if ( adc>=threshold )
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

          // shower pixels have a negative track ID for some reason
          // does this always occur?
          if ( tid<0 && (seg==(int)larcv::kROIEminus || seg==(int)larcv::kROIGamma) ) {
            tid *= -1;
	          num_neg_shower_ids++;
	        }

          bool has_id_label = false;
                if ( tid>0 || aid>0 ) {
                  nabove_thresh_withlabel[p]++;
            has_id_label = true;
          }
	  
          if ( aid>0 && (seg==(int)larcv::kROIEminus || seg==(int)larcv::kROIGamma) ) {
            shower_ancestor_ids.insert( aid );
          }

          if ( !has_id_label ) {
            // can't truth associate this pixel to any node
            continue;
          }

          Node_t* node = nullptr;

          if ( tid>0 ) {
            // first we use the instance ID
	          // this implicitly uses the _shower_daughter2mother map we filled earlier.
            node = findTrackID( tid );
            if ( node==nullptr && adc>10.0 )
              LARCV_DEBUG() << "  no node for above threshold charge-pixel tid=" << tid << std::endl;
            // if ( node && node->tid!=tid )
            //   node = nullptr; // reset (what's this?)
          }
          
          // use ancestor if we could not find the node
          if ( !node && aid>0 ) {
	          // this implicitly uses the _shower_daughter2mother map we filled earlier.	    
            node = findTrackID( aid );
            if ( node && node->tid!=aid )
              node = nullptr;
          }

          if ( node ) {
	          // if the address is not null, this means we found a particle node
	          // store the pixel.
	    
            // if ( node->tid!=tid && node->aid!=aid ) {
            //   std::cout << "pixel assigned without matching tid or aid exactly: "
            //             << " pix-tid=" << tid << " pix-aid=" << aid
            //             << " node-tid=" << node->tid << " node-mid=" << node->mtid << " node-aid=" << node->aid
            //             << std::endl;
            // }
            
            nassigned[p]++;
            node->pix_vv[p].push_back( tick );
            node->pix_vv[p].push_back( wire );
            node->pixval_vv[p].push_back( adc );
            node->pixsum_v[p] += adc;
          }
          else {
            _unassigned_pixels_vv[p].push_back( tick );
            _unassigned_pixels_vv[p].push_back( wire );
          }
          
        }//end of loop over columns
      }//end of loop over rows
      
    }//end of loop over planes

    // now make bounding boxes
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
   * @brief get pixels associated with node and its descendents
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
                                   larutil::SpaceChargeMicroBooNE& sce,
				   bool apply_sce )
  {

    imgpos4.resize(4,0);
    
    // apparent position according to image
    std::vector<double> dpos(3,0);
    for (int i=0; i<3; i++) {
      dpos[i]   = realpos4[i];
    }

    std::vector<float>  txyz(4,0);
    if ( apply_sce ) {
      std::vector<double> offset = sce.GetPosOffsets( dpos[0], dpos[1], dpos[2] );        
      dpos[0] = dpos[0] - offset[0] + 0.7;
      dpos[1] = dpos[1] + offset[1];
      dpos[2] = dpos[2] + offset[2];
    }
    
    for (int i=0; i<3; i++) {
      txyz[1+i] = realpos4[i];
    }
    txyz[0] = realpos4[3];
    float tick = 0.;
    if (apply_sce)
      tick = CrossingPointsAnaMethods::getTick( txyz, 4050.0, &sce );
    else
      tick = CrossingPointsAnaMethods::getTick( txyz, 4050.0, NULL );
    
    for (int i=0; i<3; i++) {
      imgpos4[i] = dpos[i];
    }
    imgpos4[3] = tick;

    // now make x an apparent x
    imgpos4[0] = (tick-3200)*0.5*larutil::LArProperties::GetME()->DriftVelocity();
    
  }

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
        //LARCV_DEBUG() << "  " << daughterid << " -> " << showerid << std::endl;
      }
    }
    LARCV_INFO() << "Num entries in daughter2mother map: " << _shower_daughter2mother.size() << std::endl;
  }

  int MCPixelPGraph::_define_neutrino_interaction_nodes( larlite::storage_manager& ioll )
  {
    larlite::event_mctrack* ev_mctrack
      = (larlite::event_mctrack*)ioll.get_data(larlite::data::kMCTrack,"mcreco");
    larlite::event_mcshower* ev_mcshower
      = (larlite::event_mcshower*)ioll.get_data(larlite::data::kMCShower,"mcreco");
    
    return _define_neutrino_interaction_nodes( *ev_mctrack, *ev_mcshower );
  }
  
  int MCPixelPGraph::_define_neutrino_interaction_nodes( const larlite::event_mctrack& ev_mctrack,
							 const larlite::event_mcshower& ev_mcshower )
  {
    // first we collect nodes with neutrino origin
    bool exclude_neutrons = false;
    std::vector< Node_t* > _nu_primary_v = getNeutrinoPrimaryParticles( exclude_neutrons );
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
    std::vector< Node_t* > nu_pnode_v;

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
      Node_t nu_node( nu_node_idx, type_id, fake_trackid, inuvtx,
		      pid, _eventRootNode,
		      mtid, energy, proc );

      nu_node.origin = 1;
      nu_node.aid = fake_trackid;
      nu_node.mtid = -1;

      // insert into node vector
      node_v.emplace_back( std::move(nu_node) );

      // get pointer
      Node_t* pnu_node = &(node_v.at(nu_node_idx));
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
	Node_t* pnode_nuprim = &(node_v.at(nodeidx));
	pnode_nuprim->mother = pnu_node;

	// collect daughters
	std::vector<Node_t*> prim_daughters = getNodeAndDescendentsFromTrackID( pnode_nuprim->tid );
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
    std::vector<Node_t*> all_prim_v = getPrimaryParticles(exclude_neutrons);

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

  void MCPixelPGraph::clear()
  {

    _unassigned_pixels_vv.clear();    
    node_v.clear();
    
    _eventRootNode = nullptr;    
    _nu_vertices_v.clear();
    _shower_daughter2mother.clear();
    //_map_trackid_to_nu_ancestor_v.clear();

    _nodeidx_to_photonlist_index_v.clear();
    _true_photon_v.clear();
    _true_photon_plane_trunkimg_vv.clear();
    _true_photon_plane_pixset_vv.clear();
    _true_photon_plane_pixsum_vv.clear();

    _empty_photon_pointlist_v.clear();
    _empty_pixsum_v.clear();

    // ev_adc = nullptr;
    // ev_seg = nullptr;
    // ev_ins = nullptr;
    // ev_anc = nullptr;
    
  }

  void MCPixelPGraph::_adoptNeutrinoOrphans( const larlite::event_mctruth* ev_mctruth )
  {
    // we scan our nodes without mother nodes that come from neutrinos (origin=1)
    // if we have mctruth, we look to find mother. then we make new mother node and finish graph

    std::map<int,std::vector<int> > final_state_map;
    int ifs = 0;
    if ( ev_mctruth ) {
      int imctruth=0;
      for ( auto const& mctruth : *ev_mctruth ) {
	//LARCV_DEBUG() << "MCTRUTH[" << imctruth << "] --------------" << std::endl;
	int imcpart = 0;
	for ( auto& part : mctruth.GetParticles() ) {
	  if ( part.StatusCode()==1 ) {
	    int geant_trackid = ifs+1;
	    ifs++;
	    std::vector<int> index = {imctruth,imcpart};
	    final_state_map[geant_trackid] = index;
	  }
	  imcpart++;	  
	}
	imctruth++;
      }//end of mctruth map
    }//end of if have mctruth

    std::set<int> tid_list;
    for ( auto& node : node_v )
      tid_list.insert(node.tid);

    std::map< int, int > newmom_tid_to_nodeidx;
    
    for ( auto& node : node_v ) {
      // has mother or is good neutrino mother
      if ( node.origin!=1 ) {
	continue;
      }

      // check if the mother node is not already in the node list
      auto it_tid = tid_list.find( node.mtid );
      if ( it_tid!=tid_list.end() ) {
	// has a mother in the list, no need to gather one from the mctruth or create one
	continue;
      }

      // check if we might have a mother now
      auto it_newmom = newmom_tid_to_nodeidx.find( node.mtid );
      if ( it_newmom==newmom_tid_to_nodeidx.end() ) {
	// try ancestor id
	it_newmom = newmom_tid_to_nodeidx.find( node.aid );
      }

      if ( it_newmom!=newmom_tid_to_nodeidx.end() ) {
	// connected to a new mom!
	//int momtid = it_newmom->first;	
	int momidx = it_newmom->second;
	auto& momnode = node_v.at(momidx);
	//momnode.daughter_v.push_back( &node );
	//momnode.daughter_idx_v.push_back( node.nodeidx );
	//node.mother = &momnode;
	LARCV_DEBUG() << "Connected primary orphan to a newly created mom! New mom tid=" << momnode.tid << std::endl;
	continue;
      }

      LARCV_DEBUG() << "node[" << node.nodeidx << "] with mtid=" << node.mtid << " looking to find mom or create it" << std::endl;

      bool created_mother = false;
      if ( ev_mctruth ) {
	// if we have the mctruth, try to look for missing mother id
	std::vector<int> momindex;
	int mom_tid = -1;
	
	auto it_fs = final_state_map.find( node.mtid );
	if ( it_fs!=final_state_map.end() ) {
	  // has mother, create a node
	  momindex = it_fs->second;
	  mom_tid = node.mid;
	}
	else {
	  // try using ancestor
	  auto it_fsa = final_state_map.find( node.aid );
	  auto it_tid2 = tid_list.find( node.aid );
	  // create a mother if (1) found in final state list AND (2) not found in tid_list (i.e. existing node list)
	  if ( it_fsa!=final_state_map.end() && it_tid2==tid_list.end() ) {
	    momindex = it_fsa->second;
	    mom_tid = node.aid;
	  }
	}

	if ( mom_tid<1 )
	  continue;

	LARCV_DEBUG() << "Creating mother node from Genie FSI info. mom_tid=" << mom_tid << std::endl;
	std::cin.get();

	if ( momindex.size()==2 ) {
	  // create a mother node from mcparticle info
	  const larlite::mcpart& part = ev_mctruth->at(momindex[0]).GetParticle( momindex[1] );
	  int nodeidx = node_v.size();
	  int type = 3; // genie-final-state
	  int tid  = mom_tid;	  
	  int vidx = momindex[1];
	  int pid = part.PdgCode();
	  Node_t* mother = &(node_v[0]); // the root node
	  int mid = 0;
	  int energy = part.Momentum(0)[3];
	  std::vector<float> start { (float)part.Position(0)[0],
	    (float)part.Position(0)[1],
	    (float)part.Position(0)[2],
	    (float)part.Position(0)[3]};
	  LARCV_DEBUG() << "Creating Mother node from GENIE final states: tid=" << tid << " type=" << 3 << std::endl;
	  Node_t momnode( nodeidx, type, tid, vidx, pid, mother, mid, energy, "primary" );
	  momnode.start = start;
	  momnode.imgpos4 = start;
	  momnode.mtid = tid;
	  momnode.aid  = tid;
	  momnode.origin = 1; // neutrino origin (from genie)
	  
	  //momnode.daughter_v.push_back( &node );
	  //momnode.daughter_idx_v.push_back( node.nodeidx );

	  newmom_tid_to_nodeidx[momnode.tid] = momnode.nodeidx;
	  tid_list.insert(tid);
				
	  node_v.push_back( std::move(momnode) );
	  
	  //node.mother = &(node_v.back());

	  // add new mom to root node
	  //Node_t& rootnode = node_v.front();
	  //rootnode.daughter_v.push_back( &(node_v.back())  );
	  //rootnode.daughter_idx_v.push_back( node_v.back().nodeidx );
	  
	  created_mother = true;
	  
	}
	else {
	  LARCV_DEBUG() << "No mother created for this." << std::endl;
	}
      }//end of if has ev_mctruth

      if ( created_mother )
	continue;
      
    }//end of loop over nodes
    
  }//end of adopt orphans

  // std::vector<TH2D> MCPixelPGraph::makeTH2D( std::string hist_stem_name )
  // {
  //   // first get wire image into th2d
  //   std::vector< TH2D > hist_v;
  //   if ( ev_adc==nullptr )
  //     return hist_v;

  //   for ( auto& pmarker : vis_markers_v ) {
  //     delete pmarker;
  //   };
  //   vis_markers_v.clear();

  //   for ( auto& plabel : vis_label_v ) {
  //     delete plabel;
  //   }
  //   vis_label_v.clear();
    
  //   if ( !cvis ) {
  //     cvis = new TCanvas("cmcpg","MCPixel PGraph Canvas",800,2400);
  //     cvis->Divide(1,3);
  //   }
    
  //   hist_v = larcv::rootutils::as_th2d_v( ev_adc->as_vector(), hist_stem_name );

  //   // label the nodes
  //   for (auto& node : node_v ) {

  //     for (int p=0; p<(int)hist_v.size(); p++) {
  // 	cvis->cd(p+1);
  // 	auto& hist = hist_v.at(p);
  // 	//auto& meta = ev_adc->as_vector().at(p).meta();
  // 	std::cout << "Draw marker for node[" << node.nodeidx << "]-plane[" << p << "] wire=" << node.imgpos4_edep[p] << " tick=" << node.imgpos4_edep[3] << std::endl;
  // 	TMarker* m = new TMarker( node.imgpos4_edep[p], node.imgpos4_edep[3], 20 );
  // 	//m->SetMarkerSize(3);
  // 	m->SetMarkerColor(kMagenta);	
  // 	hist.Draw("colz");
  // 	m->Draw();
  // 	vis_markers_v.push_back( m );
  //     }
      
  //   }
  //   cvis->Update();
  //   //std::cout << "[enter] to continue" << std::endl;
  //   //std::cin.get();
  //   return hist_v;
  // }

  /**
   * @brief use larflow truth to make 3D spacepoints for a given particle
   *
   * we return all points formed by the intersections
   * 
   * @return a vector of 3d positions
   * 
   */
  std::vector< std::vector<float> >
  MCPixelPGraph::makeNode3DpointsFromLArFlowTruth( Node_t& node,
						   const std::vector<larcv::Image2D>& adc_v,						   
						   const std::vector<larcv::Image2D>& larflow_v )
  {

    // we want to make 3d points from the pixels we have for the shower.    
    std::vector< std::vector<float> > pos_vv;

    // we use the plane-to-plane "flow" map to define 3d spacepoints
    struct flowpt {
      int source_plane;
      int source_wire;
      int target_plane;
      int target_wire;
      int source_trackid;
      int target_trackid;
      float src_pixval;
      float pos[3];
      int posid[3];
      bool operator<( const flowpt& rhs ) const {
	      if (posid[0]<rhs.posid[0])
	        return true;
      	else if (posid[0]>rhs.posid[0])
	        return false;
	      // neither so [0] must be equal

	      if (posid[1]<rhs.posid[1] )
	        return true;
	      else if (posid[1]>rhs.posid[1])
	        return false;

	      // neither [1] worked, so must be equal
	      if (posid[2]<rhs.posid[2])
	        return true;
	      return false;
      };
    };

    // note: this is a hardcoded map of indices relevant only to the
    // way microboone encoded the 2d wire intersections in the larflow images
    // example of how to read these indices
    // targetmap[0][0] = 1 --> larfow image with index-0 (out of 6) contains the map from U pixels [0] mapped to pixels in the V plane (=1)
    // targetmap[0][1] = 2 --> larfow image with index-1 (out of 6) contains the map from U pixels [0] mapped to pixels in the Y plane (=2)
    // targetmap[1][0] = 0 --> larfow image with index-2 (out of 6) contains the map from V pixels [1] mapped to pixels in the U plane (=0)
    // .. etc ..
    int targetmap[3][2] = { {1,2},
			    {0,2},
			    {0,1} };

    // std::cout << "adc images: " << adc_v.size() << std::endl;
    // std::cout << "larflow images: " << larflow_v.size() << std::endl;

    std::set< flowpt > pt_v;
    
    // we loop over the pixels in the wire plane image that has been 
    // associated to this particle node
    for (size_t p=0; p<node.pix_vv.size(); p++) {
      auto const& pix_v = node.pix_vv.at(p);

      // the associated pixels are stored as a pair of (tick,wire) indices for each pixel
      // the information is unrolled as a single 1D vector, so the total number of pixels 
      // are half the length of the vector.
      size_t npix = (size_t)pix_v.size()/2;
      
      // a given wire plane is mapped to the two other wire planes,
      // so we have to handle two mappings, encoding 2D wire-intersections
      for (int iflowdir=0; iflowdir<2; iflowdir++) {
        // get the index of the larflow image holding the 2d wire-intersection information
        int ilarflow_image_index = 2*p+iflowdir;
        // use the index to get the image
	      auto& flowimg = larflow_v.at(ilarflow_image_index);

      
        // loop over the pixels for the given particle
	      for (size_t ipix=0; ipix<npix; ipix++) {
	        float tick = pix_v[2*ipix]; // the row dimension variable
	        float wire = pix_v[2*ipix+1]; // the column dimension variable
	        int src_plane = p; // index of the current plane
	        int tar_plane = targetmap[src_plane][iflowdir]; // index of the wire plane we map to

          // convert (tick,wire) to the (row,col) in the image array
	        int row = 0;
	        int col = 0;
          // out of bounds calls will trigger a throw, we will just catch it and move on
	        try {
	          row = flowimg.meta().row( tick );
	          col = flowimg.meta().col( wire );
	        }
	        catch(...) {
            LARCV_WARNING() << "Tried to access larflow pixel out-of-bounds: "
                            << "(wire,tick)=(" << wire << "," << tick << ")" 
                            << std::endl;
		  LARCV_WARNING() << "FlowImage[" << ilarflow_image_index << "] meta: " << flowimg.meta().dump() << std::endl;	    
	          continue;
	        }
	  
          // if we get here, this is a valid pixel inside the flow image
          // the map from the current plane's wire to the target plane's wire
          // is stored as the shift required to get to the correct wire.
	        float pixflow = -999;
		try {
		  pixflow = flowimg.pixel( row, col );
		}
		catch(...) {
		  pixflow = -999;
		  LARCV_WARNING() << "Tried to access larflow pixel out-of-bounds: "
				  << "(wire,tick)=(" << wire << "," << tick << ")" 
				  << std::endl;
		}
	        // std::cout << "building flow pixel: plane[" << p << "] (" << wire << "," << tick << ") "
	        // 	    << "(" << col << "," << row << ") : flow=" << pixflow << std::endl;
	  
          // the sentinal value of -999 means that this pixel did not contain a true mapping.
          // if we see it, we move on.
          // in retrospect, storing the column index of the pixel would have been better. 
          // would have allowed the sentinal value to be -999 that might have allowed for easier zero
          // removal compression. but maybe root's compression recognizes repeated values.
	        if ( pixflow<=-999 )
	          continue;

          // for non-sentinal values, we get the column index of the pixel in the target plane
          // that we map to.
	        int tar_wire = col + (int)pixflow;

          // check that the wire is not out-of-bounds for some reason
          // if truth properly made, this should be unnecessary
	        if ( tar_wire<0 || tar_wire>=(int)larutil::Geometry::GetME()->Nwires(tar_plane) )
	          continue;
	  
          // now we use larutil's geometry class to calcuate the (y,z) position of the wire intersection
	        UInt_t src_ch = larutil::Geometry::GetME()->PlaneWireToChannel( (UInt_t)src_plane, (UInt_t)wire );
	        UInt_t tar_ch = larutil::Geometry::GetME()->PlaneWireToChannel( (UInt_t)tar_plane, (UInt_t)tar_wire );
	        Double_t y,z;
	        bool crosses  = larutil::Geometry::GetME()->ChannelsIntersect( src_ch, tar_ch, y, z );
          // check that the intersection point provided is inside the TPC.
          // if inside, crosses==true
	        if ( crosses ) {
	          // create a spacepoint object, inside the TPC
	          flowpt pt;
	          pt.source_plane = src_plane;
	          pt.source_wire  = wire;
	          pt.target_plane = tar_plane;
	          pt.target_wire  = tar_wire;
            pt.source_trackid = node.tid;
            pt.target_trackid = 0; // not filled for now
            pt.src_pixval = adc_v.at(src_plane).pixel( row, (unsigned int)col );
            pt.pos[1] = y;
            pt.pos[2] = z;
            pt.pos[0] = (tick-3200.0)*0.5*larutil::LArProperties::GetME()->DriftVelocity(); //< to do, need larutil function to make experiment agnostic

            // convert the flow into an integer 'posid' set of indices.
            // posid is effectively defining a voxelized grid
            // we do this to avoid close duplicate 3d pts
            pt.posid[0] = (int)(pt.pos[0]*1000);
            pt.posid[1] = (int)(pt.pos[1]*1000);
            pt.posid[2] = (int)(pt.pos[2]*1000);
            // look for the posid triplet in our set
            auto it_pt = pt_v.find( pt );
            if ( it_pt==pt_v.end() ) {
              // a new index, insert into the set
              //std::cout << "  insert intersection: (" << pt.pos[0] << "," << pt.pos[1] << "," << pt.pos[2] << ")" << std::endl;
              pt_v.insert(pt);
            }
            else {
              // if we had a previous intersection at this 3d point, we use the one with the largest
              // energy deposition.
              // (note: this might actual bias towards suprious intersections due to overlapping ionization due to tomographic projection)
              // (note: but, since we restrict to node-only pixels, probably ok. no issues seen so far. so probably minor.)
              if ( (*it_pt).src_pixval < pt.src_pixval ) {
                // std::cout << " replace src plane [" << (*it_pt).source_plane << " to " << pt.source_plane << "] "
                // 	  << " with higher pixval " << (*it_pt).src_pixval << " vs. " << pt.src_pixval
                // 	  << std::endl;
                pt_v.insert(pt); // replaces?
              }
            }
          } //end of if possible wire intersection found (and calculated)
          else {
            //std::cout << "  failed intersection: (" << y << ", " << z << ")" << std::endl;
          }
	      } // end of loop over pixels in pix_vv list
      }// end of loop over flow directions (2 of them)
    }//end of loop over planes

    LARCV_INFO() << "Node[" << node.nodeidx << "] tid=" << node.tid << " pid=" << node.pid << std::endl;
    LARCV_INFO() << "  Number of spacepoints found from using larflow: " << pt_v.size() << std::endl;
    
    // store the points into our storage vector, pos_vv, which we return.
    for (auto& pt_info : pt_v ) {
      pos_vv.push_back( std::vector<float>{ pt_info.pos[0], pt_info.pos[1], pt_info.pos[2] } );
    }

    return pos_vv;
  }
  
  /**
   * @brief Fix the location of where the start starts depositing energy 
   *
   */
  std::vector<float> MCPixelPGraph::fixingPhotonStartPoints( Node_t& node,
							     const std::vector<larcv::Image2D>& instance_v,
							     const std::vector<larcv::Image2D>& ancestor_v,
							     const std::vector<larcv::Image2D>& adc_v,
							     const std::vector<larcv::Image2D>& larflow_v )
  {

    // start of photon shower is not something accurately specified in the mcreco
    // products. detprofile is provided in mcshower -- but i think its basically trash.
    // instead we will use the instanceid + larflow truth to
    // (1) build 3d spacepoints
    // (2) find the closest spacepoint with evidence for ~MIP level deposits
    // (3) place the position there
    // the profile also provides a momentum, but I bet its garbage

    // this is only for photons (electrons have a good start point based on node.start)
    std::vector<float> shower_start_pt;
    
    if ( node.pid!=22 && abs(node.pid)!=11) {
      // we only apply the following to energy deposits made by photons (pid=22) and electrons (pid=11)
      // i.e. those producing EM showers
      return shower_start_pt; 
    }
    LARCV_INFO() << "================================================================" << std::endl;
    LARCV_INFO() << " START " << std::endl;
    
    // we must collect and scan the instance IDs related to the shower
    // int shower_aid = node.aid;
    // int shower_mid = node.mid;
    // int shower_tid = node.tid;

    // created 3d points of energy deposition from truth information saved in
    // the larflow images. The larflow images encode 2D wire-intersections
    // which provide (Y,Z) positions in the detector. The X position is
    // inferred by the drift time relative to the beam trigger (at tick 3200).
    // We only build 3d points coming from pixels associated with the node.
    pointList pos_vv = makeNode3DpointsFromLArFlowTruth( node, adc_v, larflow_v );


    LARCV_INFO() << "  Node[" << node.nodeidx << "] tid=" << node.tid << " pid=" << node.pid << std::endl;
    LARCV_INFO() << "  Number of spacepoints found from using larflow: " << pos_vv.size() << std::endl;
    
    // get the position where the photon was created
    // note: not the same as the location where it first interacted
    std::vector<float> start_reco_pt = MCPos2ImageUtils::Get()->truepos_to_recopos( node.start[0],
										    node.start[1],
										    node.start[2],
										    node.start[3],
										    true, true );
    LARCV_INFO() << "  Find closest cluster to start point: (" << start_reco_pt[0] << ","
	      << start_reco_pt[1] << ", "
	      << start_reco_pt[2] << ", "
	      << start_reco_pt[3] << " usec)"
	      << std::endl;

    // ensure the shower_start_pt vector has 4 components
    // it will store (x,y,z,t)
    shower_start_pt.resize(4);

    // we want to define the "trunk" cluster of the shower.
    // first, convert the information into the form our clustering algorithm needs.
    std::vector< std::vector<float> > data_v;
    // loop over the 3d points
    for ( auto& pt : pos_vv ) {
      //std::cout << "3d pt: (" << pt[0] << ", " << pt[1] << ", " << pt[2] << ") " << std::endl;

      // // calculate distance of 3d point to first energy deposition point
      // float dist_edep = 0.;
      // float dx = 0.;
      // for (int v=0; v<3; v++) {
	    //   dx = pt[v] - node.first_edep_pos[v];
	    //   dist_edep += dx*dx;
      // }
      // dist_edep = sqrt(dist_edep);
      
      // if ( dist_edep>photon_start_edep_radius_cm || pt.src_pixval<photon_start_pixval_threshold )
      // 	continue;

      std::vector<float> xpt(4,0);
      for (int v=0; v<3; v++)
	      xpt[v] = pt[v];
      // add time dim, but put in a useless value.
      xpt[3] = 0.0;
      data_v.push_back( xpt );
      
    }

    // use our interface to DBScan
    auto dbcluster_v = ublarcvapp::dbscan::DBScan::makeCluster3f( 0.3, 3, 50, data_v );
    LARCV_INFO() << "Clustering spacepoints. Num clusters=" << dbcluster_v.size() << std::endl;    

    // loop through the cluster and find the closest one
    // that qualifies in terms of size
    int icluster = -1;
    float closest_qualified_cluster_dist = 1.0e9;
    std::vector<float> closest_pixsum_v;
    for (int i=0; i<(int)dbcluster_v.size()-1; i++) {
      // (we skip the last noise cluster)

      // we calculate the pixel sum in each plane of each cluster.
      // (this is costly -- if it starts to be too slow -- we will need some thresholding with easy measures)
      int nhits = dbcluster_v.at(i).size();
      
      // ignore small clusters
      if ( nhits<5 )
	      continue;

      // find closest distance to photon start/creation point
      float cluster_min_dist_to_source = 1e9;
      std::vector<float> closest_pt(3,0.0);

      // extract the points for the given cluster (with index i)
      pointList cluster_pt_v;
      cluster_pt_v.reserve(nhits);
      for (int ipt=0; ipt<nhits; ipt++) {
	      auto& xpt = data_v.at( dbcluster_v.at(i).at(ipt) );
	      cluster_pt_v.push_back( xpt );

	      // calculate distance to true point
	      float dist = 0;
	      for (int v=0; v<3; v++) {
	        float dx = xpt[v]-start_reco_pt[v];
	        dist += dx*dx;
	      }
	      dist = sqrt(dist);
	      if ( dist < cluster_min_dist_to_source ) {
	        cluster_min_dist_to_source = dist;
		      closest_pt = xpt;
	      }
      }//end of loop over hits
      
      // get the pixel sum of the cluster
      std::vector<float> pixsum_v = getPlanePixelSumsFromPointList( cluster_pt_v, adc_v );
      
      // another threshold: one plane must pass 5 MeV threshol
      // using MeV = 0.00162*pixsum, from matt's conversion formula
      int planes_passing = 0;
      //std::cout << "cluster[" << i << "] dist-to-source=" << cluster_min_dist_to_source << " cm; plane MeV: (";
      for ( int p=0; p<(int)pixsum_v.size(); p++ ) {
	      float MeV = pixsum_v[p]*0.0162;
	      //std::cout << MeV;
	      // if ( p+1<nhits )
	      //   std::cout << ", ";
        // enforce an energy threshold
        // (note: is this too high? especially since we must see that amount in 2 of 3 planes)
	      if ( MeV>5.0 ) {
	        planes_passing++;
	      }
      }
      if (planes_passing>=1 && cluster_min_dist_to_source<100.0) {
        LARCV_INFO() << "cluster[" << i << "] npoints=" << nhits 
              << " dist-to-source=" << cluster_min_dist_to_source << " cm;"
              << " plane MeV: (" << pixsum_v[0]*0.0162 << ", "
              << pixsum_v[1]*0.0162 << ", "
              << pixsum_v[2]*0.0162 << ")"
              << " pt=(" << closest_pt[0] << "," << closest_pt[1] << "," << closest_pt[2] << ")"
              << "nplanes passing: " << planes_passing << std::endl;
      }
      
      // check if cluster passes 'detectability criterion'
      // and, if so, is the closest cluster to the creation point of the photon
      if ( planes_passing>=2 && cluster_min_dist_to_source < closest_qualified_cluster_dist ) {
	      closest_qualified_cluster_dist = cluster_min_dist_to_source;
	      icluster = i;
        closest_pixsum_v = pixsum_v;
      }
    }//end of loop over clusters
    
    // container to store our trunk points
    pointList trunk_pt_v;

    if ( icluster<0 ) {
      // if we did not find a valid cluster,
      // we count this shower as unreconstructable.
      // we define a sentinal value for the starting energy deposition point    
      shower_start_pt = std::vector<float>{ -1.0, -1.0, -1.0 };
      LARCV_DEBUG() << "No valid trunk cluster found." << std::endl;
    } 
    else  {
      // a valid cluster is found. store the points for this cluster
      auto& largest_cluster = dbcluster_v.at(icluster);
      // we also find the min distance between a hit in the cluster
      // and the starting reco point.
      // this 3d point will serve as the true first energy deposit point
      float mindist = 1.0e9;
      for (int ipt=0; ipt<(int)largest_cluster.size(); ipt++) {
	      auto& xpt = data_v.at( largest_cluster.at(ipt) );
	      float dist = 0;
	      for (int v=0; v<3; v++) {
	        float dx = xpt[v]-start_reco_pt[v];
	        dist += dx*dx;
	      }
	      dist = sqrt(dist);
	      if ( dist<mindist ) {
          // if the current closest point, use it to define the photon interaction point
	        for (int v=0; v<3; v++)
	          shower_start_pt[v] = xpt[v];
	        shower_start_pt[3] = xpt[0]/(0.5*larutil::LArProperties::GetME()->DriftVelocity()) + 3200;
	        mindist = dist;
	      }
	      trunk_pt_v.push_back( xpt );
      }//end of loop over largest cluster points (by index)
      LARCV_DEBUG() << "Defining observable trunk using cluster[" << icluster << "] " << std::endl;
      LARCV_DEBUG() << "  npoints=" << largest_cluster.size() 
                    << "  plane MeV: (" << closest_pixsum_v[0]*0.0162 << ", "
                    << closest_pixsum_v[1]*0.0162 << ", "
                    << closest_pixsum_v[2]*0.0162 << ")" 
                    << std::endl;
      LARCV_DEBUG() << "  minimum dist to shower startpt: " << mindist << " cm" << std::endl;
    }//end of else a valid shower cluster found
    
    LARCV_INFO() << "  startpt: (" << shower_start_pt[0] << ","
	      << shower_start_pt[1] << ","
	      << shower_start_pt[2] << ","
	      << shower_start_pt[3] << " ticks)"
	      << std::endl;

    // convert the first interaction point (x,y,z,t) position into 
    // image coordinates which are (u,v,y,tick)
    float t_ns = (shower_start_pt[3]-3200)*0.5*1000.0+3050.0; // ticks to ns
    std::vector<float> shower_imgpos4 = MCPos2ImageUtils::Get()->to_imagepos( shower_start_pt[0],
									      shower_start_pt[1],
									      shower_start_pt[2],
									      t_ns );
    LARCV_INFO() << "  shower imgpos: (" << shower_imgpos4[0] << ", "
	      << shower_imgpos4[1] << ", "
	      << shower_imgpos4[2] << ", "
	      << shower_imgpos4[3] << " ticks )"
	      << std::endl;
    
    // append the imgcoordantes to the out-going vector
    for (int v=0; v<4; v++)
      shower_start_pt.push_back( shower_imgpos4[v] );

    // get an index for this true photon pointList
    // we will use it to store information about this trunk
    // and then retrieve it later.
    int photon_point_list_index = (int)_true_photon_v.size();
    _nodeidx_to_photonlist_index_v[ node.nodeidx ] = photon_point_list_index;

    // making pixellist
    const int nplanes = adc_v.size();
    bool valid_imagemask = true;
    std::vector<float> plane_pixsum_v( nplanes, 0.0 );
    std::vector< PixelSet_t > plane_pixsets_v = _getPlanePixelSetsAndPixelSums( trunk_pt_v, adc_v, plane_pixsum_v );
    if ( plane_pixsum_v.size()>=3 ) {
      LARCV_INFO() << "  pixel sums MeV: (" << plane_pixsum_v[0]*0.0162 << ", "
                  << plane_pixsum_v[1]*0.0162 << ", "
                  << plane_pixsum_v[2]*0.0162 << ") "
                  << std::endl;
    }

    // making image masks
    std::vector< larcv::Image2D > trunk_planeimg_v;
    for (int p=0; p<nplanes; p++) {

      auto const& img = adc_v.at(p);
      auto& pixset = plane_pixsets_v.at(p);
      
      // make an image
      // we also keep track of the min/max for the
      // rows and columns. this defines 
      // a bounding box we need to define the image mask
      // crop bounds.
      int min_row = 999999;
      int max_row = 0;
      int min_col = 999999;
      int max_col = 0;

      // loop over the pixels in order to find axis-aligned bounding box
      for ( auto& pix : pixset ) {
	      int pix_row = pix.first;
	      int pix_col = pix.second;
	      if ( min_row > pix_row )  min_row = pix_row;
	      if ( max_row < pix_row )  max_row = pix_row;
	      if ( min_col > pix_col )  min_col = pix_col;
	      if ( max_col < pix_col )  max_col = pix_col;
      }
      // define the imagemeta to make the cropped larcv::Image2d object
      auto const& meta = img.meta();
      int colcount = max_col-min_col+1;
      int rowcount = max_row-min_row+1;
      if (colcount<=0 or rowcount<=0) {
        // invalid bounding box
        colcount = 1;
        rowcount = 1;
        min_col = 0;
        min_row = 0;
        max_col = 1;
        max_row = 1;
      }
      float width  = float( colcount );
      float height = meta.pixel_height()*rowcount;
      float origin_x = float( meta.pos_x( min_col ) );
      float origin_y = float( meta.pos_y( min_row ) );
	
      LARCV_DEBUG() << "min (col,row)=(" << min_col << "," << min_row << ")" << std::endl;
      LARCV_DEBUG() << "min (x,y)=(" << origin_x << "," << origin_y << ")" << std::endl;      
      larcv::ImageMeta trunk_meta( width, height, rowcount, colcount, origin_x, origin_y, meta.plane() );
      larcv::Image2D trunk_img( trunk_meta );
      trunk_img.paint(0.0);

      LARCV_DEBUG() << "Defined plane bounding box for trunk: \n" << trunk_meta.dump() << std::endl;
      LARCV_DEBUG() << "original meta: \n " << meta.dump() << std::endl;

      // loop over the pixels and fill the mask
      for (auto& pix : pixset ) {
	      int pix_row = pix.first;
	      int pix_col = pix.second;	
	      //LARCV_DEBUG() << " set pixel (row,col)=(" << pix_row << ", " << pix_col << ") " << std::endl;
	      try {
	        trunk_img.set_pixel( pix_row-min_row, pix_col-min_col, 1.0 );
	      }
	      catch (...) {
	        LARCV_ERROR() << "pixel (row,col)=(" << pix_row << "," << pix_col << ") is outside crop." << std::endl;
	      }
      }

      // store the image in the plane container
      trunk_planeimg_v.emplace_back( std::move(trunk_img) );
      
    }//end of plane loop

    // store all of this stuff
    _true_photon_v.emplace_back( std::move(trunk_pt_v) );
    _true_photon_plane_pixset_vv.emplace_back( std::move( plane_pixsets_v ) );
    _true_photon_plane_pixsum_vv.emplace_back( std::move( plane_pixsum_v ) );
    _true_photon_plane_trunkimg_vv.emplace_back( std::move( trunk_planeimg_v ) );

    // make sure we didnt screw something up
    if ( (int)_true_photon_v.size()!=(photon_point_list_index+1) ) {
      LARCV_ERROR() << "index of new true photon 3d point list doesn't match the size of the container" << std::endl;
    }
    if ( _true_photon_v.size()!=_true_photon_plane_pixset_vv.size() )
      LARCV_ERROR() << "size of the true photon 3d point container is not in the pixset container" << std::endl;
    if ( _true_photon_v.size()!=_true_photon_plane_pixsum_vv.size() )
      LARCV_ERROR() << "size of the true photon 3d point container is not in the pixel sum container" << std::endl;
    if ( _true_photon_v.size()!=_true_photon_plane_trunkimg_vv.size() )
      LARCV_ERROR() << "size of the true photon 3d point container is not the same as the trunk image mask container" << std::endl;
    
    //std::cout << "========================================================END======" << std::endl;
    
    return shower_start_pt;
    
  }

  /**
   * @brief Get the 3d points representing the true photon trunk
   *
   */
  const std::vector< std::vector<float> >& MCPixelPGraph::getTruePhotonTrunk3DPoints( Node_t& node )
  {

    auto it_pointlist = _nodeidx_to_photonlist_index_v.find( node.nodeidx );
    if ( it_pointlist==_nodeidx_to_photonlist_index_v.end() ) {
      std::cout << "[MCPixelPGraph::getTruePhotonTrunk3DPoints] [WARNING]: did not find photonlist for "
		            << " nodeidx=" << node.nodeidx
		            << " pdg=" << node.pid
		            << " trackid=" << node.tid << ". "
		            << "Returning emptry point list." 
		            << std::endl;
      _empty_photon_pointlist_v.clear();
      return _empty_photon_pointlist_v;
    }

    size_t index = it_pointlist->second;
    if ((int)index<0 || index>=_true_photon_v.size() ) {
      std::stringstream msg;
      msg << "Missing photon list index! node.trackid=" << node.tid
          << " map.trackid=" << it_pointlist->first
          << " map.index=" << it_pointlist->second
          << std::endl;
      throw std::runtime_error( msg.str() );
    }

    return _true_photon_v.at( it_pointlist->second );
  }

  /**
   * @brief Get the 3d points representing the true photon trunk by Track ID
   *
   */
  const std::vector< std::vector<float> >& MCPixelPGraph::getTruePhotonTrunk3DPoints( int trackid )
  {
    Node_t* pnode = findTrackID( trackid );
    if ( pnode==nullptr ) {
      std::cout << "[MCPixelPGraph::getTruePhotonTrunk3DPoints] [WARNING]: did not particle Node for "
          << " trackid=" << pnode->tid << ". "
          << "Returning emptry point list." 
          << std::endl;
      _empty_photon_pointlist_v.clear();
      return _empty_photon_pointlist_v;
    }

    return getTruePhotonTrunk3DPoints( *pnode );
  }  

  /**
   * @brief Get the sum of pixel values around the trunk plane
   */
  std::vector<float> MCPixelPGraph::getTruePhotonTrunkPlanePixelSums( int trackid )
  {
    // we want to sum using a kernel (3x3) pixels. We want to avoid double counting,
    // so we must mark which pixels we used
    Node_t* pnode = findTrackID( trackid );
    if ( pnode==nullptr ) {
      LARCV_ERROR() << ": did not find particle Node pointer for "
          << " trackid=" << pnode->tid << ". "
          << "Returning emptry point list." 
          << std::endl;
      _empty_photon_pointlist_v.clear();
      return _empty_pixsum_v;
    }
    
    auto it_pixsum = _nodeidx_to_photonlist_index_v.find( pnode->nodeidx );
    if ( it_pixsum==_nodeidx_to_photonlist_index_v.end() ) {
      LARCV_NORMAL() << "Did not find true photon trunk info for trackid=" << trackid << std::endl;
      for ( auto& it_x : _nodeidx_to_photonlist_index_v ) {
        pnode = &(node_v.at( it_x.first ));
        LARCV_NORMAL() << "  [nodeidx=" << it_x.first << ", trackid=" << pnode->tid << "] photonlist index=" << it_x.second << std::endl;
      }
      LARCV_ERROR() << "Stopping on error" << std::endl;
    }

    std::vector<float> pixsum_v = _true_photon_plane_pixsum_vv.at( it_pixsum->second );
    return pixsum_v;
    
  }// end of MCPixelPGraph::getTruePhotonTrunkPlanePixelSums( int trackid )

  /**
   * @brief Get the sum of pixel values around the trunk plane
   */
  const std::vector< MCPixelPGraph::PixelSet_t >&
  MCPixelPGraph::getTruePhotonTrunkPlanePixelSets( int trackid )
  {
    // we want to sum using a kernel (3x3) pixels. We want to avoid double counting,
    // so we must mark which pixels we used
    Node_t* pnode = findTrackID( trackid );
    if ( pnode==nullptr ) {
      std::cout << "[MCPixelPGraph::getTruePhotonTrunk3DPoints] [WARNING]: did not find particle Node pointer for "
        << " trackid=" << pnode->tid << ". "
        << "Returning emptry point list." 
        << std::endl;
      _empty_pixelset_v.clear();
      return _empty_pixelset_v;
    }
    
    auto it_pixsum = _nodeidx_to_photonlist_index_v.find( pnode->nodeidx );
    if ( it_pixsum==_nodeidx_to_photonlist_index_v.end() ) {
      LARCV_ERROR() << "Did not find true photon trunk info for trackid=" << trackid << std::endl;
    }
    
    return _true_photon_plane_pixset_vv.at( it_pixsum->second );
  }
  
  /**
   * @brief Get larcv image containing trunk pixels
   */
  const std::vector< larcv::Image2D >&
  MCPixelPGraph::getTruePhotonTrunkPlaneImage2DMasks( int trackid )
  {
    
    Node_t* pnode = findTrackID( trackid );
    if ( pnode==nullptr ) {
      LARCV_WARNING() << "[MCPixelPGraph::getTruePhotonTrunkPlanePixelSets] "
		      << " [WARNING]: did not find particle Node pointer for "
		      << " trackid=" << pnode->tid << ". "
		      << "Returning emptry point list." 
		      << std::endl;
      _empty_imagemask_v.clear();
      return _empty_imagemask_v;
    }
    
    auto it_pixsum = _nodeidx_to_photonlist_index_v.find( pnode->nodeidx );
    if ( it_pixsum==_nodeidx_to_photonlist_index_v.end() ) {
      LARCV_ERROR() << "Did not find true photon trunk info for trackid=" << trackid << std::endl;
    }

    return _true_photon_plane_trunkimg_vv.at( it_pixsum->second );
    
  }

  /**
   * @brief Does the calculation to get the pixels in each plane corresponding to 3D points of true shower trunk
   */
  std::vector<float>
  MCPixelPGraph::getPlanePixelSumsFromPointList( const std::vector< std::vector<float> >& pointlist,
						 const std::vector< larcv::Image2D >& adc_v )
  {
    std::vector<float> pixsum_v;
      
    const int nplanes = adc_v.size();
    const int dpix = 1;
    std::vector<int> out_of_imgpix(nplanes,0);
    // allocate space to the pixsum
    pixsum_v.resize(nplanes,0);
    
    for (int p=0; p<nplanes; p++) {
      
      auto const& img = adc_v.at(p);

      // this lists (row,col) combinations we gathered pixels from
      std::set< std::pair<int,int> > _pix_used;
      
      // loop over 3d points
      for (int ipt=0; ipt<(int)pointlist.size(); ipt++) {
        auto& pt = pointlist.at(ipt);
        // get the point in the image
        std::vector<float> imgpos4 = MCPos2ImageUtils::Get()->to_imagepos( pt[0],
                          pt[1],
                          pt[2],
                          0.0 );

	      //std::cout << "  pt[" << ipt << "] imgpos4: (" << imgpos4[0] << ", " << imgpos4[1] << ", " << imgpos4[2] << ", " << imgpos4[3] << ")" << std::endl;
	
        // sum over 3x3 kernel
        // first find center pixel
        int row_center = 0;
        int col_center = 0;
        try {
          row_center = img.meta().row( imgpos4[3] );
          col_center = img.meta().col( imgpos4[p] );
        }
        catch (...) {
          // out of image tick
          out_of_imgpix[p]++;
          continue;
        }
	
        for (int dr=-dpix; dr<=dpix; dr++) {
          for (int dc=-dpix; dc<=dpix; dc++) {
            int r = row_center + dr;
            int c = col_center + dc;

            auto it_pix = _pix_used.find( std::pair<int,int>(r,c) );
            if ( it_pix==_pix_used.end() ) {
              // did not find it in the set, so add it as part of the sum
              float pixval = 0.0;
              try {
                pixval = img.pixel( r, c );	      
                _pix_used.insert( std::pair<int,int>(r,c) );
                pixsum_v[p] += pixval;
              }
              catch (...) {
                // probably out of image
                out_of_imgpix[p]++;
                continue;
              }
            }
        
          }//end of dc loop
        }//end of dr loop
	  
      }//end of point loop
      
      //std::cout << "PixelSum[" << p << "]: sum=" << pixelsum_v[p] << "  npixels=" << _pix_used.size() << " out_of_img=" << out_of_imgpix[p] << std::endl;
    }//end of plane loop

    return pixsum_v;
  }// end of MCPixelPGraph::getPlanePixelSumsFromPointList

  /**
   * @brief Get set of pixels in each plane from a 3D point list
   *
   */
  std::vector< MCPixelPGraph::PixelSet_t >
  MCPixelPGraph::_getPlanePixelSetsAndPixelSums( const MCPixelPGraph::pointList& pt_v,
						 const std::vector<larcv::Image2D>& adc_v,
						 std::vector<float>& pixsum_v )
  {
						 

    std::vector< MCPixelPGraph::PixelSet_t > plane_pixsets_v;
    const int nplanes = adc_v.size();
    const int dpix = 2;

    std::vector<int> out_of_imgpix(nplanes,0);
    pixsum_v.resize(nplanes,0.0);    
    for (int p=0; p<nplanes; p++)
      pixsum_v[p] = 0.0;

    for (int p=0; p<nplanes; p++) {

      auto const& img = adc_v.at(p);

      // this lists (row,col) combinations we gathered pixels from
      PixelSet_t _pix_used;
      
      // loop over 3d points
      for (int ipt=0; ipt<(int)pt_v.size(); ipt++) {
        auto& pt = pt_v.at(ipt);
        // get the point in the image
        std::vector<float> imgpos4 = MCPos2ImageUtils::Get()->to_imagepos( pt[0],
                          pt[1],
                          pt[2],
                          0.0 );
        //std::cout << "  pt[" << ipt << "] imgpos4: (" << imgpos4[0] << ", " << imgpos4[1] << ", " << imgpos4[2] << ", " << imgpos4[3] << ")" << std::endl;
        
        // sum over 3x3 kernel
        // first find center pixel
        int row_center = 0;
        int col_center = 0;
        try {
          row_center = img.meta().row( imgpos4[3] );
          col_center = img.meta().col( imgpos4[p] );
        }
        catch (...) {
          // out of image tick
          out_of_imgpix[p]++;
          continue;
        }
      
        for (int dr=-dpix; dr<=dpix; dr++) {
          for (int dc=-dpix; dc<=dpix; dc++) {
            int r = row_center + dr;
            int c = col_center + dc;

            auto it_pix = _pix_used.find( std::pair<int,int>(r,c) );
            if ( it_pix==_pix_used.end() ) {
              // did not find it in the set, so add it as part of the sum
              float pixval = 0.0;
              try {
                pixval = img.pixel( r, c );	      
                _pix_used.insert( std::pair<int,int>(r,c) );
                pixsum_v[p] += pixval;
              }
              catch (...) {
                // probably out of image
                out_of_imgpix[p]++;
                continue;
              }
            }
	        }//end of dc loop
	      }//end of dr loop
      }//end of point loop 

      //std::cout << "PixelSum[" << p << "]: sum=" << pixelsum_v[p] << "  npixels=" << _pix_used.size() << " out_of_img=" << out_of_imgpix[p] << std::endl;
      plane_pixsets_v.emplace_back( std::move( _pix_used ) );   
      
    }//end of plane loop

    return plane_pixsets_v;
  }
  
  /**
  * @brief Get the position the particle first deposited energy
  *
  */ 
  std::vector<float> MCPixelPGraph::getParticleEDepPos( const int& trackid )
  {
    std::vector<float> edep_pos;
    Node_t* pnode = findTrackID( trackid );
    if (pnode==nullptr) {
      std::stringstream msg;
      msg << "[MCPixelPGraph::getParticleEDepPos] [ERROR] : Could not find a particle node for trackid=" << trackid << std::endl;
      throw std::runtime_error( msg.str() );
    }
    edep_pos.resize(3,0);
    for (int v=0; v<3; v++) {
      edep_pos[v] = pnode->first_edep_pos[v];
    }
    return edep_pos;
  }
  
  /**
   * @brief get ancestor ID for track ID by querying info in Nodes
   */
  int MCPixelPGraph::getAncestorID( int trackid ) {
    Node_t* pnode = findTrackID(trackid);
    if (pnode==nullptr)
      return -1;
    return pnode->aid;
  }

  /**
   * @brief get ancestor ID for track ID by querying info in Nodes
   */
  int MCPixelPGraph::getShowerMotherID( int trackid ) {
    auto it = _shower_daughter2mother.find( trackid );
    if ( it!=_shower_daughter2mother.end() ) {
      return it->second;
    }
    return -1;
  }

  /**
   * @brief get particle ID for track ID by querying info in Nodes
   */
  int MCPixelPGraph::getParticleID( int trackid ) {
    Node_t* pnode = findTrackID(trackid);
    if (pnode==nullptr)
      return -1;
    return pnode->pid;
  }


  
}
}
