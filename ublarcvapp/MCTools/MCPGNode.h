#ifndef __UBLARCVAPP_MCTOOLS_MCPGNODE_H__
#define __UBLARCVAPP_MCTOOLS_MCPGNODE_H__

#include <vector>
#include <string>

namespace ublarcvapp {
namespace mctools {

    class MCPGNode {

    public:

      int nodeidx;      // book-keeping index
      int type;         // track=0, shower=1, nu-vertex=2, genie_fs=3
      int vidx;         // position in mcshower or mctrack vector
      int tid;          // geant4 track-ID
      int aid;          // ancestor geant4 trackid
      int mtid;         // mother geant4 trackid
      int pid;          // particle ID
      MCPGNode* mother; // pointer to Mother Node
      int mid;          // mother nodeidx
      float E_MeV;      // energy
      std::string process;    // creating process
      std::string mother_process; // mother process
      std::string ancestor_process; // ancestor process
      std::vector<MCPGNode*>  daughter_v; // pointer to daughters 
      std::vector<int>      daughter_idx_v; // daughter node indices in node_v
      std::vector<float> start;          //< (x,y,z,t) before sce, true start of particle
      std::vector<float> first_edep_pos; //< (x,y,z,t) before sce, first step that leaves edep in cryostat
      std::vector<float> first_tpc_pos;  //< (x,y,z,t) before sce, first step inside the TPC, visible in the image
    
      int origin; // 1=neutrino, 2=cosmic, 0=unassigned, -1=unassigned

      MCPGNode()
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
        origin(-1)
      {
	      daughter_v.clear();
        daughter_idx_v.clear();
      };
        
      MCPGNode(int _nodeidx, int _type, int _tid, int _vidx,
	     int _pid,
	     MCPGNode* _mother=nullptr,
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
        origin(-1)
      {
	      daughter_v.clear();
        daughter_idx_v.clear();
      };

      bool operator<( const MCPGNode& rhs ) const {
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

}
}

#endif