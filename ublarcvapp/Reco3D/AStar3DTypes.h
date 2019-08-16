#ifndef __ASTAR3D_TYPES_H__
#define __ASTAR3D_TYPES_H__

#include <iostream>
#include <vector>
#include <array>
#include <sstream>
#include <map>
#include <algorithm>

namespace ublarcvapp {
namespace reco3d {
  
  // LATTICE POSITION OBJECT
  // stores the 3D position of a lattice point
  class A3DPixPos_t : public std::array<int,3> {
  public:
    A3DPixPos_t() : std::array<int,3>{{-1,-1,-1}}
    {};
  A3DPixPos_t(int u, int v, int w) : std::array<int,3>{{u,v,w}}
    {};
    virtual ~A3DPixPos_t() {};
  };

  // NODE DEFINITION
  // associate a point on our lattice to a 3D position
  // the node also stores scores and status (open/close)
  class AStar3DNode  {
  public:
    AStar3DNode() {
      u = v = w = 0;
      fscore=gscore=kscore=0.0;
      closed = false;
      inopenset = false;
      dir3d.resize(3,0.0);
      tyz.resize(3,0.0);
      prev = NULL;
      row = 0;
      within_image=false;
      pixval.resize(3,0.0);
      badchnode = false;
    };
    AStar3DNode( int u_, int v_, int w_, std::vector<float> tyz_ ) { 
      u = u_; 
      v = v_;
      w = w_;
      nodeid[0] = u;
      nodeid[1] = v;
      nodeid[2] = w;
      fscore=0.0;
      gscore=0.0;
      kscore=0.0;
      within_image=false;
      dir3d.resize(3,0);
      tyz = tyz_;
      closed = false;
      inopenset = false;
      prev = NULL;
      row = 0;
      pixval.resize(3,0.0);
      badchnode = false;
    };
    AStar3DNode( const AStar3DNode& src ) {
      u = src.u;
      v = src.v;
      w = src.w;
      nodeid = src.nodeid;
      fscore = src.fscore;
      gscore = src.gscore;
      kscore = src.kscore;
      dir3d = src.dir3d;
      tyz = src.tyz;
      closed = src.closed;
      within_image = src.within_image;
      inopenset = src.inopenset;
      prev = src.prev;
      row = src.row;
      cols = src.cols;
      pixval = src.pixval;
      badchnode = src.badchnode;
    };
    virtual ~AStar3DNode() {};

    int u; // lattice position corresponding to time or x
    int v; // lattice position corresponding to y
    int w; // latrice position corresponding to z
    A3DPixPos_t nodeid; // latice position key
    float fscore; // score of path: past path + to goal heuristic
    float gscore; // past path score
    float kscore; // curvature penalty

    std::vector<float> dir3d;
    std::vector<float> tyz; // (tick, y, z) in detector coordinates
    int row;
    std::vector<int> cols;
    std::vector<float> pixval;
    bool within_image;
    bool inopenset;
    bool closed;
    bool badchnode; // tag if its a badch node
    AStar3DNode* prev;

    bool operator <(const AStar3DNode& rhs ) const {
      //std::cout << "astar< " << fscore << " < " << rhs.fscore << std::endl;
      if ( fscore > rhs.fscore ) // lower scores get priority
        return true;
      return false;
    };
    bool operator==(const AStar3DNode& rhs ) const {
      if ( u==rhs.u && v==rhs.v && w==rhs.w ) return true;
      return false;
    };

    bool operator!=(const AStar3DNode& rhs ) const {
      if ( u!=rhs.u || v!=rhs.v || w!=rhs.w ) return true;
      return false;
    };

    std::string str() const { 
      std::stringstream ss;
      ss << "[node=(" << u << "," << v << "," << w << ") "
      	 << "tyz=(" << tyz[0] << "," << tyz[1] << "," << tyz[2] << ") "
      	 << "f=" << fscore << " "
      	 << "g=" << gscore << " "
      	 << "k=" << kscore << " "
      	 << "]";
      return ss.str();
    }

  };

  // NODE CONTAINER
  class AStar3DNodePtrList : public std::vector< AStar3DNode* > {

    public:
      AStar3DNodePtrList() {};
      virtual ~AStar3DNodePtrList() {};

    protected:

    // we want to sort from largest to smallest, because we will pop the nodes off the back of the various sets
    struct dirnodeptr_compare_t {
      bool operator() (const AStar3DNode* lhs, const AStar3DNode* rhs ) const {
        if ( lhs->fscore > rhs->fscore )
          return true;
        return false;
      };
    } m_thecomparator;

    struct dirnodeptr_compare_hscore_t {
      bool operator() (const AStar3DNode* lhs, const AStar3DNode* rhs ) const {
        if ( lhs->fscore-lhs->gscore > rhs->fscore-rhs->gscore )
          return true;
        return false;
      };
    } m_the_hscore_comparator;

    public:
    void sort() {
      std::sort(this->begin(),this->end(), m_thecomparator);
    }
    void sort_by_hscore() {
      std::sort(this->begin(),this->end(), m_the_hscore_comparator);
    }    

    void addnode( AStar3DNode* elem ) {
    	elem->inopenset = true;
      push_back( elem );
      //sort();
    }

  };
  
  // A PATH THROUGH THE LATTICE
  typedef std::vector< std::vector<AStar3DNode> > PathList_t;      

  //LATTICE DEFINITION
  typedef std::pair< A3DPixPos_t, AStar3DNode* > a3dpos_pair_t;
  
  
}
}

#endif
