#include "Lattice.h"

// larlite
#include "larlite/LArUtil/Geometry.h"
#include "larlite/LArUtil/LArProperties.h"

namespace ublarcvapp {
namespace reco3d {

  AStar3DNode* Lattice::getNode( const A3DPixPos_t& nodeid ) {

    // check bounds
    for (int i=0; i<3; i++) {
      if ( nodeid[i]<0 || nodeid[i]>=m_widths[i]) return nullptr;
    }

    // search map
    auto it_node = find( nodeid );

    AStar3DNode* node = nullptr;

    if ( it_node==end() ) {
      // [[ create a node ]]
      // get the detector coordinates (tick, y, z)
      std::vector<float> tyz = getPos( nodeid );
      // create the node
      node = new AStar3DNode( nodeid[0], nodeid[1], nodeid[2], tyz );
      // get the image coordinates
      getImageCoordinates( nodeid, node->row, node->cols, node->within_image );
      // put it into the map
      insert( a3dpos_pair_t(nodeid,node) );
    }
    else {
      node = (*it_node).second;
    }

    return node;
  }

  AStar3DNode* Lattice::getNode( const int u, const int v, const int w ) {

    // check for node in the map
    A3DPixPos_t nodeid(u,v,w);
    return getNode( nodeid );
  }



  A3DPixPos_t Lattice::getNodePos( const std::vector<float>& pos ) {
    A3DPixPos_t nodeid(0,0,0);
    for (int i=0; i<3; i++) {
      nodeid[i] = (int) (( pos[i]-m_origin[i])/m_cm_per_pixel[i]);
      if ( nodeid[i]<0 || nodeid[i]>=m_widths[i])
        return A3DPixPos_t(); // outside lattice, so send an empty one
    }
    return nodeid;
  }


  AStar3DNode* Lattice::getNode( const std::vector<float>& pos ) {
    return getNode( getNodePos(pos) );
  }

  std::vector<float> Lattice::getPos( const int u, const int v, const int w ) {
    A3DPixPos_t nodeid( u, v, w );
    return getPos( nodeid );
  }

  std::vector<float> Lattice::getPos( const A3DPixPos_t& nodeid ) {
    std::vector<float> pos(3,0.0);
    for (int i=0; i<3; i++) {
      // check if within lattice?
      pos[i] = m_origin[i] + nodeid[i]*m_cm_per_pixel[i];
    }
    return pos;
  }

  void Lattice::getImageCoordinates( const A3DPixPos_t& nodeid, int& row, std::vector<int>& cols, bool& within_image ) {
    // turn this lattice point into a position in the images
    row = 0;
    cols.resize(3,0);

    std::vector<float> tyz = getPos( nodeid );

    if ( tyz[0]<=m_meta_v.front()->min_y() || tyz[0]>=m_meta_v.front()->max_y() ) {
      within_image = false;
      return;
    }

    std::vector<float> wid(3,0);
    cols.resize(3,0);
    within_image = true;
    for ( int p=0; p<3; p++ ) {
      double xyz[3] { (tyz[0]-3200)*0.5*::larutil::LArProperties::GetME()->DriftVelocity(), tyz[1], tyz[2] }; // (x doesn't matter really)
      wid[p] = larutil::Geometry::GetME()->WireCoordinate( xyz, p );
      if ( wid[p]<=m_meta_v.at(p)->min_x() || wid[p]>=m_meta_v.at(p)->max_x() ) {
        within_image = false;
        //std::cout << "WARNING: sampled wire=" << wid[p] << " on plane=" << p << " outside of image" << std::endl;
        return;
      }
      cols[p] = m_meta_v.at(p)->col(wid[p]);
    }
    row = m_meta_v.front()->row(tyz[0]);
    return;
  }

  void Lattice::cleanup () {

    //if ( verbose>0 )
    //std::cout << "Lattice::cleanup";

    for ( auto &it : *this ) {
      delete it.second;
      it.second = nullptr;
    }
    //if ( verbose>0 )
    //std::cout << "... done." << std::endl;
    clear();

  }
  
  
}
}
