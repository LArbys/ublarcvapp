#include "dwall.h"
#include "larlite/LArUtil/Geometry.h"

#include <cmath>
#include <stdexcept>

namespace ublarcvapp {

  double dwall( const std::vector<double>& pos, int& boundary_type ) {
    std::vector<float> fpos(3);
    for (int i=0; i<3; i++)
      fpos[i] = (float)pos[i];
    return dwall( fpos, boundary_type);
  }
  
  float dwall( const std::vector<float>& pos, int& boundary_type ) {

    auto const geom = larlite::larutil::Geometry::GetME();
    
    TVector3 vpos( pos[0], pos[1], pos[2] );
    std::vector<int> ct = geom->GetContainingCryoAndTPCIDs( vpos );
    int tpcid  = ct[1];
    int cryoid = ct[0];

    auto const& tpcgeo   = geom->GetTPC( tpcid, cryoid );
    TVector3 tpcdriftdir = geom->TPCDriftDir( tpcid, cryoid );

    float dedge[3][2] = {0};
    int minboundary = -1;
    int mindim    = -1;
    float minboundary_dist = 1e9;
    float dwall = 1.0e9;    
    for (int v=0; v<3; v++) {
      dedge[v][0] = vpos[v]-tpcgeo.fBounds[0][v];   // dist to low bound
      dedge[v][1] = tpcgeo.fBounds[1][v] - vpos[v]; // dist to high bound
      for (int ibound=0; ibound<2; ibound++) {
	float dist = fabs(dedge[v][ibound]);
	if ( dist<minboundary_dist ) {
	  minboundary_dist = dist;
	  mindim = v;
	  minboundary = ibound;
	  dwall = dedge[v][ibound];
	}
      }
    }
    

    
    if ( mindim==1 ) {
      if (minboundary==1)
	boundary_type = 0; // top
      else
	boundary_type = 1; // botom
    }
    else if ( mindim==2 ) {
      if (minboundary==0)
	boundary_type = 2; // upstream
      else
	boundary_type = 3; // downstream
    }
    else {
      // X-direction: or dirft dir
      if ( tpcdriftdir[0]>0 ) {
	// anode is high-bound
	if ( minboundary==1 )
	  boundary_type = 4; // anode
	else
	  boundary_type = 5; // cathode
      }
      else {
	// anode is low-bound
	if ( minboundary==0 )
	  boundary_type = 4; // anode
	else
	  boundary_type = 5; // cathode
      }
    }
    
    return dwall;
    
  }

  double dspecificwall( const std::vector<double>& pos, const int boundary_type ) {
    std::vector<float> fpos(3);
    for (int i=0; i<3; i++)
      fpos[i] = (float)pos[i];
    return dspecificwall( fpos, boundary_type);
  }
  
  float dspecificwall( const std::vector<float>& pos, const int boundary_type ) {

    auto const geom = larlite::larutil::Geometry::GetME();
    
    TVector3 vpos( pos[0], pos[1], pos[2] );
    std::vector<int> ct = geom->GetContainingCryoAndTPCIDs( vpos );
    int tpcid  = ct[1];
    int cryoid = ct[0];

    auto const& tpcgeo   = geom->GetTPC( tpcid, cryoid );
    TVector3 tpcdriftdir = geom->TPCDriftDir( tpcid, cryoid );

    int dim   = -1;
    int bound = -1;
    switch ( boundary_type ) {
    case 0://top
      dim = 1;
      bound = 1;
      break;
    case 1: // bottom
      dim = 1;
      bound = 0;
      break;
    case 2: // upstream
      dim = 2;
      bound = 0;
      break;
    case 3: // downstream
      dim = 2;
      bound = 1;
      break;
    case 4: // anode (assuming -1 drift dir)
      dim = 1;
      bound = (tpcdriftdir[0]<0) ? 0 : 1;
      break;
    case 5:
      dim = 1;
      bound = (tpcdriftdir[0]<0) ? 1 : 0;
      break;
    default:
      throw std::runtime_error("dspecifcwall: imageend boundary points undefined");
    }
    
    float fdwall = 1e9;
    if ( bound==0 )
      fdwall = vpos[dim]-tpcgeo.fBounds[0][dim];   // dist to low bound
    else
      fdwall = tpcgeo.fBounds[1][dim]-vpos[dim];

    return fdwall;
    
  }

  float dwall_noAC( const std::vector<float>& pos, int& boundary_type, int tpcid, int cryoid ) {

    auto const geom = larlite::larutil::Geometry::GetME();
    
    TVector3 vpos( pos[0], pos[1], pos[2] );

    auto const& tpcgeo   = geom->GetTPC( tpcid, cryoid );

    float dedge[3][2] = {0};
    int minboundary = -1;
    int mindim    = -1;
    float minboundary_dist = 1e9;
    float dwall = 1.0e9;    
    for (int v=1; v<3; v++) {
      dedge[v][0] = vpos[v]-tpcgeo.fBounds[0][v];   // dist to low bound
      dedge[v][1] = tpcgeo.fBounds[1][v] - vpos[v]; // dist to high bound
      for (int ibound=0; ibound<2; ibound++) {
	float dist = fabs(dedge[v][ibound]);
	if ( dist<minboundary_dist ) {
	  minboundary_dist = dist;
	  mindim = v;
	  minboundary = ibound;
	  dwall = dedge[v][ibound];
	}
      }
    }
    
    if ( mindim==1 ) {
      if (minboundary==1)
	boundary_type = 0; // top
      else
	boundary_type = 1; // botom
    }
    else if ( mindim==2 ) {
      if (minboundary==0)
	boundary_type = 2; // upstream
      else
	boundary_type = 3; // downstream
    }
    else {
      throw std::runtime_error("dwall_noAC: unexpected dimension with the min distance");
    }
    
    return dwall;
    
  }

  double dwall_noAC( const std::vector<double>& pos, int& boundary_type, int tpcid, int cryoid ) {
    std::vector<float> fpos(3);
    for (int i=0; i<3; i++)
      fpos[i] = (float)pos[i];
    return (double)dwall_noAC( fpos, boundary_type, tpcid, cryoid);
  }

  /**
   * @brief a python-friendly version of dwall
   *
   */
  float  pydwall::dwall( const float x, const float y, const float z ) {
    std::vector<float> pos = { x, y, z };
    int boundary = 0;
    return (float)ublarcvapp::dwall( pos, boundary );
  }

  /**
   * @brief a python-friendly version of dwall_noAC
   *
   */
  float  pydwall::dwall_noAC( const float x, const float y, const float z, int tpcid, int cryoid ) {
    std::vector<float> pos = { x, y, z };
    int boundary = 0;
    return (float)ublarcvapp::dwall_noAC( pos, boundary, tpcid, cryoid );
  }
  
}
