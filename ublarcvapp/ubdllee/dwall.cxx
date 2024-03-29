#include "dwall.h"

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

    float dx1 = pos[0];
    float dx2 = 258-pos[0];
    float dy1 = 117.0-pos[1];
    float dy2 = pos[1]+117.0;
    float dz1 = pos[2];
    float dz2 = 1036.0-pos[2];

    float dwall = 1.0e9;

    if ( dy1<dwall ) {
      dwall = dy1;
      boundary_type = 0; // top
    }
    if ( dy2<dwall ) {
      dwall = dy2;
      boundary_type = 1; // bottom
    }
    if ( dz1<dwall ) {
      dwall = dz1;
      boundary_type = 2; // upstream
    }
    if ( dz2<dwall ) {
      dwall = dz2;
      boundary_type = 3; // downstream
    }
    if ( dx1<dwall ) {
      dwall = dx1;
      boundary_type = 4; // anode
    }
    if ( dx2<dwall ) {
      dwall = dx2;
      boundary_type = 5; // cathode
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

    float dx1 = pos[0];
    float dx2 = 258-pos[0];
    float dy1 = 117.0-pos[1];
    float dy2 = pos[1]+117.0;
    float dz1 = pos[2];
    float dz2 = 1036.0-pos[2];

    float fdwall = 1.0e9;

    switch ( boundary_type ) {
    case 0:
      fdwall = dy1;
      break;
    case 1:
      fdwall = dy2;
      break;
    case 2:
      fdwall = dz1;
      break;
    case 3:
      fdwall = dz2;
      break;
    case 4:
      fdwall = dx1;
      break;
    case 5:
      fdwall = dx2;
      break;
    default:
      throw std::runtime_error("dspecifcwall: imageend boundary points undefined");
    }

    return fdwall;
    
  }

  float dwall_noAC( const std::vector<float>& pos, int& boundary_type ) {

    float dy1 = 117.0-pos[1];
    float dy2 = pos[1]+117.0;
    float dz1 = pos[2];
    float dz2 = 1036.0-pos[2];

    float dwall = 1.0e9;

    if ( dy1<dwall ) {
      dwall = dy1;
      boundary_type = 0; // top
    }
    if ( dy2<dwall ) {
      dwall = dy2;
      boundary_type = 1; // bottom
    }
    if ( dz1<dwall ) {
      dwall = dz1;
      boundary_type = 2; // upstream
    }
    if ( dz2<dwall ) {
      dwall = dz2;
      boundary_type = 3; // downstream
    }

    return dwall;
    
  }

  double dwall_noAC( const std::vector<double>& pos, int& boundary_type ) {
    std::vector<float> fpos(3);
    for (int i=0; i<3; i++)
      fpos[i] = (float)pos[i];
    return (double)dwall_noAC( fpos, boundary_type);
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
  float  pydwall::dwall_noAC( const float x, const float y, const float z ) {
    std::vector<float> pos = { x, y, z };
    int boundary = 0;
    return (float)ublarcvapp::dwall_noAC( pos, boundary );
  }
  
}
