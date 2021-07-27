#ifndef __DWALL_H__
#define __DWALL_H__

#include <vector>

namespace ublarcvapp {
  
  float  dwall( const std::vector<float>& pos, int& boundary_type );
  double dwall( const std::vector<double>& pos, int& boundary_type );

  float  dspecificwall( const std::vector<float>& pos, const int boundary_type );
  double dspecificwall( const std::vector<double>& pos, const int boundary_type );

  float  dwall_noAC( const std::vector<float>& pos, int& boundary_type );
  double dwall_noAC( const std::vector<double>& pos, int& boundary_type );  

  /**
   * @brief class to provide python binding to dwall functions
   *
   */
  class pydwall {

  public:   
    pydwall() {};
    virtual ~pydwall() {};
    static float  dwall( const float x, const float y, const float z );
  };
  
}

#endif
