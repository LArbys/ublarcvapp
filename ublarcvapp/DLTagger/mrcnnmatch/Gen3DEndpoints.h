#ifndef __GEN_3D_ENDPOINTS_H__
#define __GEN_3D_ENDPOINTS_H__

#include <vector>
#include "FeaturesMaskCombo.h"

namespace ublarcvapp {
namespace dltagger {

  class Gen3DEndpoints {

  public:

    typedef std::vector<float> EndPoint_t;
    
    Gen3DEndpoints()
      : pfeatures(nullptr)
    {}; ///< empty constructor. avoid using. only here for ROOT dictonary making
    virtual ~Gen3DEndpoints() {}; ///< destructor

    Gen3DEndpoints( const FeaturesMaskCombo& cropdata );

    std::vector< std::vector< EndPoint_t > > mask_endpoints_vv; // a list of 2d points (min,max) for each plane
    std::vector< std::vector< float > >      mask_projdist_vv;  // projection onto PCA1 of the endpoints for each plane

    std::vector< std::vector<float> > endpt_tyz_v; // for each end point, the (tick, y, z ) coordinates
    std::vector< std::vector<int> >   endpt_wid_v; // for each end point, the wires used on each plane
    std::vector< float >              endpt_tri_v; // for each end point, a measure of how consistent the end point was
    std::vector< int >                endpt_tpc_v; // for each end point, flag indicating if wires intersect inside TPC

    const FeaturesMaskCombo* pfeatures; ///< holds data we use to make this output
    
  protected:

    void _gather_plane_endpoints( const FeaturesMaskCombo& features );
    void _gen_3dendpoints();

  };
  
}
}

#endif
