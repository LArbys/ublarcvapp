#ifndef __GEN_GRAPH_POINTS_H__
#define __GEN_GRAPH_POINTS_H__

#include <vector>

#include "larcv/core/Base/larcv_base.h"
#include "FeaturesMaskCombo.h"
#include "Gen3DEndpoints.h"

#include <Eigen/Dense>

namespace ublarcvapp {
namespace dltagger {

  class GenGraphPoints : public larcv::larcv_base {

  public:

    GenGraphPoints()
      : larcv::larcv_base("GenGraphPoints"),
      pfeatures (nullptr)
    {};
    virtual ~GenGraphPoints() {};

    GenGraphPoints( const FeaturesMaskCombo& featuredata,
                    const Gen3DEndpoints& endpointdata,
                    larcv::msg::Level_t msglevel=larcv::msg::kNORMAL );

    const FeaturesMaskCombo* pfeatures;

    struct MaskExtrema_t {

      float bounds[4]; // [minx,maxx][miny,maxy]
      float points[4][2]; // 
      MaskExtrema_t() {
        bounds[0] = 1.0e9; // minx
        bounds[1] = 0;     // maxx
        bounds[2] = 1.0e9; // miny
        bounds[3] = 0;     // maxy
      };        
      
    };

    void DefineBoundPointsWithPCAEndpoints( const Gen3DEndpoints& endpointdata,
                                            const float maxstepsize,
                                            const int pixradius,
                                            std::vector< std::vector<float> >& good_twid_v,
                                            std::vector< std::vector<float> >& good_xyz_v,                                            
                                            std::vector< int >& is_good_pt_v );
    
    
    void DefineBoundPointsWithContourMeta( const FeaturesMaskCombo& features,
                                           std::vector< MaskExtrema_t >& maskbounds_v );

    void makePointsFixedZ( const FeaturesMaskCombo& features,
                           const float detz, const int plane, const float tick,
                           std::vector< std::vector<float> >& wiretickpos_v );

    void scanTickDim( const FeaturesMaskCombo& features,
                      const float wireco, const int plane, 
                      std::vector< std::vector<float> >& wiretickpos_v );
    

    void make3Dpoints( const float ywire,
                       const float tick,
                       const std::vector< std::vector<float> >& wiretick_p0_v,
                       const std::vector< std::vector<float> >& wiretick_p1_v,
                       std::vector< std::vector<float> >& points3d_v,
                       std::vector< std::vector<float> >& twid_v,
                       bool allow_2plane_matches );

    void makeGraph( const std::vector< std::vector<float> >& points_tyz,  // input list of (tick,y,z) points
                    std::vector< Eigen::Vector3f >& points,               // node list: position of (x,y,z)
                    std::map< std::pair<int,int>, float >& distmap,       // distance between nodes
                    std::map< std::pair<int,int>, float >& pixgapdist );  // pixel gap between vertices

    float get_max_pixelgap( const std::vector<larcv::Image2D>& adcimg,
                            const std::vector<larcv::Image2D>& badch,
                            const Eigen::Vector3f& start,
                            const Eigen::Vector3f& end,
                            const float max_steplen,
                            const int pixel_search_width,
                            bool dump=false );

    void shortestpath( const std::vector< Eigen::Vector3f >& points,             // node list: position of (x,y,z)
                       const std::map< std::pair<int,int>, float >& distmap,     // distance between nodes
                       const std::map< std::pair<int,int>, float >& pixgapdist,
                       float& maxgapdist,
                       std::vector< std::vector<float> >& path_tyz ); // pixel gap between vertices
    

    // 3d points from spatial scans
    std::vector< std::vector<float> > m_points3d_v;
    std::vector< std::vector<float> > m_twid_v;

    // graph definitions
    std::vector<Eigen::Vector3f > m_graph_nodes;
    std::map< std::pair<int,int>, float > m_distmap;
    std::map< std::pair<int,int>, float > m_pixgapmap;

    float m_maxgapdist;
    std::vector< std::vector<float> > m_path_xyz;    // output path
    std::vector< std::vector<float> > m_path_twid_v; // output path in image coordinates   
    
  };
  
}
}

#endif
