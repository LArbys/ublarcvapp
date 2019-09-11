#ifndef __ContourClusterAlgo_H__
#define __ContourClusterAlgo_H__

/* -------------------------------------------------------------------------
 * ContourClusterAlgo
 * This class does 1-plane segment shape analysis using the contour tools
 * from opencv. Origionally developed by Vic Genty, Rui An, Kazu Terao
 * Ported to ublarcvapp context by Taritree
 * -----------------------------------------------------------------------*/

#include <vector>

#include "larcv/core/DataFormat/Image2D.h"

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>

#include "ContourShapeMeta.h"

namespace ublarcvapp {

  typedef std::vector<cv::Point> Contour_t;
  typedef std::vector< Contour_t > ContourList_t;
  typedef std::vector< int > ContourIndices_t;
  typedef std::vector< cv::Vec4i > Defects_t;

  class ContourClusterAlgo {
    
  public:
    ContourClusterAlgo();
    virtual ~ContourClusterAlgo() {};

    // main function
    void analyzeImages( const std::vector<larcv::Image2D>& img_v,
                        const float threshold=10.0, const int iterations=2,
                        const int min_defect_size=5,
                        const int hull_edge_pts_split=50,
                        const int n_allowed_breaks=10,
                        const int verbosity=2 );

    void analyzeImage( const larcv::Image2D& img,
                       ContourList_t& contour_v,
                       std::vector<ContourIndices_t>& hull_v,
                       std::vector<Defects_t>& defects_v,
                       ContourList_t& atomics_v,
                       std::vector< ContourShapeMeta >& atomicmeta_v,
                        const float threshold=10.0, const int iterations=2,
                        const int min_defect_size=5,
                        const int hull_edge_pts_split=50,
                        const int n_allowed_breaks=10,
                        const int verbosity=2 );

    
    void clear_intermediate_images();

  protected:
    
    void _splitContour( const std::vector<larcv::Image2D>& img_v,
                        const int min_defect_dize,
                        const int hull_edge_pts_split,
                        const int n_allowed_breaks );

    void _clear();

  public:
    
    std::vector<cv::Mat> cvimg_stage0_v; // unchanged images
    std::vector<cv::Mat> cvimg_stage1_v; // contour points over time scan
    std::vector<cv::Mat> cvimg_stage2_v; // 3D-matched contour points
    std::vector<cv::Mat> cvimg_stage3_v; // 3D-matched spacepointso

    std::vector< ContourList_t >                 m_plane_contours_v;
    std::vector< std::vector<ContourIndices_t> > m_plane_hulls_v;
    std::vector< std::vector<Defects_t> >        m_plane_defects_v;
    std::vector< ContourList_t >                 m_plane_atomics_v;
    std::vector< std::vector< ContourShapeMeta > > m_plane_atomicmeta_v;

  };


}

#endif
