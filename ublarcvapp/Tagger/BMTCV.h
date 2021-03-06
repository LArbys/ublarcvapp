#ifndef __BMTCV_H__
#define __BMTCV_H__

/* -------------------------------------------------------------------------
 * BMTCV: BoundaryMuonTagger CV
 * This class does 1-plane segment shape analysis using the contour tools
 * from opencv.
 * -----------------------------------------------------------------------*/

#include <vector>

#include "larcv/core/DataFormat/Image2D.h"

#include "TaggerCROITypes.h"

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>

// ublarcvapp
#include "ublarcvapp/ContourTools/ContourShapeMeta.h"

namespace ublarcvapp {
namespace tagger {

  typedef std::vector<cv::Point> Contour_t;
  typedef std::vector< Contour_t > ContourList_t;
  typedef std::vector< int > ContourIndices_t;
  typedef std::vector< cv::Vec4i > Defects_t;

  class BMTCV {
  public:
    BMTCV(){};
    virtual ~BMTCV() {};


    std::vector<BoundarySpacePoint> findBoundarySpacePoints( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v );

    void analyzeImages( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const float threshold, const int iterations );

    void splitContour( const std::vector<larcv::Image2D>& img_v );    

    void clear();

    std::vector<cv::Mat> cvimg_stage0_v; // unchanged images
    std::vector<cv::Mat> cvimg_stage1_v; // contour points over time scan
    std::vector<cv::Mat> cvimg_stage2_v; // 3D-matched contour points
    std::vector<cv::Mat> cvimg_stage3_v; // 3D-matched spacepointso

    std::vector< ContourList_t >                 m_plane_contours_v;
    std::vector< std::vector<ContourIndices_t> > m_plane_hulls_v;
    std::vector< std::vector<Defects_t> >        m_plane_defects_v;
    std::vector< ContourList_t >                 m_plane_atomics_v;
    std::vector< std::vector< ublarcvapp::ContourShapeMeta > > m_plane_atomicmeta_v;

  };


}
}

#endif
