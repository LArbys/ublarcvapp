#ifndef __ContourShapeMeta__
#define __ContourShapeMeta__

#include <vector>

#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/DataFormat/Image2D.h"

#include <opencv2/core.hpp>


namespace ublarcvapp {

 class ContourShapeMeta : public std::vector<cv::Point> {
   // Wrapper around OpenCV contour
   // Stores meta data for contours
   
 public:
   ContourShapeMeta();   
   ContourShapeMeta( const std::vector<cv::Point>& contour, const larcv::Image2D& img, float threshold=10.0 );   
   virtual ~ContourShapeMeta() {};

   const larcv::ImageMeta& meta() const { return m_meta; };    
   const cv::Point& getFitSegmentStart() const { return m_start; };
   const cv::Point& getFitSegmentEnd() const { return m_end; };   
   const cv::Rect&  getBBox() const  { return m_bbox; };

   std::vector<float> getEndDir() const { return m_dir; };
   std::vector<float> getStartDir() const {
     std::vector<float> reverse_dir(m_dir.size(),0);
     for (size_t i=0; i<m_dir.size(); i++) reverse_dir[i] = -1.0*m_dir[i];
     return reverse_dir;
   };
   
   float getMinX() const { return xbounds[0]; };
   float getMaxX() const { return xbounds[1]; };
   float getMinY() const { return ybounds[0]; };
   float getMaxY() const { return ybounds[1]; };

   bool hasValidPCA() const { return m_valid_pca; };
   std::vector<float> getPCAdir( int axis=0 ) const;
   double getPCAeigenvalue( int axis ) const { return eigen_val[axis]; };
   std::vector<float> getPCAStartdir() const;
   std::vector<float> getPCAEnddir() const;
   std::vector<float> getPCAStartPos() const { return m_pca_startpt; };
   std::vector<float> getPCAEndPos() const   { return m_pca_endpt; };

   const std::vector<cv::Point>& getChargePixels() const { return qpixels; };
   
 protected:
   
   // ImageMeta
   const larcv::ImageMeta m_meta;
   
   // Line Fit/Projected End
   std::vector<float> m_dir;
   cv::Point m_start;
   cv::Point m_end;
   void _fill_linefit_members();

   // Bounding Box (for collision detection)
   cv::Rect m_bbox;
   void _build_bbox();

   // Bounds
   std::vector<float> ybounds;
   std::vector<float> xbounds;
   void _get_tick_range();

   // Charge core PCA
   void _charge_core_pca( const larcv::Image2D& img, float );
   bool m_valid_pca;
   cv::Point center;
   std::vector<cv::Point2d> eigen_vecs;
   std::vector<double> eigen_val;
   // start and end points determined by radius from center and either neg or pos on major axis
   std::vector<float> m_pca_startpt; 
   std::vector<float> m_pca_endpt;
   // points within the contour with charge
   float m_threshold;
   std::vector< cv::Point > qpixels;   
   
 };
 

}

#endif
