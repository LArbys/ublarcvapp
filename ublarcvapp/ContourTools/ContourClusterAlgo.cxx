#include "ContourClusterAlgo.h"


// LArCV OpenCV
#include "larcv/core/CVUtil/CVUtil.h"

// LArOpenCV
#include "LArOpenCV/ImageCluster/AlgoClass/DefectBreaker.h"
#include "LArOpenCV/ImageCluster/AlgoData/TrackClusterCompound.h"

#include "TRandom3.h"

namespace ublarcvapp {

  ContourClusterAlgo::ContourClusterAlgo() 
  {
    _clear();
  }


  void ContourClusterAlgo::_clear() {

    cvimg_stage0_v.clear(); // unchanged images
    cvimg_stage1_v.clear(); // contour points over time scan
    cvimg_stage2_v.clear(); // 3D-matched contour points
    cvimg_stage3_v.clear(); // 3D-matched spacepointso

    m_plane_contours_v.clear();
    m_plane_hulls_v.clear();
    m_plane_defects_v.clear();
    m_plane_atomics_v.clear();
    m_plane_atomicmeta_v.clear();
  }

  /**
   * provide images and produce clusters using contour defect analysis.
   * aim is to produce 'straight' contours.
   * these are mostly the result of the analsys, but technically the algo
   * will produce any convex contour
   *
   * @param[in] img_v vector of images
   * @param[in] threshold threshold to apply to images
   * @param[in] dilate_iterations dilation iterations before contouring
   */
  void ContourClusterAlgo::analyzeImages( const std::vector<larcv::Image2D>& img_v,
                                          const float threshold,
                                          const int dilate_iterations,
                                          const int min_defect_size,
                                          const int hull_edge_pts_split,
                                          const int n_allowed_breaks,
                                          const int verbosity )
  {
    
    TRandom3 rand(1983);
    _clear();
    
    // first convert the images into cv and binarize
    for ( auto const& img : img_v ) {

      // check if image is empty
      if ( img.meta().cols()==0 || img.meta().rows()==0 ) {
        // fill blank
        cv::Mat cvimg;
        cv::Mat thresh;
        cvimg_stage0_v.emplace_back( std::move(cvimg) );
        cvimg_stage1_v.emplace_back( std::move(thresh) );
        continue;
      }
      
      cv::Mat cvimg = larcv::as_gray_mat( img, threshold, 256.0, 1.0 );
      cv::Mat cvrgb = larcv::as_mat_greyscale2bgr( img, threshold, 100.0 );
      cv::Mat thresh( cvimg );
      cv::threshold( cvimg, thresh, 0, 255, cv::THRESH_BINARY );
      cvimg_stage0_v.emplace_back( std::move(cvrgb) );
      cvimg_stage1_v.emplace_back( std::move(thresh) );
    }


    for (int p=0; p<3; p++) {

      if ( img_v[p].meta().cols()==0 || img_v[p].meta().rows()==0 ) {
        // fill blank
        cv::Mat stage2;
        cv::Mat stage3;
        cvimg_stage2_v.emplace_back( std::move(stage2) );
        cvimg_stage3_v.emplace_back( std::move(stage3) );

        m_plane_contours_v.push_back( ContourList_t() );
        m_plane_hulls_v.push_back( std::vector<ContourIndices_t>() );
        m_plane_defects_v.push_back( std::vector<Defects_t>() );
        m_plane_atomics_v.push_back( ContourList_t() );
        m_plane_atomicmeta_v.push_back( std::vector< ContourShapeMeta >() );
        continue;
      }
      
      // dilate image first
      cv::Mat& cvimg = cvimg_stage1_v[p];
      cv::dilate( cvimg, cvimg, cv::Mat(), cv::Point(-1,-1), dilate_iterations, 1, 1 );
      
      // find contours
      ContourList_t contour_v;
      cv::findContours( cvimg_stage1_v[p], contour_v, cv::RETR_LIST, cv::CHAIN_APPROX_SIMPLE );

      //std::cout << "Plane " << p << " number of contours: " << contour_v.size() << std::endl;
      
      // for each contour, find convex hull, find defect points
      std::vector< ContourIndices_t > hull_v( contour_v.size() );
      std::vector< Defects_t > defects_v( contour_v.size() );
      for ( int idx=0; idx<(int)contour_v.size(); idx++ ) {

	Contour_t& contour = contour_v[idx];
	if ( contour.size()<10 )
	  continue;
	
	// draw contours
	//cv::drawContours( cvimg_stage0_v[p], contour_v, idx, cv::Scalar( rand.Uniform(10,255),rand.Uniform(10,255),rand.Uniform(10,255),255), 1 );	
	
	// convex hull
	cv::convexHull( cv::Mat( contour ), hull_v[idx], false );

	if ( hull_v[idx].size()<=3 ) {
	  // store 
	  continue; // no defects can be found
	}

	// defects
	cv::convexityDefects( contour, hull_v[idx], defects_v[idx] );

	// plot defect point information
	for ( auto& defectpt : defects_v[idx] ) {
	  float depth = defectpt[3]/256;
	  if ( depth>3 ) {
	    int faridx = defectpt[2];
	    cv::Point ptFar( contour[faridx] );
	    cv::circle( cvimg_stage0_v[p], ptFar, 1, cv::Scalar(0,255,0,255), -1 );
	  }
	}

      }
      m_plane_contours_v.emplace_back( std::move(contour_v) );
      m_plane_hulls_v.emplace_back( std::move(hull_v) );
      m_plane_defects_v.emplace_back( std::move(defects_v) );
    }

    // split contours based on their defects, leaving only convex contours
    _splitContour( img_v, min_defect_size, hull_edge_pts_split, n_allowed_breaks );

    // mark image so we can look up contour via position in image
    //_makeContourIndexImage();

  }// end of findBoundarySpacePoints


  /**
   * given a contour with at least one defect point, 
   * we split it with the aim of producing the straightest daughter contours
   * setup the contour breaker
   *
   * @param[in] img_v original larcv images
   * @param[in] min_defect_size
   * @param[in] hull_edge_pts_split
   * @param[in] n_allowed_breaks
   *
   */
  void ContourClusterAlgo::_splitContour( const std::vector<larcv::Image2D>& img_v,
                                          const int min_defect_size,
                                          const int hull_edge_pts_split,
                                          const int n_allowed_breaks )
  {
    
    larocv::DefectBreaker defectb;
    defectb.Configure( min_defect_size,
                       hull_edge_pts_split,
                       n_allowed_breaks );

    TRandom3 rand(time(NULL));    

    m_plane_atomics_v.clear();
    m_plane_atomics_v.resize(3);
    m_plane_atomicmeta_v.clear();
    m_plane_atomicmeta_v.resize(3);
    
    for (int p=0; p<3; p++) {
      // atomic_contours;
      auto& contour_v = m_plane_contours_v[p];
      for ( auto& contour : contour_v ) {
	larocv::data::TrackClusterCompound atomics = defectb.BreakContour( contour ); // returns a vector<AtomicContour>
	for ( auto& ctr : atomics ) {
	  if ( ctr.size()<4 )
	    continue;

	  ContourShapeMeta ctrinfo( ctr, img_v[p] );	  
	  m_plane_atomics_v[p].push_back( std::move(ctr) );
	  m_plane_atomicmeta_v[p].emplace_back( std::move(ctrinfo) );
          
	}
      }

      // for visualization
      for (int idx=0; idx<(int)m_plane_atomics_v[p].size(); idx++) {
	cv::drawContours( cvimg_stage0_v[p], m_plane_atomics_v[p],
                          idx, cv::Scalar( rand.Uniform(10,255),rand.Uniform(10,255),
                                           rand.Uniform(10,255),255), 1 );
	cv::line( cvimg_stage0_v[p], m_plane_atomicmeta_v[p].at(idx).getFitSegmentStart(),
                  m_plane_atomicmeta_v[p].at(idx).getFitSegmentEnd(), cv::Scalar(255,255,255,255), 1 );
      }
    }    

    
  }//splitContour

  /**
   * clear out intermediate products to save memory
   */
  void ContourClusterAlgo::clear_intermediate_images() {
    cvimg_stage0_v.clear();
    cvimg_stage1_v.clear();
    cvimg_stage2_v.clear();
    cvimg_stage3_v.clear();
  }

  /**
   * analyze images
   *
   */
  void ContourClusterAlgo::analyzeImage( const larcv::Image2D& img,
                                         ContourList_t& contour_v,
                                         std::vector<ContourIndices_t>& hull_v,
                                         std::vector<Defects_t>& defects_v,
                                         ContourList_t& atomics_v,
                                         std::vector< ContourShapeMeta >& atomicmeta_v,                                         
                                         const float threshold,
                                         const int dilate_iterations,
                                         const int min_defect_size,
                                         const int hull_edge_pts_split,
                                         const int n_allowed_breaks,
                                         const int verbosity ) {

    if ( img.meta().cols()==0 || img.meta().rows()==0 )
      return;
    
    cv::Mat cvimg = larcv::as_gray_mat( img, threshold, 256.0, 1.0 );
    cv::Mat cvrgb = larcv::as_mat_greyscale2bgr( img, threshold, 100.0 );
    cv::Mat thresh( cvimg );
    cv::threshold( cvimg, thresh, 0, 255, cv::THRESH_BINARY );

    // dilate image first
    cv::dilate( thresh, thresh, cv::Mat(), cv::Point(-1,-1), dilate_iterations, 1, 1 );
      
    // find contours
    //ContourList_t contour_v;
    contour_v.clear();
    cv::findContours( thresh, contour_v, cv::RETR_LIST, cv::CHAIN_APPROX_SIMPLE );

    // for each contour, find convex hull, find defect points
    hull_v.clear();
    defects_v.clear();    
    hull_v.resize( contour_v.size() );
    defects_v.resize( contour_v.size() );
    
    for ( int idx=0; idx<(int)contour_v.size(); idx++ ) {

      Contour_t& contour = contour_v[idx];
      if ( contour.size()<10 )
        continue;
	
      // draw contours
      //cv::drawContours( cvimg_stage0_v[p], contour_v, idx, cv::Scalar( rand.Uniform(10,255),rand.Uniform(10,255),rand.Uniform(10,255),255), 1 );	
	
      // convex hull
      cv::convexHull( cv::Mat( contour ), hull_v[idx], false );
      
      if ( hull_v[idx].size()<=3 ) {
        // store 
        continue; // no defects can be found
      }
      
      // defects
      cv::convexityDefects( contour, hull_v[idx], defects_v[idx] );
    }

    // split contours
    larocv::DefectBreaker defectb;
    defectb.Configure( min_defect_size,
                       hull_edge_pts_split,
                       n_allowed_breaks );


    atomics_v.clear();
    atomicmeta_v.clear();
    for ( auto& contour : contour_v ) {
      larocv::data::TrackClusterCompound atomics = defectb.BreakContour( contour ); // returns a vector<AtomicContour>
      for ( auto& ctr : atomics ) {
        if ( ctr.size()<4 )
          continue;
        
        ContourShapeMeta ctrinfo( ctr, img );
        atomics_v.push_back( std::move(ctr) );
	atomicmeta_v.emplace_back( std::move(ctrinfo) );
      }
    }    
    
  }

  /**
   * make images with index stored in cluster
   *
   * @param[in] img_v original larcv images
   * @param[in] threshold pixel value threshold
   * 
   * output stored in m_plane_indeximg_v
   *
   */
  /*
  void ContourClusterAlgo::_makeContourIndexImage( const std::vector<larcv::Image2D>& img_v,
                                                   const float threshold )
  {
    
    m_plane_indeximg_v.clear();

    for ( auto const& img : img_v ) {
      larcv::Image2D out(img.meta());
      out.paint(0);

      for ( size_t r=0; r<img.meta().rows(); r++ ) {
        for ( size_t c=0; c<img.meta().cols(); c++ ) {
          if ( img.pixel(r,c)>threshold ) {
            cv::Point pt( c,r );
            for ( auto& ctr : m_plane_atomics_v.at((size_t)img.plane()) ) {
              double result =  cv::pointPolygonTest( ctr, pt, false );
              
            }
            
          }
        }
      }
    
  }//makeContourIndexImage
  */
  
  /*
  bool BMTCV::isPointInContour( const std::vector<float>& pos3d, const float maxdist, const int minsize ) {
    // provides a simple test
    // project into planes
    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    // HAS OPENCV

    //std::vector<int> imgcoords = larcv::UBWireTool::
    
    // ------------------------------------------------------------------------    

    return true;
  }
  */
}
