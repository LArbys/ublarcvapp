/**
 * \file UBSplitDetector.h
 *
 * \ingroup UBImageMod
 *
 * \brief Class def header for a class UBSplitDetector
 *
 * @author twongjirad
 *
 * We carve the detector into 3D regions and define vectors of bounding boxes to
 * then produce cropped images.
 */

/** \addtogroup UBImageMod

    @{*/
#ifndef __UBSPLITDETECTOR_H__
#define __UBSPLITDETECTOR_H__

#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

#include "larcv/core/DataFormat/EventImage2D.h"

namespace ublarcvapp {

  /**
     \class ProcessBase
     User defined class UBSplitDetector ... these comments are used to generate
     doxygen documentation!
  */
  class UBSplitDetector : public larcv::ProcessBase {

  public:

    /// Default constructor
    UBSplitDetector(const std::string name = "UBSplitDetector");

    /// Default destructor
    ~UBSplitDetector() {}

    void configure(const larcv::PSet&);

    void initialize();

    bool process( larcv::IOManager& mgr);

    void finalize();

    // algo functions

    bool process( const std::vector<larcv::Image2D>& img_v,
		  std::vector<larcv::Image2D>& outimg_v, 
		  std::vector<larcv::ROI>& outbbox_v);

    // static functions are defined to allow them to be reused in a stand-alone manner
    static larcv::ROI defineBoundingBoxFromCropCoords( const std::vector<larcv::Image2D>& img_v,
						       const int box_pixel_width, const int box_pixel_height,
						       const int t1, const int t2,
						       const int u1, const int u2,
						       const int v1, const int v2,
						       const int y1, const int y2,
						       const bool tick_forward);

    static bool cropUsingBBox2D( const larcv::ROI& bbox_vec,
				 const std::vector<larcv::Image2D>& img_v,
				 const int y1, const int y2, bool fill_y_image,
				 const float minpixfrac,
				 const int first_outidx,
				 const bool copy_imgs,
				 larcv::EventImage2D& output_imgs );

    static bool cropUsingBBox2D( const larcv::ROI& bbox_vec,
				 const std::vector<larcv::Image2D>& img_v,
				 const int y1, const int y2, bool fill_y_image,
				 const float minpixfrac,
				 const int first_outidx,
				 const bool copy_imgs,
				 std::vector<larcv::Image2D>& outimg_v );


    static std::vector<int> defineImageBoundsFromPosZT( const float zwire, const float tmid,
							const float zwidth, const float dtick,
							const int box_pixel_width, const int box_pixel_height,
							const std::vector<larcv::Image2D>& img_v,
							const bool tick_forward );


    void printElapsedTime();
    void clearElapsedTime();
    bool willRandomize() { return _randomize_crops; };

  private:

    // config parameters
    std::string _input_producer;
    std::string _output_bbox_producer;
    std::string _output_img_producer;
    bool _enable_img_crop;
    int _box_pixel_height;
    int _box_pixel_width;
    int _covered_z_width;
    bool _complete_y_crop;
    bool _debug_img;
    int _max_images;
    // randomize cropping parameters
    bool _randomize_crops;
    int  _randomize_attempts;
    float _randomize_minfracpix;
    int  _num_expected_crops;
    bool _numcrops_changed;
    bool _tick_forward;

    // cropping variables
    std::vector< std::vector<int> > m_lattice;    //< defines (t,u,v,y) wire coordinates that serve as center of cropped subimages
    std::vector<larcv::Image2D>     m_coverage_v; //< images used to mark how many times pixels are a part of subimage

    // timetracking
    static float elapsed_genbbox;
    static float elapsed_crop;
    static float elapsed_alloc;
    static float elapsed_fraccheck;
    static float elapsed_save;
    static int   num_calls;
  };

  /**
     \class larcv::UBSplitDetectorFactory
     \brief A concrete factory class for ublarcvapp::UBSplitDetector
  */
  class UBSplitDetectorProcessFactory : public larcv::ProcessFactoryBase {
  public:
    /// ctor
    UBSplitDetectorProcessFactory() { larcv::ProcessFactory::get().add_factory("UBSplitDetector", this); }
    /// dtor
    ~UBSplitDetectorProcessFactory() {}
    /// creation method
    larcv::ProcessBase* create(const std::string instance_name) { return new UBSplitDetector(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group
