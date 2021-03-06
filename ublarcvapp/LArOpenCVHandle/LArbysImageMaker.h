#ifndef LARBYSIMAGEMAKER_H
#define LARBYSIMAGEMAKER_H

#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/Base/PSet.h"
#include "LArOpenCV/Core/ImageManager.h"
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventPixel2D.h"
#include <tuple>

#ifndef __CLING__
#ifndef __CINT__
#include <opencv2/core.hpp>
#endif
#endif

namespace larcv {
  class LArbysImageMaker : public larcv_base {

  public:
  LArbysImageMaker() :
    _charge_max(255),
      _charge_min(0),
      _charge_to_gray_scale(1)
	{}
    
    ~LArbysImageMaker(){}

    void
      Configure(const PSet& pset);

    std::vector<cv::Mat>
      ExtractMat(const std::vector<Image2D>& image_v);

    cv::Mat
      ExtractMat(const Image2D& image);
    
    std::tuple<cv::Mat,larocv::ImageMeta>
      ExtractImage(const Image2D& image, size_t plane=0);
    
    std::vector<std::tuple<cv::Mat,larocv::ImageMeta> >
      ExtractImage(const std::vector<Image2D>& image_v);
    

    Image2D ConstructCosmicImage(const EventPixel2D& ev_pixel2d,
				 const Image2D& adc_image,
				 const size_t plane,
				 float value=100);

    Image2D ConstructCosmicImage(const EventPixel2D* ev_pixel2d,
				 const Image2D& adc_image,
				 const size_t plane,
				 float value);
    
    void ConstructCosmicImage(IOManager& mgr,
			      std::string producer,
			      ProductType_t datatype, 
			      const std::vector<larcv::Image2D>& adc_image_v,
			      std::vector<larcv::Image2D>& mu_image_v);
    
  public:

    float _charge_max;
    float _charge_min;
    float _charge_to_gray_scale;
    
  };
}
#endif
/** @} */ // end of doxygen group 

