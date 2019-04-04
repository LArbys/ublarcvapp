/**
 * \file InfillDataCropper.h
 *
 * \ingroup UBImageMod
 *
 * \brief Class def header for a class InfillDataCropper
 *
 * @author kmason
 *
 * We crop out data to make necessary infill images.
 *
 */

/** \addtogroup UBImageMod

    @{*/
#ifndef __INFILLDATACROPPER_H__
#define __INFILLDATACROPPER_H__

#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"
namespace ublarcvapp {

  /**
     \class ProcessBase
     User defined class InfillDataCropper ... these comments are used to generate
     doxygen documentation!
  */
  class InfillDataCropper : public larcv::ProcessBase {

  public:

    /// Default constructor
    InfillDataCropper(const std::string name = "InfillDataCropper");

    /// Default destructor
    ~InfillDataCropper() {}

    void configure(const larcv::PSet&);

    void initialize();

    bool process(larcv::IOManager& mgr);

    void finalize();

    // ----------------------------------------------------------------------
    // functions
    static std::vector<larcv::Image2D> ChStatusToLabels(
                          std::vector<larcv::Image2D>& image_label_v,
                          larcv::EventChStatus* ev_chstatus);

    static std::map<int, std::vector<int>> DeadChannelAnalyzer(
                          const std::vector<larcv::Image2D>& image_label_v,
                          const int cropwidth);

    static std::vector<larcv::Image2D> cropLabelsImage(
                          const int cropwidth,
                          const int cropheight,
                          const std::vector<larcv::Image2D>& img_v,
                          std::vector<larcv::Image2D>& cropped_labels);

    static std::vector<larcv::ImageMeta> defineBoundingBoxFromCropCoords(
                           const std::vector<larcv::Image2D>& img_v,
  									       const int box_pixel_width,
                           const int box_pixel_height,
                           const double t2, const double u1,
                           const double v1, const double y1);

    static  std::vector<larcv::Image2D> cropImagefromBox(
                           std::vector<larcv::ImageMeta> bbox_vec,
                					 const std::vector<larcv::Image2D>& img_v,
                					 std::vector<larcv::Image2D>& output_imgs);

    static bool checkLabelsCrop(std::vector<larcv::Image2D>& cropped_labels);

    static std::vector<larcv::Image2D> cropLabelsImage(
                          const int cropwidth,
                          const int cropheight,
                          const std::vector<larcv::Image2D>& img_v,
                          std::vector<larcv::Image2D>& cropped_adc,
                          std::map<int, std::vector<int>> startdict);

    static std::vector<larcv::Image2D> cropADCImage(
                          const int cropwidth,
                          const int cropheight,
                          const std::vector<larcv::Image2D>& img_v,
                          std::vector<larcv::Image2D>& cropped_adc,
                          std::map<int, std::vector<int>> startdict);

    static std::vector<larcv::Image2D> cropADCMaskImage(
                          const int cropwidth,
                          const int cropheight,
                          const std::vector<larcv::Image2D>& img_v,
                          const std::vector<larcv::Image2D>& img_labels_v,
                          const std::vector<larcv::Image2D>& img_adc_v,
                          std::vector<larcv::Image2D>& cropped_adcmasked);

    static  bool checkADCCrop(
                              const int cropwidth,
                              const int cropheight,
                              std::vector<larcv::Image2D>& cropped_adc,
                              std::vector<larcv::Image2D>& cropped_labels);

    static std::vector<larcv::Image2D> cropWeightsImage(
                                            const int cropwidth,
                                            const int cropheight,
                                            std::vector<larcv::Image2D>&weights_v,
                                            const std::vector<larcv::Image2D>& cropped_labels,
                                            const std::vector<larcv::Image2D>& cropped_adc);

    // ----------------------------------------------------------------------
    // we save data ourselves

  public:

    const larcv::IOManager& getOutputIOMan() const { return *foutIO; };

  protected:

    larcv::IOManager* foutIO;

    // ----------------------------------------------------------------------

  private:

    // config parameters
    std::string _input_chstatus_producer;
    std::string _input_adc_producer;

    //std::string _output_labels_producer;
    std::string _output_adc_producer;
    std::string _output_labels_producer;
    std::string _output_adcmasked_producer;
    std::string _output_weights_producer;
    std::string _output_meta_producer;
    std::string _output_filename;

    int _max_images;
    std::vector<float> _thresholds_v;
    bool _limit_overlap;
    float _max_overlap_fraction;
    int _verbosity_;

  };

  /**
     \class larcv::InfillDataCropperFactory
     \brief A concrete factory class for larcv::InfillDataCropper
  */
  class InfillDataCropperProcessFactory : public larcv::ProcessFactoryBase {
  public:
    /// ctor
    InfillDataCropperProcessFactory() { larcv::ProcessFactory::get().add_factory("InfillDataCropper", this); }
    /// dtor
    ~InfillDataCropperProcessFactory() {}
    /// creation method
    larcv::ProcessBase* create(const std::string instance_name) { return new InfillDataCropper(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group
