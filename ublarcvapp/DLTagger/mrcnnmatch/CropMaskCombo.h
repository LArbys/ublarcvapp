#ifndef __CROP_MASK_COMBO_H__
#define __CROP_MASK_COMBO_H__

#include <vector>
#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "MRCNNMatchTypes.h"

namespace ublarcvapp {
namespace dltagger {

  /**
   * class that stores crops across the planes based on 
   *  mask matches represented in a MaskCombo class
   *
   * the charge in the crops are masked by the MRCNN masks.
   * the masks are stored.
   * the class also stores the features we collect on those crops.
   */
  class CropMaskCombo : larcv::larcv_base {

  public:

    CropMaskCombo()
      : larcv::larcv_base("CropMaskCombo"),
      _pcombo(nullptr)
    {};
    
    CropMaskCombo( const MaskCombo& combo,
                   const std::vector<larcv::Image2D>& adc,
                   const std::vector<larcv::Image2D>& wholeview_badch_v,
                   float pixel_threshold=10.0,
                   int tick_padding=6,
                   int downsample_factor=8 );
    
    virtual ~CropMaskCombo() {};
    
    std::vector<larcv::Image2D> crops_v;
    std::vector<larcv::Image2D> mask_v;
    std::vector<larcv::Image2D> missing_v; ///< [deprecated] adds non-empty image for plane that is missing in crops_v and mask_v
    std::vector<larcv::Image2D> badch_v;   ///< crop from bad channel image
    const MaskCombo& getCombo() const { return *_pcombo; };
    bool istwoplane() const { return _twoplane_mode; };
    int  badPlane() const   { return _badplane; };
    float getThreshold() const { return _threshold; };
    
  protected:
    
    const MaskCombo* _pcombo;
    void _crop_and_mask_image( const std::vector<larcv::Image2D>& wholeview_v );
    void _make_missing_crop( const std::vector<larcv::Image2D>& wholeview_v,
                             std::vector<larcv::Image2D>& out_v,
                             std::vector<larcv::Image2D>& miss_out_v );
    void _prep_badch_crop( const std::vector<larcv::Image2D>& whole_badch_v,
                           const std::vector<larcv::Image2D>& crop_v,
                           const std::vector<larcv::Image2D>& missing_v,
                           std::vector<larcv::Image2D>& badch_crop_v );
    

    // parameters
    float _threshold;         //< pixel value threshold
    int   _tick_padding;      //< padding added to tick dimensions
    int   _downsample_factor; //< crops will be a multple of this factor in order to allow for downsampling
    bool  _twoplane_mode;     //< if true, one plane crop is inferred from the two other crops
    int   _badplane;          //< index of inferred missing plane
    
  };
  
}
}

#endif
