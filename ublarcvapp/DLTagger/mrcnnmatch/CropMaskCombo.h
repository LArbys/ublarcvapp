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
    
    CropMaskCombo( const MaskCombo& combo, const std::vector<larcv::Image2D>& adc,
                   float pixel_threshold=10.0,
                   int tick_padding=6, int downsample_factor=8 );
    
    virtual ~CropMaskCombo() {};
    
    std::vector<larcv::Image2D> crops_v;
    std::vector<larcv::Image2D> mask_v;
    std::vector<larcv::Image2D> missing_v; ///< adds non-empty image for plane that is missing in crops_v and mask_v
    int ngoodplanes;
    int badplane;
    const MaskCombo& getCombo() const { return *_pcombo; };

    
  protected:
    
    const MaskCombo* _pcombo;
    void _crop_and_mask_image( const std::vector<larcv::Image2D>& wholeview_v );
    void _make_missing_crop( const std::vector<larcv::Image2D>& wholeview_v,
                             std::vector<larcv::Image2D>& out_v );

    // parameters
    float _threshold;         //< pixel value threshold
    int   _tick_padding;      //< padding added to tick dimensions
    int   _downsample_factor; //< crops will be a multple of this factor in order to allow for downsampling

    
  };
  
}
}

#endif
