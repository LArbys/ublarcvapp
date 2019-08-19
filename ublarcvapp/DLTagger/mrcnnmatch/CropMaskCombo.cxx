#include "CropMaskCombo.h"

namespace ublarcvapp {
namespace dltagger {
  
  CropMaskCombo::CropMaskCombo( const MaskCombo& xcombo, const std::vector<larcv::Image2D>& adc )
    : _pcombo(&xcombo)
  {
    _crop_and_mask_image( adc );
  }

  void CropMaskCombo::_crop_and_mask_image( const std::vector<larcv::Image2D>& wholeview_v ) {
    
    for ( size_t p=0; p<wholeview_v.size(); p++ ) {
      
      const larcv::ClusterMask& mask = *getCombo().pmasks.at(p);
      const MaskMatchData& data      = *getCombo().pdata.at(p);
      const larcv::ImageMeta& meta   = wholeview_v.at(p).meta();
      
      // make a meta to crop
      float wire_min = meta.pos_x( (int)mask.box.min_x() );
      float tick_min = meta.pos_y( (int)mask.box.min_y() );
      float wire_max = meta.pos_x( (int)mask.box.max_x() );
      float tick_max = meta.pos_y( (int)mask.box.max_y() );
      float width    = fabs( wire_max-wire_min ) + meta.pixel_width();
      float height   = fabs( tick_max-tick_min ) + meta.pixel_height();

      larcv::ImageMeta cropmeta( width, height,
                                 int(mask.box.max_y()-mask.box.min_y())+1,
                                 int(mask.box.max_x()-mask.box.min_x())+1,
                                 wire_min, tick_min, meta.plane() );
      larcv::Image2D crop = wholeview_v.at(p).crop( cropmeta );
      larcv::Image2D imgmask( cropmeta );
      imgmask.paint(0.0);

      //std::cout << "meta: " << cropmeta.dump() << std::endl;

      // make mask
      for ( size_t ipt=0; ipt<mask.points_v.size(); ipt++ ) {
        int row = mask.points_v[ipt].y;
        int col = mask.points_v[ipt].x;
        //std::cout << " mask pt(" << row << "," << col << ")" << std::endl;        
        imgmask.set_pixel( row, col, 1.0 );
      }

      const std::vector<float>& vec_mask = imgmask.as_vector();
      std::vector<float>& vec_pix = crop.as_mod_vector();
      for ( size_t idx=0; idx<vec_pix.size(); idx++ ) {
        vec_pix[idx] *= vec_mask[idx];
      }
      for ( size_t idx=0; idx<vec_pix.size(); idx++ )
        if ( vec_pix[idx]<10.0 ) vec_pix[idx] = 0.0;
      
      crops_v.emplace_back( std::move(crop) );
      mask_v.emplace_back( std::move(imgmask) );
    }
  }
  
}
}
