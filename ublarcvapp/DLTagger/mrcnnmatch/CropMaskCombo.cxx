#include "CropMaskCombo.h"

// larlite
#include "LArUtil/Geometry.h"

namespace ublarcvapp {
namespace dltagger {
  
  CropMaskCombo::CropMaskCombo( const MaskCombo& xcombo, const std::vector<larcv::Image2D>& adc )
    : _pcombo(&xcombo)
  {
    _crop_and_mask_image( adc );
  }

  void CropMaskCombo::_crop_and_mask_image( const std::vector<larcv::Image2D>& wholeview_v ) {

    // we make crops based on the union of the mask-combos
    float tick_union_min = getCombo().union_tick[0]-6;
    float tick_union_max = getCombo().union_tick[1]+6;
    float detz_union_min = getCombo().union_detz[0];
    float detz_union_max = getCombo().union_detz[1];

    double testpts[4][3] = { { 0,  117.0, detz_union_min  },  // minz top
                             { 0, -117.0, detz_union_min  },  // minz bottom
                             { 0,  117.0, detz_union_max  },  // maxz top
                             { 0, -117.0, detz_union_max  }}; // maxz bot

    // adjust tick bounds to be divisible by pixel_height
    int tickdiff = int(tick_union_max-tick_union_min);
    int nrows    = tickdiff/int(wholeview_v.front().meta().pixel_height());
    if ( tickdiff%int(wholeview_v.front().meta().pixel_height())!=0 )
      nrows++;
    if ( nrows%8!=0 )
      nrows += 8-nrows%8;
    if ( nrows%8!=0 )
      throw std::runtime_error( "nrows not divisible by 8" );
    tick_union_max = tick_union_min + nrows*wholeview_v.front().meta().pixel_height();
    // adjust tick bounds to stay within min/max of full image
    if ( tick_union_max>=wholeview_v.front().meta().max_y() ) {
      tick_union_max = wholeview_v.front().meta().max_y();
      tick_union_min = tick_union_max - nrows*wholeview_v.front().meta().pixel_height();
    }
    if ( tick_union_min<wholeview_v.front().meta().min_y() ) {
      tick_union_min = wholeview_v.front().meta().min_y();
      tick_union_max = tick_union_min + nrows*wholeview_v.front().meta().pixel_height();
    }
    
    for ( size_t p=0; p<wholeview_v.size(); p++ ) {

      // need to get the min and max wire ID for detz position
      float near_wireco[4];
      int   near_wire[4];
      for ( int i=0; i<4; i++ ) {
        near_wireco[i] = larutil::Geometry::GetME()->WireCoordinate( testpts[i], p );
        // if ( near_wireco[i]<0 ) near_wireco[i] = 0;
        // if ( near_wireco[i]>=float(larutil::Geometry::GetME()->Nwires(p)) )
        //   near_wireco[i] = float(larutil::Geometry::GetME()->Nwires(p)-1);
          
        near_wire[i] = (int)(near_wireco[i]);
        if ( near_wire[i]<0 ) near_wire[i] = 0;
        if ( near_wire[i]>=(int)larutil::Geometry::GetME()->Nwires(p) )
          near_wire[i] = (int)larutil::Geometry::GetME()->Nwires(p) - 1;
      }
      
      int wire_min = near_wire[0];
      int wire_max = near_wire[0];
      for ( int i=1; i<4; i++ ) {
        if ( wire_min>near_wire[i] ) wire_min = near_wire[i];
        if ( wire_max<near_wire[i] ) wire_max = near_wire[i];
      }
      int ncols = wire_max-wire_min+1;
      if (ncols%8!=0)
        ncols += 8-ncols%8;
      if ( ncols%8!=0 )
        throw std::runtime_error( "nrows not divisible by 8" );
      
      wire_max = wire_min + ncols*wholeview_v[p].meta().pixel_width();
      if ( wire_max>larutil::Geometry::GetME()->Nwires(p) ) {
        wire_max = larutil::Geometry::GetME()->Nwires(p);
        wire_min = wire_max-ncols*wholeview_v[p].meta().pixel_width();
      }
      
      const larcv::ClusterMask& mask = *getCombo().pmasks.at(p);
      const MaskMatchData& data      = *getCombo().pdata.at(p);
      const larcv::ImageMeta& meta   = wholeview_v.at(p).meta();
      float width  = wire_max-wire_min;
      float height = tick_union_max-tick_union_min;
      
      // make a meta to crop
      larcv::ImageMeta cropmeta( width, height,
                                 nrows,
                                 ncols,
                                 wire_min,
                                 tick_union_min, meta.plane() );
      larcv::Image2D crop = wholeview_v.at(p).crop( cropmeta );
      larcv::Image2D imgmask( cropmeta );
      imgmask.paint(0.0);

      //std::cout << "meta: " << cropmeta.dump() << std::endl;

      // make mask
      for ( size_t ipt=0; ipt<mask.points_v.size(); ipt++ ) {
        int row = mask.points_v[ipt].y;
        int col = mask.points_v[ipt].x;
        if ( row>=mask.meta.rows() ) row--;
        // translate from ClusterMask coordinate system to crop coordinates

        try {
          float ptwire = mask.meta.pos_x( mask.box.min_x() + col );
          float pttick = mask.meta.pos_y( mask.box.min_y() + row );
          
          int xrow = cropmeta.row( pttick );
          int xcol = cropmeta.col( ptwire );
        
          imgmask.set_pixel( xrow, xcol, 50.0 );
        }
        catch (std::exception& e) {
          std::cout << "warning mask out of bounds: " << e.what() << std::endl;
        };
      }

      std::vector<float>& vec_mask = imgmask.as_mod_vector();
      std::vector<float>& vec_pix  = crop.as_mod_vector();
      // apply mask to charge image
      // for ( size_t idx=0; idx<vec_pix.size(); idx++ ) {
      //   vec_pix[idx] *= vec_mask[idx];
      // }
      // apply threshold to both charge and mask image
      for ( size_t idx=0; idx<vec_pix.size(); idx++ ) {
        if ( vec_pix[idx]<10.0 ) vec_pix[idx]  = 0.0;
        if ( vec_pix[idx]<10.0 ) vec_mask[idx] = 0.0;
      }
      
      crops_v.emplace_back( std::move(crop) );
      mask_v.emplace_back( std::move(imgmask) );
    }
  }
  
}
}
