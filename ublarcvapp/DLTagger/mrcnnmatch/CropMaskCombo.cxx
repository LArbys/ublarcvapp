#include "CropMaskCombo.h"

// larlite
#include "LArUtil/Geometry.h"

namespace ublarcvapp {
namespace dltagger {

  /** 
   * class containing image crops made from Mask-RCNN bounding boxes matched across planes
   *
   * populate the following data members:
   *  crops_v: adc image crops. tick dimension cropped is the union of the time spanned by the masks across planes.
   *  badch_v: cropped image with bad channels marked
   *  missing_v: [deprecated] contains crop whose bounding box is inferred from the other two planes
   *  mask_v:  adc image crops. However, only ADC values inside the Mask-RCNN mask is saved.
   *
   * @param[in] xcombo Produced by MRCNNMatch class. Represents Mask-RCNN match across 3 or 2 planes. 
   *                   Contains pointer to Mask-RCNN masks.
   * @param[in] adc    Whole image ADC from which we will crop.
   * @param[in] wholeview_badch_v Whole image labeling pixels that are a part of bad channels. We will crop from this.
   * @param[in] pixel_threshold Threshold ADC image with this value. [Default 10]
   * @param[in] tick_padding  [deprecated]
   * @param[in] downsample_factor [deprecated]
   */
  CropMaskCombo::CropMaskCombo( const MaskCombo& xcombo,
                                const std::vector<larcv::Image2D>& adc,
                                const std::vector<larcv::Image2D>& wholeview_badch_v,
                                float pixel_threshold,
                                int tick_padding,
                                int downsample_factor )
    : larcv::larcv_base("CropMaskCombo"),
      _pcombo(&xcombo),
      _threshold(pixel_threshold),
      _tick_padding(tick_padding),
      _downsample_factor(downsample_factor),
      _twoplane_mode(false),
      _badplane(-1)
  {
    _crop_and_mask_image( adc );
    _make_missing_crop( adc, crops_v, missing_v );
    _prep_badch_crop( wholeview_badch_v, crops_v, missing_v, badch_v );
    LARCV_INFO() << "Cropped images prepared. TwoPlaneMode=" << _twoplane_mode << std::endl;
  }
  
  /**
   * make crops_v and mask_v data members
   *
   * @param[in] wholeview_v Whole view ADC image.
   *
   */
  void CropMaskCombo::_crop_and_mask_image( const std::vector<larcv::Image2D>& wholeview_v ) {

    // we make crops based on the union of the mask-combos
    float tick_union_min = getCombo().union_tick[0]-_tick_padding;
    float tick_union_max = getCombo().union_tick[1]+_tick_padding;
    float detz_union_min = getCombo().union_detz[0];
    float detz_union_max = getCombo().union_detz[1];

    LARCV_DEBUG() << "tick_union_max = " << tick_union_max << std::endl;

    double testpts[4][3] = { { 0,  117.0, detz_union_min  },  // minz top
                             { 0, -117.0, detz_union_min  },  // minz bottom
                             { 0,  117.0, detz_union_max  },  // maxz top
                             { 0, -117.0, detz_union_max  }}; // maxz bot

    // adjust tick bounds to be divisible by pixel_height
    int tickdiff = int(tick_union_max-tick_union_min);
    int nrows    = tickdiff/int(wholeview_v.front().meta().pixel_height());
    if ( tickdiff%int(wholeview_v.front().meta().pixel_height())!=0 )
      nrows++;
    // ensure we are a multiple of some factor we can use to downsample later (for AStar)
    if ( nrows%_downsample_factor!=0 )
      nrows += _downsample_factor-int(nrows)%_downsample_factor;
    if ( nrows%_downsample_factor!=0 )
      throw std::runtime_error( "nrows not divisible by downsampling factor" );
    tick_union_max = tick_union_min + int(nrows)*wholeview_v.front().meta().pixel_height();
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

      if ( getCombo().indices[p]==-1 ) {
        // create empty image for now
        _twoplane_mode = true;
        _badplane = (int)p;
        larcv::Image2D blank;
        crops_v.push_back(blank);
        mask_v.push_back(blank);
        continue;
      }
      
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
      if (ncols%_downsample_factor!=0)
        ncols += _downsample_factor-ncols%_downsample_factor;
      if ( ncols%_downsample_factor!=0 )
        throw std::runtime_error( "nrows not divisible by _downsample_factor" );
      
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

      if ( height<0 ) { 
        LARCV_CRITICAL() << "height=" << height << " from tick_union_max=" << tick_union_max << " tick_union_min=" << tick_union_min << std::endl;
      }
      
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
          float pttick = mask.meta.pos_y( mask.box.min_y() + row ) + mask.meta.pixel_height(); // offset
          
          int xrow = cropmeta.row( pttick );
          int xcol = cropmeta.col( ptwire );

          // set to arbitrary value, but must be big enough
          // to not be scaled down to zero when converting to cv::Mat later
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
        if ( vec_pix[idx]<_threshold ) vec_pix[idx]  = 0.0;
        if ( vec_pix[idx]<_threshold ) vec_mask[idx] = 0.0;
      }
      
      crops_v.emplace_back( std::move(crop) );
      mask_v.emplace_back( std::move(imgmask) );
    }
  }

  /**
   * if one plane does not have mask bounding box, infer from other planes and replace blank image in crops_v.
   *
   * also populates missing_v member which is deprecated
   *
   */
  void CropMaskCombo::_make_missing_crop( const std::vector<larcv::Image2D>& wholeview_v,
                                          std::vector<larcv::Image2D>& out_v,
                                          std::vector<larcv::Image2D>& miss_out_v ) {
    // we do this by using the union detz of the 2 good planes to build us the bad plane crop

    // count number of good planes
    _badplane    = -1;
    int nrows = 0;
    int tick_min = 0;
    
    for ( size_t p=0; p<crops_v.size(); p++ ) {
      auto const& crop = crops_v[p];
      if ( crop.meta().cols()==0 || crop.meta().rows()==0 ) {
        _badplane = (int)p;
        _twoplane_mode = true;
      }
      else {
        nrows = crop.meta().rows();
        tick_min = crop.meta().min_y();
      }
    }

    if ( !_twoplane_mode ) {
      // nothing to do: do not fill missing_v
      return;
    }

    // use union of detz to make crop
    auto const& range_detz = getCombo().intersection_detz;    
    //auto const& range_detz = getCombo().union_detz;
    std::vector<int> plane_wire;

    double testpts[4][3] = { { 0,  117.0, range_detz[0] },  // minz top
                             { 0, -117.0, range_detz[0] },  // minz bottom
                             { 0,  117.0, range_detz[1] },  // maxz top
                             { 0, -117.0, range_detz[1] }}; // maxz bot
    
    // find possible missing-wire bounds
    float near_wireco[4];
    int   near_wire[4];
    for ( int i=0; i<4; i++ ) {
      if ( testpts[i][2]>=1036 ) // REPLACE ME
        testpts[i][2] = 1036.0;  // REPLACE ME
      if ( testpts[i][2]<0.5 )   // REPLACE ME
        testpts[i][2] = 0.5;     // REPLACE ME
      
      near_wireco[i] = larutil::Geometry::GetME()->WireCoordinate( testpts[i], _badplane );
      near_wire[i] = (int)(near_wireco[i]+0.5);
      if ( near_wire[i]<0 ) near_wire[i] = 0;
      if ( near_wire[i]>=(int)larutil::Geometry::GetME()->Nwires(_badplane) )
        near_wire[i] = (int)larutil::Geometry::GetME()->Nwires(_badplane) - 1;
    }
    
    float start_wire = larutil::Geometry::GetME()->Nwires(_badplane)-1;
    float end_wire   = 0;
    for ( int i=0; i<4; i++ ) {
      if ( start_wire>near_wire[i] ) start_wire = near_wire[i];
      if ( end_wire<near_wire[i] )   end_wire   = near_wire[i];
    }
    
    LARCV_INFO() << "make missing crop: plane=" << _badplane
                 << " wire range=[" << start_wire << "," << end_wire << "]" << std::endl;
    
    int ncols = (end_wire-start_wire)*wholeview_v.at(_badplane).meta().pixel_width();
    if ( ncols%_downsample_factor!=0 )
      ncols += (_downsample_factor-ncols%_downsample_factor);
    end_wire = start_wire + wholeview_v.at(_badplane).meta().pixel_width()*ncols;
    if ( end_wire>wholeview_v.at(_badplane).meta().max_x() ) {
      end_wire = wholeview_v.at(_badplane).meta().max_x();
      start_wire = end_wire - ncols*wholeview_v.at(_badplane).meta().pixel_width();
    }
    if ( start_wire<0 ) {
      start_wire = 0;
      end_wire   = start_wire + wholeview_v.at(_badplane).meta().pixel_width()*ncols;
    }


    // make meta
    larcv::ImageMeta cropmeta( end_wire-start_wire, nrows*wholeview_v.at(_badplane).meta().pixel_height(),
                               nrows, ncols,
                               start_wire, tick_min,
                               wholeview_v.at(_badplane).meta().plane() );
    larcv::Image2D crop = wholeview_v.at(_badplane).crop(cropmeta);
    
    // threshold
    for ( auto& pixval : crop.as_mod_vector() ) {
      if ( pixval<_threshold ) pixval = 0.;
    }
    
    for ( size_t p=0; p<crops_v.size(); p++ ) {
      if ( (int)p==_badplane ) {
        std::swap(out_v[p],crop);
        miss_out_v.push_back( crop ); // deprecate
      }
      else {
        miss_out_v.emplace_back(larcv::Image2D()); // empty. deprecate.
      }
    }
    
  }

  /**
   * make a bad channel cropped image to match size of ADC cropped images.
   *
   * @param[in] whole_badch_v Bad channel image of the whole view for each plane.
   * @param[in] input_crop_v    The ADC cropped images. We crop bad channel image to match size.
   * @param[in] input_missing_v [deprecated]
   * @param[out] badch_crop_v   Cropped bad channel images
   */
  void CropMaskCombo::_prep_badch_crop( const std::vector<larcv::Image2D>& whole_badch_v,
                                        const std::vector<larcv::Image2D>& input_crop_v,
                                        const std::vector<larcv::Image2D>& input_missing_v,
                                        std::vector<larcv::Image2D>& badch_crop_v ) {
    
    badch_crop_v.clear();

    for ( size_t p=0; p<input_crop_v.size(); p++ ) {
      auto const& cropimg = input_crop_v[p];
      const larcv::ImageMeta* cropmeta = &cropimg.meta();
      if ( cropmeta->cols()==0 || cropmeta->rows()==0 )
        cropmeta = &input_missing_v[p].meta();
      larcv::Image2D badchcrop = whole_badch_v.at(p).crop( *cropmeta );
      badch_crop_v.emplace_back( std::move(badchcrop) );
      
    }//end of plane loop
    
  }
  
  
}
}
