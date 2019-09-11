#include "AStarMaskCombo.h"

#include "ublarcvapp/UBImageMod/EmptyChannelAlgo.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

namespace ublarcvapp {
namespace dltagger {

  AStarMaskCombo::AStarMaskCombo( const Gen3DEndpoints& endptdata,
                                  const larcv::EventChStatus& evchstatus,
                                  bool run, int max_downsample_factor, int store_score_image,
                                  int astar_verbosity )
    : larcv::larcv_base("AStarMaskCombo"),
    pEndpoint( &endptdata ),
    astar_completed(0),
    used_astar(false),    
    _max_downsample_factor(max_downsample_factor),
    _store_score_image(store_score_image),
    _astar_verbosity(astar_verbosity)
  {

    set_verbosity( (larcv::msg::Level_t)_astar_verbosity );
    
    // make chstatus image
    ublarcvapp::EmptyChannelAlgo badchalgo;
    std::vector<larcv::Image2D> badch_v =
      badchalgo.makeBadChImage( 4, 3, 2400, 1008*6, 3456, 6, 1, evchstatus );

    // make a copy of the crop images    
    std::vector<larcv::Image2D> input_v;
    _prep_crop_image(input_v);
    _prep_badch_crop( badch_v, input_v );

    LARCV_INFO() << "Run combo tests: starting with linear test" << std::endl;      
    float max_disc_cm = 0;
    _run_linear_test( input_v, badch_crop_v, max_disc_cm );
    if ( run ) {    
      LARCV_INFO() << "results of linear test: complete=" << astar_completed << " max_disc_cm=" << max_disc_cm << std::endl;      
      if ( astar_completed==0 && max_disc_cm < 100.0 ) {
        astar_path.clear();
        LARCV_INFO() << "running astar" << std::endl;
        _run_astar( input_v );
        LARCV_INFO() << "results of astar test: complete=" << astar_completed << std::endl;
        used_astar = true;
      }
    }
  }


  AStarMaskCombo::AStarMaskCombo( const Gen3DEndpoints& endptdata,
                                  const std::vector<larcv::Image2D>& whole_badch_v,
                                  bool run, int max_downsample_factor, int store_score_image,
                                  int astar_verbosity )
    : larcv::larcv_base("AStarMaskCombo"),
    pEndpoint( &endptdata ),
    astar_completed(0),
    used_astar(false),    
    _max_downsample_factor(max_downsample_factor),
    _store_score_image(store_score_image),
    _astar_verbosity(astar_verbosity)    
  {

    set_verbosity( (larcv::msg::Level_t)_astar_verbosity );
    
    // make a copy of the crop images    
    std::vector<larcv::Image2D> input_v;
    _prep_crop_image(input_v);
    _prep_badch_crop( whole_badch_v, input_v );

    LARCV_INFO() << "Run combo tests: starting with linear test" << std::endl;
    float max_disc_cm = 0;
    _run_linear_test( input_v, badch_crop_v, max_disc_cm );
    if ( run ) {    
      LARCV_INFO() << "results of linear test: complete=" << astar_completed << " max_disc_cm=" << max_disc_cm << std::endl;            
      if ( astar_completed==0 && max_disc_cm < 10000.0 ) {
        astar_path.clear();
        LARCV_INFO() << "running astar" << std::endl;        
        _run_astar( input_v );
        LARCV_INFO() << "results of astar test: complete=" << astar_completed << std::endl;
        used_astar = true;
      }
    }
  }

  /**
   * prepare crops to give to AStar algorithm.
   *
   * if an empty image was given for one plane, we replace it with the crop prepared in CropMaskCombo.
   *
   * @param[inout] input_v container for images we prepare. these will be passed to AStar.
   *
   */
  void AStarMaskCombo::_prep_crop_image( std::vector<larcv::Image2D>& input_v ) {
    // sources
    const std::vector<larcv::Image2D>& crop_v    = pEndpoint->pfeatures->pcropdata->crops_v;
    const std::vector<larcv::Image2D>& missing_v = pEndpoint->pfeatures->pcropdata->missing_v;

    // make a copy of the crop images    
    input_v.clear();
    for ( size_t p=0; p<crop_v.size(); p++ ) {
      if ( crop_v[p].meta().cols()==0 || crop_v[p].meta().rows()==0 ) {
        LARCV_INFO() << "[AStarMaskCombo] replace empty image with missing plane crop" << std::endl;
        auto const& missingcrop = missing_v[p];
        input_v.push_back(missingcrop);
      }
      else {
        input_v.push_back(crop_v[p]);
      }
    }

    // dialate
    // int kernelsize=8;
    // for ( size_t p=0; p<crop_v.size(); p++ ) {
    //   auto& img  = input_v[p];
    //   int nrows = (int)img.meta().rows();
    //   int cols  = (int)img.meta().cols();
    // }
    
  }

  /**
   * make a bad channel cropped image to match size of ADC cropped images.
   *
   * @param[in] whole_badch_v Bad channel image of the whole view for each plane.
   * @param[in] crop_v        The ADC cropped images. We crop bad channel image to match size.
   *
   */
  void AStarMaskCombo::_prep_badch_crop( const std::vector<larcv::Image2D>& whole_badch_v,
                                         const std::vector<larcv::Image2D>& crop_v ) {

    badch_crop_v.clear();

    // get cropped mask images as well
    const std::vector<larcv::Image2D>& mask_v = pEndpoint->pfeatures->pcropdata->mask_v;

    for ( size_t p=0; p<crop_v.size(); p++ ) {
      auto const& cropimg = crop_v[p];

      LARCV_DEBUG() << "make badch crop (for ASTAR): " << cropimg.meta().dump() << std::endl;
      larcv::Image2D badchcrop = whole_badch_v.at(p).crop( cropimg.meta() );

      // attempt to reduce badchannels using mask
      /*
      for ( size_t c=0; c<badchcrop.meta().cols(); c++ ) {
        int status = badchcrop.pixel( 0, c );

        if ( status==0 ) continue; // good, do nothing

        // if bad, we try to flip some pixels to good outside the mask
        // this is to help astar algo pass through these regions only
        for ( size_t r=0; r<badchcrop.meta().rows(); r++ ) {

          // not in mask, we label 'good' again
          if ( cropimg.pixel(r,c)!=0 )
            badchcrop.set_pixel(r,c,0.0);
          
        }
      }
      */
      // we also want to provide some help along the extension path

      badch_crop_v.emplace_back( std::move(badchcrop) );
      
    }//end of plane loop
    
  }

  /**
   * run AStar algorithm
   *
   * @param[in] crop_v Input cropped ADC images to run on.
   *
   */
  void AStarMaskCombo::_run_astar( const std::vector<larcv::Image2D>& crop_v ) {

    ublarcvapp::reco3d::AStar3DAlgoConfig config;
    config.store_score_image = _store_score_image;
    config.astar_neighborhood.resize(10,1); // tight path (2N+1)^3 nodes checked
    config.max_steps = 2000;

    ublarcvapp::reco3d::AStar3DAlgo algo( config );
    algo.setVerbose(_astar_verbosity);

    int downsamplefactor = 1;

    // need to check what factor we can use
    while (true) {
      int test_factor = downsamplefactor*2;
      bool is_divisible = true;
      for ( auto const& crop : crop_v ) {
        if ( crop.meta().cols()%test_factor!=0 )
          is_divisible = false;
        if ( crop.meta().rows()%test_factor!=0 )
          is_divisible = false;
        if ( !is_divisible )
          break;
      }
      if ( is_divisible )
        downsamplefactor = test_factor;
      else
        break;

      // stop if downsample factor high
      if ( downsamplefactor>=_max_downsample_factor )
        break;
    }

    LARCV_DEBUG() << "downsample factor: " << downsamplefactor << std::endl;
    
    LARCV_DEBUG() << "start tick=" << pEndpoint->endpt_tyz_v[0][0] << " end tick=" << pEndpoint->endpt_tyz_v[1][0] << std::endl;
    
    int start_row = crop_v.front().meta().row( pEndpoint->endpt_tyz_v[0][0] );
    int end_row   = crop_v.front().meta().row( pEndpoint->endpt_tyz_v[1][0] );
    std::vector<int> start_col_v;
    std::vector<int> end_col_v;
    bool endpoints_in_crops = true;
    for ( size_t p=0; p<crop_v.size(); p++ ) {
      // check if wire coordinate of start and end inside the crop
      int start_wid = pEndpoint->endpt_wid_v[0][p];
      int end_wid   = pEndpoint->endpt_wid_v[1][p];

      LARCV_DEBUG() << "astar within crop image plane[" << p << "]: " << crop_v[p].meta().dump() << std::endl;      
      LARCV_DEBUG() << "  (start) wire=" << start_wid << std::endl;
      LARCV_DEBUG() << "  (end)   wire=" << end_wid << std::endl;
      if ( !crop_v[p].meta().contains( start_wid, pEndpoint->endpt_tyz_v[0][0] )
           || !crop_v[p].meta().contains( end_wid, pEndpoint->endpt_tyz_v[1][0] ) ) {
        endpoints_in_crops = false;
        break;
      }
      start_col_v.push_back( crop_v[p].meta().col( start_wid ) );
      end_col_v.push_back(   crop_v[p].meta().col( end_wid ) );
    }

    if ( endpoints_in_crops ) {
      int goal_reached = 0;
      astar_path = algo.downsampleAndFindPath( downsamplefactor,
                                               crop_v, badch_crop_v, badch_crop_v,
                                               start_row, end_row,
                                               start_col_v, end_col_v,
                                               goal_reached, 2 );
      
      astar_completed = goal_reached;
      LARCV_INFO() << "astar complete: " << astar_completed << " len(astar_path)=" << astar_path.size() << std::endl;
      for ( auto const& node : astar_path ) {
        LARCV_INFO() << "  " << node.str()
                     << " wires=(" << crop_v[0].meta().pos_x( node.cols[0] ) << ","
                     << crop_v[1].meta().pos_x( node.cols[1] ) << ","
                     << crop_v[2].meta().pos_x( node.cols[2] ) << ")"
                     << std::endl;
      }
    }
    else{
      LARCV_WARNING() << "astar not run, because start/end points not in crops" << std::endl;
    }

    // if goal reached, we need to go and find the charge point the path ran on
    // if ( astar_completed ) {
    //   // this doesnt seem to work
    //   // we limit our search to the masks
    //   auto const& cropmask_v = pEndpoint->pfeatures->pcropdata->mask_v;
      
    //   for ( auto& node : astar_path ) {

    //     for ( size_t p=0; p<3; p++ ) {
    //       int pixel_radius =  downsamplefactor;
    //       int row = node.row;
    //       int col = node.cols[p];
    //       // we know we downsampled, so we have to search in the downsample x downample window and find the max pixel
    //       float maxpix = 0.;
    //       int maxcol = -1;
    //       int maxrow = -1;
    //       for ( int dr=-pixel_radius; dr<=pixel_radius; dr++ ) {
    //         int r = row+dr;
    //         if ( r<0 || r>=(int)cropmask_v[p].meta().rows() ) continue;
    //         for ( int dc=-pixel_radius; dc<=pixel_radius; dc++ ) {
    //           int c = col+dc;
    //           if ( c<0 || c>=(int)cropmask_v[p].meta().cols() ) continue;
    //           float pixval  = crop_v[p].pixel(r,c);
    //           float maskval = cropmask_v[p].pixel(r,c);
    //           if ( (pixval>10.0 && maskval>0 && pixval<100.0) && pixval>maxpix) {
    //             maxpix = pixval;
    //             maxcol = c;
    //             maxrow = r;
    //           }
    //         }//end of col loop
    //       }//end of row loop
    //       if ( maxcol>=0 && maxrow>=0 ) {
    //         node.row = maxrow;
    //         node.cols[p] = maxcol;
    //       }
    //     }// end of plane loop
        
    //   }
    // }
    
    if ( config.store_score_image && algo.getScoreImages().size()>0 ) {
      for ( size_t p=0; p<3; p++ ) {
        //score_crop_v.push_back( algo.getScoreImages().at(3+p) );
        score_crop_v.push_back( algo.getScoreImages().at(p) );
      }
    }
  }

  /**
   * we walk along a 3D line between the end points and test if it is completable
   *
   */
  void AStarMaskCombo::_run_linear_test( const std::vector<larcv::Image2D>& input_v,
                                         const std::vector<larcv::Image2D>& badch_v,
                                         float& max_discontinuity_cm ) {

    float steplen = 0.3; // cm
    int pixel_search_width = 3;
    const size_t nplanes = input_v.size();

    // we fill the astar_path container with nodes along this line
    std::map< std::vector<int>, int > pixel_map;
    std::vector< std::vector<int> >   pixel_list;
    std::vector< std::vector<float> > pixel_tyz;

    std::vector<float> dir(3,0);
    float pathlen = 0.;
    for ( int i=0; i<3; i++ ) {
      dir[i] = pEndpoint->endpt_tyz_v[1][i] - pEndpoint->endpt_tyz_v[0][i];
      if ( i==0 ) {
        // change to cm
        dir[i] = (dir[i]*0.5)*larutil::LArProperties::GetME()->DriftVelocity(); // [tick][usec/tick][cm/usec]
      }
      pathlen += dir[i]*dir[i];
    }
    pathlen = sqrt(pathlen);
    for ( int i=0; i<3; i++ ) dir[i] /= pathlen;

    int nsteps = pathlen/steplen;
    if ( fabs(pathlen-steplen*nsteps)>0.0001 )
      nsteps++;
    steplen = pathlen/float(nsteps);

    double pos[3] = {0};
    float start_x = (pEndpoint->endpt_tyz_v[0][0]-3200)*0.5*larutil::LArProperties::GetME()->DriftVelocity();
    float end_x   = (pEndpoint->endpt_tyz_v[1][0]-3200)*0.5*larutil::LArProperties::GetME()->DriftVelocity();

    // as we step through, we want to track if the segment is good (on charge) or bad (not on charge).
    // so we have to track the state and the distance. we do this with the following struct and function
    struct SegmentState_t {
      int state; // 0=bad; 1=good;
      float curr_seg_len; // current segment length
      float max_good_seg; // max good segment seen
      float max_bad_seg;  // max bad segment seen
      std::vector<float> seg_v; // list of segs we've seen
      std::vector<int>   state_v; // list of states of the seg's we've seen
      SegmentState_t()
        : state(-1), // start in good state
          curr_seg_len(0.0),
          max_good_seg(0.0),
          max_bad_seg(0.0)
      {};
      void update( int currentstate, float steplen ) {
        if ( state==-1 ) {
          // first state: set and move on
          curr_seg_len = 0;
          state = currentstate;
          return;
        }
        
        if ( currentstate==state ) {
          // continuation of state
          curr_seg_len += steplen;
          return;
        }
        else {
            
          if ( state==1 ) {
            // good to bad
            // -----------
            // transition in the middle
            curr_seg_len += 0.5*steplen;
            // was good
            if ( curr_seg_len>max_good_seg )
              max_good_seg = curr_seg_len;
          }
          else {
            // was bad
            if ( curr_seg_len>max_bad_seg )
              max_bad_seg = curr_seg_len;
          }
          // store the seg
          seg_v.push_back( curr_seg_len );
          state_v.push_back( state );
          // set the state and start segment halfway through step
          state = currentstate;
          curr_seg_len = 0.5*steplen;
        };
      };
      void finish() {
        seg_v.push_back( curr_seg_len );
        state_v.push_back( state );
      }
      float get_max_bad_seg_notfirst() {
        float maxbad = 0.;
        // skip the first, because the end point might not be very good at first
        for ( size_t s=1; s<seg_v.size(); s++ ) {
          if ( state_v[s]==0 && seg_v[s]>maxbad )
            maxbad = seg_v[s];
        }
        return maxbad;
      };
      float get_max_bad_seg() {
        float maxbad = 0.;
        // skip the first, because the end point might not very good at first
        for ( size_t s=0; s<seg_v.size(); s++ ) {
          if ( state_v[s]==0 && seg_v[s]>maxbad )
            maxbad = seg_v[s];
        }
        return maxbad;
      };
    } segstate;
    
    for (int istep=0; istep<nsteps; istep++ ) {
      pos[0] = start_x + steplen*istep*dir[0];
      pos[1] = pEndpoint->endpt_tyz_v[0][1] + steplen*istep*dir[1];
      pos[2] = pEndpoint->endpt_tyz_v[0][2] + steplen*istep*dir[2];

      // back to ticks
      float tick = pos[0]/larutil::LArProperties::GetME()->DriftVelocity()/0.5 + 3200;
      if ( tick<input_v[0].meta().min_y() || tick>=input_v[0].meta().max_y() ) {
        // not in the image. seems bad.
        segstate.update(0,steplen);
        continue;
      }
      
      // get wires
      std::vector<int> rowcol_v; // 4 coordinates (row,col on each plane...)
      rowcol_v.reserve(4);
      rowcol_v.push_back( input_v[0].meta().row(tick) );
      for ( size_t pl=0; pl<nplanes; pl++ ) {
        float wirecoord = larutil::Geometry::GetME()->WireCoordinate( pos, pl );
        if ( wirecoord<input_v[pl].meta().min_x() || wirecoord>=input_v[pl].meta().max_x() ) {
          break;
        }
        rowcol_v.push_back( (int)input_v[pl].meta().col( wirecoord ) );
      }

      if ( rowcol_v.size()!=(1+nplanes) ) {
        // bad point by out of bounds step
        segstate.update(0,steplen);
        continue;
      }

      // we have a good point in the image
      auto it = pixel_map.find( rowcol_v );
      if ( it!=pixel_map.end() ) {
        // already in the set, no need to revaluate
        // all points in the set are good (good idea?)
        segstate.update(it->second,steplen);
        continue;
      }

      // new location in the image
      // get the pixel values
      int nplanes_w_charge = 0;
      int nplanes_w_badch  = 0;
      std::vector<float> step_tyz = { tick, (float)pos[1], (float)pos[2] };

      for ( size_t pl=0; pl<nplanes; pl++ ) {

        float maxpixval = 0.;
        float maxbadval = 0.;
        for ( int dc=-pixel_search_width; dc<=pixel_search_width; dc++ ) {
          int c = rowcol_v[pl+1]+dc;
          if ( c<0 || c>=(int)input_v[pl].meta().cols() ) continue; // skip it
          float pixval = input_v[pl].pixel( rowcol_v[0], c );
          float badval = badch_v[pl].pixel( rowcol_v[0], c );
          if ( pixval>maxpixval )
            maxpixval = pixval;
          if ( badval>maxbadval )
            maxbadval = badval;
        }
        
        if ( maxpixval>10.0 )
          nplanes_w_charge++;
        if ( maxbadval>0 ) 
          nplanes_w_badch++;
      }
      
      if ( nplanes_w_charge==3 ) {
        segstate.update(1,steplen);
        pixel_map[rowcol_v] = 1;
        pixel_list.push_back( rowcol_v );
        pixel_tyz.push_back( step_tyz );
      }
      else if ( (nplanes_w_charge + nplanes_w_badch) >= 3 && nplanes_w_badch<2 ) {
        segstate.update(1,steplen);
        pixel_map[rowcol_v] = 1;
        pixel_list.push_back( rowcol_v );
        pixel_tyz.push_back( step_tyz );          
      }
      else {
        // everything else is bad
        pixel_map[rowcol_v] = 0;
        segstate.update(0,steplen);
        pixel_list.push_back( rowcol_v );
        pixel_tyz.push_back( step_tyz );
      }
    }//end of step loop
    segstate.finish();

    // make astar node list
    astar_path.resize( pixel_list.size()+2 );

    // start node
    astar_path.front().row = input_v[0].meta().row(  pEndpoint->endpt_tyz_v[0][0] );
    astar_path.front().tyz = pEndpoint->endpt_tyz_v[0];
    astar_path.back().row = input_v[0].meta().row(  pEndpoint->endpt_tyz_v[1][0] );    
    astar_path.back().tyz  = pEndpoint->endpt_tyz_v[1];
    double start_pos[3] = { start_x, pEndpoint->endpt_tyz_v[0][1], pEndpoint->endpt_tyz_v[0][2] };
    double end_pos[3]   = { end_x,   pEndpoint->endpt_tyz_v[1][1], pEndpoint->endpt_tyz_v[1][2] };
    for ( size_t pl=0; pl<nplanes; pl++ ) {

      // start
      float start_wirecoord = larutil::Geometry::GetME()->WireCoordinate( start_pos, pl );
      if (start_wirecoord<input_v[pl].meta().min_x() )
        start_wirecoord = input_v[pl].meta().min_x();
      if ( start_wirecoord>=input_v[pl].meta().max_x() )
        start_wirecoord = input_v[pl].meta().max_x()-1;
      astar_path.front().cols[pl] = input_v[pl].meta().col( start_wirecoord );

      // end
      float end_wirecoord = larutil::Geometry::GetME()->WireCoordinate( end_pos, pl );
      if (end_wirecoord<input_v[pl].meta().min_x() )
        end_wirecoord = input_v[pl].meta().min_x();
      if ( end_wirecoord>=input_v[pl].meta().max_x() )
        end_wirecoord = input_v[pl].meta().max_x()-1;
      astar_path.back().cols[pl] = input_v[pl].meta().col( end_wirecoord );      
      
    }

    for ( size_t ipt=0; ipt<pixel_list.size(); ipt++ ) {
      for ( size_t p=0; p<3; p++ )
        astar_path[ipt+1].cols[p] = pixel_list[ipt][p+1];
      astar_path[ipt+1].row = pixel_list[ipt][0];
      astar_path[ipt+1].tyz = pixel_tyz[ipt];
    }

    max_discontinuity_cm = segstate.get_max_bad_seg_notfirst();
    if ( max_discontinuity_cm>20.0  ) 
      astar_completed = 0;
    else
      astar_completed = 1;
    
    // return this
    LARCV_DEBUG() << "-----------------------------------------" << std::endl;    
    LARCV_INFO() << "linear test. max_discontinuity_cm=" << max_discontinuity_cm << "  nsteps=" << astar_path.size() << " completed=" << astar_completed << std::endl;
    if ( logger().debug() ) {
      for ( auto const& node : astar_path )
        LARCV_DEBUG() << node.str() << std::endl;
      LARCV_DEBUG() << "-----------------------------------------" << std::endl;
    }
  }
  
}
}
