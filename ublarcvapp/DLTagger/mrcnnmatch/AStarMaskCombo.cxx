#include "AStarMaskCombo.h"

#include "ublarcvapp/UBImageMod/EmptyChannelAlgo.h"


namespace ublarcvapp {
namespace dltagger {

  AStarMaskCombo::AStarMaskCombo( const Gen3DEndpoints& endptdata,
                                  const larcv::EventChStatus& evchstatus,
                                  bool run, int max_downsample_factor, int store_score_image,
                                  int astar_verbosity )
    : larcv::larcv_base("AStarMaskCombo"),
    pEndpoint( &endptdata ),
    astar_completed(0),
    _max_downsample_factor(max_downsample_factor),
    _store_score_image(store_score_image),
    _astar_verbosity(astar_verbosity)
  {

    // make chstatus image
    ublarcvapp::EmptyChannelAlgo badchalgo;
    std::vector<larcv::Image2D> badch_v =
      badchalgo.makeBadChImage( 4, 3, 2400, 1008*6, 3456, 6, 1, evchstatus );

    // make a copy of the crop images    
    std::vector<larcv::Image2D> input_v;
    _prep_crop_image(input_v);
    _prep_badch_crop( badch_v, input_v );
    if ( run )
      _run_astar( input_v );
  }


  AStarMaskCombo::AStarMaskCombo( const Gen3DEndpoints& endptdata,
                                  const std::vector<larcv::Image2D>& whole_badch_v,
                                  bool run, int max_downsample_factor, int store_score_image,
                                  int astar_verbosity )
    : larcv::larcv_base("AStarMaskCombo"),
    pEndpoint( &endptdata ),
    astar_completed(0),
    _max_downsample_factor(max_downsample_factor),
    _store_score_image(store_score_image),
    _astar_verbosity(astar_verbosity)    
  {

    // make a copy of the crop images    
    std::vector<larcv::Image2D> input_v;
    _prep_crop_image(input_v);
    _prep_badch_crop( whole_badch_v, input_v );
    if ( run )
      _run_astar( input_v );
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

      std::cout << "make badch crop (for ASTAR): " << cropimg.meta().dump() << std::endl;
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

    if ( config.store_score_image && algo.getScoreImages().size()>0 ) {
      for ( size_t p=0; p<3; p++ ) {
        //score_crop_v.push_back( algo.getScoreImages().at(3+p) );
        score_crop_v.push_back( algo.getScoreImages().at(p) );
      }
    }
  }

}
}
