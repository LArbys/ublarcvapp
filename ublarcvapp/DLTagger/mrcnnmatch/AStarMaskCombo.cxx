#include "AStarMaskCombo.h"

#include "ublarcvapp/UBImageMod/EmptyChannelAlgo.h"


namespace ublarcvapp {
namespace dltagger {

  AStarMaskCombo::AStarMaskCombo( const Gen3DEndpoints& endptdata,
                                  const larcv::EventChStatus& evchstatus,
                                  bool run )
    : pEndpoint( &endptdata ),
      astar_completed(0)
  {

    // make chstatus image
    ublarcvapp::EmptyChannelAlgo badchalgo;
    std::vector<larcv::Image2D> badch_v =
      badchalgo.makeBadChImage( 4, 3, 2400, 1008*6, 3456, 6, 1, evchstatus );

    // sources
    const std::vector<larcv::Image2D>& crop_v    = pEndpoint->pfeatures->pcropdata->crops_v;
    const std::vector<larcv::Image2D>& missing_v = pEndpoint->pfeatures->pcropdata->crops_v;

    // make a copy of the crop images    
    std::vector<larcv::Image2D> input_v;
    for ( size_t p=0; p<crop_v.size(); p++ ) {
      if ( crop_v[p].meta().cols()==0 || crop_v[p].meta().rows()==0 ) {
        // replace empty image with missing
        std::cout << "[AStarMaskCombo] replace empty image with missing plane crop"<< std::endl;        
        auto const& missingcrop = missing_v[p];
        input_v.push_back(missingcrop);
      }
      else {
        input_v.push_back(crop_v[p]);
      }
    }
    
    _prep_badch_crop( badch_v, input_v );
    if ( run )
      _run_astar( input_v );
  }


  AStarMaskCombo::AStarMaskCombo( const Gen3DEndpoints& endptdata,
                                  const std::vector<larcv::Image2D>& whole_badch_v,
                                  bool run )
    : pEndpoint( &endptdata ),
      astar_completed(0)
  {

    // sources
    const std::vector<larcv::Image2D>& crop_v    = pEndpoint->pfeatures->pcropdata->crops_v;
    const std::vector<larcv::Image2D>& missing_v = pEndpoint->pfeatures->pcropdata->missing_v;

    // make a copy of the crop images    
    std::vector<larcv::Image2D> input_v;
    for ( size_t p=0; p<crop_v.size(); p++ ) {
      if ( crop_v[p].meta().cols()==0 || crop_v[p].meta().rows()==0 ) {
        std::cout << "[AStarMaskCombo] replace empty image with missing plane crop" << std::endl;
        auto const& missingcrop = missing_v[p];
        input_v.push_back(missingcrop);
      }
      else {
        input_v.push_back(crop_v[p]);
      }
    }
    
    _prep_badch_crop( whole_badch_v, input_v );
    if ( run )
      _run_astar( input_v );
  }

  

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

  void AStarMaskCombo::_run_astar( const std::vector<larcv::Image2D>& crop_v ) {

    ublarcvapp::reco3d::AStar3DAlgoConfig config;
    config.store_score_image = true;

    ublarcvapp::reco3d::AStar3DAlgo algo( config );
    algo.setVerbose(1);

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
      if ( downsamplefactor==16 )
        break;
    }

    std::cout << "downsample factor: " << downsamplefactor << std::endl;
    std::cout << "start tick=" << pEndpoint->endpt_tyz_v[0][0] << " end tick=" << pEndpoint->endpt_tyz_v[1][0] << std::endl;
    int start_row = crop_v.front().meta().row( pEndpoint->endpt_tyz_v[0][0] );
    int end_row   = crop_v.front().meta().row( pEndpoint->endpt_tyz_v[1][0] );
    std::vector<int> start_col_v;
    std::vector<int> end_col_v;
    bool endpoints_in_crops = true;
    for ( size_t p=0; p<crop_v.size(); p++ ) {
      // check if wire coordinate of start and end inside the crop
      int start_wid = pEndpoint->endpt_wid_v[0][p];
      int end_wid   = pEndpoint->endpt_wid_v[1][p];

      std::cout << "astar within crop image plane[" << p << "]: " << crop_v[p].meta().dump() << std::endl;      
      std::cout << "  (start) wire=" << start_wid << std::endl;
      std::cout << "  (end)   wire=" << end_wid << std::endl;
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
      std::cout << "astar complete: " << astar_completed << " len(astar_path)=" << astar_path.size() << std::endl;
      for ( auto const& node : astar_path ) {
        std::cout << "  " << node.str()
                  << " wires=(" << crop_v[0].meta().pos_x( node.cols[0] ) << ","
                  << crop_v[1].meta().pos_x( node.cols[1] ) << ","
                  << crop_v[2].meta().pos_x( node.cols[2] ) << ")"
                  << std::endl;
      }
    }
    else{
      std::cout << "astar not run, because start/end points not in crops" << std::endl;
    }

    if ( config.store_score_image && algo.getScoreImages().size()>0 ) {
      for ( size_t p=0; p<3; p++ ) {
        score_crop_v.push_back( algo.getScoreImages().at(3+p) );
      }
    }
  }

}
}
