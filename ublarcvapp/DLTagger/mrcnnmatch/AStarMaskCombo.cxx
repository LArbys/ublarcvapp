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

    _prep_badch_crop( badch_v );
    if ( run )
      _run_astar();
  }


  AStarMaskCombo::AStarMaskCombo( const Gen3DEndpoints& endptdata,
                                  const std::vector<larcv::Image2D>& whole_badch_v,
                                  bool run )
    : pEndpoint( &endptdata ),
      astar_completed(0)
  {
    _prep_badch_crop( whole_badch_v );
    if ( run )
      _run_astar();    
  }

  

  void AStarMaskCombo::_prep_badch_crop( const std::vector<larcv::Image2D>& whole_badch_v ) {

    badch_crop_v.clear();

    const std::vector<larcv::Image2D>& crop_v = pEndpoint->pfeatures->pcropdata->crops_v;
    const std::vector<larcv::Image2D>& mask_v = pEndpoint->pfeatures->pcropdata->mask_v;

    for ( size_t p=0; p<crop_v.size(); p++ ) {
      auto const& cropimg = crop_v[p];
      larcv::Image2D badchcrop = whole_badch_v.at(p).crop( cropimg.meta() );
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

      // we also want to provide some help along the extension path
      std::cout << "make badch crop: " << badchcrop.meta().dump() << std::endl;
      badch_crop_v.emplace_back( std::move(badchcrop) );
      
    }//end of plane loop
    
  }

  // void AStarMaskCombo::_prep_astar_crop( const std::vector<larcv::Image2D>& whole_adc_v ) {
  //   // recover charge
  // }

  void AStarMaskCombo::_run_astar() {

    ublarcvapp::reco3d::AStar3DAlgoConfig config;
    config.store_score_image = true;

    ublarcvapp::reco3d::AStar3DAlgo algo( config );
    algo.setVerbose(1);

    const std::vector<larcv::Image2D>& crop_v = pEndpoint->pfeatures->pcropdata->crops_v;

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
    for ( size_t p=0; p<crop_v.size(); p++ ) {
      std::cout << "astar within crop image plane[" << p << "]: " << crop_v[p].meta().dump() << std::endl;
      start_col_v.push_back( crop_v[p].meta().col( pEndpoint->endpt_wid_v[0][p] ) );
      end_col_v.push_back(   crop_v[p].meta().col( pEndpoint->endpt_wid_v[1][p] ) );
    }

    int goal_reached = 0;
    astar_path = algo.downsampleAndFindPath( downsamplefactor,
                                             crop_v, badch_crop_v, badch_crop_v,
                                             start_row, end_row,
                                             start_col_v, end_col_v,
                                             goal_reached, 2 );
    
    astar_completed = goal_reached;
    std::cout << "astar complete: " << astar_completed << " len(astar_path)=" << astar_path.size() << std::endl;
  }

}
}
