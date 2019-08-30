#ifndef __ASTAR3DALGOCONFIG_CXX__
#define __ASTAR3DALGOCONFIG_CXX__

#include "AStar3DAlgoConfig.h"

namespace ublarcvapp {
namespace reco3d {

  AStar3DAlgoConfig::AStar3DAlgoConfig()
    : astar_threshold(3,10.0),
      astar_neighborhood(3,10)
  {
    astar_start_padding = 5;
    astar_end_padding   = 5;
    lattice_padding = 5;
    accept_badch_nodes = true;
    min_nplanes_w_hitpixel = 3;
    min_nplanes_w_charge   = 2;
    restrict_path = true;
    path_restriction_radius = 50.0;
    compression_mode  = 2;
    verbosity = 0;
    store_score_image = false;
    max_steps = 1000;
  };

  AStar3DAlgoConfig AStar3DAlgoConfig::MakeFromPSet( const larcv::PSet& pset ) {
    AStar3DAlgoConfig cfg;

    cfg.astar_threshold        = pset.get< std::vector<float> >( "PixelThresholds" );
    cfg.astar_neighborhood     = pset.get< std::vector<int> >( "NeighborhoodSize" );
    cfg.astar_start_padding    = pset.get< int >( "StartPadding" );
    cfg.astar_end_padding      = pset.get< int >( "EndPadding" );
    cfg.lattice_padding        = pset.get< int >( "LatticePadding" );
    cfg.accept_badch_nodes     = pset.get< bool >( "AcceptBadChannelNodes" );
    cfg.min_nplanes_w_hitpixel = pset.get< int >( "MinNumPlanesWithHitPixel" );
    cfg.min_nplanes_w_charge   = pset.get< int >( "MinNumPlanesWithCharge" );
    cfg.compression_mode       = pset.get< int >( "CompressionMode" );
    cfg.restrict_path          = pset.get< bool >( "RestrictPath", false );
    cfg.verbosity              = pset.get< int >( "Verbosity" );
    if ( cfg.restrict_path ) {
      cfg.path_restriction_radius = pset.get<float>("PathRestrictionRadius");
    }
    else
      cfg.path_restriction_radius = pset.get<float>("PathRestrictionRadius",0.0);
    cfg.store_score_image = pset.get<bool>("StoreScoreImage",false);
    cfg.max_steps              = pset.get< int >( "MaxSteps", 1000 );

    return cfg;
  }

  void AStar3DAlgoConfig::dump() {
    std::cout << "--------------------------------" << std::endl;
    std::cout << " AStar3DAlgo Config" << std::endl;
    std::cout << " Thresholds/Neighborhoods" << std::endl;
    for (int p=0; p<(int)astar_threshold.size(); p++) {
      std::cout << "  plane " << p << ": Threshold=" << astar_threshold[p] << " Neighborhood=" << astar_neighborhood[p] << std::endl;
    }
    std::cout << " Start Padding: " << astar_start_padding << std::endl;
    std::cout << " End Padding: " << astar_end_padding << std::endl;
    std::cout << " Lattice Padding: " << lattice_padding << std::endl;
    std::cout << " Accept BadCh Nodes: " << accept_badch_nodes << std::endl;
    std::cout << " Min Nplanes w/ Valid Pixel: "   << min_nplanes_w_hitpixel << std::endl;
    std::cout << " Min Nplanes w/ Charged Pixel: " << min_nplanes_w_charge << std::endl;    
    std::cout << " Restrict Path: " << restrict_path << std::endl;
    std::cout << " Restrict Radius: " << path_restriction_radius << std::endl;
    std::cout << " Compression Mode: " << compression_mode << std::endl;
    std::cout << " Verbosity: " << verbosity << std::endl;
    std::cout << " Store Score Image: " << store_score_image << std::endl;
    std::cout << " Max Steps: " << max_steps << std::endl;
    std::cout << "--------------------------------" << std::endl;
  }

}
}
#endif
