
#include "MCPixelPMap.h"

#include "larcv/core/DataFormat/EventImage2D.h"


namespace ublarcvapp {
namespace mctools {

  void MCPixelPMap::buildmap( larcv::IOManager& iolcv, larlite::storage_manager& ioll ){
    ublarcvapp::mctools::MCPixelPGraph mcpg;
    mcpg.set_adc_treename(adc_tree);
    mcpg.buildgraph(iolcv, ioll);
    buildmap(iolcv, mcpg);
  }

  void MCPixelPMap::buildmap( larcv::IOManager& iolcv, ublarcvapp::mctools::MCPixelPGraph& pixGraph ){

    larcv::EventImage2D* ev_adc = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, adc_tree );
    const auto& adc_v = ev_adc->Image2DArray();

    for( const auto& node: pixGraph.node_v ){
      for( size_t p=0; p < adc_v.size(); p++ ){
        for( unsigned int iP = 0; iP < node.pix_vv[p].size()/2; iP++ ) {
          int row = (node.pix_vv[p][2*iP] - 2400)/6;
          int col = node.pix_vv[p][2*iP+1];
          PixLoc_t pixLoc( (int)p, row, col );
          PixCont_t pixCont( node.pid, adc_v[p].pixel(row, col) );
          if(pixMap.count(pixLoc) == 0){
            std::vector<PixCont_t> pixCont_v{pixCont};
            pixMap[pixLoc] = pixCont_v;
          }
          else{ pixMap[pixLoc].push_back(pixCont); }
        }
      }
    }

  }

  std::vector<MCPixelPMap::PixCont_t> MCPixelPMap::getPixContent(int plane, int row, int col){
    PixLoc_t pixLoc(plane, row, col);
    if(pixMap.count(pixLoc) > 0){ return pixMap.at(pixLoc); }
    else{
      std::vector<PixCont_t> pixCont_v;
      return pixCont_v;
    }
  }

  MCPixelPMap::PixCont_t MCPixelPMap::getPixParticle(int plane, int row, int col){
    auto pixCont_v = getPixContent(plane, row, col);
    PixCont_t maxPartPixCont;
    for(const auto& pixCont: pixCont_v){
      if(pixCont > maxPartPixCont) maxPartPixCont = pixCont;
    }
    return maxPartPixCont;
  }

  int MCPixelPMap::getPixParticlePDG(int plane, int row, int col){
    return getPixParticle(plane, row, col).pdg;
  }

}
}

