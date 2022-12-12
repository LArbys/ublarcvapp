
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
          PixPart_t pixPart( node.nodeidx, node.tid, node.pid );
          if(pixMap.count(pixLoc) == 0){
            PixCont_t pixCont( adc_v[p].pixel(row, col) );
            pixCont.particles.push_back(pixPart);
            pixMap[pixLoc] = pixCont;
          }
          else{ pixMap[pixLoc].particles.push_back(pixPart); }
        }
      }
    }

  }

  MCPixelPMap::PixCont_t MCPixelPMap::getPixContent(int plane, int row, int col){
    PixLoc_t pixLoc(plane, row, col);
    if(pixMap.count(pixLoc) > 0){ return pixMap.at(pixLoc); }
    else{
      PixCont_t pixCont;
      return pixCont;
    }
  }

}
}

