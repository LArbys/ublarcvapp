#include "UBSparseFlowStitcher.h"

namespace ublarcvapp {

  UBSparseFlowStitcher::UBSparseFlowStitcher( const std::vector<larcv::Image2D>& embedding_images )
  {
    // We create two images. One to hold vaues, another to hold a metric
    for ( auto const& img : embedding_images ) {
      larcv::Image2D out(img.meta());
      out.paint(0.0);
      _outimg_v.emplace_back( std::move(out) );

      larcv::Image2D metric(img.meta());
      metric.paint(-1.0);
      _metric_v.emplace_back( std::move(metric) );
    }
  }

  UBSparseFlowStitcher::~UBSparseFlowStitcher()
  {
  }


  void UBSparseFlowStitcher::addSparseData( const larcv::SparseImage& sparse ) {
    // loop through features
    unsigned int nfeats = sparse.nfeatures();
    unsigned int nvals  = nfeats-2;
    const std::vector<float>& pixlist = sparse.pixellist();

    // calc mean of meta (for metric)
    std::vector<float> coldist_v;
    for ( auto const& outmeta : sparse.meta_v() ) {
      coldist_v.push_back( (float)outmeta.cols()/2.0 );
    }
      
    
    for ( size_t ipt=0; ipt<sparse.pixellist().size()/nfeats; ipt++ ) {
      int row = pixlist.at( ipt*nfeats + 0 );
      int col = pixlist.at( ipt*nfeats + 1 ); 
      for ( size_t ifeat=0; ifeat<nvals; ifeat++ ) {
        const larcv::ImageMeta& meta = sparse.meta(ifeat);
        float wire = meta.pos_x(col);
        float tick = meta.pos_y(row);
        float met  = fabs((float)col-coldist_v[ifeat]);
        int outrow,outcol;
        try {
          outrow = _outimg_v[ifeat].meta().row( tick );
          outcol = _outimg_v[ifeat].meta().col( wire );
        }
        catch (std::exception& e) {
          continue;
        }

        float current_metric = _metric_v[ifeat].pixel( outrow, outcol );
        if ( current_metric<0 || met<current_metric ) {
          // update metric image
          _outimg_v[ifeat].set_pixel( outrow, outcol, pixlist.at(ipt*nfeats+2+ifeat) );
          _metric_v[ifeat].set_pixel( outrow, outcol, met );
        }
        
      }//end of loop over pixel features
    }//end of pixel list
  }
  
  
}
