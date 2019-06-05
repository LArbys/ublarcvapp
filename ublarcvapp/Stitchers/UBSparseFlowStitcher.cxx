#include "UBSparseFlowStitcher.h"

namespace ublarcvapp {

  UBSparseFlowStitcher::UBSparseFlowStitcher( const std::vector<larcv::Image2D>& embedding_images )
    : larcv_base("UBSparseFlowStitcher")
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

  void UBSparseFlowStitcher::clear() {
    _outimg_v.clear();
    _metric_v.clear();
  }

  void UBSparseFlowStitcher::addSparseData( const larcv::SparseImage& sparse,
                                            const larcv::ImageMeta& target1_subimg,
                                            const larcv::ImageMeta& target2_subimg ) {

    // loop through features
    unsigned int nfeats = sparse.nfeatures();
    LARCV_NORMAL() << "number of features: " << nfeats << std::endl;
    unsigned int nvals  = nfeats+2;
    const std::vector<float>& pixlist = sparse.pixellist();

    // calc mean of meta (for metric)
    std::vector<float> coldist_v;
    for ( auto const& outmeta : sparse.meta_v() ) {
      coldist_v.push_back( (float)outmeta.cols()/2.0 );
    }
      
    
    for ( size_t ipt=0; ipt<sparse.pixellist().size()/nvals; ipt++ ) {
      int row = pixlist.at( ipt*nvals + 0 );
      int col = pixlist.at( ipt*nvals + 1 ); 
      for ( size_t ifeat=0; ifeat<nfeats; ifeat++ ) {
        const larcv::ImageMeta& meta = sparse.meta(ifeat);
        float src_wire = meta.pos_x(col); // source wire
        float tick     = meta.pos_y(row);
        float met  = fabs((float)col-coldist_v[ifeat]);
        int outrow,outcol;
        try {
          outrow = _outimg_v[ifeat].meta().row( tick );
          outcol = _outimg_v[ifeat].meta().col( src_wire );
        }
        catch (std::exception& e) {
          LARCV_WARNING() << "Error filling pt["  << ipt << "] (" << src_wire << ", " << tick << ")" << std::endl;
          continue;
        }

        float current_metric = _metric_v[ifeat].pixel( outrow, outcol );
        if ( current_metric<0 || met<current_metric ) {
          // update metric image and output value
          float flowout = pixlist.at(ipt*nvals+2+ifeat);
          if ( flowout>-1000.0 ) {
            const larcv::ImageMeta& target_meta = _outimg_v.at(ifeat).meta();
            float tar_subimg_col  = col + flowout;
            float tar_subimg_wire = 0;
            if ( ifeat==0 )
              tar_subimg_wire = target1_subimg.min_x() + tar_subimg_col*target1_subimg.pixel_width();
            else
              tar_subimg_wire = target2_subimg.min_x() + tar_subimg_col*target2_subimg.pixel_width();
            try {
              flowout = (float)_outimg_v[ifeat].meta().col( tar_subimg_wire ) - (float)outcol;
              _outimg_v[ifeat].set_pixel( outrow, outcol, flowout );
              _metric_v[ifeat].set_pixel( outrow, outcol, met );
            }
            catch (std::exception& e) {
              LARCV_WARNING() << "flow out of wholeview target img" << std::endl;
            }
          }
        }
        
      }//end of loop over pixel features
    }//end of pixel list
  }
  
  
}
