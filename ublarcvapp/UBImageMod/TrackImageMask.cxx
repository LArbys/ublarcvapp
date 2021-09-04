#include "TrackImageMask.h"

#include "larlite/LArUtil/Geometry.h"
#include "larlite/LArUtil/LArProperties.h"
#include "larcv/core/DataFormat/ImageMeta.h"

#include <chrono>

namespace ublarcvapp {
namespace ubimagemod {

  void TrackImageMask::clear()
  {
    pixel_v.clear();
    pixel_map.clear();
  }
  
  int TrackImageMask::maskTrack( const larlite::track& track,
                                 const larcv::Image2D& img,
                                 larcv::Image2D& mask,
                                 const float threshold,
                                 const int dcol,
                                 const int drow,
                                 const float maxstepsize )
  {

    // clear
    clear();
    
    int npts = track.NumberTrajectoryPoints();    
    LARCV_DEBUG() << "start: track length=" << npts << std::endl;
    
    if ( npts<=1 ) {
      LARCV_WARNING() << "No mask generated for track with only 1 point" << std::endl;
      // no mask can be generated in this case
      return 0;
    }

    auto const& meta = img.meta();

    int npixels = makePixelList( track, img, threshold, maxstepsize, false );
    LARCV_DEBUG() << "Number of pixels: " << pixel_v.size() << std::endl;    
    if ( npixels==0 )
      return 0;
        
    // Now mask

    // crop width with padding added
    int colwidth = (max_col-min_col+1) + 2*dcol;
    int rowwidth = (max_row-min_row+1) + 2*drow;
    int rowlen   = (max_row-min_row+1);
    int col_origin = min_col-dcol;
    int row_origin = min_row-drow;
    LARCV_DEBUG() << "colbounds: [" << min_col << "," << max_col << "]" << std::endl;
    LARCV_DEBUG() << "rowbounds: [" << min_row << "," << max_row << "]" << std::endl;      
    LARCV_DEBUG() << "Cropped meta. Origin=(" << col_origin << "," << row_origin << ")"
                  << " ncol=" << colwidth << " nrow=" << rowwidth
                  << std::endl;

    if ( colwidth<=0 || rowwidth<=0 ) {
      LARCV_WARNING() << "Bad cropped window size. colwidth=" << colwidth << " rowwidth=" << rowwidth << std::endl;
      return 0;
    }

    // make data array
    std::vector<float> crop(colwidth*rowwidth,0);
    std::vector<float> cropmask(colwidth*rowwidth,0);

    // copy the data into the cropped matrix
    // preserve the row-ordered data
    for (int c=0; c<=(int)(max_col-min_col); c++) {
      // for the larcv image, the row-dimension is the contiguous element
      size_t origin_index = meta.index( min_row, min_col+c, __FILE__, __LINE__ );
      const float* source_ptr = img.as_vector().data() + origin_index;
      float* dest_ptr   = crop.data() + (c+dcol)*rowwidth + drow;
      memcpy( dest_ptr, source_ptr, sizeof(float)*rowlen );
    }

    // kernal loop, making mask
    // trying to write in a way that is easy to accelerate
    int N = (2*dcol+1)*(2*drow+1); // number of elements in kernel

    // parallelize masking with block kernel
    int nmasked = 0;
    std::chrono::steady_clock::time_point start_mask = std::chrono::steady_clock::now();    
    
    #pragma omp parallel for         
    for (int ikernel=0; ikernel<N; ikernel++) {

      int dr = ikernel%(2*drow+1) - drow;
      int dc = ikernel/(2*drow+1) - dcol;          
      
      for (size_t itp=0; itp<pixel_v.size(); itp++) {
        auto const& pix = pixel_v[itp];
        
        // change (col,row) from original image to mask image (with padding)
        // then add kernal shift
        int icol = (int)pix[0]-col_origin + dc;
        int irow = (int)pix[1]-row_origin + dr;
        if ( icol>=colwidth || irow>=rowwidth || icol<0 || irow<0) {
          std::stringstream ss;
          ss << "OUT OF CROP BOUNDS (" << icol << "," << irow << ") TP[" << itp << "] orig=(" << pix[0] << "," << pix[1] << ")" << std::endl;
          throw std::runtime_error(ss.str());
        }

        int cropindex = icol*rowwidth+irow;
        float pixval = crop[ cropindex ];

        if ( pixval>threshold ) {
          #pragma omp atomic
          cropmask[ cropindex ] += 1.0;
        }
      }
    }  
    #pragma omp barrier

    for (size_t ii=0; ii<cropmask.size(); ii++ ) {
      if ( cropmask[ii]>0 ) {
        nmasked++;
        cropmask[ii] = 1.0;

        int orig_row = (ii%rowwidth)-drow + row_origin;
        int orig_col = (ii/rowwidth)-dcol + col_origin;
        mask.set_pixel( orig_row, orig_col, 1.0 );
      }
    }
    LARCV_INFO() << "num pixels masked: " << nmasked << std::endl;
    std::chrono::steady_clock::time_point end_mask = std::chrono::steady_clock::now();
    if ( show_timing )
      LARCV_INFO() << "mask: " << std::chrono::duration_cast<std::chrono::microseconds>(end_mask - start_mask).count() << std::endl;
    return nmasked;
  }

  /**
   * @brief label the pixels with the minimum and maximum path along the track
   *
   */
  int TrackImageMask::labelTrackPath( const larlite::track& track,
                                      const larcv::Image2D& img, 
                                      larcv::Image2D& smin_img,
                                      larcv::Image2D& smax_img,
                                      const float threshold,
                                      const float maxstepsize )
  {

    int npts = track.NumberTrajectoryPoints();    
    LARCV_DEBUG() << "num of trajectory points=" << npts << std::endl;
    
    if ( npts<=1 ) {
      LARCV_WARNING() << "No mask generated for track with only 1 point" << std::endl;
      // no mask can be generated in this case
      return 0;
    }

    int npixels = makePixelList( track, img, threshold, maxstepsize, true );
    LARCV_DEBUG() << "Number of pixels in pixel_map: " << pixel_map.size() << std::endl;
    
    if ( npixels==0 )
      return 0;

    int nlabeled = 0;
    for (auto it=pixel_map.begin(); it!=pixel_map.end();  it++) {
      smin_img.set_pixel( it->second.row, it->second.col, it->second.smin );
      smax_img.set_pixel( it->second.row, it->second.col, it->second.smax );
      nlabeled++;
    }

    return nlabeled;
  }

  int TrackImageMask::makePixelList( const larlite::track& track,
                                     const larcv::Image2D& img,
                                     const float threshold,
                                     const float maxstepsize,                                     
                                     bool fill_map )
  {
    
    pixel_v.clear();
    pixel_map.clear();

    if ( maxstepsize<0 || std::isnan(maxstepsize) || std::isinf(maxstepsize) ) {
      LARCV_CRITICAL() << "Bad maxstepsize value: " << maxstepsize << std::endl;
    }

    std::chrono::steady_clock::time_point start_pixlist = std::chrono::steady_clock::now();
    
    int npts = track.NumberTrajectoryPoints();    
    LARCV_DEBUG() << "track length=" << npts << std::endl;
    if ( npts<=1 ) {
      LARCV_WARNING() << "No mask generated for track with only 1 point" << std::endl;
      // no mask can be generated in this case
      return 0;
    }

    const float driftv = larutil::LArProperties::GetME()->DriftVelocity();
    const float usec_per_tick = 0.5;
    auto const& meta = img.meta();
    int plane = meta.plane();
    if ( plane<0 || plane>=(int)larutil::Geometry::GetME()->Nplanes() ) {
      LARCV_WARNING() << "Invalid plane: " << plane << std::endl;
      return 0;
    }

    // microboone only
    const std::vector<Double_t>& firstwireproj = larutil::Geometry::GetME()->GetFirstWireProj(); 
    std::vector<double> orthovect = { 0,
                                      larutil::Geometry::GetME()->GetOrthVectorsY().at(plane),
                                      larutil::Geometry::GetME()->GetOrthVectorsZ().at(plane) };
    
    int nrows = meta.rows();
    int ncols = meta.cols();

    min_col = (int)meta.cols()-1;
    min_row = (int)meta.rows()-1;
    max_col = 0;
    max_row = 0;    

    std::set< std::pair<int,int> > pixel_set;
    pixel_v.reserve( 2*npts );

    float len_traveled = 0.;    
    for (int ipt=0; ipt<npts-1; ipt++) {

      TVector3 start = track.LocationAtPoint(ipt);
      TVector3 end   = track.LocationAtPoint(ipt+1);
      TVector3 dir   = end-start;
      
      double segsize = dir.Mag();

      int nsteps = 1;
      if ( segsize>maxstepsize ) {
        nsteps = segsize/maxstepsize + 1;
      }
      
      float stepsize = segsize/float(nsteps);

      for (int istep=0; istep<=nsteps; istep++) {
        // get 3d position along track
        TVector3 pos = start + istep*(stepsize/segsize)*dir;
        // project into image        
        float tick = pos[0]/driftv/usec_per_tick + 3200; // tick
        float rowcoord = (tick-meta.min_y())/meta.pixel_height();

        //  length on track
        float s = len_traveled + istep*stepsize;
        
        // from larlite Geometry::WireCoordinate(...)
        float wirecoord = pos[1]*orthovect[1] + pos[2]*orthovect[2] - firstwireproj.at(plane);
        float colcoord = (wirecoord-meta.min_x())/meta.pixel_width();
        
        int row = (int)rowcoord;
        int col = (int)colcoord;
        
        if ( row<0 || row>=nrows ) continue;
        if ( col<0 || col>=ncols ) continue;

        if ( min_row>row ) min_row = row;
        if ( max_row<row ) max_row = row;
        if ( min_col>col ) min_col = col;
        if ( max_col<col ) max_col = col;

        float pixval = img.pixel(row,col);
        if ( pixval<threshold )
          continue;
        
        std::pair<int,int> pixcoord(col,row);

        if ( fill_map ) {
          auto it = pixel_map.find( pixcoord );
          if ( it==pixel_map.end() ) {
            // new entry
            pixel_map[pixcoord] = Pix_t( col, row, s, s );
            pixel_v.push_back( std::vector<int>{col,row} );            
            it = pixel_map.find( pixcoord );
          }
          it->second.smax = s;
        }
        else {          
          auto it = pixel_set.find( pixcoord );
          if ( it==pixel_set.end() ) {
            pixel_set.insert( std::pair<int,int>(col,row) );
            pixel_v.push_back( std::vector<int>{col,row} );
          }
        }

      }//end of pixel steps
      
      len_traveled += segsize;
      
    }//end of trajectory point list
        
    LARCV_DEBUG() << "Size of pixel_v: " << pixel_v.size() << " len-traveled=" << len_traveled << " cm" << std::endl;

    std::chrono::steady_clock::time_point end_pixlist = std::chrono::steady_clock::now();
    LARCV_DEBUG() << "time to make pixlist: " << std::chrono::duration_cast<std::chrono::microseconds>(end_pixlist - start_pixlist).count() << " usec" << std::endl;
    return pixel_v.size();
  }
  
  float TrackImageMask::sumOverPixelList( const std::vector< std::vector<int> >& pixlist_v,
                                          const larcv::Image2D& img,
                                          const int dcol,
                                          const int drow,
                                          const float thresh ) const
  {

    auto const& meta = img.meta();
    int nrows = meta.rows();
    int ncols = meta.cols();
    int rowbound[2] = { nrows, 0 };
    int colbound[2] = { ncols, 0 };

    for (auto const& pix : pixlist_v ) {
      const int& col = pix[0];
      const int& row = pix[1];
    
      if ( row<0 || row>=nrows ) continue;
      if ( col<0 || col>=ncols ) continue;

      if ( rowbound[0]>row ) rowbound[0] = row;
      if ( rowbound[1]<row ) rowbound[1] = row;
      if ( colbound[0]>col ) colbound[0] = col;
      if ( colbound[1]<col ) colbound[1] = col;
    }
    rowbound[0] -= drow;
    rowbound[1] += drow;
    colbound[0] -= dcol;
    colbound[1] += dcol;
    if ( rowbound[0]<0 ) rowbound[0] = 0;
    if ( rowbound[1]>=nrows ) rowbound[1] = nrows-1;
    if ( colbound[0]<0 ) colbound[0] = 0;
    if ( colbound[1]>=ncols ) colbound[1] = ncols-1;
    
    // we make a crop with padding to avoid bound checks
    // crop width with padding added
    int colwidth = (colbound[1]-colbound[0]+1) + 2*dcol;
    int rowwidth = (rowbound[1]-rowbound[0]+1) + 2*drow;
    int rowlen   = (rowbound[1]-rowbound[0]+1);
    int col_origin = colbound[0]-dcol;
    int row_origin = rowbound[0]-drow;

    LARCV_DEBUG() << "colbounds: [" << colbound[0] << "," << colbound[1] << "]" << std::endl;
    LARCV_DEBUG() << "rowbounds: [" << rowbound[0] << "," << rowbound[1] << "]" << std::endl;      
    LARCV_DEBUG() << "Cropped meta. Origin=(" << col_origin << "," << row_origin << ")"
                  << " ncol=" << colwidth << " nrow=" << rowwidth << " rowlen=" << rowlen
                  << std::endl;
    
    if ( colwidth<=0 || rowwidth<=0 ) {
      LARCV_WARNING() << "Bad cropped window size. colwidth=" << colwidth << " rowwidth=" << rowwidth << std::endl;
      return 0;
    }

    std::vector<float> crop(colwidth*rowwidth,0);
    std::vector<float> cropmask(colwidth*rowwidth,0);

    // copy the data into the cropped matrix
    // preserve the row-ordered data
    for (int c=0; c<=(int)(colbound[1]-colbound[0]); c++) {
      // for the larcv image, the row-dimension is the contiguous element
      size_t origin_index = meta.index( rowbound[0], colbound[0]+c, __FILE__, __LINE__ );
      const float* source_ptr = img.as_vector().data() + origin_index;
      float* dest_ptr   = crop.data() + (c+dcol)*rowwidth + drow;
      memcpy( dest_ptr, source_ptr, sizeof(float)*rowlen );
      // for (int r=0;  r<rowlen; r++)      
      //   std::cout << "img(" << rowbound[0]+r << "," << colbound[0]+c << ") " << img.pixel( rowbound[0]+r, colbound[0]+c ) << std::endl;
    }
        
    // kernal loop, making mask
    // trying to write in a way that is easy to accelerate
    int N = (2*dcol+1)*(2*drow+1); // number of elements in kernel

    // parallelize masking with block kernel
    std::chrono::steady_clock::time_point start_mask = std::chrono::steady_clock::now();    

    for (int ikernel=0; ikernel<N; ikernel++) {

      int dr = ikernel%(2*drow+1) - drow;
      int dc = ikernel/(2*drow+1) - dcol;          

      //#pragma omp parallel for
      for (size_t itp=0; itp<pixlist_v.size(); itp++) {

        auto const& col = pixlist_v[itp][0];        
        auto const& row = pixlist_v[itp][1];
        
        // change (col,row) from original image to mask image (with padding)
        // then add kernal shift 
        int icol = (int)col-col_origin + dc;
        int irow = (int)row-row_origin + dr;
        if ( icol>=colwidth || irow>=rowwidth || icol<0 || irow<0) {
          std::stringstream ss;
          ss << "OUT OF CROP BOUNDS (" << icol << "," << irow << ") TP[" << itp << "] orig=(" << col << "," << row << ")" << std::endl;
          throw std::runtime_error(ss.str());
        }

        int cropindex = icol*rowwidth+irow;
        float pixval = crop[ cropindex ];
        //std::cout << "(" << icol << "," << irow << ") d=(" << dc << "," << dr << "): " << pixval << std::endl;

        if ( pixval>thresh ) {
          //#pragma omp atomic
          cropmask[ cropindex ] += 1.0;
        }
      }
      //#pragma omp barrier      
    }

    float sum = 0.;
    for (size_t i=0; i<cropmask.size(); i++) {

      if ( cropmask[i]>0 ) {
        sum += crop[i];
        //std::cout << " [*] " << crop[i] << std::endl;
      }
      else {
        //std::cout << crop[i] << std::endl;
      }
    }

    std::chrono::steady_clock::time_point end_mask = std::chrono::steady_clock::now();
    if ( show_timing )
      LARCV_INFO() << "sum and mask: " << std::chrono::duration_cast<std::chrono::microseconds>(end_mask - start_mask).count() << std::endl;
    
    return sum;
    
  }

  bool TrackImageMask::getExtremaPixelValues( const std::vector< std::vector<int> >& pixlist_v,
                                              const larcv::Image2D& img,
                                              const int dcol,
                                              const int drow,
                                              const float thresh,
                                              float& pixval_min,
                                              float& pixval_max ) const
  {

    auto const& meta = img.meta();
    int nrows = meta.rows();
    int ncols = meta.cols();
    int rowbound[2] = { nrows, 0 };
    int colbound[2] = { ncols, 0 };

    for (auto const& pix : pixlist_v ) {
      const int& col = pix[0];
      const int& row = pix[1];
    
      if ( row<0 || row>=nrows ) continue;
      if ( col<0 || col>=ncols ) continue;

      if ( rowbound[0]>row ) rowbound[0] = row;
      if ( rowbound[1]<row ) rowbound[1] = row;
      if ( colbound[0]>col ) colbound[0] = col;
      if ( colbound[1]<col ) colbound[1] = col;
    }
    rowbound[0] -= drow;
    rowbound[1] += drow;
    colbound[0] -= dcol;
    colbound[1] += dcol;
    if ( rowbound[0]<0 ) rowbound[0] = 0;
    if ( rowbound[1]>=nrows ) rowbound[1] = nrows-1;
    if ( colbound[0]<0 ) colbound[0] = 0;
    if ( colbound[1]>=ncols ) colbound[1] = ncols-1;
    
    // we make a crop with padding to avoid bound checks
    // crop width with padding added
    int colwidth = (colbound[1]-colbound[0]+1) + 2*dcol;
    int rowwidth = (rowbound[1]-rowbound[0]+1) + 2*drow;
    int rowlen   = (rowbound[1]-rowbound[0]+1);
    int col_origin = colbound[0]-dcol;
    int row_origin = rowbound[0]-drow;
    if ( colwidth<=0 || rowwidth<=0 ) {
      LARCV_WARNING() << "Bad cropped window size. colwidth=" << colwidth << " rowwidth=" << rowwidth << std::endl;
      return false;
    }

    std::vector<float> crop(colwidth*rowwidth,0);
    std::vector<float> cropmask(colwidth*rowwidth,0);

    // copy the data into the cropped matrix
    // preserve the row-ordered data
    for (int c=0; c<=(int)(colbound[1]-colbound[0]); c++) {
      // for the larcv image, the row-dimension is the contiguous element
      size_t origin_index = meta.index( rowbound[0], colbound[0]+c, __FILE__, __LINE__ );
      const float* source_ptr = img.as_vector().data() + origin_index;
      float* dest_ptr   = crop.data() + (c+dcol)*rowwidth + drow;
      memcpy( dest_ptr, source_ptr, sizeof(float)*rowlen );
    }
        
    // kernal loop, making mask
    // trying to write in a way that is easy to accelerate
    int N = (2*dcol+1)*(2*drow+1); // number of elements in kernel

    // parallelize masking with block kernel
    std::chrono::steady_clock::time_point start_mask = std::chrono::steady_clock::now();    

    for (int ikernel=0; ikernel<N; ikernel++) {

      int dr = ikernel%(2*drow+1) - drow;
      int dc = ikernel/(2*drow+1) - dcol;          

      #pragma omp parallel for
      for (size_t itp=0; itp<pixlist_v.size(); itp++) {

        auto const& col = pixlist_v[itp][0];        
        auto const& row = pixlist_v[itp][1];
        
        // change (col,row) from original image to mask image (with padding)
        // then add kernal shift
        int icol = (int)col-col_origin + dc;
        int irow = (int)row-row_origin + dr;
        if ( icol>=colwidth || irow>=rowwidth || icol<0 || irow<0) {
          std::stringstream ss;
          ss << "OUT OF CROP BOUNDS (" << icol << "," << irow << ") TP[" << itp << "] orig=(" << col << "," << row << ")" << std::endl;
          throw std::runtime_error(ss.str());
        }

        int cropindex = icol*rowwidth+irow;
        float pixval = crop[ cropindex ];

        if ( pixval>thresh ) {
          #pragma omp atomic
          cropmask[ cropindex ] += 1.0;
        }
      }
      #pragma omp barrier      
    }

    pixval_min = 1e9;
    pixval_max = -1e9;
    for (size_t i=0; i<cropmask.size(); i++) {
      if ( cropmask[i]>0 ) {
        if ( crop[i]<pixval_min ) pixval_min = crop[i];
        if ( crop[i]>pixval_max ) pixval_max = crop[i];
      }
    }

    std::chrono::steady_clock::time_point end_mask = std::chrono::steady_clock::now();
    if ( show_timing )
      LARCV_INFO() << "mask and find extrema: " << std::chrono::duration_cast<std::chrono::microseconds>(end_mask - start_mask).count() << std::endl;
    
    return true;    
  }
  
  
}
}
