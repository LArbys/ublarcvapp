#include "TrackImageMask.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"
#include "larcv/core/DataFormat/ImageMeta.h"

namespace ublarcvapp {
namespace ubimagemod {

  void TrackImageMask::maskTrack( const larlite::track& track,
                                  const larcv::Image2D& img,
                                  larcv::Image2D& mask,
                                  const float threshold,
                                  const int dcol,
                                  const int drow,
                                  const float maxstepsize )
  {

    int npts = track.NumberTrajectoryPoints();    
    LARCV_DEBUG() << "track length=" << npts << std::endl;
    
    if ( npts<=1 ) {
      LARCV_WARNING() << "No mask generated for track with only 1 point" << std::endl;
      // no mask can be generated in this case
      return;
    }

    const float driftv = larutil::LArProperties::GetME()->DriftVelocity();
    const float usec_per_tick = 0.5;
    auto const& meta = img.meta();
    int plane = meta.plane();
    if ( plane<0 || plane>=(int)larutil::Geometry::GetME()->Nplanes() ) {
      LARCV_WARNING() << "Invalid plane: " << plane << std::endl;
      return;
    }

    // microboone only
    const std::vector<Double_t>& firstwireproj = larutil::Geometry::GetME()->GetFirstWireProj(); 
    std::vector<double> orthovect = { 0,
                                      larutil::Geometry::GetME()->GetOrthVectorsY().at(plane),
                                      larutil::Geometry::GetME()->GetOrthVectorsZ().at(plane) };
    
    int nrows = meta.rows();
    int ncols = meta.cols();

    // first we make a list of pixels covered by the track
    // we also get the min and max col bounds
    std::set< std::pair<int,int> > pixel_list;
    int min_col = (int)nrows-1;
    int min_row = (int)ncols-1;
    int max_col = 0;
    int max_row = 0;
         
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
        
        // from larlite Geometry::WireCoordinate(...)
        float wirecoord = pos[1]*orthovect[1] + pos[2]*orthovect[2] - firstwireproj.at(plane);
        float colcoord = (wirecoord-meta.min_x())/meta.pixel_width();
        
        int row = (int)rowcoord;
        int col = (int)colcoord;
        if ( row<0 || row>=nrows ) continue;
        if ( col<0 || col>=ncols ) continue;

        pixel_list.insert( std::pair<int,int>(col,row) );
        if ( min_row>row ) min_row = row;
        if ( max_row<row ) max_row = row;
        if ( min_col>col ) min_col = col;
        if ( max_col<col ) max_col = col;
      }//end of pixel steps
    }//end of trajectory point list

    // convert set into list
    std::vector< std::vector<int> > pixel_v;
    pixel_v.reserve( pixel_list.size() );
    for ( auto& it : pixel_list ) {
      pixel_v.push_back( std::vector<int>{ it.first, it.second} );
    }
    LARCV_DEBUG() << "Number of pixels: " << pixel_v.size() << std::endl;
    
    // ok now for the masking

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
      return;
    }

    // make data array
    std::vector<float> crop(colwidth*rowwidth,0);
    std::vector<float> cropmask(colwidth*rowwidth,0);

    // copy the data into the cropped matrix
    // preserve the row-ordered data
    for (int c=0; c<=(int)(max_col-min_col); c++) {
      // for the larcv image, the row-dimension is the contiguous element
      size_t origin_index = img.meta().index( min_row, min_col+c, __FILE__, __LINE__ );
      const float* source_ptr = img.as_vector().data() + origin_index;
      float* dest_ptr   = crop.data() + (c+dcol)*rowwidth + drow;
      memcpy( dest_ptr, source_ptr, sizeof(float)*rowlen );
    }

    // kernal loop, making mask
    // trying to write in a way that is easy to accelerate
    int N = (2*dcol+1)*(2*drow+1); // number of elements in kernel

    // parallelize masking with block kernel
    int nmasked = 0;
    
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
    
  }
  
}
}
