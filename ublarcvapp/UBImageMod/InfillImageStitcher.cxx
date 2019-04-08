#include "InfillImageStitcher.h"

namespace ublarcvapp {

    void InfillImageStitcher::PixelScaling(larcv::Image2D& inputimg,
                    float scalefactor){

        // function to scale pixel values for mcc8/mcc9 conversion
        larcv::ImageMeta inputimgmeta=inputimg.meta();
        double x_min = inputimgmeta.min_x();
        double x_max = inputimgmeta.max_x();
        double y_min = inputimgmeta.min_y();
        double y_max = inputimgmeta.max_y();

        size_t row_min = inputimgmeta.row(y_min);
        size_t col_min = inputimgmeta.col(x_min);

        size_t nrows = (y_max - y_min) / inputimgmeta.pixel_height();
        size_t ncols = (x_max - x_min) / inputimgmeta.pixel_width();
        for(size_t col_index=0; col_index < ncols; ++col_index) {
          for(size_t row_index=0; row_index < nrows; ++row_index) {
            //change each pixel by scalefactor
            double original = inputimg.pixel(row_index+row_min, col_index+col_min);
            double newvalue = original*scalefactor;
            inputimg.set_pixel(row_index+row_min, col_index+col_min,newvalue);
          }
        }//end of loops

    }//end of function

    void InfillImageStitcher::Croploop(larcv::ImageMeta& output_meta,
                                        larcv::Image2D& outimg,
                                        larcv::Image2D& outputimg,
                                        larcv::Image2D& overlapcountimg) {
      // function that loops through the crops to stitch them together
      // also creates the output overlay image
      //Inputs:
        //p: int specifying current plane
        //outputmeta: meta of the final image
        // outimg_v: the output crops from the network
        // outputimg: the final image
        // overlapcountimg: an image to keep track of overlaps for averaging
                            // of the final pixel values


      //loop through all crops to stitch onto overlay image
      larcv::ImageMeta crop_meta=outimg.meta();
      double x_min = std::max(output_meta.min_x(),crop_meta.min_x());
      double x_max = std::min(output_meta.max_x(),crop_meta.max_x());
      double y_min = std::max(output_meta.min_y(),crop_meta.min_y());
      double y_max = std::min(output_meta.max_y(),crop_meta.max_y());

      size_t row_min1 = output_meta.row(y_min);
      size_t col_min1 = output_meta.col(x_min);

      size_t row_min2 = crop_meta.row(y_min);
      size_t col_min2 = crop_meta.col(x_min);

      size_t nrows = (y_max - y_min) / output_meta.pixel_height();
      size_t ncols = (x_max - x_min) / output_meta.pixel_width();

      // first pixel loop
      for(size_t col_index=0; col_index < ncols; ++col_index) {
        for(size_t row_index=0; row_index < nrows; ++row_index) {
          double newvalue = outimg.pixel(row_index+row_min2, col_index+col_min2);
          double original = outputimg.pixel(row_index+row_min1, col_index+col_min1);
          double overlapcount = overlapcountimg.pixel(row_index+row_min1, col_index+col_min1);

          outputimg.set_pixel(row_index+row_min1, col_index+col_min1, newvalue+original);
          overlapcountimg.set_pixel(row_index+row_min1, col_index+col_min1, overlapcount+1.0);
        }
      }

    }//end of function

    void InfillImageStitcher::Overlayloop(int p,larcv::ImageMeta& output_meta,
                                        larcv::Image2D& outputimg,
                                        larcv::Image2D& overlapcountimg,
                                        const std::vector<larcv::Image2D>& wholeview_v,
                                        larcv::EventChStatus& ev_chstatus){

      // loop through pixels of whole img to take average and overlay
      double x_min = output_meta.min_x();
      double x_max = output_meta.max_x();
      double y_min = output_meta.min_y();
      double y_max = output_meta.max_y();

      size_t row_min = output_meta.row(y_min);
      size_t col_min = output_meta.col(x_min);

      size_t nrows = (y_max - y_min) / output_meta.pixel_height();
      size_t ncols = (x_max - x_min) / output_meta.pixel_width();

      for(size_t col_index=0; col_index < ncols; ++col_index) {
        for(size_t row_index=0; row_index < nrows; ++row_index) {
          double original = outputimg.pixel(row_index+row_min, col_index+col_min);
          double overlapcount = overlapcountimg.pixel(row_index+row_min, col_index+col_min);
          double truevalue= wholeview_v.at(p).pixel(row_index+row_min, col_index+col_min);
          if (overlapcount > 0){
            outputimg.set_pixel(row_index+row_min, col_index+col_min, original/overlapcount);
          }
          if (p != 2){
            if (col_index<2400){
                //if not dead, set as original value
                if (ev_chstatus.Status(p).as_vector()[col_index]==4){
                    outputimg.set_pixel(row_index+row_min, col_index+col_min,truevalue);}
            }
          }
          else{
            if (ev_chstatus.Status(p).as_vector()[col_index]==4){
                outputimg.set_pixel(row_index+row_min, col_index+col_min,truevalue);}
          }
        }
      }//end of pixel loop

    }//end of funtion

}//end of ublarcvapp namespace
