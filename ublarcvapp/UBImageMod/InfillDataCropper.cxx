#ifndef __INFILLDATACROPPER_CXX__
#define __INFILLDATACROPPER_CXX__

#include <sstream>
#include <random>

#include "InfillDataCropper.h"

//larcv
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"
#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/DataFormat/IOManager.h"
// #include "../../core/ROOTUtil/ROOTUtils.h"

// //larlite
// #include "larlite/LArUtil/Geometry.h"
// #include "larlite/LArUtil/LArProperties.h"

//root
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"

namespace ublarcvapp {

  static InfillDataCropperProcessFactory __global_InfillDataCropperProcessFactory__;

  InfillDataCropper::InfillDataCropper(const std::string name)
    : larcv::ProcessBase(name)
  {}

  void InfillDataCropper::configure(const larcv::PSet& cfg)
  {

    _verbosity_             = cfg.get<int>("Verbosity");

    _input_adc_producer     = cfg.get<std::string>("InputADCProducer");
    _input_chstatus_producer     = cfg.get<std::string>("InputChStatusProducer");

    _output_adc_producer     = cfg.get<std::string>("OutputADCProducer");
    _output_adcmasked_producer     = cfg.get<std::string>("OutputADCMaskedProducer");
    _output_labels_producer     = cfg.get<std::string>("OutputLabelsProducer");
    _output_weights_producer     = cfg.get<std::string>("OutputWeightsProducer");

    _output_filename        = cfg.get<std::string>("OutputFilename");

    _max_images             = cfg.get<int>("MaxImages",-1);
    _thresholds_v           = cfg.get< std::vector<float> >("Thresholds",std::vector<float>(3,10.0) );


    // sparsity requirement: prevent cropping over the same regions
    _limit_overlap        = cfg.get<bool>("LimitOverlap",false);
    _max_overlap_fraction = cfg.get<float>("MaxOverlapFraction", 0.5 );

    // output file
    foutIO = new larcv::IOManager( larcv::IOManager::kWRITE );
    foutIO->set_out_file( _output_filename );
    foutIO->initialize();

  }

  void InfillDataCropper::initialize()
  {}

  bool InfillDataCropper::process(larcv::IOManager& mgr)
  {
    // ---------------------------------------------------------------
    // get data

    // input ADC "truth"
    auto ev_in_adc  = (larcv::EventImage2D*)(mgr.get_data(larcv::kProductImage2D, _input_adc_producer));
    if (!ev_in_adc) {
      LARCV_CRITICAL() << "No Input ADC Image2D found with a name: " << _input_adc_producer << std::endl;
      throw larcv::larbys();
    }
    std::vector< larcv::Image2D > img_adc_v = ev_in_adc->Image2DArray();

    // input chstatus
    larcv::EventChStatus* ev_chstatus = nullptr;
    if(!_input_chstatus_producer.empty()) {
      ev_chstatus = (larcv::EventChStatus*)(mgr.get_data(larcv::kProductChStatus,_input_chstatus_producer));
      if(!ev_chstatus) {
        LARCV_CRITICAL() << "ChStatus by " << _input_chstatus_producer << " not found!" << std::endl;
        throw larcv::larbys();
      }
    }

    // ----------------------------------------------------------------

    // Output containers
    larcv::EventImage2D* ev_out_adc  = (larcv::EventImage2D*)foutIO->get_data(larcv::kProductImage2D, _output_adc_producer);
    larcv::EventImage2D* ev_out_labels  = (larcv::EventImage2D*)foutIO->get_data(larcv::kProductImage2D, _output_labels_producer);
    larcv::EventImage2D* ev_out_adcmasked  = (larcv::EventImage2D*)foutIO->get_data(larcv::kProductImage2D,_output_adcmasked_producer);
    larcv::EventImage2D* ev_out_weights  = (larcv::EventImage2D*)foutIO->get_data(larcv::kProductImage2D,_output_weights_producer);

    ev_out_adc->clear();
    ev_out_labels->clear();
    ev_out_adcmasked->clear();
    ev_out_weights->clear();

    // ----------------------------------------------------------------

    // Copies to work with
    // Make a copy to store labels
    // 1 = dead wires
    // 0 = not dead wires
    std::vector<larcv::Image2D> wire_image_v = ev_in_adc->Image2DArray();
    std::vector<larcv::Image2D> image_label_v = {wire_image_v[0].meta() , wire_image_v[1].meta() , wire_image_v[2].meta()};

    // Make a copy to store masked ADC
    auto ev_in_adc_masked  = (larcv::EventImage2D*)(mgr.get_data(larcv::kProductImage2D, _input_adc_producer));
    std::vector< larcv::Image2D > img_adcmasked_v = ev_in_adc_masked->Image2DArray();

    // ----------------------------------------------------------------

    // Run, subrun, event
    int run    = ev_in_adc->run();
    int subrun = ev_in_adc->subrun();
    int event  = ev_in_adc->event();

    int n=0;
    int saved = 0;
    int nboxes = 100;
    int cropwidth = 496; //width of crops
    int cropheight = 512; //width of crops

    image_label_v = ChStatusToLabels(image_label_v,ev_chstatus);
    std::map<int, std::vector<int>> CropStartDict;
    CropStartDict = DeadChannelAnalyzer(image_label_v,cropwidth);
    int nimages = 0;

    // store an image used to keep track of previously cropped pixels
    // std::vector<larcv::Image2D> overlap_img;
    // const std::vector<larcv::ImageMeta>& src_meta = (wire_image_v[0].meta() , wire_image_v[1].meta() , wire_image_v[2].meta());
    // if ( _limit_overlap ) {
    //   for (int plane = 0; plane < 3; plane++ ){
    //     larcv::Image2D overlap( wire_image_v[plane].meta() );
    //     overlap.paint(0.0);
    //     overlap_img.emplace_back( std::move(overlap) );
    //   }
    // }

    while (nimages < nboxes ){
      std::vector<larcv::Image2D> cropped_adc;
      std::vector<larcv::Image2D> cropped_labels;
      std::vector<larcv::Image2D> cropped_adcmasked;
      std::vector<larcv::Image2D> cropped_weights;

      cropped_labels = cropLabelsImage(
                    cropwidth,
                    cropheight,
                    image_label_v,
                    cropped_labels);

      bool passed_labelcheck;
      bool passed_adccheck;

      passed_labelcheck = checkLabelsCrop(cropped_labels);

      if (passed_labelcheck){
        cropped_adc = cropADCImage(
                          cropwidth,
                          cropheight,
                          wire_image_v,
                          cropped_adc,
                          CropStartDict);


        // if limiting overlap, check that we are not overlapping too many pixels
      //   if ( _limit_overlap ) {
      //     for (int plane = 0; plane<3; plane++){
      //       const larcv::ImageMeta& cropped_src_meta = cropped_adc[plane].meta();
      //       float npix_overlap = 0.0;
      //       int crop_r_start = src_meta[plane].row( cropped_src_meta.min_y() );
      //       int crop_c_start = src_meta[plane].col( cropped_src_meta.min_x() );
      //       std::cout << "crop meta: " << cropped_src_meta.min_y() << " , " <<cropped_src_meta.min_x() <<std::endl;
      //       std::cout << "src meta: " << src_meta[plane].min_y() << " , " << src_meta[plane].min_x() <<std::endl;
      //       for (int r=crop_r_start; r<crop_r_start+cropheight; r++) {
      //         for (int c=crop_c_start; c<crop_c_start+cropwidth; c++) {
      //           if ( overlap_img[plane].pixel(r,c)>0 )
      //             npix_overlap+=1.0;
      //         }
      //       }
      //       float frac_overlap = npix_overlap/float(cropped_src_meta.rows()*cropped_src_meta.cols());
      //       if ( frac_overlap>_max_overlap_fraction ) {
      //         LARCV_NORMAL() << "Skipping overlapping image. Frac overlap=" << frac_overlap << "." << std::endl;
      //         continue;
      //       }
      //       std::cout << "Overlap fraction: " << frac_overlap << std::endl;
      //     }
      // }

        cropped_adcmasked = cropADCMaskImage(
                          cropwidth,
                          cropheight,
                          wire_image_v,
                          cropped_labels,
                          cropped_adc,
                          cropped_adcmasked);

        passed_adccheck = checkADCCrop(
                          cropwidth,
                          cropheight,
                          cropped_adc,
                          cropped_labels);

        if (passed_adccheck) {
             nimages += 1;
             // std::cout << "saved Images: " << nimages <<std::endl;
            //make weights image

            std::vector< larcv::Image2D > weights_v;
            for (auto const& img : cropped_labels){
              larcv::Image2D copy = img;
              weights_v.emplace_back(std::move(copy));
            }

            weights_v = cropWeightsImage(
                          cropwidth,
                          cropheight,
                          weights_v,
                          cropped_labels,
                          cropped_adc);
          //
          // if ( _limit_overlap ) {
          //   for (int plane = 0; plane<3; plane++){
          // 	  const larcv::ImageMeta& cropped_src_meta = cropped_adc[plane].meta();
          // 	  int crop_r_start = src_meta[plane].row( cropped_src_meta.min_y() );
          // 	  int crop_c_start = src_meta[plane].col( cropped_src_meta.min_x() );
          // 	  for (int r=crop_r_start; r<crop_r_start+(int)cropped_src_meta.rows(); r++) {
          // 	    for (int c=crop_c_start; c<crop_c_start+(int)cropped_src_meta.cols(); c++) {
          // 	      overlap_img[plane].set_pixel(r,c,1.0);
          // 	    }
          // 	  }
          //   }
        	// }


          ev_out_weights->Emplace( std::move(weights_v) );
          ev_out_labels->Emplace( std::move(cropped_labels) );
          ev_out_adcmasked->Emplace( std::move(cropped_adcmasked) );
          ev_out_adc->Emplace( std::move(cropped_adc) );

          std::cout<<"made crops" <<std::endl;

        	foutIO->set_id( run, subrun, 100*event+nimages );
        	foutIO->save_entry();
          saved++;
          std::cout << "saved Images: " << saved<<std::endl;
          if ( _max_images>0 && saved >= _max_images) break;

        }
      }

      }
    n++;
   //    if ( _max_images>0 && saved >= _max_images) break;
   //
   // }
    return true;
  }

// =============================================================================
  std::vector<larcv::Image2D> InfillDataCropper::ChStatusToLabels(
              std::vector<larcv::Image2D>& image_label_v,
              larcv::EventChStatus* ev_chstatus){

    //function to convert chstatus object to a labels image
    //In order to label, first set all pixels in image to zero

    for(int plane =0; plane<3; plane++)
      {
        image_label_v[plane].paint(0.0);
      }

    //loop that: sets dead channels to one for labels image
    //loop through planes

    for(int plane=0; plane<(int)(image_label_v.size()); plane++) {

      std::vector<size_t> wire_v;
      if (plane ==2){
        std::vector<size_t> wire_v(3456);
        std::iota(wire_v.begin(), wire_v.end(), 0);
      }
      else{
        std::vector<size_t> wire_v(2400);
        std::iota(wire_v.begin(), wire_v.end(), 0);
      }

      // construct wire numbers to mask
      std::set<size_t> wire_s;
      for(auto const& w : wire_v) wire_s.insert(w);

      if(ev_chstatus) {
        auto const& status_v = ev_chstatus->Status(plane).as_vector();

        for(size_t w=0; w<status_v.size(); ++w) {
          if(status_v[w] < 4) wire_s.insert(w);
        }
      }


       //create an image for basic labels
       auto& img_label = image_label_v.at(plane);
       auto const& meta_label = img_label.meta();
       std::vector<float> empty_column(img_label.meta().rows(),1);

       // Loop over wires, set dead channels in labels to one
       for(auto const& ch : wire_s) {
        if(ch < meta_label.min_x() || ch >= meta_label.max_x()) continue;
        auto col = meta_label.col((double)ch);
        //set columns in label image to one
        img_label.copy(1,col,empty_column);
       }
    }
    return image_label_v;
  }

  std::vector<larcv::Image2D> InfillDataCropper::ChStatusToLabelsX(larcv::EventChStatus& ev_chstatus) {
    std::vector<larcv::Image2D> label_v;
    InfillDataCropper::ChStatusToLabels( label_v, &ev_chstatus );
    return label_v;
  }
								  


// =============================================================================
  std::map<int, std::vector<int>> InfillDataCropper::DeadChannelAnalyzer(
        const std::vector<larcv::Image2D>& image_label_v,
        const int cropwidth){
    // function to characterize dead channels
    // returns dictionary
    // 0,1,2 are the possible starts of crops in each plane
    std::map<int, std::vector<int>> CropStartDict;
    for (int plane=0; plane<3; plane++){
      size_t cols;
      // hard code num of cols for uboone
      if (plane == 2) { cols = 3456;}
      else {cols = 2400;}
      int numlive = 0;
      int numpossiblecrops = 0;
      std::vector<int> startchannels;
      for (unsigned int col=0; col<cols; col++) {
        // if live channel
        if (image_label_v[plane].pixel(0,col) == 0) {
          numlive += 1;
          if (numlive >= cropwidth){
            startchannels.push_back(col - cropwidth);
            numpossiblecrops+=1;}
        }
        // if dead channel
        else {
          numlive = 0;
          }
      }
      CropStartDict[plane] = startchannels;
    }
    return CropStartDict;
  }


// =============================================================================
std::vector<larcv::Image2D> InfillDataCropper::cropLabelsImage(
    const int cropwidth,
    const int cropheight,
    const std::vector<larcv::Image2D>& img_v,
    std::vector<larcv::Image2D>& cropped_labels){

  // make vectors of possible indices and pick a random entry for
  // do loop until we've saved enough crops
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> distr_t(0, 1008-cropheight);
    double t2 = (double) distr_t(eng);
    std::uniform_int_distribution<> distr_u(0, 2400-cropwidth);
    double u1 = (double) distr_u(eng);
    std::uniform_int_distribution<> distr_v(0, 2400-cropwidth);
    double v1 = (double) distr_v(eng);
    std::uniform_int_distribution<> distr_y(0, 3456-cropwidth);
    double y1 = (double) distr_y(eng);

    std::vector<larcv::ImageMeta> meta_labels;
    meta_labels = defineBoundingBoxFromCropCoords(
        img_v, cropwidth, cropheight,
        t2,u1,v1,y1);

    // crop the image
    cropped_labels = cropImagefromBox(
                    meta_labels,
                    img_v,
                    cropped_labels);

    return cropped_labels;
}

// =============================================================================

 bool InfillDataCropper::checkLabelsCrop(std::vector<larcv::Image2D>& cropped_labels){
      // function to make sure that there are dead regions in the labels crop

      float dead = 0.;
      bool saveimg = true;

      for (int plane=0; plane<3; plane++){
        dead = 0;
        const larcv::Image2D& labels = cropped_labels[plane];
        for (int col=0; col<(int)labels.meta().cols(); col++) {
          if ( labels.pixel(0,col) > 0 ) {dead+=1.0;}
        }
        if ( dead < 3 ) saveimg = false;
        if ( !saveimg ){return false;}
      }

      return true;
      }

// =============================================================================


std::vector<larcv::Image2D> InfillDataCropper::cropADCImage(
    const int cropwidth,
    const int cropheight,
    const std::vector<larcv::Image2D>& img_v,
    std::vector<larcv::Image2D>& cropped_adc,
    std::map<int, std::vector<int>> startdict){

    //  start by randomly accessing values in t - same as labels
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> distr_t(0, 1008-cropheight);
    double t2 = (double) distr_t(eng);

    //  randomly access from startdict
    std::uniform_int_distribution<> distr_u(0, startdict[0].size()-1);
    double u1 = (double) startdict[0][distr_u(eng)];

    std::uniform_int_distribution<> distr_v(0, startdict[1].size()-1);
    double v1 = (double) startdict[1][distr_v(eng)];

    std::uniform_int_distribution<> distr_y(0, startdict[2].size()-1);
    double y1 = (double) startdict[2][distr_y(eng)];

    std::vector<larcv::ImageMeta> meta_ADC;
    meta_ADC = defineBoundingBoxFromCropCoords(
        img_v, cropwidth, cropheight,
        t2,u1,v1,y1);

    cropped_adc = cropImagefromBox(
                    meta_ADC,
                    img_v,
                    cropped_adc);


    return cropped_adc;
    }

// =============================================================================
std::vector<larcv::Image2D> InfillDataCropper::cropADCMaskImage(
    const int cropwidth,
    const int cropheight,
    const std::vector<larcv::Image2D>& img_v,
    const std::vector<larcv::Image2D>& img_labels_v,
    const std::vector<larcv::Image2D>& img_adc_v,
    std::vector<larcv::Image2D>& cropped_adcmasked){

    // function to mask good adc regions with cropped masks

    for(int plane=0; plane< 3; plane++) {
      const larcv::Image2D& adc = img_adc_v[plane];
      const larcv::Image2D& labels = img_labels_v[plane];

      // set new image meta based on meta of adc
      larcv::ImageMeta meta_ADC = img_adc_v[plane].meta();
      larcv::Image2D cropimg = img_v[plane].crop( meta_ADC );
      cropped_adcmasked.emplace_back( std::move(cropimg) );

      for (int col=0; col<cropwidth; col++){
          for (int row=0; row<cropheight; row++) {
            if (labels.pixel(row,col) == 1){
              cropped_adcmasked[plane].set_pixel(row,col,0);
            }
            else{
              float pixelval =adc.pixel(row,col);
              cropped_adcmasked[plane].set_pixel(row,col,pixelval);
            }
          }
       }
    }


    return cropped_adcmasked;
    }

// =============================================================================

bool InfillDataCropper::checkADCCrop(
    const int cropwidth,
    const int cropheight,
    std::vector<larcv::Image2D>& cropped_adc,
    std::vector<larcv::Image2D>& cropped_labels){

     float occupied = 0.;
     bool saveimg = true;
     int minpixfrac=10;
     if ( minpixfrac>=0 ) {
       for (int plane=0; plane<3; plane++){
         occupied=0;
         const larcv::Image2D& adc = cropped_adc[plane];
         const larcv::Image2D& labels = cropped_labels[plane];
         for (int col=0; col<cropwidth; col++) {
           for (int row=0; row<cropheight; row++) {
            float labelvalue = labels.pixel(row,col);
            float adcvalue = adc.pixel(row,col);
            if ( (adcvalue*labelvalue) > 0 ){
              occupied+=1.0;}
            }
          }
         if ( occupied <= minpixfrac ) saveimg = false;
         if ( !saveimg )return false;
       }
     }
     return true;
     }

// =============================================================================
std::vector<larcv::Image2D> InfillDataCropper::cropWeightsImage(
              const int cropwidth,
              const int cropheight,
              std::vector<larcv::Image2D>&weights_v,
              const std::vector<larcv::Image2D>& cropped_labels,
              const std::vector<larcv::Image2D>& cropped_adc){

  //first do a loop to get totals
  for (int plane=0; plane<3; plane++){
    //four different classes
    // 0: adc = 0
    // 1: adc>= 10
    // 2: 0 <= adc <10
    // 3: dead channel , no charge

    float total0=0.0;
    float total1=0.0;
    float total2=0.0;
    float total3=0.0;
    size_t cols = cropwidth;
    size_t rows = cropheight;

    for (unsigned int col=0; col<cols; col++) {
      for (unsigned int row=0; row<rows; row++) {
          float pix_id = cropped_adc[plane].pixel(row,col);
          if ( cropped_labels[plane].pixel(row,col) == 0 ) total0 += 1.0;
          // class 1
          else if ( pix_id >= 10 ) total1 += 1.0;
          //class 2
          else if ( pix_id > 0) total2 += 1.0;
          //class 3
          else {total3 +=1.0;}
       }
     }

    //loop again to set values in weights
    for (unsigned int col=0; col<cols; col++) {
      for (unsigned int row=0; row<rows; row++) {
          float pix_id = cropped_adc[plane].pixel(row,col);
          if ( cropped_labels[plane].pixel(row,col) == 0 ) weights_v[plane].set_pixel(row,col, 1.0/total0);
          // class 1
          else if ( pix_id >= 10 ) weights_v[plane].set_pixel(row,col,100.0/total1);
          //class 2
          else if ( pix_id > 0) weights_v[plane].set_pixel(row,col, 50.0/total2);
          //class 3
          else {weights_v[plane].set_pixel(row,col, 10.0/total3);}
       }
     }

   }

  return weights_v;
}


// =============================================================================
  std::vector<larcv::ImageMeta> InfillDataCropper::defineBoundingBoxFromCropCoords(
                         const std::vector<larcv::Image2D>& img_v,
									       const int box_pixel_width, const int box_pixel_height,
									       const double t2, const double u1,
                         const double v1, const double y1) {

    // takes pre-defined image bounds on all 3 planes (given in min/max row/col)
    // note, box_pixel_width and box_pixel_height are meant to be the same
    // variable as member variables _box_pixel_width and _box_pixel_height.
    // we pass them here in order to make this function static, so it can be used as a stand-alone function.

    // input
    // ------
    // img_v: source ADC images from which we are cropping
    // box_pixel_width: corresponds to _box_pixel_width
    // box_pixel_height: corresponds to _box_pixel_height
    // (x)1, (x)2, where x=t,u,v,y
    // row bounds (t) and col bounds (u,v,y) max-exclusive [x1,x2)
    //
    // output
    // -------
    // (return) vector of bounding boxes defined for (u,v,y) plane

    std::vector< larcv::ImageMeta > meta_vec; // we create one for each plane
    const larcv::ImageMeta& meta = img_v.front().meta();

    double maxt = meta.pos_y(t2);

    // define tick and row bounds
    int nrows = box_pixel_height*6;
    int ncols = box_pixel_width;

    larcv::PlaneID_t plane = 0;
    larcv::ImageMeta metacropu( box_pixel_width, box_pixel_height*6, nrows, ncols, u1, maxt, plane);

    plane = 1;
    larcv::ImageMeta metacropv( box_pixel_width, box_pixel_height*6, nrows, ncols, v1, maxt, plane);

    plane = 2;
    larcv::ImageMeta metacropy( box_pixel_width, box_pixel_height*6, nrows, ncols, y1, maxt, plane);

    meta_vec.emplace_back( std::move(metacropu) );
    meta_vec.emplace_back( std::move(metacropv) );
    meta_vec.emplace_back( std::move(metacropy) );

    return meta_vec;

  }

// =============================================================================
  std::vector<larcv::Image2D> InfillDataCropper::cropImagefromBox(
           std::vector<larcv::ImageMeta> meta_vec,
					 const std::vector<larcv::Image2D>& img_v,
					 std::vector<larcv::Image2D>& output_imgs) {
    // inputs
    // ------
    // bbox_v, vector of bounding boxes for (u,v,y)
    // img_v, source adc images
    // nimages, number of the bbox to use
    // outputs
    // --------
    // output_imgs, cropped output image2d instances filled into eventimage2d container


    // get bounding boxes
    const larcv::ImageMeta metacrop_u = meta_vec[0];
    const larcv::ImageMeta metacrop_v = meta_vec[1];
    const larcv::ImageMeta metacrop_y = meta_vec[2];


    //crop
    larcv::Image2D crop_up = img_v[0].crop( metacrop_u );
    output_imgs.emplace_back( std::move(crop_up) );

    larcv::Image2D crop_vp = img_v[1].crop( metacrop_v );
    output_imgs.emplace_back( std::move(crop_vp) );
    //std::cout<< "nimages " <<nimages <<std::endl;
    larcv::Image2D crop_yp = img_v[2].crop( metacrop_y );
    output_imgs.emplace_back( std::move(crop_yp) );
    
    return output_imgs;
  }

// =============================================================================
  void InfillDataCropper::finalize()
  {
    foutIO->finalize();
  }

}
#endif
