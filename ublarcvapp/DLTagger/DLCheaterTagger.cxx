#include "DLCheaterTagger.h"

#include "larcv/core/DataFormat/EventImage2D.h"

namespace ublarcvapp {
namespace dltagger {

  void DLCheaterTagger::configure( const larcv::PSet& pset ) {
    // run on MC input
    _input_adc_producer      = "wiremc";
    _input_chstatus_producer = "wiremc";
    _input_instance_image_producer = "instance";

    _output_tagger_image_producer  = "cheatertagger";
    _output_tagger_pixcluster_producer = "cheatertagger";    
  }

  void DLCheaterTagger::initialize() {
    
  }

  /**
   * process one event.
   *
   * make a cosmic tag image by simply checking which pixels are from MC.
   * this, of course, only works for overlay MC samples.
   * this is the concrete implementation of the abstract interface method from larcv::ProcessBase
   * 
   * @param[inout] mgr The IOManager from which we get data and to which we store our output.
   *
   */  
  bool DLCheaterTagger::process( larcv::IOManager& mgr ) {

    // the event container for wholeview ADC images
    auto const ev_adc = (larcv::EventImage2D*)mgr.get_data( larcv::kProductImage2D, _input_adc_producer );

    // the event container for instance ID images 
    auto const ev_instance = (larcv::EventImage2D*)mgr.get_data( larcv::kProductImage2D, _input_instance_image_producer );


    // make container of output images
    std::vector<larcv::Image2D> out_v;
    for ( auto const& img : ev_adc->Image2DArray() ) {
      larcv::Image2D out( img.meta() );
      out.paint(0.0);

      auto const& instance = ev_instance->Image2DArray().at((int)img.meta().plane());

      for ( size_t c=0; c<img.meta().cols(); c++ ) {
        for ( size_t r=0; r<img.meta().rows(); r++ ) {
          if ( instance.pixel(r,c)>0 ) {
            out.set_pixel(r,c, img.pixel(r,c));
          }
        }
      }

      out_v.emplace_back( std::move(out) );
    }

    auto ev_outtag = (larcv::EventImage2D*)mgr.get_data( larcv::kProductImage2D, _output_tagger_image_producer );
    ev_outtag->Emplace( std::move(out_v) );

    return true;
    
  }

  void DLCheaterTagger::finalize() {

  }
  
}
}
