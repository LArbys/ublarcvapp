#include "DLSSNetOverlayCheater.h"

#include <sstream>
#include "larcv/core/DataFormat/EventImage2D.h"

namespace ublarcvapp {
namespace dltagger {

  static DLSSNetOverlayCheaterFactory __global_DLSSNetOverlayCheaterFactory__;
  
  void DLSSNetOverlayCheater::configure( const larcv::PSet& pset ) {
    _input_adc_name         = "wire";
    _input_ssnet_input_stem = "uburn";
    _input_segment_name     = "segment";
    _output_stem_treename   = "cheaterssnet";
    _mask_cosmic_pixels     = false;
    _apply_threshold        = true;
    _adc_threshold_v.resize(3,10.0);
  }

  void DLSSNetOverlayCheater::initialize() {
  }

  
  bool DLSSNetOverlayCheater::process( larcv::IOManager& mgr ) {

    larcv::EventImage2D* ev_uburn[3] = { nullptr };
    for ( int p=0; p<3; p++ ) {
      std::stringstream treename;
      treename << _input_ssnet_input_stem << "_plane" << p;      
      ev_uburn[p] = (larcv::EventImage2D*)mgr.get_data(larcv::kProductImage2D,treename.str());
    }

    larcv::EventImage2D* ev_adc     = (larcv::EventImage2D*)mgr.get_data(larcv::kProductImage2D,_input_adc_name);
    larcv::EventImage2D* ev_segment = (larcv::EventImage2D*)mgr.get_data(larcv::kProductImage2D,_input_segment_name);

    
    for ( int p=0; p<3; p++ ) {
      std::vector<larcv::Image2D> ssnetout_v;

      auto const& segment = ev_segment->Image2DArray()[p];
      auto const& adc     = ev_adc->Image2DArray()[p];
        
      // now modify image with truth
      for (int i=0; i<2; i++) {
        // i=0: shower
        // i=1: track
        auto const& origssnet = ev_uburn[p]->Image2DArray()[i];
        larcv::Image2D ssnetout( origssnet );
        
        for (int r=0; r<(int)origssnet.meta().rows(); r++) {
          float tick = origssnet.meta().pos_y(r);
          
          for ( int c=0; c<(int)origssnet.meta().cols(); c++ ) {
            
            // translate ssnet image pixel location to whole-view segment image
            float wire = origssnet.meta().pos_x(c);
            if ( !segment.meta().contains( wire, tick ) ) continue;
            
            int xrow = segment.meta().row(tick);
            int xcol = segment.meta().col(wire);
          
            int segid = segment.pixel(xrow,xcol);

            if ( _apply_threshold ) {            
              float pixadc = adc.pixel(xrow,xcol);
              if ( pixadc < _adc_threshold_v.at(p) ) continue;
            }

            if ( segid>=(int)larcv::kROIEminus && segid<=(int)larcv::kROIPizero ) {
              if ( i==0 ) ssnetout.set_pixel(r,c,1.0); // shower
              else ssnetout.set_pixel(r,c,0.0);        // track
            }
            else if ( segid>=(int)larcv::kROIMuminus && segid<=(int)larcv::kROIProton ) {
              if (i==0) ssnetout.set_pixel(r,c,0.0); // shower
              else ssnetout.set_pixel(r,c,1.0);      // track
            }
            else {
              if ( _mask_cosmic_pixels && segid==0 ) {
                // mask pixels with no mc truth
                ssnetout.set_pixel(r,c,0.0);
              }
              // otherwise, we leave it alone
            }
          }
        }

        ssnetout_v.emplace_back( std::move(ssnetout) );
      }//end of ssnet image loop {shower,track}
      
      std::stringstream outname;
      outname << _output_stem_treename << "_plane" << p;

      larcv::EventImage2D* ev_out = (larcv::EventImage2D*)mgr.get_data(larcv::kProductImage2D, outname.str());
      ev_out->Emplace( std::move(ssnetout_v) );
      
    }//end of plane loop
    
    return true;
  }

  
  void DLSSNetOverlayCheater::finalize() {

  }


  
  
}
}
