#include "NeutrinoPixelFilter.h"

#include <iostream>

#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventImage2D.h"
// #include "DataFormat/storage_manager.h"
// #include "DataFormat/mctrack.h"
// #include "DataFormat/mcshower.h"

namespace ublarcvapp {
namespace mctools {

  /**
   * @brief make images with only nu pixels using event data in larcv and larlite io manager
   */
  void NeutrinoPixelFilter::process( larcv::IOManager& iolarcv )
  {
    larcv::EventImage2D* ev_adc =
      (larcv::EventImage2D*)iolarcv.get_data(larcv::kProductImage2D,"wiremc"); // wire data
    // larcv::EventImage2D* ev_instance =
    //   (larcv::EventImage2D*)iolarcv.get_data(larcv::kProductImage2D,"instance"); // geant trackid that made charge in pixel
    // larcv::EventImage2D* ev_ancestor =
    //   (larcv::EventImage2D*)iolarcv.get_data(larcv::kProductImage2D,"ancestor"); // geant ancestor id that made charge in pixel
    larcv::EventImage2D* ev_segment =
      (larcv::EventImage2D*)iolarcv.get_data(larcv::kProductImage2D,"segment"); // geant segment id that made charge in pixel

    // larlite::event_mctrack* ev_mctrack
    //   = ( larlite::event_mctrack*)ioll.get_data(larlite::data::kMCTrack,"mcreco"); // geant track ids
    // larlite::event_mcshower* ev_mcshower
    //   = ( larlite::event_mcshower*)ioll.get_data(larlite::data::kMCShower,"mcreco"); // geant shower ids

    // we extract the track ids, anccestor ids for particles with origin=1 (neutrino)
    // std::set<int> neutrino_ids;
    // for (auto const& track : *ev_mctrack ) {
    //   if ( track.Origin()==1 ) {
    //     neutrino_ids.insert( track.TrackID() );
    //     neutrino_ids.insert( track.MotherTrackID() );
    //     neutrino_ids.insert( track.AncestorTrackID() );
    //   }
    // }
    // for (auto const& shower : *ev_mcshower ) {
    //   if ( shower.Origin()==1 ) {
    //     neutrino_ids.insert( shower.TrackID() );
    //     neutrino_ids.insert( shower.MotherTrackID() );
    //     neutrino_ids.insert( shower.AncestorTrackID() );
    //   }
    // }

    // std::stringstream ss;
    // for ( auto& iid : neutrino_ids )
    //   ss << iid << " ";
    // std::cout << "neutrino-origin track ids: " << ss.str() << std::endl;

    // now we loop over the ancestor id
    std::vector< larcv::Image2D > nuimg_v;
    for ( size_t p=0; p<ev_adc->as_vector().size(); p++) {
      auto const& adc = ev_adc->as_vector().at(p);
      // auto const& instance = ev_instance->as_vector().at(p);
      // auto const& ancestor = ev_ancestor->as_vector().at(p);
      auto const& segment  = ev_segment->as_vector().at(p);
      larcv::Image2D nuimg( adc.meta() );
      nuimg.paint(0.0);

      for (size_t r=0; r<adc.meta().rows(); r++) {
        for (size_t c=0; c<adc.meta().cols(); c++) {
          // int iid = (int)instance.pixel(r,c);
          // int aid = (int)ancestor.pixel(r,c);
          int sid = (int)segment.pixel(r,c);

          // segment ID for the win!
          if ( sid>0 ) {
            nuimg.set_pixel(r,c,adc.pixel(r,c));
          }
          else {
            //std::cout << "pixel not found: iid=" << iid << " aid=" << aid << std::endl;
          }
        }
      }
      nuimg_v.emplace_back( std::move(nuimg) );
    }

    larcv::EventImage2D* evout_nupix =
      (larcv::EventImage2D*)iolarcv.get_data(larcv::kProductImage2D,"nupix");
    evout_nupix->Emplace( std::move(nuimg_v) );
    
  }
  
}
}
