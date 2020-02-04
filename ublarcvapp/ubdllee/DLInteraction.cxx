#include "DLInteraction.h"

namespace ublarcvapp {
namespace ubdllee {

  void DLInteraction::clear() {
    run = -1;
    subrun = -1;
    event = -1;
    vertex_index = -1;

    croi_filled   = 0;
    tagger_filled = 0;
    vertex_filled = 0;
    shower_filled = 0;

    vertex_nprongs = 0;
    nthrumu_tracks = 0;
    nshowers  = 0;
    ntracks   = 0;

    croi_v.clear();
    track_v.clear();
    shower_v.clear();
    vertex_prong_contour_v.clear();
    track_pixels_v.clear();
    track_goodflag_v.clear();
    thrumu_pixels_v.clear();
    thrumu_tracks_v.clear();
    
  }
  
}
}
