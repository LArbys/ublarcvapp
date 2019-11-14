#ifndef __DLInteraction_h__
#define __DLInteraction_h__

#include <vector>

// larlite
#include "core/DataFormat/track.h"
#include "core/DataFormat/shower.h"

// larcv
#include "larcv/core/DataFormat/ROI.h"
#include "larcv/core/DataFormat/Pixel2DCluster.h"

namespace ublarcvapp {
namespace ubdllee {

  class DLInteraction {

  public:
    
    DLInteraction() { clear(); };
    virtual ~DLInteraction() {};

    // RUN, SUBRUN, EVENT
    int run;
    int subrun;
    int event;
    
    // CROI
    int croi_filled;
    std::vector< larcv::ROI >      croi_v; // merged

    // thrumu cosmic tagger
    int tagger_filled;    
    int nthrumu_tracks;
    std::vector< larlite::track > thrumu_tracks_v;
    std::vector< std::vector< larcv::Pixel2DCluster > > thrumu_pixels_v; // (wire,tick) in full image
    
    // vertex
    int vertex_filled;
    int vertex_index;
    std::vector<float> vertex_pos;
    std::vector<int>   vertex_wires;
    int                vertex_tick;
    int                vertex_nprongs;
    std::vector< std::vector<larcv::Pixel2DCluster> > vertex_prong_contour_v;
    std::vector< std::vector<larcv::ImageMeta> >      vertex_prong_meta_v;

    // Track Reco
    int tracker_filled;
    int ntracks;
    std::vector< larlite::track > track_v;
    std::vector< std::vector< larcv::Pixel2DCluster > > track_pixels_v; // (wire,tick) in full image
    std::vector< int >  track_goodflag_v;
    // -- track variables --

    // shower reco
    int shower_filled;
    int nshowers;
    std::vector< larlite::shower > shower_v;
    // -- shower variables --

    void clear();
    
  };

}
}

#endif
