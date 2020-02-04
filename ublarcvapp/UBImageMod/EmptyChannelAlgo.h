#ifndef __EMPTY_CHANNEL_ALGO__
#define __EMPTY_CHANNEL_ALGO__

#include <vector>

// larlite
#include "larlite/core/Base/DataFormatConstants.h"
#include "larlite/core/DataFormat/chstatus.h" // needed because of duplicate name!

// larcv
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"


namespace ublarcvapp {

  class EmptyChannelAlgo {
  public:
    
    EmptyChannelAlgo() {};
    virtual ~EmptyChannelAlgo() {};

    larcv::Image2D labelEmptyChannels( float threshold, const larcv::Image2D& tpcimg, const float max_value=-1.0 );

    std::vector<int> findEmptyChannels( float threshold, const larcv::Image2D& tpcimg, const float max_value=-1.0 );

    std::vector<larcv::Image2D> makeBadChImage( int minstatus, int nplanes, int start_tick, int nticks, int nchannels, 
						int time_downsample_factor, int wire_downsample_factor,
						const larlite::event_chstatus& ev_status );
    
    std::vector<larcv::Image2D> makeBadChImage( int minstatus, int nplanes, int start_tick, int nticks, int nchannels, 
						int time_downsample_factor, int wire_downsample_factor,
						const larcv::EventChStatus& ev_status );
    
    std::vector<larcv::Image2D> findMissingBadChs( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badchimgs_v,
      const float empty_ch_threshold, const int max_empty_gap, const float max_vaue=-1.0 );

    std::vector<larcv::Image2D> makeGapChannelImage( const std::vector<larcv::Image2D>& img_v,
                                                     const larcv::EventChStatus& ev_status, int minstatus,
                                                     int nplanes, int start_tick, int nticks, int nchannels, 
                                                     int time_downsample_factor, int wire_downsample_factor,                                                
                                                     const float empty_ch_threshold,
                                                     const int max_empty_gap,
                                                     const float empty_ch_max );
    
    
  };
  
}

#endif

