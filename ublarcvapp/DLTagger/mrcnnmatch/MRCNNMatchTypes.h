#ifndef __MRCNN_MATCH_TYPES_H__
#define __MRCNN_MATCH_TYPES_H__

#include <iostream>
#include "larcv/core/DataFormat/ClusterMask.h"

namespace ublarcvapp {
namespace dltagger {

  /**
   * Class used to calculate and store tick and detector-z bounds from larcv::ClusterMask
   */
  class MaskMatchData {
  public:

    MaskMatchData( int plane, int index, const larcv::ClusterMask& mask );
    virtual ~MaskMatchData() {};

    int plane;
    int index; // in clustermask vector
    float tick_min;
    float tick_max;
    float wire_min;
    float wire_max;
    float detz_min;
    float detz_max;
    bool  used;

    // comparison operator needed for sorting
    bool operator < (const MaskMatchData& b) const {
      if ( tick_min<b.tick_min ) return true;
      else if ( tick_min==b.tick_min ) {
        if ( detz_min<b.detz_min ) return true;
      }
      return false;
    };

    // stdout streamer
    friend std::ostream& operator<<(std::ostream &os,const MaskMatchData& m);
    
  };

  class MaskCombo {
  public:
    MaskCombo()
      : pmasks(3,nullptr),
      pdata(3,nullptr),
      indices(3,-1),
      maskdata_indices(3,-1),
      intersection_tick(2,-1),
      union_tick(2,-1),
      intersection_detz(2,-1),
      union_detz(2,-1),
      IOU(0)
        {};
    virtual ~MaskCombo() {};
    
    void addMask( const larcv::ClusterMask& mask, const MaskMatchData& data );
    bool iscompatible( const MaskMatchData& ) const;
    float calc_iou( bool wdetz_correction=false ) const;
    float iou() const { return IOU; };

    std::vector<const larcv::ClusterMask*> pmasks;
    std::vector<const MaskMatchData*> pdata;    
    std::vector<int>   indices; //< clustermask input vector index for each plane
    std::vector<int>   maskdata_indices; //< source maskdata input vector index for each plane
    std::vector<float> intersection_tick;
    std::vector<float> union_tick;
    std::vector<float> intersection_detz;
    std::vector<float> union_detz;
    float IOU;

    // stdout streamer
    friend std::ostream& operator<<(std::ostream &os,const MaskCombo& m);
    
    // comparison operator for sorting: we sort high to low
    bool operator < (const MaskCombo& rhs ) const {
      if ( iou() > rhs.iou() ) return true;
      return false;
    };
  };


  
}
}

#endif
