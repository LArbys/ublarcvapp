#ifndef __FOX_TROT_TRACKER_ALGO_H__
#define __FOX_TROT_TRACKER_ALGO_H__

/* ==================================================
 *
 * FoxTrotTrackerAlgo
 *
 * Tracker using box steps.
 *
 ==================================================*/

#include <vector>

//larcv
#include "larcv/core/DataFormat/Image2D.h"

// ublarcvapp
#include "TaggerCROITypes.h"
#include "BMTrackCluster3D.h"

#include "RadialSegmentSearch.h"
#include "FoxTrotTrackerAlgoConfig.h"
#include "FoxTrotTrackerAlgoTypes.h"
#include "FoxTrotLead.h"

namespace ublarcvapp {
namespace tagger {

  class FoxTrotLeadStraight : public FoxTrotLead {
    // This is default FoxTrotLead instance
  public:
    FoxTrotLeadStraight() { min_dcos = 0.0; };
    virtual ~FoxTrotLeadStraight() {};

    void setMinCos( float mc ) { min_dcos = mc; };

  protected:

    bool _chooseBestSegment_( const FoxStep& current, std::vector<Segment3D_t>& seg3d_v, const FoxTrack& past,
			      const FoxTrotTrackerAlgoConfig& config, const std::vector<larcv::Image2D>& tagged_v, int& best_seg_idx );
    
    float min_dcos;
    
  };
  
  class FoxTrotTrackerAlgo {

    FoxTrotTrackerAlgo() {};

  public:

    FoxTrotTrackerAlgo( const FoxTrotTrackerAlgoConfig& config );
    virtual ~FoxTrotTrackerAlgo();

    FoxTrack followTrack( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v,
			  const BoundarySpacePoint& start );

    FoxTrack followTrack( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, const std::vector<larcv::Image2D>& tagged_v,
			  const BoundarySpacePoint& start, const std::vector<float>& startdir );

    void setUserLead( FoxTrotLead* userlead ) { m_user_lead = userlead; }; // is not responsible for this pointer

  protected:

    FoxTrotTrackerAlgoConfig m_config;

    RadialSegmentSearch m_radialalgo;

    FoxStep getNextStep( const FoxStep& current, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
			 const std::vector<larcv::Image2D>& stepped_v, const FoxTrack& path );
    
    void tagStep( const FoxStep& current, const FoxStep& last, std::vector<larcv::Image2D>& stepped_v, float maxstepsize );    

    FoxTrotLead* m_user_lead;
    FoxTrotLeadStraight m_default_lead;
    
  };

}
}

#endif
