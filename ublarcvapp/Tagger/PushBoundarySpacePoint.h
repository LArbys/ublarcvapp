#ifndef __PushBoundaryPoints_h__

#include <vector>
#include "larcv/core/DataFormat/Image2D.h"

#include "BoundaryMuonTaggerTypes.h"
#include "BoundarySpacePoint.h"
#include "FoxTrotTrackerAlgoTypes.h"
#include "FoxTrotTrackerAlgoConfig.h"
#include "FoxTrotTrackerAlgo.h"

namespace ublarcvapp {
namespace tagger {
  
    class PushBoundarySpacePoint {
    public:
        PushBoundarySpacePoint();
        virtual ~PushBoundarySpacePoint() {};

    BoundarySpacePoint pushPoint( const BoundarySpacePoint& boundarypoint,
                                  const std::vector<larcv::Image2D>& img_v,
                                  const std::vector<larcv::Image2D>& badch_v );

    bool isPixelWithinImage( const std::vector<larcv::Image2D>& img_v,
                             std::vector<int> imgcoords); // Check to ensure that a pixel is located within the image.


    void clear();

    protected:
        // submethods
        FoxTrack runFoxTrot( const BoundarySpacePoint& sp,
                             const std::vector<larcv::Image2D>& img_v,
                             const std::vector<larcv::Image2D>& badch_v ); //< runs foxtrottrackeralgo to try and go from spacepoint to end of track

        BoundarySpacePoint scanTrackForEndPoint( const BoundarySpacePoint& original,
                                                 const FoxTrack& track,
                                                 const std::vector<larcv::Image2D>& img_v,
                                                 const std::vector<larcv::Image2D>& badch_v ); //< scans along foxtrot track to look for new end point

        BoundarySpacePoint evalEndPoint( const BoundarySpacePoint& sp ); //< evaluate end point quality -- and reclassify based on location

        // algo-variables
        std::vector<FoxTrack>           m_tracklist_v;
        std::vector<BoundarySpacePoint> m_endpoints_v;

        // used algos
        FoxTrotTrackerAlgoConfig m_foxalgo_cfg;
        FoxTrotTrackerAlgo       m_foxalgo;
        Segment3DAlgo            m_segment_algo;

    };
}
}

#endif
