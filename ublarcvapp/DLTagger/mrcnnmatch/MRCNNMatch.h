#ifndef __MRCNN_MATCH_H__
#define __MRCNN_MATCH_H__

#include <vector>

#include "larcv/core/DataFormat/ClusterMask.h"

#include "MRCNNMatchTypes.h"
#include "CropMaskCombo.h"
#include "FeaturesMaskCombo.h"
#include "Gen3DEndpoints.h"
#include "AStarMaskCombo.h"

namespace ublarcvapp {
namespace dltagger {

  class MRCNNMatch {
  public:
    MRCNNMatch() {};
    virtual ~MRCNNMatch() {};

    void matchMasksAcrossPlanes( const std::vector< std::vector<larcv::ClusterMask> >& clustermask_vv,
                                 const std::vector<larcv::Image2D>& wholeview_v,
                                 const larcv::EventChStatus& ev_chstatus,                                 
                                 std::vector< std::vector<int> >& match_indices );

    std::vector< MaskCombo >         m_combo_3plane_v;   // matches across planes
    std::vector< CropMaskCombo >     m_combo_crops_v;    // cropped images from matches
    std::vector< FeaturesMaskCombo > m_combo_features_v; // features extracted from the images
    std::vector< Gen3DEndpoints >    m_combo_endpt3d_v;  // 3D endpoints created for mask-combo
    std::vector< AStarMaskCombo >    m_combo_astar_v;    // AStar path


  protected:
    
    void run3PlanePass( const std::vector<std::vector<larcv::ClusterMask>>& clustermask_vv,
                        const std::vector<larcv::Image2D>& wholeview_v,
                        const std::vector<larcv::Image2D>& badch_v,
                        std::vector< std::vector<MaskMatchData> >& matchdata_vv,
                        std::vector< std::vector<int> >& match_indices );

    void run2PlanePass( const std::vector<std::vector<larcv::ClusterMask>>& clustermask_vv,
                        const std::vector<larcv::Image2D>& wholeview_v,
                        const std::vector<larcv::Image2D>& badch_v,
                        std::vector< std::vector<MaskMatchData> >& matchdata_vv,
                        std::vector< std::vector<int> >& match_indices );
  
  };
}
}

#endif
