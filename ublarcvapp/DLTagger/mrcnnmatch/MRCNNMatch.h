#ifndef __MRCNN_MATCH_H__
#define __MRCNN_MATCH_H__

#include <vector>

#include "larcv/core/DataFormat/ClusterMask.h"

#include "ublarcvapp/ContourTools/ContourClusterAlgo.h"
#include "MRCNNMatchTypes.h"
#include "CropMaskCombo.h"

namespace ublarcvapp {
namespace dltagger {

  class MRCNNMatch {
  public:
    MRCNNMatch() {};
    virtual ~MRCNNMatch() {};

    void matchMasksAcrossPlanes( const std::vector<std::vector<larcv::ClusterMask>>& clustermask_vv,
                                 const std::vector<larcv::Image2D>& wholeview_v,
                                 std::vector< std::vector<int> >& match_indices );

    std::vector< MaskCombo >     m_combo_3plane_v; // matches across planes
    std::vector< CropMaskCombo > m_combo_crops_vv; // cropped images from matches

    std::vector< ublarcvapp::ContourClusterAlgo > m_combo_charge_contour_v; // contours on the charge
    std::vector< ublarcvapp::ContourClusterAlgo > m_combo_mask_contour_v;   // contours on the mask

  };
}
}

#endif
