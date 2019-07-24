#ifndef __LINEAR_3D_CHARGETAGGER_H__
#define __LINEAR_3D_CHARGETAGGER_H__

#include <vector>

// LArCV
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/Base/PSet.h"

// ublarcvapp/tagger
#include "TaggerCROITypes.h"
#include "Linear3DChargeTaggerTypes.h"
#include "Linear3DChargeTaggerConfig.h"

namespace ublarcvapp {
namespace tagger {


  class Linear3DChargeTagger  {
  public:
    Linear3DChargeTagger( const Linear3DChargeTaggerConfig& cfg ) : m_config(cfg) {}; 
    virtual ~Linear3DChargeTagger() {};

    Linear3DChargeTaggerConfig m_config;

    PointInfoList findpath( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
      const int start_row, const int goal_row, const std::vector<int>& start_cols, const std::vector<int>& goal_cols  );

    PointInfoList pointsOnTrack(const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
      const std::vector<float>& initial_point_coords, const std::vector<float>& final_point_coords, 
     const float step_size, const float min_ADC_value, const int neighborhood_size);

    bool doesNeighboringPixelHaveCharge(const larcv::Image2D& img, int central_row, int central_col, int neighborhood_size, float min_ADC_value);

    bool isNeighboringPixelDeadPixel(const larcv::Image2D& badch_v, int central_row, int central_col, int neighborhood_size);

    BMTrackCluster3D makeTrackCluster3D( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, 
      const BoundarySpacePoint& start_pt, const BoundarySpacePoint& end_pt, const PointInfoList& infolist );

    void getTrackExtension( const PointInfoList& infolist, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
      const float max_extension_length, PointInfoList& start_extension, PointInfoList& end_extension );


  };

}
}

#endif
