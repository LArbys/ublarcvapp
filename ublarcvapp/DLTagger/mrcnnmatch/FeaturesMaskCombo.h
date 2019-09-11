#ifndef __FEATURES_MASK_COMBO_H__
#define __FEATURES_MASK_COMBO_H__

#include <vector>
#include "ublarcvapp/ContourTools/ContourClusterAlgo.h"
#include "CropMaskCombo.h"

namespace ublarcvapp {
namespace dltagger {

  /**
   * class stores features extracted for a combination of 
   * of MRCNN masks
   *
   */
  class FeaturesMaskCombo {
  public:
    
    FeaturesMaskCombo()
      : pcropdata(nullptr)
    {}; ///< default constructor. do not use. for ROOT dictionary building.

    FeaturesMaskCombo( CropMaskCombo& cropdata,
                       bool extend_mask_wpca=true,
                       float max_dist_2_pcaline=10.0 );
    virtual ~FeaturesMaskCombo() {};

    CropMaskCombo* pcropdata; ///< inputs used to generate information from this class

    // outputs: contours
    ublarcvapp::ContourClusterAlgo combo_charge_contour; // contours on the charge
    ublarcvapp::ContourClusterAlgo combo_mask_contour;   // contours on the mask
    
    // outputs: global pca
    std::vector< std::vector<float> > pca_mean_vv;        // mean position of cluster pixels for each plane
    std::vector< std::vector<float> > pca1_dir_vv;        // 1st PCA-axis direction for each plane
    std::vector< std::vector<float> > pca2_dir_vv;        // 2nd PCA-axis direction for each plane
    std::vector< std::vector<float> > pca_eigenvalues_vv; // eigenvalues for each plane

    void extendMaskWithPCAregion( const float max_dist_2_pcaline=10.0 );
    
  protected:

    void _make_contours( const CropMaskCombo& cropdata );
    void _calc_mask_pca( const CropMaskCombo& cropdata );


  };

}
}

#endif
