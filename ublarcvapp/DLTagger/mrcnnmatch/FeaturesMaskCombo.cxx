#include "FeaturesMaskCombo.h"

#include <cilantro/principal_component_analysis.hpp>

namespace ublarcvapp {
namespace dltagger {

  /**
   * class containing contour and PCA information for cropped images
   *
   */
  FeaturesMaskCombo::FeaturesMaskCombo( CropMaskCombo& cropdata )
    : pcropdata(&cropdata)
  {
    _make_contours(cropdata);
    _calc_mask_pca(cropdata);
  }

  /**
   * create contour clusters for both charge and mask pixels
   *
   */
  void FeaturesMaskCombo::_make_contours( const CropMaskCombo& cropdata) {

    // contour clusters on charge crops
    combo_charge_contour.analyzeImages( cropdata.crops_v );
    combo_charge_contour.clear_intermediate_images();

    // contour on mask clusters
    combo_mask_contour.analyzeImages( cropdata.mask_v, 0.5, 1 );
    combo_mask_contour.clear_intermediate_images();
    
  }

  /**
   * using the mask pixels, we perform PCA to get a rough direction of the entire mask
   *
   * we use cilantro to do this.
   */
  void FeaturesMaskCombo::_calc_mask_pca( const CropMaskCombo& cropdata ) {

    
    for ( size_t p=0; p< cropdata.mask_v.size(); p++ ) {
      auto const& mask = cropdata.mask_v[p];
      if ( mask.meta().cols()==0 || mask.meta().rows()==0 ) {
        // empty. fill blank
        std::vector<float> mean(2,0);
        std::vector<float> pca1(2,0);
        std::vector<float> pca2(2,0);
        std::vector<float> eigenvals(2,0);
        
        pca_mean_vv.push_back( mean );
        pca1_dir_vv.push_back( pca1 );
        pca2_dir_vv.push_back( pca2 );
        pca_eigenvalues_vv.push_back(eigenvals);
        continue;
      }

      std::vector< Eigen::Vector3f > pt_v;
      pt_v.reserve( int(mask.meta().cols()*mask.meta().rows()*0.05) );
      for ( size_t c=0; c<mask.meta().cols(); c++ ) {
        for ( size_t r=0; r<mask.meta().rows(); r++ ) {
          if ( mask.pixel(r,c)>0.5 ) {
            pt_v.push_back( Eigen::Vector3f((float)c,(float)r,0.0));
          }
        }
      }

      if ( pt_v.size()==0 ) {
        // empty. fill blank
        std::vector<float> mean(2,0);
        std::vector<float> pca1(2,0);
        std::vector<float> pca2(2,0);
        std::vector<float> eigenvals(2,0);
        
        pca_mean_vv.push_back( mean );
        pca1_dir_vv.push_back( pca1 );
        pca2_dir_vv.push_back( pca2 );
        pca_eigenvalues_vv.push_back(eigenvals);
        continue;
      }
      
      cilantro::PrincipalComponentAnalysis3f pca(pt_v);

      std::vector<float> mean(2);
      std::vector<float> pca1(2);
      std::vector<float> pca2(2);
      std::vector<float> eigenvals(2);

      float pca1len = 0.;
      float pca2len = 0.;
      for (int i=0; i<2; i++) {
        mean[i] = pca.getDataMean()(i);
        eigenvals[i] = pca.getEigenValues()(i);
        pca1[i] = pca.getEigenVectors()(i,0); // first axis
        pca2[i] = pca.getEigenVectors()(i,1); // second axis
        pca1len += pca1[i]*pca1[i];
        pca2len += pca2[i]*pca2[i];
      }
      pca1len = sqrt(pca1len);
      pca2len = sqrt(pca2len);
      for (int i=0; i<2; i++ ) {
        pca1[i] /= pca1len;
        pca2[i] /= pca2len;
      }
      
      pca_mean_vv.push_back( mean );
      pca1_dir_vv.push_back( pca1 );
      pca2_dir_vv.push_back( pca2 );
      pca_eigenvalues_vv.push_back(eigenvals);
      
      
    }//end of loop over planes
    
  }

}
}
