#ifndef __DLTAGGER_H__
#define __DLTAGGER_H__

#include <vector>
#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/Pixel2DCluster.h"
#include "larcv/core/DataFormat/EventPixel2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"
#include "larcv/core/DataFormat/EventROI.h"
#include "larcv/core/DataFormat/ClusterMask.h"

// larlite
#include "DataFormat/opflash.h"

// mrcnnmatch
#include "MRCNNMatch.h"


namespace ublarcvapp {
namespace dltagger {

  class DLTagger : public larcv::larcv_base {

  public:
    DLTagger()
      : larcv::larcv_base("DLTagger"),
      hasRun(false)
    {};
    virtual ~DLTagger() {};

    void runTagger( const std::vector<larcv::Image2D>& wholeview_v,
                    const larcv::EventChStatus& ev_chstatus,
                    const std::vector<larlite::opflash>& intime_opflash_v,                    
                    const std::vector< std::vector<larcv::ClusterMask> >& clustermask_vv );

    void transferImages( std::vector<larcv::Image2D>& cosmic,
                         std::vector<larcv::Image2D>& notcosmic_v );
    
    void transferPixelClusters( larcv::EventPixel2D& cosmic_clusters,
                                larcv::EventPixel2D& notcosmic_clusters );

    void transferCROI( larcv::EventROI& ev_croi_v, larcv::EventROI& ev_croi_merged );
    
    bool hasData() { return hasRun; };

    struct CosmicSelectVars_t {
      float dwall_outermost;
      float dwall_innermost;
      float dtick_outoftime;
      std::vector<float> frac_per_plane;
      float total_frac;
    };

    void reset();
      
    
  protected:

    // Mask-RCNN Matching Across Planes
    // ---------------------------------
    MRCNNMatch _mask_match_algo;
    std::vector< std::vector<larcv::Pixel2DCluster> > m_pixel_cluster_vv;
    std::vector< larcv::ImageMeta > m_pixel_cluster_meta_v;
    std::vector<larcv::Image2D> m_tagged_v;    ///< wholeview image tagged with cosmic clusters
    std::vector<larcv::Image2D> m_notcosmic_v; ///< wholeview image tagged with not-cosmic clusters
    bool hasRun;

    void _makePixelClusters( const std::vector<larcv::Image2D>& wholeview_v,
                             const MRCNNMatch& matchdata,
                             std::vector< std::vector<larcv::Pixel2DCluster> >& pixel_cluster_vv,
                             std::vector< larcv::ImageMeta >& pixel_cluster_meta_v );

    void _tagPixels( const std::vector<larcv::Image2D>& wholeview_v,
                     const std::vector< std::vector<larcv::Pixel2DCluster> >& pixel_cluster_vv,
                     const std::vector<int>& iscosmic_v,
                     std::vector<larcv::Image2D>& whole_cosmic_v,
                     std::vector<larcv::Image2D>& whole_notcosmic_v );
    
    larcv::ROI _mergeROI( const int nplanes,
                          const std::vector<larcv::ROI>& roi_v );
    
    void _evaluateClusters( const std::vector<larcv::Image2D>& wholeview_v,
                            const MRCNNMatch& matchdata,
                            const larcv::ROI& croimerged,
                            const std::vector< std::vector<larcv::Pixel2DCluster> >& pixel_cluster_vv,
                            std::vector<int>& iscosmic_v );
    
    // Flash CROI
    // -----------
    std::vector<larcv::ROI> m_croi_v;
    larcv::ROI              m_croi_merged;

    // Cosmic Selection
    std::vector<int> m_iscosmic_v;
    std::vector<CosmicSelectVars_t> m_select_vars_v;
    void _calcDwall( const Gen3DEndpoints& endpts, float& outermost_dwall, float& innermost_dwall );
    float _calcOutOfTimeTicks( const Gen3DEndpoints& endpts );
    void _calcOutOfROIfrac( const std::vector<larcv::Image2D>& wholeview_v,
                            const larcv::ROI& croi,
                            const std::vector<const larcv::Pixel2DCluster*>& ppixel_cluster_v,
                            std::vector<float>& frac_per_plane,
                            float& total_frac );
    
    
    
  };

}
}


#endif
