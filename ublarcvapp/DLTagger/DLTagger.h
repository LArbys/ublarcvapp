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
      hasRun(false),
      _cut_dtick_outoftime(-100.0),
      _cut_frac_out_of_croi(0.025),
      _cut_frac_dwall_innermost(-10.0),
      _cut_frac_dwall_outermost(-10.0),
      _max_pixdist_from_path(8.0)
    {};
    virtual ~DLTagger() {};

    void runTagger( const std::vector<larcv::Image2D>& wholeview_v,
                    const larcv::EventChStatus& ev_chstatus,
                    const std::vector<larlite::opflash>& intime_opflash_v,                    
                    const std::vector< std::vector<larcv::ClusterMask> >& clustermask_vv );

    void recoTruthMatching( const std::vector<larcv::Image2D>& mcinstance_v );
    
    void transferImages( std::vector<larcv::Image2D>& cosmic,
                         std::vector<larcv::Image2D>& notcosmic_v );
    
    void transferPixelClusters( larcv::EventPixel2D& cosmic_clusters,
                                larcv::EventPixel2D& notcosmic_clusters );

    void transferCROI( larcv::EventROI& ev_croi_v, larcv::EventROI& ev_croi_merged );
    
    bool hasData() { return hasRun; };

    struct CosmicSelectVars_t {
      int astar_complete;
      std::vector<int>   numpixels;
      float dwall_outermost;
      float dwall_innermost;
      float dtick_outoftime;
      float total_frac;
      std::vector<float> frac_per_plane;
      float total_nufrac;
      std::vector<float> nufrac_per_plane;
      std::vector<float> outermost_endpt_tyz;
      std::vector<float> innermost_endpt_tyz;      
      
      CosmicSelectVars_t()
      : astar_complete(0),
        dwall_outermost(-100),
        dwall_innermost(-100),
        dtick_outoftime(0),
        total_frac(0),
        total_nufrac(0)
      {
        numpixels.resize(3,0);
        frac_per_plane.resize(3,0);
        nufrac_per_plane.resize(3,0);
        outermost_endpt_tyz.resize(3,0);
        innermost_endpt_tyz.resize(3,0);
      };
        
    };
    const std::vector<CosmicSelectVars_t>& getSelectionVars() const { return m_select_vars_v; };

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
    void _calcDwall( const Gen3DEndpoints& endpts, float& outermost_dwall, float& innermost_dwall,
                     std::vector<float>& outermost_endpt_tyz,
                     std::vector<float>& innermost_endpt_tyz );
    float _calcOutOfTimeTicks( const Gen3DEndpoints& endpts );
    void _calcOutOfROIfrac( const std::vector<larcv::Image2D>& wholeview_v,
                            const larcv::ROI& croi,
                            const std::vector<const larcv::Pixel2DCluster*>& ppixel_cluster_v,
                            std::vector<float>& frac_per_plane,
                            std::vector<int>& npixs,
                            float& total_frac );
    
    // Cut Variable Values
    // -------------------
    float _cut_dtick_outoftime;
    float _cut_frac_out_of_croi;
    float _cut_frac_dwall_innermost;
    float _cut_frac_dwall_outermost;
    float _max_pixdist_from_path;
    
  };

}
}


#endif
