#ifndef __UB_SPARSE_FLOW_STITCHER_H__
#define __UB_SPARSE_FLOW_STITCHER_H__

#include <vector>
#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/SparseImage.h"

namespace ublarcvapp {

  class UBSparseFlowStitcher : public larcv::larcv_base {
  public:

    UBSparseFlowStitcher( const std::vector<larcv::Image2D>& embedding_images );
    ~UBSparseFlowStitcher();
    
    void addSparseData( const larcv::SparseImage& data,
                        const larcv::ImageMeta& tar1_subimg,
                        const larcv::ImageMeta& tar2_subimg );
                        
    void clear();

    std::vector<larcv::Image2D> _outimg_v;
    std::vector<larcv::Image2D> _metric_v;

  };


}

#endif
