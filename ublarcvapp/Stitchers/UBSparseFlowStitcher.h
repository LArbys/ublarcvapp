#ifndef __UB_SPARSE_FLOW_STITCHER_H__
#define __UB_SPARSE_FLOW_STITCHER_H__

#include <vector>
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/SparseImage.h"

namespace ublarcvapp {

  class UBSparseFlowStitcher {
  public:

    UBSparseFlowStitcher( const std::vector<larcv::Image2D>& embedding_images );
    ~UBSparseFlowStitcher();
    
    void addSparseData( const larcv::SparseImage& data );

    std::vector<larcv::Image2D> _outimg_v;
    std::vector<larcv::Image2D> _metric_v;

  };


}

#endif
