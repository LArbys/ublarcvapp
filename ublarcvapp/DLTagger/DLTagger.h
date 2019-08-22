#ifndef __DLTAGGER_H__
#define __DLTAGGER_H__

#include <vector>
#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"
#include "larcv/core/DataFormat/ClusterMask.h"

// mrcnnmatch
#include "MRCNNMatch.h"

namespace ublarcvapp {
namespace dltagger {

  class DLTagger : larcv::larcv_base {

  public:
    DLTagger()
      : larcv::larcv_base("DLTagger"),
      hasRun(false)
    {};
    virtual ~DLTagger() {};

    void runTagger( const std::vector<larcv::Image2D>& wholeview_v,
                    const larcv::EventChStatus& ev_chstatus,
                    const std::vector< std::vector<larcv::ClusterMask> >& clustermask_vv );

    void transferImages( std::vector<larcv::Image2D>& container );
    bool hasData() { return hasRun; };
    
  protected:

    //bool hasRun;
    MRCNNMatch _mask_match_algo;
    std::vector<larcv::Image2D> m_tagged_v;
    bool hasRun;
    
    void _tagPixels( const std::vector<larcv::Image2D>& wholeview_v,
                     std::vector<larcv::Image2D>& tagged_v );
    
  };

}
}


#endif
