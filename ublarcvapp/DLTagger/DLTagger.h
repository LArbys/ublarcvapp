#ifndef __DLTAGGER_H__
#define __DLTAGGER_H__

// mrcnnmatch

#include "MRCNNMatch.h"

namespace ublarcvapp {
namespace dltagger {

  class DLTagger {
  public:
    DLTagger() {};
    virtual ~DLTagger() {};

    MRCNNMatch _mask_match_algo;
    
  };

}
}


#endif
