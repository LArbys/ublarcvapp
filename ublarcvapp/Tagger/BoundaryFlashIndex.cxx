#include "BoundaryFlashIndex.h"

#include <sstream>

namespace ublarcvapp {
namespace tagger {

  BoundaryFlashIndex::BoundaryFlashIndex()
    : ivec(-1),idx(-1),popflash(NULL)
  {}
  
  BoundaryFlashIndex::BoundaryFlashIndex(int iivec, int iidx, const larlite::opflash* popf )
    : ivec(iivec),idx(iidx),popflash(popf)
  {}

  bool BoundaryFlashIndex::operator< (const BoundaryFlashIndex& rhs) const {
    if ( ivec<rhs.ivec )
      return true;
    else if ( ivec>rhs.ivec )
      return false;
    else {

      if ( idx<rhs.idx ) {
	return true;
      }
      else if ( idx>=rhs.idx )
	return false;
    }
    return false;
  }


  std::string BoundaryFlashIndex::getInfo() const {
    std::stringstream ss;
    ss << "[ivec=" << ivec << ", idx=" << idx << ", opflash_ptr=" << popflash << "]";
    return ss.str();
  }
}
}
