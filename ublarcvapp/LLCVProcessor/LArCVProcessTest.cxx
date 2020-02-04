#include "LArCVProcessTest.h"

namespace ublarcvapp {
namespace llcv {

  static LArCVProcessTestFactory __global_LArCVProcessTestFactory__;

  bool LArCVProcessTest::process( larcv::IOManager& mgr ) {
    LARCV_NORMAL() << "called process(IOManager)" << std::endl;
    return true;
  }

}
}
  
