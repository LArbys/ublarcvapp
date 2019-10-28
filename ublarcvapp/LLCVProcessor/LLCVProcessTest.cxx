#include "LLCVProcessTest.h"

namespace ublarcvapp {
namespace llcv {

  static LLCVProcessTestFactory __global_LLCVProcessTestFactory__;

  void LLCVProcessTest::configure( const larcv::PSet& pset ) {
    LARCV_NORMAL() << "called configure(pset)" << std::endl;
  }

  void LLCVProcessTest::initialize() {
    LARCV_NORMAL() << "called initialize" << std::endl;
  }

  bool LLCVProcessTest::process( larcv::IOManager& mgr, larlite::storage_manager& llio ) {
    LARCV_NORMAL() << "called process(IOManager,storage_manager)" << std::endl;    
    return true;
  }

  void LLCVProcessTest::finalize() {
    LARCV_NORMAL() << "called finalize()" << std::endl;
  }

}
}
