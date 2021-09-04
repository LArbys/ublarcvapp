#ifndef __LLCV_PROCESS_TEST_H__
#define __LLCV_PROCESS_TEST_H__

#include "LLCVProcessBase.h"

#include <string>
#include "larcv/core/Processor/ProcessFactory.h"
#include "larcv/core/DataFormat/IOManager.h"
#include "larlite/DataFormat/storage_manager.h"

namespace ublarcvapp {
namespace llcv {

  class LLCVProcessTest : public LLCVProcessBase {
  public:
    LLCVProcessTest( const std::string name )
      : LLCVProcessBase(name)
    {};

    void initialize();
    void configure( const larcv::PSet& pset );
    bool process( larcv::IOManager& mgr, larlite::storage_manager& llio );
    void finalize();
    
  };

  class LLCVProcessTestFactory : public larcv::ProcessFactoryBase {
  public:
    LLCVProcessTestFactory() { larcv::ProcessFactory::get().add_factory("LLCVProcessTest",this); };
    ~LLCVProcessTestFactory() {};
    larcv::ProcessBase* create(const std::string instance_name) { return new LLCVProcessTest(instance_name); };
  };
  
  
}
}



#endif
