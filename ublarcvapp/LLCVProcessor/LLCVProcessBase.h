#ifndef __LLCV_PROCESS_BASE__
#define __LLCV_PROCESS_BASE__

#include <string>
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/DataFormat/IOManager.h"

// larlite
#include "DataFormat/storage_manager.h"

namespace larcv {
  class ProcessDriver;
  class ProcessFactory;
}

namespace ublarcvapp {
namespace llcv {

  class LLCVProcessBase : public larcv::ProcessBase {
    friend class larcv::ProcessDriver;
    friend class larcv::ProcessFactory;
    friend class LLCVProcessDriver;
    
  public:

    LLCVProcessBase( const std::string name="LLCVProcessBase" );

    virtual ~LLCVProcessBase() {};

    virtual bool process( larcv::IOManager& mgr );
    virtual bool process( larcv::IOManager& mgr, larlite::storage_manager& ) = 0;
    

  };
  
}
}

#endif
