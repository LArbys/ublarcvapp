#include "LLCVProcessBase.h"


namespace ublarcvapp {
namespace llcv {

  LLCVProcessBase::LLCVProcessBase( const std::string name )
    : larcv::ProcessBase(name)
  {}
  
  bool LLCVProcessBase::process( larcv::IOManager& mgr ) {
    std::string errmsg = "Called process(IOManager&)\
. This process (name=" + name() + ") is meant to be used with LLCVProcessDriver, \
which calles process(larcv::IOManager, larlite::storage_manager).";
    LARCV_CRITICAL() << errmsg << std::endl;
    throw std::runtime_error(errmsg.c_str());
  }
  
}
}
