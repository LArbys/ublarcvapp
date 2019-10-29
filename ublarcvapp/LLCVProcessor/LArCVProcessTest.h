#ifndef __LARCV_PROCESS_TEST_H__
#define __LARCV_PROCESS_TEST_H__


/**
 * 
 * \brief A test module that inherits from larcv::ProcessBase only
 *
 */

#include <string>
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

namespace ublarcvapp {
namespace llcv {

  class LArCVProcessTest : public larcv::ProcessBase {

  public:

    LArCVProcessTest( std::string name )
      : larcv::ProcessBase(name)
      {};

    virtual ~LArCVProcessTest() {};

    void configure( const larcv::PSet& pset ) {};
    void initialize() {};
    bool process( larcv::IOManager& io );
    void finalize() {};
    

  };


  class LArCVProcessTestFactory : public larcv::ProcessFactoryBase {
  public:
    LArCVProcessTestFactory() { larcv::ProcessFactory::get().add_factory("LArCVProcessTest",this); };
    ~LArCVProcessTestFactory() {};
    larcv::ProcessBase* create(const std::string instance_name) { return new LArCVProcessTest(instance_name); };
  };
  
  
}
}

#endif
