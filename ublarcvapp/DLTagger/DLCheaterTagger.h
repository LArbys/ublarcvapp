#ifndef __DL_CHEATER_TAGGER_H__
#define __DL_CHEATER_TAGGER_H__

#include <string>
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

namespace ublarcvapp {
namespace dltagger {

  class DLCheaterTagger : public larcv::ProcessBase {
  public:

    DLCheaterTagger( std::string instance_name )
      : larcv::ProcessBase(instance_name) {};
    virtual ~DLCheaterTagger() {};

    // required functions
    void configure(const larcv::PSet& );
    void initialize();
    bool process( larcv::IOManager& mgr );
    void finalize();

  protected:
      
    std::string _input_adc_producer;
    std::string _input_chstatus_producer;
    std::string _input_instance_image_producer;
    std::string _output_tagger_image_producer;
    std::string _output_tagger_pixcluster_producer;
    
  };

  class DLCheaterTaggerFactory : public larcv::ProcessFactoryBase {
  public:
    DLCheaterTaggerFactory() { larcv::ProcessFactory::get().add_factory("DLCheaterTagger",this); };
    ~DLCheaterTaggerFactory() {};
    larcv::ProcessBase* create(const std::string instance_name) { return new ublarcvapp::dltagger::DLCheaterTagger(instance_name); };
  };
  
  
}
}

#endif
