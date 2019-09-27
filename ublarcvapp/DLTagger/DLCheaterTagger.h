#ifndef __DL_CHEATER_TAGGER_H__
#define __DL_CHEATER_TAGGER_H__

#include <string>
#include "larcv/core/Processor/ProcessBase.h"

namespace ublarcvapp {
namespace dltagger {

  class DLCheaterTagger : public larcv::ProcessBase {
  public:

    DLCheaterTagger() {};
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
  
}
}

#endif
