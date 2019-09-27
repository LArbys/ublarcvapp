#ifndef __DLSSNetOverlayCheater_H__
#define __DLSSNetOverlayCheater_H__

/**
 * \class DLSSNetOverlayCheater
 * \brief This module makes ssnet output using the truth info from the MC 'segment' image.
 *
 * This is meant to work with overlay samples. 
 * The  module can be run in two modes. The first only produces SSNet values where 
 * there are segment image values, i.e. only for MC-generated pixels.
 * The second more keeps ssnet values on the cosmic muon pixels, but only modifies the values
 * on the MC neutrino.
 *
 */

#include <string>
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

namespace ublarcvapp {
namespace dltagger {

  class DLSSNetOverlayCheater : public larcv::ProcessBase {
  public:

    DLSSNetOverlayCheater( std::string instance_name )
      : larcv::ProcessBase(instance_name)
      {};
    virtual ~DLSSNetOverlayCheater() {};

    // required functions
    void configure(const larcv::PSet& );
    void initialize();
    bool process( larcv::IOManager& mgr );
    void finalize();    

  protected:

    std::string _input_adc_name;
    std::string _input_ssnet_input_stem;
    std::string _input_segment_name;
    std::string _output_stem_treename;
    bool        _mask_cosmic_pixels;
    bool        _apply_threshold;
    std::vector<float> _adc_threshold_v;
    
  };

  class DLSSNetOverlayCheaterFactory : public larcv::ProcessFactoryBase {
  public:
    DLSSNetOverlayCheaterFactory() {
      // insert into ProcessorFactor
      larcv::ProcessFactory::get().add_factory("DLSSNetOverlayCheater",this);
    };
    ~DLSSNetOverlayCheaterFactory() {};
    
    larcv::ProcessBase* create(const std::string instance_name) {
      return new ublarcvapp::dltagger::DLSSNetOverlayCheater(instance_name);
    };
    
  };
  
}
}

#endif
