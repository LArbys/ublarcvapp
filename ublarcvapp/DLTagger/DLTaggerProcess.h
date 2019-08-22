#ifndef __DLTAGGER_PROCESS_H__
#define __DLTAGGER_PROCESS_H__

#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

#include "DLTagger.h"

namespace ublarcvapp {
namespace dltagger {

  class DLTaggerProcess : public larcv::ProcessBase {
  public:

    DLTaggerProcess(std::string instance_name )
      : larcv::ProcessBase(instance_name)
    {};
    virtual ~DLTaggerProcess() {};

    // require functions
    void configure(const larcv::PSet& );
    void initialize();
    bool process( larcv::IOManager& mgr );
    void finalize();

  protected:

    DLTagger m_tagger;

    // processor configuration parameters
    std::vector<std::string> _larlite_files; //< list of input larlite files
    
    std::string _input_adc_producer;      //< name of producer for input ADC images
    std::string _input_chstatus_producer; //< name of producer for ChStatus
    std::string _input_opflash_producer;  //< name of producer for in-beam window opflash
    std::string _input_ophit_producer;    //< name of producer for in-beam window ophit
    std::string _input_mask_producer;     //< name of producer containing Mask R-CNN output masks
    std::string _output_tagged_image;     //< name to give event container storing tagged whole-view images
    std::string _output_pixel_clusters;   //< name to give event container storing pixels from found tracks
    std::string _output_larlite_file;     //< name of output larlite file
    std::string _output_tracks;           //< name to give event container storing found 3D tracks
    
  };

  class DLTaggerProcessFactory : public larcv::ProcessFactoryBase {
  public:
    DLTaggerProcessFactory() { larcv::ProcessFactory::get().add_factory("DLTaggerProcess",this); };
    ~DLTaggerProcessFactory() {};
    larcv::ProcessBase* create(const std::string instance_name) { return new DLTaggerProcess(instance_name); };
  };
  
}
}

#endif
