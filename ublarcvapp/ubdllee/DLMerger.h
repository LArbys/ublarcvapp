#ifndef __DL_MERGER_H__
#define __DL_MERGER_H__

#include <string>
#include <vector>
#include <map>

#include "larcv/core/Processor/ProcessDriver.h"

namespace ublarcvapp {
namespace ubdllee {

  class DLMerger : public larcv::ProcessDriver {

  public:

    DLMerger(std::string drivername="DLMerger")
      : larcv::ProcessDriver(drivername),
      _file_dir("")
        {};

    virtual ~DLMerger() {};

    void configure( const std::string process_driver_cfg,
                    const std::string json_filelist,
                    const std::string file_dir="" );

    void initialize();

  protected:

    std::string _file_dir;
    std::map< std::string, std::vector<std::string> > _filelist_map;
    
  };
  
}
}

#endif
