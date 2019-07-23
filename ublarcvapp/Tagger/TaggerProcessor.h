#ifndef __TAGGER_PROCESSOR_H__
#define __TAGGER_PROCESSOR_H__

#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/Base/PSet.h"

#include "ublarcvapp/LArliteHandler/LArliteManager.h"

#include "TaggerCROITypes.h"
#include "TaggerCROIAlgoConfig.h"

namespace ublarcvapp {
  namespace tagger {


    class TaggerProcessor : public larcv::ProcessBase {
    public:
      TaggerProcessor( std::string cfg_path="" );
      
      virtual ~TaggerProcessor();
      
      // required implementations of pure virtual functions from larcv::ProcessBase
      virtual void configure(const larcv::PSet&);
      virtual bool process(larcv::IOManager& io);
      virtual void initialize();      
      virtual void finalize();
      
      std::string _cfg_path; //< path to configuration file (if done directly)
      
      ublarcvapp::LArliteManager* _larlite_io; //< handles interface/sync to larlite info
      void add_larlite_input( std::string llinput );


      // primary functions
      void loadInput( larcv::IOManager& io_larcv,
                      ublarcvapp::LArliteManager& io_larlite,
                      InputPayload& input );
      
      
      // config parameters
      TaggerCROIAlgoConfig m_config;
      
      
    };
    
  }
}

#endif
