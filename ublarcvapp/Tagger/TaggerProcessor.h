#ifndef __TAGGER_PROCESSOR_H__
#define __TAGGER_PROCESSOR_H__

//larlite
#include "DataFormat/storage_manager.h"
#include "SelectionTool/LEEPreCuts/LEEPreCut.h"

// larcv
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/Base/PSet.h"

// ublarcvapp
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
      virtual void configure(const larcv::PSet&); //< configure from a PSet object
      virtual void configure(const std::string);  //< configure from a file
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
      void runPMTPrecuts( larlite::storage_manager* llio );
      void runBoundaryTagger( const InputPayload& input, ThruMuPayload& output );
      
      
      // config parameters
      TaggerCROIAlgoConfig m_config; //< holds configuration parameters
      bool _RunPMTPrecuts;
      bool _ApplyPMTPrecuts;

      // timing monitor
      std::vector<float> m_time_tracker;
      enum Stages_t { kInputTotal=0, kPMTPrecuts,
                      kThruMuConfig, kThruMuContour, kThruMuBMT,
                      kThruMuFlash, kThruMuFilter, kThruMuTracker,
                      kStopMuTracker, kRecluster,
                      kUntagged, kPCAmerge, kCROI, kNumStages };

      
      // sub-algorithms
      //BMTCV m_bmtcv_algo;
      
    };
    
  }
}

#endif
