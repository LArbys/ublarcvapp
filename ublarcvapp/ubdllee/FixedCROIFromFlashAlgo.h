#ifndef __FIXED_CROI_FROM_FLASH_ALGO__
#define __FIXED_CROI_FROM_FLASH_ALGO__

#include <vector>

#include "larcv/core/DataFormat/ROI.h"
#include "larcv/core/DataFormat/ImageMeta.h"

#include "DataFormat/opflash.h"

#include "FixedCROIFromFlashConfig.h"

namespace ublarcvapp {
  namespace ubdllee {

    class FixedCROIFromFlashAlgo {
    public:
      FixedCROIFromFlashAlgo();
      FixedCROIFromFlashAlgo( const FixedCROIFromFlashConfig& config );
      virtual ~FixedCROIFromFlashAlgo() {};
      
      std::vector<larcv::ROI> findCROIfromFlash( const larlite::opflash& intimeflash );

    protected:
      FixedCROIFromFlashConfig m_config;

      float getFlashMeanZ( const larlite::opflash& flash );
    };

  }
}

#endif
