#ifndef __INFILL_IMAGE_STITCHER__
#define __INFILL_IMAGE_STITCHER__

#include <vector>
#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"



namespace ublarcvapp {

    class InfillImageStitcher {
    public:
      InfillImageStitcher(){};
      virtual ~InfillImageStitcher(){};

      void Croploop(larcv::ImageMeta& output_meta,
                      larcv::Image2D& outimg,
                      larcv::Image2D& outputimg,
                      larcv::Image2D& overlapcountimg);

      void Overlayloop(int p,larcv::ImageMeta& output_meta,
                      larcv::Image2D& outputimg,
                      larcv::Image2D& overlapcountimg,
                      const std::vector<larcv::Image2D>& wholeview_v,
                      larcv::EventChStatus& ev_chstatus);
    };


}

#endif
