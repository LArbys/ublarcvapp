#ifndef __UBLARCVAPP_MCTOOLS_MCPIXEL_LABEL_MAKER_H__
#define __UBLARCVAPP_MCTOOLS_MCPIXEL_LABEL_MAKER_H__

#include <array>
#include <map>

#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/IOManager.h"
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"

#include "ublarcvapp/MCTools/MCParticleGraph.h"

#include "EventMCPixelLabels.h"

namespace ublarcvapp {
namespace mctools {

class MCPixelLabelMaker : public larcv::larcv_base {

public:

  MCPixelLabelMaker()
  : larcv::larcv_base("MCPixelLabelMaker")
  {};

  virtual ~MCPixelLabelMaker() {};

  void process( larlite::storage_manager& ioll, 
                larcv::IOManager& iolcv,
                std::string image2d_tree_name );

  void make_truthlabels_fromsimch(
      std::string image2d_tree_name,
      larlite::storage_manager& ioll, 
      larcv::IOManager& iolcv,
      ublarcvapp::mctools::MCParticleGraph& mcpg,
      larutil::SpaceChargeMicroBooNE* psce );

  void export_as_hdf(std::string hdf_outfile);

  EventMCPixelLabels _pixels_v;


};

}
}


#endif