#ifndef __UBLARCVAPP_MCTOOLS_MCPIXEL_LABEL_MAKER_H__
#define __UBLARCVAPP_MCTOOLS_MCPIXEL_LABEL_MAKER_H__

#include <array>
#include <map>
#include <cmath>

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
  : larcv::larcv_base("MCPixelLabelMaker"),
  dwire(1),
  drow(1),
  source("largeant"),
  preverse_sce(nullptr),
  psce(nullptr)
  {};

  virtual ~MCPixelLabelMaker();

  void process( larlite::storage_manager& ioll, 
                larcv::IOManager& iolcv,
                std::string image2d_tree_name );

  void make_truthlabels_fromsimch(
      std::string image2d_tree_name,
      larlite::storage_manager& ioll, 
      larcv::IOManager& iolcv,
      ublarcvapp::mctools::MCParticleGraph& mcpg,
      larutil::SpaceChargeMicroBooNE* psce );

  void set_dwire( int dwiremod ) { dwire=std::abs(dwiremod); };

  void set_drow( int drowmod ) { drow=std::abs(drowmod); };

  void set_largeant_source() { source="largeant"; };

  void set_driftwc_source() { source="driftWC:simpleSC:Detsim"; };

  void export_as_hdf(std::string hdf_outfile);

  void clear() {
    _pixels_v.clear();
  };

  EventMCPixelLabels _pixels_v;

  int dwire; ///< create copies of MCPixelLabels with triplet indices with modified wire index, wire+dwire
  int drow;  ///< create copies of MCPixelLabels with triplet indices with modified row index, row+drow
  std::string source;
  larutil::SpaceChargeMicroBooNE* preverse_sce;
  larutil::SpaceChargeMicroBooNE* psce;


};

}
}


#endif