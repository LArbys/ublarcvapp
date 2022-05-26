#include "NeutrinoVertex.h"

#include "larlite/LArUtil/Geometry.h"
#include "larlite/DataFormat/mctruth.h"

#include "crossingPointsAnaMethods.h"

namespace ublarcvapp {
namespace mctools {

  /*
   * get image coordinate for neutrino vertex
   * 
   * @param[in] ioll The larlite storage manager that contains mctruth class
   * @return Image coordinate vector with entries: (u,v,y,tick)
   */
  std::vector<int> NeutrinoVertex::getImageCoords( larlite::storage_manager& ioll )
  {
    std::vector<float> pos3dt_wsce = getPos3DwSCE( ioll );
    TVector3 dpos(0,0,0);
    for (int i=0; i<3; i++) dpos[i] = (double)pos3dt_wsce[i];
    std::vector<int> imgcoord(4,0); // U,V,Y,tick
    for (int p=0; p<3; p++ )
      imgcoord[p] = (int)(larlite::larutil::Geometry::GetME()->WireCoordinate(dpos,p)+0.5);//0.5 for rounding
    imgcoord[3] = pos3dt_wsce[3];
    return imgcoord;
  }
  
  /**
   * get the neutrino vertex in detector coordinates, UNcorrected by space-charge effect
   *
   * @param[in] ioll The larlite storage manager that contains mctruth class
   * @return         length 3 vector with position in detector and tick: (x,y,z,tick)
   *
   */
  std::vector<float> NeutrinoVertex::getPos3D( larlite::storage_manager& ioll )
  {

    larlite::event_mctruth* ev_mctruth
      = (larlite::event_mctruth*)ioll.get_data(larlite::data::kMCTruth,"generator");

    auto const& mctruth = ev_mctruth->at(0);
    const larlite::mcstep& start = mctruth.GetNeutrino().Nu().Trajectory().front();
    float tick = CrossingPointsAnaMethods::getTick(start, 4050.0, nullptr);
    std::vector<float> pos3d;
    pos3d[0] = start.X();
    pos3d[1] = start.Y();
    pos3d[2] = start.Z();
    pos3d[3] = tick;
    return pos3d;
  }
  
  /**
   * get the neutrino vertex in detector coordinates, corrected by space-charge effect
   *
   * @param[in] ioll The larlite storage manager that contains mctruth class
   * @param[in] psce Pointer to a SpaceChargeMicroBooNE instance. If NULL provided, 
   *                 creates and destroys instance inside function.
   * @return         length 3 vector with position in detector and tick: (x,y,z,tick)
   */
  std::vector<float> NeutrinoVertex::getPos3DwSCE( larlite::storage_manager& ioll,
                                                   larutil::SpaceChargeMicroBooNE* psce )
  {

    larutil::SpaceChargeMicroBooNE* _sce = nullptr;
    if ( psce ) _sce = psce;
    else {
      _sce = new larutil::SpaceChargeMicroBooNE;
    }
    
    larlite::event_mctruth* ev_mctruth
      = (larlite::event_mctruth*)ioll.get_data(larlite::data::kMCTruth,"generator");

    auto const& mctruth = ev_mctruth->at(0);
    const larlite::mcstep& start = mctruth.GetNeutrino().Nu().Trajectory().front();
    float tick = CrossingPointsAnaMethods::getTick(start, 4050.0, _sce);
    std::vector<float> pos3d_wtick(4,0);
    pos3d_wtick[0] = start.X();
    pos3d_wtick[1] = start.Y();
    pos3d_wtick[2] = start.Z();

    std::vector<double> offset = _sce->GetPosOffsets( start.X(), start.Y(), start.Z() );
    pos3d_wtick[0] += -offset[0] + 0.7;
    pos3d_wtick[1] += offset[1];
    pos3d_wtick[2] += offset[2];
    pos3d_wtick[3] = tick;

    // clearn up space charge
    if ( !psce )
      delete _sce;

    return pos3d_wtick;
  }
  
}
}
    
