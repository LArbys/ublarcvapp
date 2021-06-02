#include "TruthShowerTrunkSCE.h"

namespace ublarcvapp {
namespace mctools {

    /**
   * @brief default constructor
   *
   *  This instance will make its own SpaceChargeMicroBooNE instance
   *
   */
  TruthShowerTrunkSCE::TruthShowerTrunkSCE()
    : larcv::larcv_base("TruthShowerTrunkSCE"),
    _kown_sce(true),
    _p_sce(nullptr)
  {
    _p_sce = new larutil::SpaceChargeMicroBooNE();
  }

  /**
   * @brief constructor using external spacechargemicroboone instance
   *
   * @param psce Pointer to instance of SpaceChargeMicroBooNE
   */  
  TruthShowerTrunkSCE::TruthShowerTrunkSCE( larutil::SpaceChargeMicroBooNE* psce )
    : larcv::larcv_base("TruthShowerTrunkSCE"),
    _kown_sce(false),
    _p_sce(psce)
  {
  }

  TruthShowerTrunkSCE::~TruthShowerTrunkSCE()
  {
    if ( _kown_sce && _p_sce )
      delete _p_sce;
    _p_sce = nullptr;
  }

  /**
   * @brief define shower trunk that takes space-charge effect into account
   *
   * @param[in] shower Truth info about simulated electron or photon
   * @return track with 2 points
   */
  larlite::track TruthShowerTrunkSCE::applySCE( const larlite::mcshower& shower )
  {
    
    larlite::track trunk;

    auto const& profile = shower.DetProfile();
    TVector3 dir = shower.Start().Momentum().Vect();
    float pmom = dir.Mag();
    TVector3 vstart = shower.Start().Position().Vect();
    TVector3 pstart = shower.DetProfile().Position().Vect();

    if ( shower.PdgCode()==22 ) {
      // for gamma use shower profile
      vstart = pstart;
      dir = profile.Momentum().Vect();
      pmom = dir.Mag();
    }

    if ( dir.Mag()<0.1 )
      return trunk;
      
    std::vector<float> fdir(3,0);
    std::vector<float> fstart(3,0);
    std::vector<float> fend(3,0);
    TVector3 vend;
    for (int i=0; i<3; i++) {
      dir[i] /= pmom;
      fdir[i] = (float)dir[i];
      fstart[i] = vstart[i];
      fend[i] = fstart[i] + 10.0*fdir[i];
      vend[i] = fend[i];
    }
    
    // space charge correction if not gamma
    TVector3 v3end = { vend[0], vend[1], vend[2] };
      
    if ( shower.PdgCode()!=22 ) {
      std::vector<double> s_offset = _p_sce->GetPosOffsets(vstart[0],vstart[1],vstart[2]);
      vstart[0] = fstart[0] - s_offset[0] + 0.7;
      vstart[1] = fstart[1] + s_offset[1];
      vstart[2] = fstart[2] + s_offset[2];
      
      std::vector<double> e_offset = _p_sce->GetPosOffsets(v3end[0],v3end[1],v3end[2]);
      v3end[0] = vend[0] - e_offset[0] + 0.7;
      v3end[1] = vend[1] + e_offset[1];
      v3end[2] = vend[2] + e_offset[2];
    }
    
    TVector3 sce_dir = v3end-vstart;
    float sce_dir_len = sce_dir.Mag();
    if ( sce_dir_len>0 ) {
      for (int i=0; i<3; i++)
        sce_dir[i] /= sce_dir_len;
    }      
    
    // make shower trunk object
    trunk.reserve(2);
    trunk.add_vertex( vstart );
    trunk.add_vertex( v3end );
    trunk.add_direction( sce_dir );
    trunk.add_direction( sce_dir );      

    return trunk;
    
  }

  
}
}
