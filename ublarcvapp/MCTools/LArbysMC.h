#ifndef __LARLITECV_LARBYSMC_H__
#define __LARLITECV_LARBYSMC_H__

#include <vector>
#include <string>

#include "TTree.h"

#include "DataFormat/storage_manager.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/IOManager.h"

namespace larutil {
  class SpaceChargeMicroBooNE;
}

namespace ublarcvapp {
namespace mctools {

  /**
   * @class LArbysMC
   * @brief Calculate truth variables useful for low energy excess analysis
   */
  class LArbysMC  {

  public:

    /// Default constructor
    LArbysMC();

    /// Default destructor
    ~LArbysMC();

    //void configure(const PSet&);
    void initialize();
    bool process(larlite::storage_manager& mgr);
    bool process(larcv::IOManager& iolcv, larlite::storage_manager& mgr);
    void finalize();
    void Clear();

    void bindAnaVariables( TTree* );
    void printInteractionInfo();
    void setWireProducer( std::string name ) { _producer_wireimage=name; };

    void calculateVertexPixsum( larcv::IOManager& iolcv,
                                float vtx_tick, float vtx_wireu, float vtx_wirev, float vtx_wirey,
                                int pix_radius,
                                std::vector<float>& plane_vtx_pixsum  );
    
  protected:

    std::string _name;
    
    std::string _producer_mctruth;
    std::string _producer_mctrack;
    std::string _producer_mcshower;    

    TTree* _mc_tree;

    larutil::SpaceChargeMicroBooNE* _psce;

  public:
    
    /// Event ID
    int _run;    ///< run number
    int _subrun; ///< subrun number
    int _event;  ///< event number
    int _entry;  ///< entry number

    std::string _rse_producer; ///< larlite data product to use to get (run,subrun,event)

    // neutrino interaction info
    bool _neutrino_present; ///< event has neutrino
    int _current_type;      ///< CC=0, NC=1
    int _interaction_type;  ///< Nuance interaction mode
    int _genie_mode;        ///< GENIE interaction mode

    float _Enu_true;    ///< true neutrino energy
    float _vtx_t;       ///< true neutrino vertex time
    float _vtx_x;       ///< true neutrino vertex x-position
    float _vtx_y;       ///< true neutrino vertex y-position
    float _vtx_z;       ///< true neutrino vertex z-position
    float _vtx_sce_x;   ///< true neutrino vertex x-position, space charge applied
    float _vtx_sce_y;   ///< true neutrino vertex y-position, space charge applied
    float _vtx_sce_z;   ///< true neutrino vertex z-position, space charge applied
    float _vtx_detx;    ///< true neutrino vertex non-t0-corrected x-position, space charge applied
    float _vtx_tick;    ///< true neutrino vertex tick, space charge applied
    float _vtx_wire[3]; ///< true neutrino vertex [U,V,Y] plane, space charge applied

    float _evis;        ///< total visible energy: KE sum of leptons, protons, charged pions, photon
    float _evis_had;    ///< total hadronic visible energy: KE sum of proton and charged pions
    float _evis_vtx;    ///< total KE of protons with KE<60 MeV, charged pions<35 MeV
    float _evis_lep;    ///< total KE of either primary muon or electron

    int _hi_lep_pdg;    ///< pdg of highest KE muon or electron
    float _hi_lep_e;    ///< highest lepton KE
    
    int _nprimary;      ///< number of primary final state particles
    int _ntotal;        ///< not used

    int _nproton;        ///< number of final state protons 
    int _nproton_60mev;  ///< number of final state protons above 60 MeV
    int _nlepton;        ///< number of final state primary leptons
    int _nlepton_35mev;  ///< number of final state primary leptons above 35 MeV
    int _nmeson;         ///< number of charged pions or pi0
    int _nmeson_35mev;   ///< number of charged pions > 35 MeV
    int _npi0;           ///< number of pi0
    int _nshower;        ///< number of showers
    int _nneutron;       ///< number of neutrons
    int _1l1p0pi;        ///< final state is 1 lepton>35 MeV, 1 proton>60 MeV, 0 charged pion>35 MeV, 0 pi0
    int _1l0p0pi;        ///< final state is 1 lepton>35 MeV, 0 proton>60 MeV, 0 charged pion>35 MeV, 0 pi0

    std::string  _producer_wireimage;  ///< producer name for wire image
    float        _plane_vtx_pixsum[3]; ///< vertex pixel sum
    float        _vtx_med_pixsum;      ///< vertex median pixel sum

    float _vtx_dwall;    ///< dwall using the true (no SCE) vertex
    int   _vtx_boundary; ///< nearest boundary

    struct entry_info{
      int run;
      int subrun;
      int event;
    };

    entry_info _entry_info;
    bool _is_signal;


  };

  /**
     \class larcv::LArbysMCFactory
     \brief A concrete factory class for larcv::LArbysMC
  */
  /* class LArbysMCProcessFactory : public ProcessFactoryBase { */
  /* public: */
  /*   /// ctor */
  /*   LArbysMCProcessFactory() { ProcessFactory::get().add_factory("LArbysMC",this); } */
  /*   /// dtor */
  /*   ~LArbysMCProcessFactory() {} */
  /*   /// creation method */
  /*   ProcessBase* create(const std::string instance_name) { return new LArbysMC(instance_name); } */

  /* }; */

}
}

#endif
/** @} */ // end of doxygen group
