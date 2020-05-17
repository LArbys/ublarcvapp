#ifndef __LARLITECV_LARBYSMC_H__
#define __LARLITECV_LARBYSMC_H__

#include <vector>
#include <string>

#include "TTree.h"

#include "DataFormat/storage_manager.h"

namespace larutil {
  class SpaceChargeMicroBooNE;
}

namespace ublarcvapp {
namespace mctools {

  /**
     \class ProcessBase
     User defined class LArbysImageMC ... these comments are used to generate
     doxygen documentation!
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
    void finalize();
    void Clear();

    void bindAnaVariables( TTree* );
    void printInteractionInfo();
    
  protected:

    std::string _name;
    
    std::string _producer_mctruth;
    std::string _producer_mctrack;
    std::string _producer_mcshower;    

    TTree* _mc_tree;

    larutil::SpaceChargeMicroBooNE* _psce;

  public:
    
    /// Event ID
    int _run;
    int _subrun;
    int _event;
    int _entry;

    std::string _rse_producer;

    // neutrino interaction info
    bool _neutrino_present;    
    int _current_type;
    int _interaction_type;
    int _genie_mode;

    float _Enu_true;
    float _vtx_t;
    float _vtx_x;
    float _vtx_y;
    float _vtx_z;
    float _vtx_sce_x;
    float _vtx_sce_y;
    float _vtx_sce_z;
    float _vtx_tick;
    float _vtx_wire[3];

    float _evis;
    float _evis_had;
    float _evis_vtx;
    float _evis_lep;

    int _hi_lep_pdg;
    float _hi_lep_e;
    
    int _nprimary;
    int _ntotal;

    int _nproton;
    int _nproton_60mev;    
    int _nlepton;
    int _nlepton_35mev;
    int _nmeson;
    int _nmeson_35mev;
    int _npi0;
    int _nshower;
    int _nneutron;
    int _1l1p0pi;

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
