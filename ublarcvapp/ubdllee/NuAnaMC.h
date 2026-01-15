/**
 * \file NuAnaMC.h
 *
 * \ingroup Filter
 * 
 * \brief Defines 1 lepton 1 proton signal events using larlite info.
 *
 * @author twongjirad
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __NU_ANA_MC_H__
#define __NU_ANA_MC_H__

#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"
#include "larcv/core/DataFormat/EventROI.h"
#include "larcv/core/DataFormat/EventChStatus.h"
#include "larcv/core/DataFormat/IOManager.h"

#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"
#include "TTree.h"

#include "ublarcvapp/LArliteHandler/LArliteManager.h"

namespace ublarcvapp {
namespace ubdllee {
  
  /**
     \class NuAnaMC
  */
  class NuAnaMC : public larcv::ProcessBase {

  public:
    
    /// Default constructor
    NuAnaMC(const std::string name="NuAnaMC");
    
    /// Default destructor
    virtual ~NuAnaMC();

    void configure(const larcv::PSet&);
    void initialize();
    bool process(larcv::IOManager& mgr);
    void finalize();

    void addBranchesToTree( TTree* ttree );
    void setLarliteInputFile( const std::string inputfilename );
    void clearVariables();
    
  protected:

    // configuration parameters
    std::string _mctruth_producer;
    std::string _mctrack_producer;
    std::string _mcshower_producer;
    std::string _chstatus_producer;

    // larlite iomanager
    ublarcvapp::LArliteManager* _ll;
    std::vector<std::string>    _larlite_file_v;

    // output ana tree
    TTree* _event_tree;

    // space charge
    larutil::SpaceChargeMicroBooNE* _sce;


  public:
    
    // analysis variables
    // -------------------

    // event info
    int _run;                 ///< run number
    int _subrun;              ///< subrun number
    int _event;               ///< event
    int _entry;               ///< entry (in file)

    // interaction info
    int _nu_pdg;              ///< neutrino pdg; -1 if no neutrino
    int _ccnc;                ///< CC or NC
    int _genie_mode;          ///< GENIE interaction mode
    int _nuance_mode;         ///< NUANCE Mode
    int _is1l1p;              ///< is 1e1p or 1mu1p + 0 visible gamma + 0 charged pion
    int _isfv;                ///< is within 10 cm from TPC wall
    int _inbadch[3];          ///< 1 if in badch on plane
    float _dwall;             ///< distance to closest TPC wall
    float _enu;               ///< true neutrino energy
    float _q2;                ///< q-squared
    float _bjorken_x;         ///< bjorken-X
    float _W;                 ///< W: energy transfer to lepton
    float _t;                 ///< geant4 time of event
    float _vtx[3];            ///< x,y,z of true vertex
    float _vtx_sce[3];        ///< space-charge corrected true vertex location
    float _wid[3];            ///< wires vertex occurred on in images
    float _tick;              ///< tick vertex occurs on in images

    
    /* float _dep_sum_lepton;    ///< total deposited energy by the primary lepton */
    /* float _dep_sum_proton;    ///< total deposited energy by the primary protons */
    /* float _evis;              ///< total visible energy: charged leptons, charged pions, gamma showers */
    /* float _evis_had;          ///< total visible hadronic energy: all the above, no lepton */
    /* float _evis_vertex;       ///< total visible energy for below threshold particles at vertex */
    /* float _adctot_vertex[3];  ///< total pixel ADC sum around true vertex */
    /* float _neutron_energy;    ///< total energy carried off by neutrons */
    /* float _min_gamma_gap_cm;  ///< shortest distance from visible shower start to vertex */

    /* // particle counts: at the vertex */
    /* int    _nproton_60mev;    ///< number of protons above 60 MeV threshold */
    /* int    _nproton;          ///< number of total final state protons */
    /* int    _npion_35mev;      ///< number of charged pions above 35 MeV threshold */
    /* int    _npion;            ///< number of pions */
    /* int    _nneutron;         ///< number of neutrons from vertex */
    /* int    _ngamma;           ///< number of visible gammas */
    
    /* struct aparticle{ */
    /*   int pdg; */
    /*   int trackid; */
    /*   int ptrackid; */
    /*   bool primary; */
    /*   float depeng; */
    /* }; */

  protected:
    
    void fillInteractionVariables( const larlite::event_mctruth&  ev_mctruth,
                                   const larlite::event_mctrack&  ev_mctrack,
                                   const larlite::event_mcshower& ev_mcshower,
                                   const larcv::EventChStatus& ev_chstatus );

  };

  /**
     \class larcv::NuAnaMCFactory
     \brief A concrete factory class for larcv::NuAnaMC
  */
  class NuAnaMCProcessFactory : public larcv::ProcessFactoryBase {
  public:
    /// ctor
    NuAnaMCProcessFactory() { larcv::ProcessFactory::get().add_factory("NuAnaMC",this); }
    /// dtor
    ~NuAnaMCProcessFactory() {}
    /// creation method
    larcv::ProcessBase* create(const std::string instance_name) {
      return new ublarcvapp::ubdllee::NuAnaMC(instance_name);
    };
    
  };

}
}
#endif
/** @} */ // end of doxygen group 

