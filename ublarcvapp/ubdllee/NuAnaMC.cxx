#ifndef __NU_ANA_MC_CXX__
#define __NU_ANA_MC_CXX__

#include "NuAnaMC.h"

#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventROI.h"
#include "larcv/core/DataFormat/EventChStatus.h"

#include "larlite/LArUtil/LArProperties.h"
#include "larlite/LArUtil/Geometry.h"
#include "larlite/LArUtil/DetectorProperties.h"
#include "larlite/LArUtil/TimeService.h"
#include "larlite/DataFormat/mctruth.h"
#include "larlite/DataFormat/mctrack.h"
#include "larlite/DataFormat/mcshower.h"

#include "ublarcvapp/ubdllee/dwall.h"

namespace ublarcvapp {
namespace ubdllee {

  static NuAnaMCProcessFactory __global_NuAnaMCProcessFactory__;

  NuAnaMC::NuAnaMC(const std::string name)
    : larcv::ProcessBase(name), _event_tree(nullptr)
  {
  }

  NuAnaMC::~NuAnaMC() {
    delete _sce;
    delete _ll;
  }
    
  void NuAnaMC::configure(const larcv::PSet& cfg)
  {

    // only create own tree if not defined yet.
    // this means user can override tree, by calling addBranchestoTree before configure.
    if ( !_event_tree ) {
      _event_tree = new TTree("NuAnaMCTree","Neutrino Truth Analysis");
      addBranchesToTree( _event_tree );
    }

    _mctruth_producer  = cfg.get< std::string >( "MCTruthProducer",  "generator" );
    _mctrack_producer  = cfg.get< std::string >( "MCTrackProducer",  "mcreco" );
    _mcshower_producer = cfg.get< std::string >( "MCShowerProducer", "mcreco" );
    _chstatus_producer = cfg.get< std::string >( "ChStatusProducer", "wire" );
    _larlite_file_v    = cfg.get< std::vector<std::string> >( "LArliteInputFiles" );
    
  }

  void NuAnaMC::addBranchesToTree( TTree* tree ) {
    
    _event_tree = tree;

    // event info
    tree->Branch("run"      , &_run      , "run/I");
    tree->Branch("subrun"   , &_subrun   , "subrun/I");
    tree->Branch("event"    , &_event    , "event/I");
    tree->Branch("entry"    , &_entry    , "entry/I");

    // interaction info
    tree->Branch("nupdg",   &_nu_pdg,           "nu_pdg/I");
    tree->Branch("iscc",    &_ccnc,             "iscc/I" );
    tree->Branch("mode",    &_genie_mode,       "mode/I");
    tree->Branch("nuanceid",&_nuance_mode,      "nuanceid/I");    
    tree->Branch("is1l1p",  &_is1l1p,           "is1l1p/I");
    tree->Branch("isfv",    &_isfv,             "isfv/I" );
    tree->Branch("inbad",   _inbadch,           "inbad[3]/I");
    tree->Branch("dwall",   &_dwall,            "dwall/F");
    tree->Branch("enu_mev", &_enu,              "enu_mev/F");
    tree->Branch("q2",      &_q2,               "q2/F");
    tree->Branch("x_b",     &_bjorken_x,        "x_b/F");
    tree->Branch("W",       &_W,                "W/F");    
    tree->Branch("vtx",     _vtx,               "vtx[3]/F");
    tree->Branch("vtx_sce", _vtx_sce,           "vtx_sce[3]/F");
    tree->Branch("wireid",  _wid,               "wireid[3]/F");
    tree->Branch("tick",    &_tick,             "tick/F");
    
    
  }

  void NuAnaMC::clearVariables() {

    // -----------------------
    // interaction variables
    // -----------------------
    _nu_pdg = -1;             ///< neutrino pdg; -1 if no neutrino
    _ccnc   = -1;             ///< CC or NC
    _genie_mode  = -1;        ///< GENIE interaction mode
    _nuance_mode = -1;        ///< NUANCE interaction mode    
    _is1l1p = 0;              ///< is 1e1p or 1mu1p + 0 visible gamma + 0 charged pion
    _isfv = 0;                ///< is within 10 cm from TPC wall
    _dwall = 0;               ///< distance to closest TPC wall
    _enu   = 0;               ///< true neutrino energy
    _q2    = -1.0;            ///< q-squared
    _bjorken_x = -1.0;        ///< bjorken-X
    _W     = -1.0;            ///< W
    _tick = -1.0;             ///< tick vertex occurs on in images
    _t    = -1.0;             ///< geant4 time of the event
    for (int i=0; i<3; i++ ) {
      _inbadch[i] = 0;        ///< 1 if in badch on plane
      _vtx[i] = -2000;        ///< x,y,z of true vertex
      _vtx_sce[i] = -2000;    ///< space-charge corrected true vertex location
      _wid[i] = -1;           ///< wires vertex occurred on in images
    }
    
  }

  void NuAnaMC::fillInteractionVariables( const larlite::event_mctruth&  ev_mctruth,
                                          const larlite::event_mctrack&  ev_mctrack,
                                          const larlite::event_mcshower& ev_mcshower,
                                          const larcv::EventChStatus& ev_chstatus ) {

    // loop over interactions
    for ( auto const& mctruth : ev_mctruth) {
      LARCV_DEBUG() << "==== [interaction] ======================================" << std::endl;
      LARCV_DEBUG() << "mctruth origin: " << mctruth.Origin() << std::endl;
      bool nufilled = false;
      for ( auto const& mcpart : mctruth.GetParticles() ) {
	LARCV_DEBUG() << "  mcpart: "
                      << " status=" << mcpart.StatusCode()
                      << " pdg=" << mcpart.PdgCode() 
                      << " E=" << mcpart.Trajectory().front().E() 
                      << " vtx=(" << mcpart.Trajectory().front().X() << "," << mcpart.Trajectory().front().Y() << "," << mcpart.Trajectory().front().Z() << ")"
                      << " t=" << mcpart.Trajectory().front().T()
                      << std::endl;
	if ( !nufilled && mcpart.StatusCode()==0 && 
	     ( abs(mcpart.PdgCode())==12 || abs(mcpart.PdgCode())==14 || abs(mcpart.PdgCode())==16 ) ) {

          // vertex
	  _vtx[0] =  mcpart.Trajectory().front().X();
	  _vtx[1] =  mcpart.Trajectory().front().Y();
	  _vtx[2] =  mcpart.Trajectory().front().Z();
	  _t      =  mcpart.Trajectory().front().T();
          std::vector<double> vtx_v = { (double)_vtx[0], (double)_vtx[1], (double)_vtx[2] };

          //double g4Ticks = detClocks->TPCG4Time2Tick( _t ) + detProperties->GetXTicksOffset(0,0,0) - detProperties->TriggerOffset();
          //double xtimeoffset = detProperties->ConvertTicksToX(g4Ticks,0,0,0);
          double g4ticks = larutil::TimeService::GetME()->TPCG4Time2Tick(_t);
          double xticksoffset = larutil::DetectorProperties::GetME()->GetXTicksOffset(0);
          double trigoffset   =  larutil::DetectorProperties::GetME()->TriggerOffset();
          std::cout << "g4ticks=" << g4ticks << " xticksoffset=" << xticksoffset << " trigoffset=" << trigoffset << std::endl;
          double tick_offset = g4ticks;
          double x_offset = larutil::DetectorProperties::GetME()->ConvertTicksToX( tick_offset, 0 );
          std::cout << "tick_offset=" << tick_offset << " x_offset=" << x_offset << std::endl;
          
          std::vector<double> offsets = _sce->GetPosOffsets( vtx_v[0], vtx_v[1], vtx_v[2] );
          _vtx_sce[0] = _vtx[0] - offsets[0] + 0.6;
          _vtx_sce[1] = _vtx[1] + offsets[1];
          _vtx_sce[2] = _vtx[2] + offsets[2];

          _tick = 3200 + _vtx_sce[0]/larutil::LArProperties::GetME()->DriftVelocity()/0.5;

          // ---------------------------------------
          // hack: 
          _tick += _sce->tickoffset_forward_hack( _tick );
          // ---------------------------------------
          
          for ( size_t p=0; p<3; p++ ) {
            _wid[p] = larutil::Geometry::GetME()->WireCoordinate( _vtx_sce, p );
            auto const& chs = ev_chstatus.Status( (larcv::PlaneID_t)p );
            const std::vector<short>& status_v = chs.as_vector();
            if ( (int)_wid[p]>=0 && (int)_wid[p]<status_v.size() ) {
              if ( status_v[ (int)_wid[p] ]>=4 )
                _inbadch[p] = 0;
              else
                _inbadch[p] = 1;
            }
          }
          
          int boundary_type;
          _dwall = ublarcvapp::dwall( vtx_v, boundary_type );
          
	  // input neutrino
	  _nu_pdg = mcpart.PdgCode();
	  _ccnc = mctruth.GetNeutrino().CCNC();
	  _genie_mode = mctruth.GetNeutrino().Mode();
	  _nuance_mode = mctruth.GetNeutrino().InteractionType();
          _isfv = ( _dwall>0 && _dwall<10.0 ) ? 1 : 0;
	  _enu    = mcpart.Trajectory().front().E()*1000.0; // save in mev
	  _q2   = mctruth.GetNeutrino().QSqr();
          _bjorken_x = mctruth.GetNeutrino().X();
          _W    = mctruth.GetNeutrino().W();

	  nufilled = true;
          
	}
      }
    
  }

  /*

    for ( auto const& mctruth : *ev_mctruth) {
      std::cout << "mctruth origin: " << mctruth.Origin() << std::endl;
      for ( auto const& mcpart : mctruth.GetParticles() ) {
	std::cout << "  mcpart: "
		  << " status=" << mcpart.StatusCode()
		  << " pdg=" << mcpart.PdgCode() 
		  << " E=" << mcpart.Trajectory().front().E() 
		  << " vtx=(" << mcpart.Trajectory().front().X() << "," << mcpart.Trajectory().front().Y() << "," << mcpart.Trajectory().front().Z() << ")"
		  << " t=" << mcpart.Trajectory().front().T()
		  << std::endl;
	if ( !nufilled && mcpart.StatusCode()==0 && 
	     ( abs(mcpart.PdgCode())==12 || abs(mcpart.PdgCode())==14 || abs(mcpart.PdgCode())==16 ) ) {
	  // input neutrino
	  nupdg = mcpart.PdgCode();
	  nuE = mcpart.Trajectory().front().E()*1000.0; // save in mev
	  x =  mcpart.Trajectory().front().X();
	  y =  mcpart.Trajectory().front().Y();
	  z =  mcpart.Trajectory().front().Z();
	  t =  mcpart.Trajectory().front().T();

	  ccnc = mctruth.GetNeutrino().CCNC();
	  mode = mctruth.GetNeutrino().Mode();
	  nuancecode = mctruth.GetNeutrino().InteractionType();
	  q2   = mctruth.GetNeutrino().QSqr();
	  nufilled = true;
	}

	float ke = (mcpart.Trajectory().front().E() - mcpart.Mass())*1000.0;
	if ( mcpart.StatusCode()==1 ) {
	  // final state particles
	  float ke = (mcpart.Trajectory().front().E() - mcpart.Mass())*1000.0;
	  if ( abs(mcpart.PdgCode())==11 ) {
	    nlepton++;
	    if ( ke>35.0 )
	      nlepton_abovethresh++;
	    if ( ke>maxleptonKEmev ) {
	      maxleptonKEmev = ke;
	    }

	  }
	  else if ( abs(mcpart.PdgCode())==2212 ) {
	    nproton++;
	    if ( ke>35.0 )
	      nproton_abovethresh++;
	    if ( ke>maxprotonKEmev ) {
	      maxprotonKEmev = ke;
	    }

	  }
	  else {
	    if(ke>35.0)
	      nother_abovethresh++;
	  }
	}//end of if final state
      }//end of particle loop
    }//end of mctruth loop

    if ( nlepton_abovethresh==1 && nproton_abovethresh==1 && nother_abovethresh==0)
      l1p1 = 1;
    else
      l1p1 = 0;


  
  //from rui an.
  bool NuAnaMC::MCSelect(const EventROI* ev_roi) {

    bool pdg_b = false;
    bool interaction_mode_b = false;
    bool energy_init_b = false;
    bool vis_lepton_b = false;
    bool vis_one_proton_b = false;
    bool vis_no_unknowns_b = false;
    
    int nlepton  = 0;
    int nproton  = 0;
    
    int unknown_ctr = 0.0;
    
    // Get neutrino ROI
    const auto& roi = ev_roi->at(0);
    
    int pdgcode          = roi.PdgCode();
    int interaction_mode = roi.NuInteractionType();
    double energy_init   = roi.EnergyInit();

    std::vector<aparticle> protons_v;
    std::vector<aparticle> leptons_v;
    std::vector<aparticle> unknown_v;

    bool first = false;
    
    for(const auto& roi : ev_roi->ROIArray()) {

      if(!first) { first=true; continue; }
      
      int pdgcode = std::abs(roi.PdgCode());
      
      if (pdgcode==2212) {
	aparticle thispro;
	thispro.pdg           = pdgcode;
	thispro.trackid       = roi.TrackID();
	thispro.ptrackid      = roi.ParentTrackID();
	thispro.depeng        = roi.EnergyDeposit();
	thispro.primary       = ( roi.TrackID() == roi.ParentTrackID() );
	protons_v.push_back(std::move(thispro));
      }
      
      else if (pdgcode==11 or pdgcode==13) {
	aparticle thislep;
	thislep.pdg           = pdgcode;
	thislep.trackid       = roi.TrackID();
	thislep.ptrackid      = roi.ParentTrackID();
	thislep.depeng        = roi.EnergyDeposit();
	thislep.primary       = ( roi.TrackID() == roi.ParentTrackID() );
	leptons_v.push_back(std::move(thislep));
      }

      else {
	aparticle thisunknown;
	thisunknown.pdg           = pdgcode;
	thisunknown.trackid       = roi.TrackID();
	thisunknown.ptrackid      = roi.ParentTrackID();
	thisunknown.depeng        = roi.EnergyDeposit();
	thisunknown.primary       = ( roi.TrackID() == roi.ParentTrackID() );
	unknown_v.push_back(std::move(thisunknown));
      }
      
    }
    
    //calculate the visible lepton energy
    int lepton_engs_ctr = 0 ;
    for (size_t p1=0;p1 < leptons_v.size() ; p1++ ){
      const auto& lepton1 = leptons_v[p1];
      if (! lepton1.primary ) continue;
      nlepton+=1;
      float this_lepton_eng = lepton1.depeng;
      for (size_t p2=0;p2 < leptons_v.size() ; p2++ ){
      	if (p1==p2) continue;
      	const auto& lepton2 = leptons_v[p2];
      	if (lepton2.ptrackid != lepton1.trackid) continue;
      	this_lepton_eng+=lepton2.depeng;
      }
      if ( this_lepton_eng > _dep_sum_lepton ) lepton_engs_ctr ++;
    }

    //calculate the visible proton energy
    int proton_engs_ctr = 0 ;
    for (size_t p1=0;p1 < protons_v.size() ; p1++ ){
      const auto& proton1 = protons_v[p1];
      if (! proton1.primary ) continue;
      nproton+=1;
      float this_proton_eng = proton1.depeng;
      for (size_t p2=0;p2 < protons_v.size() ; p2++ ){
	if (p1==p2) continue;
	const auto& proton2 = protons_v[p2];
	if (proton2.ptrackid != proton1.trackid) continue;
	this_proton_eng+=proton2.depeng;
      }
      if ( this_proton_eng > _dep_sum_proton ) proton_engs_ctr ++;
    }


    for(const auto& unknown : unknown_v) {
      if (!unknown.primary) continue;
      unknown_ctr+=1;
    }

    
    // requested nu PDG code
    if (pdgcode == _nu_pdg or _nu_pdg == 0) {
      pdg_b = true;
    }

    // requested nu interaction mode
    if ((interaction_mode == _genie_mode) or (_genie_mode == 0)) {
      interaction_mode_b = true;
    } 

    // neutrino energy range
    if (energy_init >= _min_nu_init_e && energy_init <= _max_nu_init_e) {
      energy_init_b = true;
    }

    // interaction must have visible leptons
    if (lepton_engs_ctr == 1) {
      vis_lepton_b = true;
    }

    // There must be 1 visible proton
    if (proton_engs_ctr == 1) {
      vis_one_proton_b = true;
    }

    // There must be no primary unknowns
    if (unknown_ctr==0) {
      vis_no_unknowns_b = true;
    }
    
    bool selected =
      pdg_b               &&
      interaction_mode_b  &&
      energy_init_b       &&
      vis_lepton_b        &&
      vis_one_proton_b    &&
      vis_no_unknowns_b;
    
    if (! pdg_b)              _n_fail_nupdg      += 1;
    if (! interaction_mode_b) _n_fail_inter      += 1;
    if (! energy_init_b)      _n_fail_nuE        += 1;
    if (! vis_lepton_b)       _n_fail_lepton_dep += 1;
    if (! vis_one_proton_b)   _n_fail_proton_dep += 1;
    if (! vis_no_unknowns_b)  _n_fail_unknowns   += 1;
    
    return selected;
  }
  */
  }
  
  void NuAnaMC::initialize()
  {
    _ll = new ublarcvapp::LArliteManager( larlite::storage_manager::kREAD, "NuAnaMCLarlite" );
    for ( auto const& input_larlite : _larlite_file_v ) {
      _ll->add_in_filename( input_larlite );
    }
    _ll->open();

    _sce = new larutil::SpaceChargeMicroBooNE( larutil::SpaceChargeMicroBooNE::kMCC9_Forward, "" );
  }

  bool NuAnaMC::process( larcv::IOManager& mgr )
  {
    LARCV_DEBUG() << "start" << std::endl;

    auto chstatus_v = (larcv::EventChStatus*)   mgr.get_data( larcv::kProductChStatus,  _chstatus_producer );    
    _ll->syncEntry(mgr);
    
    const int kINVALID_INT=-1;
    _run      = kINVALID_INT;
    _subrun   = kINVALID_INT;
    _event    = kINVALID_INT;
    _entry    = kINVALID_INT;

    auto mctruth_v  = (larlite::event_mctruth*) _ll->get_data(larlite::data::kMCTruth,  _mctruth_producer );
    auto mctrack_v  = (larlite::event_mctrack*) _ll->get_data(larlite::data::kMCTrack,  _mctrack_producer );
    auto mcshower_v = (larlite::event_mcshower*)_ll->get_data(larlite::data::kMCShower, _mcshower_producer );

    _run    = _ll->run_id();
    _subrun = _ll->subrun_id();
    _event  = _ll->event_id();
    _entry  = _ll->get_index();

    fillInteractionVariables( *mctruth_v, *mctrack_v, *mcshower_v, *chstatus_v );

    _event_tree->Fill();

    return true;
  }
  
  void NuAnaMC::finalize()
  {
    _ll->close();
    _event_tree->Write();
  }

}
}
#endif
