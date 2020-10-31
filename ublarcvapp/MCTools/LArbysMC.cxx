#ifndef __LARLITECV_LARBYSMC_CXX__
#define __LARLITECV_LARBYSMC_CXX__

#include <algorithm>

#include "LArbysMC.h"

#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"
#include "Base/MCConstants.h"
#include "DataFormat/storage_manager.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctruth.h"

namespace ublarcvapp {
namespace mctools {
  
  //static LArbysImageMCProcessFactory __global_LArbysImageMCProcessFactory__;

  LArbysMC::LArbysMC()
    : _name("LArbysMC"),
      _mc_tree(nullptr)
  {
    _producer_mctruth  = "generator";
    _producer_mctrack  = "mcreco";
    _producer_mcshower = "mcreco";
    _producer_wireimage = "wire";
    _psce = new larutil::SpaceChargeMicroBooNE();
    for (int p=0; p<3; p++) _plane_vtx_pixsum[p] = 0.;
      
    
  }

  LArbysMC::~LArbysMC()
  {
    if (_mc_tree)
      delete _mc_tree;
    _mc_tree = nullptr;
  }

  // void LArbysMC::configure(const PSet& cfg)
  // {
  //   // _rse_producer       = cfg.get<std::string>("RSEProducer");
  //   // _producer_roi       = cfg.get<std::string>("MCProducer");
  //   // _producer_image2d   = cfg.get<std::string>("Image2DProducer");
  //   // _neutrino_present   = cfg.get<bool>("NeutrinoPresent",true);
  // }

  /**
   * initialize interal storage tree
   *
   */
  void LArbysMC::initialize()
  {
    _mc_tree = new TTree("LArbysMCTree","MC infomation");
    bindAnaVariables(_mc_tree);
  }

  /**
   * @brief bind variables to the TTree we will fill to
   *
   * @param[in] mc_tree ROOT TTree to which we add branches
   *
   */
  void LArbysMC::bindAnaVariables( TTree* mc_tree )
  {

    // event indexing
    mc_tree->Branch("larbysmc_run",     &_run,     "larbysmc_run/I");
    mc_tree->Branch("larbysmc_subrun",  &_subrun,  "larbysmc_subrun/I");
    mc_tree->Branch("larbysmc_event",   &_event,   "larbysmc_event/I");
    mc_tree->Branch("larbysmc_entry"  , &_entry,   "larbysmc_entry/I");

    mc_tree->Branch("nuPresent",       &_neutrino_present, "nuPresent/O" );
    mc_tree->Branch("currentType",     &_current_type,     "currentType/I");
    mc_tree->Branch("interactionType", &_interaction_type, "interactionType/I");
    mc_tree->Branch("genieMode",       &_genie_mode,       "genieMode/I");
    mc_tree->Branch("Enu_true",        &_Enu_true,         "Enu_true/F");
    mc_tree->Branch("vtx_x",           &_vtx_x,            "vtx_x/F");
    mc_tree->Branch("vtx_y",           &_vtx_y,            "vtx_y/F");
    mc_tree->Branch("vtx_z",           &_vtx_z,            "vtx_z/F");        
    mc_tree->Branch("vtx_sce_x",       &_vtx_sce_x,        "vtx_sce_x/F");
    mc_tree->Branch("vtx_sce_y",       &_vtx_sce_y,        "vtx_sce_y/F");
    mc_tree->Branch("vtx_sce_z",       &_vtx_sce_z,        "vtx_sce_z/F");
    mc_tree->Branch("vtx_tick",        &_vtx_tick,         "vtx_tick/F");
    mc_tree->Branch("vtx_wire",        _vtx_wire,          "vtx_wire[3]/F");
    mc_tree->Branch("vtx_pixsum",      _plane_vtx_pixsum,  "vtx_pixsum[3]/F");
    mc_tree->Branch("vtx_med_pixsum",  &_vtx_med_pixsum,   "vtx_med_pixsum/F");
    
    mc_tree->Branch("evis",            &_evis,             "evis/F");
    mc_tree->Branch("evis_had",        &_evis_had,         "evis_had/F");
    mc_tree->Branch("evis_vtx",        &_evis_vtx,         "evis_vtx/F");
    mc_tree->Branch("evis_lep",        &_evis_lep,         "evis_lep/F");

    mc_tree->Branch("hi_lep_pdg",      &_hi_lep_pdg,       "hi_lep_pdg/I");
    
    mc_tree->Branch("nprimary", &_nprimary,"nprimary/I");

    mc_tree->Branch("nproton", &_nproton,"nproton/I");
    mc_tree->Branch("nproton_60mev", &_nproton_60mev,"nproton_60mev/I");    
    mc_tree->Branch("nlepton", &_nlepton,"nlepton/I");
    mc_tree->Branch("nlepton_35mev", &_nlepton_35mev,"nlepton_35mev/I");    
    mc_tree->Branch("nmeson", &_nmeson,"nmeson/I");
    mc_tree->Branch("nmeson_35mev", &_nmeson_35mev,"nmeson_35mev/I");
    mc_tree->Branch("nshower", &_nshower,"nshower/I");
    mc_tree->Branch("nneutron", &_nneutron,"nneutron/I");
    mc_tree->Branch("npi0",     &_npi0, "npi0/I" );
    mc_tree->Branch("is1l1p0pi", &_1l1p0pi, "is1l1p0pi/I");
    mc_tree->Branch("is1l0p0pi", &_1l0p0pi, "is1l0p0pi/I");

  }

  bool LArbysMC::process(larlite::storage_manager& mgr)
  {
    
    Clear();

    larlite::event_mctruth* ev_mctruth =
      (larlite::event_mctruth*)mgr.get_data(larlite::data::kMCTruth,_producer_mctruth);

    larlite::event_mctrack* ev_mctrack =
      (larlite::event_mctrack*)mgr.get_data(larlite::data::kMCTrack,_producer_mctrack);

    larlite::event_mcshower* ev_mcshower =
      (larlite::event_mcshower*)mgr.get_data(larlite::data::kMCShower,_producer_mcshower);
    
    _run    = mgr.run_id();
    _subrun = mgr.subrun_id();
    _event  = mgr.event_id();
    _entry  = (int)mgr.get_index();

    // If we've got a neutrino, sample that
    auto const& mct = ev_mctruth->at(0);
    _neutrino_present = mct.Origin()==larlite::simb::kBeamNeutrino;
    
    // way to fill an empty entry?
    if ( !_neutrino_present ) {
      if ( _mc_tree )
        _mc_tree->Fill();
      return false;
    }
        
    auto const& nu = mct.GetNeutrino();

    _current_type     = nu.CCNC();
    _interaction_type = nu.InteractionType();
    _genie_mode       = nu.Mode();
    _Enu_true         = nu.Nu().Trajectory().front().E()*1000.0;
    _vtx_t            = nu.Nu().Trajectory().front().T();
    _vtx_x            = nu.Nu().Trajectory().front().X();
    _vtx_y            = nu.Nu().Trajectory().front().Y();
    _vtx_z            = nu.Nu().Trajectory().front().Z();    

    // SCE correction
    std::vector<double> pos_offset = _psce->GetPosOffsets( _vtx_x, _vtx_y, _vtx_z );
    _vtx_sce_x = _vtx_x - pos_offset[0] + 0.7;
    _vtx_sce_y = _vtx_y + pos_offset[1];
    _vtx_sce_z = _vtx_z + pos_offset[2];
    
    const float trig_time = 4050.0;
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    _vtx_tick = ( _vtx_t*1.0e-3 - (trig_time-4050.0) )/0.5 + _vtx_sce_x/cm_per_tick + 3200.0;

    if ( _vtx_y<-117.0 || _vtx_y>117.0
         || _vtx_z<0 || _vtx_z>1036.0 ) {
      // out of tpc
      for (int p=0; p<3; p++ ) _vtx_wire[p] = -1;
    }
    else {
      // in TPC
      double dpos[3] = { _vtx_sce_x, _vtx_sce_y, _vtx_sce_z };
      for (int p=0;p<3; p++) {
        _vtx_wire[p] = larutil::Geometry::GetME()->WireCoordinate( dpos, p );
      }
    }

    for (int p=0; p<3; p++) _plane_vtx_pixsum[p] = 0.;
    _vtx_med_pixsum = 0.;
        
    // final state tally
    _nprimary = 0;
    _ntotal   = 0;
    _nproton  = 0;
    _nproton_60mev  = 0;    
    _nneutron = 0;
    _nlepton  = 0;
    _nlepton_35mev  = 0;    
    _nmeson   = 0;
    _nmeson_35mev   = 0;    
    _nshower  = 0;
    _npi0 = 0;
    _1l1p0pi = 0;
    _1l0p0pi = 0;    
    
    _evis = 0;
    _evis_had = 0;
    _evis_vtx = 0;
    _evis_lep = 0;

    _hi_lep_pdg      = -1;
    _hi_lep_e        = 0;


    int vidx = -1;
    for(auto const& mct : *ev_mctrack ) {
      vidx++;
      if ( mct.Origin()==1 ) {
        // neutrino interaction primary
        int pid  = mct.PdgCode();

        std::cout << "track[" << vidx << "] origin=" << mct.Origin()
                  << " tid=" << mct.TrackID()
                  << " [ mid=" << mct.MotherTrackID()
                  << " mpdg=" << mct.MotherPdgCode() << " ]"
                  << " [ aid=" << mct.AncestorTrackID()
                  << " apdg=" << mct.AncestorPdgCode() << "]"
                  << " pid=" << mct.PdgCode()
                  << " E=" << mct.Start().E()          
                  << std::endl;
        
        
        if ( abs(pid)==12 || abs(pid)==14 || abs(pid)==16 )
          continue; // neutrino: skip

        if ( mct.TrackID()!=mct.MotherTrackID() )
          continue; // skip secondaries
        
        _nprimary++;
        
        float pE = mct.Start().E();
        if ( abs(pid)==13 ){
          // muon
          float ke=pE-105.0;
          if ( ke>35.0 ) _nlepton_35mev++;
          _nlepton++;
          _evis_lep += ke;
          _evis += ke;
          if ( ke>_hi_lep_e ) {
            _hi_lep_e = ke;
            _hi_lep_pdg =pid;
          }
        }        
        else if ( abs(pid)==11 ) {
          // electron
          float ke = pE-0.511;
          if ( ke>35.0 ) _nlepton_35mev++;
          _nlepton++;
          _evis += ke;
          _evis_lep += ke;
          if ( ke>_hi_lep_e ) {
            _hi_lep_e = ke;
            _hi_lep_pdg =pid;
          }          
        }
        else if ( pid==2212 ) {
          // proton
          float ke = pE-938.0;
          if( ke>60.0 ) {
            _nproton_60mev++;
          }
          else {
            _evis_vtx += ke;
          }
          _nproton++;
          _evis += ke;
          _evis_had += ke;
        }
        else if ( abs(pid)==211 ) {
          float ke = pE-135.0;
          if ( ke>35.0 ) {
            _nmeson_35mev++;
          }
          else {
            _evis_vtx += ke;
          }
          _nmeson++;          
          _evis_had += ke;
          _evis += ke;
        }
        else if ( pid==2112 ) {
          _nneutron++;
        }
        
      }
    }

    vidx = -1;
    std::set<int> pi0_ids; // list of pi0 id's we've already used
    for(auto const& mcs : *ev_mcshower ) {
      vidx++;
      if ( mcs.Origin()==1 ) {
        // neutrino interaction primary
        int pid  = mcs.PdgCode();

        std::cout << "shower[" << vidx << "] origin=" << mcs.Origin()
                  << " tid=" << mcs.TrackID()
                  << " [ mid=" << mcs.MotherTrackID()
                  << " mpd=" << mcs.MotherPdgCode() << " ]"
                  << " [ aid=" << mcs.AncestorTrackID()
                  << " apdg=" << mcs.AncestorPdgCode() << " ]"
                  << " pid=" << mcs.PdgCode()
                  << " E=" << mcs.Start().E()
                  << std::endl;
        
        if ( abs(pid)==12 || abs(pid)==14 || abs(pid)==16 )
          continue; // neutrino: skip

        if ( abs(pid)!=22 && mcs.TrackID()!=mcs.MotherTrackID() )
          continue; // skip secondaries (non photons)
        if ( abs(pid)==22 && mcs.MotherTrackID()!=mcs.AncestorTrackID() )
          continue; // skip secondaries (photons)
        
        float pE = mcs.Start().E();
        if ( abs(pid)==11 ) {
          // electron
          float ke = pE-0.511;
          if ( ke>35.0 ) _nlepton_35mev++;
          _nlepton++;
          _nshower++;          
          _evis += ke;
          _evis_lep += ke;
          if ( ke>_hi_lep_e ) {
            _hi_lep_e = ke;
            _hi_lep_pdg =pid;
          }          
        }
        else if ( abs(pid)==22 ) {

          if ( mcs.MotherPdgCode()==111 ) {
            // this is a pi0 photon
            int mid = mcs.MotherTrackID();
            if (  pi0_ids.find(mid)==pi0_ids.end() ) {
              // new id
              pi0_ids.insert(mid);
              _npi0++;
              _nmeson++;
            }
          }
          
          // photon
          float ke = pE;
          _nshower++;
          _evis += ke;
        }
        
      }
    }

    if ( _nproton_60mev==1 && _nlepton_35mev==1 && _nmeson_35mev==0 && _npi0==0 )
      _1l1p0pi = 1;
    if ( _nproton_60mev==0 && _nlepton_35mev==1 && _nmeson_35mev==0 && _npi0==0 )
      _1l0p0pi = 1;
    

    if ( _mc_tree )
      _mc_tree->Fill();

    return true;
  }

  /**
   * @brief calculate MC info that includes image-based metrics
   *
   * @param[inout] iolcv LArCV data manager
   * @param[inout] mgr   larlite data manager
   */
  bool LArbysMC::process(larcv::IOManager& iolcv, larlite::storage_manager& mgr)
  {
    
    process( mgr );
    std::vector<float> plane_vtx_pixsum(3,0);
    calculateVertexPixsum( iolcv, _vtx_tick, _vtx_wire[0], _vtx_wire[1], _vtx_wire[2], 3, plane_vtx_pixsum );
    for (int p=0; p<3; p++)
      _plane_vtx_pixsum[p] = plane_vtx_pixsum[p];
    // sort plane sum to get median
    struct plane_pixsum_t {
      int plane;
      float pixsum;
      bool operator<(const plane_pixsum_t& rhs ) {
        if ( pixsum<rhs.pixsum ) return true;
        return false;
      };
    };
    std::vector< plane_pixsum_t > ranked_pixsum_v;
    for (int p=0; p<3; p++) {
      plane_pixsum_t pp;
      pp.plane = p;
      pp.pixsum = _plane_vtx_pixsum[p];
      ranked_pixsum_v.push_back( pp );
    }
    std::sort( ranked_pixsum_v.begin(), ranked_pixsum_v.end() );
    _vtx_med_pixsum = ranked_pixsum_v[1].pixsum;
    
    return true;
  }

  void LArbysMC::calculateVertexPixsum( larcv::IOManager& iolcv,
                                        float vtx_tick, float vtx_wireu, float vtx_wirev, float vtx_wirey,
                                        int pix_radius,
                                        std::vector<float>& plane_vtx_pixsum  )
  {

    larcv::EventImage2D* ev_adc =
      (larcv::EventImage2D*)iolcv.get_data(larcv::kProductImage2D,_producer_wireimage);

    auto const& adc_v = ev_adc->as_vector();
    int vtx_row = adc_v.front().meta().row( vtx_tick );
    int row_max = (int)adc_v.front().meta().rows();
    int col_max = (int)adc_v.front().meta().cols();
    std::vector<int> vtx_col(adc_v.size(),0);
    vtx_col[0] = adc_v[0].meta().col( vtx_wireu );
    vtx_col[1] = adc_v[1].meta().col( vtx_wirev );
    vtx_col[2] = adc_v[2].meta().col( vtx_wirey );

    plane_vtx_pixsum.resize(3,0);
    for (int p=0; p<3; p++)
      plane_vtx_pixsum[p] = 0.;
    
    int rad = abs(pix_radius);
    for (int r=-rad; r<=rad; r++) {
      int row = vtx_row+r;
      if ( row<0 || row>=row_max )continue;
      for (int c=-rad; c<=rad; c++) {
        for (int p=0; p<3; p++) {
          int col = vtx_col[p]+c;
          if ( col<0 || col>=col_max ) continue;
          float pixval = adc_v[p].pixel( row, col );
          if ( pixval>10.0 )
            plane_vtx_pixsum[p] += pixval;
        }//end of plane loop
      }//end of c loop
    }//end of r loop
    
  }


  void LArbysMC::finalize()
  {
    if ( _mc_tree )
      _mc_tree->Write();
  }

  void LArbysMC::Clear() {

    _neutrino_present = false;
    _current_type = -1;
    _interaction_type = -1;
    _genie_mode = -1;
    
  }


  void LArbysMC::printInteractionInfo() {
    std::cout << "[LArbysMC::printInteractionInfo] ==============================================" << std::endl;
    std::cout << " rse: (" << _run << "," << _subrun << "," << _event << ") entry=" << _entry << std::endl;
    std::cout << " True Nu E: " << _Enu_true << " GeV" << std::endl;
    std::cout << " genie mode: " << _genie_mode << std::endl;
    std::cout << " is CC: " << (_current_type==0) << std::endl;
    std::cout << " interaction type: " << _interaction_type << std::endl;
    std::cout << " nlepton: " << _nlepton << "; above 35 MeV " << _nlepton_35mev << std::endl;
    std::cout << " nproton: " << _nproton << "; above 60 MeV " << _nproton_60mev << std::endl;
    std::cout << " nmeson: "  << _nmeson  << "; above 35 MeV " << _nmeson_35mev << std::endl;
    std::cout << " nneutron: " << _nneutron << std::endl;
    std::cout << " npi0: " << _npi0 << std::endl;
    std::cout << " is 1L1P (+ 0 pi + 0 pi0): " << _1l1p0pi << std::endl;
    std::cout << " evis: " << _evis << "; lepton: " << _evis_lep << "; hadronic: " << _evis_had << "; vertex " << _evis_vtx << std::endl;
    std::cout << " (x,y,z) true: (" << _vtx_x << "," << _vtx_y << "," << _vtx_z << ")" << std::endl;
    std::cout << " (x,y,z) sce: (" << _vtx_sce_x << "," << _vtx_sce_y << "," << _vtx_sce_z << ")" << std::endl;
    std::cout << " tick: " << _vtx_tick << std::endl;
    std::cout << " wires: (" << _vtx_wire[0] << "," << _vtx_wire[1] << "," << _vtx_wire[2] << ")" << std::endl;
    std::cout << "==============================================================================" << std::endl;    
  }
  
}
}
#endif
