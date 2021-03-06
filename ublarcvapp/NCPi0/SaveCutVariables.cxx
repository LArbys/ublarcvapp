#include "SaveCutVariables.h"

namespace ublarcvapp {
namespace ncpi0{

void SaveCutVariables::configure( const larcv::PSet& pset ) {
    // all inputs needed to get cut variables
    // larcv truth inputs
    _input_adc_producer          = pset.get<std::string>("InputADCProducer");
    _input_partroi_producer      = pset.get<std::string>("InputPartROIProducer");
    _input_segment_producer      = pset.get<std::string>("InputSegmentProducer");
    _input_instance_producer     = pset.get<std::string>("InputInstanceProducer");
    //ssnet inputs
    _input_ssnet_uplane_producer = pset.get<std::string>("InputSSNetUProducer");
    _input_ssnet_vplane_producer = pset.get<std::string>("InputSSNetVProducer");
    _input_ssnet_yplane_producer = pset.get<std::string>("InputSSNetYProducer");
    // larlite (mcinfo) inputs
    _input_mctrack_producer      = pset.get<std::string>("InputMCTrackProducer");
    _input_mcshower_producer     = pset.get<std::string>("InputMCShowerProducer");
    _input_flux_producer         = pset.get<std::string>("InputFluxProducer");
    _input_mctruth_producer      = pset.get<std::string>("InputMCTruthProducer");
    //larlite (tracker) inputs
    _input_recotrack_producer    = pset.get<std::string>("InputRecoTrackProducer");
    _input_vtxtracker_producer   = pset.get<std::string>("InputVtxTrackerProducer");

  }

  void SaveCutVariables::initialize() {
    OutFile = new TFile("NCPi0CutVariables.root","RECREATE");
    setupAnaTree();
    fin = new TFile("/cluster/tufts/wongjiradlab/kmason03/testdir/ubdl/ublarcvapp/ublarcvapp/NCPi0/test/hadd_NCPi0Probabilities.root","read");
    LoadProbHists();
  }

bool SaveCutVariables::process( larcv::IOManager& io,larlite::storage_manager& ioll,larcv::IOManager& ioforward ){
  //set utils object
  Utils Utils;
  SaveProbabilities Probs;
  // inputs
  const auto ev_img            = (larcv::EventImage2D*)io.get_data( larcv::kProductImage2D, _input_adc_producer );
  const auto ev_instance       = (larcv::EventImage2D*)io.get_data( larcv::kProductImage2D, _input_instance_producer );
  const auto ev_segment        = (larcv::EventImage2D*)io.get_data( larcv::kProductImage2D, _input_segment_producer );
  const auto ev_partroi        = (larcv::EventROI*)(io.get_data( larcv::kProductROI,_input_partroi_producer));
  const auto& ev_mctrack       = *((larlite::event_mctrack*)ioll.get_data(larlite::data::kMCTrack,  _input_mctrack_producer ));
  const auto& ev_mcshower      = *((larlite::event_mcshower*)ioll.get_data(larlite::data::kMCShower,  _input_mcshower_producer));
  const auto ev_flux           = ((larlite::event_mcflux*)ioll.get_data(larlite::data::kMCFlux,  _input_flux_producer ));
  const auto ev_mctruth        = (larlite::event_mctruth*)ioll.get_data(larlite::data::kMCTruth,  _input_mctruth_producer );
  const auto ev_recotrack      = (larlite::event_track*)ioll.get_data(larlite::data::kTrack,  _input_recotrack_producer );
  const auto ev_vtxtracker     = (larlite::event_vertex*)ioll.get_data(larlite::data::kVertex,  _input_vtxtracker_producer );
  const auto ssnet_uplane      = (larcv::EventImage2D*)ioforward.get_data( larcv::kProductImage2D, _input_ssnet_uplane_producer );
  const auto ssnet_vplane      = (larcv::EventImage2D*)ioforward.get_data( larcv::kProductImage2D, _input_ssnet_vplane_producer );
  const auto ssnet_yplane      = (larcv::EventImage2D*)ioforward.get_data( larcv::kProductImage2D, _input_ssnet_yplane_producer );
  run    = ioll.run_id();
  subrun = ioll.subrun_id();
  event  = ioll.event_id();
  //grab array wire image
  auto const& wire_img = ev_img->Image2DArray();
  auto const& wireu_meta = wire_img.at(0).meta();
  auto const& wirev_meta = wire_img.at(1).meta();
  auto const& wirey_meta = wire_img.at(2).meta();
  auto const& instance_img = ev_instance->Image2DArray();
  auto const& segment_img = ev_segment->Image2DArray();
  auto const& ssnetu_img = ssnet_uplane->Image2DArray();
  auto const& ssnetv_img = ssnet_vplane->Image2DArray();
  auto const& ssnety_img = ssnet_yplane->Image2DArray();
  // MCC9 SCE correction
  TFile* newSCEFile = new TFile("/cluster/tufts/wongjiradlab/rshara01/ubdl/SCEoffsets_dataDriven_combined_fwd_Jan18.root");
  TH3F* sceDx = (TH3F*) newSCEFile->Get("hDx");
  TH3F* sceDy = (TH3F*) newSCEFile->Get("hDy");
  TH3F* sceDz = (TH3F*) newSCEFile->Get("hDz");

  //get sce-corrected true Vertex
  //first get 3d version
  true_vtx_location_3D = GetTrueVtxLoc(ev_partroi, sceDx, sceDy, sceDz);
  newSCEFile->Close();
  //get 2d projection of location
  true_vtx_location_2D= Utils.getProjectedPixel(true_vtx_location_3D, wirey_meta, 3);
  vtx_true_fid = Utils.InsideFiducial(true_vtx_location_3D[0],true_vtx_location_3D[1],true_vtx_location_3D[2]);

  //get reco vtx location
  //first get 3d version
  reco_vtx_location_3D_v = GetRecoVtxLocs(ev_vtxtracker);
  //get 2d location Projection
  for (int ii = 0; ii < reco_vtx_location_3D_v.size(); ii++){
    reco_vtx_location_2D_v.push_back(Utils.getProjectedPixel(reco_vtx_location_3D_v[ii], wirey_meta, 3));
  }

  // mark each vtx as good or bad
  vtx_status_v = IsVtxGood(reco_vtx_location_3D_v,true_vtx_location_3D);
  //is each vtx in fiducial volume
  vtx_reco_fid_v = IsVtxInFid(reco_vtx_location_3D_v);
  //truth match each track to particle type
  //format: top vector = vertex level, for each vtx there is a vector of track id
  //see function for id definitions
  true_track_id_v = Probs.TruthMatchTracks(reco_vtx_location_3D_v,ev_recotrack,instance_img,segment_img,wirey_meta);
  //get track length and containment
  TrackLengthAndContainment(reco_vtx_location_3D_v,ev_recotrack);
  // get avg and max dqdx of track
  AvgMaxDQDX(reco_vtx_location_3D_v,ev_recotrack);
  //Get ssnet frac for each plane
  SaveSSNetFrac(reco_vtx_location_3D_v,ev_recotrack,ssnetu_img, wireu_meta, 0);
  SaveSSNetFrac(reco_vtx_location_3D_v,ev_recotrack,ssnetv_img, wirev_meta, 1);
  SaveSSNetFrac(reco_vtx_location_3D_v,ev_recotrack,ssnety_img, wirey_meta, 2);
  //Calculate R Values
  CalculateR_Proton(reco_vtx_location_3D_v,ev_recotrack);
  CalculateR_Gamma(reco_vtx_location_3D_v,ev_recotrack);

  for (int ii = 0; ii< tracklength_v.size(); ii++){
    std::cout<< "\n For vertex ("<<ii<<")"<<std::endl;
    for (int tracknum = 0; tracknum<true_track_id_v[ii].size(); tracknum++){
      std::cout<<" Track ID: "<<true_track_id_v[ii][tracknum];
      // std::cout<<" Track Length: "<<tracklength_v[ii][tracknum];
      // std::cout<<" Is contained? "<<vtx_cont_v[ii][tracknum];
      // std::cout<<" MaxDQDX: "<<max_dqdx_v[ii][tracknum];
      // std::cout<<" AvgDQDX: "<<avg_dqdx_v[ii][tracknum]<<std::endl;
      // std::cout<<"  SSNetShowerU: "<<SSNet_shower_frac_u_v[ii][tracknum];
      // std::cout<<" SSNetTrackU: "<<SSNet_track_frac_u_v[ii][tracknum];
      // std::cout<<" SSNetShowerV: "<<SSNet_shower_frac_v_v[ii][tracknum];
      // std::cout<<" SSNetTrackV: "<<SSNet_track_frac_v_v[ii][tracknum];
      // std::cout<<" SSNetShowerY: "<<SSNet_shower_frac_y_v[ii][tracknum];
      // std::cout<<" SSNetTrackY: "<<SSNet_track_frac_y_v[ii][tracknum]<<std::endl;
      std::cout<<" R(proton): "<<R_proton_v[ii][tracknum];
      // std::cout<<" R(gamma): "<<R_gamma_v[ii][tracknum]<<std::endl;
      std::cout<<std::endl;

    }
  }

  _ana_tree->Fill();
  ClearBranches();
  return true;
}
// -----------------------------------------------------------------------------
void SaveCutVariables::setupAnaTree() {

  std::cout << "Setup analysis tree" << std::endl;
  // ana_file().cd();
  OutFile->cd();
  _ana_tree = new TTree("ncpi0cutvariable", "NCPi0 Cut Variables");

  // CLUSTER RECO VARIABLES
  //create output file and branches
  _ana_tree->Branch("run",&run);
  _ana_tree->Branch("subrun",&subrun);
  _ana_tree->Branch("event",&event);
  _ana_tree->Branch("true_track_id_v",&true_track_id_v);
  _ana_tree->Branch("R_proton_v",&R_proton_v);
  _ana_tree->Branch("R_gamma_v",&R_gamma_v);
  _ana_tree->Branch("SSNet_shower_frac_u_v",&SSNet_shower_frac_u_v);
  _ana_tree->Branch("SSNet_track_frac_u_v",&SSNet_track_frac_u_v);
  _ana_tree->Branch("SSNet_shower_frac_v_v",&SSNet_shower_frac_v_v);
  _ana_tree->Branch("SSNet_track_frac_v_v",&SSNet_track_frac_v_v);
  _ana_tree->Branch("SSNet_shower_frac_y_v",&SSNet_shower_frac_y_v);
  _ana_tree->Branch("SSNet_track_frac_y_v",&SSNet_track_frac_y_v);
  _ana_tree->Branch("tracklength_v",&tracklength_v);
  _ana_tree->Branch("max_dqdx_v",&max_dqdx_v);
  _ana_tree->Branch("avg_dqdx_v",&avg_dqdx_v);
  _ana_tree->Branch("vtx_cont_v",&vtx_cont_v);
  _ana_tree->Branch("vtx_true_fid",&vtx_true_fid);
  _ana_tree->Branch("vtx_reco_fid_v",&vtx_reco_fid_v);
  _ana_tree->Branch("vtx_status_v",&vtx_status_v);
  _ana_tree->Branch("true_vtx_location_3D",&true_vtx_location_3D);
  _ana_tree->Branch("true_vtx_location_2D",&true_vtx_location_2D);
  _ana_tree->Branch("reco_vtx_location_3D_v",&reco_vtx_location_3D_v);
  _ana_tree->Branch("reco_vtx_location_2D_v",&reco_vtx_location_2D_v);

}
//------------------------------------------------------------------------------
void SaveCutVariables::LoadProbHists(){
  fin->cd();
  Utils Utils;
  Muon_length_prob_h = (TH1F*)fin->Get("MuonLength");
	Gamma_length_prob_h = (TH1F*)fin->Get("GammaLength");
	Gamma_ssnetu_h = (TH1F*)fin->Get("GammaSSNetU");
	Gamma_ssnetv_h = (TH1F*)fin->Get("GammaSSNetV");
	Gamma_ssnety_h = (TH1F*)fin->Get("GammaSSNetY");
	Muon_ssnetu_h = (TH1F*)fin->Get("MuonSSNetU");
	Muon_ssnetv_h = (TH1F*)fin->Get("MuonSSNetV");
	Muon_ssnety_h = (TH1F*)fin->Get("MuonSSNetY");
  //scale to probability
	Muon_length_prob_h->Scale(1/Muon_length_prob_h->Integral());
	Gamma_length_prob_h->Scale(1/Gamma_length_prob_h->Integral());
	Muon_ssnetu_h->Scale(1/Muon_ssnetu_h->Integral());
	Muon_ssnetv_h->Scale(1/Muon_ssnetv_h->Integral());
	Muon_ssnety_h->Scale(1/Muon_ssnety_h->Integral());
	Gamma_ssnetu_h->Scale(1/Gamma_ssnetu_h->Integral());
	Gamma_ssnetv_h->Scale(1/Gamma_ssnetv_h->Integral());
	Gamma_ssnety_h->Scale(1/Gamma_ssnety_h->Integral());

  //get and weight prob slices
  for (int ii = 0;ii<20;ii++){
    Muon_resrange_slice_prob_v[ii] = (TH1D*)fin->Get(Form("muon_%d",ii));
    Proton_resrange_slice_prob_v[ii] = (TH1D*)fin->Get(Form("proton_%d",ii));
    if(Muon_resrange_slice_prob_v[ii]->Integral()> 0)Muon_resrange_slice_prob_v[ii]->Scale(1/Muon_resrange_slice_prob_v[ii]->Integral());
    if(Proton_resrange_slice_prob_v[ii]->Integral()> 0) Proton_resrange_slice_prob_v[ii]->Scale(1/Proton_resrange_slice_prob_v[ii]->Integral());
    Muon_resrange_slice_prob_v[ii]->Rebin(4);
    Proton_resrange_slice_prob_v[ii]->Rebin(4);
    Rweight_v[ii] = Utils.CalculateSliceWeight(Muon_resrange_slice_prob_v[ii],Proton_resrange_slice_prob_v[ii]);
  }
}

// -----------------------------------------------------------------------------
std::vector<double> SaveCutVariables::GetTrueVtxLoc(larcv::EventROI* ev_partroi, TH3F* sceDx,
      TH3F* sceDy, TH3F* sceDz){
  // grab true neutrino vtx_new
  //...fill into a vtx to compare to trackervtx...
  std::vector<double> _scex(3,0);
  std::vector<double> _tx(3,0);
  Utils Utils;
  bool isTrueFid;
  for(auto const& roi : ev_partroi->ROIArray()){
    if(std::abs(roi.PdgCode()) == 12 || std::abs(roi.PdgCode()) == 14) {
      _tx.resize(3,0);
      _tx[0] = roi.X();
      _tx[1] = roi.Y();
      _tx[2] = roi.Z();
      auto const offset = Utils.GetPosOffsets(_tx,sceDx,sceDy,sceDz);
      _scex = Utils.MakeSCECorrection(_scex,_tx,offset);
    }
  }
  return _scex;

}//end of function

//------------------------------------------------------------------------------

std::vector<std::vector<double>> SaveCutVariables::GetRecoVtxLocs(larlite::event_vertex* ev_vtxtracker){
  //get 3d locations of vtx points and save them all
  std::vector<std::vector<double>> all_vtx_locations;
  for(int ii = 0; ii<ev_vtxtracker->size(); ii++){
    std::vector<double> vtx_loc;
    vtx_loc.push_back(ev_vtxtracker->at(ii).X());
    vtx_loc.push_back(ev_vtxtracker->at(ii).Y());
    vtx_loc.push_back(ev_vtxtracker->at(ii).Z());
    all_vtx_locations.push_back(vtx_loc);
  }//end of vtx loop
  return all_vtx_locations;
}//end of function

//------------------------------------------------------------------------------
std::vector<bool> SaveCutVariables::IsVtxGood(
    std::vector<std::vector<double>> reco_vtx_v, std::vector<double> true_vtx){
    //function that marks each vtx variable as good or bad depending on
    //distance from true vtx (5cm). output vector of bools, 1 = good, 0 = bad
    std::vector<bool>vtxstatus;
    double truex = true_vtx[0];
    double truey = true_vtx[1];
    double truez = true_vtx[2];
    for (int iii = 0; iii<reco_vtx_v.size(); iii++){
      double recox = reco_vtx_v[iii][0];
      double recoy = reco_vtx_v[iii][1];
      double recoz = reco_vtx_v[iii][2];
      double distance;
      distance = std::sqrt(std::pow(recox-truex,2)+std::pow(recoy-truey,2)+std::pow(recoz-truez,2));
      if (distance < 5) vtxstatus.push_back(true);
      else vtxstatus.push_back(false);
    }//end of vtx loop
    return vtxstatus;
}//end of function

//------------------------------------------------------------------------------

std::vector<bool> SaveCutVariables::IsVtxInFid(std::vector<std::vector<double>> reco_vtx_v){
  //function that marks whether or not each vtx is in fiducial volume
  Utils Utils;
  std::vector<bool> vtxfid_v;
  for (int iii = 0; iii<reco_vtx_v.size(); iii++){
    bool vtxfid = Utils.InsideFiducial(reco_vtx_v[iii][0],reco_vtx_v[iii][1],reco_vtx_v[iii][2]);
    vtxfid_v.push_back(vtxfid);
  }//end of vtx loop
  return vtxfid_v;
}//end of function

//------------------------------------------------------------------------------

void SaveCutVariables::TrackLengthAndContainment(
  std::vector<std::vector<double>> reco_vtx_v, larlite::event_track*ev_recotrack){
  //function to calculate track length and if contained.

  Utils Utils;
  SaveProbabilities Probs;
  for (int vtx = 0; vtx<reco_vtx_v.size(); vtx++){
    double xcoor = reco_vtx_v[vtx][0];
    double ycoor = reco_vtx_v[vtx][1];
    double zcoor = reco_vtx_v[vtx][2];
    std::vector<float> tracklength;
    std::vector<bool> contained;
    for(int ii = 0; ii<ev_recotrack->size(); ii++){
      TVector3 trackvtx = ev_recotrack->at(ii).Vertex();
      float xvtx= trackvtx.X();
      float yvtx= trackvtx.Y();
      float zvtx= trackvtx.Z();
      TVector3 trackvtxend = ev_recotrack->at(ii).End();
      float xvtxend = trackvtxend.X();
      float yvtxend = trackvtxend.Y();
      float zvtxend = trackvtxend.Z();
      //loop through vertices to find closest match
      int matchindex = Probs.MatchVtxToTrack(reco_vtx_v,xvtx,yvtx,zvtx);
      //check if matches vtx id
      if(matchindex == vtx){
        //calculate and save tracklength and containment
        float length =std::sqrt(std::pow(xvtxend-xvtx,2)+std::pow(yvtxend-yvtx,2)+std::pow(zvtxend-zvtx,2));
        tracklength.push_back(length);
        bool cont = Utils.InsideFiducial(xvtxend,yvtxend,zvtxend);
        contained.push_back(cont);
      }//end of if vtx match statement
    }//end of track loop
    tracklength_v.push_back(tracklength);

    vtx_cont_v.push_back(contained);
  }//end of vtx loop
  return;

}//end of function

//------------------------------------------------------------------------------
void SaveCutVariables::AvgMaxDQDX(
  std::vector<std::vector<double>> reco_vtx_v, larlite::event_track*ev_recotrack){
  //function to calculate avg and max dqdx and save to tree variable
  Utils Utils;
  SaveProbabilities Probs;
  for (int vtx = 0; vtx<reco_vtx_v.size(); vtx++){
    std::vector<float> avg_dqdx;
    std::vector<float> max_dqdx;
    double xcoor = reco_vtx_v[vtx][0];
    double ycoor = reco_vtx_v[vtx][1];
    double zcoor = reco_vtx_v[vtx][2];
    for(int ii = 0; ii<ev_recotrack->size(); ii++){
      TVector3 trackvtx = ev_recotrack->at(ii).Vertex();
      float xvtx= trackvtx.X();
      float yvtx= trackvtx.Y();
      float zvtx= trackvtx.Z();
      TVector3 trackvtxend = ev_recotrack->at(ii).End();
      float xvtxend = trackvtxend.X();
      float yvtxend = trackvtxend.Y();
      float zvtxend = trackvtxend.Z();
      //loop through vertices to find closest match
      int matchindex = Probs.MatchVtxToTrack(reco_vtx_v,xvtx,yvtx,zvtx);
      //check if matches vtx id
      if(matchindex == vtx){
        float totalcharge = 0;
        float max = 0;
        for (int pt = 0; pt<ev_recotrack->at(ii).NumberdQdx(larlite::geo::View_t::kW);pt++){
          totalcharge += ev_recotrack->at(ii).DQdxAtPoint(pt,larlite::geo::View_t::kW);
          if (ev_recotrack->at(ii).DQdxAtPoint(pt,larlite::geo::View_t::kW)>max){
            max = ev_recotrack->at(ii).DQdxAtPoint(pt,larlite::geo::View_t::kW);
          }
        }
        float avg = totalcharge/(float)ev_recotrack->at(ii).NumberdQdx(larlite::geo::View_t::kW);
        avg_dqdx.push_back(avg);
        max_dqdx.push_back(max);
      }//end of if match vtx
    }//end of loop through tracks
    avg_dqdx_v.push_back(avg_dqdx);
    max_dqdx_v.push_back(max_dqdx);
  }//end of loop through verticies
  return;
}//end of function

//------------------------------------------------------------------------------

void SaveCutVariables::SaveSSNetFrac(std::vector<std::vector<double>> reco_vtx_v,
  larlite::event_track*ev_recotrack,std::vector<larcv::Image2D> ssnet_img,
  larcv::ImageMeta wire_meta, int plane){
  //function to go through tracks and save ssnet track and shower fraction on a plane
  Utils Utils;
  SaveProbabilities Probs;
  for (int vtx = 0; vtx<reco_vtx_v.size(); vtx++){
    std::vector<float> showerfrac_v;
    std::vector<float> trackfrac_v;
    double xcoor = reco_vtx_v[vtx][0];
    double ycoor = reco_vtx_v[vtx][1];
    double zcoor = reco_vtx_v[vtx][2];
    for(int ii = 0; ii<ev_recotrack->size(); ii++){
      TVector3 trackvtx = ev_recotrack->at(ii).Vertex();
      float xvtx= trackvtx.X();
      float yvtx= trackvtx.Y();
      float zvtx= trackvtx.Z();
      TVector3 trackvtxend = ev_recotrack->at(ii).End();
      float xvtxend = trackvtxend.X();
      float yvtxend = trackvtxend.Y();
      float zvtxend = trackvtxend.Z();
      //loop through vertices to find closest match
      int matchindex = Probs.MatchVtxToTrack(reco_vtx_v,xvtx,yvtx,zvtx);
      //check if matches vtx id
      if(matchindex == vtx){
        //ssnet loop
        float showercount = 0;
        float trackcount = 0;
        float numwssnet = 0;
        bool nossnet = true;
        int totalpts = ev_recotrack->at(ii).NumberTrajectoryPoints();

        auto const& ssnet_meta = ssnet_img.at(0).meta();

        for(int pt = 0; pt<totalpts ;pt++){
          TVector3 pt_3vec=ev_recotrack->at(ii).LocationAtPoint(pt);
          std::vector<double> pt_vec;
          pt_vec.push_back(pt_3vec.X());
          pt_vec.push_back(pt_3vec.Y());
          pt_vec.push_back(pt_3vec.Z());
          std::vector<int> ssnet_pixel = Utils.getProjectedPixel(pt_vec,ssnet_meta,3);
          float trackval = ssnet_img.at(1).pixel(ssnet_pixel[0],ssnet_pixel[plane+1]);
          float showerval = ssnet_img.at(0).pixel(ssnet_pixel[0],ssnet_pixel[plane+1]);

          if(showerval>0 || trackval>0){
            if (showerval>.5)showercount++;
            else trackcount++;
            numwssnet++;
          }
        }//end of loop over points
        float showerfrac = 0;
        float trackfrac = 0;
        if(numwssnet>0){
          showerfrac = (float)showercount/(float)numwssnet;
          trackfrac = (float)trackcount/(float)numwssnet;
        }
        showerfrac_v.push_back(showerfrac);
        trackfrac_v.push_back(trackfrac);
      }//end of if match vtx
    }//end of loop through tracks
    if (plane == 0){
      SSNet_shower_frac_u_v.push_back(showerfrac_v);
      SSNet_track_frac_u_v.push_back(trackfrac_v);
    }
    else if (plane == 1){
      SSNet_shower_frac_v_v.push_back(showerfrac_v);
      SSNet_track_frac_v_v.push_back(trackfrac_v);
    }
    else if (plane == 2){
      SSNet_shower_frac_y_v.push_back(showerfrac_v);
      SSNet_track_frac_y_v.push_back(trackfrac_v);
    }
    else std::cout<<"BAD PLANE ID"<<std::endl;
  }//end of loop through verticies
  return;

}//end of function
//-----------------------------------------------------------------------------
void SaveCutVariables::CalculateR_Proton(std::vector<std::vector<double>> reco_vtx_v,
  larlite::event_track*ev_recotrack){
  Utils Utils;
  SaveProbabilities Probs;
  for (int vtx = 0; vtx<reco_vtx_v.size(); vtx++){
    double xcoor = reco_vtx_v[vtx][0];
    double ycoor = reco_vtx_v[vtx][1];
    double zcoor = reco_vtx_v[vtx][2];
    std::vector<float> R_tracks_proton;
    for(int ii = 0; ii<ev_recotrack->size(); ii++){
      TVector3 trackvtx = ev_recotrack->at(ii).Vertex();
      float xvtx= trackvtx.X();
      float yvtx= trackvtx.Y();
      float zvtx= trackvtx.Z();
      TVector3 trackvtxend = ev_recotrack->at(ii).End();
      float xvtxend = trackvtxend.X();
      float yvtxend = trackvtxend.Y();
      float zvtxend = trackvtxend.Z();
      //loop through vertices to find closest match
      int matchindex = Probs.MatchVtxToTrack(reco_vtx_v,xvtx,yvtx,zvtx);
      //check if matches vtx id
      std::vector<float> RVal_proton_v;
      std::vector<float> weights_v;
      if(matchindex == vtx){
        for (int iii = 0; iii<ev_recotrack->at(ii).NumberTrajectoryPoints();iii++){
          float dqdx_val = -1;
          dqdx_val = ev_recotrack->at(ii).DQdxAtPoint(iii,larlite::geo::View_t::kW);
          if (dqdx_val == 0){
            float uval = ev_recotrack->at(ii).DQdxAtPoint(iii,larlite::geo::View_t::kU);
            float vval = ev_recotrack->at(ii).DQdxAtPoint(iii,larlite::geo::View_t::kV);
            dqdx_val = (uval+vval)/2.0;
          }
          TVector3 trackpt = ev_recotrack->at(ii).LocationAtPoint(iii);
          float trackpt_x = trackpt.X();
          float trackpt_y = trackpt.Y();
          float trackpt_z = trackpt.Z();
          float resrange = std::sqrt(std::pow(trackpt_x - xvtxend,2)+std::pow(trackpt_y - yvtxend,2)+std::pow(trackpt_z - zvtxend,2));
          //get prob slice number
          int slice;
          float muonprob;
          float protonprob;
          int rebin =4;
          if (resrange<10 &&dqdx_val > 4){
            slice = Utils.GetSliceNumber(resrange);
            muonprob = Muon_resrange_slice_prob_v[slice]->GetBinContent(dqdx_val/rebin+1);
            protonprob = Proton_resrange_slice_prob_v[slice]->GetBinContent(dqdx_val/rebin+1);
            if (muonprob == 0) muonprob = 1e-5;
            if (protonprob == 0) protonprob = 1e-5;
            float R_value;
            R_value = log(protonprob)-log(muonprob);
            RVal_proton_v.push_back(R_value*Rweight_v[slice]);
            weights_v.push_back(Rweight_v[slice]);
          }
        }
        //get total R Proton
        float totalR_proton = 0;
        float total_weights = 0;
        for (int rr =0 ;rr<RVal_proton_v.size();rr++){
          totalR_proton = totalR_proton+RVal_proton_v[rr];
          total_weights = total_weights + weights_v[rr];
        }

        R_tracks_proton.push_back(totalR_proton/total_weights);
      }//end of if match
    }//end of loop through tracks
    R_proton_v.push_back(R_tracks_proton);
  }//end of loop through verticies
  return;
}//end of function
//-----------------------------------------------------------------------------
void SaveCutVariables::CalculateR_Gamma(std::vector<std::vector<double>> reco_vtx_v,
  larlite::event_track*ev_recotrack){
  Utils Utils;
  SaveProbabilities Probs;
  for (int vtx = 0; vtx<reco_vtx_v.size(); vtx++){
    double xcoor = reco_vtx_v[vtx][0];
    double ycoor = reco_vtx_v[vtx][1];
    double zcoor = reco_vtx_v[vtx][2];
    int trackid = 0;
    std::vector<float> R_tracks_gamma;
    for(int ii = 0; ii<ev_recotrack->size(); ii++){
      TVector3 trackvtx = ev_recotrack->at(ii).Vertex();
      float xvtx= trackvtx.X();
      float yvtx= trackvtx.Y();
      float zvtx= trackvtx.Z();
      TVector3 trackvtxend = ev_recotrack->at(ii).End();
      float xvtxend = trackvtxend.X();
      float yvtxend = trackvtxend.Y();
      float zvtxend = trackvtxend.Z();
      //loop through vertices to find closest match
      int matchindex = Probs.MatchVtxToTrack(reco_vtx_v,xvtx,yvtx,zvtx);
      //check if matches vtx id
      if(matchindex == vtx){
        float muon_length_prob;
        float muon_ssnetu_prob;
        float muon_ssnetv_prob;
        float muon_ssnety_prob;
        float gamma_length_prob;
        float gamma_ssnetu_prob;
        float gamma_ssnetv_prob;
        float gamma_ssnety_prob;
        float tracklength = std::sqrt(std::pow(xvtxend - xvtx,2)+std::pow(yvtxend - yvtx,2)+std::pow(zvtxend - zvtx,2));


        gamma_length_prob=Gamma_length_prob_h->GetBinContent(tracklength/2+1);
        if (gamma_length_prob < 10e-5)gamma_length_prob = 10e-5;
        muon_length_prob=Muon_length_prob_h->GetBinContent(tracklength/2+1);
        if (muon_length_prob < 10e-5)muon_length_prob = 10e-5;
        muon_ssnetu_prob=Muon_ssnetu_h->GetBinContent(SSNet_shower_frac_u_v[vtx][trackid]*20+1);
        if (muon_ssnetu_prob < 10e-5)muon_ssnetu_prob = 10e-5;
        muon_ssnetv_prob=Muon_ssnetv_h->GetBinContent(SSNet_shower_frac_v_v[vtx][trackid]*20+1);
        if (muon_ssnetv_prob < 10e-5)muon_ssnetv_prob = 10e-5;
        muon_ssnety_prob=Muon_ssnety_h->GetBinContent(SSNet_shower_frac_y_v[vtx][trackid]*20+1);
        if (muon_ssnety_prob < 10e-5)muon_ssnety_prob = 10e-5;
        gamma_ssnetu_prob=Gamma_ssnetu_h->GetBinContent(SSNet_shower_frac_u_v[vtx][trackid]*20+1);
        if (gamma_ssnetu_prob < 10e-5)gamma_ssnetu_prob = 10e-5;
        gamma_ssnetv_prob=Gamma_ssnetv_h->GetBinContent(SSNet_shower_frac_v_v[vtx][trackid]*20+1);
        if (gamma_ssnetv_prob < 10e-5)gamma_ssnetv_prob = 10e-5;
        gamma_ssnety_prob=Gamma_ssnety_h->GetBinContent(SSNet_shower_frac_y_v[vtx][trackid]*20+1);
        if (gamma_ssnety_prob < 10e-5)gamma_ssnety_prob = 10e-5;

        trackid++;

        float ll_gamma_length = log(gamma_length_prob)-log(muon_length_prob);
        float ll_gamma_ssnetu = log(gamma_ssnetu_prob)-log(muon_ssnetu_prob);
        float ll_gamma_ssnetv = log(gamma_ssnetv_prob)-log(muon_ssnetv_prob);
        float ll_gamma_ssnety = log(gamma_ssnety_prob)-log(muon_ssnety_prob);
        float totalR_gamma = (ll_gamma_length+ll_gamma_ssnetu+ll_gamma_ssnetv+ll_gamma_ssnety);

        R_tracks_gamma.push_back(totalR_gamma);



      }//end of if match
    }//end of loop through tracks
    R_gamma_v.push_back(R_tracks_gamma);
  }//end of loop through verticies
  return;
}//end of function
//------------------------------------------------------------------------------
void SaveCutVariables::ClearBranches(){
  //clear all branches
  true_track_id_v.clear();
  R_proton_v.clear();
  R_gamma_v.clear();
  SSNet_shower_frac_u_v.clear();
  SSNet_track_frac_u_v.clear();
  SSNet_shower_frac_v_v.clear();
  SSNet_track_frac_v_v.clear();
  SSNet_shower_frac_y_v.clear();
  SSNet_track_frac_y_v.clear();
  tracklength_v.clear();
  max_dqdx_v.clear();
  avg_dqdx_v.clear();
  vtx_cont_v.clear();
  vtx_reco_fid_v.clear();
  vtx_status_v.clear();
  true_vtx_location_3D.clear();
  true_vtx_location_2D.clear();
  reco_vtx_location_3D_v.clear();
  reco_vtx_location_2D_v.clear();
}//end of function

//------------------------------------------------------------------------------

void SaveCutVariables::finalize() {
  // if ( has_ana_file() && _ana_tree ) {
  //   _ana_tree->Write();
  //   OutFile->Write();
  // }
  OutFile->cd();
  _ana_tree->Write();
  OutFile->Write();

}
}//end of ncpi0 namespace
}//end of ublarcvapp namespace
