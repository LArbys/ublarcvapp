#include "SaveProbabilities.h"

namespace ublarcvapp {
namespace ncpi0{

void SaveProbabilities::configure( const larcv::PSet& pset ) {
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

  void SaveProbabilities::initialize() {
    OutFile = new TFile("NCPi0Probablities.root","RECREATE");
    setupAnaTree();
  }

bool SaveProbabilities::process(larcv::IOManager& io,larlite::storage_manager& ioll,larcv::IOManager& ioforward){
  //set object
  Utils Utils;
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

  //get reco vtx location
  //first get 3d version
  std::vector<std::vector<double>> reco_vtx_location_3D_v = GetRecoVtxLocs(ev_vtxtracker);
  //truth match each track to particle type
  //format: top vector = vertex level, for each vtx there is a vector of track id
  //see function for id definitions
  std::vector<std::vector<int>> true_track_id_v = TruthMatchTracks(reco_vtx_location_3D_v,ev_recotrack,instance_img,segment_img,wirey_meta);
  // calculate dqdx slices
  DQDXProbs(true_track_id_v, reco_vtx_location_3D_v, ev_recotrack);
  // calculate track length
  TrackLengthProbs(true_track_id_v, reco_vtx_location_3D_v, ev_recotrack);
  // calculate ssnet shower frac
  SSNetShowerProbs(reco_vtx_location_3D_v,ev_recotrack,ssnetu_img,wire_img, wireu_meta, 0,true_track_id_v);
  SSNetShowerProbs(reco_vtx_location_3D_v,ev_recotrack,ssnetv_img,wire_img, wirev_meta, 1,true_track_id_v);
  SSNetShowerProbs(reco_vtx_location_3D_v,ev_recotrack,ssnety_img,wire_img, wirey_meta, 2,true_track_id_v);
  //Fill slice histograms
  GetResRangeSlices();


	// // loop through tracker vertices to get characteristics
	// for(int ii = 0; ii<ev_recotrack->size(); ii++){
	// 	int totalpts = ev_recotrack->at(ii).NumberTrajectoryPoints();
	// 	for(int pt = 0; pt<totalpts ;pt++){
	// 		TVector3 pt_3vec=ev_recotrack->at(ii).LocationAtPoint(pt);
	// 		std::vector<double> pt_vec;
	// 		pt_vec.push_back(pt_3vec.X());
	// 		pt_vec.push_back(pt_3vec.Y());
	// 		pt_vec.push_back(pt_3vec.Z());
	// 		std::vector<int> track_pixel =Utils.getProjectedPixel(pt_vec,wirey_meta,3);
	// 		tracker_image->SetBinContent(track_pixel[3],track_pixel[0],1);
	// 	}
	// }
  //
	// TCanvas can3("can", "histograms ", 3456, 1008);
	// tracker_image->SetMarkerStyle(kCircle);
  // tracker_image->Draw("");
  // adc_image->Draw("SAME COLZ");
  // can3.SaveAs("Tracks.png");

  return true;
}
//------------------------------------------------------------------------------

std::vector<std::vector<double>> SaveProbabilities::GetRecoVtxLocs(larlite::event_vertex* ev_vtxtracker){
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
int SaveProbabilities::MatchVtxToTrack(std::vector<std::vector<double>> reco_vtx_v,
  float xvtx,float yvtx,float zvtx){
  //returns the index of the closest vtx to the track
  int matchindex = -1;
  float mindistance = 1000;
  for(int vtxii = 0; vtxii<reco_vtx_v.size();vtxii++ ){
    double xcoortmp = reco_vtx_v[vtxii][0];
    double ycoortmp = reco_vtx_v[vtxii][1];
    double zcoortmp = reco_vtx_v[vtxii][2];
    float distance = std::sqrt(std::pow(xcoortmp-xvtx,2)+std::pow(ycoortmp-yvtx,2)+std::pow(zcoortmp-zvtx,2));
    if (distance<mindistance){
      mindistance = distance;
      matchindex = vtxii;
    }
  }
  return matchindex;
}//end of function
// -----------------------------------------------------------------------------

std::vector<std::vector<int>> SaveProbabilities::TruthMatchTracks(
  std::vector<std::vector<double>> reco_vtx_v, larlite::event_track* ev_recotrack,
  std::vector<larcv::Image2D> instance_img, std::vector<larcv::Image2D> segment_img,
  larcv::ImageMeta wirey_meta){
  //function that loops through tracks and truth matches them to recovtx and particle type
  // choose id of largest fraction of pixels in track
  //ids (based on segment img):
  //  0: unknown
  //  1: cosmic muon --- not applicable
  //  2: BNB ---not applicable
  //  3: Electron
  //  4: Gamma
  //  5: Pi0
  //  6: Muon
  //  7: Kaon
  //  8: Charged Pion
  //  9: Proton
  Utils Utils;
  std::vector<std::vector<int>> trackid_v;
  for (int vtx = 0; vtx<reco_vtx_v.size(); vtx++){
    std::vector<int> trackid;
    double xcoor = reco_vtx_v[vtx][0];
    double ycoor = reco_vtx_v[vtx][1];
    double zcoor = reco_vtx_v[vtx][2];
    int NumMuonCosmic   = 0;
    int NumMuonNeutrino = 0;
    int NumElectron     = 0;
    int NumGamma        = 0;
    int NumPi0          = 0;
    int NumPiCharged    = 0;
    int NumKaon         = 0;
    int NumProton       = 0;
    // std::cout<<"Number of Reco tracks: "<<ev_recotrack->size()<<std::endl;
    for(int ii = 0; ii<ev_recotrack->size(); ii++){
      TVector3 trackvtx = ev_recotrack->at(ii).Vertex();
      float xvtx= trackvtx.X();
      float yvtx= trackvtx.Y();
      float zvtx= trackvtx.Z();
      //loop through vertices to find closest match
      int matchindex = MatchVtxToTrack(reco_vtx_v,xvtx,yvtx,zvtx);
      //check if matches vtx id
      if(matchindex == vtx){
        // to do truth match to particle type
        // -- loop through pixels, get instance image of each pixel
        // -- type is whatever has highest number of instance image
				int NumMuonCosmicPix   = 0;
        int NumMuonNeutrinoPix = 0;
        int NumElectronPix     = 0;
        int NumGammaPix        = 0;
        int NumPi0Pix          = 0;
        int NumPiChargedPix    = 0;
        int NumKaonPix         = 0;
        int NumProtonPix       = 0;
        int numtrackpts = 0;
        std::vector<double> prev_pix (2,0.0);
				for (int iii = 0; iii<ev_recotrack->at(ii).NumberTrajectoryPoints(); iii++){
					TVector3 trackpt = ev_recotrack->at(ii).LocationAtPoint(iii);
					std::vector<double> trackpt_v;
					trackpt_v.push_back(trackpt.X());
					trackpt_v.push_back(trackpt.Y());
					trackpt_v.push_back(trackpt.Z());
					std::vector<int> trackpt_pix = Utils.getProjectedPixel(trackpt_v,wirey_meta,3);
          // std::cout<<"trackpt: "<<trackpt_pix[0]<< " , "<<trackpt_pix[3]<<std::endl;
          if(prev_pix[0]!=trackpt_pix[0] || prev_pix[1]!=trackpt_pix[3]){
            // std::cout<<" @ point: "<<trackpt_pix[0]<<" , "<<trackpt_pix[3]<<std::endl;
            float instance_val = instance_img[2].pixel(trackpt_pix[0],trackpt_pix[3]);
						if (instance_val > 0){
							//get segement value
							float segment_val = segment_img[2].pixel(trackpt_pix[0],trackpt_pix[3]);
							// segment types in larcv/core/dataformat/DataFormatTypes.h
							if (segment_val==3) NumElectronPix++;
              else if (segment_val==4) NumGammaPix++;
              else if (segment_val==5) NumPi0Pix++;
              else if (segment_val==6) NumMuonNeutrinoPix++;
              else if (segment_val==7) NumKaonPix++;
              else if (segment_val==8) NumPiChargedPix++;
              else if (segment_val==9) NumProtonPix++;
						}
						else NumMuonCosmicPix++;
            prev_pix[0]=trackpt_pix[0];
            prev_pix[1]=trackpt_pix[3];
            numtrackpts++;
          }

				}//end of loop through track points

        float FracProton       = (float)NumProtonPix/((float)numtrackpts);
				float FracMuonCosmic   = (float)NumMuonCosmicPix/((float)numtrackpts);
        float FracMuonNeutrino = (float)NumMuonNeutrinoPix/((float)numtrackpts);
        float FracGamma        = (float)NumGammaPix/((float)numtrackpts);
        float FracPi0          = (float)NumPi0Pix/((float)numtrackpts);
        float FracPiCharged    = (float)NumPiChargedPix/((float)numtrackpts);
        float FracKaon         = (float)NumKaonPix/((float)numtrackpts);
        float FracElectron     = (float)NumElectronPix/((float)numtrackpts);

        //order of adding: gamma,electron,pi0,proton,picharged,kaon,muon,cosmicmuon;
        int id = 0;
        if (FracGamma>=.5)             id = 4;
        else if (FracElectron>=.5)     id = 3;
        else if (FracPi0>=.5)          id = 5;
        else if (FracProton>=.5)       id = 9;
        else if (FracPiCharged>=.5)    id = 8;
        else if (FracKaon>=.5)         id = 7;
        else if (FracMuonNeutrino>=.5) id = 6;
        else if (FracMuonCosmic>=.5)   id = 1;

        if (id ==0){
          std::cout<<"Unknown Particle"<<std::endl;
          std::cout<<"FracProton: "<<FracProton<<std::endl;
          std::cout<<"FracGamma: "<<FracGamma<<std::endl;
          std::cout<<"FracCosmic: "<<FracMuonCosmic<<std::endl;
          std::cout<<"FracMuon: "<<FracMuonNeutrino<<std::endl;
          std::cout<<"FracElectron: "<<FracElectron<<std::endl;
          std::cout<<"FracKaon: "<<FracKaon<<std::endl;
          std::cout<<"FracPi0: "<<FracPi0<<std::endl;
          std::cout<<"FracPiCharged: "<<FracPiCharged<<std::endl;
        }


        trackid.push_back(id);
      }//end of if vtx match statement
    }//end of track loop
    trackid_v.push_back(trackid);
  }//end of vtx loop
  return trackid_v;
}//end of function

//------------------------------------------------------------------------------
void SaveProbabilities::SSNetShowerProbs(std::vector<std::vector<double>> reco_vtx_v,
  larlite::event_track*ev_recotrack,std::vector<larcv::Image2D> ssnet_img,std::vector<larcv::Image2D> wire_img,
  larcv::ImageMeta wire_meta, int plane,std::vector<std::vector<int>> trueid){
  //function to go through tracks and save ssnet track and shower fraction on a plane
  Utils Utils;
  SaveProbabilities Probs;
  for (int vtx = 0; vtx<reco_vtx_v.size(); vtx++){
    double xcoor = reco_vtx_v[vtx][0];
    double ycoor = reco_vtx_v[vtx][1];
    double zcoor = reco_vtx_v[vtx][2];
    int trackid = 0;
    for(int ii = 0; ii<ev_recotrack->size(); ii++){
      TVector3 trackvtx = ev_recotrack->at(ii).Vertex();
      float xvtx= trackvtx.X();
      float yvtx= trackvtx.Y();
      float zvtx= trackvtx.Z();
      //loop through vertices to find closest match
      int matchindex = Probs.MatchVtxToTrack(reco_vtx_v,xvtx,yvtx,zvtx);
      //check if matches vtx id
      if(matchindex == vtx){
        //ssnet loop
        float showercount = 0;
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
            numwssnet++;
          }
          pt_vec.clear();
        }//end of loop over points
        float showerfrac = 0;
        if(numwssnet>0){
          showerfrac = (float)showercount/(float)numwssnet;
        }
        if (plane == 0 && trueid[vtx][trackid]==1) Muon_ssnetu_h->Fill(showerfrac,1);
        else if (plane == 1 && trueid[vtx][trackid]==1) Muon_ssnetv_h->Fill(showerfrac,1);
        else if (plane == 2 && trueid[vtx][trackid]==1) Muon_ssnety_h->Fill(showerfrac,1);
        else if (plane == 0 && trueid[vtx][trackid]==4) Gamma_ssnetu_h->Fill(showerfrac,1);
        else if (plane == 1 && trueid[vtx][trackid]==4) Gamma_ssnetv_h->Fill(showerfrac,1);
        else if (plane == 2 && trueid[vtx][trackid]==4) Gamma_ssnety_h->Fill(showerfrac,1);
        trackid++;
      }//end of if match vtx
    }//end of loop through tracks
  }//end of loop through verticies
  return;

}//end of function
//------------------------------------------------------------------------------

void SaveProbabilities::TrackLengthProbs(std::vector<std::vector<int>> true_track_id_v,
  std::vector<std::vector<double>> reco_vtx_v, larlite::event_track*ev_recotrack){
  //function to calculate track length and if contained.

  Utils Utils;
  for (int vtx = 0; vtx<reco_vtx_v.size(); vtx++){
    double xcoor = reco_vtx_v[vtx][0];
    double ycoor = reco_vtx_v[vtx][1];
    double zcoor = reco_vtx_v[vtx][2];
    int trackid = 0;
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
      int matchindex = MatchVtxToTrack(reco_vtx_v,xvtx,yvtx,zvtx);
      //check if matches vtx id
      if(matchindex == vtx){
        //calculate and save tracklength and containment
        float length =std::sqrt(std::pow(xvtxend-xvtx,2)+std::pow(yvtxend-yvtx,2)+std::pow(zvtxend-zvtx,2));
        if(length>199.9)length =199.9;
        if (true_track_id_v[vtx][trackid]==1)Muon_length_prob_h->Fill(length,1);
        else if(true_track_id_v[vtx][trackid]==4)Gamma_length_prob_h->Fill(length,1);
        trackid++;
      }//end of if vtx match statement
    }//end of track loop
  }//end of vtx loop
  return;

}//end of function

//------------------------------------------------------------------------------
void SaveProbabilities::DQDXProbs(std::vector<std::vector<int>> true_track_id_v,
    std::vector<std::vector<double>> reco_vtx_v, larlite::event_track*ev_recotrack){
  //function to calculate avg and max dqdx and save to tree variable
  Utils Utils;
  for (int vtx = 0; vtx<reco_vtx_v.size(); vtx++){
    double xcoor = reco_vtx_v[vtx][0];
    double ycoor = reco_vtx_v[vtx][1];
    double zcoor = reco_vtx_v[vtx][2];
    int trackid = 0;
    for(int ii = 0; ii<ev_recotrack->size(); ii++){
      TVector3 trackvtx = ev_recotrack->at(ii).Vertex();
      float xvtx= trackvtx.X();
      float yvtx= trackvtx.Y();
      float zvtx= trackvtx.Z();
      //loop through vertices to find closest match
      int matchindex = MatchVtxToTrack(reco_vtx_v,xvtx,yvtx,zvtx);
      //check if matches vtx id
      if(matchindex == vtx){
        for (int iii = 0; iii<ev_recotrack->at(ii).NumberTrajectoryPoints();iii++){
          float dqdx_val = ev_recotrack->at(ii).DQdxAtPoint(iii,larlite::geo::View_t::kW);
          if (dqdx_val == 0){
            float uval = ev_recotrack->at(ii).DQdxAtPoint(iii,larlite::geo::View_t::kU);
            float vval = ev_recotrack->at(ii).DQdxAtPoint(iii,larlite::geo::View_t::kV);
            dqdx_val = (uval+vval)/2.0;
          }
          TVector3 trackpt = ev_recotrack->at(ii).LocationAtPoint(iii);
          float trackpt_x = trackpt.X();
          float trackpt_y = trackpt.Y();
          float trackpt_z = trackpt.Z();
          TVector3 trackendpt = ev_recotrack->at(ii).End();
					float trackendpt_x = trackendpt.X();
					float trackendpt_y = trackendpt.Y();
					float trackendpt_z = trackendpt.Z();
          float resrange = std::sqrt(std::pow(trackpt_x - trackendpt_x,2)+std::pow(trackpt_y - trackendpt_y,2)+std::pow(trackpt_z - trackendpt_z,2));
          //fill dq/dx plots
          if (true_track_id_v[vtx][trackid] == 9 && dqdx_val>0){
            Proton_dqdx_resrange_h->Fill(resrange,dqdx_val,1);
          }
          else if (true_track_id_v[vtx][trackid] == 1 && dqdx_val>0){
            Muon_dqdx_resrange_h->Fill(resrange,dqdx_val,1);
          }
          trackid++;
        }//end of dqdx loop
      }//end of if match vtx
    }//end of loop through tracks
  }//end of loop through verticies
  return;
}//end of function
//------------------------------------------------------------------------------
void SaveProbabilities::GetResRangeSlices(){
  for (int ii = 0;ii<numslices;ii++){
		Muon_resrange_slice_v[ii] = Muon_dqdx_resrange_h->ProjectionY(Form("muon_%d",ii),ii,ii+1);
		Proton_resrange_slice_v[ii] = Proton_dqdx_resrange_h->ProjectionY(Form("proton_%d",ii),ii,ii+1);
  }//end of loop over slices
}//end of function
//------------------------------------------------------------------------------
void SaveProbabilities::setupAnaTree() {

  std::cout << "Setup analysis tree" << std::endl;
  // ana_file().cd();
  OutFile->cd();
  // setup histograms
  //for R(Gamma)
  Muon_dqdx_resrange_h = new TH2F("MuonDQDX","MuonDQDX",200,0,100,200,0,200);
	Proton_dqdx_resrange_h = new TH2F("ProtonDQDX","ProtonDQDX",200,0,100,200,0,200);
  Muon_length_prob_h = new TH1F("MuonLength","MuonLength",100,0,200);
	Gamma_length_prob_h = new TH1F("GammaLength","GammaLength",100,0,200);
	Gamma_ssnetu_h = new TH1F("GammaSSNetU","GammaSSNetU",20,0,1);
	Gamma_ssnetv_h = new TH1F("GammaSSNetV","GammaSSNetV",20,0,1);
	Gamma_ssnety_h = new TH1F("GammaSSNetY","GammaSSNetY",20,0,1);
	Muon_ssnetu_h = new TH1F("MuonSSNetU","MuonSSNetU",20,0,1);
	Muon_ssnetv_h = new TH1F("MuonSSNetV","MuonSSNetV",20,0,1);
	Muon_ssnety_h = new TH1F("MuonSSNetY","MuonSSNetY",20,0,1);


}


//------------------------------------------------------------------------------

void SaveProbabilities::finalize() {

  OutFile->cd();
  // _ana_tree->Write();
  OutFile->Write();

}
}//end of ncpi0 namespace
}//end of ublarcvapp namespace
