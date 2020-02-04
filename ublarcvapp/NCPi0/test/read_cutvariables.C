#include <iostream>
#include <map>
#include <utility>
// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TPaveLabel.h"

#include <vector>

//Macro to read in info from root files, save histograms, and determine Values
//made to be used with ncpi0

TFile* fin;
TTree* tree;
int run =-1;
int subrun =-1;
int event =-1;

// // Initalize values from file
std::vector<std::vector<int> >* true_track_id_v = NULL; //vector storing track id for reco tracks for each vtx
std::vector<std::vector<float> >* R_proton_v = NULL; //vector of R(Proton) for recotracks for each vtx
std::vector<std::vector<float> >* R_gamma_v = NULL; //vector of R(Proton) for recotracks for each vtx
std::vector<std::vector<float>>* SSNet_shower_frac_u_v = NULL;//vector of SSNet(uplane) shower fraction for recotracks for each vtx
std::vector<std::vector<float>>* SSNet_track_frac_u_v = NULL;//vector of SSNet(uplane) shower fraction for recotracks for each vtx
std::vector<std::vector<float>>* SSNet_shower_frac_v_v = NULL;//vector of SSNet(vplane) shower fraction for recotracks for each vtx
std::vector<std::vector<float>>* SSNet_track_frac_v_v = NULL;//vector of SSNet(vplane) shower fraction for recotracks for each vtx
std::vector<std::vector<float>>* SSNet_shower_frac_y_v = NULL;//vector of SSNet(yplane) shower fraction for recotracks for each vtx
std::vector<std::vector<float>>* SSNet_track_frac_y_v = NULL;//vector of SSNet(yplane) shower fraction for recotracks for each vtx
std::vector<std::vector<float>>* tracklength_v = NULL;//vector of track length of recotracks for eachvtx;
std::vector<std::vector<float>>* avg_dqdx_v = NULL;// vector of avg dqdx of recotracks for each vtx;
std::vector<std::vector<float>>* max_dqdx_v = NULL;// vector of max dqdx of recotracks for each vtx;
std::vector<std::vector<bool>>* vtx_cont_v = NULL;//vector of boolians of is track contained?
bool vtx_true_fid;//vector of boolians on whether vtx is in fid
std::vector<bool>* vtx_reco_fid_v = NULL;//vector of boolians on whether vtx is in fid
std::vector<bool>* vtx_status_v = NULL;// 0 if bad vertex, 1 if good
// for location objects:
//  2D: [time, uwire, vwire, ywire]
//  3D: [X,Y,Z]
std::vector<double>* true_vtx_location_3D = NULL; //location of true vtx in event;
std::vector<int>* true_vtx_location_2D = NULL; //location of true vtx in event;
std::vector<std::vector<double>>* reco_vtx_location_3D_v = NULL; //vector of locations of reco vtx
std::vector<std::vector<int>>* reco_vtx_location_2D_v = NULL; //vector of locations of reco vtx

void read_cutvariables(){
  // load file, values, and histograms
  fin = new TFile("hadd_NCPi0CutVariables.root","read");
  tree = (TTree*)fin->Get("ncpi0cutvariable");

  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("subrun",&subrun);
  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("true_track_id_v",&true_track_id_v);
  tree->SetBranchAddress("R_proton_v",&R_proton_v);
  tree->SetBranchAddress("R_gamma_v",&R_gamma_v);
  tree->SetBranchAddress("SSNet_shower_frac_u_v",&SSNet_shower_frac_u_v);
  tree->SetBranchAddress("SSNet_track_frac_u_v",&SSNet_track_frac_u_v);
  tree->SetBranchAddress("SSNet_shower_frac_v_v",&SSNet_shower_frac_v_v);
  tree->SetBranchAddress("SSNet_track_frac_v_v",&SSNet_track_frac_v_v);
  tree->SetBranchAddress("SSNet_shower_frac_y_v",&SSNet_shower_frac_y_v);
  tree->SetBranchAddress("SSNet_track_frac_y_v",&SSNet_track_frac_y_v);
  tree->SetBranchAddress("tracklength_v",&tracklength_v);
  tree->SetBranchAddress("max_dqdx_v",&max_dqdx_v);
  tree->SetBranchAddress("avg_dqdx_v",&avg_dqdx_v);
  tree->SetBranchAddress("vtx_cont_v",&vtx_cont_v);
  tree->SetBranchAddress("vtx_true_fid",&vtx_true_fid);
  tree->SetBranchAddress("vtx_reco_fid_v",&vtx_reco_fid_v);
  tree->SetBranchAddress("vtx_status_v",&vtx_status_v);
  tree->SetBranchAddress("true_vtx_location_3D",&true_vtx_location_3D);
  tree->SetBranchAddress("true_vtx_location_2D",&true_vtx_location_2D);
  tree->SetBranchAddress("reco_vtx_location_3D_v",&reco_vtx_location_3D_v);
  tree->SetBranchAddress("reco_vtx_location_2D_v",&reco_vtx_location_2D_v);

  //initialize histograms to fill with variables
  // just start with R values
  TH1F* R_value_proton = new TH1F("RValProton","RValProton",150,-2,2);
  TH1F* R_value_muon_p = new TH1F("RValMuon_p","RValMuon_p",150,-2,2);
	TH1F* R_value_gamma = new TH1F("RValGamma","RValGamma",150,-15,15);
  TH1F* R_value_muon_g = new TH1F("RValMuon_g","RValMuon_g",150,-15,15);
  TH2F* R_2D_proton = new TH2F("R2DProton","R2DProton",150,-12,12,150,-3,3);
  TH2F* R_2D_gamma = new TH2F("R2DGamma","R2DGamma",150,-12,12,150,-3,3);
  TH2F* R_2D_muon = new TH2F("R2DMuon","R2DMuon",150,-12,12,150,-3,3);
  TH1F* ParticleID = new TH1F("ParticleID","ParticleID",10,0,10);
  TH1F* UnknownTrackLength = new TH1F("UnknownTrackLength","UnknownTrackLength",100,0,100);

  int Nentries = tree->GetEntries();
  std::cout<<"NENTRIES "<<Nentries<<std::endl;

  for (int i =0; i<Nentries;i++){
    tree->GetEntry(i);
    std::cout<<"ENTRY: "<<i<<std::endl;

    for (int vtx = 0; vtx < R_proton_v->size();vtx++){
      bool tracks_contained = true;
      for (int track = 0; track<R_proton_v->at(vtx).size(); track++){
        //particle id histogram
        ParticleID->Fill(true_track_id_v->at(vtx).at(track));
        if (true_track_id_v->at(vtx).at(track) ==0){
          UnknownTrackLength->Fill(tracklength_v->at(vtx).at(track));
        }
        //cosmic muon
        if (true_track_id_v->at(vtx).at(track) == 1){
          R_value_muon_g->Fill(R_gamma_v->at(vtx).at(track));
          R_value_muon_p->Fill(R_proton_v->at(vtx).at(track));
          R_2D_muon->Fill(R_gamma_v->at(vtx).at(track),R_proton_v->at(vtx).at(track));
        }
        //gamma
        else if (true_track_id_v->at(vtx).at(track) == 4){
          R_value_gamma->Fill(R_gamma_v->at(vtx).at(track));
          R_2D_gamma->Fill(R_gamma_v->at(vtx).at(track),R_proton_v->at(vtx).at(track));
        }
        //proton
        else if (true_track_id_v->at(vtx).at(track) == 9){
          R_value_proton->Fill(R_proton_v->at(vtx).at(track));
          R_2D_proton->Fill(R_gamma_v->at(vtx).at(track),R_proton_v->at(vtx).at(track));
        }

        // //get values for cuts
        // if (vtx_cont_v->vtx())
      }
    }
  }//end of entry loop


  // //calculate cut point (intersection)
  // int nbins = R_value_muon->GetSize();
  // float difference = 1000000;
  // int int_bin;
  // float max_muon = 0;
  // int max_muon_bin = 0;
  // float max_proton = 0;
  // int max_proton_bin = 0;
  // for (int ii =1;ii<nbins-1;ii++){
  //   float muonval = R_value_muon->GetBinContent(ii);
  //   float protonval = R_value_proton->GetBinContent(ii);
  //   if (muonval>max_muon){
  //     max_muon = muonval;
  //     max_muon_bin = ii;
  //   }
  //   if (protonval>max_proton){
  //     max_proton = protonval;
  //     max_proton_bin = ii;
  //   }
  // }
  // std::cout<<"Max muon value: "<<max_muon <<" @ bin: "<<max_muon_bin<<std::endl;
  // std::cout<<"Max Proton value: "<<max_proton <<" @ bin: "<<max_proton_bin<<std::endl;
  // int minbin;
  // int maxbin;
  // if (max_muon_bin > max_proton_bin){
  //   maxbin = max_muon_bin;
  //   minbin = max_proton_bin;
  // }
  // else if (max_proton_bin > max_muon_bin){
  //   maxbin = max_proton_bin;
  //   minbin = max_muon_bin;
  // }
  // else {
  //   maxbin=max_muon_bin;
  //   minbin=max_muon_bin;
  // }
  //
  // // get intersection point
  // for(int ii = minbin; ii<=maxbin;ii++){
  //   if(fabs(R_value_muon->GetBinContent(ii) - R_value_proton->GetBinContent(ii))<difference){
  //     difference = fabs(R_value_muon->GetBinContent(ii) - R_value_proton->GetBinContent(ii));
  //     int_bin = ii;
  //   }
  // }
  // std::cout<< "Intersection @ bin: "<<int_bin<< " with difference: "<< difference <<std::endl;
  // std::cout<<"Cut value: "<<int_bin*(6.0/200.0)-3<<std::endl;

  gStyle->SetOptStat(0);
  //particle id histogram
  TCanvas can("can","histograms",600,600);
  ParticleID->SetXTitle("True Particle ID");
  can.SetLogy();
  ParticleID->Draw("HIST");
  can.SaveAs("Images/ParticleID.png");

  TCanvas can1("can","histograms",600,600);
  UnknownTrackLength->Draw("HIST");
  can1.SaveAs("Images/UnknownTrackLength.png");

  //1d R histograms
  TCanvas can3("can","histograms",600,600);
  R_value_muon_p ->SetLineColor(kBlue);
  R_value_proton->SetLineColor(kRed);
  R_value_muon_p->Scale(1/R_value_muon_p->Integral());
  R_value_proton->Scale(1/R_value_proton->Integral());
  R_value_muon_p->Draw("HIST");
  R_value_proton->Draw("HIST SAME");
  can3.SaveAs("Images/LogLikelihoodProton.png");

  TCanvas can4("can","histograms",600,600);
  R_value_muon_g ->SetLineColor(kBlue);
  R_value_gamma->SetLineColor(kGreen);
  R_value_muon_g->Scale(1/R_value_muon_g->Integral());
  R_value_gamma->Scale(1/R_value_gamma->Integral());
  R_value_muon_g->Draw("HIST");
  R_value_gamma->Draw("HIST SAME");
  can4.SaveAs("Images/LogLikelihoodGamma.png");

  //2d R plots
  R_2D_muon->Rebin2D(2,2);
  R_2D_proton->Rebin2D(2,2);
  R_2D_gamma->Rebin2D(2,2);

  TCanvas can5("can","histograms",600,600);
  R_2D_muon->SetXTitle("R(Gamma)");
  R_2D_muon->SetYTitle("R(Proton)");
  R_2D_muon->Draw("COLZ");
  can5.SaveAs("Images/2DMuon.png");

  TCanvas can6("can","histograms",600,600);
  R_2D_gamma->SetXTitle("R(Gamma)");
  R_2D_gamma->SetYTitle("R(Proton)");
  R_2D_gamma->Draw("COLZ");
  can6.SaveAs("Images/2DGamma.png");

  TCanvas can7("can","histograms",600,600);
  R_2D_proton->SetXTitle("R(Gamma)");
  R_2D_proton->SetYTitle("R(Proton)");
  R_2D_proton->Draw("COLZ");
  can7.SaveAs("Images/2DProton.png");

  // //save Aval plot
  // // a_val_muon->Rebin(15);
  // // a_val_proton->Rebin(15);
  // // a_val_gamma->Rebin(15);
  // TCanvas can9("can","histograms",600,600);
  // a_val_muon ->SetLineColor(kBlue);
  // a_val_gamma->SetLineColor(kGreen);
  // a_val_proton->SetLineColor(kRed);
  // a_val_muon->Draw("HIST");
  // a_val_gamma->Draw("HIST SAME");
  // a_val_proton->Draw("HIST SAME");
  // can9.SaveAs("Images/AVal.png");
  //
  // //loop to make eff vs rej plot
  // float Efficiency_Proton;
  // float Efficiency_Gamma;
  // float Rejection;
  // float total_muon = a_val_muon->Integral();
  // float total_gamma = a_val_gamma->Integral();
  // float total_proton = a_val_proton->Integral();
  // float gamma_sofar = 0;
  // float proton_sofar = 0;
  // float muon_sofar = 0;
  // TH2F* Eff_Rej_gamma = new TH2F("EffRejGamma","EffRejGamma",100,0,1,100,0,1);
  // TH2F* Eff_Rej_proton = new TH2F("EffRejProton","EffRejProton",100,0,1,100,0,1);
  // for(int ii = 1; ii<=150;ii++){
  //   gamma_sofar+=a_val_gamma->GetBinContent(ii);
  //   muon_sofar+=a_val_muon->GetBinContent(ii);
  //   proton_sofar+=a_val_proton->GetBinContent(ii);
  //   Efficiency_Proton = ((total_proton - proton_sofar))/(total_proton);
  //   Efficiency_Gamma = ((total_gamma - gamma_sofar))/(total_gamma);
  //   Rejection = (muon_sofar)/total_muon;
  //   // std::cout<<"Efficiency: "<<Efficiency<<std::endl;
  //   // std::cout<<"Rejection: "<<Rejection<<std::endl;
  //   Eff_Rej_proton->Fill(Efficiency_Proton,Rejection,ii);
  //   Eff_Rej_gamma->Fill(Efficiency_Gamma,Rejection,ii);
  // }
  // TCanvas can10("can","histograms",600,600);
  // Eff_Rej_proton->SetXTitle("Efficiency");
  // Eff_Rej_proton->SetYTitle("Rejection");
  // Eff_Rej_proton->SetMarkerStyle(kFullCrossX);
  // Eff_Rej_proton->SetMarkerColor(kRed);
  // Eff_Rej_gamma->SetMarkerStyle(kFullCrossX);
  // Eff_Rej_gamma->SetMarkerColor(kGreen);
  // Eff_Rej_proton->Draw("COLZ");
  // // Eff_Rej_gamma->Draw("COLZ SAME");
  // can10.SaveAs("Images/EffRej.png");

  //clear pointers
  delete true_track_id_v; //vector storing track id for reco tracks for each vtx
  delete R_proton_v; //vector of R(Proton) for recotracks for each vtx
  delete R_gamma_v; //vector of R(Proton) for recotracks for each vtx
  delete SSNet_shower_frac_u_v;//vector of SSNet(uplane) shower fraction for recotracks for each vtx
  delete SSNet_track_frac_u_v;//vector of SSNet(uplane) shower fraction for recotracks for each vtx
  delete SSNet_shower_frac_v_v;//vector of SSNet(vplane) shower fraction for recotracks for each vtx
  delete SSNet_track_frac_v_v;//vector of SSNet(vplane) shower fraction for recotracks for each vtx
  delete SSNet_shower_frac_y_v;//vector of SSNet(yplane) shower fraction for recotracks for each vtx
  delete SSNet_track_frac_y_v;//vector of SSNet(yplane) shower fraction for recotracks for each vtx
  delete tracklength_v;//vector of track length of recotracks for eachvtx;
  delete avg_dqdx_v;// vector of avg dqdx of recotracks for each vtx;
  delete max_dqdx_v;// vector of max dqdx of recotracks for each vtx;
  delete vtx_cont_v;//vector of boolians of is track contained?
  delete vtx_reco_fid_v;//vector of boolians on whether vtx is in fid
  delete vtx_status_v;// 0 if bad vertex, 1 if good
  delete true_vtx_location_3D; //location of true vtx in event;
  delete true_vtx_location_2D; //location of true vtx in event;
  delete reco_vtx_location_3D_v; //vector of locations of reco vtx
  delete reco_vtx_location_2D_v; //vector of locations of reco vtx




}
