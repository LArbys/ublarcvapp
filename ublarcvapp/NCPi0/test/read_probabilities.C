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

//Macro to read in info from root files, save histograms, and determine Values
//made to be used with ncpi0

TFile* fin;
TTree* tree;

// Initalize intger values from file


void read_probabilities(){
  // load file, values, and histograms
  fin = new TFile("hadd2_NCPi0Probabilities.root","read");

  //Set Histograms
  TH2F* Muon_dqdx_resrange_h   = (TH2F*)fin->Get("MuonDQDX");
  TH2F* Proton_dqdx_resrange_h = (TH2F*)fin->Get("ProtonDQDX");
  TH1F* Muon_length_prob_h = (TH1F*)fin->Get("MuonLength");
  TH1F* Gamma_length_prob_h = (TH1F*)fin->Get("GammaLength");

  TH1F* Gamma_ssnetu_h = (TH1F*)fin->Get("GammaSSNetU");
  TH1F* Gamma_ssnetv_h = (TH1F*)fin->Get("GammaSSNetV");
  TH1F* Gamma_ssnety_h = (TH1F*)fin->Get("GammaSSNetY");
  TH1F* Muon_ssnetu_h = (TH1F*)fin->Get("MuonSSNetU");
  TH1F* Muon_ssnetv_h = (TH1F*)fin->Get("MuonSSNetV");
  TH1F* Muon_ssnety_h = (TH1F*)fin->Get("MuonSSNetY");

  Muon_dqdx_resrange_h->SetMarkerColor(kBlue);
  Proton_dqdx_resrange_h->SetMarkerColor(kRed);
  Muon_dqdx_resrange_h->SetMarkerStyle(kDot);
  Proton_dqdx_resrange_h->SetMarkerStyle(kDot);
  Proton_dqdx_resrange_h->Rebin2D(3,3);
  Muon_dqdx_resrange_h->Rebin2D(3,3);

  // //make resrange
  TCanvas can("can","histograms",600,600);
  Proton_dqdx_resrange_h->Draw("COLZ");
  Proton_dqdx_resrange_h->GetYaxis()->SetRangeUser(1,100);
  can.SaveAs("Images/ResRange_Proton.png");

  TCanvas can2("can","histograms",600,600);
  Muon_dqdx_resrange_h->Draw("COLZ");
  Muon_dqdx_resrange_h->GetYaxis()->SetRangeUser(1,100);
  can2.SaveAs("Images/ResRange_Muon.png");

  gStyle->SetOptStat(0);

  TCanvas can8("can","histograms",600,600);
  Muon_length_prob_h->SetLineColor(kBlue);
  Muon_length_prob_h->Scale(1/Muon_length_prob_h->Integral());
  Muon_length_prob_h->Draw("HIST");
  Gamma_length_prob_h->SetLineColor(kRed);
  Gamma_length_prob_h->Scale(1/Gamma_length_prob_h->Integral());
  Gamma_length_prob_h->Draw("HIST");
  Muon_length_prob_h->Draw("HIST SAME");
  can8.SaveAs("Images/Tracklength.png");

  TCanvas can9("can","histograms",600,600);
  Muon_ssnetu_h->SetLineColor(kBlue);
  Muon_ssnetu_h->Scale(1/Muon_ssnetu_h->Integral());
  Muon_ssnetu_h->Draw("HIST");
  Gamma_ssnetu_h->SetLineColor(kRed);
  Gamma_ssnetu_h->Scale(1/Gamma_ssnetu_h->Integral());
  Gamma_ssnetu_h->Draw("HIST");
  Muon_ssnetu_h->Draw("HIST SAME");
  can9.SaveAs("Images/SSnetU.png");

  TCanvas can10("can","histograms",600,600);
  Muon_ssnetv_h->SetLineColor(kBlue);
  Muon_ssnetv_h->Scale(1/Muon_ssnetv_h->Integral());
  Muon_ssnetv_h->Draw("HIST");
  Gamma_ssnetv_h->SetLineColor(kRed);
  Gamma_ssnetv_h->Scale(1/Gamma_ssnetv_h->Integral());
  Gamma_ssnetv_h->Draw("HIST");
  Muon_ssnetv_h->Draw("HIST SAME");
  can10.SaveAs("Images/SSnetV.png");

  TCanvas can11("can","histograms",600,600);
  Muon_ssnety_h->SetLineColor(kBlue);
  Muon_ssnety_h->Scale(1/Muon_ssnety_h->Integral());
  Muon_ssnety_h->Draw("HIST");
  Gamma_ssnety_h->SetLineColor(kRed);
  Gamma_ssnety_h->Scale(1/Gamma_ssnety_h->Integral());
  Gamma_ssnety_h->Draw("HIST");
  Muon_ssnety_h->Draw("HIST SAME");
  can11.SaveAs("Images/SSnetY.png");

  //making a for loop for slices
  int numslices = 20;
  bool save_prob = true;
  TH1D* Muon_resrange_slice_v[numslices];
  TH1D* Proton_resrange_slice_v[numslices];
  for (int ii = 0;ii<numslices;ii++){
    Muon_resrange_slice_v[ii] = Muon_dqdx_resrange_h->ProjectionY(Form("muon_%d",ii),ii,ii+1);
    Proton_resrange_slice_v[ii] = Proton_dqdx_resrange_h->ProjectionY(Form("proton_%d",ii),ii,ii+1);
    // change each to a probability by normalizing
    double muon_integral = Muon_resrange_slice_v[ii]->Integral();
    double proton_integral = Proton_resrange_slice_v[ii]->Integral();
    Muon_resrange_slice_v[ii]->Scale(1/muon_integral);
    Proton_resrange_slice_v[ii]->Scale(1/proton_integral);
    if (save_prob){
      TCanvas can("can","histograms",600,600);
      can.SetTitle(Form("slice_%d",ii));
      Muon_resrange_slice_v[ii]->SetLineColor(kBlue);
      Proton_resrange_slice_v[ii]->SetLineColor(kRed);
      Muon_resrange_slice_v[ii]->Draw("HIST");
      Proton_resrange_slice_v[ii]->Draw("HIST SAME");
      TLegend legend = TLegend(.6,.9,.6,.9);
      legend.AddEntry(Muon_resrange_slice_v[ii],"Muon","l");
      legend.AddEntry(Proton_resrange_slice_v[ii],"Proton","l");
      legend.Draw();
      can.SaveAs(Form("Images/ResRange_Slice%d.png",ii));
    }

  }



}
