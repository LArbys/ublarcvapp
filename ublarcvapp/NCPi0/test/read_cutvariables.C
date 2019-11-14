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
std::vector<std::vector<float> >* SSNet_shower_frac_u_v = NULL;//vector of SSNet(uplane) shower fraction for recotracks for each vtx
std::vector<std::vector<float> >* SSNet_track_frac_u_v = NULL;//vector of SSNet(uplane) shower fraction for recotracks for each vtx
std::vector<std::vector<float> >* SSNet_shower_frac_v_v = NULL;//vector of SSNet(vplane) shower fraction for recotracks for each vtx
std::vector<std::vector<float> >* SSNet_track_frac_v_v = NULL;//vector of SSNet(vplane) shower fraction for recotracks for each vtx
std::vector<std::vector<float> >* SSNet_shower_frac_y_v = NULL;//vector of SSNet(yplane) shower fraction for recotracks for each vtx
std::vector<std::vector<float> >* SSNet_track_frac_y_v = NULL;//vector of SSNet(yplane) shower fraction for recotracks for each vtx
std::vector<std::vector<float> >* tracklength_v = NULL;//vector of track length of recotracks for eachvtx;
std::vector<std::vector<float> >* avg_dqdx_v = NULL;// vector of avg dqdx of recotracks for each vtx;
std::vector<std::vector<float> >* max_dqdx_v = NULL;// vector of max dqdx of recotracks for each vtx;
std::vector<std::vector<bool> >* vtx_cont_v = NULL;//vector of boolians of is track contained?
bool vtx_true_fid;//vector of boolians on whether vtx is in fid
std::vector<bool>* vtx_reco_fid_v = NULL;//vector of boolians on whether vtx is in fid
std::vector<bool>* vtx_status_v = NULL;// 0 if bad vertex, 1 if good
double pot;
// for location objects:
//  2D: [time, uwire, vwire, ywire]
//  3D: [X,Y,Z]
std::vector<double>* true_vtx_location_3D = NULL; //location of true vtx in event;
std::vector<int>* true_vtx_location_2D = NULL; //location of true vtx in event;
std::vector<std::vector<double> >* reco_vtx_location_3D_v = NULL; //vector of locations of reco vtx
std::vector<std::vector<int> >* reco_vtx_location_2D_v = NULL; //vector of locations of reco vtx
// shower variables
int showerreco_nshowers;
std::vector<std::vector<std::vector<int> > >* ShowerAssociation_vvv = NULL;//object that contains all of the shower associations
std::vector<std::vector<int> >* showerid_v = NULL; //vector of ids of all showers for each vertex
std::vector<std::vector<double> >* shower_reco_totalE_v = NULL;//total energy of each reco shower
std::vector<std::vector<double> >* shower_reco_length_v = NULL;//length of shower
std::vector<std::vector<double> >* shower_reco_openingangle_v = NULL;//opening angle of shower
std::vector<std::vector<double> >* shower_reco_ssnetshower_v = NULL; //shower fraction of hits
std::vector<std::vector<std::vector<double> > >* shower_reco_dedx_vv = NULL; //vector of dedx for each reco shower object

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
  tree->SetBranchAddress("showerreco_nshowers",&showerreco_nshowers);
  tree->SetBranchAddress("pot",&pot);
  tree->SetBranchAddress("ShowerAssociation_vvv",&ShowerAssociation_vvv);
  tree->SetBranchAddress("showerid_v",&showerid_v);
  tree->SetBranchAddress("shower_reco_totalE_v",&shower_reco_totalE_v);
  tree->SetBranchAddress("shower_reco_length_v",&shower_reco_length_v);
  tree->SetBranchAddress("shower_reco_openingangle_v",&shower_reco_openingangle_v);
  tree->SetBranchAddress("shower_reco_ssnetshower_v",&shower_reco_ssnetshower_v);
  tree->SetBranchAddress("shower_reco_dedx_vv",&shower_reco_dedx_vv);

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

  TH1F* AvgDQDX_proton = new TH1F("AvgDQDX_proton","AvgDQDX_proton",100,0,100);
  TH1F* AvgDQDX_muon = new TH1F("AvgDQDX_muon","AvgDQDX_muon",100,0,100);
  TH1F* AvgDQDX_gamma = new TH1F("AvgDQDX_gamma","AvgDQDX_gamma",100,0,100);
  TH1F* MaxDQDX_proton = new TH1F("MaxDQDX_proton","MaxDQDX_proton",100,0,200);
  TH1F* MaxDQDX_muon = new TH1F("MaxDQDX_muon","MaxDQDX_muon",100,0,200);
  TH1F* MaxDQDX_gamma = new TH1F("MaxDQDX_gamma","MaxDQDX_gamma",100,0,200);

  TH2F* Efficiency = new TH2F("Efficiency","Efficiency",60,-15,15,40,-2,2);
  TH2F* Rejection = new TH2F("Rejection","Rejection",40,-5,15,30,-1,2);

  //shower histograms
  TH1F* ShowerID_good_h = new TH1F("ShowerIDGood","ShowerIDGood",10,0,10);
  TH1F* ShowerID_h = new TH1F("ShowerID","ShowerID",10,0,10);
  //oh boy here we go......
  TH1F* showerlength_unknown_h = new TH1F("showerlengthunknown","showerlengthunknown",100,0,150);
  TH1F* showerlength_cosmic_h = new TH1F("showerlengthcosmic","showerlengthcosmic",100,0,150);
  TH1F* showerlength_electron_h = new TH1F("showerlengthelectron","showerlengthelectron",100,0,150);
  TH1F* showerlength_gamma_h = new TH1F("showerlengthgamma","showerlengthgamma",100,0,150);
  TH1F* showerlength_pion_h = new TH1F("showerlengthpion","showerlengpion",100,0,150);
  TH1F* showerlength_proton_h = new TH1F("showerlengthproton","showerlengthproton",100,0,150);

  TH1F* showerenergy_unknown_h = new TH1F("showerenergyunknown","showerenergyunknown",100,0,150);
  TH1F* showerenergy_cosmic_h = new TH1F("showerenergycosmic","showerenergycosmic",100,0,150);
  TH1F* showerenergy_electron_h = new TH1F("showerenergyelectron","showerenergyelectron",100,0,150);
  TH1F* showerenergy_gamma_h = new TH1F("showerenergygamma","showerenergygamma",100,0,150);
  TH1F* showerenergy_pion_h = new TH1F("showerenergypion","showerlengpion",100,0,150);
  TH1F* showerenergy_proton_h = new TH1F("showerenergyproton","showerenergyproton",100,0,150);

  TH1F* showeropeningangle_unknown_h = new TH1F("showeropeningangleunknown","showeropeningangleunknown",100,0,3.14);
  TH1F* showeropeningangle_cosmic_h = new TH1F("showeropeninganglecosmic","showeropeninganglecosmic",100,0,3.14);
  TH1F* showeropeningangle_electron_h = new TH1F("showeropeningangleelectron","showeropeningangleelectron",100,0,3.14);
  TH1F* showeropeningangle_gamma_h = new TH1F("showeropeninganglegamma","showeropeninganglegamma",100,0,3.14);
  TH1F* showeropeningangle_pion_h = new TH1F("showeropeninganglepion","showerlengpion",100,0,3.14);
  TH1F* showeropeningangle_proton_h = new TH1F("showeropeningangleproton","showeropeningangleproton",100,0,3.14);

  TH1F* showerssnetshower_unknown_h = new TH1F("showerssnetshowerunknown","showerssnetshowerunknown",25,0,1.00001);
  TH1F* showerssnetshower_cosmic_h = new TH1F("showerssnetshowercosmic","showerssnetshowercosmic",25,0,1.00001);
  TH1F* showerssnetshower_electron_h = new TH1F("showerssnetshowerelectron","showerssnetshowerelectron",25,0,1.00001);
  TH1F* showerssnetshower_gamma_h = new TH1F("showerssnetshowergamma","showerssnetshowergamma",25,0,1.00001);
  TH1F* showerssnetshower_pion_h = new TH1F("showerssnetshowerpion","showerlengpion",25,0,1.00001);
  TH1F* showerssnetshower_proton_h = new TH1F("showerssnetshowerproton","showerssnetshowerproton",25,0,1.00001);

  TH2F* showerdqdx_muon_h = new TH2F("MuonDQDX","MuonDQDX",20,0,10,100,0,100);

  int Nentries = tree->GetEntries();
  std::cout<<"NENTRIES "<<Nentries<<std::endl;

  //initialize totals from cuts
  int total_good = 0;
  int total_bad = 0;
  int total_good_fid = 0;
  int total_bad_fid = 0;
  int total_good_contained = 0;
  int total_bad_contained = 0;
  int total_good_showerqualitycut = 0;
  int total_bad_showerqualitycut = 0;
  int total_good_showerssnet = 0;
  int total_bad_showerssnet = 0;
  int total_good_Rcut = 0;
  int total_bad_Rcut = 0;
  int total_truefid = 0;
  double totpot = 0;

  for (int i =0; i<Nentries;i++){
    tree->GetEntry(i);
    totpot +=pot;
    std::cout<<"ENTRY: "<<i<<std::endl;
    if (vtx_true_fid) total_truefid++;
    for (int vtx = 0; vtx < R_proton_v->size();vtx++){
      //add to good and bad total
      if (vtx_status_v->at(vtx)==1) total_good ++;
      else total_bad++;
      //initialize boolians for cuts
      bool AllTracksContained = true;
      bool VtxInFiducial = vtx_reco_fid_v->at(vtx);
      bool Passes2DRCut = false;
      bool PassesGammaCut = false;
      bool PassesProtonCut = false;
      bool PassesShowerQualCut = false;
      bool PassesShowerSSnetCut = false;

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
          AvgDQDX_muon->Fill(avg_dqdx_v->at(vtx).at(track));
          MaxDQDX_muon->Fill(max_dqdx_v->at(vtx).at(track));
        }
        //gamma
        else if (true_track_id_v->at(vtx).at(track) == 4){
          R_value_gamma->Fill(R_gamma_v->at(vtx).at(track));
          R_2D_gamma->Fill(R_gamma_v->at(vtx).at(track),R_proton_v->at(vtx).at(track));
          AvgDQDX_gamma->Fill(avg_dqdx_v->at(vtx).at(track));
          MaxDQDX_gamma->Fill(max_dqdx_v->at(vtx).at(track));
        }
        //proton
        else if (true_track_id_v->at(vtx).at(track) == 9){
          R_value_proton->Fill(R_proton_v->at(vtx).at(track));
          R_2D_proton->Fill(R_gamma_v->at(vtx).at(track),R_proton_v->at(vtx).at(track));
          AvgDQDX_proton->Fill(avg_dqdx_v->at(vtx).at(track));
          MaxDQDX_proton->Fill(max_dqdx_v->at(vtx).at(track));
        }

        //Does it pass R Cut?
        if (R_proton_v->at(vtx).at(track)>1) PassesProtonCut =true;
        if (R_gamma_v->at(vtx).at(track)>2) PassesGammaCut = true;


      }//end of track loop

      //loop through showers
      for(int shower = 0; shower<showerid_v->at(vtx).size();shower++){
        int showerid = showerid_v->at(vtx).at(shower);
        double length = shower_reco_length_v->at(vtx).at(shower);
        double energy = shower_reco_totalE_v->at(vtx).at(shower);
        double ssnetshower = shower_reco_ssnetshower_v->at(vtx).at(shower);
        //make sure not garbage shower
        if(ssnetshower>0&&length>0&&energy>0){
          PassesShowerQualCut = true;
        }
        if (ssnetshower>.5) PassesShowerSSnetCut = true;
      }//end of shower loop

      //figure out if cut
      if(PassesProtonCut||PassesGammaCut)Passes2DRCut=true;
      if (vtx_status_v->at(vtx) ==1){
        if(VtxInFiducial){
          total_good_fid++;
          if(AllTracksContained){
            total_good_contained++;
            if(Passes2DRCut){
              total_good_Rcut++;
              if(PassesShowerQualCut){
                total_good_showerqualitycut++;
                if(PassesShowerSSnetCut){
                  total_good_showerssnet++;
                }
              }
            }
          }
        }
      }
      else{
        if(VtxInFiducial){
          total_bad_fid++;
          if(AllTracksContained){
            total_bad_contained++;
            if(Passes2DRCut){
              total_bad_Rcut++;
              if(PassesShowerQualCut){
                total_bad_showerqualitycut++;
                if(PassesShowerSSnetCut){
                  total_bad_showerssnet++;
                }
              }
            }
          }
        }
      }
    }//end of vertex loop
    //showers loop
    // std::cout<<shower_reco_length_v->size()<<std::endl;
    for(int vtx = 0; vtx<showerid_v->size(); vtx++){
      //loop through Showers
      for(int shower = 0; shower<showerid_v->at(vtx).size();shower++){
        int showerid = showerid_v->at(vtx).at(shower);
        double length = shower_reco_length_v->at(vtx).at(shower);
        double energy = shower_reco_totalE_v->at(vtx).at(shower);
        double ssnetshower = shower_reco_ssnetshower_v->at(vtx).at(shower);
        double openingangle = shower_reco_openingangle_v->at(vtx).at(shower);

        //make sure not garbage shower
        if(ssnetshower>0&&length>0&&energy>0){
          ShowerID_h->Fill(showerid);
          if(vtx_status_v->at(vtx)==1){
            ShowerID_good_h->Fill(showerid);
          }
          //sort by particle type:
          if(showerid ==0){//unknown
            showerlength_unknown_h->Fill(length);
            showerenergy_unknown_h->Fill(energy);
            showeropeningangle_unknown_h->Fill(openingangle);
            showerssnetshower_unknown_h->Fill(ssnetshower);
          }
          else if(showerid ==1){//cosmic
            showerlength_cosmic_h->Fill(length);
            showerenergy_cosmic_h->Fill(energy);
            showeropeningangle_cosmic_h->Fill(openingangle);
            showerssnetshower_cosmic_h->Fill(ssnetshower);
            //dqdx values
            for(int ii = 0; ii< shower_reco_dedx_vv->at(vtx).at(shower).size();ii++ ){
              float dqdx =  shower_reco_dedx_vv->at(vtx).at(shower).at(ii);
              showerdqdx_muon_h->Fill(ii,dqdx,1);
            }

          }
          else if(showerid ==3){//electron
            showerlength_electron_h->Fill(length);
            showerenergy_electron_h->Fill(energy);
            showeropeningangle_electron_h->Fill(openingangle);
            showerssnetshower_electron_h->Fill(ssnetshower);
          }
          else if(showerid ==4){//gamma
            showerlength_gamma_h->Fill(length);
            showerenergy_gamma_h->Fill(energy);
            showeropeningangle_gamma_h->Fill(openingangle);
            showerssnetshower_gamma_h->Fill(ssnetshower);
          }
          else if(showerid ==8){//charged pion
            showerlength_pion_h->Fill(length);
            showerenergy_pion_h->Fill(energy);
            showeropeningangle_pion_h->Fill(openingangle);
            showerssnetshower_pion_h->Fill(ssnetshower);
          }
          else if(showerid ==9){//proton
            showerlength_proton_h->Fill(length);
            showerenergy_proton_h->Fill(energy);
            showeropeningangle_proton_h->Fill(openingangle);
            showerssnetshower_proton_h->Fill(ssnetshower);
          }

        }
      }//end of shower loop
    }//end of shower vtx loop
  }//end of entry loop

  // loop through RValues
  // for(float rgamma = -5.0; rgamma<15;rgamma+=.5){
  //   for (float rproton = -1; rproton< 2;rproton+=.1){
  //     std::cout<<"r gamma: "<<rgamma<<" r proton: "<<rproton<< std::endl;
  //     float eff_val;
  //     float rej_val;
  //     float total_good_pass = 0;
  //     float total_bad_pass = 0;
  //     //loop through entries... very time consuming....
  //     for (int i =0; i<Nentries;i++){
  //       tree->GetEntry(i);
  //       //loop through vertices
  //       for (int vtx = 0; vtx < R_proton_v->size();vtx++){
  //         bool VtxPasses = false;
  //         //turns to true if >0 tracks passes;
  //         //loop through tracks
  //         for (int track = 0; track<R_proton_v->at(vtx).size(); track++){
  //           if(R_proton_v->at(vtx).at(track)>rproton) VtxPasses =true;
  //           if (R_gamma_v->at(vtx).at(track)>rgamma) VtxPasses = true;
  //         }
  //         if (vtx_status_v->at(vtx)==1){
  //           if (VtxPasses) total_good_pass++;
  //         }
  //         else{
  //           if (VtxPasses) total_bad_pass++;
  //         }
  //       }
  //     }
  //     eff_val = (float)total_good_pass/(float)total_good;
  //     rej_val = (float)(total_bad-total_bad_pass)/(float)total_bad;
  //     Efficiency->Fill(rgamma,rproton,eff_val);
  //     Rejection->Fill(rgamma,rproton,rej_val);
  //   }
  // }
  gStyle->SetOptStat(0);
  //save Eff and rej plots
  // TCanvas caneff("can","histograms",600,600);
  // Efficiency->SetXTitle("R(Gamma)");
  // Efficiency->SetYTitle("R(Proton)");
  // Efficiency->Draw("COLZ");
  // caneff.SaveAs("Images/Efficiency.png");
  //
  // TCanvas canrej("can","histograms",600,600);
  // Rejection->SetXTitle("R(Gamma)");
  // Rejection->SetYTitle("R(Proton)");
  // Rejection->Draw("COLZ");
  // canrej.SaveAs("Images/Rejection.png");

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

  //particle id histogram
  // TCanvas can("can","histograms",600,600);
  // ParticleID->SetXTitle("True Particle ID");
  // can.SetLogy();
  // ParticleID->Draw("HIST");
  // can.SaveAs("Images/ParticleID.png");
  //
  // TCanvas can1("can","histograms",600,600);
  // UnknownTrackLength->Draw("HIST");
  // can1.SaveAs("Images/UnknownTrackLength.png");
  //
  // //1d R histograms
  // TCanvas can3("can","histograms",600,600);
  // R_value_muon_p ->SetLineColor(kBlue);
  // R_value_proton->SetLineColor(kRed);
  // R_value_muon_p->Scale(1/R_value_muon_p->Integral());
  // R_value_proton->Scale(1/R_value_proton->Integral());
  // R_value_muon_p->Draw("HIST");
  // R_value_proton->Draw("HIST SAME");
  // can3.SaveAs("Images/LogLikelihoodProton.png");
  //
  // TCanvas can4("can","histograms",600,600);
  // R_value_muon_g ->SetLineColor(kBlue);
  // R_value_gamma->SetLineColor(kGreen);
  // R_value_muon_g->Scale(1/R_value_muon_g->Integral());
  // R_value_gamma->Scale(1/R_value_gamma->Integral());
  // R_value_muon_g->Draw("HIST");
  // R_value_gamma->Draw("HIST SAME");
  // can4.SaveAs("Images/LogLikelihoodGamma.png");
  //
  // //2d R plots
  // R_2D_muon->Rebin2D(2,2);
  // R_2D_proton->Rebin2D(2,2);
  // R_2D_gamma->Rebin2D(2,2);
  //
  // TCanvas can5("can","histograms",600,600);
  // R_2D_muon->SetXTitle("R(Gamma)");
  // R_2D_muon->SetYTitle("R(Proton)");
  // R_2D_muon->Draw("COLZ");
  // can5.SaveAs("Images/2DMuon.png");
  //
  // TCanvas can6("can","histograms",600,600);
  // R_2D_gamma->SetXTitle("R(Gamma)");
  // R_2D_gamma->SetYTitle("R(Proton)");
  // R_2D_gamma->Draw("COLZ");
  // can6.SaveAs("Images/2DGamma.png");
  //
  // TCanvas can7("can","histograms",600,600);
  // R_2D_proton->SetXTitle("R(Gamma)");
  // R_2D_proton->SetYTitle("R(Proton)");
  // R_2D_proton->Draw("COLZ");
  // can7.SaveAs("Images/2DProton.png");
  //
  // //save avg and max dqdx plots
  // TCanvas can8("can","histograms",600,600);
  // AvgDQDX_muon->SetXTitle("Avg Dqdx");
  // AvgDQDX_muon->SetLineColor(kBlue);
  // AvgDQDX_proton->SetLineColor(kRed);
  // AvgDQDX_gamma->SetLineColor(kGreen);
  // AvgDQDX_muon->Scale(1.0/AvgDQDX_muon->Integral());
  // AvgDQDX_proton->Scale(1.0/AvgDQDX_proton->Integral());
  // AvgDQDX_gamma->Scale(1.0/AvgDQDX_gamma->Integral());
  // AvgDQDX_muon->Draw("HIST");
  // AvgDQDX_proton->Draw("HIST SAME");
  // AvgDQDX_gamma->Draw("HIST SAME");
  // can8.SaveAs("Images/AvgDQDX.png");
  //
  // TCanvas can9("can","histograms",600,600);
  // MaxDQDX_muon->SetXTitle("Max Dqdx");
  // MaxDQDX_muon->SetLineColor(kBlue);
  // MaxDQDX_proton->SetLineColor(kRed);
  // MaxDQDX_gamma->SetLineColor(kGreen);
  // MaxDQDX_muon->Scale(1.0/MaxDQDX_muon->Integral());
  // MaxDQDX_proton->Scale(1.0/MaxDQDX_proton->Integral());
  // MaxDQDX_gamma->Scale(1.0/MaxDQDX_gamma->Integral());
  // MaxDQDX_muon->Draw("HIST");
  // MaxDQDX_proton->Draw("HIST SAME");
  // MaxDQDX_gamma->Draw("HIST SAME");
  // can9.SaveAs("Images/MaxDQDX.png");
  //
  // //shower Images
  // TCanvas can10("can","histograms",600,600);
  // ShowerID_h->SetXTitle("True Shower ID");
  // can10.SetLogy();
  // ShowerID_h->Draw("HIST");
  // can10.SaveAs("Images/ShowerID.png");
  //
  // TCanvas can11("can","histograms",600,600);
  // ShowerID_good_h->SetXTitle("True Shower ID (and with good vtx)");
  // can10.SetLogy();
  // ShowerID_good_h->Draw("HIST");
  // can11.SaveAs("Images/ShowerID_good.png");
  //
  // TCanvas can12("can","histograms",600,600);
  // showerlength_unknown_h->SetXTitle("Length of shower (cm)");
  // //set colors
  // showerlength_unknown_h->SetLineColor(kBlack);
  // showerlength_cosmic_h->SetLineColor(kBlue);
  // showerlength_electron_h->SetLineColor(kYellow);
  // showerlength_gamma_h->SetLineColor(kGreen);
  // showerlength_pion_h->SetLineColor(kOrange);
  // showerlength_proton_h->SetLineColor(kRed);
  // // normalize to 1
  // showerlength_unknown_h->Scale(1/showerlength_unknown_h->Integral());
  // showerlength_cosmic_h->Scale(1/showerlength_cosmic_h->Integral());
  // showerlength_electron_h->Scale(1/showerlength_electron_h->Integral());
  // showerlength_gamma_h->Scale(1/showerlength_gamma_h->Integral());
  // showerlength_pion_h->Scale(1/showerlength_pion_h->Integral());
  // showerlength_proton_h->Scale(1/showerlength_proton_h->Integral());
  // //draw
  // showerlength_unknown_h->Draw("HIST");
  // showerlength_cosmic_h->Draw("HIST SAME");
  // showerlength_electron_h->Draw("HIST SAME");
  // showerlength_gamma_h->Draw("HIST SAME");
  // showerlength_pion_h->Draw("HIST SAME");
  // showerlength_proton_h->Draw("HIST SAME");
  // can12.SaveAs("Images/ShowerLength.png");
  // TCanvas can13("can","histograms",600,600);
  // showerlength_cosmic_h->SetXTitle("Length of shower (cm)");
  // showerlength_cosmic_h->Draw("HIST SAME");
  // showerlength_gamma_h->Draw("HIST SAME");
  // showerlength_proton_h->Draw("HIST SAME");
  // can13.SaveAs("Images/ShowerLength_simple.png");
  //
  // TCanvas can14("can","histograms",600,600);
  // showerenergy_unknown_h->SetXTitle("Total reco energy of shower (MeV)");
  // //set colors
  // showerenergy_unknown_h->SetLineColor(kBlack);
  // showerenergy_cosmic_h->SetLineColor(kBlue);
  // showerenergy_electron_h->SetLineColor(kYellow);
  // showerenergy_gamma_h->SetLineColor(kGreen);
  // showerenergy_pion_h->SetLineColor(kOrange);
  // showerenergy_proton_h->SetLineColor(kRed);
  // // normalize to 1
  // showerenergy_unknown_h->Scale(1/showerenergy_unknown_h->Integral());
  // showerenergy_cosmic_h->Scale(1/showerenergy_cosmic_h->Integral());
  // showerenergy_electron_h->Scale(1/showerenergy_electron_h->Integral());
  // showerenergy_gamma_h->Scale(1/showerenergy_gamma_h->Integral());
  // showerenergy_pion_h->Scale(1/showerenergy_pion_h->Integral());
  // showerenergy_proton_h->Scale(1/showerenergy_proton_h->Integral());
  // //draw
  // showerenergy_unknown_h->Draw("HIST");
  // showerenergy_cosmic_h->Draw("HIST SAME");
  // showerenergy_electron_h->Draw("HIST SAME");
  // showerenergy_gamma_h->Draw("HIST SAME");
  // showerenergy_pion_h->Draw("HIST SAME");
  // showerenergy_proton_h->Draw("HIST SAME");
  // can14.SaveAs("Images/Showerenergy.png");
  // TCanvas can15("can","histograms",600,600);
  // showerenergy_cosmic_h->SetXTitle("Total reco energy of shower (MeV)");
  // showerenergy_cosmic_h->Draw("HIST SAME");
  // showerenergy_gamma_h->Draw("HIST SAME");
  // showerenergy_proton_h->Draw("HIST SAME");
  // can15.SaveAs("Images/Showerenergy_simple.png");
  //
  // TCanvas can16("can","histograms",600,600);
  // showeropeningangle_unknown_h->SetXTitle("Opening Angle (radians)");
  // //set colors
  // showeropeningangle_unknown_h->SetLineColor(kBlack);
  // showeropeningangle_cosmic_h->SetLineColor(kBlue);
  // showeropeningangle_electron_h->SetLineColor(kYellow);
  // showeropeningangle_gamma_h->SetLineColor(kGreen);
  // showeropeningangle_pion_h->SetLineColor(kOrange);
  // showeropeningangle_proton_h->SetLineColor(kRed);
  // // normalize to 1
  // showeropeningangle_unknown_h->Scale(1/showeropeningangle_unknown_h->Integral());
  // showeropeningangle_cosmic_h->Scale(1/showeropeningangle_cosmic_h->Integral());
  // showeropeningangle_electron_h->Scale(1/showeropeningangle_electron_h->Integral());
  // showeropeningangle_gamma_h->Scale(1/showeropeningangle_gamma_h->Integral());
  // showeropeningangle_pion_h->Scale(1/showeropeningangle_pion_h->Integral());
  // showeropeningangle_proton_h->Scale(1/showeropeningangle_proton_h->Integral());
  // //draw
  // showeropeningangle_unknown_h->Draw("HIST");
  // showeropeningangle_cosmic_h->Draw("HIST SAME");
  // showeropeningangle_electron_h->Draw("HIST SAME");
  // showeropeningangle_gamma_h->Draw("HIST SAME");
  // showeropeningangle_pion_h->Draw("HIST SAME");
  // showeropeningangle_proton_h->Draw("HIST SAME");
  // can16.SaveAs("Images/Showeropeningangle.png");
  // TCanvas can17("can","histograms",600,600);
  // showeropeningangle_cosmic_h->SetXTitle("Opening Angle (radians)");
  // showeropeningangle_cosmic_h->Draw("HIST SAME");
  // showeropeningangle_gamma_h->Draw("HIST SAME");
  // showeropeningangle_proton_h->Draw("HIST SAME");
  // can17.SaveAs("Images/Showeropeningangle_simple.png");
  //
  // TCanvas can18("can","histograms",600,600);
  // showerssnetshower_unknown_h->SetXTitle("SSNet shower fraction");
  // //set colors
  // showerssnetshower_unknown_h->SetLineColor(kBlack);
  // showerssnetshower_cosmic_h->SetLineColor(kBlue);
  // showerssnetshower_electron_h->SetLineColor(kYellow);
  // showerssnetshower_gamma_h->SetLineColor(kGreen);
  // showerssnetshower_pion_h->SetLineColor(kOrange);
  // showerssnetshower_proton_h->SetLineColor(kRed);
  // // normalize to 1
  // showerssnetshower_unknown_h->Scale(1/showerssnetshower_unknown_h->Integral());
  // showerssnetshower_cosmic_h->Scale(1/showerssnetshower_cosmic_h->Integral());
  // showerssnetshower_electron_h->Scale(1/showerssnetshower_electron_h->Integral());
  // showerssnetshower_gamma_h->Scale(1/showerssnetshower_gamma_h->Integral());
  // showerssnetshower_pion_h->Scale(1/showerssnetshower_pion_h->Integral());
  // showerssnetshower_proton_h->Scale(1/showerssnetshower_proton_h->Integral());
  // //draw
  // showerssnetshower_unknown_h->Draw("HIST");
  // showerssnetshower_cosmic_h->Draw("HIST SAME");
  // showerssnetshower_electron_h->Draw("HIST SAME");
  // showerssnetshower_gamma_h->Draw("HIST SAME");
  // showerssnetshower_pion_h->Draw("HIST SAME");
  // showerssnetshower_proton_h->Draw("HIST SAME");
  // can18.SaveAs("Images/Showerssnetshower.png");
  // TCanvas can19("can","histograms",600,600);
  // showerssnetshower_cosmic_h->SetXTitle("SSNet shower fraction");
  // showerssnetshower_cosmic_h->Draw("HIST SAME");
  // showerssnetshower_gamma_h->Draw("HIST SAME");
  // showerssnetshower_proton_h->Draw("HIST SAME");
  // can19.SaveAs("Images/Showerssnetshower_simple.png");
  //
  // TCanvas can20("can","histograms",600,600);
  // showerdqdx_muon_h->Draw("COLZ");
  // can20.SaveAs("Images/showerdqdx_Muon.png");

  //print out results of cuts....
  float scaleratio = (13.0*10e20)/totpot;

  std::cout<<"---------------------------------------------"<<std::endl;
  std::cout<<"totpot: "<<totpot<<std::endl;
  std::cout<<"scalerato: "<<scaleratio<<std::endl;
  std::cout<<"-----------------------------------------------------"<<std::endl;
  std::cout<<"SCALING COUNTS BY POT"<<std::endl;
  std::cout<<"Total Number of Events with true in Fid: "<<total_truefid*scaleratio<<std::endl;
  std::cout<<"Total Number of Good: "<<total_good*scaleratio<<std::endl;
  std::cout<<"Total Number of Bad: "<<total_bad*scaleratio<<std::endl;
  std::cout<<std::endl;
  std::cout<<"AFTER FIDUCIAL VOLUME CUT"<<std::endl;
  std::cout<<"Number Good: "<<total_good_fid*scaleratio<<", Number Bad: "<<total_bad_fid*scaleratio<<std::endl;
  std::cout<<"Efficiency: "<<(float)total_good_fid/(float)total_good<<std::endl;
  std::cout<<"Rejection: "<<(float)(total_bad-total_bad_fid)/(float)total_bad<<std::endl;
  std::cout<<std::endl;
  std::cout<<"AFTER TRACK CONTAINMENT CUT"<<std::endl;
  std::cout<<"Number Good: "<<total_good_contained*scaleratio<<", Number Bad: "<<total_bad_contained*scaleratio<<std::endl;
  std::cout<<"Efficiency: "<<(float)total_good_contained/(float)total_good<<std::endl;
  std::cout<<"Rejection: "<<(float)(total_bad-total_bad_contained)/(float)total_bad<<std::endl;
  std::cout<<std::endl;
  std::cout<<"AFTER 2D R-TRACK CUT"<<std::endl;
  std::cout<<"Number Good: "<<total_good_Rcut*scaleratio<<", Number Bad: "<<total_bad_Rcut*scaleratio<<std::endl;
  std::cout<<"Efficiency: "<<(float)total_good_Rcut/(float)total_good<<std::endl;
  std::cout<<"Rejection: "<<(float)(total_bad-total_bad_Rcut)/(float)total_bad<<std::endl;
  std::cout<<std::endl;
  std::cout<<"AFTER SHOWER QUALITY CUT"<<std::endl;
  std::cout<<"Number Good: "<<total_good_showerqualitycut*scaleratio<<", Number Bad: "<<total_bad_showerqualitycut*scaleratio<<std::endl;
  std::cout<<"Efficiency: "<<(float)total_good_showerqualitycut/(float)total_good<<std::endl;
  std::cout<<"Rejection: "<<(float)(total_bad-total_bad_showerqualitycut)/(float)total_bad<<std::endl;
  std::cout<<std::endl;
  std::cout<<"AFTER SHOWER SSNET CUT"<<std::endl;
  std::cout<<"Number Good: "<<total_good_showerssnet*scaleratio<<", Number Bad: "<<total_bad_showerssnet*scaleratio<<std::endl;
  std::cout<<"Efficiency: "<<(float)total_good_showerssnet/(float)total_good<<std::endl;
  std::cout<<"Rejection: "<<(float)(total_bad-total_bad_showerssnet)/(float)total_bad<<std::endl;

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
  delete ShowerAssociation_vvv;//object that contains all of the shower associations
  delete showerid_v; //vector of ids of all showers for each vertex
  delete shower_reco_totalE_v;//total energy of each reco shower
  delete shower_reco_length_v;//length of shower
  delete shower_reco_openingangle_v;//opening angle of shower
  delete shower_reco_ssnetshower_v; //shower fraction of hits
  delete shower_reco_dedx_vv ; //vector of dedx for each reco shower object





}
