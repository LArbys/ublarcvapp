/**
 * \file SaveCutVariables
 *
 * \brief Class def header for a class SaveCutVariables
 * An overarching class that reads in all necessary root files for a run,subrun,event:
 *  Calculate all needed cut variables
 *  Save to an output root file that can be hadded
 *
 * @author katie
 */

/** \addtogroup core_DataFormat

    @{*/
#ifndef __SAVECUTVARIABLES_H__
#define __SAVECUTVARIABLES_H__

// includes
#include <iostream>
#include <map>
#include <utility>
#include <string>
#include <cstring>
#include <bits/stdc++.h>
// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TVector3.h"
// larlite
#include "DataFormat/storage_manager.h"
#include "DataFormat/larflow3dhit.h"
#include "DataFormat/larflowcluster.h"
#include "DataFormat/opflash.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcnu.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/event_ass.h"
// larutil
#include "LArUtil/LArProperties.h"
#include "LArUtil/DetectorProperties.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/ClockConstants.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"

// larcv
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"
#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"
#include "larcv/core/DataFormat/EventClusterMask.h"
#include "larcv/core/DataFormat/EventPGraph.h"
#include "larcv/core/DataFormat/EventROI.h"

#include <cstdlib>
#include <math.h>

#include "Utils.h"
#include "SaveProbabilities.h"

namespace ublarcvapp {
namespace ncpi0 {

  class SaveCutVariables{

  public:
    SaveCutVariables(){};
    virtual ~SaveCutVariables() {};

    //main code to run.
    void configure(const larcv::PSet& );
    void initialize();
    bool process(  larcv::IOManager& io,larlite::storage_manager& ioll,larcv::IOManager& ioforward  );
    void finalize();  // function to get true_vtx_location
    void setupAnaTree();
    void ClearBranches();
    void LoadProbHists();

    // side functions
    std::vector<double> GetTrueVtxLoc(larcv::EventROI* ev_partroi, TH3F* sceDx,
          TH3F* sceDy, TH3F* sceDz);
    std::vector<std::vector<double>> GetRecoVtxLocs(larlite::event_vertex* ev_vtxtracker);
    std::vector<bool> IsVtxGood(std::vector<std::vector<double>> reco_vtx_v,
          std::vector<double> true_vtx);
    std::vector<bool> IsVtxInFid(std::vector<std::vector<double>> reco_vtx_v);
    void TrackLengthAndContainment(std::vector<std::vector<double>> reco_vtx_v,
          larlite::event_track*ev_recotrack);
    void AvgMaxDQDX( std::vector<std::vector<double>> reco_vtx_v,
          larlite::event_track*ev_recotrack);
    void SaveSSNetFrac(std::vector<std::vector<double>> reco_vtx_v,
      larlite::event_track*ev_recotrack,std::vector<larcv::Image2D> ssnet_img,
      larcv::ImageMeta wire_meta, int plane);
    void CalculateR_Proton(std::vector<std::vector<double>> reco_vtx_v,
      larlite::event_track*ev_recotrack);
    void CalculateR_Gamma(std::vector<std::vector<double>> reco_vtx_v,
      larlite::event_track*ev_recotrack);

  protected:

    // initialize all variables here...
    //output file
    TFile* fin;
    TFile* OutFile;
    TTree* _ana_tree;

    int run;
    int subrun;
    int event;
    //input config variables
    std::string _input_adc_producer;
    std::string _input_partroi_producer;
    std::string _input_segment_producer;
    std::string _input_instance_producer;
    //ssnet inputs
    std::string _input_ssnet_uplane_producer;
    std::string _input_ssnet_vplane_producer;
    std::string _input_ssnet_yplane_producer;
    // larlite (mcinfo) inputs
    std::string _input_mctrack_producer;
    std::string _input_mcshower_producer;
    std::string _input_flux_producer;
    std::string _input_mctruth_producer;
    //larlite (tracker) inputs
    std::string _input_recotrack_producer;
    std::string _input_vtxtracker_producer;

    //input Histograms
    std::vector<float> Rweight_v = std::vector<float> (20,1.0);
    TH1D* Muon_resrange_slice_prob_v[20];
    TH1D* Proton_resrange_slice_prob_v[20];
    TH1F* Muon_length_prob_h;
    TH1F* Gamma_length_prob_h;
    TH1F* Gamma_ssnetu_h;
    TH1F* Gamma_ssnetv_h;
    TH1F* Gamma_ssnety_h;
    TH1F* Muon_ssnetu_h;
    TH1F* Muon_ssnetv_h;
    TH1F* Muon_ssnety_h;

    //output variables
    std::vector<std::vector<int>> true_track_id_v; //vector storing track id for reco tracks for each vtx
    std::vector<std::vector<float>> R_proton_v; //vector of R(Proton) for recotracks for each vtx
    std::vector<std::vector<float>> R_gamma_v; //vector of R(Proton) for recotracks for each vtx
    std::vector<std::vector<float>> SSNet_shower_frac_u_v;//vector of SSNet(uplane) shower fraction for recotracks for each vtx
    std::vector<std::vector<float>> SSNet_track_frac_u_v;//vector of SSNet(uplane) shower fraction for recotracks for each vtx
    std::vector<std::vector<float>> SSNet_shower_frac_v_v;//vector of SSNet(vplane) shower fraction for recotracks for each vtx
    std::vector<std::vector<float>> SSNet_track_frac_v_v;//vector of SSNet(vplane) shower fraction for recotracks for each vtx
    std::vector<std::vector<float>> SSNet_shower_frac_y_v;//vector of SSNet(yplane) shower fraction for recotracks for each vtx
    std::vector<std::vector<float>> SSNet_track_frac_y_v;//vector of SSNet(yplane) shower fraction for recotracks for each vtx
    std::vector<std::vector<float>> tracklength_v;//vector of track length of recotracks for eachvtx;
    std::vector<std::vector<float>> max_dqdx_v;// vector of max dqdx of recotracks for each vtx;
    std::vector<std::vector<float>> avg_dqdx_v;// vector of avg dqdx of recotracks for each vtx;
    std::vector<std::vector<bool>> vtx_cont_v;//vector of boolians of is track contained?
    bool vtx_true_fid;//vector of boolians on whether vtx is in fid
    std::vector<bool> vtx_reco_fid_v;//vector of boolians on whether vtx is in fid
    std::vector<bool> vtx_status_v;// 0 if bad vertex, 1 if good
    // for location objects:
    //  2D: [time, uwire, vwire, ywire]
    //  3D: [X,Y,Z]
    std::vector<double> true_vtx_location_3D; //location of true vtx in event;
    std::vector<int> true_vtx_location_2D; //location of true vtx in event;
    std::vector<std::vector<double>> reco_vtx_location_3D_v; //vector of locations of reco vtx
    std::vector<std::vector<int>> reco_vtx_location_2D_v; //vector of locations of reco vtx

  };

}
}


#endif
