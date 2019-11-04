/**
 * \file SaveProbabilities
 *
 * \brief Class def header for a class SaveCutVariables
 * An overarching class that reads in all necessary root files for a run,subrun,event:
 *  Calculate all needed probabilities for likelihood
 *
 * @author katie
 */

#ifndef __SAVEPROBABILITIES_H__
#define __SAVEPROBABILITIES_H__

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

namespace ublarcvapp {
namespace ncpi0 {

  class SaveProbabilities{

  public:
    SaveProbabilities(){};
    virtual ~SaveProbabilities() {};

    //main code to run.
    void configure(const larcv::PSet& );
    void initialize();
    bool process( larcv::IOManager& io,larlite::storage_manager& ioll,larcv::IOManager& ioforward );
    void finalize();  // function to get true_vtx_location
    void setupAnaTree();

    // side functions
    std::vector<std::vector<double>> GetRecoVtxLocs(larlite::event_vertex* ev_vtxtracker);
    std::vector<std::vector<int>> TruthMatchTracks(std::vector<std::vector<double>> reco_vtx_v,
      larlite::event_track* ev_recotrack, std::vector<larcv::Image2D> instance_img,
      std::vector<larcv::Image2D> segment_img,larcv::ImageMeta wirey_meta);
    int MatchVtxToTrack(std::vector<std::vector<double>> reco_vtx_v,
      float xvtx,float yvtx,float zvtx);
    void DQDXProbs(std::vector<std::vector<int>> true_track_id_v,
        std::vector<std::vector<double>> reco_vtx_v, larlite::event_track*ev_recotrack);
    void TrackLengthProbs(std::vector<std::vector<int>> true_track_id_v,
      std::vector<std::vector<double>> reco_vtx_v, larlite::event_track*ev_recotrack);
    void SSNetShowerProbs(std::vector<std::vector<double>> reco_vtx_v,
      larlite::event_track*ev_recotrack,std::vector<larcv::Image2D> ssnet_img,std::vector<larcv::Image2D> wire_img,
      larcv::ImageMeta wire_meta, int plane,std::vector<std::vector<int>> trueid);
    void GetResRangeSlices();

  protected:

    // initialize all variables here...
    //output file
    TFile* OutFile;
    TFile* fin;

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

    //output variables
    int numslices = 20;
    // for R(Proton)
    TH2F* Muon_dqdx_resrange_h;
  	TH2F* Proton_dqdx_resrange_h;
    TH1D* Muon_resrange_slice_v[20];
    TH1D* Proton_resrange_slice_v[20];
    //for R(Gamma)
    TH1F* Muon_length_prob_h;
    TH1F* Gamma_length_prob_h;
    TH1F* Gamma_ssnetu_h;
    TH1F* Gamma_ssnetv_h;
    TH1F* Gamma_ssnety_h;
    TH1F* Muon_ssnetu_h;
    TH1F* Muon_ssnetv_h;
    TH1F* Muon_ssnety_h;
  };

}
}


#endif
