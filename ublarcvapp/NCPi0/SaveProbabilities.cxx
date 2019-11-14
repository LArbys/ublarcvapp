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
    //larlite (shower) inputs
    _input_recoshower_producer   = pset.get<std::string>("InputRecoShowerProducer");
    _input_pfpartshower_producer = pset.get<std::string>("InputPfPartShowerProducer");
    _input_hitsshower_producer   = pset.get<std::string>("InputHitsShowerProducer");
    _input_clustershower_producer= pset.get<std::string>("InputClusterShowerProducer");
    _input_assshower_producer    = pset.get<std::string>("InputAssShowerProducer");
    _input_assdlshower_producer  = pset.get<std::string>("InputAssDLShowerProducer");
    _input_vtxshower_producer  = pset.get<std::string>("InputVtxShowerProducer");

  }

  void SaveProbabilities::initialize(std::string showerrecoananame) {
    ShowerRecoFile = new TFile(showerrecoananame.c_str(),"read");
    // GetShowerRecoVals();
    //hard code in CalibrationFile
    CalibrationFile = new TFile("/cluster/tufts/wongjiradlab/kmason03/testdir/ubdl/ublarcvapp/ublarcvapp/NCPi0/test/CalibrationMaps_MCC9.root","read");
    hImageCalibrationMap_00 = (TH3D*)CalibrationFile->Get("hImageCalibrationMap_00");
    hImageCalibrationMap_01 = (TH3D*)CalibrationFile->Get("hImageCalibrationMap_01");
    hImageCalibrationMap_02 = (TH3D*)CalibrationFile->Get("hImageCalibrationMap_02");
    OutFile = new TFile("NCPi0Probablities.root","RECREATE");
    setupAnaTree();
  }

bool SaveProbabilities::process(larcv::IOManager& io,larlite::storage_manager& ioll,larcv::IOManager& ioforward,int ientry){
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
  const auto ev_recoshower     = (larlite::event_shower*)ioll.get_data(larlite::data::kShower,  _input_recoshower_producer );
  const auto ev_pfpartshower   = (larlite::event_pfpart*)ioll.get_data(larlite::data::kPFParticle,  _input_pfpartshower_producer );
  const auto ev_hitsshower     = (larlite::event_hit*)ioll.get_data(larlite::data::kHit, _input_hitsshower_producer);
  const auto ev_assshower      = (larlite::event_ass*)ioll.get_data(larlite::data::kAssociation, _input_assshower_producer);
  const auto ev_assdlshower      = (larlite::event_ass*)ioll.get_data(larlite::data::kAssociation, _input_assdlshower_producer);
  const auto ev_clustershower  = (larlite::event_cluster*)ioll.get_data(larlite::data::kCluster, _input_clustershower_producer);
  const auto ev_vtxshower     = (larlite::event_vertex*)ioll.get_data(larlite::data::kVertex,  _input_vtxtracker_producer );

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

  //shower variables
  std::cout<<"size of event_shower: "<<ev_recoshower->size()<<std::endl;
  std::cout<<"size of vtx from shower: "<<ev_vtxshower->size()<<std::endl;
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
  //get hit association vectors
  std::vector<std::vector<int>> hitsid = HitsVtxAssociation(ev_hitsshower,ev_clustershower,ev_assdlshower);
  //also get the opposite (vtx to hits)
  std::vector<std::vector<int>> vtxhits_v = VtxHitsAssociation( hitsid, ev_vtxshower->size());
  std::vector<std::vector<int>> showerid = ShowerVtxAssociation(ev_hitsshower,
    ev_recoshower,hitsid,wireu_meta,wirev_meta,wirey_meta);

  // make one large association object: vertex[shower][hit]
  std::vector<std::vector<std::vector<int>>> ShowerAssociation_vvv = TotalShowerAssociation(ev_recoshower,
    ev_hitsshower,vtxhits_v,showerid,ev_vtxshower->size(),wireu_meta,wirev_meta,wirey_meta);

  std::vector<std::vector<int>> showerid_v = TruthMatchShowers(ShowerAssociation_vvv,
        ev_hitsshower,instance_img, segment_img, wireu_meta,wirev_meta,wirey_meta);

  std::cout<<"Number of Vertices: "<< ShowerAssociation_vvv.size()<<std::endl;
  for (int ii =0;ii< ShowerAssociation_vvv.size();ii++){
    std::cout<< "-- Vertex: "<<ii<<" has "<< ShowerAssociation_vvv.at(ii).size()<<" associated showers"<<std::endl;
    for (int iii = 0;iii< ShowerAssociation_vvv.at(ii).size();iii++){
      std::cout<<"----Shower: "<<iii<< " has "<< ShowerAssociation_vvv.at(ii).at(iii).size()<<" associated hits"<<std::endl;
      std::cout<<"--------It's shower number is: "<< ShowerAssociation_vvv.at(ii).at(iii).at(0)<<std::endl;
      std::cout<<"--------It's truth ID is: " << showerid_v.at(ii).at(iii)<<std::endl;
      if(ShowerAssociation_vvv.at(ii).at(iii).size()>1){
        showerid_h->Fill(showerid_v.at(ii).at(iii));
      }
    }
  }
  return true;
}
//------------------------------------------------------------------------------
void SaveProbabilities::GetShowerRecoVals(){
  ShowerRecoFile->cd();
  showerrecotree = (TTree*)ShowerRecoFile->Get("ShowerQuality_DL");
  showerrecotree->SetBranchAddress("nshowers",&nshowers);
}//end of function
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
//------------------------------------------------------------------------------
std::vector<std::vector<int>> SaveProbabilities::HitsVtxAssociation(larlite::event_hit* ev_hitsshower,
    larlite::event_cluster* ev_clustershower,larlite::event_ass* ev_assdlshower){
  //function to use the association tree to match hits with vtx
  //output is a vtx with size  = number of hits, each entry is list of associated vtx
  std::vector<std::vector<int>> hitid_v (ev_hitsshower->size());
  //middle step vectors (hits->cluster->pfpart->vertex)
  std::vector<std::vector<int>> hits_to_cluster (ev_hitsshower->size());
  std::vector<std::vector<int>> cluster_to_pfpart (ev_hitsshower->size());
  //loop through hit to get association
  if (ev_assdlshower->size()>0){
    //get partid
    larlite::product_id hitid = larlite::product_id(larlite::data::kHit,"dl");
    larlite::product_id clusterid = larlite::product_id(larlite::data::kCluster,"dl");
    larlite::product_id pfpartid = larlite::product_id(larlite::data::kPFParticle,"dl");
    larlite::product_id vertexid = larlite::product_id(larlite::data::kVertex, "dl");

    //get index of associations
    auto const& index0 = ev_assdlshower->assid(hitid,clusterid);
    auto const& index1 = ev_assdlshower->assid(clusterid,pfpartid);
    auto const& index2 = ev_assdlshower->assid(pfpartid,vertexid);
    //get association data
    std::vector<std::vector<unsigned int> > hitclusterdata = ev_assdlshower->association(index0);
    std::vector<std::vector<unsigned int> > clusterpfpartdata = ev_assdlshower->association(index1);
    std::vector<std::vector<unsigned int> > pfpartvertexdata = ev_assdlshower->association(index2);

    //loop through hits to match to clusters
    for (int hit = 0; hit<ev_hitsshower->size(); hit++){
      auto& cluster_idx = hitclusterdata.at(hit);
      if ( cluster_idx.size()>0 ) {
        for ( auto const& cidx : cluster_idx ){
           hits_to_cluster.at(hit).push_back(cidx);
           // std::cout<< " cidx: "<< cidx <<std::endl;
         }
      }
    }

    //loop through clusters to match to pfpart
    for(int cluster = 0; cluster<hits_to_cluster.size(); cluster++){
      // loop through all associated clusters to find pfpart
      for (int ii = 0; ii<hits_to_cluster.at(cluster).size();ii++){
        auto& pfpart_idx = clusterpfpartdata.at(hits_to_cluster.at(cluster).at(ii));
        cluster_to_pfpart.at(cluster).push_back(pfpart_idx.at(0));
      }
    }

    //loop through pfpart to match to vertex
    for(int hit =0;hit<cluster_to_pfpart.size();hit++){
      for(int cluster = 0; cluster<cluster_to_pfpart.at(hit).size();cluster++ ){
        int vtx_idx =  pfpartvertexdata.at(cluster_to_pfpart.at(hit).at(cluster)).at(0);
        hitid_v.at(hit).push_back(vtx_idx);
      }
    }
  }

  return hitid_v;
}//end of function
//-----------------------------------------------------------------------------
std::vector<std::vector<int>> SaveProbabilities::VtxHitsAssociation(std::vector<std::vector<int>> hitsid, int nvertices){
  //function to switch from hits->vtx to vtx->hits
  std::vector<std::vector<int>> vtx_to_hitsid_v (nvertices);
  //loop through hitsid
  for(int hit = 0; hit<hitsid.size(); hit++){
    if(hitsid.at(hit).size()<=nvertices){
      for (int vtx = 0; vtx<hitsid.at(hit).size(); vtx++){
        vtx_to_hitsid_v.at(vtx).push_back(hit);
      }
    }
  }
  return vtx_to_hitsid_v;
}//end of function
//------------------------------------------------------------------------------
std::vector<std::vector<int>> SaveProbabilities::ShowerVtxAssociation(larlite::event_hit* ev_hitsshower,
      larlite::event_shower* ev_recoshower,std::vector<std::vector<int>>hitsid,
      larcv::ImageMeta wireu_meta,larcv::ImageMeta wirev_meta,larcv::ImageMeta wirey_meta){
   //function to match reco showers to vtx (sort of  hacky :( )
   Utils Utils;
   std::vector<std::vector<int>> showerid (ev_recoshower->size());
   //loop through showers
   for(int shower = 0; shower<ev_recoshower->size(); shower++){
     // get shower start point
     TVector3 shower_start = ev_recoshower->at(shower).ShowerStart();
     std::vector<double> startpt_v;
     startpt_v.push_back(shower_start.X());
     startpt_v.push_back(shower_start.Y());
     startpt_v.push_back(shower_start.Z());
     std::vector<int> showerpt_pix = Utils.getProjectedPixel(startpt_v,wirey_meta,3);
     //which hit is closest to start point - loop through hits
     //(only using y plane for simplicity)
     float mindistance = 1000000;
     int closesthit = -1;
     for (int hit =0; hit<ev_hitsshower->size();hit++){
       //get distance to hit
       if(ev_hitsshower->at(hit).View()==2){
         int row = wirey_meta.row(ev_hitsshower->at(hit).PeakTime()+2400);
         int col = wirey_meta.col(ev_hitsshower->at(hit).Channel()-2*2400);
         float distance = std::sqrt(std::pow(float(row)-showerpt_pix[0],2)+std::pow(float(col)-showerpt_pix[3],2));
         if(distance<mindistance){
           mindistance = distance;
           closesthit = hit;
         }
       }
       if (ev_hitsshower->at(hit).View()==0){
         int row = wireu_meta.row(ev_hitsshower->at(hit).PeakTime()+2400);
         int col = wireu_meta.col(ev_hitsshower->at(hit).Channel());
         float distance = std::sqrt(std::pow(float(row)-showerpt_pix[0],2)+std::pow(float(col)-showerpt_pix[1],2));
         if(distance<mindistance){
           mindistance = distance;
           closesthit = hit;
         }
       }
       if (ev_hitsshower->at(hit).View()==1){
         int row = wirev_meta.row(ev_hitsshower->at(hit).PeakTime()+2400);
         int col = wirev_meta.col(ev_hitsshower->at(hit).Channel()-2400);
         float distance = std::sqrt(std::pow(float(row)-showerpt_pix[0],2)+std::pow(float(col)-showerpt_pix[2],2));
         if(distance<mindistance){
           mindistance = distance;
           closesthit = hit;
         }
       }

     }//end of loop though hits
     // vtxid of hit?
     std::vector<int> vtxlist = hitsid.at(closesthit);
     //take all ids -- overly cautious
     showerid.at(shower)=vtxlist;
   }
   return showerid;
 }//end of function
//------------------------------------------------------------------------------
std::vector<std::vector<std::vector<int>>> SaveProbabilities::TotalShowerAssociation(
    larlite::event_shower* ev_recoshower,larlite::event_hit* ev_hitsshower,
    std::vector<std::vector<int>> vtxhits_v,std::vector<std::vector<int>> showerid,
    int nvertices,larcv::ImageMeta wireu_meta,larcv::ImageMeta wirev_meta,
    larcv::ImageMeta wirey_meta){
  //function to create one giant association object...vertex[shower[hit]]
  Utils Utils;
  std::vector<std::vector<std::vector<int>>> total_ass_vvv (nvertices);
  //loop through Vertices
  for(int vtx = 0; vtx<nvertices; vtx++){
    //loop through shower id to see how many match
    int nshowers = 0;
    for(int shower = 0; shower<showerid.size(); shower++){
      bool showermatch = false;
      for(int ii =0;ii<showerid.at(shower).size();ii++){
        if( showerid.at(shower).at(ii) == vtx) showermatch = true;
      }//end of loop trhough vertices associated with shower
      if (showermatch)nshowers++;
    }//end of loop through shower id
    //loop again to save ids and locations
    std::vector<int> matchshowerid_v;
    //vectors to save positions of showers
    std::vector<float> wireu_loc;
    std::vector<float> wirev_loc;
    std::vector<float> wirey_loc;
    std::vector<float> time_loc;
    for(int shower = 0; shower<showerid.size(); shower++){
      bool showermatch = false;
      for(int ii =0;ii<showerid.at(shower).size();ii++){
        if( showerid.at(shower).at(ii) == vtx) showermatch = true;
      }//end of loop trhough vertices associated with shower
      //if there is a map, save locations...
      if (showermatch){
        matchshowerid_v.push_back(shower);
        TVector3 shower_start = ev_recoshower->at(shower).ShowerStart();
        std::vector<double> startpt_v;
        startpt_v.push_back(shower_start.X());
        startpt_v.push_back(shower_start.Y());
        startpt_v.push_back(shower_start.Z());
        std::vector<int> showerpt_pix = Utils.getProjectedPixel(startpt_v,wirey_meta,3);
        time_loc.push_back(showerpt_pix[0]);
        wireu_loc.push_back(showerpt_pix[1]);
        wirev_loc.push_back(showerpt_pix[2]);
        wirey_loc.push_back(showerpt_pix[3]);
      }
    }//end of loop through shower id
    std::vector<std::vector<int>> shower_vv (nshowers);

    //loop through hits associated with vtx, figure out closest shower
    std::vector<int> vtx_shower_v;
    for(int ii = 0;ii<vtxhits_v.at(vtx).size();ii++){
      //get hit location
      int row;
      int col;
      int plane;
      if(ev_hitsshower->at(ii).View()==2){
        plane = 2;
        row = wirey_meta.row(ev_hitsshower->at(ii).PeakTime()+2400);
        col = wirey_meta.col(ev_hitsshower->at(ii).Channel()-2*2400);
      }
      else if (ev_hitsshower->at(ii).View()==0){
        row = wireu_meta.row(ev_hitsshower->at(ii).PeakTime()+2400);
        col = wireu_meta.col(ev_hitsshower->at(ii).Channel());
        plane = 0;
      }
      if (ev_hitsshower->at(ii).View()==1){
        row = wirev_meta.row(ev_hitsshower->at(ii).PeakTime()+2400);
        col = wirev_meta.col(ev_hitsshower->at(ii).Channel()-2400);
        plane = 1;
      }
      int closestshower = -1;
      int mindistance = 100000;

      for(int shower = 0; shower<nshowers;shower++){
        float distance;
        if (plane ==0) distance = std::sqrt(std::pow(float(row)-time_loc[shower],2)+std::pow(float(col)-wireu_loc[shower],2));
        else if (plane == 1) distance = std::sqrt(std::pow(float(row)-time_loc[shower],2)+std::pow(float(col)-wirev_loc[shower],2));
        else distance = std::sqrt(std::pow(float(row)-time_loc[shower],2)+std::pow(float(col)-wirey_loc[shower],2));
        if(distance<mindistance){
          mindistance = distance;
          closestshower  = shower;
        }
      }//end of shower loop
    vtx_shower_v.push_back(closestshower);
    }//end of loop through hits

    for (int shower = 0;shower<matchshowerid_v.size();shower++){
      //FIRST INDEX IS SHOWER ID!!
      std::vector<int> hits_v;
      hits_v.push_back(matchshowerid_v.at(shower));
      for(int ii = 0;ii<vtxhits_v.at(vtx).size();ii++){
        if (vtx_shower_v.at(ii) == shower){
          hits_v.push_back(vtxhits_v.at(vtx).at(ii));
        }
      }

      shower_vv.at(shower) = hits_v;
    }

    total_ass_vvv.at(vtx)= shower_vv;
  }//end of loop through vertices - first pass

  //now going through again to add hits

  return total_ass_vvv;
}//end of funciton
//------------------------------------------------------------------------------
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
std::vector<std::vector<int>> SaveProbabilities::TruthMatchShowers(
      std::vector<std::vector<std::vector<int>>> ShowerAssociation_vvv,
      larlite::event_hit* ev_hitsshower,std::vector<larcv::Image2D> instance_img,
      std::vector<larcv::Image2D> segment_img,larcv::ImageMeta wireu_meta,larcv::ImageMeta wirev_meta,
      larcv::ImageMeta wirey_meta){
  //function to truth match showers
  //output is a vector of truthids for each vertex
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
  std::vector<std::vector<int>> showerid_v;
  for (int vtx = 0; vtx<ShowerAssociation_vvv.size(); vtx++){
    std::vector<int> showerid;
    int NumMuonCosmic   = 0;
    int NumMuonNeutrino = 0;
    int NumElectron     = 0;
    int NumGamma        = 0;
    int NumPi0          = 0;
    int NumPiCharged    = 0;
    int NumKaon         = 0;
    int NumProton       = 0;
    //loop through associated showers
    for(int ii = 0; ii<ShowerAssociation_vvv.at(vtx).size(); ii++){

      // to do truth match to particle type
      // -- loop through hits, get instance image of each pixel
      // -- type is whatever has highest number of instance image
      int NumMuonCosmicPix   = 0;
      int NumMuonNeutrinoPix = 0;
      int NumElectronPix     = 0;
      int NumGammaPix        = 0;
      int NumPi0Pix          = 0;
      int NumPiChargedPix    = 0;
      int NumKaonPix         = 0;
      int NumProtonPix       = 0;
      int numshowerpts = 0;
      //loop through hits
      for (int iii = 1; iii<ShowerAssociation_vvv.at(vtx).at(ii).size(); iii++){
        //get location:
        int hitindex = ShowerAssociation_vvv.at(vtx).at(ii).at(iii);
        int plane;
        int row;
        int col;
        if(ev_hitsshower->at(hitindex).View()==2){
          plane = 2;
          row = wirey_meta.row(ev_hitsshower->at(hitindex).PeakTime()+2400);
          col = wirey_meta.col(ev_hitsshower->at(hitindex).Channel()-2*2400);
        }
        else if (ev_hitsshower->at(hitindex).View()==0){
          row = wireu_meta.row(ev_hitsshower->at(hitindex).PeakTime()+2400);
          col = wireu_meta.col(ev_hitsshower->at(hitindex).Channel());
          plane = 0;
        }
        else if (ev_hitsshower->at(hitindex).View()==1){
          row = wirev_meta.row(ev_hitsshower->at(hitindex).PeakTime()+2400);
          col = wirev_meta.col(ev_hitsshower->at(hitindex).Channel()-2400);
          plane = 1;
        }
        //get image values
        float instance_val;
        float segment_val;
        if(plane ==0){
          instance_val = instance_img[0].pixel(row,col);
          segment_val = segment_img[0].pixel(row,col);
        }
        else if(plane ==1){
          instance_val = instance_img[1].pixel(row,col);
          segment_val = segment_img[1].pixel(row,col);
        }
        else if(plane ==2){
          instance_val = instance_img[2].pixel(row,col);
          segment_val = segment_img[2].pixel(row,col);
        }
        if (instance_val>0){
          if (segment_val==3) NumElectronPix++;
          else if (segment_val==4) NumGammaPix++;
          else if (segment_val==5) NumPi0Pix++;
          else if (segment_val==6) NumMuonNeutrinoPix++;
          else if (segment_val==7) NumKaonPix++;
          else if (segment_val==8) NumPiChargedPix++;
          else if (segment_val==9) NumProtonPix++;
        }
        else NumMuonCosmicPix++;
        numshowerpts++;

      }//end of loop through hits
      if (numshowerpts == 0)numshowerpts = 1;
      float FracProton       = (float)NumProtonPix/((float)numshowerpts);
      float FracMuonCosmic   = (float)NumMuonCosmicPix/((float)numshowerpts);
      float FracMuonNeutrino = (float)NumMuonNeutrinoPix/((float)numshowerpts);
      float FracGamma        = (float)NumGammaPix/((float)numshowerpts);
      float FracPi0          = (float)NumPi0Pix/((float)numshowerpts);
      float FracPiCharged    = (float)NumPiChargedPix/((float)numshowerpts);
      float FracKaon         = (float)NumKaonPix/((float)numshowerpts);
      float FracElectron     = (float)NumElectronPix/((float)numshowerpts);

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

      // if (id ==0 && numshowerpts>1){
      //   std::cout<<"Unknown Particle: "<<ShowerAssociation_vvv.at(vtx).at(ii).at(0)<<std::endl;
      //   std::cout<<"FracProton: "<<FracProton<<std::endl;
      //   std::cout<<"FracGamma: "<<FracGamma<<std::endl;
      //   std::cout<<"FracCosmic: "<<FracMuonCosmic<<std::endl;
      //   std::cout<<"FracMuon: "<<FracMuonNeutrino<<std::endl;
      //   std::cout<<"FracElectron: "<<FracElectron<<std::endl;
      //   std::cout<<"FracKaon: "<<FracKaon<<std::endl;
      //   std::cout<<"FracPi0: "<<FracPi0<<std::endl;
      //   std::cout<<"FracPiCharged: "<<FracPiCharged<<std::endl;
      // }
      showerid.push_back(id);
    }//end of shower loop
    showerid_v.push_back(showerid);
  }//end of vtx loop
  return showerid_v;
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
          TVector3 trackpt = ev_recotrack->at(ii).LocationAtPoint(iii);
          //get calibration factors..
          double bin = hImageCalibrationMap_00->FindBin(trackpt.X(),trackpt.Y(),trackpt.Z());
          double calib_u = hImageCalibrationMap_00->GetBinContent(bin);
          double calib_v = hImageCalibrationMap_01->GetBinContent(bin);
          double calib_y = hImageCalibrationMap_02->GetBinContent(bin);
          //use y plane, if dead, average the others.
          if (dqdx_val == 0){
            float uval = ev_recotrack->at(ii).DQdxAtPoint(iii,larlite::geo::View_t::kU)*calib_u;
            float vval = ev_recotrack->at(ii).DQdxAtPoint(iii,larlite::geo::View_t::kV)*calib_v;
            dqdx_val = (uval+vval)/2.0;
          }
          else dqdx_val = dqdx_val*calib_y;

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
  // OutFile->cd();
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
  showerid_h = new TH1F("ShowerID","ShowerID",10,0,10);


}


//------------------------------------------------------------------------------

void SaveProbabilities::finalize() {

  OutFile->cd();
  // _ana_tree->Write();
  OutFile->Write();

}
}//end of ncpi0 namespace
}//end of ublarcvapp namespace
