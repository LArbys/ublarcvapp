#include "MergeDLInteraction.h"

#include <string>

#include "larcv/core/DataFormat/EventPGraph.h"
#include "larcv/core/DataFormat/EventPixel2D.h"

namespace ublarcvapp {
namespace ubdllee {

  static MergeDLInteractionFactory __global_MergeDLInteractionFactory__;  

  void MergeDLInteraction::configure( const larcv::PSet& pset) {

    _treename = pset.get<std::string>("TreeName");
    
  }

  void MergeDLInteraction::initialize() {

    char treename[254];
    sprintf( treename, "dlinteraction_%s_tree", _treename.c_str() );

    _event_interaction_v = new std::vector<DLInteraction>;
    _tree = new TTree( treename, "DL interaction tree: merges reco for each candidate vertex" );
    _tree->Branch( "event", _event_interaction_v );
    
  }

  bool MergeDLInteraction::process( larcv::IOManager& mgr ) {
    _event_interaction_v->clear();

    // get vertices
    larcv::EventPGraph*  event_vertex = (larcv::EventPGraph*)mgr.get_data(larcv::kProductPGraph, "test");
    larcv::EventPixel2D* event_ctrpix = (larcv::EventPixel2D*)mgr.get_data(larcv::kProductPixel2D, "test_ctor" );

    int nvertices = event_vertex->PGraphArray().size();
    LARCV_INFO() << "number of vertices=" << nvertices << std::endl;
    
    if ( nvertices==0 ) {
      _tree->Fill();
      return true;
    }

    _larliteio->syncEntry(mgr);
    
    // tracks
    larlite::event_track* track_v     = (larlite::event_track*)_larliteio->get_data(larlite::data::kTrack, "trackReco" );
    larlite::event_track* track_sce_v = (larlite::event_track*)_larliteio->get_data(larlite::data::kTrack, "trackReco_sceAdded" );

    larlite::event_shower* shower_v   = (larlite::event_shower*)_larliteio->get_data(larlite::data::kShower, "showerreco" );
    larlite::event_shower* shwr_dl_v  = (larlite::event_shower*)_larliteio->get_data(larlite::data::kShower, "dl" );    

    for (int ivtxidx=0; ivtxidx<nvertices; ivtxidx++ ) {

      DLInteraction dlvertex;
      dlvertex.run    = mgr.event_id().run();
      dlvertex.subrun = mgr.event_id().subrun();
      dlvertex.event  = mgr.event_id().event();    
      dlvertex.vertex_index = ivtxidx;

      dlvertex.vertex_pos.resize(3,0);
      dlvertex.vertex_wires.resize(3,0);      

      const larcv::PGraph& pgraph = event_vertex->PGraphArray().at(ivtxidx);
      int nparticles = pgraph.NumParticles();
      auto const& roi = pgraph.ParticleArray().at(0);
      dlvertex.vertex_pos[0] = roi.X();
      dlvertex.vertex_pos[1] = roi.Y();
      dlvertex.vertex_pos[2] = roi.Z();

      // contours
      for ( size_t ictr=0; ictr<pgraph.ClusterIndexArray().size(); ictr++ ) {
        int cluster_index = pgraph.ClusterIndexArray().at(ictr);
        std::vector< larcv::Pixel2DCluster > pixcluster_v;
        std::vector< larcv::ImageMeta >      meta_v;
        for (size_t p=0; p<3; p++ ) {
          pixcluster_v.push_back( event_ctrpix->Pixel2DClusterArray( p ).at(cluster_index) );
          meta_v.push_back( event_ctrpix->ClusterMetaArray( p ).at(cluster_index) );
        }
        dlvertex.vertex_prong_contour_v.emplace_back( std::move(pixcluster_v) );
        dlvertex.vertex_prong_meta_v.emplace_back( std::move(meta_v) );
        dlvertex.vertex_nprongs++;
      }

      // tracks
      for ( size_t itrk=0; itrk<track_v->size(); itrk++ ) {
        float dist=0;
        const larlite::track& track = track_v->at(itrk);
        for ( int v=0; v<3; v++ ) {
          float dx = dlvertex.vertex_pos[v]-track.LocationAtPoint(0)[v];
          dist += dx*dx;
        }
        //std::cout << "  track[" << itrk << "] dist2vertex=" << dist << std::endl;
        if ( dist<1.0e-6 ) {
          dlvertex.track_v.push_back( track );
          dlvertex.ntracks++;
        }
      }
      
      // showers
      for (size_t ishr=0; ishr<shower_v->size(); ishr++) {
        float dist=0;
        const larlite::shower& shower = shower_v->at(ishr);
        for ( int v=0; v<3; v++ ) {
          float dx = dlvertex.vertex_pos[v]-shower.ShowerStart()[v];
          dist += dx*dx;
        }
        //std::cout << "  shower[" << ishr << "] dist2vertex=" << dist << std::endl;
        if ( dist<1.0e-6 ) {
          dlvertex.shower_v.push_back( shower );
          dlvertex.nshowers++;
        }
      }

      // std::cout << "Vertex [" << dlvertex.vertex_index << "]" << std::endl;
      // std::cout << "  nprongs=" << dlvertex.vertex_nprongs << std::endl;
      // std::cout << "  ntracks=" << dlvertex.ntracks << std::endl;
      // std::cout << "  nshowers=" << dlvertex.nshowers << std::endl;
      
      _event_interaction_v->emplace_back( std::move(dlvertex) );
    }

    _tree->Fill();    
    return true;
  }

  void MergeDLInteraction::finalize() {
    _larliteio->close();
    _tree->Write();
  }


  void MergeDLInteraction::setup_larlite_io( const std::vector< std::string >& input_larlite_files ) {
    
    _larliteio = new ublarcvapp::LArliteManager( larlite::storage_manager::kREAD );
    for ( auto const& llfile : input_larlite_files )
      _larliteio->add_in_filename( llfile );
    _larliteio->open();

  }

}
}
