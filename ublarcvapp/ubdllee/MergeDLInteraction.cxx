#include "MergeDLInteraction.h"

#include <string>

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

    
    
    _tree->Fill();
    return true;
  }

  void MergeDLInteraction::finalize() {
    _tree->Write();
  }
  
}
}
