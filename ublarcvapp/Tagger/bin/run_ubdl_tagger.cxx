#include <iostream>

#include "larcv/core/Processor/ProcessDriver.h"
#include "larcv/core/DataFormat/IOManager.h"

#include "ublarcvapp/Tagger/TaggerProcessor.h"

#include "ublarcvapp/Tagger/RadialSegmentSearch.h"

int main( int nargs, char** argv ) {

  std::cout << "=======================================" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "  UB DL LEE 'classic' cosmic tagger " << std::endl;
  std::cout << "" << std::endl;
  std::cout << "=======================================" << std::endl;

  std::string supera = argv[1];
  std::string opreco = argv[2];

  larcv::IOManager io(larcv::IOManager::kBOTH,"tagger_larcv_input",larcv::IOManager::kTickBackward);
  io.set_out_file( "output_tagger_larcv.root" );
  io.add_in_file( supera );
  io.add_in_file( "saved_endpts.root" );
  io.initialize();
  
  ublarcvapp::tagger::TaggerProcessor tagger;
  tagger.configure( "example_cfgs/tagger_data_v2_splity.cfg" );
  tagger.add_larlite_input( opreco );
  tagger.set_larlite_output( "output_tagger_larlite.root" );
  tagger.initialize();

  int nentries = io.get_n_entries();
  std::cout << "Number of entries: " << nentries << std::endl;

  // for debugging one of the many subalgos
  //ublarcvapp::tagger::RadialSegmentSearch::_global_verbosity = 2;
  
  for (int ientry=0; ientry<nentries; ientry++ ) {
    std::cout << "ENTRY [" << ientry << " of " << nentries << "]" << std::endl;
    io.read_entry(ientry);
    tagger.process(io);
    io.save_entry();
    break;
  }

  std::cout << "----------" << std::endl;
  std::cout << "cleaning up" << std::endl;
  tagger.finalize();
  io.finalize();
  std::cout << "done" << std::endl;
  
  return 0;
}
