#include <iostream>

#include "larcv/core/Processor/ProcessDriver.h"
#include "larcv/core/DataFormat/IOManager.h"

#include "ublarcvapp/Tagger/TaggerProcessor.h"

int main( int nargs, char** argv ) {

  std::cout << "=======================================" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "  UB DL LEE 'classic' cosmic tagger " << std::endl;
  std::cout << "" << std::endl;
  std::cout << "=======================================" << std::endl;

  std::string supera = argv[1];
  std::string opreco = argv[2];

  larcv::IOManager io(larcv::IOManager::kBOTH,"tagger_larcv_input");
  io.set_out_file( "output_tagger.root" );
  io.add_in_file( supera );
  io.initialize();
  
  ublarcvapp::tagger::TaggerProcessor tagger;
  tagger.configure( "example_cfgs/tagger_data_v2_splity.cfg" );
  tagger.add_larlite_input( opreco );
  tagger.initialize();

  int nentries = io.get_n_entries();
  std::cout << "Number of entries: " << nentries << std::endl;

  for (int ientry=0; ientry<nentries; ientry++ ) {
    std::cout << "ENTRY [" << ientry << " of " << nentries << "]" << std::endl;
    io.read_entry(ientry);
    tagger.process(io);
  }

  std::cout << "----------" << std::endl;
  std::cout << "cleaning up" << std::endl;
  tagger.finalize();
  io.finalize();
  std::cout << "done" << std::endl;
  
  return 0;
}
