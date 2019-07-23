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


  larcv::IOManager io(larcv::IOManager::kBOTH,"tagger");
  io.set_out_file( "output_tagger.root" );
  io.add_in_file( argv[1] );
  
  ublarcvapp::tagger::TaggerProcessor tagger;
  
  
  return 0;
}
