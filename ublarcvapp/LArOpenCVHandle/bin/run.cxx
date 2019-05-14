#include "larcv/core/Processor/ProcessDriver.h"

#include "larcv/app/ImageMod/ROIClipper.h"
#include "ublarcvapp/LArOpenCVHandle/BlankImage.h"
#include "ublarcvapp/Filter/NuFilter.h"
#include "LArOpenCV/ImageCluster/AlgoClass/SuperClusterer.h"
#include "LArOpenCV/ImageCluster/AlgoModule/SuperClusterMaker.h"

int main(int argc, char** argv) {

  //std::string CFG     = "../cfg/for_pressnet.cfg";
  std::string CFG     = "../cfg/prod_fullchain_ssnet_combined_newtag_base_c10_union_server.cfg";
  std::string ANAFILE = "test_ana.root";
  std::string PGRFILE = "test_pgraph.root";
  std::string INFILE1 = "../../../../testdata/mcc9v12_intrinsicoverlay/supera-Run004955-SubRun000079.root";
  std::string INFILE2 = "../../../../testdata/mcc9v12_intrinsicoverlay/taggeroutv2-larcv-Run004955-SubRun000079.root";
  std::string INFILE3 = "../../../../testdata/mcc9v12_intrinsicoverlay/ssnetserveroutv2-larcv-Run004955-SubRun000079.root";
  std::string MCFILE  = "../../../../testdata/mcc9v12_intrinsicoverlay/larcvtruth-Run004955-SubRun000079.root";
  std::string OUTDIR  = ".";

  larcv::ROIClipper test;
  larcv::BlankImage test2;
  larocv::SuperClusterMaker test3;
  larcv::NuFilter   test4;
  
  larcv::ProcessDriver proc("ProcessDriver");
  std::cout << "Got configuration: " << CFG << std::endl;
  proc.configure(CFG);
  
  std::vector<std::string> input_v;
  input_v.push_back( INFILE1 );
  input_v.push_back( INFILE2 );
  input_v.push_back( INFILE3 );
  input_v.push_back( MCFILE  );  
  proc.override_input_file(input_v);

  proc.override_ana_file(    OUTDIR+"/"+ANAFILE );
  proc.override_output_file( OUTDIR+"/"+PGRFILE );

  proc.initialize();

  proc.batch_process(0,10);

  proc.finalize();

  return 0;
}
