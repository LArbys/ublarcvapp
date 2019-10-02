import os,sys
import ROOT as rt
from ROOT import std
from larlite import larlite
from larcv import larcv
from ublarcvapp import ublarcvapp
print ublarcvapp.ubdllee.NuAnaMC

larlite_inputfiles = ["../../../../testdata/mcc9_v13_nueintrinsic_overlay_run1/reco2d-Run00499-SubRun000006.root"]
larcv_files = ["../../../../testdata/mcc9_v13_nueintrinsic_overlay_run1/supera-Run004999-SubRun000006.root"]
config = "driver_nuanamc.cfg"

out_larcv = "temp.root"

inputfiles_v = std.vector("std::string")()
for lcv in larcv_files:
    inputfiles_v.push_back(lcv)

driver = larcv.ProcessDriver("NuAnaMCProcessor")
driver.configure( config )
driver.override_input_file( inputfiles_v )
driver.override_output_file( out_larcv )

# get processor, add larlite files
processors = driver.process_map()
it_process = processors.find("DLTaggerProcess")
nuanamc = driver.process_ptr(it_process.second)

driver.initialize()

nentries = driver.io().get_n_entries()

for ientry in xrange(0,nentries):
    driver.process_entry(ientry)

driver.finalize()

