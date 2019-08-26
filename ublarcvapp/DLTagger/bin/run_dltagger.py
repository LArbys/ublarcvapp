import os,sys

import ROOT as rt
from ROOT import std
from larcv import larcv
from ublarcvapp import ublarcvapp
factory = ublarcvapp.dltagger.DLTaggerProcessFactory() # load the library, which registers the factory
print factory

inputfiles = std.vector("std::string")()
inputfiles.push_back( "../test/testset1/out_larcv_test.root" )

driver = larcv.ProcessDriver("DLTagger")
driver.configure( "dltagger.cfg" )
driver.override_input_file( inputfiles )
driver.override_output_file( "output_dltagger_test.root" )
driver.initialize()

nentries = driver.io().get_n_entries()

for ientry in xrange(nentries):
    driver.process_entry(ientry)
#driver.process_entry(2);

driver.finalize()
