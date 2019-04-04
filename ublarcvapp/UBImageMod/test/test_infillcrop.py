import ROOT as rt
from larcv import larcv
from ublarcvapp import ublarcvapp

import sys
sys.argv.append("-b")

datafile = "larcv_test.root"

io = larcv.IOManager(larcv.IOManager.kBOTH,"IOManager")
io.add_in_file( datafile )
io.set_out_file( "cropinfill_data.root" )
io.initialize()

# -------------------------------------
# UBCropInfill

lfcrop_cfg="""Verbosity:0

InputADCProducer: \"wire\"
InputChStatusProducer: \"wire\"

OutputADCProducer:  \"ADC\"
OutputLabelsProducer:  \"Labels\"
OutputADCMaskedProducer:  \"ADCMasked\"
OutputWeightsProducer:  \"Weights\"
OutputCroppedMetaProducer: \"meta\"
OutputFilename: \"cropinfill_data.root\"

MaxImages: 1
LimitOverlap: false
MaxOverlapFraction: 0.25
"""

lfcfg = open("InfillDataCropper.cfg",'w')
print >>lfcfg,lfcrop_cfg
lfcfg.close()
lfpset = larcv.CreatePSetFromFile( "InfillDataCropper.cfg", "UBCropInfill" )

# -------------------------------------
# ALGOS
lfcrop_algo = ublarcvapp.InfillDataCropper()
lfcrop_algo.configure(lfpset)
lfcrop_algo.initialize()

# -------------------------------------

nentries = io.get_n_entries()
print "Num Entries: ",nentries

for n in range(0,nentries):
    print "on entry: " , n
    io.read_entry(n)
    lfcrop_algo.process( io );

lfcrop_algo.finalize()
io.finalize()


print "FIN"
