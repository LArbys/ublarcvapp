import ROOT as rt
from larcv import larcv
from ublarcvapp import ublarcvapp

import sys
sys.argv.append("-b")

datafile = "cropinfill_data.root"

io = larcv.IOManager(larcv.IOManager.kBOTH,"IOManager")
io.add_in_file( datafile )
io.set_out_file( "sparseinfill_data.root" )
io.initialize()

# -------------------------------------
# UBCropInfill

lfcrop_cfg="""Verbosity:0

InputADCProducer: \"ADC\"
InputLabelsProducer:  \"Labels\"
InputADCMaskedProducer:  \"ADCMasked\"

OutputADCProducer:  \"ADC\"
OutputADCMaskedProducer:  \"ADCMasked\"

OutputCroppedMetaProducer: \"meta\"
OutputFilename: \"sparseinfill_data.root\"
"""

lfcfg = open("InfillSparsifyImage.cfg",'w')
print >>lfcfg,lfcrop_cfg
lfcfg.close()
lfpset = larcv.CreatePSetFromFile( "InfillSparsifyImage.cfg", "InfillSparsifyImage" )

# -------------------------------------
# ALGOS
lfcrop_algo = ublarcvapp.InfillSparsifyImage()
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
