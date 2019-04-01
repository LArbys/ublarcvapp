import sys
import ROOT as rt
from larcv import larcv
larcv.load_rootutil()
from ublarcvapp import ublarcvapp
from ROOT import std

rt.gStyle.SetOptStat(0)

superafile = sys.argv[1]
input_adc_producer = sys.argv[2]

io = larcv.IOManager(larcv.IOManager.kREAD,"",larcv.IOManager.kTickBackward)
io.add_in_file( superafile )
io.initialize()

out = larcv.IOManager(larcv.IOManager.kWRITE)
out.set_out_file( "baka.root" )
out.initialize()

# -------------------------------------
# UBSplitDetector

scfg="""Verbosity: 2
InputProducer: \"%s\"
OutputBBox2DProducer: \"detsplit\"
CropInModule: true
OutputCroppedProducer: \"detsplit\"
BBoxPixelHeight: 512
BBoxPixelWidth: 832
CoveredZWidth: 310
FillCroppedYImageCompletely: true
DebugImage: false
MaxImages: -1
RandomizeCrops: false
MaxRandomAttempts: 50
MinFracPixelsInCrop: -0.0001
"""

fcfg = open("ubsplit.cfg",'w')
print >>fcfg,scfg%(input_adc_producer)
fcfg.close()


cfg = larcv.CreatePSetFromFile( "ubsplit.cfg", "UBSplitDetector" )
algo = ublarcvapp.UBSplitDetector()
algo.initialize()
algo.configure(cfg)
algo.set_verbosity(0)

# -------------------------------------

cfgcrop="""
"""


nentries = io.get_n_entries()
print "Num Entries: ",nentries
nentries = 1

for n in range(nentries):
    io.read_entry(n)

    ev_adc = io.get_data(larcv.kProductImage2D,"wiremc")
    adc_v  = ev_adc.Image2DArray()

    roi_v = std.vector("larcv::ROI")()
    out_v = std.vector("larcv::Image2D")()
    algo.process( adc_v, out_v, roi_v  )

    detsplit = out.get_data( larcv.kProductImage2D, "detsplit" )
    print "num rois: ",roi_v.size()
    print "num crops: ",out_v.size()

    if True:
        # visualize
        h_v = {}
        for i in xrange(out_v.size()):
            h_v[i] = larcv.as_th2d( out_v.at(i), "test%d"%(i) )
            h_v[i].GetZaxis().SetRangeUser(-10,100)
        #print h_v

        c = rt.TCanvas("c","c",1500,400)
        c.Divide(3,1)
        for i in range(len(h_v)/3):
            for j in range(3):
                c.cd(j+1)
                h_v[3*i+j].Draw("COLZ")
                print out_v.at(3*i+j).meta().dump()

            c.Update()
            raw_input()
    for i in xrange(out_v.size()):
        detsplit.Append(out_v.at(i))
    out.set_id( io.event_id().run(), io.event_id().subrun(), io.event_id().event() )
    print "save entry"
    out.save_entry()

algo.printElapsedTime()
io.finalize()
out.finalize()

print "=========================================="

io2 = larcv.IOManager(larcv.IOManager.kREAD)
io2.add_in_file( "baka.root" )
io2.initialize()

for i in xrange(io2.get_n_entries()):
    io2.read_entry(i)
    evdetsplit = io2.get_data(larcv.kProductImage2D,"detsplit")
    print "entry ",i,": desplit images = ",evdetsplit.Image2DArray().size()

print "FIN"
