from __future__ import print_function

import os,sys
from math import fabs

import ROOT as rt
from ROOT import std
from larlite import larlite
from larcv import larcv
from ublarcvapp import ublarcvapp
larcv.load_rootutil()

rt.gStyle.SetOptStat(0)


#inputfile = "../../../../out_larcv_test.root"
#inputfile = "testset1/out_larcv_test.root"
#inputfiles = ["out_larcv_test.root"] # forward example
inputfiles = ["../bin/testset2/supera-Run007704-SubRun000023.root", "../bin/testset2/mrcnn-mcc9_v13_nueintrinsics_overlay_run1-Run007704-SubRun000023.root"] # backward example

llinputfiles = ["../bin/testset2/opreco-Run007704-SubRun000023.root"]

dltagger = ublarcvapp.dltagger.DLTagger()
dltagger.set_verbosity(0)

io = larcv.IOManager(larcv.IOManager.kREAD, "IO", larcv.IOManager.kTickBackward )
#io = larcv.IOManager(larcv.IOManager.kREAD, "IO" )
for inputfile in inputfiles:
    io.add_in_file( inputfile )
io.initialize()

ioll = larlite.storage_manager(larlite.storage_manager.kREAD)
for l in llinputfiles:
    ioll.add_in_filename( l )
ioll.open()

nentries = io.get_n_entries()

for ientry in xrange(nentries):

    io.read_entry(ientry)
    ioll.go_to(ientry)
    
    ev_wire  = io.get_data(larcv.kProductImage2D, "wire" )
    ev_masks = io.get_data(larcv.kProductClusterMask, "mrcnn_masks" )
    ev_chstatus = io.get_data(larcv.kProductChStatus, "wire" )
    ev_opflash  = ioll.get_data(larlite.data.kOpFlash, "simpleFlashBeam")

    mask_vv = ev_masks.as_vector()
    print("Number of masks: ",[mask_vv.at(x).size() for x in range(3)])
    hwire_v = [ larcv.as_th2d( ev_wire.as_vector().at(p), "hwire_p%d"%(p) ) for p in xrange(3) ]
    meta_v  = [ ev_wire.as_vector().at(p).meta() for p in xrange(3) ]

    dltagger.runTagger( ev_wire.Image2DArray(),
                        ev_chstatus,
                        ev_opflash,
                        mask_vv )

    tagged_v = std.vector("larcv::Image2D")()
    nottagged_v = std.vector("larcv::Image2D")()    
    dltagger.transferImages(tagged_v,nottagged_v)

    htagged_v = [ larcv.as_th2d( tagged_v.at(x), "htagged_p%d"%(x) ) for x in range(3) ]
    hnottagged_v = [ larcv.as_th2d( nottagged_v.at(x), "hnottagged_p%d"%(x) ) for x in range(3) ]    

    c = rt.TCanvas("c","Tagged", 1500,1200)
    c.Divide(1,3)
    c.Draw()
    for p in xrange(3):
        c.cd(p+1)
        htagged_v[p].Draw("colz")
    c.Update()
    c.SaveAs("testset1/thrumuimg_ientry%02d.png"%(ientry))

io.finalize()

raw_input()
