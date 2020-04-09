import os,sys,argparse

parser = argparse.ArgumentParser("Test MCPixelPGraph")
parser.add_argument("-ill", "--input-larlite",required=True,type=str,help="Input larlite file")
parser.add_argument("-ilcv","--input-larcv",required=True,type=str,help="Input LArCV file")
parser.add_argument("-adc", "--adc",type=str,default="wire",help="Name of tree with Wire ADC values [default: wire]")
parser.add_argument("-tb",  "--tick-backward",action='store_true',default=False,help="Input LArCV data is tick-backward [default: false]")
    
args = parser.parse_args()

import ROOT as rt
from larlite import larlite
from larcv import larcv
from ublarcvapp import ublarcvapp

rt.gROOT.ProcessLine( "gErrorIgnoreLevel = 3002;" )

"""
test script that demos the Neutrino Vertex class.
"""

rt.gStyle.SetOptStat(0)

ioll = larlite.storage_manager( larlite.storage_manager.kREAD )
ioll.add_in_filename(  args.input_larlite )
ioll.open()

if args.tick_backward:
    iolcv = larcv.IOManager( larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickBackward )
else:
    iolcv = larcv.IOManager( larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickForward )
iolcv.add_in_file( args.input_larcv )
iolcv.reverse_all_products()
iolcv.initialize()

nentries = iolcv.get_n_entries()
print "Number of entries: ",nentries

print "Start loop."
vtxutil = ublarcvapp.mctools.NeutrinoVertex
print vtxutil

c = rt.TCanvas("c","c",1200,1800)
c.Divide(1,3)

for ientry in xrange( nentries ):

    print 
    print "=========================="
    print "===[ EVENT ",ientry," ]==="
    ioll.go_to(ientry)
    iolcv.read_entry(ientry)


    ev_adc = iolcv.get_data( larcv.kProductImage2D, args.adc )
    print "number of images: ",ev_adc.Image2DArray().size()
    adc_v = ev_adc.Image2DArray()
    for p in xrange(adc_v.size()):
        print " image[",p,"] ",adc_v[p].meta().dump()

    # vertex in image coords
    vtxcoord = vtxutil.getImageCoords( ioll )
    vtxpos   = vtxutil.getPos3DwSCE( ioll )
    print "vertex pos: (",vtxpos[0],",",vtxpos[1],",",vtxpos[2],") tick=",vtxpos[3]
    print "image coords: (",vtxcoord[0],",",vtxcoord[1],",",vtxcoord[2],",",vtxcoord[3],")"

    intpc = True
    if vtxpos[0]<0 or vtxpos[0]>256 or vtxpos[1]>117 or vtxpos[1]<-117 or vtxpos[2]<0 or vtxpos[2]>1036:
        intpc = False
    print "intpc: ",intpc
    
    # make histogram
    hist_v = larcv.rootutils.as_th2d_v( adc_v, "hentry%d"%(ientry) )
    for ih in xrange(adc_v.size()):
        h = hist_v[ih]
        h.GetZaxis().SetRangeUser(0,100)

    # make vertex graph
    g_v = []
    for p in xrange(3):
        g = rt.TGraph(1)
        if intpc:
            g.SetMarkerStyle(20)
        else:
            g.SetMarkerStyle(24)
        g.SetPoint(0,vtxcoord[p],vtxcoord[3])        
        g.SetMarkerColor(rt.kRed)
        g_v.append(g)
        

    #draw canvas
    for p in xrange(3):
        c.cd(p+1)
        hist_v[p].Draw("colz")
        g_v[p].Draw("P")
    c.Update()

    print "[enter to continue]"
    raw_input()



print "=== FIN =="
