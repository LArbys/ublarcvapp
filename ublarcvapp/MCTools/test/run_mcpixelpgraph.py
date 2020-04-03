import os,sys,argparse

parser = argparse.ArgumentParser("Test MCPixelPGraph")
parser.add_argument("-ill", "--input-larlite",required=True,type=str,help="Input larlite file")
parser.add_argument("-ilcv","--input-larcv",required=True,type=str,help="Input LArCV file")
parser.add_argument("-adc", "--adc",type=str,default="wire",help="Name of tree with Wire ADC values [default: wire]")
parser.add_argument("-tb",  "--tick-backward",action='store_true',default=False,help="Input LArCV data is tick-backward [default: false]")
args = parser.parse_args()

import ROOT as rt
from larcv import larcv
from larlite import larlite
from ublarcvapp import ublarcvapp

"""
test script that demos the MCPixelPGraph class.
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
nentries = 10

print "Start loop."

mcpg = ublarcvapp.mctools.MCPixelPGraph()
mcpg.set_adc_treename( args.adc )

tmp = rt.TFile("temp.root","recreate")

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
    
    # make histogram
    hist_v = larcv.rootutils.as_th2d_v( adc_v, "hentry%d"%(ientry) )
    for ih in xrange(adc_v.size()):
        h = hist_v[ih]
        h.GetZaxis().SetRangeUser(0,100)

    mcpg.buildgraph( iolcv, ioll )
    #mcpg.printAllNodeInfo()
    mcpg.printGraph()

    primaries = mcpg.getPrimaryParticles()    

    # get primary electron, make tgraph of pixels
    graph_v = []
    for i in xrange(primaries.size()):
        node = primaries.at(i)
        #print "primary pid[",node.pid,"]"
        if node.pid in [11,2212,13,-13]:
            #print "making tgraph for pid=",node.pid
            e_v = []
            for p in xrange(3):
                if node.pix_vv[p].size()==0:
                    e_v.append(None)
                    continue
                g = rt.TGraph( node.pix_vv[p].size()/2 )
                for j in xrange( node.pix_vv[p].size()/2 ):
                    g.SetPoint(j, node.pix_vv[p][2*j+1], node.pix_vv[p][2*j] ) # wire, tick
                g.SetMarkerStyle(20)
                g.SetMarkerSize(0.5)                
                if node.pid==11:
                    if node.origin==1:
                        g.SetMarkerColor(rt.kRed)
                elif node.pid in [13,-13]:
                    if node.origin==2:
                        g.SetMarkerColor(rt.kGreen)
                    elif node.origin==1:
                        g.SetMarkerColor(rt.kMagenta)
                elif node.pid in [2212]:
                    if node.origin==1:                    
                        g.SetMarkerColor(rt.kBlue)
                e_v.append(g)
            graph_v.append(e_v)

    #draw canvas
    for p in xrange(3):
        c.cd(p+1)
        hist_v[p].Draw("colz")
        for e_v in graph_v:
            if e_v[p] is not None:
                e_v[p].Draw("P")
    c.Update()

    print "[enter to continue]"
    raw_input()    


print "=== FIN =="
