from __future__ import print_function
import os,sys,argparse

parser = argparse.ArgumentParser("Test MCPixelPGraph")
parser.add_argument("-ill", "--input-larlite",required=True,type=str,help="Input larlite file")
parser.add_argument("-ilcv","--input-larcv",required=False,default=None,type=str,help="Input LArCV file")
parser.add_argument("-adc", "--adc",type=str,default="wire",help="Name of tree with Wire ADC values [default: wire]")
parser.add_argument("-tb",  "--tick-backward",action='store_true',default=False,help="Input LArCV data is tick-backward [default: false]")
parser.add_argument("-d",   "--debug", action='store_true', default=False, help="Run in debug mode")
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
ioll.set_data_to_read( larlite.data.kMCTrack,  "mcreco" )
ioll.set_data_to_read( larlite.data.kMCShower, "mcreco" )
ioll.set_data_to_read( larlite.data.kMCTruth,  "generator" )
ioll.add_in_filename(  args.input_larlite )
ioll.open()

if args.input_larcv is not None:
    # has larcv images to get pixels
    HAS_LARCV = True
    if args.tick_backward:
        iolcv = larcv.IOManager( larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickBackward )
    else:
        iolcv = larcv.IOManager( larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickForward )
    iolcv.add_in_file( args.input_larcv )
    iolcv.reverse_all_products()
    iolcv.initialize()
else:
    HAS_LARCV = False

nentries = ioll.get_entries()
print("Number of entries: ",nentries)
#nentries = 10

print("Start loop.")

mcpg = ublarcvapp.mctools.MCPixelPGraph()
mcpg_nu = ublarcvapp.mctools.MCPixelPGraph()
if args.debug:
    mcpg.set_verbosity( larcv.msg.kDEBUG )
    mcpg_nu.set_verbosity( larcv.msg.kDEBUG )    
else:
    mcpg.set_verbosity( larcv.msg.kINFO )
    mcpg_nu.set_verbosity( larcv.msg.kINFO )    

if HAS_LARCV:
    mcpg.set_adc_treename( args.adc )

if HAS_LARCV:
    tmp = rt.TFile("temp.root","recreate")
    c = rt.TCanvas("c","c",1200,1800)
    c.Divide(1,3)

print("[ENTER to Start")
input()

for ientry in range( nentries ):

    print() 
    print("==========================")
    print("===[ EVENT ",ientry," ]===")
    ioll.go_to(ientry)

    mcpg.clear()
    mcpg_nu.clear()
    mcpg.set_cluster_neutrino_particles(False)
    mcpg_nu.set_cluster_neutrino_particles(True)
    
    if HAS_LARCV:
        iolcv.read_entry(ientry)
        ev_adc = iolcv.get_data( larcv.kProductImage2D, args.adc )
        print("number of images: ",ev_adc.Image2DArray().size())
        adc_v = ev_adc.Image2DArray()
        for p in range(adc_v.size()):
            print(" image[",p,"] ",adc_v[p].meta().dump())
    
        # make histogram
        hist_v = larcv.rootutils.as_th2d_v( adc_v, "hentry%d"%(ientry) )
        for ih in range(adc_v.size()):
            h = hist_v[ih]
            h.GetZaxis().SetRangeUser(0,100)

        mcpg.buildgraph( iolcv, ioll )
    else:
        mcpg.buildgraphonly( ioll )
        mcpg_nu.buildgraphonly( ioll )

    if args.debug:
        print("ALL NODE INFO [NO NU]======================================")
        mcpg.printAllNodeInfo()
        print("====================================================")
        print("CONSTRUCTED PARTICLE GRAPH [NO NU]")
        mcpg.printGraph(0,False)
        print("====================================================")
        print("====================================================")
        print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
        print("====================================================")
    
    print("====================================================")
    print("ALL NODE INFO [WITH NU VERTEX GROUPING]======================================")
    mcpg_nu.printAllNodeInfo()
    print("====================================================")
    print("CONSTRUCTED PARTICLE GRAPH [WITH NU VERTEX GROUPING]")
    mcpg_nu.printGraph(0,False)
    print("====================================================")
    
    print("[ENTER] to continue")
    input()
    if True:
        continue

    #primaries = mcpg.getPrimaryParticles()
    primaries = mcpg.node_v

    # get primary electron, make tgraph of pixels
    graph_v = []
    bbox_v  = []
    for i in range(primaries.size()):
        node = primaries.at(i)
        print("primary pid[",node.pid,"]")
        if node.pid in [11,2212,13,-13,22,211,-211]:
            print("  making tgraph for pid=",node.pid)
            e_v = []
            bb_v = []
            for p in range(3):
                if node.pix_vv[p].size()==0:
                    e_v.append(None)
                    bb_v.append(None)
                    continue
                g = rt.TGraph( node.pix_vv[p].size()/2 )

                bb = rt.TBox( node.plane_bbox_twHW_vv[p][1]-node.plane_bbox_twHW_vv[p][3],
                              node.plane_bbox_twHW_vv[p][0]-node.plane_bbox_twHW_vv[p][2],
                              node.plane_bbox_twHW_vv[p][1]+node.plane_bbox_twHW_vv[p][3],
                              node.plane_bbox_twHW_vv[p][0]+node.plane_bbox_twHW_vv[p][2] )
                bb.SetFillStyle(0)
                bb.SetLineWidth(2)
                
                for j in range( node.pix_vv[p].size()/2 ):
                    g.SetPoint(j, node.pix_vv[p][2*j+1], node.pix_vv[p][2*j] ) # wire, tick
                g.SetMarkerStyle(20)
                g.SetMarkerSize(0.5)                
                if node.pid==11:
                    if node.origin==1:
                        g.SetMarkerColor(rt.kRed)
                        bb.SetLineColor(rt.kRed)
                elif node.pid in [13,-13]:
                    if node.origin==2:
                        g.SetMarkerColor(rt.kGreen)
                        bb.SetLineColor(rt.kGreen)                        
                    elif node.origin==1:
                        g.SetMarkerColor(rt.kMagenta)
                        bb.SetLineColor(rt.kMagenta)                        
                elif node.pid in [2212]:
                    if node.origin==1:                    
                        g.SetMarkerColor(rt.kBlue)
                        bb.SetLineColor(rt.kBlue)                        
                elif node.pid in [22]:
                    g.SetMarkerColor(rt.kOrange)
                    bb.SetLineColor(rt.kOrange)                    
                elif node.pid in [211,-211]:
                    g.SetMarkerColor(rt.kViolet)
                    bb.SetLineColor(rt.kViolet)                    
                e_v.append(g)
                bb_v.append(bb)
            graph_v.append(e_v)
            bbox_v.append(bb_v)

    print("num graphs: ",len(graph_v))
    
    #draw canvas
    for p in range(3):
        c.cd(p+1)
        hist_v[p].Draw("colz")
        for e_v in graph_v:
            if e_v[p] is not None:
                e_v[p].Draw("P")
        for bb_v in bbox_v:
            if bb_v[p] is not None:
                bb_v[p].Draw()
    c.Update()

    print("[enter to continue]")
    raw_input()    


print("=== FIN ==")
