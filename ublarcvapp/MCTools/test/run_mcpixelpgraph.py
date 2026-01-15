#!/bin/env python3
from __future__ import print_function
import os,sys,argparse

parser = argparse.ArgumentParser("Test MCPixelPGraph")
parser.add_argument("-ill", "--input-larlite",required=True,type=str,help="Input larlite file")
parser.add_argument("-ilcv","--input-larcv",required=False,default=None,type=str,help="Input LArCV file")
parser.add_argument("-adc", "--adc",type=str,default="wire",help="Name of tree with Wire ADC values [default: wire]")
parser.add_argument("-tb",  "--tick-backward",action='store_true',default=False,help="Input LArCV data is tick-backward [default: false]")
parser.add_argument("-d",   "--debug", action='store_true', default=False, help="Run in debug mode")
parser.add_argument('-n', "--nentries", required=False, type=int, default=-1, help="number of entries to run")
parser.add_argument('-e', "--entry", required=False, type=int, default=-1, help="start at given entry number")
parser.add_argument('-v', "--vis", required=False, default=False, action='store_true', help="if flag provided, will visualize event")
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
ioll.set_data_to_read( "mctrack",  "mcreco" ) # larlite.data.kMCTrack
ioll.set_data_to_read( "mcshower", "mcreco" ) # larlite.data.kMCShower
ioll.set_data_to_read( "mctruth", "generator" ) # larlite.data.kMCTruth
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

print("HAS LARCV: ",HAS_LARCV)
tot_nentries = ioll.get_entries()
start_entry = 0
print("Number of entries: ",tot_nentries)
if args.entry > 0:
    start_entry = args.entry
if args.nentries>0:
    nentries = args.nentries
else:
    nentries = tot_nentries
end_entry = start_entry + nentries
if end_entry>tot_nentries:
    end_entry = tot_nentries
    


print("Start loop.")

mcpg = ublarcvapp.mctools.MCPixelPGraph()
mcpg_nu = ublarcvapp.mctools.MCPixelPGraph()
if args.debug:
    mcpg.set_verbosity( "debug" )
    mcpg_nu.set_verbosity( "debug" )
else:
    mcpg.set_verbosity( "info" )
    mcpg_nu.set_verbosity( "info" )

if HAS_LARCV:
    mcpg.set_adc_treename( args.adc )

if HAS_LARCV:
    tmp = rt.TFile("temp.root","recreate")
    c = rt.TCanvas("c","c",2100,1500)
    c.Divide(3,3)

#print("[ENTER] to Start")
#input()

for ientry in range( start_entry, end_entry ):

    print() 
    print("==========================")
    print("===[ EVENT ",ientry," ]===")
    ioll.go_to(ientry)

    mcpg.clear()
    mcpg_nu.clear()
    mcpg.set_cluster_neutrino_particles(False)
    mcpg_nu.set_cluster_neutrino_particles(True)

    if not HAS_LARCV:
        mcpg.buildgraphonly( ioll )
        mcpg_nu.buildgraphonly( ioll )
    else:
        print("HAS_LARCV: get images needed for pixel matching to particles")
        #mcpg.buildgraphonly( ioll )        
        iolcv.read_entry(ientry)
        ev_adc = iolcv.get_data( "image2d", args.adc )
        ev_instance = iolcv.get_data( "image2d", "instance" )
        ev_ancestor = iolcv.get_data( "image2d", "ancestor" )
        ev_segment  = iolcv.get_data( "image2d", "segment" )
        ev_larflow  = iolcv.get_data( "image2d", "larflow" )
        print("number of adc images: ",ev_adc.Image2DArray().size())
        print("number of larflow images: ",ev_larflow.Image2DArray().size())
        adc_v = ev_adc.Image2DArray()
        for p in range(adc_v.size()):
            print(" image[",p,"] ",adc_v[p].meta().dump())        
        mcpg_nu.buildgraph( iolcv, ioll )        
        

    photon_starts = {}            
    if HAS_LARCV and False:
        for i in range(mcpg_nu.node_v.size()):
            node = mcpg_nu.node_v.at(i)
            if node.pid!=22:
                continue;
            print("running fixingPhotonStartPoints(...) for photons")
            gamma_start = mcpg_nu.fixingPhotonStartPoints( node,
                                                           ev_instance.as_vector(),
                                                           ev_ancestor.as_vector(),
                                                           ev_adc.as_vector(),
                                                           ev_larflow.as_vector() )
            photon_starts[node.tid] = gamma_start
            
    ## isolate events
    stop = False
    for inode in range(mcpg_nu.node_v.size()):
        node = mcpg_nu.node_v.at(inode)
        if node.pid in [22]:
            npixsum = 0
            for p in range(node.pix_vv.size()):
                npixsum += node.pix_vv[p].size()
            if npixsum>1:
                stop=True
    print("STOP TO VISUALIZE/DUMP  MCPG INFO: ",stop)
    if not stop:
        continue
    
    if args.debug or args.vis:
        if not HAS_LARCV:
            print("================================================")
            print("PARSING PARTICLE GRAPH ONLY: No pixel matching  ")
            print("================================================")        
            print("ALL NODE INFO [NO NU] --------------------------")
            mcpg.printAllNodeInfo()
            print(" -----------------------------------------------")
            print("CONSTRUCTED PARTICLE GRAPH [NO NU]")
            mcpg.printGraph(0,False)
            print("====================================================")
        else:
            #print("================================================")
            #print("PARSING PARTICLE GRAPH ONLY: No pixel matching  ")
            #print("================================================")        
            #print("CONSTRUCTED PARTICLE GRAPH [NO NU]")
            #mcpg.printGraph(0,False)
            #print("====================================================")
            
            print("====================================================")
            print("CONSTRUCTED PARTICLE GRAPH [WITH NU VERTEX GROUPING]")
            mcpg_nu.printGraph(0,False)
            print("====================================================")
            
    
        # make histogram
        marker_v = []
        #hist_v = mcpg_nu.makeTH2D( "hentry%d"%(ientry) )
        hist_v = larcv.rootutils.as_th2d_v( ev_adc.as_vector(), "hentry%d"%(ientry) )
        segment_v = larcv.rootutils.as_th2d_v( ev_segment.as_vector(), "hsegment%d"%(ientry))
        instance_v = larcv.rootutils.as_th2d_v( ev_instance.as_vector(), "hinstance%d"%(ientry))
        ancestor_v = larcv.rootutils.as_th2d_v( ev_ancestor.as_vector(), "hancestor%d"%(ientry))

        c.cd()
        c.Draw()
        c.Update()
        for ih in range(hist_v.size()):
            c.cd(3*ih+0+1)
            h = hist_v[ih]
            h.Draw("colz")            
            h.GetZaxis().SetRangeUser(0,300)

            for i in range(mcpg_nu.node_v.size()):
                node = mcpg_nu.node_v.at(i)
                if node.pid not in [22,11,-11]:
                    continue
                
                x = node.imgpos4_edep[ih]
                y = node.imgpos4_edep[3]
                #print("nodeidx=",node.nodeidx," pid=",node.pid," (x,y)=",(x,y))
                m = rt.TMarker( x, y, 4 )
                m.SetNDC(False)
                m.SetMarkerColor(rt.kMagenta)
                m.Draw()
                marker_v.append(m)

                if node.tid in photon_starts:
                    x2 = photon_starts[node.tid][4+ih]
                    y2 = photon_starts[node.tid][4+3]
                    m2 = rt.TMarker( x2, y2, 4 )
                    m2.SetNDC(False)
                    m2.SetMarkerColor(rt.kRed)
                    m2.Draw()
                    marker_v.append(m2)
                
                x4 = node.imgpos4_start[ih]
                y4 = node.imgpos4_start[3]
                m4 = rt.TMarker( x4, y4, 5 )
                m4.SetNDC(False)
                m4.SetMarkerSize(3)
                m4.SetMarkerColor(rt.kBlack)
                #print("x4: ",(x4,y4))
                m4.Draw()
                marker_v.append(m4)

                #if ih==2:
                #    print("DUMP PIXELS [plane=2] for node.nodeidx=",node.nodeidx," trackid=",node.tid," pid=",node.pid)
                pix_v = node.pix_vv.at(ih)
                npix = int(pix_v.size()/2)
                for ipix in range(npix):
                    xbin = int(pix_v.at(2*ipix+1))
                    ybin = int((pix_v.at(2*ipix)-2400)/6.0)
                    #if ih==2:
                    #    print(" ipix[",ipix,"]: (xbin,ybin)=",(xbin,ybin)," (wire,tick)=",(pix_v.at(2*ipix+1),pix_v.at(2*ipix)))
                    ancestor_v[ih].SetBinContent( xbin+1, ybin+1, node.tid )
                
            c.cd(3*ih+1+1)
            #instance_v[ih].SetTitle("instance plane[%d]"%(ih))            
            #instance_v[ih].Draw("colz")
            segment_v[ih].SetTitle("segment plane[%d]"%(ih))            
            segment_v[ih].Draw("colz")
            c.cd(3*ih+2+1)
            ancestor_v[ih].SetTitle("ancestor plane[%d]"%(ih))
            ancestor_v[ih].Draw("colz")

            c.Update()
        # end of loop
        c.Update()

        # draw individual showers
        cpart = []
        for inode in range(mcpg_nu.node_v.size()):
            node = mcpg_nu.node_v.at(inode)
            if node.pid not in [22]:
                continue
            print("Make canvas of particle, tid=",node.tid)
            pixsum_v = mcpg_nu.getTruePhotonTrunkPlanePixelSums( node.tid )
            if pixsum_v.size()>=3:
                print(" pixelsum (",pixsum_v[0],", ",pixsum_v[1],", ",pixsum_v[2],")")
            meta = ev_adc.as_vector().at(2).meta()
            hpart = rt.TH2D("hnode%d"%(inode),"",3456,0,3456,1008,2400,2400+6*1008)
            pix_v = node.pix_vv.at(ih)
            npix = int(pix_v.size()/2)
            for ipix in range(npix):
                xbin = int(pix_v.at(2*ipix+1))
                ybin = int((pix_v.at(2*ipix)-2400)/6.0)
                hpart.SetBinContent( xbin+1, ybin+1, 1 )
            cpt = rt.TCanvas("cnode%d"%(inode),"node[%d] pid[%d] tid[%d]"%(inode,node.pid,node.tid),800,600)
            hpart.GetZaxis().SetRangeUser(0,10.0)
            hpart.Draw("colz")
            cpt.Update()
            cpart.append(cpt)
            cpart.append(hpart)
        
    
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
