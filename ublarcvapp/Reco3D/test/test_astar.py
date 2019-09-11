import os,sys,time

import ROOT as rt
from larcv import larcv
from ublarcvapp import ublarcvapp
larcv.load_rootutil()
from ROOT import std

from ctypes import c_int

rt.gStyle.SetOptStat(0)

infile = "out_larcv_test.root"

tickdir = larcv.IOManager.kTickForward
#tickdir = larcv.IOManager.kTickBackward

# CONFIG
algoconfig = ublarcvapp.reco3d.AStar3DAlgoConfig()
algoconfig.store_score_image = True

algo = ublarcvapp.reco3d.AStar3DAlgo(algoconfig)
algo.setVerbose(1)

io = larcv.IOManager(larcv.IOManager.kREAD,"io",tickdir)
io.add_in_file( infile )
io.initialize()

badchmaker = ublarcvapp.EmptyChannelAlgo()

entry = 0

io.read_entry( entry )

ev_wire = io.get_data( larcv.kProductImage2D, "wire" )
wire_v  = ev_wire.Image2DArray()
meta    = wire_v.at(0).meta()
print "num images: ",wire_v.size()

ev_chstatus = io.get_data( larcv.kProductChStatus, "wire" )
badch_v   = badchmaker.makeBadChImage(4,3,2400,1008*6,3456,6,1,ev_chstatus)
for p in xrange(3):
    print "badch meta: ",badch_v.at(p).meta().dump()
    

# test points
# entry zero
start_pts = [ (324,8088),
              (953,8088),
              (607,8088) ] # wire,tick
end_pts = [ (755,6348),
            (103,6348),
            (192,6348)]

start_v = std.vector("int")()
end_v   = std.vector("int")()
start_tick = start_pts[0][1]
end_tick   = end_pts[0][1]
for p in xrange(3):
    start_v.push_back( start_pts[p][0] )
    end_v.push_back( end_pts[p][0] )

# call algo
dtalgo = time.time()
reached = c_int()
path = algo.downsampleAndFindPath( 16, wire_v, badch_v, badch_v,
                                   meta.row(start_tick), meta.row(end_tick),
                                   start_v, end_v, reached, 2 );
dtalgo = time.time()-dtalgo
print "Algo completed in ",dtalgo," seconds"


tgraph_v = []
for p in xrange(3):
    g = rt.TGraph(path.size())
    for i in xrange(path.size()):
        node = path.at(i)
        tick = node.tyz.at(0)
        wire = node.cols.at(p)
        g.SetPoint(i,wire,tick)
    g.SetMarkerStyle(20)
    g.SetMarkerColor(rt.kRed)
    tgraph_v.append(g)


scoreimg_v = algo.getScoreImages()
print "Astar visualization images: ",scoreimg_v.size()

hwire = [ larcv.as_th2d( wire_v[x], "hwire_p%d"%(x) ) for x in xrange(3) ]
    
cscore = rt.TCanvas("cscore","cscore",1800,1500)
cscore.Divide(3,3)
hscore = [ larcv.as_th2d( scoreimg_v.at(x), "hscore_p%d"%(x) ) for x in xrange(scoreimg_v.size()) ]
for p in xrange(3):
    g = tgraph_v[p]
    
    cscore.cd(3*p+1)
    #hscore[p].Draw("colz") # ADC image
    hwire[p].Draw("colz")
    hwire[p].SetTitle("Plane %d: ADC;wire;tick"%(p))
    hscore[p].GetXaxis().SetRangeUser( min( start_v[p], end_v[p])-50, max(start_v[p], end_v[p])+50 )
    hscore[p].GetYaxis().SetRangeUser( min( start_tick, end_tick)-50, max(start_tick, end_tick)+50 )    
    g.Draw("LP")
    
    cscore.cd(3*p+2).SetLogz(1)
    hscore[3+p].Draw("colz") # ADC image
    hscore[3+p].SetTitle("Plane %d: F-score;wire;tick"%(p))    
    hscore[3+p].GetXaxis().SetRangeUser( min( start_v[p], end_v[p])-50, max(start_v[p], end_v[p])+50 )
    hscore[p].GetYaxis().SetRangeUser( min( start_tick, end_tick)-50, max(start_tick, end_tick)+50 )        
    g.Draw("LP")    

    cscore.cd(3*p+3).SetLogz(1)
    hscore[6+p].Draw("colz")
    hscore[6+p].SetTitle("Plane %d: G-score;wire;tick"%(p))        
    hscore[6+p].GetXaxis().SetRangeUser( min( start_v[p], end_v[p])-50, max(start_v[p], end_v[p])+50 )
    hscore[6+p].GetYaxis().SetRangeUser( min( start_tick, end_tick)-50, max(start_tick, end_tick)+50 ) 
    g.Draw("LP")
    
cscore.Update()

#c = rt.TCanvas("c","c",800,600)
#c.Draw()

#hist_v[2].Draw("colz")

#c.Update()

raw_input()
