import os,sys

import ROOT as rt
from ROOT import std

from ublarcvapp import ublarcvapp


""" Test DBScan implementation """

# Create gaussians
gaus_clusters = [ (2.0,2.0,0.5),
                  (5.0,1.0,0.2),
                  (0.5,6.0,1.0) ]
npts_per_gaus = [50, 100, 100]
color_v = [ rt.kRed,
            rt.kBlue,
            rt.kGreen ]

rand = rt.TRandom3()

truth_graphs = []
data_v = std.vector("std::vector<float>")()
for npts,pars,color in zip(npts_per_gaus,gaus_clusters,color_v):
    g = rt.TGraph(npts)
    for ipt in xrange(npts):
        x = rand.Gaus( pars[0], pars[2] )
        y = rand.Gaus( pars[1], pars[2] )
        g.SetPoint(ipt,x,y)
        pt = std.vector("float")(3,0)
        pt[0] = x
        pt[1] = y
        data_v.push_back( pt )
    g.SetMarkerStyle(20)
    g.SetMarkerColor( color )
    truth_graphs.append(g)

print "generated {} points".format( data_v.size() )

# algo
dbscan = ublarcvapp.dbscan.DBScan()
print dbscan

cluster_v = dbscan.makeCluster3f( 0.2, 3, 5, data_v )

reco_graphs = []
for iclust in xrange(cluster_v.size()):
    clust = cluster_v.at(iclust)
    print "cluster [{}] has {} points".format(iclust,clust.size())
    g = rt.TGraph( clust.size() )
    for n,ihit in enumerate(xrange(clust.size())):
        index = clust[ihit]
        x = data_v[index][0]
        y = data_v[index][1]
        g.SetPoint(n,x,y)
    if iclust+1<cluster_v.size():
        g.SetMarkerColor( int(rand.Uniform(50))+1 )
    else:
        g.SetMarkerColor(rt.kBlack)
    g.SetMarkerStyle(20)
    reco_graphs.append(g)
    

print "generated {} clusters from data points".format(cluster_v.size())
    
h = rt.TH2D("h","",100,0,10, 100, 0, 10 )

c = rt.TCanvas("c","c",1600,800)
c.Draw()
c.Divide(2,1)

c.cd(1)
h.Draw()
for g in truth_graphs:
    g.Draw("p")

c.cd(2)
h.Draw()
for g in reco_graphs:
    g.Draw("p")

c.Update()


print "[enter] to exit"
raw_input()
