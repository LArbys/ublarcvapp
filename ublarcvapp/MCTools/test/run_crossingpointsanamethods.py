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
from ROOT import larutil
from ublarcvapp import ublarcvapp


rt.gROOT.ProcessLine( "gErrorIgnoreLevel = 3002;" )

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
xptana = ublarcvapp.mctools.CrossingPointsAnaMethods
sce = larutil.SpaceChargeMicroBooNE()

tmp = rt.TFile("temp.root","recreate")

c = rt.TCanvas("c","c",1200,1800)
c.Divide(1,3)

def tyz_from_mcstep( mcstep ):
    tick = xptana.getTick( mcstep, 4050.0, sce )
    tyz = []
    for p in xrange(0,3):
        vec = rt.TVector3()
        vec[0] = mcstep.X()
        vec[1] = mcstep.Y()
        vec[2] = mcstep.Z()
        w = larutil.Geometry.GetME().WireCoordinate( vec, p )
        tyz.append(w)
    tyz.append(tick)
    return tyz

    

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

    # larlite
    ev_mctrack = ioll.get_data( larlite.data.kMCTrack, "mcreco" )
    print "num of mctracks: ",ev_mctrack.size()

    xpt_list = []
    for itrack in xrange(ev_mctrack.size()):
        print "[mctrack ",itrack,"] "        
        mct = ev_mctrack.at(itrack)
        nsteps = mct.size()
        print " mcsteps: ",nsteps
        if nsteps==0:
            continue
        print " start: (",mct.front().X(),",",mct.front().Y(),",",mct.front().Z(),"t=",mct.front().T(),")"," : ",
        tyz_start = tyz_from_mcstep( mct.front() )
        print tyz_start
        print " end: (",mct.back().X(),",",mct.back().Y(),",",mct.back().Z(),"t=",mct.back().T(),")"
        tyz_end   = tyz_from_mcstep( mct.back() )

        if False:
            s_v = []
            for p in xrange(0,3):
                g = rt.TGraph(1)
                g.SetPoint(0,tyz_start[p],tyz_start[3])
                g.SetMarkerStyle(20)
                g.SetMarkerColor(rt.kRed)
                s_v.append(g)
            xpt_list.append(s_v)

            e_v = []
            for p in xrange(0,3):
                g = rt.TGraph(1)
                g.SetPoint(0,tyz_end[p],tyz_end[3])
                g.SetMarkerStyle(20)
                g.SetMarkerColor(rt.kBlue)
                e_v.append(g)
            xpt_list.append(e_v)
        
        start_imgcoord = xptana.getFirstStepPosInsideImage( mct, adc_v[0].meta(), 4050, True, 0.3, 0.0, sce )
        if start_imgcoord.size()>0:
            print " len(start_imgcoord): ",start_imgcoord.size()," (",start_imgcoord[0],",",start_imgcoord[1],",",start_imgcoord[2],",",start_imgcoord[3],")"
            g_v = []            
            for p in xrange(3):
                g = rt.TGraph(1)
                g.SetPoint( 0, start_imgcoord[p+1], adc_v[p].meta().pos_y(start_imgcoord[0]) )
                g.SetMarkerStyle(20)
                g.SetMarkerColor(rt.kMagenta)
                g_v.append(g)
            xpt_list.append( g_v )

        end_imgcoord = xptana.getFirstStepPosInsideImage( mct, adc_v[0].meta(), 4050, False, 0.3, 0.0, sce )
        if end_imgcoord.size()>0:
            print " len(end_imgcoord): ",end_imgcoord.size()," (",end_imgcoord[0],",",end_imgcoord[1],",",end_imgcoord[2],",",end_imgcoord[3],")"
            g_v = []            
            for p in xrange(3):
                g = rt.TGraph(1)
                g.SetPoint( 0, end_imgcoord[p+1], adc_v[p].meta().pos_y(end_imgcoord[0]) )
                g.SetMarkerStyle(20)
                g.SetMarkerColor(rt.kCyan)
                g_v.append(g)
            xpt_list.append( g_v )
            
        
    # make histogram
    hist_v = larcv.rootutils.as_th2d_v( adc_v, "hentry%d"%(ientry) )
    for ih in xrange(adc_v.size()):
        h = hist_v[ih]
        h.GetZaxis().SetRangeUser(0,100)

    #draw canvas
    for p in xrange(3):
        c.cd(p+1)
        hist_v[p].Draw("colz")
        for g_v in xpt_list:
            g_v[p].Draw("P")
    c.Update()

    print "[enter to continue]"
    raw_input()
    sys.exit(0)    


print "=== FIN =="
