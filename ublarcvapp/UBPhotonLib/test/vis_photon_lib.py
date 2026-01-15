import os,sys
import ROOT as rt
from ROOT import std
from ublarcvapp import ublarcvapp
import lardly
from lardly.ubdl.pmtpos import getPMTPosByOpChannel,getPMTPosByOpDet

photlib = ublarcvapp.ubphotonlib.UBPhotonLib.getPhotonLib()
print(photlib)

rt.gStyle.SetOptStat(0)

def make_yz_hist( xslice, opch, histname="hvis2d_sliceyz", nbinsz=100, nbinsy=20 ):

    h2d = rt.TH2D( histname, "", nbinsz, 0, 1040.0, nbinsy, -120, 120 )
    pos = std.vector("float")(3)    
    for i in range(nbinsz):
        print("get vis for z-bin[",i,"]")
        for j in range(nbinsy):
            pos[0] = xslice
            pos[1] = h2d.GetYaxis().GetBinCenter(j+1)
            pos[2] = h2d.GetXaxis().GetBinCenter(i+1)
            vis = photlib.getVisibility( pos, opch )
            voxcoord = photlib.getVoxelCoords( pos )
            h2d.SetBinContent( i+1, j+1, vis )
            #print("get vis at pos=",pos," for opch=",opch," vis=",vis," voxcoord=",voxcoord)
    return h2d


if __name__=="__main__":

    iopch = 15
    xpos = 50.0
    temp = rt.TFile("test.root","recreate")
    h2d = make_yz_hist(xpos, iopch, nbinsz=200, nbinsy=40)
    c = rt.TCanvas("c","c",2000,600)
    c.Draw()
    h2d.Draw("colz")
    h2d.SetTitle("Photon Visibility vs. (y,z) with x=%.1f cm for OpDet=%d; z position (cm); y position (cm)"%(xpos,iopch))

    # drop opchannel pos
    opch_v = []
    opch_label_v = []
    for ich in range(32):
        pos = getPMTPosByOpDet( ich, use_v4_geom=False )
        #pos = getPMTPosByOpChannel( ich, use_v4_geom=False )        
        circ = rt.TEllipse( pos[2], pos[1], 15.2, 15.2 )
        circ.SetFillStyle(0)
        oplabel = rt.TText(pos[2]-10.0,pos[1]-5.0,"%02d"%(ich))        
        if ich==iopch:
            circ.SetLineColor(rt.kRed)
            oplabel.SetTextColor(rt.kRed)
        circ.Draw()

        oplabel.Draw()
        opch_v.append( circ )
        opch_label_v.append(oplabel)
    tt = rt.TText(826.0, 99.0, "Labels are the OpDetID")
    tt.Draw()
    opch_label_v.append(tt)
    circ.Draw()
    
    c.Update()
    print("[Enter] to continue")
    input()
    
    
