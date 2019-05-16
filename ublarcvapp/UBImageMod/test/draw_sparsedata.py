import os,sys
import ROOT as rt
from larcv import larcv

rt.gStyle.SetOptStat(0)

io = larcv.IOManager(larcv.IOManager.kREAD)
io.add_in_file(sys.argv[1])
io.initialize()

for ientry in xrange(io.get_n_entries()):
    io.read_entry(ientry)
    ev_sparse = io.get_data(larcv.kProductSparseImage,sys.argv[2])

    c = {}
    h = {}
    for ifeat in xrange(3):
        c[ifeat] = rt.TCanvas("c%d"%(ifeat),"Image %d"%(ifeat),1200,1000)
        h[ifeat] = rt.TH2D("h%d"%(ifeat),"",496,0,496,512,2400,2400+512)
        sparseimg=ev_sparse.at(ifeat)

        nvalues = sparseimg.pixellist().size()
        print nvalues
        npts    = nvalues/3
        for ipt in xrange(npts):
            start = ipt*(3)
            row = int(sparseimg.pixellist()[start])
            col = int(sparseimg.pixellist()[start+1])

            h[ifeat].SetBinContent( col+1 , row+1,
                                    sparseimg.pixellist()[start+2] )
            if sparseimg.pixellist()[start+2] ==0:
                h[ifeat].SetBinContent( col+1, row+1, .1 )

    for ifeat in xrange(3):
        c[ifeat].cd()
        c[ifeat].Draw()
        h[ifeat].Draw("COLZ")
        c[ifeat].Update()
            #c[ifeat].SaveAs("dumpedimages/test_sparseimg_img%d.pdf"%(ifeat))
    print "[ENTER] for next entry"
    raw_input()
