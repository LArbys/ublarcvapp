from __future__ import print_function
import os,sys
import ROOT as rt
from larlite import larlite
from larcv import larcv
from ublarcvapp import ublarcvapp



"""
Test the class TrackImageMask
"""

PLOT_ME = True

# make a dummy larcv image
# first define a meta, looks like a UB image
meta = larcv.ImageMeta( 3456, 1008*6.0, 1008, 3456, 0, 2400.0, 0 )
img = larcv.Image2D( meta )
mask = larcv.Image2D( meta )
smin = larcv.Image2D( meta )
smax = larcv.Image2D( meta )

# define a larlite track
start = rt.TVector3( 20, -50, 250 )
end   = rt.TVector3( 200, 80, 820 )
vdir = end-start
lltrack = larlite.track()
lltrack.reserve(2)
lltrack.add_vertex(start)
lltrack.add_vertex(end)
lltrack.add_direction(vdir)
lltrack.add_direction(vdir)

algo = ublarcvapp.ubimagemod.TrackImageMask()
algo.set_verbosity( larcv.msg.kDEBUG )
algo.show_timing = True
print(algo)


algo.maskTrack( lltrack, img, mask, -1.0, 2, 2, 0.2 )
algo.labelTrackPath( lltrack, img, smin, smax )

if PLOT_ME:

    rt.gStyle.SetOptStat(0)
    c = rt.TCanvas("c","c",1200,2400)
    c.Divide(1,4)
    himg = larcv.rootutils.as_th2d(img,"img")
    hmask = larcv.rootutils.as_th2d(mask,"mask")
    hpath_min = larcv.rootutils.as_th2d(smin,"smin")
    hpath_max = larcv.rootutils.as_th2d(smax,"smax")
    c.cd(1)
    himg.Draw("colz")
    c.cd(2)
    hmask.Draw("colz")
    c.cd(3)
    hpath_min.Draw("colz")
    c.cd(4)
    hpath_max.Draw("colz")
    print("[ENTER] to continue")
    raw_input()

