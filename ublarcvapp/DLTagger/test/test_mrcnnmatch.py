from __future__ import print_function

import os,sys

import ROOT as rt
from ROOT import std
from larcv import larcv
from ublarcvapp import ublarcvapp
larcv.load_rootutil()

rt.gStyle.SetOptStat(0)

inputfile = "out_larcv_test.root"

matchalgo = ublarcvapp.dltagger.MRCNNMatch()

io = larcv.IOManager(larcv.IOManager.kREAD, "IO" )
io.add_in_file( inputfile )
io.initialize()

io.read_entry(0)

ev_wire  = io.get_data(larcv.kProductImage2D, "wire" )
ev_masks = io.get_data(larcv.kProductClusterMask, "mrcnn_masks" )

mask_vv = ev_masks.as_vector()
print("Number of masks: ",[mask_vv.at(x).size() for x in range(3)])
hwire_v = [ larcv.as_th2d( ev_wire.as_vector().at(p), "hwire_p%d"%(p) ) for p in xrange(3) ]
meta_v  = [ ev_wire.as_vector().at(p).meta() for p in xrange(3) ]

#matchdata = ublarcvapp.dltagger.MaskMatchData( 0, 0, mask_vv.at(0).at(0) )
indices = std.vector("vector<int>")()
matchalgo.matchMasksAcrossPlanes( mask_vv, indices )

# visualize individual matches
ccombos = rt.TCanvas("ccombos","Combos", 1500, 500)
ccombos.Divide(3,1)

for icombo in xrange( matchalgo.m_combo_3plane_v.size() ):
    ccombos.Clear()
    ccombos.Divide(3,1)
    
    combo = matchalgo.m_combo_3plane_v.at(icombo)

    # get masks
    combo_indices = [ combo.indices.at(p) for p in xrange(3) ]
    print("indicies: {}".format(combo_indices))
    
    combo_masks = []
    for p,idx in enumerate(combo_indices):
        if idx<0:
            combo_masks.append( None )
        else:
            combo_masks.append( mask_vv.at(p).at( idx ) )
    print( "combo: {}".format(combo_masks) )

    # draw canvas and boxes
    tbox_v = []
    for p in xrange(3):
        ccombos.cd(1+p)
        if combo_masks[p] is None:
            continue
        hwire_v[p].Draw("colz")
        if p in [0,1]:
            hwire_v[p].GetXaxis().SetRangeUser(0,2400)
        print("bbox: ",combo_masks[p].box.min_x(), combo_masks[p].box.min_y(),combo_masks[p].box.max_x(), combo_masks[p].box.max_y() )
        bbox = rt.TBox( combo_masks[p].box.min_x(), meta_v[p].pos_y( int(combo_masks[p].box.min_y()) ),
                        combo_masks[p].box.max_x(), meta_v[p].pos_y( int(combo_masks[p].box.max_y()) ) )
        bbox.SetLineColor(rt.kRed)
        bbox.SetLineWidth(3)
        bbox.SetFillStyle(0)
        bbox.Draw("same")
        tbox_v.append(bbox)

    ccombos.Draw()
    ccombos.Update()
    ccombos.SaveAs("example_combos/combo_%02d.png"%(icombo))
    
    
    if combo.iou()<0.5:
        break

    
