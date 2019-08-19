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
matchalgo.matchMasksAcrossPlanes( mask_vv, ev_wire.Image2DArray(), indices )

# visualize individual matches
ccombos = rt.TCanvas("ccombos","Combos", 1500, 1200)
ccombos.Divide(3,2)

for icombo in xrange( matchalgo.m_combo_3plane_v.size() ):
    ccombos.Clear()
    ccombos.Divide(3,2)
    
    combo = matchalgo.m_combo_3plane_v.at(icombo)
    mask_contours = matchalgo.m_combo_mask_contour_v.at(icombo)

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

    # crops
    combocrop = ublarcvapp.dltagger.CropMaskCombo( combo, ev_wire.Image2DArray() )

    # draw canvas and boxes
    tbox_v = []
    hcrop  = [ larcv.as_th2d( combocrop.crops_v.at(p), "hcrop_combo%d_p%d"%(idx,p) ) for p in xrange(3) ]
    tmarkers = []
    tcontours = []
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

        # crop
        ccombos.cd(4+p)
        hcrop[p].Draw("colz")

        # draw contours, pca, end points and stuff
        mask_meta = combocrop.mask_v.at(p).meta()
        ncontours = int(mask_contours.m_plane_atomicmeta_v.at(p).size())
        gmarker = rt.TGraph( ncontours*2 )
        for ictr in xrange( ncontours ):
            contour = mask_contours.m_plane_atomicmeta_v.at(p).at(ictr)
            # contour end points
            startpt = contour.getFitSegmentStart()
            endpt   = contour.getFitSegmentEnd()
            gmarker.SetPoint( 2*ictr,   mask_meta.pos_x(startpt.x), mask_meta.pos_y(startpt.y) )
            gmarker.SetPoint( 2*ictr+1, mask_meta.pos_x(endpt.x),   mask_meta.pos_y(endpt.y) )
            # contour graph
            gcontour = rt.TGraph( contour.size()+1 )
            for ipt in xrange(contour.size()):
                gcontour.SetPoint( ipt, mask_meta.pos_x(contour.at(ipt).x), mask_meta.pos_y(contour.at(ipt).y) )
                if ipt==0:
                    gcontour.SetPoint( contour.size(), mask_meta.pos_x(contour.at(ipt).x), mask_meta.pos_y(contour.at(ipt).y) )
            gcontour.SetLineColor(rt.kBlack)
            gcontour.SetLineWidth(1)
            gcontour.Draw("L")
            tcontours.append(gcontour)
        gmarker.SetMarkerStyle(20)
        gmarker.SetMarkerColor(rt.kRed)
        gmarker.Draw("P")
        tmarkers.append(gmarker)
        
    ccombos.Draw()
    ccombos.Update()
    ccombos.SaveAs("example_combos/combo_%02d.png"%(icombo))
    if True:
        break
    
    
    if combo.iou()<0.5:
        break

raw_input()
