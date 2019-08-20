from __future__ import print_function

import os,sys
from math import fabs

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
ev_chstatus = io.get_data(larcv.kProductChStatus, "wire" )

mask_vv = ev_masks.as_vector()
print("Number of masks: ",[mask_vv.at(x).size() for x in range(3)])
hwire_v = [ larcv.as_th2d( ev_wire.as_vector().at(p), "hwire_p%d"%(p) ) for p in xrange(3) ]
meta_v  = [ ev_wire.as_vector().at(p).meta() for p in xrange(3) ]

indices = std.vector("vector<int>")()
matchalgo.matchMasksAcrossPlanes( mask_vv, ev_wire.Image2DArray(), ev_chstatus, indices )

# visualize individual matches
ccombos = rt.TCanvas("ccombos","Combos", 1500, 1200)
ccombos.Divide(3,2)

for icombo in xrange( matchalgo.m_combo_3plane_v.size() ):
    ccombos.Clear()
    ccombos.Divide(3,2)

    # get algo outputs
    combo     = matchalgo.m_combo_3plane_v.at(icombo)

    if icombo>=matchalgo.m_combo_crops_v.size():
        break
    
    combocrop = matchalgo.m_combo_crops_v.at(icombo)
    features  = matchalgo.m_combo_features_v.at(icombo)
    endpt3d   = matchalgo.m_combo_endpt3d_v.at(icombo)
    astarout  = matchalgo.m_combo_astar_v.at(icombo)
    
    mask_contours = features.combo_mask_contour

    # get masks
    combo_indices = [ combo.indices.at(p) for p in xrange(3) ]
    print("indicies: {}".format(combo_indices))

    # get clustermask objects
    combo_masks = []
    for p,idx in enumerate(combo_indices):
        if idx<0:
            combo_masks.append( None )
        else:
            combo_masks.append( mask_vv.at(p).at( idx ) )
    print( "combo: {}".format(combo_masks) )


    # draw canvas and boxes
    tbox_v = []
    hcrop  = [ larcv.as_th2d( combocrop.crops_v.at(p), "hcrop_combo%d_p%d"%(idx,p) ) for p in xrange(3) ]
    tmarkers = []
    tcontours = []
    gpca_v = []
    tendpt_v = []
    tastar_v = []
    for p in xrange(3):
        ccombos.cd(1+p)
        if combo_masks[p] is None:
            continue
        hwire_v[p].Draw("colz")
        if p in [0,1]:
            hwire_v[p].GetXaxis().SetRangeUser(0,2400)
        hwire_v[p].SetTitle("Plane %d: IOU %.3f;wire;tick"%(p,combo.iou()))
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

        # markers and contours
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

        # pca lines
        pca_mean = [ mask_meta.pos_x( int(features.pca_mean_vv.at(p).at(0)) ),
                     mask_meta.pos_y( int(features.pca_mean_vv.at(p).at(1)) ) ]
        pca1_dir = [ features.pca1_dir_vv.at(p).at(x) for x in xrange(2) ]
        pca1_dir[1] *= 6.0 # change y-axis scale from row pixels to ticks

        # forward dir
        if pca1_dir[0]>0:
            dw = (mask_meta.max_x()-pca_mean[0])/pca1_dir[0]
        else:
            dw = (mask_meta.min_x()-pca_mean[0])/pca1_dir[0]
        if pca1_dir[1]>0:
            dt = (mask_meta.max_y()-pca_mean[1])/pca1_dir[1]
        else:
            dt = (mask_meta.min_y()-pca_mean[1])/pca1_dir[1]
        
        if fabs(dt)<fabs(dw):
            dd = dt
        else:
            dd = dw
        ptmax = [ pca_mean[x] + dd*pca1_dir[x] for x in xrange(2) ]

        # backwards
        if pca1_dir[0]<0:
            dw = (mask_meta.max_x()-pca_mean[0])/pca1_dir[0]
        else:
            dw = (mask_meta.min_x()-pca_mean[0])/pca1_dir[0]
        if pca1_dir[1]<0:
            dt = (mask_meta.max_y()-pca_mean[1])/pca1_dir[1]
        else:
            dt = (mask_meta.min_y()-pca_mean[1])/pca1_dir[1]
        if fabs(dt)<fabs(dw):
            dd = dt
        else:
            dd = dw
        ptmin = [ pca_mean[x] + dd*pca1_dir[x] for x in xrange(2) ]
        
        gpca = rt.TGraph(3)
        print("plane[%d] pca mean: ",pca_mean)
        print("plane[%d] pca ends: min="%(p),ptmin," max=",ptmax)
        print("plane[%d] pca1-dir: "%(p),pca1_dir)        
        gpca.SetPoint(0,ptmin[0],ptmin[1])
        gpca.SetPoint(1,pca_mean[0],pca_mean[1])
        gpca.SetPoint(2,ptmax[0],ptmax[1])
        gpca.SetMarkerStyle(24)
        gpca.SetMarkerColor(rt.kMagenta)
        #gpca.SetMarkerSize(2)
        gpca.SetLineColor(rt.kMagenta)
        gpca.SetLineWidth(2)
        gpca.Draw("LP")
        gpca_v.append(gpca)


        # 3d points (projected of course)
        tendpt = rt.TGraph(2)
        print("plane[%d] endpoint: "%(p), endpt3d.endpt_wid_v[0][p], ",", endpt3d.endpt_tyz_v[0][0] )
        print("plane[%d] endpoint: "%(p), endpt3d.endpt_wid_v[1][p], ",", endpt3d.endpt_tyz_v[1][0] )
        tendpt.SetPoint(0, endpt3d.endpt_wid_v[0][p], endpt3d.endpt_tyz_v[0][0] )
        tendpt.SetPoint(1, endpt3d.endpt_wid_v[1][p], endpt3d.endpt_tyz_v[1][0] )
        tendpt.SetMarkerStyle(23)
        tendpt.SetMarkerSize(2)
        tendpt.SetMarkerColor( rt.kGreen )
        tendpt.Draw("P")
        tendpt_v.append(tendpt)
        hcrop[p].SetTitle("Plane %d: tri=(%.3f,%.3f) x=(%d,%d);wire;tick"%(p,endpt3d.endpt_tri_v[0],endpt3d.endpt_tri_v[1],endpt3d.endpt_tpc_v[0],endpt3d.endpt_tpc_v[1]))

        
        # astar path
        astar_path = astarout.astar_path
        tastar = rt.TGraph( astar_path.size() )
        for ipt in xrange( astar_path.size() ):
            astar_node = astar_path[ipt]
            node_wire = mask_meta.pos_x( int(astar_node.cols[p]) )
            print("Astar node[%d] tick=%d wire=%d"%(p,astar_node.tyz[0],node_wire))
            tastar.SetPoint( ipt, node_wire, astar_node.tyz[0] )
        tastar.SetMarkerStyle(24)
        if astarout.astar_completed==1:
            tastar.SetMarkerColor( rt.kBlue )
            tastar.SetLineColor( rt.kBlue )
        else:
            tastar.SetMarkerColor( rt.kBlack )
            tastar.SetLineColor( rt.kBlack )
        tastar.Draw("LP")
        tastar_v.append( tastar )
        
    ccombos.Draw()
    ccombos.Update()
    ccombos.SaveAs("example_combos/combo_%02d.png"%(icombo))
    
    if combo.iou()<0.5:
        break

raw_input()
