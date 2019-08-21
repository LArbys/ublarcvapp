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

# collect graphical items
tbox_v = {0:[],1:[],2:[]}
tmarkers = {0:[],1:[],2:[]}
tcontours = {0:[],1:[],2:[]}
gpca_v = {0:[],1:[],2:[]}
tendpt_v = {0:[],1:[],2:[]}
tastar_v = {0:[],1:[],2:[]}

for icombo in xrange( matchalgo.m_combo_3plane_v.size() ):
    ccombos.Clear()
    ccombos.Divide(3,2)

    # get algo outputs
    combo     = matchalgo.m_combo_3plane_v.at(icombo)

    if icombo>=matchalgo.m_combo_crops_v.size():
        break
    
    combocrop = matchalgo.m_combo_crops_v.at(icombo)
    features  = matchalgo.m_combo_features_v.at(icombo)
    if icombo<matchalgo.m_combo_endpt3d_v.size():
        endpt3d   = matchalgo.m_combo_endpt3d_v.at(icombo)
    else:
        endpt3d   = None
    if icombo<matchalgo.m_combo_astar_v.size():
        astarout  = matchalgo.m_combo_astar_v.at(icombo)
    else:
        astarout  = None
    
    mask_contours = features.combo_mask_contour

    # get masks
    combo_indices = [ combo.indices.at(p) for p in xrange(3) ]
    print("COMBO[{}] indicies={}".format(icombo,combo_indices))

    # get clustermask objects
    combo_masks = []
    for p,idx in enumerate(combo_indices):
        if idx<0:
            combo_masks.append( None )
        else:
            combo_masks.append( mask_vv.at(p).at( idx ) )
    print( "combo: {}".format(combo_masks) )


    # draw canvas and boxes
    hcrop  = [ None if combo_indices[p]<0 else larcv.as_th2d( combocrop.crops_v.at(p), "hcrop_combo%d_p%d"%(icombo,p) ) for p in xrange(3) ]
    is3plane = True    
    for ip,h in enumerate(hcrop):
        if h is None:
            is3plane = False
            hcrop[ip] = larcv.as_th2d( combocrop.missing_v[ip], "hmissingcrop_combo%d_p%d"%(icombo,p))

    castar = rt.TCanvas("caster","ASTAR", 1500, 500 )
    castar.Divide(3,1)
    hastar_v = []
    for p in xrange(3):
        ccombos.cd(1+p)

        hwire_v[p].Draw("colz")
        if p in [0,1]:
            hwire_v[p].GetXaxis().SetRangeUser(0,2400)
        hwire_v[p].SetTitle("Plane %d: IOU %.3f;wire;tick"%(p,combo.iou()))
        if combo_masks[p] is not None:
            print("bbox: ",combo_masks[p].box.min_x(), combo_masks[p].box.min_y(),combo_masks[p].box.max_x(), combo_masks[p].box.max_y() )
            bbox = rt.TBox( combo_masks[p].box.min_x(), meta_v[p].pos_y( int(combo_masks[p].box.min_y()) ),
                            combo_masks[p].box.max_x(), meta_v[p].pos_y( int(combo_masks[p].box.max_y()) ) )
            if is3plane:
                bbox.SetLineColor(rt.kRed)
            else:
                bbox.SetLineColor(rt.kBlue)
                
            bbox.SetLineWidth(3)
            bbox.SetFillStyle(0)
            bbox.Draw("same")
            tbox_v[p].append(bbox)

        # crop
        ccombos.cd(4+p)
        hcrop[p].Draw("colz")

        # get the meta for the crop
        if combo_masks[p] is not None:
            mask_meta = combocrop.mask_v.at(p).meta()
        else:
            mask_meta = combocrop.missing_v.at(p).meta()

        # draw contours, pca, end points and stuff

        # markers and contours
        if combo_masks[p] is not None:
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
                tcontours[p].append(gcontour)
                
            gmarker.SetMarkerStyle(20)
            gmarker.SetMarkerColor(rt.kRed)
            gmarker.Draw("P")
            tmarkers[p].append(gmarker)

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
            # backward dir
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
            #print("plane[%d] pca mean: "%(p),pca_mean)
            #print("plane[%d] pca ends: min="%(p),ptmin," max=",ptmax)
            #print("plane[%d] pca1-dir: "%(p),pca1_dir)        
            gpca.SetPoint(0,ptmin[0],ptmin[1])
            gpca.SetPoint(1,pca_mean[0],pca_mean[1])
            gpca.SetPoint(2,ptmax[0],ptmax[1])
            gpca.SetMarkerStyle(24)
            gpca.SetMarkerColor(rt.kMagenta)
            #gpca.SetMarkerSize(2)
            gpca.SetLineColor(rt.kMagenta)
            gpca.SetLineWidth(2)
            gpca.Draw("LP")
            gpca_v[p].append(gpca)


        # 3d points (projected of course)
        if endpt3d is not None:
            tendpt = rt.TGraph(2)
            #print("plane[%d] endpoint: "%(p), endpt3d.endpt_wid_v[0][p], ",", endpt3d.endpt_tyz_v[0][0] )
            #print("plane[%d] endpoint: "%(p), endpt3d.endpt_wid_v[1][p], ",", endpt3d.endpt_tyz_v[1][0] )
            tendpt.SetPoint(0, endpt3d.endpt_wid_v[0][p], endpt3d.endpt_tyz_v[0][0] )
            tendpt.SetPoint(1, endpt3d.endpt_wid_v[1][p], endpt3d.endpt_tyz_v[1][0] )
            tendpt.SetMarkerStyle(23)
            tendpt.SetMarkerSize(2)
            tendpt.SetMarkerColor( rt.kGreen )
            tendpt.Draw("P")
            tendpt_v[p].append(tendpt)
            hcrop[p].SetTitle("Plane %d: tri=(%.3f,%.3f) x=(%d,%d);wire;tick"%(p,endpt3d.endpt_tri_v[0],endpt3d.endpt_tri_v[1],endpt3d.endpt_tpc_v[0],endpt3d.endpt_tpc_v[1]))

        
        # astar path
        if astarout is not None:
            astar_path = astarout.astar_path
            tastar = rt.TGraph( astar_path.size() )
            tastar.SetMarkerStyle(24)            
            for ipt in xrange( astar_path.size() ):
                astar_node = astar_path[ipt]
                node_wire = mask_meta.pos_x( int(astar_node.cols[p]) )
                #print("Astar node[%d] tick=%d wire=%d"%(p,astar_node.tyz[0],node_wire))
                tastar.SetPoint( ipt, node_wire, astar_node.tyz[0] )

            if astarout.astar_completed==1:
                tastar.SetMarkerColor( rt.kBlue )
                tastar.SetLineColor( rt.kBlue )
            else:
                tastar.SetMarkerColor( rt.kBlack )
                tastar.SetLineColor( rt.kBlack )
            tastar.Draw("LP")
            tastar_v[p].append( tastar )

            # astar canvas
            castar.cd(1+p)
            if astarout.score_crop_v.size()>0:
                hastar = larcv.as_th2d( astarout.score_crop_v.at( p ), "hastar_combo%d_p%d"%(icombo,p) )
                hastar.Draw("colz")
                tastar.Draw("LP")
                hastar_v.append(hastar)
            castar.Update()
    
    ccombos.Draw()
    ccombos.Update()
    ccombos.SaveAs("example_combos/combo_%02d.png"%(icombo))

    castar.Draw()
    castar.Update()
    castar.SaveAs("example_combos/astar_%02d.png"%(icombo))

call = rt.TCanvas("call", "All", 1200,1500 )
call.Divide(1,3)
for p in xrange(3):
    call.cd(1+p)
    hwire_v[p].Draw("colz")
    for g in tastar_v[p]:
        g.Draw("LP")
    for g in gpca_v[p]:
        g.Draw("LP")
    for g in tbox_v[p]:
        g.Draw()
call.Update()
call.Draw()
call.SaveAs("example_combos/all_combos_pass1.png")

raw_input()
