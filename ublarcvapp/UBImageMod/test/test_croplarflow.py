from __future__ import print_function
import ROOT as rt
from ROOT import std
from larcv import larcv
from ublarcvapp import ublarcvapp
import sys
sys.argv.append("-b")
larcv.load_rootutil()

superafile = sys.argv[1]
#superafile = "../../../../testdata/mcc9mar_bnbcorsika/larcv_mctruth_ee881c25-aeca-4c92-9622-4c21f492db41.root"

VIS_SPLITS = False
VIS_FLOWS  = True

io = larcv.IOManager(larcv.IOManager.kBOTH,"",larcv.IOManager.kTickBackward)
io.add_in_file( superafile )
io.set_out_file( "baka.root" )
io.initialize()

# -------------------------------------
# UBSplitDetector: setup for
#  -- whole-view splits
#  -- no cropping of values

split_cfg="""Verbosity:2
InputProducer: \"wiremc\"
OutputBBox2DProducer: \"detsplit\"
CropInModule: true
OutputCroppedProducer: \"detsplit\"
BBoxPixelHeight: 512
BBoxPixelWidth: 832
CoveredZWidth: 310
FillCroppedYImageCompletely: true
DebugImage: false
MaxImages: -1
RandomizeCrops: false
MaxRandomAttempts: 1000
MinFracPixelsInCrop: 0.0
"""
# create config file
print(split_cfg,file=open("ubsplit.cfg",'w'))
split_pset = larcv.CreatePSetFromFile( "ubsplit.cfg", "UBSplitDetector" )

# -------------------------------------
# UBLArFlowCropDetector

lfcrop_cfg="""Verbosity:0
InputBBoxProducer: \"detsplit\"
InputCroppedADCProducer: \"detsplit\"
InputADCProducer: \"wire\"
InputVisiProducer: \"pixvisi\"
InputFlowProducer: \"pixflow\"
OutputCroppedADCProducer:  \"adc\"
OutputCroppedVisiProducer: \"visi\"
OutputCroppedFlowProducer: \"flow\"
OutputCroppedMetaProducer: \"flowmeta\"
OutputFilename: \"baka_lf.root\"
SaveOutput: false
CheckFlow:  true
MakeCheckImage: true
DoMaxPool: false
RowDownsampleFactor: 2
ColDownsampleFactor: 2
MaxImages: -1
LimitOverlap: false
RequireMinGoodPixels: false
MaxOverlapFraction: 0.2
IsMC: true
UseVectorizedCode: true
"""

print(lfcrop_cfg,file=open("ublarflowcrop.cfg",'w'))
lfpset = larcv.CreatePSetFromFile( "ublarflowcrop.cfg", "UBCropLArFlow" )

# -------------------------------------
# ALGOS

split_algo = ublarcvapp.UBSplitDetector()
split_algo.configure(split_pset)
split_algo.initialize()

lfcrop_algo = ublarcvapp.UBCropLArFlow()
lfcrop_algo.configure(lfpset)
lfcrop_algo.initialize()

print("Configured ready to go. [ENTER] to continue.")
raw_input()

# -------------------------------------

nentries = io.get_n_entries()
print("Num Entries: ",nentries)
nentries = 1

thresholds_v = std.vector("float")(3,10.0)

if VIS_FLOWS:
    cflow = rt.TCanvas("cflow","cflow",1500,1200)
    cflow.Divide(3,3)
    cflow.Draw()

for n in range(nentries):
    io.read_entry(n)

    ev_adc = io.get_data(larcv.kProductImage2D,"wiremc")
    adc_v  = ev_adc.Image2DArray()

    ev_flow = io.get_data(larcv.kProductImage2D,"larflow")
    flow_v  = ev_flow.Image2DArray()

    ev_status = io.get_data(larcv.kProductChStatus,"wiremc")
    
    roi_v = std.vector("larcv::ROI")()
    out_v = std.vector("larcv::Image2D")()
    split_algo.process( adc_v, out_v, roi_v  )

    #detsplit = out.get_data( larcv.kProductImage2D, "detsplit" )
    print( "num rois: ",roi_v.size())
    print( "num crops: ",out_v.size())

    if VIS_SPLITS:
        # visualize
        h_v = {}
        for i in xrange(out_v.size()):
            h_v[i] = larcv.as_th2d( out_v.at(i), "test%d"%(i) )
            h_v[i].GetZaxis().SetRangeUser(-10,100)
        #print h_v

        c = rt.TCanvas("c","c",1500,400)
        c.Divide(3,1)
        for i in range(len(h_v)/3):
            for j in range(3):
                c.cd(j+1)
                h_v[3*i+j].Draw("COLZ")
                print( out_v.at(3*i+j).meta().dump())

            c.Update()
            raw_input()

    flow_vv = []
    for i in xrange(roi_v.size()):
        # get roi, pass to larflow cropper
        roi = roi_v.at(i)
        cropped_adc_v = std.vector("larcv::Image2D")()
        for j in xrange(3):
            cropped_adc_v.push_back( out_v.at(3*i+j) )
        cropped_flow_v = std.vector("larcv::Image2D")()
        cropped_visi_v = std.vector("larcv::Image2D")() # dummy for now
        print("Crop Flow ROI[%d] of %d"%(i,roi_v.size()))
        lfcrop_algo.make_cropped_flow_images( 2, roi, adc_v, ev_status, flow_v, thresholds_v, cropped_flow_v )

        hvis = std.vector("TH2D")()
        has_visi  = False
        vis_check = True
        lfcrop_algo.check_cropped_images( 2, cropped_adc_v, ev_status,
                                          thresholds_v, cropped_flow_v, cropped_visi_v, hvis, has_visi, vis_check )
        
        flow_vv.append( cropped_flow_v )
        if VIS_FLOWS:
            # first row: images
            h_v = {}
            canvspot = [1,2,0]
            htitle = ["TARGET1","TARGET2","SOURCE"]
            for j in range(3):
                cflow.cd(canvspot[j]+1)
                h_v[j] = larcv.as_th2d( out_v.at(3*i+j), "test%d_%d"%(i,j) )
                h_v[j].GetZaxis().SetRangeUser(-10,250)
                h_v[j].SetTitle("%s"%(htitle[j]))
                h_v[j].Draw("COLZ")                
            # first col as well for matrix
            cflow.cd(3+1)
            h_v[2].Draw("COLZ")
            cflow.cd(6+1)
            h_v[2].Draw("COLZ")

            # second row, the flows
            for j in xrange(0,2):
                cflow.cd(3+j+2)
                h_v[3+j] = larcv.as_th2d( cropped_flow_v.at(j), "flow%d_%d"%(i,j) )
                h_v[3+j].GetZaxis().SetRangeUser(-832,832)                
                h_v[3+j].Draw("COLZ")

            # third row, errors
            if vis_check:
                print("number of hists: %d"%(hvis.size()))
                for j in xrange(0,2):
                    cflow.cd(6+j+2)
                    hvis.at(2*j+1).Draw("COLZ")
                
            cflow.Update()
            print("Visualized adc and flow images")
            raw_input()
            
        
    for i in xrange(out_v.size()):
        detsplit.Append(out_v.at(i))
    out.set_id( io.event_id().run(), io.event_id().subrun(), io.event_id().event() )
    #print "save entry"
    #out.save_entry()
    break

algo.printElapsedTime()
io.finalize()
out.finalize()

print("==========================================")


print("FIN")
