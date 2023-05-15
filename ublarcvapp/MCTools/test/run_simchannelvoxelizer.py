from __future__ import print_function
import os,sys,argparse

parser = argparse.ArgumentParser("Run SimChannelVoxelizer")
parser.add_argument("-ill", "--input-larlite",required=True,type=str,help="Input larlite file")
parser.add_argument("-t","--tree-name",type=str,default="driftWC simpleSC",help="Set name of tree with SimCh objects")
parser.add_argument("-d","--detector",type=str,default="uboone",help="Set detector name: uboone, sbnd, icarus. default: uboone")
parser.add_argument("--debug", action='store_true', default=False, help="Run in debug mode")
args = parser.parse_args()

import ROOT as rt
from ROOT import std
from larlite import larlite
from larlite import larutil
from ublarcvapp import ublarcvapp

"""
test script that runs the SimChannelVoxelizer

This class takes simulation truth metadata from the SimChannel object which stores individual energy deposits
assigned to each wire and stores it into a voxel grid as various labels.
"""

if args.detector=="uboone":
    larutil.LArUtilConfig.SetDetector( larlite.geo.kMicroBooNE )
elif args.detector=="sbnd":    
    larutil.LArUtilConfig.SetDetector( larlite.geo.kSBND )
else:
    raise ValueError("unrecognized detector name: ",args.detector)


rt.gStyle.SetOptStat(0)

ioll = larlite.storage_manager( larlite.storage_manager.kREAD )
ioll.add_in_filename(  args.input_larlite )
ioll.open()

nentries = ioll.get_entries()
print("Number of entries: ",nentries)
#nentries = 1 # to override for debugging

tmp = rt.TFile("testout_simchannelvoxelizer.root","recreate")
out = rt.TTree("simchannelvoxelizer","Output of simchannel voxelizer")
index_v  = std.vector("larcv::NumpyArrayInt")()
trackid_v  = std.vector("larcv::NumpyArrayInt")()
ancestor_v  = std.vector("larcv::NumpyArrayInt")()
pdg_v  = std.vector("larcv::NumpyArrayInt")()
charge_v = std.vector("larcv::NumpyArrayFloat")()
truepos_v = std.vector("larcv::NumpyArrayFloat")()
simchimg_v = std.vector("larcv::NumpyArrayFloat")()
out.Branch( "coordindex_v", index_v )
out.Branch( "index_v", index_v )
out.Branch( "trackid_v", trackid_v )
out.Branch( "ancestor_v", ancestor_v )
out.Branch( "pdg_v", pdg_v )
out.Branch( "charge_v", charge_v )
out.Branch( "truepos_v", truepos_v )
out.Branch( "simchimg_v", simchimg_v )


print("Create the sim channel voxelizer.")
#input()

simchvoxer = ublarcvapp.mctools.SimChannelVoxelizer()
simchvoxer.set_simch_treename( args.tree_name )

print("run the event loop")
#input()

for ientry in range( nentries ):

    charge_v.clear()
    index_v.clear()
    ancestor_v.clear()
    trackid_v.clear()
    pdg_v.clear()
    truepos_v.clear()
    simchimg_v.clear()
    
    print() 
    print("==========================")
    print("===[ EVENT ",ientry," ]===")
    ioll.go_to(ientry)

    simchvoxer.process( ioll )

    for i in range( simchvoxer._tpcdata_v.size() ):
        tpcinfo = simchvoxer._tpcdata_v.at(i);
        print("tpc origin: ",tpcinfo._origin_cm_v[0]," ",tpcinfo._origin_cm_v[1]," ",tpcinfo._origin_cm_v[2])
        charge_v.push_back( tpcinfo._charge_v )
        index_v.push_back( tpcinfo._coordindex_v )
        trackid_v.push_back( tpcinfo._trackid_v )
        ancestor_v.push_back( tpcinfo._ancestorid_v )
        pdg_v.push_back( tpcinfo._pdg_v )
        truepos_v.push_back( tpcinfo._coordpos_v )
        for ii in range( tpcinfo._simch_img_v.size() ):
            simchimg_v.push_back( tpcinfo._simch_img_v.at(ii) )
    
    out.Fill()
    print("[enter to continue]")
    #input()
    if ientry==0:
        break

tmp.cd()
out.Write()
tmp.Close()
print("=== FIN ==")
