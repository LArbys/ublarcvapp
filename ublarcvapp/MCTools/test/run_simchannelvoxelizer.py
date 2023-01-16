from __future__ import print_function
import os,sys,argparse

parser = argparse.ArgumentParser("Run SimChannelVoxelizer")
parser.add_argument("-ill", "--input-larlite",required=True,type=str,help="Input larlite file")
parser.add_argument("-d",   "--debug", action='store_true', default=False, help="Run in debug mode")
args = parser.parse_args()

import ROOT as rt
from ROOT import std
from larlite import larlite
from larlite import larutil
from ublarcvapp import ublarcvapp

"""
test script that demos the MCPixelPGraph class.
"""

larutil.LArUtilConfig.SetDetector( larlite.geo.kSBND )

rt.gStyle.SetOptStat(0)

ioll = larlite.storage_manager( larlite.storage_manager.kREAD )
ioll.add_in_filename(  args.input_larlite )
ioll.open()

nentries = ioll.get_entries()
print("Number of entries: ",nentries)
nentries = 1 # to override for debugging

tmp = rt.TFile("testout_simchannelvoxelizer.root","recreate")
out = rt.TTree("simchannelvoxelizer","Output of simchannel voxelizer")
index_v  = std.vector("larcv::NumpyArrayInt")()
trackid_v  = std.vector("larcv::NumpyArrayInt")()
ancestor_v  = std.vector("larcv::NumpyArrayInt")()
pdg_v  = std.vector("larcv::NumpyArrayInt")()
charge_v = std.vector("larcv::NumpyArrayFloat")()
out.Branch( "coordindex_v", index_v )
out.Branch( "index_v", index_v )
out.Branch( "trackid_v", trackid_v )
out.Branch( "ancestor_v", ancestor_v )
out.Branch( "pdg_v", pdg_v )
out.Branch( "charge_v", charge_v )


print("Create the sim channel voxelizer.")
input()

simchvoxer = ublarcvapp.mctools.SimChannelVoxelizer()

print("run the event loop")
input()

for ientry in range( nentries ):

    charge_v.clear()
    index_v.clear()
    ancestor_v.clear()
    trackid_v.clear()
    pdg_v.clear()
    
    print() 
    print("==========================")
    print("===[ EVENT ",ientry," ]===")
    ioll.go_to(ientry)

    #ev_simch = ioll.get_data( larlite.data.kSimChannel, "simdrift" )
    #print("Number of elements in event simdrift: ",ev_simch.size())    
    #simchvoxer.process( ev_simch )
    simchvoxer.process( ioll )

    for i in range( simchvoxer._tpcdata_v.size() ):
        tpcinfo = simchvoxer._tpcdata_v.at(i);
        charge_v.push_back( tpcinfo._charge_v )
        index_v.push_back( tpcinfo._coordindex_v )
        trackid_v.push_back( tpcinfo._trackid_v )
        ancestor_v.push_back( tpcinfo._ancestorid_v )
        pdg_v.push_back( tpcinfo._pdg_v )        
    
    out.Fill()
    print("[enter to continue]")
    input()    

tmp.cd()
out.Write()
tmp.Close()
print("=== FIN ==")
