#!/bin/env python3
from __future__ import print_function
import os,sys,argparse

parser = argparse.ArgumentParser("Test MCPixelPGraph")
parser.add_argument("-ill", "--input-larlite",required=True,type=str,help="Input larlite file")
#parser.add_argument("-adc", "--adc",type=str,default="wire",help="Name of tree with Wire ADC values [default: wire]")
#parser.add_argument("-tb",  "--tick-backward",action='store_true',default=False,help="Input LArCV data is tick-backward [default: false]")
parser.add_argument("-d",   "--debug", action='store_true', default=False, help="Run in debug mode")
parser.add_argument('-n', "--nentries", required=False, type=int, default=-1, help="number of entries to run")
parser.add_argument('-e', "--entry", required=False, type=int, default=-1, help="start at given entry number")
#parser.add_argument('-v', "--vis", required=False, default=False, action='store_true', help="if flag provided, will visualize event")
parser.add_argument("--cluster-nu",action='store_true',default=False, help="Cluster Nu particles together.")
parser.add_argument('-p','--pause',action='store_true',default=False, help="If provided, pause after each entry")
args = parser.parse_args()

import ROOT as rt
#from larcv import larcv
from larlite import larlite
from ublarcvapp import ublarcvapp

"""
test script that demos the MCPixelPGraph class.
"""

rt.gStyle.SetOptStat(0)

ioll = larlite.storage_manager( larlite.storage_manager.kREAD )
ioll.set_data_to_read( "mctrack",  "mcreco" ) # larlite.data.kMCTrack
ioll.set_data_to_read( "mcshower", "mcreco" ) # larlite.data.kMCShower
ioll.set_data_to_read( "mctruth", "generator" ) # larlite.data.kMCTruth
ioll.add_in_filename(  args.input_larlite )
ioll.open()

tot_nentries = ioll.get_entries()
start_entry = 0
print("Number of entries: ",tot_nentries)
if args.entry > 0:
    start_entry = args.entry
if args.nentries>0:
    nentries = args.nentries
else:
    nentries = tot_nentries
end_entry = start_entry + nentries
if end_entry>tot_nentries:
    end_entry = tot_nentries
    
print("Start loop.")

mcpg    = ublarcvapp.mctools.MCParticleGraph()
if args.debug:
    mcpg.set_verbosity( "debug" )
else:
    mcpg.set_verbosity( "info" )


for ientry in range( start_entry, end_entry ):

    print() 
    print("="*80)
    print("===[ EVENT ",ientry," ]===")
    ioll.go_to(ientry)

    mcpg.clear()
    if args.cluster_nu:
        mcpg.cluster_nu_particles( True )
    else:
        mcpg.cluster_nu_particles( False )

    mcpg.buildgraph( ioll )
    
    print("================================================",flush=True)
    print("PARSING PARTICLE GRAPH ONLY: No pixel matching  ",flush=True)
    print("================================================",flush=True)        
    #print("ALL NODE INFO [NO NU] --------------------------")
    #mcpg.printAllNodeInfo()
    #print(" -----------------------------------------------")
    print("CONSTRUCTED PARTICLE GRAPH",flush=True)
    mcpg.printGraph(0,False)
    #print("====================================================",flush=True)
    if args.pause:
        print("[ENTER] to continue",flush=True)
        input()

#print("=== FIN ==",flush=True)
