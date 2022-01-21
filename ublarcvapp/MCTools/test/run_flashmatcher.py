import os,sys,argparse

parser = argparse.ArgumentParser("Test FlashMatcher")
parser.add_argument("-ill", "--input-larlite",required=True,type=str,help="Input larlite mcinfo file")
# want another argument for opflash file
parser.add_argument("-opreco", "--input-opreco",required=True,type=str,help="Input larlite opreco file")
#parser.add_argument("-iv", "--input-voxelfile",required=True,type=str,help="Input voxeldata file")
parser.add_argument("-out", "--output_file",required=True,type=str,help="Output file")
#parser.add_argument("-ilcv","--input-larcv",required=True,type=str,help="Input LArCV file")
#parser.add_argument("-adc", "--adc",type=str,default="wire",help="Name of tree with Wire ADC values [default: wire]")
#parser.add_argument("-tb",  "--tick-backward",action='store_true',default=False,help="Input LArCV data is tick-backward [default: false]")

args = parser.parse_args()

import ROOT as rt
from larlite import larlite
#from larcv import larcv
from ublarcvapp import ublarcvapp

rt.gROOT.ProcessLine( "gErrorIgnoreLevel = 3002;" )

"""
test script that demos the Flash Matcher class.
"""

rt.gStyle.SetOptStat(0)

ioll = larlite.storage_manager( larlite.storage_manager.kREAD )
ioll.add_in_filename(  args.input_larlite )
ioll.open()

opio = larlite.storage_manager( larlite.storage_manager.kREAD )
opio.add_in_filename(  args.input_opreco )
opio.open()

f = TFile(args.input_voxelfile,"READ")
print("passed tfile part")

#voxio = larlite.storage_manager( larlite.storage_manager.kREAD )
#voxio.add_in_filename(  args.input_voxelfile )
#voxio.open()

outio = larlite.storage_manager( larlite.storage_manager.kWRITE )
outio.set_out_filename(  args.output_file )
outio.open()

nentries = ioll.get_entries()
print("Number of entries: ",nentries)

print("Start loop.")
fmutil = ublarcvapp.mctools.FlashMatcher()
#vtxutil = ublarcvapp.mctools.NeutrinoVertex
#print(vtxutil)
print(fmutil)

#c = rt.TCanvas("c","c",1200,1800)
#c.Divide(1,3)

listy = []

# for testing:
#match = fmutil.matchTicks( 100.0, [], 1 )
#print(match)

print("isCosmic from ctor: ",fmutil.isCosmic)

fmutil.initialize( voxio )

#for ientry in range( nentries ):
for ientry in range( 3 ):

    print()
    print("==========================")
    print("===[ EVENT ",ientry," ]===")
    ioll.go_to(ientry)
    opio.go_to(ientry)

    numTracks = fmutil.numTracks( ioll )

    for i in range( 0, numTracks, 1 ):
        track_tick = fmutil.grabTickFromMCTrack( ioll, i )

        producer = fmutil.producer
        isCosmic = fmutil.isCosmic

        print("Now isCosmic has been set to: ",fmutil.isCosmic)
        print(producer)

        op_tick = fmutil.grabTickFromOpflash( opio )
        #    fmtrack_tick = vtxutil.getImageCoords( ioll )
        print("track_tick",track_tick)
        print("op_tick",op_tick)
        #    iolcv.read_entry(ientry)
        match = fmutil.matchTicks( track_tick, op_tick )
        print("Found match: ", match )

        if match != -999.999 and match != 999.999 and track_tick != -999.997:
            fmutil.process( ioll )

        if isCosmic == 1 and match != -999.999 and track_tick != -999.997:
            listy.append( match - track_tick )

print("list: ", listy)

fmutil.finalize()

print("size of list: ", len(listy))

print("=== FIN ==")
