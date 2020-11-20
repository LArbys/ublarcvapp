import os,sys,argparse

parser = argparse.ArgumentParser("Test MCPixelPGraph")
parser.add_argument("-ill", "--input-larlite",required=True,type=str,help="Input larlite file")
parser.add_argument("-ilcv","--input-larcv",required=True,type=str,help="Input LArCV file")
parser.add_argument("-o","--output",required=True,type=str,help="Output LArCV file")
parser.add_argument("-adc", "--adc",type=str,default="wire",help="Name of tree with Wire ADC values [default: wire]")
parser.add_argument("-tb",  "--tick-backward",action='store_true',default=False,help="Input LArCV data is tick-backward [default: false]")
args = parser.parse_args()

import ROOT as rt
from larcv import larcv
from larlite import larlite
from ublarcvapp import ublarcvapp

"""
test script that demos the MCPixelPGraph class.
"""

rt.gStyle.SetOptStat(0)

ioll = larlite.storage_manager( larlite.storage_manager.kREAD )
ioll.add_in_filename(  args.input_larlite )
ioll.set_data_to_read( larlite.data.kMCTrack,  "mcreco" )
ioll.set_data_to_read( larlite.data.kMCShower, "mcreco" )
ioll.open()

if args.tick_backward:
    iolcv = larcv.IOManager( larcv.IOManager.kBOTH, "larcv", larcv.IOManager.kTickBackward )
else:
    iolcv = larcv.IOManager( larcv.IOManager.kBOTH, "larcv", larcv.IOManager.kTickForward )
iolcv.add_in_file( args.input_larcv )
iolcv.specify_data_read( larcv.kProductImage2D, "wire" );
iolcv.specify_data_read( larcv.kProductImage2D, "wiremc" );
iolcv.specify_data_read( larcv.kProductChStatus, "wire" );
iolcv.specify_data_read( larcv.kProductChStatus, "wiremc" );
iolcv.specify_data_read( larcv.kProductImage2D, "thrumu" );
iolcv.specify_data_read( larcv.kProductImage2D, "ancestor" );
iolcv.specify_data_read( larcv.kProductImage2D, "segment" );
iolcv.specify_data_read( larcv.kProductImage2D, "instance" );
iolcv.specify_data_read( larcv.kProductImage2D, "larflow" );
iolcv.reverse_all_products()
iolcv.addto_storeonly_list( larcv.kProductImage2D, "nupix" )
iolcv.set_out_file( args.output )
iolcv.initialize()

nentries = iolcv.get_n_entries()
print "Number of entries: ",nentries
#nentries = 10

print "Start loop."

nupix = ublarcvapp.mctools.NeutrinoPixelFilter()

tmp = rt.TFile("temp.root","recreate")

c = rt.TCanvas("c","c",1200,1800)
c.Divide(1,3)

for ientry in xrange( nentries ):

    print 
    print "=========================="
    print "===[ EVENT ",ientry," ]==="
    ioll.go_to(ientry)
    iolcv.read_entry(ientry)

    nupix.process(iolcv,ioll)
    iolcv.save_entry()
    print "[enter to continue]"
    raw_input()
    break

ioll.close()
iolcv.finalize()
print "=== FIN =="
