import os,sys,argparse

parser = argparse.ArgumentParser("Test FlashMatcher")
parser.add_argument("-ill", "--input-larlite",required=True,type=str,help="Input larlite file")
# want another argument for opflash file
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

nentries = ioll.get_entries()
print("Number of entries: ",nentries)

print("Start loop.")
fmutil = ublarcvapp.mctools.FlashMatcher
print(fmutil)

#c = rt.TCanvas("c","c",1200,1800)
#c.Divide(1,3)

for ientry in range( nentries ):

    print() 
    print("==========================")
    print("===[ EVENT ",ientry," ]===")
    ioll.go_to(ientry)
#    iolcv.read_entry(ientry)

 
print("=== FIN ==")
