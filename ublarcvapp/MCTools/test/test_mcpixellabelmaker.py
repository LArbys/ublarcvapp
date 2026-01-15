import os,sys
import argparse

parser = argparse.ArgumentParser("Test MCPixelLabelMaker")
parser.add_argument("-i",   "--input-dlmerged",required=True,type=str,help="Input larlite file")
parser.add_argument('-e',   "--entry", required=True, type=int, default=0, help="entry to process")
parser.add_argument("-adc", "--adc", required=False, type=str,default="wiremc",help="Name of tree with Wire ADC values [default: wiremc]")
parser.add_argument("-tb",  "--tick-backward",action='store_true',default=False,help="Input LArCV data is tick-backward [default: false]")
parser.add_argument("-d",   "--debug", action='store_true', default=False, help="Run in debug mode")
#parser.add_argument('-n', "--nentries", required=False, type=int, default=-1, help="number of entries to run")
parser.add_argument('-o',   "--output-hdf5", required=False, default="", help="if provided, will save data to HDF5 file")
args = parser.parse_args()


from larlite import larlite
from larcv import larcv
from ublarcvapp import ublarcvapp

#input_dlmerged = "/mnt/ddrive/data/ub_on_tufts/corsika_bnb_nu_pi0/dlmerged_coriska_bnb_nu_pi0_fileno000001.root"
input_dlmerged = args.input_dlmerged

ENTRY=args.entry

mcpixel_label_maker = ublarcvapp.mctools.MCPixelLabelMaker()
mcpixel_label_maker.set_verbosity(1)

ioll = larlite.storage_manager( larlite.storage_manager.kREAD )
ioll.add_in_filename( input_dlmerged )
ioll.set_verbosity(2)
ioll.open()

iolcv = larcv.IOManager( larcv.IOManager.kREAD, "larcv" )
iolcv.add_in_file( input_dlmerged )
iolcv.set_verbosity(2)
iolcv.initialize()

ioll.go_to(ENTRY)
iolcv.read_entry(ENTRY)

mcpixel_label_maker.process( ioll, iolcv, args.adc )

if args.output_hdf5!="":
    mcpixel_label_maker.export_as_hdf(args.output_hdf5)
