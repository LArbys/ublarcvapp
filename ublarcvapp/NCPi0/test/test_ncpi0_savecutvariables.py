from __future__ import print_function

import os,sys
from math import fabs
import argparse

import ROOT as rt
from ROOT import std
from larlite import larlite
from larcv import larcv
from ublarcvapp import ublarcvapp

def parse_args():
    # parse input arguments
    parser = argparse.ArgumentParser(description='GetNCPI0LikelihoodProbs')
    parser.add_argument('--larcvtruth')
    parser.add_argument('--ssnet')
    parser.add_argument('--mcinfo')
    parser.add_argument('--tracker')
    parser.add_argument('--showerreco_ana')
    parser.add_argument('--showerreco_larlite')

    return parser.parse_args();


def main():
    args = parse_args()
    print(args)
    larcv.load_rootutil()

    larcvinputfiles = [args.larcvtruth] # backward example
    larcvforwardinputfiles = [args.ssnet]
    larliteinputfiles = [args.mcinfo,args.tracker,args.showerreco_larlite]


    # #larcv io manager
    # io_cfg = """
    # Verbosity:1
    # Name: \"IOManager\"
    # IOMode: 0
    # ReadOnlyName: [\"wiremc\",\"segment\",\"segment\",\"instance\",\"uburn_plane0\",\"uburn_plane1\",\"uburn_plane2\"]
    # ReadOnlyType: [0,0,1,0,0,0,0]
    # TickBackward: true
    # ReverseImage2D: [\"wiremc\",\"segment\",\"segment\",\"instance\",\"uburn_plane0\",\"uburn_plane1\",\"uburn_plane2\"]
    # """
    # iocfg = open("larcviomanager.cfg",'w')
    # print(io_cfg, file=iocfg)
    # iocfg.close()
    # iopset = larcv.CreatePSetFromFile("larcviomanager.cfg", "LarcvIOManager")
    # io = larcv.IOManager(larcv.CreatePSetFromFile("larcviomanager.cfg", "LarcvIOManager"))
    io = larcv.IOManager(larcv.IOManager.kREAD, "IO", larcv.IOManager.kTickBackward )
    for inputfile in larcvinputfiles:
        io.add_in_file( inputfile )
    io.initialize()

    ioforward = larcv.IOManager(larcv.IOManager.kREAD, "IO", larcv.IOManager.kTickForward)
    for inputfile in larcvforwardinputfiles:
        ioforward.add_in_file( inputfile )
    ioforward.initialize()

    savecutvariables = ublarcvapp.ncpi0.SaveCutVariables()

    # larlite io manager
    ioll = larlite.storage_manager(larlite.storage_manager.kREAD)
    for l in larliteinputfiles:
        ioll.add_in_filename( l )
    ioll.open()

    # # input anafiles
    # showerrecoanafile = new TFile(args.showerreco_ana,"read")

    # ----make config file----
    savevariables_cfg = """
    InputADCProducer: \"wiremc\"
    InputPartROIProducer: \"segment\"
    InputSegmentProducer: \"segment\"
    InputInstanceProducer: \"instance\"
    InputSSNetUProducer: \"uburn_plane0\"
    InputSSNetVProducer: \"uburn_plane1\"
    InputSSNetYProducer: \"uburn_plane2\"
    InputMCTrackProducer: \"mcreco\"
    InputMCShowerProducer: \"mcreco\"
    InputFluxProducer: \"generator\"
    InputMCTruthProducer: \"generator\"
    InputRecoTrackProducer: \"trackReco\"
    InputVtxTrackerProducer: \"trackReco\"
    InputRecoShowerProducer: \"showerreco\"
    InputPotProducer: \"generator\"
    InputPfPartShowerProducer: \"dl\"
    InputHitsShowerProducer: \"dl\"
    InputClusterShowerProducer: \"dl\"
    InputAssShowerProducer: \"showerreco\"
    InputAssDLShowerProducer: \"dl\"
    InputVtxShowerProducer: \"dl\"
    """
    lfcfg = open("savecropvariables.cfg",'w')
    print(savevariables_cfg, file=lfcfg)
    lfcfg.close()
    lfpset = larcv.CreatePSetFromFile( "savecropvariables.cfg", "SaveCropVariables" )

    # ----Algos-----
    savecutvariables.configure(lfpset)
    savecutvariables.initialize(args.showerreco_ana)

    nentries = io.get_n_entries()
    print("Num Entries: ",nentries)

    for ientry in xrange(nentries):

        io.read_entry(ientry)
        ioll.go_to(ientry)
        ioforward.read_entry(ientry)
        print("ON ENTRY: ",ientry)
        savecutvariables.process(io,ioll,ioforward,ientry)
        print(" ")
    # ientry = 16
    # io.read_entry(ientry)
    # ioll.go_to(ientry)
    # ioforward.read_entry(ientry)
    # print("ON ENTRY: ",ientry)
    # savecutvariables.process(io,ioll,ioforward,ientry)
    # print(" ")

    savecutvariables.finalize()
    io.finalize()
    ioll.close()
    ioforward.finalize()

    print("FIN")

if __name__ == '__main__':
    main()
