from __future__ import print_function

import os,sys

import ROOT as rt
from larcv import larcv
from ublarcvapp import ublarcvapp

inputfile = "out_larcv_test.root"

matchalgo = ublarcvapp.dltagger.MRCNNMatch()

io = larcv.IOManager(larcv.IOManager.kREAD, "IO" )
io.add_in_file( inputfile )
io.initialize()

io.read_entry(0)

ev_masks = io.get_data(larcv.kProductClusterMask, "mrcnn_masks" )

mask_vv = ev_masks.as_vector()
print("Number of masks: ",mask_vv.at(0).size())

matchdata = ublarcvapp.dltagger.MaskMatchData( 0, 0, mask_vv.at(0).at(0) )
