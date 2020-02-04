import os,sys
import ROOT

from larcv import larcv
from ublarcvapp import ublarcvapp

print ublarcvapp.ubdllee.DLMerger

merger = ublarcvapp.ubdllee.DLMerger("DLMerger")
merger.configure( "dlmerger.cfg",
                  "fileset.json",
                  "../../../../testdata/mcc9_v13_nueintrinsic_overlay_run1/complete" )

merger.initialize()

merger.batch_process()

merger.finalize()
