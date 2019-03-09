import os,sys

'''
Goal is to test if FixedCROIFromFlashAlgo works properly.
'''

import ROOT

from larlite import larlite
from larcv import larcv
from ublarcvapp import ublarcvapp

supera="../../../../testdata/ex1/supera-Run000001-SubRun006867.root"
opreco="../../../../testdata/ex1/opreco-Run000001-SubRun006867.root"

# FILE IO
# -------
iolarcv = larcv.IOManager(larcv.IOManager.kREAD, "io", larcv.IOManager.kTickBackward)
iolarlite = ublarcvapp.LArliteManager(larlite.storage_manager.kREAD)

iolarcv.add_in_file( supera )
iolarcv.initialize()
iolarcv.set_verbosity(0)

iolarlite.add_in_filename( opreco )
iolarlite.open()
iolarlite.set_verbosity(0)

# Algo Config
# -----------
config = ublarcvapp.ubdllee.FixedCROIFromFlashConfig()
algo   = ublarcvapp.ubdllee.FixedCROIFromFlashAlgo( config )

# load LArCV Entry
iolarcv.read_entry(1)
ev_wire = iolarcv.get_data(larcv.kProductImage2D,"wire")

# load larlite entry
iolarlite.syncEntry(iolarcv)
ev_opflash = iolarlite.get_data(larlite.data.kOpFlash,"simpleFlashBeam")

# pick out the intime flash
usec_min = 190*0.015625
usec_max = 320*0.015625
intimeflash = None
for iflash in xrange(ev_opflash.size()):
    opflash = ev_opflash.at(iflash)
    tusec = opflash.Time() # usec since trigger
    
    if usec_min < tusec and tusec < usec_max:
        intimeflash = opflash

# run algo
if intimeflash is not None:
    roi_v = algo.findCROIfromFlash( intimeflash )





