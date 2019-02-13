from __future__ import print_function
import os,sys
import ROOT

from larlite import larlite
from larcv import larcv
from ublarcvapp import ublarcvapp

larcv_input = sys.argv[1]
larlite_input = sys.argv[2]

iolarcv   = larcv.IOManager( larcv.IOManager.kREAD )
iolarlite = ublarcvapp.LArliteManager( larlite.storage_manager.kREAD )

iolarcv.add_in_file( larcv_input )
iolarcv.initialize()
iolarcv.set_verbosity(0)

iolarlite.add_in_filename( larlite_input )
iolarlite.open()
iolarlite.set_verbosity(0)

print("test first entry")
iolarcv.read_entry(0)
iolarcv.get_data(larcv.kProductImage2D,"wire")
print("larcv RSE: {} {} {}".format(iolarcv.event_id().run(),iolarcv.event_id().subrun(),iolarcv.event_id().event()))
iolarlite.syncEntry( iolarcv )
print("larlite index: {}".format( int(iolarlite.get_index()) ) )
print("larlite RSE: {} {} {}".format(iolarlite.run_id(),iolarlite.subrun_id(),iolarlite.event_id()))

print("test 1")
iolarcv.read_entry(5)
iolarcv.get_data(larcv.kProductImage2D,"wire")
print("larcv RSE: {} {} {}".format(iolarcv.event_id().run(),iolarcv.event_id().subrun(),iolarcv.event_id().event()))
iolarlite.syncEntry( iolarcv )
print("larlite index: {}".format( int(iolarlite.get_index()) ) )
print("larlite RSE: {} {} {}".format(iolarlite.run_id(),iolarlite.subrun_id(),iolarlite.event_id()))

print("test 1: repeat")
iolarcv.read_entry(5)
iolarcv.get_data(larcv.kProductImage2D,"wire")
print("larcv RSE: {} {} {}".format(iolarcv.event_id().run(),iolarcv.event_id().subrun(),iolarcv.event_id().event()))
iolarlite.syncEntry( iolarcv )
print("larlite index: {}".format( int(iolarlite.get_index()) ) )
print("larlite RSE: {} {} {}".format(iolarlite.run_id(),iolarlite.subrun_id(),iolarlite.event_id()))

# second test
print("test 2")
iolarcv.read_entry(3)
iolarcv.get_data(larcv.kProductImage2D,"wire")
iolarlite.syncEntry( iolarcv )
print("larlite index: {}".format( int(iolarlite.get_index()) ) )
print("larlite RSE: {} {} {}".format(iolarlite.run_id(),iolarlite.subrun_id(),iolarlite.event_id()))

print("END")
