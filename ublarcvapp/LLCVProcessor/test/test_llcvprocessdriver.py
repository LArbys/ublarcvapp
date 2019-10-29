from __future__ import print_function
import os,sys

import ROOT

from ublarcvapp import ublarcvapp

llcvdriver = ublarcvapp.llcv.LLCVProcessDriver("LLCVProcessDriver")
llcvdriver.configure( "llcvdriver.cfg" )
llcvdriver.initialize()


llcvdriver.process_entry(0)
