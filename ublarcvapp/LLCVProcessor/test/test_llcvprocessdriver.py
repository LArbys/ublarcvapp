from __future__ import print_function
import os,sys

import ROOT

from ublarcvapp import ublarcvapp

llcvdriver = ublarcvapp.llcv.LLCVProcessDriver("LLCVProcessDriver")

print(llcvdriver)
llcvdriver.configure( "llcvdriver.cfg" )
