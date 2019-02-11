import ROOT,os
if not 'LARCV_BASEDIR' in os.environ:
    print '$LARCV_BASEDIR shell env. var. not found (run configure.sh)'
    raise ImportError

# must load dependencies first
# LARLITE
if 'LARLITE_BASEDIR' in os.environ:
    from larlite import larlite
if 'LAROPENCV_BASEDIR' in os.environ:
    from larocv import larocv
from larcv import larcv

ublarcvapp_dir = os.environ['UBLARCVAPP_BASEDIR']+"/build/lib"
# We need to load in order
for l in [x for x in os.listdir(ublarcvapp_dir) if x.endswith('.so')]:
    ROOT.gSystem.Load(l)
