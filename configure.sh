#!/bin/bash


export UBLARCVAPP_BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export UBLARCVAPP_LIBDIR=${UBLARCVAPP_BASEDIR}/build/lib
export UBLARCVAPP_INCDIR=${UBLARCVAPP_BASEDIR}/build/include

# add to python path
[[ ":$LD_LIBRARY_PATH:" != *":${UBLARCVAPP_LIBDIR}:"* ]] && export LD_LIBRARY_PATH="${UBLARCVAPP_LIBDIR}:${LD_LIBRARY_PATH}"

# add bin to path
[[ ":$PATH:" != *":${UBLARCVAPP_BASEDIR}/bin:"* ]] && export PATH="${UBLARCVAPP_BASEDIR}/bin:${PATH}"

# add to python path
[[ ":$PYTHONPATH:" != *":${UBLARCVAPP_BASEDIR}/python:"* ]] && export PYTHONPATH="${UBLARCVAPP_BASEDIR}/python:${PYTHONPATH}"

