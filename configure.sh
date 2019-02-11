#!/bin/bash


export UBLARCVAPP_BASEDIR=$PWD
export UBLARCVAPP_LIBDIR=${UBLARCVAPP_BASEDIR}/build/lib
export UBLARCVAPP_INCDIR=${UBLARCVAPP_BASEDIR}/build/include

# add to python path
[[ ":$LD_LIBRARY_PATH:" != *":${UBLARCVAPP_LIBDIR}:"* ]] && export LD_LIBRARY_PATH="${UBLARCVAPP_LIBDIR}:${LD_LIBRARY_PATH}"

# add to python path
[[ ":$PYTHONPATH:" != *":${UBLARCVAPP_BASEDIR}/python:"* ]] && export PYTHONPATH="${UBLARCVAPP_BASEDIR}/python:${PYTHONPATH}"

