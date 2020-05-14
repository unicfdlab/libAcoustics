#!/bin/bash
#
# Generates time history of acoustic pressure at predefined points:
# - (0 0 1)
# - (0 0 2)
# - (0 0 3)
# - (0 0 4)
# - (0 0 5)
#
export PYTHONPATH=$PYTHONPATH:../../Source/lib/samplingBEM/python2.7

C0=100
FREQ=100
INFILE="complexPressureData/triSurfaceSamplingFine_surfacePressure/$FREQ/dataWithMesh.msh"
python pulsosphere.py $C0 $FREQ $INFILE

#
#END-OF-FILE
#

