#!/bin/bash

source libEnv.sh

cd FoamFourierAnalysis/$FFTW_LIB

make clean

make distclean

cd ../..

wclean

#
#END-OF-FILE
#

