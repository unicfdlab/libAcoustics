#!/bin/bash

source libEnv.sh

rm -rf fftw.stamp
cd FoamFourierAnalysis/$FFTW_LIB
make clean
make distclean
cd ../..

wclean

#
#END-OF-FILE
#

