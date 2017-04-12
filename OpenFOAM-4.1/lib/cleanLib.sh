#!/bin/bash

source libEnv.sh

cd FoamFourierAnalysis/$FFTW_LIB

make clean

make distclean

rm -rf build.stamp

cd ../..

wclean

#
#END-OF-FILE
#

