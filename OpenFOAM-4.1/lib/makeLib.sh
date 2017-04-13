#!/bin/bash

source ./libEnv.sh

#
# Make FFT Library
#
if [ ! -e fftw.stamp ]
then
    THIS_DIR=`pwd`
    cd $THIS_DIR/FoamFourierAnalysis/$FFTW_LIB
    
    CFLAGS=-fPIC\\
    CXXFLAGS=-fPIC\\
    ./configure --prefix=$FOAM_USER_LIBBIN/$FFTW_LIB --enable-shared
    make
    make install
    
    cd $FOAM_USER_LIBBIN
    libs=`ls $FFTW_LIB/lib/lib*.so*`
    ls $libs
    
    for lib in $libs
    do
        ln -s $lib
    done
    
    cd $THIS_DIR
    touch fftw.stamp
fi

#
# Make libAcoustics library
#
wmake libso

#
#END-OF-FILE
#

