#!/bin/bash

source ./libEnv.sh

THIS_DIR=`pwd`

#
# Make FFT Library
#

cd $THIS_DIR/FoamFourierAnalysis/$FFTW_LIB

if [ ! -e build.stamp ]
then
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
    
    touch build.stamp
fi

#
# Make libAcoustics library
#

#
cd $THIS_DIR
wmake libso

#
#END-OF-FILE
#

