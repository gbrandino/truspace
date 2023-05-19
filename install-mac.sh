#!/bin/bash
ROOTDIR=$PWD
CLN=cln-1.3.6
GINAC=ginac-1.8.6

if [ -z "$CC" ]
then
      export CC=gcc;
      export CXX=g++;
      echo "Compiling with $CC and $CXX";
      sleep 5
else
      echo "Compiling with $CC and $CXX";
      sleep 5
fi


cd $ROOTDIR/$CLN && \
touch configure.ac aclocal.m4 configure Makefile.am Makefile.in && \
./configure --prefix=$ROOTDIR/dep/$CLN/ && make -j4 && make install && \
cd $ROOTDIR/$GINAC && \
touch configure.ac aclocal.m4 configure Makefile.am Makefile.in && \
CLN_LIBS="-L$ROOTDIR/dep/$CLN/lib -lcln" CLN_CFLAGS=-I$ROOTDIR/dep/$CLN/include ./configure --prefix=$ROOTDIR/dep/$GINAC && make -j4 && make install && \
cd $ROOTDIR/src/nrg && make -f Makefile_mac serial && \
cd $ROOTDIR/src/truspace && make -f Makefile_mac serial && \
cd $ROOTDIR/src/truspace_extension_example && make -f Makefile_mac 

