#!/bin/bash
ROOTDIR=$PWD
CLN=cln-1.3.6
GINAC=ginac-1.8.6
BLAS=OpenBLAS-0.3.17
SCALAPACK=scalapack-2.0.0

cd $ROOTDIR/$CLN
make distclean
cd $ROOTDIR/$GINAC
make distclean
cd $ROOTDIR/$BLAS
make clean
cd $ROOTDIR/$SCALAPACK
make clean
cd $ROOTDIR/src/nrg
make -f Makefile_linux clean
cd $ROOTDIR/src/truspace
make -f Makefile_linux clean
cd $ROOTDIR/src/truspace_extension_example
make -f Makefile_linux clean
cd $ROOTDIR/dep
rm -rf *
rm -f  $ROOTDIR/bin
rm -f  $ROOTDIR/lib

