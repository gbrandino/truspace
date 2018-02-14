#!/bin/bash
ROOTDIR=$PWD
CLN=cln-1.3.4
GINAC=ginac-1.7.1
BLAS=OpenBLAS-0.2.19
SCALAPACK=scalapack-2.0.0

cd $ROOTDIR/$CLN
make distclean
cd $ROOTDIR/$GINAC
make distclean
cd $ROOTDIR/$BLAS
make clean
cd $ROOTDIR/$SCALAPACK
make clean
cd $ROOTDIR/src
make clean
cd $ROOTDIR/dep
rm -rf *
rm -f  $ROOTDIR/truspace_mpi.x

