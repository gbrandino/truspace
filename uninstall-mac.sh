#!/bin/bash
ROOTDIR=$PWD
CLN=cln-1.3.6
GINAC=ginac-1.8.6

cd $ROOTDIR/$CLN
make distclean
cd $ROOTDIR/$GINAC
make distclean
cd $ROOTDIR/src/nrg
make -f Makefile_mac clean
cd $ROOTDIR/src/truspace
make -f Makefile_mac clean
cd $ROOTDIR/src/truspace_extension_example
make -f Makefile_mac clean
cd $ROOTDIR/dep
rm -rf *
rm -rf  $ROOTDIR/bin
rm -rf  $ROOTDIR/lib
rm -rf  $ROOTDIR/include
