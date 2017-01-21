#!/bin/bash
ROOTDIR=$PWD
CLN=cln-1.3.4
GINAC=ginac-1.7.1
BLAS=OpenBLAS-0.2.19

cd $ROOTDIR/$CLN
./configure --prefix=$ROOTDIR/dep/$CLN/ && make -j4 && make install
cd $ROOTDIR/$GINAC 
CLN_LIBS="-L$ROOTDIR/dep/$CLN/lib -lcln" CLN_CFLAGS=-I$ROOTDIR/dep/$CLN/include ./configure --prefix=$ROOTDIR/dep/$GINAC && make -j4 && make install 
cd $ROOTDIR/$BLAS
make  && make install PREFIX=$ROOTDIR/dep/$BLAS
cd $ROOTDIR/src
make
cp main.x $ROOTDIR/
