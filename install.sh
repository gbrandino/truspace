#!/bin/bash
ROOTDIR=$PWD
CLN=cln-1.3.4
GINAC=ginac-1.7.1
BLAS=OpenBLAS-0.2.19

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

if which gfortran > /dev/null; then
    FC=gfortran
else
    echo "gfortran not found!"
    echo "Please specify the full path to your fortran compiler"
    read FC
fi

i="0"
while [ $i -lt 1 ]
do
    if [ -x $FC ]; then
        i=$[$i+1]
    else
        echo $FC " not found! Please specify a proper path"
        read FC
    fi
done


cd $ROOTDIR/$CLN
touch configure.ac aclocal.m4 configure Makefile.am Makefile.in
./configure --prefix=$ROOTDIR/dep/$CLN/ && make -j4 && make install
cd $ROOTDIR/$GINAC
touch configure.ac aclocal.m4 configure Makefile.am Makefile.in
CLN_LIBS="-L$ROOTDIR/dep/$CLN/lib -lcln" CLN_CFLAGS=-I$ROOTDIR/dep/$CLN/include ./configure --prefix=$ROOTDIR/dep/$GINAC && make -j4 && make install
cd $ROOTDIR/$BLAS
make FC=$FC  && make install PREFIX=$ROOTDIR/dep/$BLAS
cd $ROOTDIR/src
make
cp truspace.x $ROOTDIR/

