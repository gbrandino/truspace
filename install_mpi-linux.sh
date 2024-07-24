#!/bin/bash
ROOTDIR=$PWD
CLN=cln-1.3.7
GINAC=ginac-1.8.7
BLAS=OpenBLAS-0.3.17
SCALAPACK=scalapack-2.0.0

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

if type "gfortran" > /dev/null 2>&1; then
    FC=$(which gfortran)
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


cd $ROOTDIR/$CLN  && \
touch configure.ac aclocal.m4 configure Makefile.am Makefile.in  && \
./configure --prefix=$ROOTDIR/dep/$CLN/ && make -j4 && make install  && \
cd $ROOTDIR/$GINAC  && \
touch configure.ac aclocal.m4 configure Makefile.am Makefile.in  && \
CLN_LIBS="-L$ROOTDIR/dep/$CLN/lib -lcln" CLN_CFLAGS=-I$ROOTDIR/dep/$CLN/include ./configure --prefix=$ROOTDIR/dep/$GINAC && make -j4 && make install  && \
cd $ROOTDIR/$BLAS  && \ 
make FC=$FC  && make install PREFIX=$ROOTDIR/dep/$BLAS  && \
cd $ROOTDIR/$SCALAPACK  && \
BLASDIR=$ROOTDIR/dep/$BLAS/lib make  && \
cd $ROOTDIR/src/nrg && make -f Makefile_linux mpi && \
cd $ROOTDIR/src/truspace && make -f Makefile_linux mpi && \
cd $ROOTDIR/src/truspace_extension_example && make -f Makefile_linux mpi

