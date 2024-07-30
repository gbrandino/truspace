#!/bin/bash
ROOTDIR=$PWD

cd ${ROOTDIR}/src/truspace && make -f Makefile_linux test 
