#!/bin/bash

export CC=mpicc
export CXX=mpicxx
export FC=mpif90

#export SOURCE_ROOT=/Users/pabosle/Desktop/LPM-v2/fortran

rm -rf CMakeCache.txt
rm -rf CMakeFiles/

cmake \
	-D CMAKE_BUILD_TYPE:STRING=DEBUG \
	-D CMAKE_INSTALL_PREFIX=$SOURCE_ROOT/install \
	..
	