#!/bin/bash

export CC=mpicc
export CXX=mpicxx
export FC=mpifort

rm -rf CMakeCache.txt
rm -rf CMakeFiles/

cmake \
	-D CMAKE_BUILD_TYPE:STRING=RELEASE \
	-D CMAKE_INSTALL_PREFIX:FILEPATH=/Users/pabosle/Desktop/LPM-1.0/ \
.