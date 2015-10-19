#!/bin/bash

export CC=mpicc
export CXX=mpicxx
export FC=mpif90

export NETCDF=/opt/netcdf-4.3.2

export SOURCE_ROOT=$HOME/LPPM

rm -rf CMakeCache.txt
rm -rf CMakeFiles/

cmake \
	-D CMAKE_BUILD_TYPE:STRING=DEBUG \
	-D CMAKE_INSTALL_PREFIX=$SOURCE_ROOT/install \
	$SOURCE_ROOT