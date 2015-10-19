#!/bin/bash

export CC=mpicc
export CXX=mpicxx
export FC=mpif90

rm -rf CMakeCache.txt
rm -rf CMakeFiles/

cmake \
	-D CMAKE_BUILD_TYPE:STRING=RELEASE \
	-D CMAKE_INSTALL_PREFIX:FILEPATH=/ascldap/users/pabosle/LPPM/install \
.