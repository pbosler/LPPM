#!/bin/bash

rm -f mlsqIcosTriConvTestOut.txt

for i in `seq 1 6`
do
	mpirun -np 1 sphereTrivarMLSQ.exe 3 $i 2>&1 | tee -a mlsqIcosTriConvTestOut.txt
done;

rm -f mlsqCubedSphereConvTestOut.txt

for i in `seq 1 7`
do
	mpirun -np 1 sphereTrivarMLSQ.exe 4 $i 2>&1 | tee -a mlsqCubedSphereConvTestOut.txt
done;