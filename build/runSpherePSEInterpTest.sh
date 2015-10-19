#!/bin/bash

rm -f spherePSEIcosTriInterpConvTestOut.txt

for i in `seq 1 6`
do
	mpirun -np 10 spherePSEInterpTest.exe 3 $i 2>&1 | tee -a spherePSEIcosTriInterpConvTestOut.txt
done;

rm -f spherePSECubedSphereInterpConvTestOut.txt

for i in `seq 1 7`
do 
	mpirun -np 10 spherePSEInterpTest.exe 4 $i 2>&1 | tee -a spherePSECubedSphereInterpConvTestOut.txt
done;
