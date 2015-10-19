#!/bin/bash	

echo "Runnnig PSE convergence test on spherical meshes ..."

rm -f spherePSETest*Out.txt

for i in `seq 1 5`
do
	mpirun -np 6 spherePSETest.exe 3 $i 2>&1 | tee -a spherePSETestIcosTriOut.txt
done

for i in `seq 1 6`
do 
	mpirun -np 6 spherePSETest.exe 4 $i 2>&1 | tee -a spherePSETestCubedSphereOut.txt
done
