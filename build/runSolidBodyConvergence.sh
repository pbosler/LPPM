#!/bin/bash

lpmExeDir=/ascldap/users/pabosle/lpm-v2/fortran/build
outputDir=/fscratch/pabosle/BVE
runDir=/fscratch/pabosle/BVE

cd $runDir
cp $lpmExeDir/bveSolidBody.exe .

iNests=(1 2 3 4 5 6 7)
dts=(0.01 0.01 0.005 0.0025 0.00125 0.000625 0.0003125)
nTimesteps=(100 100 200 400 800 1600 3200)

for (( i = 1; i <= 7; i++))
do
cat <<EOF > solidBodyConvTest.namelist
&meshDefine
faceKind = 3
initNest = $i
maxNest = $i
amrLimit = 0
/
&timestepping
dt = ${dts[i-1]}
tfinal = 1.0
/
&fileIO
outputDir='${outputDir}'
outputRoot='bveSolidBodyConv'
frameOut = ${nTimesteps[i-1]}
/

EOF

cat solidBodyConvTest.namelist
	
done


