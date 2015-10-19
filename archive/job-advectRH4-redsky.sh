#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --account=fy150039
#SBATCH --job-name=advectRH4
#SBATCH --partition=ec
#SBATCH --mail-type=ALL

outputLoc=/gscratch1/pabosle/rh4Wave
exeLoc=/ascldap/users/pabosle/LPPM

cp $exeLoc/advectRH4.exe $outputLoc/.
cd $outputLoc

cat <<EOF > AdvectRH4.namelist
&meshDefine
initNest=6
AMR=0
panelKind=3
amrLimit=4
tracerMassTol=0.1
tracerVarTol=0.25
/

&timestepping
dt = 0.01 ! days
tfinal = 24.0 ! days
remeshInterval = 10
resetAlphaInterval = 2000000
/

&fileIO
outputDir='${outputLoc}'
jobPrefix='advectRH4_'
frameOut=4
/

EOF

mpirun -np 8 ./advectRH4.exe
