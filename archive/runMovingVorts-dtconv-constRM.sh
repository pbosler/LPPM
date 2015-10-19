#!/bin/bash

#SBATCH --nodes=2
#SBATCH --time=1:00:00
#SBATCH --account=fy150039
#SBATCH --partition=ec
#SBATCH --job-name=mv-dtConv
#SBATCH --mail-type=ALL

outputLoc=/fscratch/pabosle

cp $HOME/LPPM/advectMovingVortices.exe $outputLoc/.

cd $outputLoc

DT=(0.48 0.24 0.12 0.06 0.03 0.015 0.0075 0.00375)
FO=(12 25 50 100 200 400 800 1600)

namelistfilename='mv_dtConv_constRM.namelist'

for i in `seq 0 7`
do
cat <<EOF > ${namelistfilename}
&meshDefine
	initNest = 4
	AMR = 0
	panelKind = 3
	amrLimit = 4
	tracerMassTol = 1.0e20
	tracerVarTol = 1.0e20
	lagVarTol = 1.0e20
	tracerRefVal = 1.0
	tracerRefTol = 0.0
/

&timestepping
	dt = ${DT[i]}
	tfinal = 12.0
	remeshInterval = 20
	resetAlphaInterval = 20000
/

&fileIO
	outputDir = '$outputLoc/'
	jobPrefix = 'mv_dtConv_${DT[i]}_rm20_'
	frameOut = ${FO[i]}
/

EOF

mpirun -np 16 advectMovingVortices.exe ${namelistfilename}

done