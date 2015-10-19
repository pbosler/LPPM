#!/bin/bash

#SBATCH --nodes=2
#SBATCH --time=0:40:00
#SBATCH --account=fy150039
#SBATCH --job-name=mvDirect
#SBATCH --partition=ec
#SBATCH --mail-type=ALL

outputLoc=/fscratch/pabosle

cp $HOME/LPPM/advectMovingVorticesDirect.exe $outputLoc/.

cd $outputLoc

for i in `seq 3 5`;
do

cat <<EOF > MovingVortices2Direct.namelist
&meshDefine
	initNest = ${i}
	AMR = 0
	panelKind = 3
	amrLimit = 0
	maxCircTol = 0.2
	vortVarTol = 1.0
	tracerMassTol = 2.0
	tracerVarTol = 1.0
	lagVarTol = 2.2
/

&timestepping
	dt = 0.03			! days (note : 1 period = 12 days)
	tfinal = 12.0		! days
	remeshInterval = 20
	resetAlphaInterval = 4000
/

&fileIO
	outputDir = '$outputLoc/'
	jobPrefix = 'mv_conv_direct_'
	frameOut = 200
/

EOF

mpirun -np 16 ./advectMovingVorticesDirect.exe

done