#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=1:30:00
#SBATCH --account=fy150039
#SBATCH --job-name=mv-sphere
#SBATCH --partition=ec
#SBATCH --mail-type=ALL

#outputLoc=/fscratch/pabosle/movingVortices

#cp $HOME/LPPM/advectMovingVortices.exe $outputLoc/.

#cd $outputLoc

cd $HOME/LPPM

cat <<EOF > MovingVortices2Direct.namelist
&meshDefine
	initNest = 7
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
	outputDir = '$HOME/modelData/'
	jobPrefix = 'mv_conv_direct_'
	frameOut = 200
/

EOF

mpirun -np 8 ./advectMovingVorticesDirect.exe