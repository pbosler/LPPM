#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=0:30:00
#SBATCH --account=fy150039
#SBATCH --job-name=adv-gaussHills
#SBATCH --partition=ec
#SBATCH --mail-type=ALL

cd $HOME/LPPM

cat <<EOF > AdvectGaussianHills.namelist
&meshDefine
	initNest = 7
	AMR = 0
	panelKind = 3
	amrLimit = 5
	tracerMassTol = 0.008
	tracerVarTol = 0.09
/

&timestepping
	dt = 0.03			! days (note : 1 period = 12 days)
	tfinal = 12.0		! days
	remeshInterval = 20
	resetAlphaInterval = 50000
/

&fileIO
	outputDir = '$HOME/modelData'
	jobPrefix = 'gHills_conv_'
	frameOut = 200
/
EOF

mpirun -np 8 advectGaussHills.exe