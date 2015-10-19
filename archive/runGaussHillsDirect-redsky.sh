#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --account=fy150039
#SBATCH --job-name=adv-gaussHillsDirect
#SBATCH --partition=ec
#SBATCH --mail-type=ALL

cd $HOME/LPPM

for i in `seq 1 7`;
do
cat <<EOF > AdvectGaussHillsDirect.namelist
&meshDefine
	initNest = ${i}
	AMR = 0
	panelKind = 3
	amrLimit = 0
	tracerMassTol = 0.1
	tracerVarTol = 1.0
/

&timestepping
	dt = 0.03			! days (note : 1 period = 12 days)
	tfinal = 12.0		! days
	remeshInterval = 20
	resetAlphaInterval = 5000
/

&fileIO
	outputDir = '$HOME/modelData/'
	jobPrefix = 'gHills_conv_DIRECT_'
	frameOut = 200
/
EOF

mpirun -np 8 advectGaussHillsDirect.exe

done