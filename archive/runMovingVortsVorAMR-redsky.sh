#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=0:40:00
#SBATCH --account=fy150039
#SBATCH --job-name=mvVorAMR-sphere
#SBATCH --partition=ec
#SBATCH --mail-type=ALL

outputLoc=/fscratch/pabosle/movingVortices

cp $HOME/LPPM/advectMovingVorticesWithVorticityRefinement.exe $outputLoc/.

cd $outputLoc

cat <<EOF > MovingVortices2.namelist
&meshDefine
	initNest = 3
	AMR = 3
	panelKind = 3
	amrLimit = 3
	maxCircTol = 0.02
	vortVarTol = 1.0e20
	tracerMassTol = 1.0e20
	tracerVarTol = 0.03
	lagVarTol = 1.0e20
/

&timestepping
	dt = 0.03			! days (note : 1 period = 12 days)
	tfinal = 12.0		! days
	remeshInterval = 20
	resetAlphaInterval = 4000
/

&fileIO
	outputDir = '${outputLoc}/'
	jobPrefix = 'mv_vorTracer_rm20_intPuzz_'
	frameOut = 100
/

EOF

mpirun -np 8 ./advectMovingVorticesWithVorticityRefinement.exe