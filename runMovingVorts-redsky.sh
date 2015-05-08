#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=0:20:00
#SBATCH --account=fy150039
#SBATCH --job-name=mvAMR-sphere
#SBATCH --partition=ec
#SBATCH --mail-type=ALL

outputLoc=/fscratch/pabosle/movingVortices

cp $HOME/LPPM/advectMovingVortices.exe $outputLoc/.

cd $outputLoc

cat <<EOF > MovingVortices.namelist
&meshDefine
	initNest = 4
	AMR = 3
	panelKind = 3
	amrLimit = 2
	tracerMassTol = 1.0e20 	! smaller values => more panels
	tracerVarTol =0.06	! smaller values => more panels
	lagVarTol = 1.0e20	! smaller values => more panels
	tracerRefVal = 1.0
	tracerRefTol = 0.00625 ! larger values => more panels
/

&timestepping
	dt = 0.03			! days (note : 1 period = 12 days)
	tfinal = 12.0		! days
	remeshInterval = 20
	resetAlphaInterval = 40000
/

&fileIO
	outputDir = '${outputLoc}/'
	jobPrefix = 'mv_intPuzz_'
	frameOut = 100
/

EOF

mpirun -np 8 ./advectMovingVortices.exe