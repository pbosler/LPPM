#!/bin/bash

#SBATCH --nodes=2
#SBATCH --time=2:40:00
#SBATCH --account=fy150039
#SBATCH --job-name=mv-conv
#SBATCH --partition=ec
#SBATCH --mail-type=ALL

outputLoc=/fscratch/pabosle

cp $HOME/LPPM/advectMovingVortices.exe $outputLoc/.

cd $outputLoc

for i in 1 5 10 20 40 100 200
do 

cat <<EOF > mvRemeshFreq.namelist
&meshDefine
	initNest = 5
	AMR = 0
	panelKind = 3
	amrLimit = 2
	tracerMassTol = 1.0e20 	! smaller values => more panels
	tracerVarTol =0.06	! smaller values => more panels
	lagVarTol = 1.0e20	! smaller values => more panels
	tracerRefVal = 1.0
	tracerRefTol = 0.00625 ! larger values => more panels
/

&timestepping
	dt = 0.015			! days (note : 1 period = 12 days)
	tfinal = 12.0		! days
	remeshInterval = ${i}
	resetAlphaInterval = 40000
/

&fileIO
	outputDir = '$outputLoc/'
	jobPrefix = 'mv_remeshFreq_halfDt_rm${i}'
	frameOut = 400
/

EOF

mpirun -np 16 ./advectMovingVortices.exe mvRemeshFreq.namelist

done