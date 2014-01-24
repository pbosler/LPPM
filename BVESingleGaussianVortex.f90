program BVESolidBody
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	This program solves the barotropic vorticity equation for the case of rigid rotation about the z-axis.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
!----------------
use NumberKindsModule
use LoggerModule
use SphereMeshModule
use RefineRemeshModule
use TracerSetupModule
use VTKOutputModule
use BVESetupModule
use BVEDirectSumModule

implicit none

include 'mpif.h'

!
!	mesh variables
!
type(SphereMesh) :: sphere
integer(kint) :: panelKind, initNest, AMR, nTracer, problemKind
!
!	time stepping variables
!
type(BVERK4Data) :: bveRK4
real(kreal) :: t, dt, tfinal
integer(kint) :: timesteps, timeJ
!
!	tracer variables
!
type(tracerSetup) :: cosBell
real(kreal) :: h0, RR, lat0, lon0
integer(kint) :: tracerID
!
!	vorticity variables
!
type(BVESetup) :: gaussVort
real(kreal) :: vortLat, vortLon, bb, maxVort
!
!	refinement variables
!
type(RefinementSetup) :: tracerRefine
real(kreal) :: tracerMaxTol, tracerVarTol
integer(kint) :: refinementLimit

type(RefinementSetup) :: vortRefine
real(kreal) :: circMaxTol, vortVarTol

type(RefinementSetup) :: flowMapRefine
real(kreal) :: lagVarTol

!
!	remeshing variables
!
integer(kint) :: remeshInterval, remeshCounter
!
!	User input
!
character(len = 128) :: namelistInputFile = 'gaussianVortex.namelist'
integer(kint) :: readStat
!
!	Output variables
!
type(VTKSource) :: vtkOut
character(len = 128) :: vtkRoot, vtkFile, outputDir, jobPrefix
character(len = 56 ) :: amrString
integer(kint) :: frameCounter, frameOut
!
!	logging
!
type(Logger) :: exeLog
character(len=28) :: logKey
character(len=128) :: logString
!
!	general and mpi variables
!
integer(kint) :: i, j, k
integer(kint) :: mpiErrCode
integer(kint), parameter  :: BROADCAST_INT_SIZE = 5, BROADCAST_REAL_SIZE = 11
integer(kint) :: broadcastIntegers(BROADCAST_INT_SIZE)
real(kreal) :: broadcastReals(BROADCAST_REAL_SIZE)

namelist /sphereDefine/ panelKind, initNest, AMR, tracerMaxTol, tracerVarTol, &
						refineMentLimit, circMaxTol, vortVarTol, lagVarTol
namelist /vorticityDefine/ 	vortLat, vortLon, bb, maxVort					
namelist /timeStepping/ tfinal, dt,  remeshInterval
namelist /fileIO/ outputDir, jobPrefix, frameOut

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	PROGRAM START
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numProcs,mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD,procRank,mpiErrCode)

call InitLogger(exeLog,procRank)
!
!	read user input from namelist file, broadcast to all processes
!
call ReadNameListFile(procRank)
call MPI_BCAST(broadcastIntegers,BROADCAST_INT_SIZE,MPI_INTEGER,0,MPI_COMM_WORLD,mpiErrCode)
panelKind = broadcastIntegers(1)
AMR = broadcastIntegers(2)
initNest = broadcastIntegers(3)
refinementLimit = broadcastIntegers(4)
remeshInterval = broadcastIntegers(5)

call MPI_BCAST(broadcastReals,BROADCAST_REAL_SIZE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiErrCode)
tracermaxTol = broadcastReals(1)
tracerVarTol = broadcastReals(2)
tfinal = broadcastReals(3)*ONE_DAY ! convert tfinal to seconds
circMaxTol = broadcastReals(4)
vortVarTol = broadcastReals(5)
dt = broadcastReals(6) * ONE_DAY
lagVarTol = broadcastReals(7)
vortLat = broadcastReals(8)
vortLon = broadcastReals(9)
bb = broadcastReals(10)
maxVort = broadcastReals(11)

!
!	define test case in accordance with Williamson et al., JCP 1992, test case 1
!
problemKind = BVE_SOLVER
nTracer = 2
h0 = 1000.0_kreal	! height of cosine bell center
lat0 = 0.0_kreal	! latitude of cosine bell center
lon0 = 3.0_kreal*PI/2.0_kreal ! longitude of cosine bell center
RR = EARTH_RADIUS/3.0_kreal ! radius of cosine bell
tracerID = 1
call New(cosBell, COS_BELL_NINT, COS_BELL_NREAL)
call InitCosineBellTracer(cosBell,lat0,lon0,RR,h0,tracerID)

call New(gaussVort,SINGLE_VORTEX_NINT,SINGLE_VORTEX_NREAL)
call InitSingleGaussianVortex(gaussVort,vortLat,vortLon,bb,maxVort)

!
!	initialize sphere
!
call New(sphere,panelKind,initNest,AMR,nTracer,problemKind)
call SetCosineBellTracerOnMesh(sphere,cosBell)
call SetInitialLatitudeTracerOnMesh(sphere,2)
call SetSingleGaussianVortexOnMesh(sphere,gaussVort)
if ( AMR > 0 ) then
	call New(tracerRefine,refinementLimit,tracermaxTol,tracerVarTol,TRACER_REFINE,tracerID)
	call New(vortRefine,refinementLimit,circMaxTol,vortVarTol,RELVORT_REFINE)
	call New(flowMapREfine,refinementLimit,100000.0_kreal,lagVarTol,FLOWMAP_REFINE)
	call InitialRefinement(sphere,tracerRefine,SetCosineBellTracerOnMesh, cosBell, &
						   vortRefine, SetSingleGaussianVortexOnMesh,gaussVort)
	call SetInitialLatitudeTracerOnMesh(sphere,2)					   
	if ( panelKind == QUAD_PANEL) then
		write(amrString,'(A,I1,A,I0.2)') 'quadAMR_',initNest,'to',initNest+refinementLimit
	endif
else
	call New(tracerRefine)
	call New(vortRefine)
	call New(flowMapRefine)
	if ( panelKind == QUAD_PANEL ) then
		write(amrString,'(A,I1,A)') 'quadUnif',initNest,'_'
	endif
endif
!
!	initialize output, output t = 0 data
!
frameCounter = 0
if ( procRank == 0 ) then
	call LogStats(sphere,exeLog,'initial mesh :')
	write(vtkRoot,'(A,A,A,A,A)') trim(outputDir), '/vtkOut/',trim(jobPrefix),trim(amrString),'_'
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),frameCounter,'.vtk'

	call New(vtkOut,sphere,vtkFile,'gaussVort')
	call VTKOutput(vtkOut,sphere)
	frameCounter = frameCounter + 1
endif
!
!	intialize timestepping
!
call New(bveRK4,sphere,numProcs)
timesteps = floor(tfinal/dt)
t = 0.0_kreal
remeshCounter = 0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	RUN THE PROBLEM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
do timeJ = 0, timesteps - 1
	!
	!	remesh if necessary
	!
	if ( mod(timeJ+1,remeshInterval) == 0 ) then
		remeshCounter = remeshCounter + 1
		!
		!	delete objects associated with old mesh
		!
		call Delete(bveRK4)
		call Delete(vtkOut)
		!
		! build new mesh
		!
		call LagrangianRemesh(sphere, SetSingleGaussianVortexOnMesh, gaussVort, vortRefine, &
									  SetCosineBellTracerOnMesh, cosBell, tracerRefine, &
									  flowMapRefine)
		!
		!	create new associated objects
		!
		call New(bveRK4,sphere,numProcs)
		call New(vtkOut,sphere,vtkFile,'gaussVort')

		write(logString,'(A,I4)') 'remesh ', remeshCounter
		if ( procRank == 0 ) call LogStats(sphere,exelog,trim(logString))
	endif
	!
	!	advance timestep
	!
	call BVERK4Timestep(bveRK4, sphere, dt, procRank, numProcs)

	t = real(timeJ+1,kreal)*dt
	!
	!	output timestep data
	!
	if ( procRank == 0 .AND. mod(timeJ+1,frameOUt) == 0  ) then
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,'day = ',t/ONE_DAY)
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
		call UpdateFileName(vtkOut,vtkFile)
		call VTKOutput(vtkOut,sphere)
		frameCounter = frameCounter + 1
	endif
enddo
!
!	output final data
!

! TO DO : output final data 



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	Clear memory and Finalize
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
call Delete(bveRK4)
if ( procRank == 0 ) call Delete(vtkOut)
call Delete(vortRefine)
call Delete(tracerRefine)
call Delete(flowmapRefine)
call Delete(cosBell)
call Delete(gaussVort)
call Delete(exeLog)
call MPI_FINALIZE(mpiErrCode)

contains

subroutine ReadNamelistfile(rank)
	integer(kint), intent(in) :: rank
	if ( rank == 0 ) then
		open(unit = READ_UNIT, file = namelistInputFile, status = 'OLD', action = 'READ',iostat= readStat)
			if ( readstat /= 0 ) then
				stop 'cannot read namelist file.'
			endif
			read(READ_UNIT,nml=sphereDefine)
			rewind(READ_UNIT)
			read(READ_UNIT,nml=timestepping)
			rewind(READ_UNIT)
			read(READ_UNIT,nml=fileIO)

			broadcastIntegers(1) = panelKind
			broadcastIntegers(2) = AMR
			broadcastIntegers(3) = initNest
			broadcastIntegers(4) = refinementLimit
			broadcastIntegers(5) = remeshInterval

			broadcastReals(1) = tracerMaxTol
			broadcastReals(2) = tracerVarTol
			broadcastReals(3) = tfinal
			broadcastReals(4) = circMaxTol
			broadcastReals(5) = vortVarTol
			broadcastReals(6) = dt
			broadcastReals(7) = lagVarTol
			broadcastReals(8) = vortLat
			broadcastReals(9) = vortLon
			broadcastReals(10) = bb
			broadcastReals(11) = maxVort
	endif
end subroutine

subroutine InitLogger(alog,rank)
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	if ( rank == 0 ) then
		call New(aLog,DEBUG_LOGGING_LEVEL)
	else
		call New(aLog,WARNING_LOGGING_LEVEL)
	endif
        write(logKey,'(A,I0.2,A)') 'EXE_LOG',procRank,' : '
end subroutine

end program
