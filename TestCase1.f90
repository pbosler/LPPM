program TestCase1
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	This program is runs Williamson et al., JCP 1992, test case 1.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!

!----------------
use NumberKindsModule
use LoggerModule
use SphereMeshModule
use AdvectionModule
use RefineRemeshModule
use TracerSetupModule
use VTKOutputModule
use BVESetupModule

implicit none

!
!	mesh variables
!
type(SphereMesh) :: sphere
integer(kint) :: panelKind, initNest, AMR, nTracer, problemKind
!
!	time stepping variables
!
type(RK4Data) :: advRK4
real(kreal) :: t, dt, tfinal
integer(kint) :: timesteps, timeJ
!
!	tracer variables
!
type(tracerSetup) :: cosBell
real(kreal) :: h0, RR, lat0, lon0
integer(kint) :: tracerID
!
!	refinement variables
!
type(RefinementSetup) :: tracerRefine
real(kreal) :: tracerMaxTol, tracerVarTol
integer(kint) :: refinementLimit
!
!	remeshing variables
!
integer(kint) :: remeshInterval, remeshCounter
!
!	User input
!
character(len = 128) :: namelistFile = 'testcase1.namelist'
integer(kint) :: readStat
!
!	Output variables
!
type(VTKSource) :: vtkOut
character(len = 128) :: vtkRoot, vtkFile
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
integer(kint), parameter  :: BROADCAST_INT_SIZE = 5, BROADCAST_REAL_SIZE = 4
integer(kint) :: broadcastIntegers(BROADCAST_INT_SIZE)
real(kreal) :: broadcastReals(BROADCAST_REAL_SIZE)
type(BVESetup) :: nullVort

namelist /sphereDefine/ panelKind, initNest, AMR, tracerMaxTol, tracerVarTol, refineMentLimit
namelist /timeStepping/ tfinal, dt, remeshInterval
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
dt = broadcastreals(4)
!
!	define test case in accordance with Williamson et al., JCP 1992
!
problemKind = ADVECTION_SOLVER
nTracer = 1
h0 = 1000.0_kreal	! height of cosine bell center
lat0 = 0.0_kreal	! latitude of cosine bell center
lon0 = 3.0_kreal*PI/2.0_kreal ! longitude of cosine bell center
RR = EARTH_RADIUS/3.0_kreal ! radius of cosine bell
tracerID = 1
call New(cosBell, COS_BELL_NINT, COS_BELL_NREAL)
call InitCosineBellTracer(cosBell,lat0,lon0,RR,h0,tracerID)
!
!	initialize sphere
!
call New(sphere,panelKind,initNest,AMR,nTracer,problemKind)
call SetCosineBellTracerOnMesh(sphere,cosBell)
if ( AMR > 0 ) then
	call New(tracerRefine,refinementLimit,tracermaxTol,tracerVarTol,REFINE_TRACER,tracerID)
	call InitialRefinement(sphere,tracerRefine,SetCosineBellTracerOnMesh, cosBel, &
						   tracerRefine, NullVorticity,nullVort)
	if ( panelKind == QUAD_PANEL) then						   
		write(amrString,'(A,I1,A,I0.2)') 'quadAMR_',initNest,'to',initNest+refinementLimit
	endif
else
	call New(tracerRefine)
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
	
	call New(vtkOut,sphere,vtkFile,'tc1')
	call VTKOutput(vtkOut,sphere)
	frameCounter = frameCounter + 1
endif
!
!	intialize timestepping
!
call New(advRK4,sphere,numProcs)
timesteps = floor(tfinal,dt)
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
		call Delete(advRK4)
		call Delete(vtkOut)
		!
		! build new mesh
		!
		call LagrangianRemesh(sphere, NullVorticity, nullVort, tracerRefine, &
									  SetCosineBellOnMesh, cosBell, tracerRefine, &
									  tracerRefine)
		!
		!	create new associated objects
		!							  
		call New(advRK4,sphere,numProcs)
		call New(vtkOut,sphere,vtkFile,'tc1')	  
		
		write(logString,'(A,I4)') 'remesh ', remeshCounter
		if ( procRank == 0 ) call LogStats(sphere,exelog,trim(logString))
	endif
	!
	!	advance timestep
	!
	call AdvectionRK4(advRK4,sphere, dt, t, procRank, numProcs, LauritzenEtAlNonDivergentWind)
	
	t = real(timeJ+1,kreal)*dt
	
	!
	!	output timestep data
	!
	if ( procRank == 0 .AND. mod(timeJ+1,frameOUt) == 0 ) then
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

! TO DO : output final data requested by Test case 1

!
!	Clean up / free memory
!
call Delete(advRK4)
if ( procRank == 0 ) call Delete(vtkOut)
call Delete(tracerRefine)
call Delete(sphere)
call Delete(cosBell)
call Delete(exeLog)

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
			broadcastReals(4) = dt		
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

end subroutine

end program