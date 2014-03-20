program NonlinearSteadySWE
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	This program solves the shallow water equations for Williamson et al (1992) test case 2.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
!	Williamson, D.L. et al, "A standard test set for numerical approximations to the shallow water equations in 
!		spherical geometry." J. Comp. Phys., 102:211-224, 1992.
!
!----------------

use NumberKindsModule
use SphereGeomModule
use OutputWriterModule
use LoggerModule
use ParticlesModule
use PanelsModule
use SphereMeshModule
use VTKOutputModule
use SWESetupModule
use TracerSetupModule
use SWEDirectSumModule

implicit none

include 'mpif.h'

!
! mesh variables
!
type(SphereMesh) :: sphere
integer(kint) :: panelKind, initNest, AMR, nTracer, problemKind
type(Particles), pointer :: sphereParticles => null()
type(Panels), pointer :: spherePanels => null()

!
! time stepping variables
!
type(SWERK4Data) :: sweRK4
real(kreal) :: t, dt, tfinal
integer(kint) :: timesteps, timeJ

!
! vorticity / divergence variables
!
type(SWESetup) :: testCase2
real(kreal) :: u0, rotAlpha, h0

!
! output variables
!
type(VTKSource) :: vtkOut
character(len=128) :: vtkRoot, vtkFile, outputDir, jobPrefix
character(len=128) :: dataFile, summaryFile
character(len=56) :: amrString
integer(kint) :: frameCounter, frameOut, readWriteStat
type(OutputWriter) :: writer

!
! user input
!
character(len=128) :: namelistInputFile = 'TestCase2.namelist'
namelist /sphereDefine/ panelKind, initNest, AMR
namelist /testCaseDefine/ u0, rotAlpha, h0
namelist /timestepping/ tfinal, dt
namelist /fileIO/ outputDir, jobPrefix, frameOut

!
! logging
!
type(Logger) :: exeLog
character(len=28) :: logKey
character(len=128) :: logString

!
! general and mpi variables
!
integer(kint) :: i, j, k, errCode
integer(kint), parameter :: BROADCAST_INT_SIZE = 3, BROADCAST_REAL_SIZE = 5
integer(kint) :: broadcastIntegers(BROADCAST_INT_SIZE)
real(kreal) :: broadcastReals(BROADCAST_REAL_SIZE)
real(kreal) :: wallClock

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	PROGRAM START
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
call MPI_INIT(errCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numProcs,errCode)
call MPI_COMM_RANK(MPI_COMM_WORLD,procRank,errCode)

call InitLogger(exeLog,procRank)

wallClock = MPI_WTIME()

!
!	read user input from namelist file, set starting state
!
call ReadNamelistInputFile(procRank)
call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logkey,' initial broadcast done.')

!
! define test case
!
call New(testCase2,TC2_NINT, TC2_NREAL)
call InitNonlinearSteadyStateTC2(testCase2, u0, rotAlpha, h0)

!
! initialize sphere
!
nTracer = 1
problemKind = SWE_SOLVER
call New(sphere,panelKind,initNest,AMR,nTracer,problemKind)
call SetFlowMapLatitudeTracerOnMesh(sphere,1)
call SetNonlinearSteadyStateTC2OnMesh(sphere,testCase2)
call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,trim(logKey),' base sphere (t = 0) done.')

if ( AMR > 0 ) then
! TO DO : AMR
else
	if ( procRank == 0 .AND. panelKind == QUAD_PANEL ) then
		write(amrString, '(A,I1,A)') 'quadUnif',initNest,'_'
	endif
endif	

!
! initialize output, output t=0 data
!
frameCounter = 0
if ( procRank == 0 ) then
	call LogStats(sphere,exeLog,'initial mesh :')
	write(vtkRoot,'(A,A,A,A)') trim(outputDir),'/vtkOut/',trim(jobPrefix),trim(amrString)
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),frameCounter,'.vtk'
	
	call New(vtkOut,sphere,vtkFile,'SWE test case 2')
	call VTKOutput(vtkOut,sphere)
	frameCounter = frameCounter + 1
endif

!
! initialize timestepping
!
call New(sweRK4,sphere,numProcs)
timesteps = floor(tfinal/dt)
t = 0.0_kreal

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	RUN THE PROBLEM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
do timeJ = 0, timesteps-1
	! TO DO : remeshing for SWE
	
	!
	! increment time
	!
	call SWERK4Timestep(sweRK4, sphere, dt, procRank, numProcs)
	
	!
	! TO DO : error calculation
	!
	
	t = real(timeJ+1,kreal)*dt
	
	!
	! output timestep data
	!
	if ( procRank == 0 .AND. mod(timeJ+1,frameOut) == 0 ) then
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,'day = ',t/ONE_DAY)
		
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
		call UpdateFileName(vtkOut,vtkFile)
		call VTKOutput(vtkOut,sphere)
		
		frameCounter = frameCounter + 1
	endif
enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	Clear memory and Finalize
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
call Delete(sweRK4)
if (procRank == 0 ) call Delete(vtkOut)
call Delete(sphere)
call Delete(testCase2)
call Delete(exeLog)
call MPI_FINALIZE(ErrCode)

contains

subroutine ReadNamelistInputFile(rank)
	integer(kint), intent(in) :: rank
	if ( rank == 0 ) then
		open(unit = READ_UNIT, file=namelistInputFile, status = 'OLD', action = 'READ', iostat = readWriteStat)
		if ( readWriteStat /= 0 ) then
			stop 'cannot read namelist file'
		endif
		read(READ_UNIT,nml=sphereDefine)
		rewind(READ_UNIT)
		read(READ_UNIT,nml=testCaseDefine)
		rewind(READ_UNIT)
		read(READ_UNIT,nml=timestepping)
		rewind(READ_UNIT)
		read(READ_UNIT,nml=fileIO)
		
		broadcastIntegers(1) = panelKind
		broadcastIntegers(2) = initNest
		broadcastIntegers(3) = AMR
		
		broadcastReals(1) = u0
		broadcastReals(2) = rotAlpha
		broadcastReals(3) = h0
		broadcastReals(4) = tfinal
		broadcastReals(5) = dt
	endif
	
	call MPI_BCAST(broadcastIntegers, BROADCAST_INT_SIZE, MPI_INTEGER, 0, MPI_COMM_WORLD, errCode)
	panelKind = broadcastIntegers(1)
	initNest = broadcastIntegers(2)
	AMR = broadcastIntegers(3)
	
	call MPI_BCAST(broadcastReals, BROADCAST_REAL_SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errCode)
	u0 = broadcastReals(1)
	rotAlpha = broadcastReals(2)
	h0 = broadcastReals(3)
	tfinal = broadcastReals(4)*ONE_DAY ! convert tfinal from days to seconds
	dt = broadcastReals(5)*ONE_DAY ! convert dt from days to seconds
end subroutine

subroutine InitLogger(aLog, rank)
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	if (rank == 0 ) then
		call New(aLog,DEBUG_LOGGING_LEVEL)
	else
		call New(aLog,DEBUG_LOGGING_LEVEL)
	endif
	write(logkey,'(A,I0.2,A)') 'EXE_LOG',rank,' : '
end subroutine


end program
