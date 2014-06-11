program GaussianHillsAdvection

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use SphereMeshModule
use AdvectionModule
use SphereMeshModule
use TracerSetupModule
use VTKOutputModule
use BVESetupModule

implicit none

include 'mpif.h'

!
! mesh variables
!
type(SphereMesh) :: sphere
integer(kint) :: panelKind, initNest, AMR, nTracer

!
! tracer variables
!
type(TracerSetup) :: gHills
integer(kint) :: tracerID
real(kreal) :: hmax, beta

!
! vorticity placeholder
!
type(BVESetup) :: nullVort

!
! remeshing / refinement variables
!
type(RemeshSetup) :: remesh
integer(kint) :: remeshInterval, resetAlphaInterval, amrLimit, remeshCounter
real(kreal) :: tracerMassTol, tracerVarTol
type(ReferenceSphere) :: reference

!
! time stepping variables
!
type(AdvRK4Data) :: timekeeper
real(kreal) :: t, tfinal, dt
integer(kint) :: timesteps, timeJ

!
! output variables
!
type(VTKSource) :: vtkOut
character(len = MAX_STRING_LENGTH) :: vtkRoot, vtkFile, outputDir, jobPrefix, dataFile, summaryFile
character(len = 56) :: amrString
integer(kint) :: frameCounter, frameOut, readWriteStat
type(OutputWriter) :: writer

!
! test case variables
!
real(kreal), allocatable(:) :: totalMassGHills(:)

!
! logging
!
type(Logger) :: exeLog
character(len=28) :: logkey
character(len=MAX_STRING_LENGTH) :: logstring

!
! mpi / computing environment / general variables
!
integer(kint) :: errCode
real(kreal) :: wallclock
integer(kint) :: j

!
! namelists and user input
!
character(len=MAX_STRING_LENGTH) :: namelistFile = 'AdvectGaussianHills.namelist'
namelist /meshDefine/ initNest, AMR, panelKind, amrLimit, tracerMassTol, tracerVarTol 
namelist /timestepping/ tfinal, dt, remeshInterval, resetAlphaInterval
namelist /fileIO/ outputDir, jobPrefix, frameOut


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	INITIALIZE COMPUTER, MESH, TEST CASE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

call MPI_INIT(errCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, errCode)
call MPI_COMM_RANK(MPI_COMM_WORLD, procRank, errCode)

call InitLogger(exeLog, procRank)

wallclock = MPI_WTIME()

nTracer = 1
tracerID = 1
hmax = 0.95_kreal
beta = 5.0_kreal

!
! get user input
!
call ReadNamelistFile(procRank)

!
! define tracer
!
call New(gHills, GAUSS_HILLS_N_INT, GAUSS_HILLS_N_REAL)
call InitGaussianHillsTracer(gHills, hmax, beta, tracerID)

!
! build initial mesh
!
call New(sphere, panelKind, initNest, AMR, nTracer, ADVECTION_SOLVER)
call SetGaussianHillsTracerOnMesh(sphere, gHills)

!
! initialize remeshing and refinement
!
call ConvertFromRelativeTolerances(sphere, tracerMassTol, tracerVarTol)
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'tracerMassTol = ', tracerMassTol )
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'tracerVarTol  = ', tracerVarTol )
call New(remesh, tracerID, tracerMassTol, tracerVarTol, amrLimit)
if ( AMR > 0 ) then
	call InitialRefinement(sphere, remesh, SetGaussianHillsTracerOnMesh, gHills, NullVorticity, nullvort)
	if ( panelKind == QUAD_PANEL ) &
		write(amrstring,'(A,I1,A,I0.2,A)') 'quadAMR_', initNest, 'to', initNest+amrLimit, '_'
else
	if ( panelKind == QUAD_PANEL ) &
		write(amrstring,'(A,I1,A)') 'quadUnif_', initNest, '_'
endif

!
! initialize output
!
if ( procrank == 0 ) then
	write(vtkRoot,'(A,A,A,A,A)') trim(outputDir), '/vtkOut/',trim(jobPrefix),trim(amrString),'_'
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),0,'.vtk'
	write(summaryFile,'(A,A,A,A)') trim(outputDir), trim(jobPrefix), trim(amrString), '_summary.txt'
	write(datafile,'(A,A,A,A)') trim(outputDir), trim(jobPrefix), trim(amrstring), '_calculatedData.txt'
	call New(vtkOut, sphere, vtkFile, 'Gaussian hills advection')
	call VTKOutput(vtkOut, sphere)
endif

!
! initialize time stepping
!
call New(timekeeper, sphere, numProcs)
timesteps = floor(tfinal / dt)
t = 0.0_kreal
remeshCounter = 0
frameCounter = 1
allocate(totalMassGHills(0:timesteps))
totalMassGHills = 0.0_kreal
totalMassGHills(0) = TotalMass(sphere, tracerID)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	RUN THE PROBLEM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do timeJ = 0, timesteps - 1
	if ( mod( timeJ+1, remeshInterval) == 0 ) then
		!
		! remesh before timestep
		!
		remeshCounter = remeshCounter + 1
		!
		! choose appropriate remeshing procedure
		!
		if ( remeshCounter < resetAlphaInterval ) then
			!
			! remesh to t = 0
			!
			call LagrangianRemeshToInitialTime(sphere, remesh, NullVorticity, nullVort, SetGaussianHillsTracerOnMesh, gHills)
			
		elseif ( remeshCounter == resetAlphaInterval ) then
			!
			! remesh to t = 0, create reference mesh to current time
			!
			call LagrangianRemeshToInitialTime(sphere, remesh, NullVorticity, nullVort, SetGaussianHillsTracerOnMesh, gHills)
			call New(reference, sphere)
			call ResetLagrangianParameter(sphere)
			
		elseif ( remeshCounter > resetAlphaInterval .AND. mod(remeshCounter, resetAlphaInterval) == 0 ) then
			!
			! remesh to existing reference, then create new reference to current time
			!
			call LagrangianRemeshToReference( sphere, reference, remesh)
			call Delete(reference)
			call New( reference, sphere)
			call ResetLagrangianParameter(sphere)
		else
			!
			! remesh to existing reference
			!
			call LagrangianRemeshToReference(sphere, reference, remesh)
			
		endif
		!
		! delete objects associated with old mesh
		!
		call Delete(timekeeper)
		if ( procrank == 0 ) call Delete(vtkOUt)
		!
		! create new associated objects for new mesh
		!
		call New(timekeeper, sphere, numProcs)
		if ( procRank == 0 ) call New(vtkOut, sphere, vtkFile, 'Gaussian hills advection')
	endif ! remesh
	
	!
	! advance time
	!
	call BVERK4Timestep(timekeeper, sphere, dt, procRank, numProcs)
	
	totalMassGHills(timeJ+1) = TotalMass(sphere, tracerID)
	
	t = real( timeJ+1, kreal) * dt
	
	if ( procRank == 0 .AND. mod( timeJ+1, frameOut) == 0 ) then
		call LogMessage(exelog, TRACE_LOGGING_LEVEL, 'day = ', t/ONE_DAY)
		
		write(vtkFile, '(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
		call UpdateFilename(vtkOut, vtkFile)
		call VTKOutput(vtkOut, sphere)
		
		frameCounter = frameCounter + 1
	endif	
enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	OUTPUT FINAL DATA
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	FREE MEMORY, CLEAN UP, FINALIZE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
call Delete(reference)
deallocate(totalMassGHills)
call Delete(timekeeper)
call Delete(remesh)
if ( procrank == 0 ) call Delete(vtkOut)
call Delete(sphere)
call Delete(gHills)
call Delete(exeLog)

call MPI_FINALIZE(errCode)

contains

subroutine ConvertFromRelativeTolerances(aMesh, tracerMassTol, tracerVarTol, tracerID)
	type(SphereMesh), intent(in) :: sphere
	real(kreal), intent(inout) :: tracerMassTol, tracerVarTol
	integer(kint), intent(in) :: tracerID
	tracerMassTol = tracerMassTol * MaximumTracerMass(aMesh, tracerID)
	tracerVarTol = tracerVarTol * MaximumTracerVariation(aMesh, tracerID)
end subroutine

subroutine ReadNamelistFile(rank)
	integer(kint), intent(in) :: rank
	integer(kint), parameter :: BCAST_INT_SIZE = 6, BCASE_REAL_SIZE= 4
	integer(kint) :: broadcastIntegers(BCAST_INT_SIZE)
	real(kreal) :: broadcastReals(BCAST_REAL_SIZE)
	
	if ( rank == 0 ) then
		open(unit=READ_UNIT, file=namelistfile, status='OLD', action='READ', iostat=readWriteStat)
			if ( readWriteStat /= 0 ) stop 'cannot read namelist file.'
			read(READ_UNIT, nml=meshDefine)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=timestepping)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=fileIO)
			rewind(READ_UNIT)
		close(READ_UNIT)
		
		broadcastIntegers(1) = panelKind
		broadcastIntegers(2) = initNest
		broadcastIntegers(3) = AMR
		broadcastIntegers(4) = amrLimit
		broadcastIntegers(5) = remeshInterval
		broadcastIntegers(6) = resetAlphaInterval
		
		broadcastReals(1) = tracerMassTol
		broadcastReals(2) = tracerVarTol
		broadcastReals(3) = dt
		broadcastReals(4) = tfinal
	endif

	call MPI_BCAST(broadcastIntegers, BCAST_INT_SIZE, MPI_INTEGER, 0, MPI_COMM_WORLD, errCode)
	panelKind = broadcastIntegers(1)
	initNest = broadcastIntegers(2)
	AMR = broadcastIntegers(3)
	amrLimit = broadcastIntegers(4)
	remeshInterval = broadcastIntegers(5)
	resetAlphaInterval = broadcastIntegers(6)
	
	call MPI_BCAST(broadcastReals, BCAST_REAL_SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errCode)
	tracerMassTol = broadcastReals(1)
	tracerVarTol = broadcastReals(2)
	dt = broadcastReals(3) * ONE_DAY		! convert time to seconds
	tfinal = broadcastReals(4) * ONE_DAY	! convert time to seconds
end subroutine

subroutine InitLogger(alog,rank)
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	if ( rank == 0 ) then
		call New(aLog,DEBUG_LOGGING_LEVEL)
	else
		call New(aLog,WARNING_LOGGING_LEVEL)
	endif
    write(logKey,'(A,I0.2,A)') 'EXE_LOG',rank,' : '
end subroutine

end program