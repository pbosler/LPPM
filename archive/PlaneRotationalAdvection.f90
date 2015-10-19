program PlaneRotateAdvection

use NumberKindsModule
use LoggerModule
use ParticlesModule
use PanelsModule
use PlaneMeshModule
use PlaneOutputModule
use PlaneDirectSumModule
use PlaneTracerModule
use PlaneVorticityModule
use PlaneRemeshModule
use BIVARInterfaceModule

implicit none

include 'mpif.h'

!
! mesh variables
!
type(PlaneMesh) :: mesh
integer(kint) :: initNest, AMR, nTracer
real(kreal) :: xmin, xmax, ymin, ymax
integer(kint) :: boundaryType = FREE_BOUNDARIES

!
! tracer variables
!
type(TracerSetup) :: tracer3
integer(kint) :: tracerID

!
! remeshing / refinement variables
!
type(RemeshSetup) :: remesh
real(kreal) :: maxTol, tracerVarTol, lagVarTol
integer(kint) :: amrLimit, remeshInterval, remeshCounter, resetAlphaInterval
type(BIVARSetup) :: reference

!
! timestepping variables
!
type(PlaneRK4DirectSum) :: timekeeper
real(kreal) :: dt, tfinal, t
integer(kint) :: timeJ, timesteps

!
! computation variables
!
real(kreal), allocatable :: totalMass(:)

!
! input / output variables
!
type(PlaneOutput) :: meshOut
character(len=MAX_STRING_LENGTH) :: vtkFilename, vtkFileroot, outputDir, jobPrefix, summaryFile
character(len=128) :: namelistFile = 'RotateAdvection.namelist'

!
! mpi / computing environment variables
!
type(Logger) :: exeLog
character(len=MAX_STRING_LENGTH) :: logstring
integer(kint) :: errCode
real(kreal) :: wallclock
integer(kint), parameter :: BCAST_INT_SIZE = 5, BCAST_REAL_SIZE = 9
integer(kint) :: broadcastIntegers(BCAST_INT_SIZE)
real(kreal) :: broadcastReals(BCAST_REAL_SIZE)
integer(kint) :: j, iostat
type(VorticitySetup) :: nullVort

!
! namelists for user input
!
namelist /meshDefine/ initNest, AMR, amrLimit, xmin, xmax, ymin, ymax
namelist /timestepping/ dt, tfinal
namelist /remeshing/ maxTol, tracerVarTol, lagVarTol, remeshInterval, resetAlphaInterval
namelist /fileIO/ outputDir, jobPrefix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	INITIALIZE COMPUTER, MESH, TEST CASE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

call MPI_INIT(errCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numProcs,errCode)
call MPI_COMM_RANK(MPI_COMM_WORLD,procRank,errCode)

call New(exeLog,DEBUG_LOGGING_LEVEL)

wallclock = MPI_WTIME()

nTracer = 3
tracerID = 1

call InitLogger(exelog, procRank)

call ReadNamelistFile(procRank)

if ( procrank == 0 ) write(vtkFilename,'(A,I0.4,A)') trim(vtkFileRoot), 0, '.vtk'

!
! build initial mesh
!
call New(mesh, initNest, AMR, nTracer)
call InitializeRectangle(mesh, xmin, xmax, ymin, ymax, boundaryType)
write(logstring,'(A,I1,A,I8,A,F16.12)') 'nest = ', initNest, ', nPanels = ', mesh%panels%N_Active, ' total area = ', TotalArea(mesh)
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'base mesh : ',trim(logstring))

!
! initialize tracer
!
call InitThreeTracers(tracer3, tracerID)
call SetThreeTracersOnMesh(mesh, tracer3)

if ( AMR > 0 ) call InitialRefinement(mesh, remesh, SetThreeTracersOnMesh, tracer3, NullVorticity, nullvort)

!
! initialize output
!
if ( procrank == 0 ) then
	call New(meshOut, mesh, vtkFilename)
	call LogStats(mesh, exeLog)
	call OutputForVTK(meshout, mesh)
endif

!
! initialize time stepping
!
call New( timekeeper, mesh, numProcs)
timesteps = floor(tfinal / dt)
remeshCounter = 0
allocate(totalMass(0:timesteps))
totalMass = 0.0_kreal
totalMass(0) = GetTotalMass(mesh, tracerID)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	RUN THE PROBLEM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t = 0.0_kreal
do timeJ = 0, timesteps - 1
	!
	! remesh if necessary
	!
	if ( mod(timeJ+1, remeshInterval ) == 0 ) then
		remeshCounter = remeshCounter + 1
		!
		! choose appropriate remeshing procedure
		!
		if ( remeshcounter < resetAlphaInterval ) then
			!
			! remesh to t = 0
			!
			call LagrangianRemeshToInitialTime( mesh, remesh, NullVorticity, nullVort, SetThreeTracersOnMesh, tracer3)
			
		elseif ( remeshCounter == resetAlphaInterval ) then
			!
			! remesh to t = 0, then build a new reference mesh
			!
			call LagrangianRemeshToInitialTime( mesh, remesh, NullVorticity, nullVort, SetThreeTracersOnMesh, tracer3)
			call New(reference, mesh)
			call ResetLagrangianParameter(reference)
			call ResetLagrangianParameter(mesh)
		elseif ( remeshCounter > resetAlphaInterval .AND. mod(remeshCounter, resetAlphaInterval) == 0 ) then
			!
			! remesh to previous reference mesh, then create new reference mesh
			!
			call LagrangianRemeshToReferenceTime(mesh, reference, remesh)
			call Delete(reference)
			call New(reference, mesh)
			call ResetLagrangianParameter(reference)
			call ResetLagrangianParameter(mesh)
		else
			!
			! remesh to reference
			!
			call LagrangianRemeshToReferenceTime(mesh, reference, remesh)
		endif

		!
		! delete objects associated with old mesh
		!
		call Delete(timekeeper)
		if ( procrank == 0 ) call Delete(meshOut)
		!
		! create new mesh-associated objects
		!
		call New(timekeeper, mesh, numProcs)
		if ( procRank == 0 ) call New(meshout, mesh, vtkFilename)
	endif! remeshing

	!
	! increment time
	!
	call RK4TimestepIncompressibleAdvection( timekeeper, mesh, dt, t, procRank, numProcs, LeVeque93)
	t = real(timeJ+1,kreal) * dt
	!
	! compute conserved integrals
	!
	totalMass(timeJ+1) = GetTotalMass(mesh,tracerID)
	!
	! output timestep data
	!
	if ( procRank == 0 ) then
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 't = ', t)
		write(vtkFilename,'(A,I0.4,A)') trim(vtkFileRoot), timeJ+1, '.vtk'
		call UpdateFilename(meshOut, vtkFilename)
		call OutputForVTK(meshOUt, mesh)
	endif
enddo

deallocate(totalMass)
call Delete(timekeeper)
call Delete(tracer3)
if (procRank == 0 ) call Delete(meshOut)
call Delete(remesh)
call Delete(mesh)
call Delete(exeLog)

call MPI_FINALIZE(errCode)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	supplemental functions and subroutines
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine ConvertFromRelativeTolerances(aMesh, maxMassTol, tracerVarTol, lagVarTol)
	type(PlaneMesh), intent(in) :: aMesh
	real(kreal), intent(inout) :: maxMassTol, tracerVarTol, lagVarTol
end subroutine

subroutine ReadNamelistfile(rank)
!	processor 0 reads user input, then broadcasts initialization variables to all other processes.
!
	integer(kint), intent(in) :: rank
	integer(kint) :: readStat
	if ( rank == 0 ) then
		open(unit=READ_UNIT, file=namelistFile, status='OLD', action='READ',iostat=readStat)
		if ( readStat /= 0 ) stop 'cannot read namelist file.'

		read(READ_UNIT, nml = meshDefine)
		rewind(READ_UNIT)
		read(READ_UNIT, nml = timestepping)
		rewind(READ_UNIT)
		read(READ_UNIT, nml = fileIO)
		rewind(READ_UNIT)
		read(READ_UNIT, nml = remeshing)
		rewind(READ_UNIT)

		close(READ_UNIT)

		broadcastIntegers(1) = initNest
		broadcastIntegers(2) = AMR
		broadcastIntegers(3) = amrLimit
		broadcastIntegers(4) = remeshInterval
		broadcastIntegers(5) = resetAlphaInterval

		broadcastReals(1) = dt
		broadcastReals(2) = tfinal
		broadcastReals(3) = maxTol
		broadcastReals(4) = tracerVarTol
		broadcastReals(5) = lagVarTol
		broadcastReals(6) = xmin
		broadcastReals(7) = xmax
		broadcastReals(8) = ymin
		broadcastReals(9) = ymax

		write(vtkFileRoot,'(4A)') trim(outputDir), 'vtkOut/', trim(jobPrefix), '_'
		write(summaryFile,'(3A)') trim(outputDir), trim(jobPrefix), '_summary.txt'
	endif

	call MPI_BCAST(broadcastIntegers, BCAST_INT_SIZE, MPI_INTEGER, 0, MPI_COMM_WORLD, errCode)
	initNest = broadcastIntegers(1)
	AMR = broadcastIntegers(2)
	amrLimit = broadcastIntegers(3)
	remeshInterval = broadcastIntegers(4)
	resetAlphaInterval = broadcastIntegers(5)

	call MPI_BCAST(broadcastReals, BCAST_REAL_SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errCode)
	dt = broadcastReals(1)
	tfinal = broadcastReals(2)
	maxTol = broadcastReals(3)
	tracerVarTol = broadcastReals(4)
	lagVarTol = broadcastReals(5)
	xmin = broadcastReals(6)
	xmax = broadcastReals(7)
	ymin = broadcastReals(8)
	ymax = broadcastReals(9)
end subroutine

subroutine InitLogger(aLog, rank)
	type(Logger), intent(inout) :: aLog
	integer(kint), intent(in) :: rank
	if ( rank == 0 ) then
		call New(aLog,DEBUG_LOGGING_LEVEL)
	else
		call New(aLog,WARNING_LOGGING_LEVEL)
	endif
end subroutine

end program
