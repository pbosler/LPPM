program MovingVorticesAdvection

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use SphereGeomModule
use SphereMeshModule
use AdvectionModule
use ParticlesModule
use PanelsModule
use SphereMeshModule
use TracerSetupModule
use VTKOutputModule
use BVESetupModule
use SphereRemeshModule

implicit none

include 'mpif.h'

!
! mesh variables
!
type(SphereMesh) :: sphere
integer(kint) :: panelKind, initNest, AMR, nTracer
type(Particles), pointer :: sphereParticles
type(Panels), pointer :: spherePanels

!
! tracer variables
!
type(TracerSetup) :: testCaseTracer
integer(kint) :: tracerID
real(kreal) :: rho0, gamma

!
! vorticity placeholder
!
type(BVESetup) :: nullVort

!
! remeshing / refinement variables
!
type(RemeshSetup) :: remesh
integer(kint) :: remeshInterval, resetAlphaInterval, amrLimit, remeshCounter
real(kreal) :: tracerMassTol, tracerVarTol, lagVarTol
type(ReferenceSphere), pointer :: reference

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
real(kreal), allocatable :: totalMasstestCaseTracer(:), sphereL2(:), sphereLinf(:), panelsLinf(:), particlesLinf(:), phiMax(:), phiMin(:), tracerVar(:)
real(kreal) :: deltaPhi, phimax0, phimin0
real(kreal) :: mass0, var0

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
character(len=MAX_STRING_LENGTH) :: namelistFile = 'MovingVortices.namelist'
namelist /meshDefine/ initNest, AMR, panelKind, amrLimit, tracerMassTol, tracerVarTol, lagVarTol
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

nTracer = 2
tracerID = 1
rho0 = 3.0_kreal
gamma = 5.0_kreal

!
! get user input
!
call ReadNamelistFile(procRank)

!
! define tracer
!
call New(testCaseTracer, 1, 1)
testCaseTracer%reals(1) = 0.0_kreal
testCaseTracer%integers(1) = tracerID

!
! build initial mesh
!
call New(sphere, panelKind, initNest, AMR, nTracer, ADVECTION_SOLVER)

call LogStats( sphere, exeLog, 'DEBUG')

call SetMovingVortsTracerOnMesh(sphere, testCaseTracer)

!
! initialize remeshing and refinement
!
call ConvertFromRelativeTolerances(sphere, tracerMassTol, tracerVarTol, tracerID, lagVarTol)
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'tracerMassTol = ', tracerMassTol )
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'tracerVarTol  = ', tracerVarTol )
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, '   lagVarTol  = ', lagVarTol )
call New(remesh, tracerID, tracerMassTol, tracerVarTol, lagVarTol, amrLimit)
nullify(reference)
if ( AMR > 0 ) then
	call InitialRefinement(sphere, remesh, SetMovingVortsTracerOnMesh, testCaseTracer, NullVorticity, nullvort)
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

	call LogStats( sphere, exeLog)

	write(vtkRoot,'(A,A,A,A,A)') trim(outputDir), '/vtkOut/',trim(jobPrefix),trim(amrString),'_'
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),0,'.vtk'
	write(summaryFile,'(A,A,A,A)') trim(outputDir), trim(jobPrefix), trim(amrString), '_summary.txt'
	write(datafile,'(A,A,A,A)') trim(outputDir), trim(jobPrefix), trim(amrstring), '_calculatedData.m'
	call New(vtkOut, sphere, vtkFile, 'moving vortices')
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
allocate(totalMasstestCaseTracer(0:timesteps))
totalMasstestCaseTracer = 0.0_kreal
mass0 = TotalMass(sphere, tracerID)
allocate(sphereL2(0:timesteps))
sphereL2 = 0.0_kreal
allocate(sphereLinf(0:timesteps))
sphereLinf = 0.0_kreal
allocate(particlesLinf(0:timesteps))
particlesLinf = 0.0_kreal
allocate(panelsLinf(0:timesteps))
panelsLinf = 0.0_kreal
allocate(phiMax(0:timesteps))
phiMax = 0.0_kreal
allocate(phiMin(0:timesteps))
phiMin = 0.0_kreal
allocate(tracerVar(0:timesteps))
tracerVar = 0.0_kreal
var0 = TracerVariance(sphere, tracerID)

sphereParticles => sphere%particles
spherePanels => sphere%panels
phimax0 = max( maxval(sphereParticles%tracer(1:sphereParticles%N,1)), maxval(spherePanels%tracer(1:spherePanels%N,1)) )
phimin0 = 0.0_kreal
deltaPhi = phimax0 - phimin0


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
			call LagrangianRemeshToInitialTime(sphere, remesh, NullVorticity, nullVort, SetMovingVortsTracerOnMesh, testCaseTracer)

		elseif ( remeshCounter == resetAlphaInterval ) then
			!
			! remesh to t = 0, create reference mesh to current time
			!
			call LagrangianRemeshToInitialTime(sphere, remesh, NullVorticity, nullVort, SetMovingVortsTracerOnMesh, testCaseTracer)
			allocate(reference)
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
		sphereParticles => sphere%particles
		spherePanels => sphere%panels
	endif ! remesh

	!
	! advance time
	!
	call AdvectionRK4Timestep(timekeeper, sphere, dt, t, procRank, numProcs, MovingVorticesVelocity)
	t = real( timeJ+1, kreal) * dt
	testCaseTracer%reals(1) = t

	!
	! calculate error
	!
	do j = 1, sphereParticles%N
		sphereParticles%tracer(j,2) = sphereParticles%tracer(j,1) - testCaseTracerExact(sphereParticles%x(:,j), t, rho0, gamma)
	enddo
	do j = 1, spherePanels%N
		if ( spherePanels%hasChildren(j) ) then
			spherePanels%tracer(j,2) = 0.0_kreal
		else
			spherePanels%tracer(j,2) = spherePanels%tracer(j,1) - testCaseTracerExact( spherePanels%x(:,j), t, rho0, gamma)
		endif
	enddo
	totalMasstestCaseTracer(timeJ+1) = ( TotalMass(sphere, tracerID) - mass0 ) / mass0
	tracerVar(timeJ+1) = ( TracerVariance(sphere, tracerID) - var0 ) / var0
	
	particlesLinf(timeJ+1) = maxval(sphereParticles%tracer(1:sphereParticles%N,2))  / maxval(sphereParticles%tracer(1:sphereParticles%N,1))
	panelsLinf(timeJ+1) = maxval( spherePanels%tracer(1:spherePanels%N,2) ) / maxval( spherePanels%tracer(1:spherePanels%N,1) )

	sphereLinf(timeJ+1) = max( particlesLinf(timeJ+1), panelsLinf(timeJ+1) )
	sphereL2(timeJ+1) = sum( spherePanels%tracer(1:spherePanels%N,2) * spherePanels%tracer(1:spherePanels%N,2) * spherePanels%area(1:spherePanels%N) )
	sphereL2(timeJ+1) = sphereL2(timeJ+1) / sum( spherePanels%tracer(1:spherePanels%N,1) * spherePanels%tracer(1:spherePanels%N,1) * spherePanels%area(1:spherePanels%N) )
	sphereL2(timeJ+1) = sqrt(sphereL2(timeJ+1))

	phimax(timeJ+1) = ( max( maxval(sphereParticles%tracer(1:sphereParticles%N,1)), maxval( spherePanels%tracer(1:spherePanels%N,1)) ) - phimax0) / deltaPhi
	phimin(timeJ+1) = ( min( minval(sphereParticles%tracer(1:sphereParticles%N,1)), minval( spherePanels%tracer(1:spherePanels%N,1)) ) - phimin0)/ deltaPhi

	!
	! output data
	!
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

	if ( procRank == 0 ) then
		open( unit = WRITE_UNIT_1, file = datafile, status = 'REPLACE', action = 'WRITE', iostat = readwritestat)
		if ( readwritestat /= 0 ) then
			call LogMessage(exeLog, ERROR_LOGGING_LEVEL, 'data file ERROR : ', ' failed to open data file.')
		else
			write(WRITE_UNIT_1,'(A,F24.15,A)') 'passiveLinf = [', particlesLinf(0), ' ;'
			do j = 1, timesteps -1
				write(WRITE_UNIT_1,'(F24.15,A)') particlesLinf(j), ' ;'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') particlesLinf(timesteps), ' ];'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'activeLinf = [', panelsLinf(0), ' ;'
			do j = 1, timesteps -1
				write(WRITE_UNIT_1,'(F24.15,A)') panelsLinf(j), ' ;'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') panelsLinf(timesteps), ' ];'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'sphereLinf = [', sphereLinf(0), ' ;'
			do j = 1, timesteps -1
				write(WRITE_UNIT_1,'(F24.15,A)') sphereLinf(j), ' ;'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') sphereLinf(timesteps), ' ];'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'sphereL2 = [', sphereL2(0), ' ;'
			do j = 1, timesteps -1
				write(WRITE_UNIT_1,'(F24.15,A)') sphereL2(j), ' ;'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') sphereL2(timesteps), ' ];'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'phi_max = [', phimax(0), ' ;'
			do j = 1, timesteps -1
				write(WRITE_UNIT_1,'(F24.15,A)') phimax(j), ' ;'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') phimax(timesteps), ' ];'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'phi_min = [', phimin(0), ' ;'
			do j = 1, timesteps -1
				write(WRITE_UNIT_1,'(F24.15,A)') phimin(j), ' ;'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') phimin(timesteps), ' ];'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'dt_day = ', dt / ONE_DAY, ' ;'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'tfinal_day = ', tfinal / ONE_DAY, ' ;'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'mass = [ ', totalMasstestCaseTracer(0), ' ; ...'
			do j = 1, timesteps-1
				write(WRITE_UNIT_1,'(F24.15,A)') totalMasstestCaseTracer(j), ' ; ...'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') totalMasstestCaseTracer(timesteps), ' ] ;'
		
			write(WRITE_UNIT_1,'(A,F24.15,A)') 'tracerVar = [ ', tracerVar(0), ' ; ...'
			do j = 1, timesteps-1
				write(WRITE_UNIT_1,'(F24.15,A)') tracerVar(j), ' ; ...'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') tracerVar(timesteps), ' ] ;'
		endif
		close(WRITE_UNIT_1)

		write(logstring,'(A, F8.2,A)') 'elapsed time = ', (MPI_WTIME() - wallClock)/60.0, ' minutes.'
		call LogMessage(exelog,TRACE_LOGGING_LEVEL,'PROGRAM COMPLETE : ',trim(logstring))

	endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	FREE MEMORY, CLEAN UP, FINALIZE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (associated(reference)) then
	call Delete(reference)
	deallocate(reference)
endif
deallocate(totalMasstestCaseTracer)
deallocate(tracerVar)
deallocate(sphereL2)
deallocate(sphereLinf)
deallocate(particlesLinf)
deallocate(panelsLinf)
deallocate(phiMax)
deallocate(phiMin)
call Delete(timekeeper)
call Delete(remesh)
if ( procrank == 0 ) call Delete(vtkOut)
call Delete(sphere)
call Delete(testCaseTracer)
call Delete(exeLog)

call MPI_FINALIZE(errCode)

contains

function testCaseTracerExact(xyz, t, rho0, gamma)
	real(kreal) :: testCaseTracerExact
	real(kreal), intent(in) :: xyz(3), rho0, gamma, t
	!
	real(kreal) :: wr, rho, u0

	rho = rho0 * cos( Latitude(xyz) )
	u0 = 2.0_kreal * PI * EARTH_RADIUS / (12.0_kreal * ONE_DAY)
	wr = u0 * 0.5_kreal * 3.0_kreal * sqrt(3.0_kreal) * tanh(rho) / cosh(rho) / cosh(rho) * ( rho / (rho*rho + ZERO_TOL*ZERO_TOL) )

	testCaseTracerExact = 1.0_kreal - tanh(rho* sin(Longitude(xyz) - wr * t /EARTH_RADIUS) / gamma )
end function

subroutine ConvertFromRelativeTolerances(aMesh, tracerMassTol, tracerVarTol, tracerID, lagVarTol)
	type(SphereMesh), intent(in) :: amesh
	real(kreal), intent(inout) :: tracerMassTol, tracerVarTol, lagVarTol
	integer(kint), intent(in) :: tracerID
	tracerMassTol = tracerMassTol * MaximumTracerMass(aMesh, tracerID)
	tracerVarTol = tracerVarTol * MaximumTracerVariation(aMesh, tracerID)
	lagVarTol = lagVarTol * MaximumLagrangianParameterVariation(aMesh)
end subroutine

subroutine ReadNamelistFile(rank)
	integer(kint), intent(in) :: rank
	integer(kint), parameter :: BCAST_INT_SIZE = 6, BCAST_REAL_SIZE= 5
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
		broadcastReals(5) = lagVarTol
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
	lagVarTol = broadcastReals(5)
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
