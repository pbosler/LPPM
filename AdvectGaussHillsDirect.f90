program GaussianHillsAdvectionDirect

use NumberKindsModule
use OutputWriterModule
use LoggerModule
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
real(kreal), allocatable :: totalMassGHills(:), tracerVar(:)
real(kreal) :: sphereL2, sphereLinf, panelsLinf, particlesLinf, phiMax, phiMin, deltaPhi, phimax0, phimin0
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
character(len=MAX_STRING_LENGTH) :: namelistFile = 'AdvectGaussHillsDirect.namelist'
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

nTracer = 2
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
!call New(gHills, 0, 4)
!call InitCosineBellTracer(gHills, 0.0_kreal, 0.0_kreal, EARTH_RADIUS/3.0_kreal, 1000.0_kreal, tracerID)

!
! build initial mesh
!
call New(sphere, panelKind, initNest, AMR, nTracer, ADVECTION_SOLVER)
call SetGaussianHillsTracerOnMesh(sphere, gHills)
!call SetCosineBellTracerOnMesh(sphere, gHills)

!
! initialize remeshing and refinement
!
call ConvertFromRelativeTolerances(sphere, tracerMassTol, tracerVarTol, tracerID)
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'tracerMassTol = ', tracerMassTol )
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'tracerVarTol  = ', tracerVarTol )
call New(remesh, tracerID, tracerMassTol, tracerVarTol, amrLimit)
nullify(reference)
if ( AMR > 0 ) then
	call InitialRefinement(sphere, remesh, SetGaussianHillsTracerOnMesh, gHills, NullVorticity, nullvort)
	if ( panelKind == QUAD_PANEL ) &
		write(amrstring,'(A,I1,A,I0.2,A)') 'quadAMR_', initNest, 'to', initNest+amrLimit, '_'
	if ( panelKind == TRI_PANEL ) &
		write(amrstring,'(A,I1,A,I0.2,A)') 'triAMR_', initNest, 'to', initNest+amrLimit, '_'
else
	if ( panelKind == QUAD_PANEL ) &
		write(amrstring,'(A,I1,A)') 'quadUnif_', initNest, '_'
	if ( panelKind == TRI_PANEL ) &
		write(amrstring,'(A,I1,A)') 'triUnif_', initNest, '_'
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
mass0 = TotalMass(sphere, tracerID)
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

		call DirectRemesh(sphere, remesh)

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
	call AdvectionRK4Timestep(timekeeper, sphere, dt, t, procRank, numProcs, LauritzenEtAlNonDivergentWind)

	totalMassGHills(timeJ+1) = ( TotalMass(sphere, tracerID) - mass0 ) / mass0
	tracerVar(timeJ+1) = ( TracerVariance(sphere, tracerID) - var0 ) / var0

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

	!
	! calculate error : exact solution at final time should equal exact tracer from initial distribution
	!
	sphereParticles => sphere%particles
	spherePanels => sphere%panels
	do j = 1, sphereParticles%N
		sphereParticles%tracer(j,2) = abs( GHillsExact(sphereParticles%x(:,j)/EARTH_RADIUS, hmax, beta) - sphereParticles%tracer(j,1))
	enddo
	do j = 1, spherePanels%N
		if ( spherePanels%hasChildren(j) ) then
			spherePanels%tracer(j,2) = 0.0_kreal
		else
			spherePanels%tracer(j,2) = abs( GHillsExact(spherePanels%x(:,j)/EARTH_RADIUS, hmax, beta) - spherePanels%tracer(j,1) )
		endif
	enddo

	particlesLinf = maxval(sphereParticles%tracer(1:sphereParticles%N,2))  / maxval(sphereParticles%tracer(1:sphereParticles%N,1))
	panelsLinf = maxval( spherePanels%tracer(1:spherePanels%N,2) ) / maxval( spherePanels%tracer(1:spherePanels%N,1) )

	sphereLinf = max( particlesLinf, panelsLinf )
	sphereL2 = sum( spherePanels%tracer(1:spherePanels%N,2) * spherePanels%tracer(1:spherePanels%N,2) * spherePanels%area(1:spherePanels%N) )
	sphereL2 = sphereL2 / sum( spherePanels%tracer(1:spherePanels%N,1) * spherePanels%tracer(1:spherePanels%N,1) * spherePanels%area(1:spherePanels%N) )
	sphereL2 = sqrt(sphereL2)

	phimax = ( max( maxval(sphereParticles%tracer(1:sphereParticles%N,1)), maxval( spherePanels%tracer(1:spherePanels%N,1)) ) - phimax0) / deltaPhi
	phimin = ( min( minval(sphereParticles%tracer(1:sphereParticles%N,1)), minval( spherePanels%tracer(1:spherePanels%N,1)) ) - phimin0)/ deltaPhi
	if ( procRank == 0 ) then
		open( unit = WRITE_UNIT_1, file = datafile, status = 'REPLACE', action = 'WRITE', iostat = readwritestat)
		if ( readwritestat /= 0 ) then
			call LogMessage(exeLog, ERROR_LOGGING_LEVEL, 'data file ERROR : ', ' failed to open data file.')
		else
			write(WRITE_UNIT_1,'(A,F24.15,A)') 'passiveLinf = ', particlesLinf, ' ;'
			write(WRITE_UNIT_1,'(A,F24.15,A)') 'activeLinf = ', panelsLinf, ' ;'
			write(WRITE_UNIT_1,'(A,F24.15,A)') 'sphereLinf = ', sphereLinf, ' ;'
			write(WRITE_UNIT_1,'(A,F24.15,A)') 'sphereL2 = ', sphereL2, ' ;'
			write(WRITE_UNIT_1,'(A,F24.15,A)') 'phi_max = ', phimax, ' ;'
			write(WRITE_UNIT_1,'(A,F24.15,A)') 'phi_min = ', phimin, ' ;'
			write(WRITE_UNIT_1,'(A,F24.15,A)') 'dt_day = ', dt / ONE_DAY, ' ;'
			write(WRITE_UNIT_1,'(A,F24.15,A)') 'tfinal_day = ', tfinal / ONE_DAY, ' ;'
			write(WRITE_UNIT_1,'(A,F24.15,A)') 'mass = [ ', totalMassGHills(0), ' ; ...'
			do j = 1, timesteps-1
				write(WRITE_UNIT_1,'(F24.15,A)') totalMassGHills(j), ' ; ...'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') totalMassGHills(timesteps), ' ] ;'

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
deallocate(totalMassGHills)
deallocate(tracerVar)
call Delete(timekeeper)
call Delete(remesh)
if ( procrank == 0 ) call Delete(vtkOut)
call Delete(sphere)
call Delete(gHills)
call Delete(exeLog)

call MPI_FINALIZE(errCode)

contains

function GHillsExact(xyz, hmax, beta)
	real(kreal) :: GHillsExact
	real(kreal), intent(in) :: xyz(3), hmax, beta
	!
	real(kreal) :: xc1(3), xc2(3), h1, h2

	xC1 = [ cos(5.0_kreal * PI / 6.0_kreal), sin( 5.0_kreal * PI / 6.0_kreal ), 0.0_kreal ]
	xC2 = [ cos(7.0_kreal * PI / 6.0_kreal), sin( 7.0_kreal * PI / 6.0_kreal ), 0.0_kreal ]

	h1 = hmax * exp( -beta * ( sum( (xyz-xc1) * (xyz-xc1) ) ) )
	h2 = hmax * exp( -beta * ( sum( (xyz-xc2) * (xyz-xc2) ) ) )
	GHillsExact = h1 + h2
end function

subroutine ConvertFromRelativeTolerances(aMesh, tracerMassTol, tracerVarTol, tracerID)
	type(SphereMesh), intent(in) :: amesh
	real(kreal), intent(inout) :: tracerMassTol, tracerVarTol
	integer(kint), intent(in) :: tracerID
	tracerMassTol = tracerMassTol * MaximumTracerMass(aMesh, tracerID)
	tracerVarTol = tracerVarTol * MaximumTracerVariation(aMesh, tracerID)
end subroutine

subroutine ReadNamelistFile(rank)
	integer(kint), intent(in) :: rank
	integer(kint), parameter :: BCAST_INT_SIZE = 6, BCAST_REAL_SIZE= 4
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
