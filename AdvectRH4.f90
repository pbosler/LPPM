program rh4advection

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use SphereMeshModule
use AdvectionModule
use ParticlesModule
use PanelsModule
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
type(TracerSetup) :: nullScalar
integer(kint) :: tracerID

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
type(VTKSource) :: vtkOut, meshOut
character(len = MAX_STRING_LENGTH) :: vtkRoot, vtkFile, vtkMeshFile, outputDir, jobPrefix, dataFile, summaryFile
character(len = 56) :: amrString
integer(kint) :: frameCounter, frameOut, readWriteStat
type(OutputWriter) :: writer

!
! test case variables
!
real(kreal), allocatable :: totalTracer(:), tracerVar(:)

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
character(len=MAX_STRING_LENGTH) :: namelistFile = 'AdvectRH4.namelist'
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

!
! get user input
!
call ReadNamelistFile(procRank)

!
! build initial mesh
!
call New(sphere, panelKind, initNest, AMR, nTracer, ADVECTION_SOLVER)
call SetFlowMapLatitudeTracerOnMesh(sphere,1)
call SetFlowMapLatitudeTracerOnMesh(sphere,2)

!
! initialize remeshing and refinement
!
call ConvertFromRelativeTolerances(sphere, tracerMassTol, tracerVarTol, tracerID)
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'tracerMassTol = ', tracerMassTol )
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'tracerVarTol  = ', tracerVarTol )
call New(remesh, tracerID, tracerMassTol, tracerVarTol, AMR)
nullify(reference)
if ( AMR > 0 ) then
	!call InitialRefinement(sphere, remesh, SetGaussianHillsTracerOnMesh, gHills, NullVorticity, nullvort)
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

sphereParticles => sphere%particles
spherePanels => sphere%panels

!
! initialize output
!
if ( procrank == 0 ) then

	call LogStats( sphere, exeLog)

	write(vtkRoot,'(A,A,A,A,A)') trim(outputDir), '/vtkOut/',trim(jobPrefix),trim(amrString),'_'
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),0,'.vtk'
	
	write(vtkMeshFile,'(A,A,I0.4,A)') trim(vtkRoot),'_mesh_',0,'.vtk'
	
	write(summaryFile,'(A,A,A,A)') trim(outputDir), trim(jobPrefix), trim(amrString), '_summary.txt'
	write(datafile,'(A,A,A,A)') trim(outputDir), trim(jobPrefix), trim(amrstring), '_calculatedData.m'

	call New(vtkOut, sphere, vtkFile, 'rh4 advection')
	call VTKOutput(vtkOut, sphere)
	
	call New(meshOut, sphere, vtkMeshFile, 'rh4 advection')
	call VTKOutputMidpointRule(meshOUt,sphere)
	
endif

!
! initialize time stepping
!
call New(timekeeper, sphere, numProcs)
timesteps = floor(tfinal / dt)
t = 0.0_kreal
remeshCounter = 0
frameCounter = 1
allocate(totalTracer(0:timesteps))
totalTracer = 0.0_kreal
totalTracer(0) = TotalMass(sphere,1)
allocate(tracerVar(0:timesteps))
tracerVar = 0.0_kreal
tracerVar(0) = TracerVariance(sphere,1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	RUN THE PROBLEM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do timeJ = 0, timesteps-1
	if ( mod(timeJ+1, remeshInterval) == 0 ) then
		!
		!	remesh before timestep
		!
		remeshCounter = remeshCounter + 1
		!
		!	choose appropriate remeshing procedure
		!
		if ( remeshCounter < resetAlphaInterval ) then
			!
			!	remesh to t = 0
			!
			call LagrangianRemeshToInitialTime( sphere, remesh, NullVorticity, nullVort, nullTracer, nullScalar)
			call SetFlowMapLatitudeTracerOnMesh(sphere, 1)
			call SetFlowMapLatitudeTracerOnMesh(sphere, 2)
		elseif ( remeshCounter == resetAlphaInterval ) then
			!
			!	remesh to t = 0, then create reference mesh to current time
			!
			call LagrangianRemeshToInitialTime( sphere, remesh, NullVorticity, nullVort, nullTracer, nullScalar)
			call SetFlowMapLatitudeTracerOnMesh(sphere, 1)
			call SetFlowMapLatitudeTracerOnMesh(sphere, 2)
			allocate(reference)
			call New(reference, sphere)
			call ResetLagrangianParameter(sphere)
		elseif ( remeshCounter > resetAlphaInterval .AND. mod(remeshCounter, resetAlphaInterval) == 0 ) then
			!
			!	remesh to existing reference time, then create new reference mesh to current time
			!
			call LagrangianRemeshToReference( sphere, reference, remesh )
			call Delete(reference)
			call New( reference, sphere )
			call SetFlowMapLatitudeTracerOnMesh(sphere, 2)
			call ResetLagrangianParameter(sphere)
		else
			!
			!	remesh to existing reference time
			!
			call LagrangianRemeshToReference(sphere, reference, remesh)
		endif
		!
		!	delete objects associated with old mesh, create new ones for new mesh
		!
		sphereParticles => sphere%particles
		spherePanels => sphere%panels
		call Delete(timekeeper)
		call New( timekeeper, sphere, numProcs )
		if ( procRank == 0 ) then
			call Delete(vtkOut)
			call Delete(meshOut)
			call New(vtkOut, sphere, vtkFile, 'rh4 advection')
			call New(meshOut, sphere, vtkMeshFile, 'rh4 advection')
		endif
	endif ! remesh
	
	!
	!	advance time
	!
	call AdvectionRK4Timestep( timekeeper, sphere, dt, t, procRank, numProcs, rh4velocity)
	
	totalTracer(timeJ+1) = TotalMass(sphere, 1)
	tracerVar(timeJ+1) = TracerVariance(sphere, 1)
	
	t = real( timeJ+1, kreal) * dt
	
	if ( procRank == 0 .AND. mod(timeJ+1,frameOut) == 0 ) then
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'day = ', t/ONE_DAY )
		
		write(vtkFile, '(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
		write(vtkMeshFile, '(A,A,I0.4,A)') trim(vtkRoot), '_mesh_', frameCounter, '.vtk'
		
		call UpdateFilename(vtkOut, vtkFile)
		call UpdateFilename(meshOut, vtkMeshFile)
		
		call VTKOutput(vtkOut, sphere)
		call VTKOutputMidpointRule(meshOut, sphere)
		
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
			write(WRITE_UNIT_1, *) 'totalTracer = [', totalTracer(0), ';...'
			do j = 1, timesteps-1
				write(WRITE_UNIT_1, *) totalTracer(j), ';...'
			enddo
			write(WRITE_UNIT_1, *) totalTracer(timesteps), '];'
			
			write(WRITE_UNIT_1, *) 'tracerVar = [', tracerVar(0), ';...'
			do j = 1, timesteps-1
				write(WRITE_UNIT_1, *) tracerVar(j), ';...'
			enddo
			write(WRITE_UNIT_1,*) tracerVar(timesteps), '];'
		endif
	close(WRITE_UNIT_1)
	write(logstring,'(A, F8.2,A)') 'elapsed time = ', (MPI_WTIME() - wallClock)/60.0, ' minutes.'
	call LogMessage(exelog,TRACE_LOGGING_LEVEL,'PROGRAM COMPLETE : ',trim(logstring))
endif

if ( associated(reference) ) then
	call Delete(reference)
	deallocate(reference)
endif
deallocate(totalTracer)
deallocate(tracerVar)
call Delete(timekeeper)
call Delete(remesh)
if ( procRank == 0 ) then 
	call Delete(vtkOut)
	call Delete(meshOut)
endif
call Delete(sphere)
call Delete(exeLog)

call MPI_Finalize(errCode)

contains

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
