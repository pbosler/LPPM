program ReversibleDipolesMain

use NumberKindsModule
use LoggerModule
use ParticlesModule
use PanelsModule
use PlaneMeshModule
use PlaneOutputModule
use PlaneVorticityModule
use PlaneDirectSumModule
use PlaneTracerModule
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
! vorticity variables
!
type(VorticitySetup) :: dipoles
real(kreal) :: xc1, yc1, rad1, u1, xc2, yc2, rad2, u2

!
! tracer variables
!
type(tracerSetup) :: dipoleIDs

!
! remeshing / refinement variables
!
type(RemeshSetup) :: remesh
real(kreal) :: maxCircTol, vortVarTol, lagVarTol
integer(kint) :: amrLimit, remeshInterval, remeshCounter, resetAlphaInterval
type(BIVARSetup) :: reference

!
! timestepping variables
!
type(PlaneRK4DirectSum) :: timekeeper
real(kreal) :: dt, tfinal
integer(kint) :: timeJ, timesteps

!
! computation variables
!
real(kreal), allocatable :: totalKE(:), totalEnstrophy(:), totalCirc(:)

!
! input / output variables
!
type(PlaneOutput) :: meshOut
character(len=MAX_STRING_LENGTH) :: vtkFilename, vtkFileroot, outputDir, jobPrefix, summaryFile
character(len=128) :: namelistFile = 'ReversibleDipoles.namelist'
integer(kint) :: frameOut, frameCounter

!
! mpi / computing environment variables
!
type(Logger) :: exeLog
character(len=MAX_STRING_LENGTH) :: logstring
integer(kint) :: errCode
real(kreal) :: wallclock
integer(kint), parameter :: BCAST_INT_SIZE = 5, BCAST_REAL_SIZE = 17
integer(kint) :: broadcastIntegers(BCAST_INT_SIZE)
real(kreal) :: broadcastReals(BCAST_REAL_SIZE)
integer(kint) :: j, iostat

!
! namelists for user input
!
namelist /meshDefine/ initNest, AMR, amrLimit, xmin, xmax, ymin, ymax
namelist /vorticityDefine/ xc1, yc1, rad1, u1, xc2, yc2, rad2, u2
namelist /timestepping/ dt, tfinal
namelist /remeshing/ maxCircTol, vortVarTol, lagVarTol, remeshInterval, resetAlphaInterval
namelist /fileIO/ outputDir, jobPrefix, frameOut

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	INITIALIZE COMPUTER, MESH, TEST CASE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

call MPI_INIT(errCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numProcs,errCode)
call MPI_COMM_RANK(MPI_COMM_WORLD,procRank,errCode)

call New(exeLog,DEBUG_LOGGING_LEVEL)

wallclock = MPI_WTIME()

nTracer = 3

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
! initialize vorticity
!
call New(dipoles, TWO_DIPOLES_N_INT, TWO_DIPOLES_N_REAL)
call InitTwoDipoles(dipoles, xc1, yc1, rad1, u1, xc2, yc2, rad2, u2)
call SetTwoDipolesOnMesh(mesh, dipoles)

!
! initialize tracer
!
call InitDipoleIDTracer(dipoleIDs, xc1, yc1, rad1, xc2, yc2, rad2, 2)
call SetDipoleIDTracerOnMesh(mesh, dipoleIDs)

!
! initialize remeshing
!
call ConvertToRelativeTolerances( mesh, maxCircTol, vortVarTol, lagVarTol)
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'maxCircTol = ', maxCircTol)
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'vortVarTol = ', vortVarTol)
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'lagVarTol  = ', lagVarTol)
call New( remesh, maxCircTol, vortVarTol, lagVarTol, amrLimit)
if ( AMR > 0 ) call InitialRefinement(mesh, remesh, SetDipoleIDTracerOnMesh, dipoleIDs, SetTwoDipolesOnMesh, dipoles)

call StoreLagrangianXCoordinateInTracer(mesh, 1)



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
frameCounter = 1
allocate(totalCirc(0:timesteps))
allocate(totalEnstrophy(0:timesteps))
allocate(totalKE(0:timesteps))
totalCirc = 0.0_kreal
totalEnstrophy = 0.0_kreal
totalKE = 0.0_kreal

totalCirc(0) = GetTotalCirculation(mesh)
totalEnstrophy(0) = GetTotalEnstrophy(mesh)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	RUN THE PROBLEM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
			call LagrangianRemeshToInitialTime( mesh, remesh, SetTwoDipolesOnMesh, dipoles, SetDipoleIDTracerOnMesh, dipoleIDs)
			call StoreLagrangianXCoordinateInTracer(mesh, 1)
			if ( timeJ > timesteps/2) call ReverseVorticity(mesh, reference)
		elseif ( remeshCounter == resetAlphaInterval ) then
			!
			! remesh to t = 0, then build a new reference mesh
			!
			call LagrangianRemeshToInitialTime( mesh, remesh, SetTwoDipolesOnMesh, dipoles, SetDipoleIDTracerOnMesh, dipoleIDs)
			call StoreLagrangianXCoordinateInTracer(mesh, 1)
			if ( timeJ > timesteps/2) call ReverseVorticity(mesh, reference)
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
	call RK4TimestepNoRotation( timekeeper, mesh, dt, procRank, numProcs)

	if ( timeJ + 1 == timesteps/2 ) then
		call ReverseVorticity(mesh, reference)
		call Delete(timekeeper)


		call New(timekeeper,mesh,numProcs)
		if ( procRank == 0 ) then
			call Delete(meshOut)
			call New(meshOut, mesh, vtkFilename)
		endif
	endif

	!
	! compute conserved integrals
	!
	totalCirc(timeJ+1) = GetTotalCirculation(mesh)
	totalEnstrophy(timeJ+1) = GetTotalEnstrophy(mesh)
	totalKE(timeJ+1) = GetTotalKE(mesh)
	if ( timeJ == 0 ) totalKE(0) = totalKE(timeJ+1)

	!
	! output timestep data
	!
	if ( procRank == 0 .AND. mod(timeJ+1,frameOut) == 0 ) then
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 't = ', real(timeJ+1,kreal) * dt)
		write(vtkFilename,'(A,I0.4,A)') trim(vtkFileRoot), frameCounter, '.vtk'
		call UpdateFilename(meshOut, vtkFilename)
		call OutputForVTK(meshOUt, mesh)
		frameCounter = frameCounter + 1
	endif
enddo


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	OUTPUT FINAL DATA
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if ( procRank == 0 ) then
	open(unit=WRITE_UNIT_1,file=summaryFile,status='REPLACE',action='WRITE', iostat=iostat)
	if (iostat == 0 ) then
		write(WRITE_UNIT_1,'(4A24)') 't ', 'totalCirc', 'totalKE', 'totalEns'
		do j = 0, timesteps
			write(WRITE_UNIT_1,'(4F24.12)') j*dt, totalCirc(j), totalKE(j), totalEnstrophy(j)
		enddo
	else
		call LogMessage(exeLog,ERROR_LOGGING_LEVEL,'summaryFile ERROR : ','cannot open datafile.')
	endif
	close(WRITE_UNIT_1)
	write(logstring,'(A, F8.2,A)') 'elapsed time = ', (MPI_WTIME() - wallClock)/60.0, ' minutes.'
	call LogMessage(exelog,TRACE_LOGGING_LEVEL,'PROGRAM COMPLETE : ',trim(logstring))
endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	FREE MEMORY, CLEAN UP, FINALIZE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
deallocate(totalCirc)
deallocate(totalEnstrophy)
deallocate(totalKE)
call Delete(timekeeper)
call Delete(dipoleIDs)
if (procRank == 0 ) call Delete(meshOut)
call Delete(remesh)
call Delete(dipoles)
call Delete(mesh)
call Delete(exeLog)

call MPI_FINALIZE(errCode)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	supplemental functions and subroutines
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ReverseVorticity(aMesh, aRef)
	type(PlaneMesh), intent(inout) :: aMesh
	type(BIVARSetup), intent(inout) :: aRef
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	do j = 1, aParticles%N
		aParticles%relVort(j) = - aParticles%relVort(j)
	enddo
	do j = 1, aPanels%N
		if ( aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = 0.0_kreal
		else
			aPanels%relVort(j) = -aPanels%relVort(j)
		endif
	enddo

	if ( associated(aRef%vort) ) then
		do j = 1, aRef%n
			aRef%vort(j) = - aRef%vort(j)
		enddo
	endif
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
		read(READ_UNIT, nml = vorticityDefine)
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
		broadcastReals(3) = maxCircTol
		broadcastReals(4) = vortVarTol
		broadcastReals(5) = lagVarTol
		broadcastReals(6) = xc1
		broadcastReals(7) = yc1
		broadcastReals(8) = u1
		broadcastReals(9) = rad1
		broadcastReals(10) = xc2
		broadcastReals(11) = yc2
		broadcastReals(12) = u2
		broadcastReals(13) = rad2
		broadcastReals(14) = xmin
		broadcastReals(15) = xmax
		broadcastReals(16) = ymin
		broadcastReals(17) = ymax

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
	maxCircTol = broadcastReals(3)
	vortVarTol = broadcastReals(4)
	lagVarTol = broadcastReals(5)
	xc1 = broadcastReals(6)
	yc1 = broadcastReals(7)
	u1 = broadcastReals(8)
	rad1 = broadcastReals(9)
	xc2 = broadcastReals(10)
	yc2 = broadcastReals(11)
	u2 = broadcastReals(12)
	rad2 = broadcastReals(13)
	xmin = broadcastReals(14)
	xmax = broadcastReals(15)
	ymin = broadcastReals(16)
	ymax = broadcastReals(17)
end subroutine

subroutine ConvertToRelativeTolerances(aMesh, maxCircTol, vortVarTol, lagVarTol)
!	converts relative AMR tolerances (input between 0 and 1) to absolute tolerances (output between 0 and max( panel circulation)
!
	type(PlaneMesh), intent(in) :: aMesh
	real(kreal), intent(inout) :: maxCircTol, vortVarTol, lagVarTol
	maxCircTol = maxCircTol * MaximumCirculation(aMesh)
	vortVarTol = vortVarTol * MaximumVorticityVariation(aMesh)
	lagVarTol = lagVarTol * MaximumLagrangianParameterVariation(aMesh)
end subroutine

subroutine StoreLagrangianXCoordinateInTracer(aMesh, tracerID)
	use ParticlesModule
	use PanelsModule
	implicit none
	type(PlaneMesh), intent(inout) :: aMesh
	integer(kint), intent(in) :: tracerID
	!
	integer(kint) :: j
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	aParticles => aMesh%particles
	aPanels => aMesh%panels
	do j = 1, aParticles%N
		aParticles%tracer(j, tracerID) = aParticles%x0(1,j)
	enddo
	do j = 1, aPanels%N
		aPanels%tracer(j, tracerID) = aPanels%x0(1,j)
	enddo
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
