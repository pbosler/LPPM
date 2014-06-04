program LambTest

use NumberKindsModule
use LoggerModule
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
type(VorticitySetup) :: lamb
real(kreal) :: xc, yc, rad, u0
!
! tracer variables
!
type(TracerSetup) :: noTracer
!
! remeshing / refinement variables
!
type(RemeshSetup) :: remesh
real(kreal) :: maxCircTol, vortVarTol, lagVarTol
integer(kint) :: amrlimit, remeshInterval, remeshCounter, resetAlphaInterval
type(BIVARSetup) :: reference
!
! timestepping variables
!
type(PlaneRK4DirectSum) :: timekeeper
real(kreal) :: dt, tfinal
integer(kint) :: timeJ, timesteps
!
! input / output variables
!
type(PlaneOutput) :: meshOut
character(len=MAX_STRING_LENGTH) :: filename, fileroot, outputDir, jobPrefix
character(len=128) :: namelistInputFile = 'LambDipole.namelist'
!
! computing environment variables
!
type(Logger) :: exeLog
character(len=MAX_STRING_LENGTH) :: logstring
integer(kint) :: mpiErrCode
real(kreal) :: wallclock
integer(kint), parameter :: BCAST_INT_SIZE = 5, BCAST_REAL_SIZE = 7
integer(kint) :: broadcastIntegers(BCAST_INT_SIZE)
real(kreal) :: broadcastReals(BCAST_REAL_SIZE)
!
! namelists
!
namelist /meshDefine/ initNest, AMR, amrLimit
namelist /vorticityDefine/ u0, rad
namelist /timestepping/ dt, tfinal
namelist /remeshing/ maxCircTol, vortVarTol, lagVarTol, remeshInterval, resetAlphaInterval
namelist /fileIO/ outputDir, jobPrefix

call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numProcs,mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD,procRank,mpiErrCode)

call New(exeLog,DEBUG_LOGGING_LEVEL)

wallclock = MPI_WTIME()

! mesh definition variables
xmin = -2.0_kreal
xmax = 2.0_kreal
ymin = -2.0_kreal
ymax = 2.0_kreal

nTracer = 2

! vorticity definition variables
xc = 0.0_kreal
yc = 0.0_kreal

call ReadNamelistFile(procRank)

! i/o variables
if ( procRank == 0 ) write(filename,'(A,I0.4,A)') trim(fileroot), 0, '.vtk'
!
! build the initial mesh
!
call New(mesh,initNest,AMR,nTracer)
call InitializeRectangle(mesh,xmin,xmax,ymin,ymax,boundaryType)
write(logstring,'(A,I1,A,F16.12)') 'nest = ', initNest, ' , total area = ',  TotalArea(mesh)
call LogMessage(exeLog,TRACE_LOGGING_LEVEL,'mesh info ', logstring)
!
! initialize vorticity
!
call New(lamb, LAMB_DIPOLE_N_INT, LAMB_DIPOLE_N_REAL)
call InitLambDipole(lamb, rad, u0, xc, yc)
call SetLambDipoleOnMesh(mesh, lamb)

!
! initialize remesher
!
call ConvertToRelativeTolerances(mesh, maxCircTol, vortVarTol, lagVarTol)
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'maxCircTol = ', maxCircTol)
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'vortVarTol = ', vortVarTol)
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'lagVarTol  = ', lagVarTol)

call New(remesh, maxCircTol, vortVarTol, lagVarTol, amrlimit)
if ( AMR > 0 ) then
	call InitialRefinement(mesh, remesh, nullTracer, noTracer, SetLambDipoleOnMesh, lamb)
endif
!
! initialize output
!
if ( procRank == 0 ) then
	call New(meshOut,mesh,filename)
	call LogStats(mesh,exeLog)
	call OutputForVTK(meshOut,mesh)
endif
!
! initialize timestepping
!
call New(timekeeper, mesh, numProcs)
timesteps = floor(tfinal / dt)
remeshCounter = 0

do timeJ = 0, timesteps - 1
	!
	! remesh if necessary
	!
	if ( mod(timeJ+1, remeshInterval) == 0 ) then
		remeshCounter = remeshCounter + 1
		!
		! perform the appropriate remeshing procedure
		!
		if ( remeshCounter < resetAlphaInterval ) then
			!
			! remesh to t = 0
			!
			call LagrangianRemeshToInitialTime( mesh, remesh, SetLambDipoleOnMesh, lamb, nullTracer, noTracer)
		elseif ( remeshCounter == resetAlphaInterval ) then
			!
			! remesh to t = 0, create new reference mesh
			!
			call LagrangianRemeshToInitialTime( mesh, remesh, SetLambDipoleOnMesh, lamb, nullTracer, noTracer)
			call New(reference, mesh)
			call LogMessage(exeLog,TRACE_LOGGING_LEVEL,'resetAlpha : ','ReferenceMesh created.')
			call ResetLagrangianParameter(reference)
			call ResetLagrangianParameter(mesh)
		elseif ( remeshCounter > resetAlphaInterval .AND. mod(remeshCounter, resetAlphaInterval) == 0 ) then
			!
			! remesh to previous reference time, create new reference mesh
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
		if ( procRank == 0 ) call Delete(meshOut)
		!
		! create new objects for new mesh
		!
		call New(timekeeper, mesh, numProcs)
		if ( procRank == 0 ) then
			call New(meshOut, mesh, filename)
		endif
	endif
	!
	! increment time
	!
	call RK4TimestepNoRotation( timekeeper, mesh, dt, procRank, numProcs)
	!
	! output timestep data
	!
	if ( procRank == 0 ) then
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, ' t = ', real(timeJ+1,kreal) * dt)
		write(filename,'(A,I0.4,A)') trim(fileroot), timeJ+1, '.vtk'
		call UpdateFilename(meshOut, filename)
		call OutputForVTK(meshout,mesh)
	endif
enddo

if ( procRank == 0 ) then
	write(logstring,'(A, F8.2,A)') 'elapsed time = ', (MPI_WTIME() - wallClock)/60.0, ' minutes.'
	call LogMessage(exelog,TRACE_LOGGING_LEVEL,'PROGRAM COMPLETE : ',trim(logstring))
endif

call Delete(lamb)
call Delete(mesh)
call Delete(remesh)
call Delete(meshOut)
call Delete(exelog)

call MPI_FINALIZE(mpiErrCode)

contains

subroutine ReadNamelistFile(rank)
	integer(kint), intent(in) :: rank
	!
	integer(kint) :: readStat
	if ( rank == 0 ) then
		open(unit = READ_UNIT, file=namelistInputFile, status='OLD', action='READ', iostat=readStat )
		if ( readstat /= 0 ) then
			stop 'cannot read namelist file.'
		endif
		read(READ_UNIT,nml=meshDefine)
		rewind(READ_UNIT)
		read(READ_UNIT,nml=timestepping)
		rewind(READ_UNIT)
		read(READ_UNIT,nml=remeshing)
		rewind(READ_UNIT)
		read(READ_UNIT,nml=fileIO)
		rewind(READ_UNIT)
		read(READ_UNIT, nml=vorticityDefine)

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
		broadcastReals(6) = u0
		broadcastReals(7) = rad

		write(fileRoot,'(4A)') trim(outputDir), 'vtkOut/', trim(jobPrefix), '_'
	endif
	call MPI_BCAST(broadcastIntegers, BCAST_INT_SIZE, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiErrCode)
	initNest = broadcastIntegers(1)
	AMR = broadcastIntegers(2)
	amrLimit = broadcastIntegers(3)
	remeshInterval = broadcastIntegers(4)
	resetAlphaInterval = broadcastIntegers(5)
	call MPI_BCAST(broadcastReals, BCAST_REAL_SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpiErrCode)
	dt = broadcastReals(1)
	tfinal = broadcastReals(2)
	maxCircTol = broadcastReals(3)
	vortVarTol = broadcastReals(4)
	lagVarTol = broadcastReals(5)
	u0 = broadcastReals(6)
	rad = broadcastReals(7)
end subroutine

subroutine ConvertToRelativeTolerances(aMesh, maxCircTol, vortVarTol, lagVarTol)
	type(PlaneMesh), intent(in) :: aMesh
	real(kreal), intent(inout) :: maxCircTol, vortVarTol, lagVarTol
	maxCircTol = maxCircTol * MaximumCirculation(aMesh)
	vortVarTol = vortVarTol * MaximumVorticityVariation(aMesh)
	lagVarTol = lagVarTol * MaximumLagrangianParameterVariation(aMesh)
end subroutine

end program
