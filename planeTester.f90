program PlaneUnitTest

use NumberKindsModule
use LoggerModule
use PlaneMeshModule
use PlaneOutputModule
use PlaneVorticityModule
use PlaneDirectSumModule

implicit none

include 'mpif.h'

type(PlaneMesh) :: mesh
integer(kint) :: initNest, AMR, nTracer
real(kreal) :: xmin, xmax, ymin, ymax
integer(kint) :: boundaryType = FREE_BOUNDARIES

type(VorticitySetup) :: rankine
real(kreal) :: xc, yc, rad, str

type(PlaneOutput) :: meshOut
character(len=MAX_STRING_LENGTH) :: filename, fileroot

type(PlaneRK4DirectSum) :: timekeeper
real(kreal) :: dt, tfinal
integer(kint) :: timeJ, timesteps

type(Logger) :: exeLog
character(len=MAX_STRING_LENGTH) :: logstring

integer(kint) :: mpiErrCode
real(kreal) :: wallclock

call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numProcs,mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD,procRank,mpiErrCode)

call New(exeLog,DEBUG_LOGGING_LEVEL)

wallclock = MPI_WTIME()

! mesh definition variables
xmin = -1.0_kreal
xmax = 1.0_kreal
ymin = -1.0_kreal
ymax = 1.0_kreal
initNest = 4
AMR  = 0
nTracer = 2

! vorticity definition variables
xc = 0.0_kreal
yc = 0.0_kreal
rad = 0.25_kreal
str = 1.0_kreal

! time definition variables
dt = 0.01
tfinal = 0.25

! i/o variables
fileroot = '~/Desktop/modelData/vtkOut/'
write(filename,'(A,A,I0.4,A)') trim(fileroot),'planeTest', 0, '.vtk'
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
call New(rankine, RANKINE_N_INT, RANKINE_N_REAL)
call InitRankineVortex(rankine, xc, yc, rad, str)
call SetRankineVortexOnMesh(mesh, rankine)
!
! initialize output
!
call New(meshOut,mesh,filename)
call LogStats(mesh,exeLog)
if ( procRank == 0 ) call OutputForVTK(meshOut,mesh)
!
! initialize timestepping
!
call New(timekeeper, mesh, numProcs)
timesteps = floor(tfinal / dt)
do timeJ = 0, timesteps - 1
	call RK4TimestepNoRotation( timekeeper, mesh, dt, procRank, numProcs)
	if ( procRank == 0 ) then
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, ' t = ', real(timeJ+1,kreal) * dt)
		write(filename,'(A,A,I0.4,A)') trim(fileroot),'planeTest', timeJ+1, '.vtk'
		call UpdateFilename(meshOut, filename)
		call OutputForVTK(meshout,mesh)
	endif
enddo

if ( procRank == 0 ) then
	write(logstring,'(A, F8.2,A)') 'elapsed time = ', (MPI_WTIME() - wallClock)/60.0, ' minutes.'
	call LogMessage(exelog,TRACE_LOGGING_LEVEL,'PROGRAM COMPLETE : ',trim(logstring))
endif

call Delete(rankine)
call Delete(mesh)
call Delete(meshOut)
call Delete(exelog)

call MPI_FINALIZE(mpiErrCode)


end program
