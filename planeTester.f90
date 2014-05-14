program PlaneUnitTest

use NumberKindsModule
use LoggerModule
use PlaneMeshModule
use PlaneOutputModule
use PlaneVorticityModule

implicit none


type(PlaneMesh) :: mesh
integer(kint) :: initNest, AMR, nTracer
real(kreal) :: xmin, xmax, ymin, ymax
integer(kint) :: boundaryType = FREE_BOUNDARIES

type(VorticitySetup) :: rankine
real(kreal) :: xc, yc, rad, str

type(PlaneOutput) :: meshOut
character(len=MAX_STRING_LENGTH) :: filename

type(Logger) :: exeLog
character(len=MAX_STRING_LENGTH) :: logstring

call New(exeLog,DEBUG_LOGGING_LEVEL)

xmin = -1.0_kreal
xmax = 1.0_kreal
ymin = -1.0_kreal
ymax = 1.0_kreal

xc = 0.0_kreal
yc = 0.0_kreal
rad = 0.25_kreal
str = 1.0_kreal

!initNest = 0
!AMR  = 0
!nTracer = 0
!filename = '~/Desktop/baseMesh.vtk'
!call New(mesh,initNest,AMR,nTracer)
!call InitializeRectangle(mesh,xmin,xmax,ymin,ymax,boundaryType)
!call New(meshOut,mesh,filename)
!call LogStats(mesh,exeLog)
!!write(logstring,'(A,F16.12)') 'total area = ', TotalArea(mesh)
!call LogMessage(exeLog,TRACE_LOGGING_LEVEL,'nest = 0 ', logstring)
!!call OutputForVTK(meshOut,mesh)
!filename = '~/Desktop/baseMesh.m'
!call UpdateFilename(meshOut,filename)
!call OutputForMatlab(meshOut,mesh)
!!call Delete(mesh)
!call Delete(meshout)

initNest = 4
AMR  = 0
nTracer = 2
filename = '~/Desktop/meshNest.vtk'

call New(mesh,initNest,AMR,nTracer)
call InitializeRectangle(mesh,xmin,xmax,ymin,ymax,boundaryType)
write(logstring,'(A,F16.12)') 'total area = ', TotalArea(mesh)
call LogMessage(exeLog,TRACE_LOGGING_LEVEL,'nest = 2 ', logstring)

call New(rankine, RANKINE_N_INT, RANKINE_N_REAL)
call InitRankineVortex(rankine, xc, yc, rad, str)
call SetRankineVortexOnMesh(mesh, rankine)


call New(meshOut,mesh,filename)
call LogStats(mesh,exeLog)


call OutputForVTK(meshOut,mesh)
filename = '~/Desktop/meshNest.m'
call UpdateFilename(meshOut,filename)
call OutputForMatlab(meshOut,mesh)

call Delete(rankine)
call Delete(mesh)
call Delete(meshOut)
call Delete(exelog)

end program
