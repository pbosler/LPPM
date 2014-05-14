program PlaneUnitTest

use NumberKindsModule
use LoggerModule
use PlaneMeshModule
use PlaneOutputModule

implicit none


type(PlaneMesh) :: mesh
integer(kint) :: initNest, AMR, nTracer
real(kreal) :: xmin, xmax, ymin, ymax
integer(kint) :: boundaryType = FREE_BOUNDARIES


type(PlaneOutput) :: meshOut
character(len=MAX_STRING_LENGTH) :: filename

type(Logger) :: exeLog

call New(exeLog,DEBUG_LOGGING_LEVEL)

xmin = -1.0_kreal
xmax = 1.0_kreal
ymin = -1.0_kreal
ymax = 1.0_kreal

initNest = 0
AMR  = 0
nTracer = 0
filename = '~/Desktop/baseMesh.vtk'

call New(mesh,initNest,AMR,nTracer)
call InitializeRectangle(mesh,xmin,xmax,ymin,ymax,boundaryType)
call New(meshOut,mesh,filename)
call LogStats(mesh,exeLog)

call OutputForVTK(meshOut,mesh)
filename = '~/Desktop/baseMesh.m'
call UpdateFilename(meshOut,filename)
call OutputForMatlab(meshOut,mesh)

call Delete(mesh)
call Delete(meshout)


initNest = 2
AMR  = 0
nTracer = 2
filename = '~/Desktop/meshNest.vtk'

call New(mesh,initNest,AMR,nTracer)
call InitializeRectangle(mesh,xmin,xmax,ymin,ymax,boundaryType)
call New(meshOut,mesh,filename)
call LogStats(mesh,exeLog)

call OutputForVTK(meshOut,mesh)
filename = '~/Desktop/meshNest.m'
call UpdateFilename(meshOut,filename)
call OutputForMatlab(meshOut,mesh)

call Delete(mesh)
call Delete(meshOut)
call Delete(exelog)

end program
