program sphereMeshTester

use NumberKindsModule
use LoggerModule
use FieldModule
use PolyMesh2dModule

implicit none

integer(kint) :: i, j

type(PolyMesh2d) :: icosTri
type(PolyMesh2d) :: cubedSphere
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: amrLimit
real(kreal) :: ampFactor

type(Logger) :: exeLog

integer(kint) :: narg
character(len=100) :: arg

call New(exeLog, DEBUG_LOGGING_LEVEL )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "SphereMeshTest: ", "build a Sphere mesh of triangles and quadrilaterals.")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "SphereMeshTest: ", "test point query algorithm in both meshes.")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "SphereMeshTest: ", "output to VTK/ParaView.")

narg = IARGC()
if (  narg == 0 ) then
	initNest = 0
else 
	call GETARG(1, arg)
	read(arg,*) initNest
endif

maxNest = initNest
amrLimit = 0
ampFactor = 1.0_kreal

call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, "PlanarMeshTest: ", "BUILDING TRIANGULAR MESH.")
call New(triMesh, ICOS_TRI_SPHERE_SEED, initnest, maxnest, amrlimit, ampFactor)
call LogStats(triMesh, exeLog)

if ( procRank == 0 .AND. initNest <= 2 ) then
	print *, "DEBUG : PRINTING ALL MESH INFO  "
	call PrintDebugInfo( triMesh )
endif

call Delete(triMesh)

call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, "PlanarMeshTest: ", "BUILDING QUADRILATERAL MESH.")
call New(quadMesh, CUBED_SPHERE_SEED, initNest, maxNest, amrLimit, ampFactor)
call LogStats(quadMesh, exeLog)

if ( procRank == 0 .AND. initNest <= 2 ) then
	print *, "DEBUG : PRINTING ALL MESH INFO  "
	call PrintDebugInfo( quadMesh )
endif


call Delete(quadMesh)


call Delete(exeLog)
end program 