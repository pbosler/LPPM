program planarMeshTester

use NumberKindsModule
use LoggerModule
use FieldModule
use PolyMesh2dModule

implicit none

integer(kint) :: i, j

type(PolyMesh2d) :: triMesh
type(PolyMesh2d) :: quadMesh
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: amrLimit
real(kreal) :: ampFactor

type(Logger) :: exeLog

integer(kint), parameter :: nn = 101
real(kreal), parameter :: dx = 0.1
real(kreal), parameter :: xmin = -5.0_kreal
real(kreal), parameter :: xmax = 5.0_kreal
real(kreal), parameter :: ymin = xmin
real(kreal), parameter :: ymax = xmax
real(kreal) :: x(nn) 
real(kreal) :: y(nn)
integer(kint) :: triFaces(nn,nn)
integer(kint) :: quadFaces(nn,nn)

integer(kint) :: narg
character(len=100) :: arg

call New(exeLog, DEBUG_LOGGING_LEVEL )
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "PlanarMeshTest: ", "build a planar mesh of triangles and quadrilaterals.")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "PlanarMeshTest: ", "test point query algorithm in both meshes.")
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "PlanarMeshTest: ", "output to Matlab.")

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
call New(triMesh, TRI_HEX_SEED, initnest, maxnest, amrlimit, ampFactor)
call LogStats(triMesh, exeLog)

if ( procRank == 0 .AND. initNest <= 2 ) then
	print *, "DEBUG : PRINTING ALL MESH INFO  "
	call PrintDebugInfo( triMesh )
endif

open(unit=WRITE_UNIT_1,file="planeTriMeshTestMatlab.m",status='REPLACE')

	call WriteMeshToMatlab(triMesh, WRITE_UNIT_1)

	write(WRITE_UNIT_1,*) "figure(1);clf;hold on;"
	write(WRITE_UNIT_1,*) "plot(x,y,'k+');"
	write(WRITE_UNIT_1,*) "for i = 1:", triMesh%edges%N
	write(WRITE_UNIT_1,*) "		eX = [ x(edgeVerts(i,1)), x(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     eY = [ y(edgeVerts(i,1)), y(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     if edgeHasChildren(i) > 0 "
	write(WRITE_UNIT_1,*) "			plot(eX,eY,'b:x','LineWidth',2);"
	write(WRITE_UNIT_1,*) "     else"
	write(WRITE_UNIT_1,*) "     	plot(eX,eY,'k-o');"
	write(WRITE_UNIT_1,*) "		end "
	write(WRITE_UNIT_1,*) "end"
	write(WRITE_UNIT_1,*) "title('edges and particles');"
	write(WRITE_UNIT_1,*) "xlim([-1.1,1.1]);"
	write(WRITE_UNIT_1,*) "ylim([-1.1,1.1]);"

	write(WRITE_UNIT_1,*) "figure(2);clf;hold on;"
	write(WRITE_UNIT_1,*) "plot(x,y,'k+');"
	write(WRITE_UNIT_1,*) "for i = 1:", triMesh%faces%N
	write(WRITE_UNIT_1,*) "		fX = [ x(faceVerts(i,1)), x(faceVerts(i,2)), x(faceVerts(i,3)), x(faceVerts(i,1))];"
	write(WRITE_UNIT_1,*) "     fY = [ y(faceVerts(i,1)), y(faceVerts(i,2)), y(faceVerts(i,3)), y(faceVerts(i,1))];"
	write(WRITE_UNIT_1,*) "     if faceHasChildren(i) > 0 "
	write(WRITE_UNIT_1,*) "			plot(fX,fY,'b:x','LineWidth',2);"
	write(WRITE_UNIT_1,*) "     else"
	write(WRITE_UNIT_1,*) "     	plot(fX,fY,'k-o');"
	write(WRITE_UNIT_1,*) "		end "
	write(WRITE_UNIT_1,*) "		if faceCenterParticle(i) > 0 "
	write(WRITE_UNIT_1,*) "			plot(x(faceCenterParticle(i)),y(faceCenterParticle(i)),'ro');"
	write(WRITE_UNIT_1,*) "     end"
	write(WRITE_UNIT_1,*) "end"
	write(WRITE_UNIT_1,*) "title('faces and particles');"
	write(WRITE_UNIT_1,*) "xlim([-1.1,1.1]);"
	write(WRITE_UNIT_1,*) "ylim([-1.1,1.1]);"
	
close(WRITE_UNIT_1)

call Delete(triMesh)

call LogMessage(exeLog, DEBUG_LOGGING_LEVEL, "PlanarMeshTest: ", "BUILDING QUADRILATERAL MESH.")
call New(quadMesh, QUAD_RECT_SEED, initNest, maxNest, amrLimit, ampFactor)
call LogStats(quadMesh, exeLog)

if ( procRank == 0 .AND. initNest <= 2 ) then
	print *, "DEBUG : PRINTING ALL MESH INFO  "
	call PrintDebugInfo( quadMesh )
endif

open(unit=WRITE_UNIT_1,file="planeQuadMeshTestMatlab.m",status='REPLACE')

	call WriteMeshToMatlab(quadMesh, WRITE_UNIT_1)

	write(WRITE_UNIT_1,*) "figure(3);clf;hold on;"
	write(WRITE_UNIT_1,*) "plot(x,y,'k+');"
	write(WRITE_UNIT_1,*) "for i = 1:", quadMesh%edges%N
	write(WRITE_UNIT_1,*) "		eX = [ x(edgeVerts(i,1)), x(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     eY = [ y(edgeVerts(i,1)), y(edgeVerts(i,2))];"
	write(WRITE_UNIT_1,*) "     if edgeHasChildren(i) > 0 "
	write(WRITE_UNIT_1,*) "			%plot(eX,eY,'b:x','LineWidth',2);"
	write(WRITE_UNIT_1,*) "     else"
	write(WRITE_UNIT_1,*) "     	plot(eX,eY,'k-o');"
	write(WRITE_UNIT_1,*) "		end "
	write(WRITE_UNIT_1,*) "end"
	write(WRITE_UNIT_1,*) "title('edges and particles');"
	write(WRITE_UNIT_1,*) "xlim([-1.1,1.1]);"
	write(WRITE_UNIT_1,*) "ylim([-1.1,1.1]);"

	write(WRITE_UNIT_1,*) "figure(4);clf;hold on;"
	write(WRITE_UNIT_1,*) "plot(x,y,'k+');"
	write(WRITE_UNIT_1,*) "for i = 1:", quadMesh%faces%N
	write(WRITE_UNIT_1,*) "		fX = [ x(faceVerts(i,1)), x(faceVerts(i,2)), x(faceVerts(i,3)), x(faceVerts(i,4)), x(faceVerts(i,1))];"
	write(WRITE_UNIT_1,*) "     fY = [ y(faceVerts(i,1)), y(faceVerts(i,2)), y(faceVerts(i,3)), y(faceVerts(i,4)), y(faceVerts(i,1))];"
	write(WRITE_UNIT_1,*) "     if faceHasChildren(i) > 0 "
	write(WRITE_UNIT_1,*) "			%plot(fX,fY,'b:x','LineWidth',2);"
	write(WRITE_UNIT_1,*) "     else"
	write(WRITE_UNIT_1,*) "     	plot(fX,fY,'k-o');"
	write(WRITE_UNIT_1,*) "		end "
	write(WRITE_UNIT_1,*) "		if faceCenterParticle(i) > 0 "
	write(WRITE_UNIT_1,*) "			plot(x(faceCenterParticle(i)),y(faceCenterParticle(i)),'ro');"
	write(WRITE_UNIT_1,*) "     end"
	write(WRITE_UNIT_1,*) "end"
	write(WRITE_UNIT_1,*) "title('faces and particles');"
	write(WRITE_UNIT_1,*) "xlim([-1.1,1.1]);"
	write(WRITE_UNIT_1,*) "ylim([-1.1,1.1]);"

	
close(WRITE_UNIT_1)

call Delete(quadMesh)


call Delete(exeLog)
end program 