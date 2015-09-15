program MorePSETests

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use PolyMesh2dModule
use ParticlesModule
use EdgesModule
use FacesModule
use STDIntVectorModule
use FieldModule
use PSEDirectSumModule
use MPISetupModule

implicit none

include 'mpif.h'
!
!	computing environment
!
type(Logger) :: exeLog
character(len=MAX_STRING_LENGTH) :: logstring
character(len=15) :: logKey = "morePSE"
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: namelistFilename
integer(kint) :: mpiErrCode
real(kreal) :: programStart, programEnd
integer(kint) :: i, j

!
!	mesh variables
!
type(PolyMesh2d) :: mesh
integer(kint) :: meshSeed
integer(kint) :: faceKind
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: amrLimit
real(kreal), parameter :: meshRadius = 8.0_kreal
character(len=MAX_STRING_LENGTH) :: meshString
integer(kint) :: meshInts(2)
type(MPISetup) :: mpiParticles
type(PSE) :: pseSetup

!
!	testing mesh
!
integer(kint), parameter :: nn = 501
real(kreal), parameter :: xmin = -8.0_kreal
real(kreal), parameter :: ymin = xmin
real(kreal), parameter :: xmax = 8.0_kreal
real(kreal), parameter :: ymax = xmax
real(kreal), parameter :: dx = (xmax - xmin) / real(nn-1,kreal)
real(kreal) :: xi(nn)
real(kreal) :: yi(nn)
type(MPISetup) :: mpiInterp
real(kreal) :: scalarInterp(nn,nn), scalarExact(nn,nn)
type(Field) :: scalarPlane
type(Field) :: vectorLinear
type(Field) :: estLap
type(Field) :: estDoubleDot
real(kreal) :: interpLoc(3)
real(kreal) :: maxInterpErr

!
!	I/O variables
!
character(len=MAX_STRING_LENGTH) :: outputDir, outputRoot, vtkFile, matlabFile, arg

!--------------------------------
!	initialize : setup computing environment
!--------------------------------
call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD, procRank, mpiErrCode)
programStart = MPI_WTIME()

outputDir = "/Users/pabosle/modelData/SWE"
outputRoot = "morePSETests"

call InitLogger(exeLog, procRank)

do i = 1, nn
	xi(i) = xmin + dx * (i-1)
	yi(i) = ymin + dx * (i-1)
enddo

if ( procRank == 0 ) then
	call GET_COMMAND_ARGUMENT( 1, arg )
	read(arg, *) faceKind
	call GET_COMMAND_ARGUMENT( 2, arg )
	read(arg, *) initNest
	if ( faceKind == 3 ) then	
		meshInts(1) = TRI_HEX_SEED
		write(meshString, '(A,I1)') '_triHex', initNest
	elseif ( faceKind == 4 ) then
		meshInts(1) = QUAD_RECT_SEED
		write(meshString, '(A,I1)') '_quadRect', initNest
	endif
	meshInts(2) = initNest
	
	write(vtkFile,'(5A)') trim(outputDir), '/', trim(outputRoot), trim(meshString), '.vtk'
	write(matlabFile,'(5A)') trim(outputDir), '/', trim(outputRoot), trim(meshString), '.m'
endif
call MPI_BCAST( meshInts, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiErrCode)
meshSeed = meshInts(1)
initNest = meshInts(2)
maxNest = initNest
amrLimit = 0

!--------------------------------
!	initialize : setup problem / build mesh
!--------------------------------

call New(mesh, meshSeed, initNest, maxNest, amrLimit, meshRadius)

call New(pseSetup, mesh)

if ( procRank == 0 ) then
	call LogStats(mesh, exeLog)
	!call OutputToVTKFile( mesh, vtkFile)
endif

call New(scalarPlane, 1, mesh%particles%N, "scalarPlane","n/a")
call New(vectorLinear, 2, mesh%particles%N, "vectorLinear", "n/a")
call New(estLap, 1, mesh%particles%N, "estLapPlane", "n/a")
call New(estDoubleDot, 1, mesh%particles%N, "estDoubleDot","n/a")

do i = 1, mesh%particles%N
	call InsertScalarToField(scalarPlane, PlaneScalar( mesh%particles%x(i), mesh%particles%y(i)) )
	call InsertVectorToField(vectorLinear, LinearVector( mesh%particles%x(i), mesh%particles%y(i)) )
enddo

call New(mpiParticles, mesh%particles%n, numProcs)
call New(mpiInterp, nn, numProcs)

!--------------------------------
!	run : compare PSE results to exact solution
!--------------------------------

interpLoc = 0.0_kreal
do j = mpiInterp%indexStart(procRank), mpiInterp%indexEnd(procRank)
	interpLoc(1) = xi(j)
	do i = 1, nn
		interpLoc(2) = yi(i)
		scalarInterp(i,j) = PSEPlaneInterpolateScalar(pseSetup, mesh, scalarPlane, interpLoc)
	enddo
enddo

do j = 1, nn
	interpLoc(1) = xi(j)
	do i = 1, nn
		interpLoc(2) = yi(i)
		scalarExact(i,j) = sum(interpLoc)
	enddo
enddo

do j = 0, numProcs - 1
	call MPI_BCAST( scalarInterp(:, mpiInterp%indexStart(j):mpiInterp%indexEnd(j)), nn * mpiInterp%messageLength(j), &
			MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, mpiErrCode)
enddo

maxInterpErr = maxval(abs(scalarInterp - scalarExact))/16.0_kreal

call PSEPlaneDoubleDotProductAtParticles(pseSetup, mesh, vectorLinear, estDoubleDot, mpiParticles)

call PSEPlaneLaplacianAtParticles(pseSetup, mesh, scalarPlane, estLap, mpiParticles)

if ( procRank == 0 ) then
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "maxInterpErr = ", maxInterpErr)
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "maxErr-EstLap = ", maxval(abs(estLap%scalar)) )
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, "maxErr-EstDoubleDot = ", maxval(abs(estDoubleDot%scalar - 2.0_kreal))/2.0_kreal )
	open(unit=WRITE_UNIT_1, file=vtkFile, action='WRITE',status='REPLACE')
		! record topology
		call WriteVTKPoints( mesh%particles, WRITE_UNIT_1)
		call WriteFacesToVTKPolygons( mesh%faces, WRITE_UNIT_1)
		! record point data
		call WriteVTKLagCoords(mesh%particles, WRITE_UNIT_1)
		call WriteFieldToVTKPointData(scalarPlane, WRITE_UNIT_1)
		call WriteFieldToVTKPointData(estLap, WRITE_UNIT_1)
		call WriteFieldToVTKPointData(vectorLinear, WRITE_UNIT_1)
		call WriteFieldToVTKPointData(estDoubleDot, WRITE_UNIT_1)
		! record cell data
		call WriteFaceAreaToVTKCellData(mesh%faces, mesh%particles, WRITE_UNIT_1)
	close(WRITE_UNIT_1)
	open(unit=WRITE_UNIT_1, file=matlabFile, action='WRITE', status='REPLACE')
		call WriteParticlesToMatlab( mesh%particles, WRITE_UNIT_1)
		call WriteFieldToMatlab(scalarPlane, WRITE_UNIT_1)
		call WriteFieldToMatlab(estLap, WRITE_UNIT_1)
		call WriteFieldToMatlab(vectorLinear, WRITE_UNIT_1)
		call WriteFieldToMatlab(estDoubleDot, WRITE_UNIT_1)
		call WriteToMatlab(xi, WRITE_UNIT_1, "xi")
		call WriteToMatlab(yi, WRITE_UNIT_1, "yi")
		call WriteToMatlab(scalarInterp, WRITE_UNIT_1,"scalarInterp")
	close(WRITE_UNIT_1)
endif

!--------------------------------
!	finalize : clean up
!--------------------------------

programEnd = MPI_WTIME()

write(logstring,'(A,F12.2,A)') "PROGRAM COMPLETE : elapsed time ", programEnd - programStart, " seconds."
call LogMessage(exelog, TRACE_LOGGING_LEVEL, trim(logkey)//" ", logString)

call Delete(mpiInterp)
call Delete(mpiParticles)
call Delete(estDoubleDot)
call Delete(estLap)
call Delete(vectorLinear)
call Delete(scalarPlane)
call Delete(pseSetup)
call Delete(mesh)

call MPI_FINALIZE(mpiErrCode)

contains

pure function PlaneScalar( x, y)
	real(kreal) :: PlaneScalar
	real(kreal), intent(in) :: x, y
	PlaneScalar = x + y
end function

pure function PlaneGradient( x, y )
	real(kreal) :: PlaneGradient(2)
	real(kreal), intent(in) :: x, y
	PlaneGradient(1) = 1.0_kreal
	PlaneGradient(2) = 1.0_kreal
end function 

pure function LinearVector( x, y )
	real(kreal) :: LinearVector(2)
	real(kreal), intent(in) :: x, y
	LinearVector(1) = x
	LinearVector(2) = y
end function 

pure function LinearVectorDoubleDot( x, y)
	real(kreal) :: LinearVectorDoubleDot
	real(kreal), intent(in) :: x, y
	LinearVectorDoubleDot = 2.0_kreal
end function

subroutine InitLogger(log, rank)
	type(Logger), intent(inout) :: log
	integer(kint), intent(in) :: rank
	if ( rank == 0 ) then
		call New(log, logLevel)
	else
		call New(log, ERROR_LOGGING_LEVEL)
	endif
	write(logKey,'(A,I0.2,A)') trim(logKey)//"_", rank, ":"
end subroutine

end program