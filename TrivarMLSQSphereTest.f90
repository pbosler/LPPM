program TrivariateMLSQTest

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use ParticlesModule
use EdgesModule
use FacesModule
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use SphereGeomModule
use TrivariateMLSQModule
use STDIntVectorModule

implicit none

include 'mpif.h'

! mesh variables
type(PolyMesh2d) :: sphere
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: meshSeed
integer(kint) :: amrLimit
real(kreal) :: radius

! fields for source data
type(Field) :: const
real(kreal), parameter :: CONST_VAL = 2.0_kreal / 3.0_kreal
type(Field) :: linear
type(Field) :: quadratic
type(Field) :: sphHarm
type(Field) :: nearbyIndicator
type(STDIntVector) :: nearFace
type(STDIntVector) :: nearVert
integer(kint) :: faceCheckIndex
integer(kint) :: vertCheckIndex

! interpolation method
type(QuadMLSQ) :: constInterp
type(QuadMLSQ) :: linearInterp
type(QuadMLSQ) :: quadInterp
type(QuadMLSQ) :: harmInterp

! uniform grid for interpolation output
real(kreal) :: lons(360)
real(kreal) :: lats(181)
integer(kint), parameter :: nLat = 181
integer(kint), parameter :: nLon = 360
real(kreal), dimension(nLat, nLon) :: constErr, linearErr, quadErr, harmErr, harmI
integer(kint), dimension(nLat, nLon) :: inFace, nearParticle

! general computing / machine variables
type(MPISetup) :: mpiParticles
type(MPISetup) :: mpiLongitude
integer(kint) :: mpiErrCode, narg
character(len=100) :: arg
type(Logger) :: exeLog
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString
character(len=20) :: logKey = "sphereMLSQ"
real(kreal) :: progStartTime, progEndTime
integer(kint) :: i, j
real(kreal) :: xVec(3), vecA(3), vecB(3), interp, exact

! output variables
character(len=100) :: matlabFilename
character(len=100) :: paraviewFilename

!
!----------------
! PROGRAM START
!----------------
!
call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD, procRank, mpiErrCode)

progStartTime = MPI_WTIME()

call InitLogger(exeLog, procRank)

if ( procRank == 0 ) then
	narg = IARGC()
	if ( narg /= 2 ) then
		stop "usage : spherePSETest.exe meshType initNest"
	else
		call GETARG(1, arg)
		read(arg, *) meshSeed
		if ( meshSeed == 3 ) then
			meshSeed = ICOS_TRI_SPHERE_SEED
		elseif (meshSeed == 4) then
			meshSeed = CUBED_SPHERE_SEED
		else
			stop "invalid arg1: meshType must be 3 or 4"
		endif
		call GETARG(2, arg)
		read(arg, *) initNest
	endif
endif
call MPI_BCAST(meshSeed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiErrCode)
call MPI_BCAST(initNest, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiErrCode)

if ( initNest == 3 ) then
	faceCheckIndex = 800

elseif ( initNest == 4) then
	faceCheckIndex = 1800
	vertCheckIndex = 1901
endif

!
!----------------
! build mesh
!----------------
!
if ( meshSeed == ICOS_TRI_SPHERE_SEED ) then
	write(logString,'(A,I1,A)') "building icos tri mesh to initNest = ", initNest, "..."
	write(paraviewFilename,'(A,I1,A)') "icosTriMLSQ_", initNest, ".vtk"
	write(matlabFilename,'(A,I1,A)') "icosTriMLSQ_", initNest, ".m"
else
	write(logString,'(A,I1,A)') "building cubed sphere mesh to initNest = ", initNest, "..."
	write(paraviewFilename,'(A,I1,A)') "cubedSphereMLSQ_", initNest, ".vtk"
	write(matlabFilename,'(A,I1,A)') "cubedSphereMLSQ_", initNest, ".m"
endif
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, logKey, logString )

maxNest = initNest
amrLimit = 0
radius = 1.0_kreal
call New( sphere, meshSeed, initNest, maxNest, amrLimit, radius )
call LogStats( sphere, exeLog )

call New( mpiParticles, sphere%particles%N, numProcs)

do j = 1, 360
	lons(j) = (j-1) * DEG_2_RAD
enddo
do i = 1, 181
	lats(i) = -PI/2.0_kreal + (i-1) * DEG_2_RAD
enddo

call New( mpiLongitude, nLon, numProcs )

!
!----------------
! set up source field data
!----------------
!
call New(const, 1, sphere%particles%N, "const", "n/a")
call New(linear, 1, sphere%particles%N, "linear", "n/a")
call New(quadratic, 1, sphere%particles%N, "quad", "n/a")
call New(sphHarm, 1, sphere%particles%N, "sphHarm54", "n/a")
call New(nearbyIndicator, 1, sphere%particles%n, "nearbyInd", "n/a")

do i = 1, sphere%particles%N
	xVec = [sphere%particles%x(i), sphere%particles%y(i), sphere%particles%z(i)]
	call InsertScalarToField( const, CONST_VAL )
	call InsertScalarToField( linear, sphere%particles%x(i) )
	call InsertScalarToField( quadratic, xVec(1) * xVec(1) )
	call InsertScalarToField( sphHarm, SphereHarmonic54(xVec) )
enddo

call SetFieldToZero(nearbyIndicator)
nearbyIndicator%N = sphere%particles%N
vertCheckIndex = sphere%faces%vertices(1,faceCheckIndex)

if ( .NOT. sphere%faces%hasChildren(faceCheckIndex) ) then
	call GetParticlesNearFace(nearFace, sphere, faceCheckIndex)
endif
if ( .NOT. sphere%particles%isActive(vertCheckIndex) ) then
	call GetParticlesNearVertex(nearVert, sphere, vertCheckIndex)
endif
do i = 1, nearFace%N
	nearbyIndicator%scalar( nearFace%int(i) ) = 1.0_kreal
enddo
do i = 1, nearVert%N
	nearbyIndicator%scalar( nearVert%int(i) ) = nearbyIndicator%scalar( nearVert%int(i) ) + 3.0_kreal
enddo

!
!----------------
! interpolate
!----------------
!
call New(constInterp, sphere, const)
call New(linearInterp, sphere, linear)
call New(quadInterp, sphere, quadratic)
call New(harmInterp, sphere, sphHarm)

if ( initNest <= 4 ) then
	call LogStats(harmInterp, exeLog)
endif

inFace = 0
nearParticle = 0
do j = 1, nLon
	do i = 1, nLat
		xVec = [ cos(lons(j))*cos(lats(i)), sin(lons(j))*cos(lats(i)), sin(lats(i)) ]
		constErr(i,j) = InterpolateScalar(constInterp, sphere, xVec) - CONST_VAL
		linearErr(i,j) = InterpolateScalar(linearInterp, sphere, xVec) - xVec(1)
		quadErr(i,j) = InterpolateScalar(quadInterp, sphere, xVec) - xVec(1) * xVec(1)
		harmI(i,j) = InterpolateScalar(harmInterp, sphere, xVec)
		harmErr(i,j) = harmI(i,j) - SphereHarmonic54(xVec)
		inFace(i,j) = LocateFaceContainingPoint(sphere, xVec)
		nearParticle(i,j) = nearestParticle(sphere, xVec)
	enddo
enddo

!
!----------------
! write output
!----------------
!
write(logString,'(A,E14.8)') " Test 1, ConstInterp : max rel. err. = ", maxval(abs(constErr))/CONST_VAL
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey), logString)

write(logString,'(A,E14.8)') " Test 2, LinearInterp : max rel. err. = ", maxval(abs(linearErr))
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey), logString)

write(logString,'(A,E14.8)') " Test 3, QuadraticInterp : max rel. err. = ", maxval(abs(quadErr))
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey), logString)

write(logString,'(A,E14.8)') " Test 4, HarmonicInterp : max rel. err. = ", maxval(abs(harmErr)) / &
	maxval(abs(sphHarm%scalar))
call LogMessage(exeLog,TRACE_LOGGING_LEVEL, trim(logKey), logString)

if ( procRank == 0 ) then
	open(unit=WRITE_UNIT_1, file=paraviewFilename, action='WRITE', status='REPLACE')
		call WriteVTKPoints( sphere%particles, WRITE_UNIT_1 )
		call WriteFacesToVTKPolygons( sphere%faces, WRITE_UNIT_1)
		! write all point data before first cell data
		call WriteVTKPointDataSectionHeader(WRITE_UNIT_1, sphere%particles%N)
		call WriteVTKLagCoords( sphere%particles, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( const, WRITE_UNIT_1 )
		call WriteFieldToVTKPointData( linear, WRITE_UNIT_1 )
		call WriteFieldToVTKPointData( quadratic, WRITE_UNIT_1 )
		call WriteFieldToVTKPointData( sphHarm, WRITE_UNIT_1 )
		call WriteFieldToVTKPointData( nearbyIndicator, WRITE_UNIT_1)
	close(WRITE_UNIT_1)
	open(unit=WRITE_UNIT_1, file=matlabFilename, action='WRITE', status='REPLACE')
		call WriteToMatlab( lons, WRITE_UNIT_1, "lons")
		call WriteToMatlab( lats, WRITE_UNIT_1, "lats")
		call WriteToMatlab( constErr, WRITE_UNIT_1, "constErr")
		call WriteToMatlab( linearErr, WRITE_UNIT_1, "linearErr")
		call WriteToMatlab( quadErr, WRITE_UNIT_1, "quadErr")
		call WriteToMatlab( harmErr, WRITE_UNIT_1, "harmErr")
		call WriteToMatlab( harmI, WRITE_UNIT_1, "harmInterp")
		call WriteToMatlab( nearParticle, WRITE_UNIT_1, "nearest")
		call WriteToMatlab( inFace, WRITE_UNIT_1, "inFace")
	close(WRITE_UNIT_1)
endif

!
!----------------
! PROGRAM END
!----------------
!
progEndTime= MPI_WTIME()

write(logstring,'(A,F12.2,A)') "PROGRAM COMPLETE : elapsed time ", progEndTime - progStartTime, " seconds."
call LogMessage(exelog, TRACE_LOGGING_LEVEL, trim(logkey)//" ", logString)

call Delete(harmInterp)
call Delete(quadInterp)
call Delete(linearInterp)
call Delete(constInterp)
call Delete(sphHarm)
call Delete(quadratic)
call Delete(linear)
call Delete(const)
call Delete(mpiLongitude)
call Delete(mpiParticles)
call Delete(sphere)

call MPI_Finalize(mpiErrCode)

contains 

pure function SphereHarmonic54(xyz)
	real(kreal) :: SphereHarmonic54
	real(kreal), intent(in) :: xyz(3)
	SphereHarmonic54 = 3.0_kreal * sqrt(35.0_kreal) * cos( 4.0_kreal * Longitude(xyz)) * sin(Latitude(xyz)) * &
					   (-1.0_kreal + sin(Latitude(xyz)) * sin(Latitude(xyz)))**2
end function

subroutine InitLogger(aLog,rank)
! Initialize a logger for this processor
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	if ( rank == 0 ) then
		call New(aLog,logLevel)
	else
		call New(aLog,ERROR_LOGGING_LEVEL)
	endif
end subroutine

end program