program DeltaKernelInterpTest

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
use PSEDirectSumModule

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

! interpolation method
type(PSE) :: pseSetup

! uniform grid for interpolation output
real(kreal) :: lons(360)
real(kreal) :: lats(181)
integer(kint), parameter :: nLat = 181
integer(kint), parameter :: nLon = 360
real(kreal), dimension(nLat, nLon) :: constErr, linearErr, quadErr, harmErr, harmI
integer(kint), dimension(nLat, nLon) :: inFace, nearParticle
real(kreal) :: constL2, linearL2, quadL2, harmL2, constDenom, linearDenom, quadDenom, harmDenom, cosLat

! general computing / machine variables
type(MPISetup) :: mpiParticles
type(MPISetup) :: mpiLongitude
integer(kint) :: mpiErrCode, narg
integer(kint) :: mpiStatus(MPI_STATUS_SIZE)
character(len=100) :: arg
type(Logger) :: exeLog
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString
character(len=20) :: logKey = "spherePSE"
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

!
!----------------
! build mesh
!----------------
!
if ( meshSeed == ICOS_TRI_SPHERE_SEED ) then
	write(logString,'(A,I1,A)') "building icos tri mesh to initNest = ", initNest, "..."
	write(paraviewFilename,'(A,I1,A)') "icosTriPSE_", initNest, ".vtk"
	write(matlabFilename,'(A,I1,A)') "icosTriPSE_", initNest, ".m"
else
	write(logString,'(A,I1,A)') "building cubed sphere mesh to initNest = ", initNest, "..."
	write(paraviewFilename,'(A,I1,A)') "cubedSpherePSE_", initNest, ".vtk"
	write(matlabFilename,'(A,I1,A)') "cubedSpherePSE_", initNest, ".m"
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

!
!----------------
! interpolate
!----------------
!
call New(pseSetup, sphere )

do j = mpiLongitude%indexStart(procRank), mpiLongitude%indexEnd(procRank)
	do i = 1, nLat
		xVec = [cos(lats(i))*cos(lons(j)), cos(lats(i))*sin(lons(j)), sin(lats(i))]
		constErr(i,j) = PSESphereInterpolateScalar(pseSetup, sphere, const, xVec) - CONST_VAL
		linearErr(i,j) = PSESphereInterpolateScalar(pseSetup, sphere, linear, xVec) - xVec(1)
		quadErr(i,j) = PSESphereInterpolateScalar(pseSetup, sphere, quadratic, xVec) - xVec(1) * xVec(1)
		harmI(i,j) = PSESphereInterpolateScalar(pseSetup, sphere, sphHarm, xVec)
		harmErr(i,j) = harmI(i,j) - SphereHarmonic54(xVec)
	enddo
enddo

do i = 1, numProcs - 1
	if ( procRank == i ) then
		call MPI_SEND( constErr(:,mpiLongitude%indexStart(i):mpiLongitude%indexEnd(i)), nLat*mpiLongitude%messageLength(i), &
			MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, mpiErrCode)
	endif
	if ( procRank == 0 ) then
		call MPI_RECV( constErr(:,mpiLongitude%indexStart(i):mpiLongitude%indexEnd(i)), nLat*mpiLongitude%messageLength(i), &
			MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, mpiStatus, mpiErrCode)
	endif
	if ( procRank == i ) then
		call MPI_SEND( linearErr(:,mpiLongitude%indexStart(i):mpiLongitude%indexEnd(i)), nLat*mpiLongitude%messageLength(i), &
			MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, mpiErrCode)
	endif
	if ( procRank == 0 ) then
		call MPI_RECV( linearErr(:,mpiLongitude%indexStart(i):mpiLongitude%indexEnd(i)), nLat*mpiLongitude%messageLength(i), &
			MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, mpiStatus, mpiErrCode)
	endif
	if ( procRank == i ) then
		call MPI_SEND( quadErr(:,mpiLongitude%indexStart(i):mpiLongitude%indexEnd(i)), nLat*mpiLongitude%messageLength(i), &
			MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, mpiErrCode)
	endif
	if ( procRank == 0 ) then
		call MPI_RECV( quadErr(:,mpiLongitude%indexStart(i):mpiLongitude%indexEnd(i)), nLat*mpiLongitude%messageLength(i), &
			MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, mpiStatus, mpiErrCode)
	endif
	if ( procRank == i ) then
		call MPI_SEND( harmErr(:,mpiLongitude%indexStart(i):mpiLongitude%indexEnd(i)), nLat*mpiLongitude%messageLength(i), &
			MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, mpiErrCode)
	endif
	if ( procRank == 0 ) then
		call MPI_RECV( harmErr(:,mpiLongitude%indexStart(i):mpiLongitude%indexEnd(i)), nLat*mpiLongitude%messageLength(i), &
			MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, mpiStatus, mpiErrCode)
	endif
	if ( procRank == i ) then
		call MPI_SEND( harmI(:,mpiLongitude%indexStart(i):mpiLongitude%indexEnd(i)), nLat*mpiLongitude%messageLength(i), &
			MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, mpiErrCode)
	endif
	if ( procRank == 0 ) then
		call MPI_RECV( harmI(:,mpiLongitude%indexStart(i):mpiLongitude%indexEnd(i)), nLat*mpiLongitude%messageLength(i), &
			MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, mpiStatus, mpiErrCode)
	endif
enddo


!
!----------------
! write output
!----------------
!
if ( procRank == 0 ) then
	write(logString,'(A,E14.8)') " Test 1, ConstInterp : max rel. err. = ", maxval(abs(constErr))/CONST_VAL
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey), logString)

	write(logString,'(A,E14.8)') " Test 2, LinearInterp : max rel. err. = ", maxval(abs(linearErr))
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey), logString)

	write(logString,'(A,E14.8)') " Test 3, QuadraticInterp : max rel. err. = ", maxval(abs(quadErr))
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey), logString)

	write(logString,'(A,E14.8)') " Test 4, HarmonicInterp : max rel. err. = ", maxval(abs(harmErr)) / &
		maxval(abs(sphHarm%scalar))
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL, trim(logKey), logString)
	
	constL2 = 0.0_kreal
	constDenom = 0.0_kreal
	linearL2 = 0.0_kreal
	linearDenom = 0.0_kreal
	quadL2 = 0.0_kreal
	quadDenom = 0.0_kreal
	harmL2 = 0.0_kreal
	harmDenom = 0.0_kreal
	do j = 1, nLon
		do i = 1, nLat
			xVec = [cos(lons(j))*cos(lats(i)), cos(lats(i))*sin(lons(j)), sin(lats(i))]
			cosLat = cos(lats(i))
			constL2 = constL2 + constErr(i,j) * constErr(i,j) * cosLat
			constDenom = constDenom + CONST_VAL * CONST_VAL * cosLat
			linearL2 = linearL2 + linearErr(i,j) * linearErr(i,j) * cosLat
			linearDenom = linearDenom + xVec(1) * xVec(1) * cosLat
			quadL2 = quadL2 + quadErr(i,j) * quadErr(i,j) * cosLat
			quadDenom = quadDenom + ( xVec(1)**4 ) * cosLat
			harmL2 = harmL2 + harmErr(i,j) * harmErr(i,j) * cosLat
			harmDenom = harmDenom + SphereHarmonic54(xVec) * SphereHarmonic54(xVec) * cosLat
		enddo
	enddo
	
		
	write(logString,'(A,E14.8)') " Test 1, ConstInterp : l2 rel. err. = ", constL2 / constDenom
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey), logString)

	write(logString,'(A,E14.8)') " Test 2, LinearInterp : l2 rel. err. = ", linearL2 / linearDenom
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey), logString)

	write(logString,'(A,E14.8)') " Test 3, QuadraticInterp : l2 rel. err. = ", quadL2 / quadDenom
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey), logString)

	write(logString,'(A,E14.8)') " Test 4, HarmonicInterp : l2 rel. err. = ", harmL2 / harmDenom
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL, trim(logKey), logString)
	
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
	close(WRITE_UNIT_1)
	
	open(unit=WRITE_UNIT_1, file=matlabFilename, action='WRITE', status='REPLACE')
		call WriteToMatlab( lons, WRITE_UNIT_1, "lons")
		call WriteToMatlab( lats, WRITE_UNIT_1, "lats")
		call WriteToMatlab( constErr, WRITE_UNIT_1, "constErr")
		call WriteToMatlab( linearErr, WRITE_UNIT_1, "linearErr")
		call WriteToMatlab( quadErr, WRITE_UNIT_1, "quadErr")
		call WriteToMatlab( harmErr, WRITE_UNIT_1, "harmErr")
		call WriteToMatlab( harmI, WRITE_UNIT_1, "harmInterp")
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
call Delete(pseSetup)
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