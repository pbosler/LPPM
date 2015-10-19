program SphereMLSQTester

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use PolyMesh2dModule
use ParticlesModule
use EdgesModule
use FacesModule
use STDIntVectorModule
use FieldModule
use SphereGeomModule
use MPISetupModule
use MLSQModule

implicit none

include 'mpif.h'

!
!	logger / console output
!
type(Logger) :: exeLog
character(len=MAX_STRING_LENGTH) :: logstring
character(len=18) :: logKey = "sphereMLSQ"
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
!
!	mesh
!
type(PolyMesh2d) :: sphere
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: amrLimit
real(kreal) :: radius
integer(kint) :: faceType
integer(kint) :: meshSeed
!
!	interpolation
!
type(MLSQ) :: scalarMlsq, vectorMlsq
real(kreal), allocatable :: lons(:)
real(kreal), allocatable :: lats(:)
real(kreal), allocatable :: scalarInterp(:,:), scalarGradX(:,:), scalarGradY(:,:), scalarGradZ(:,:), scalarLap(:,:)
real(kreal), allocatable :: exactScalar(:,:), exactScalarGradX(:,:), exactScalarGradY(:,:), exactScalarGradZ(:,:)
real(kreal), allocatable :: exactScalarLap(:,:)
integer(kint), parameter :: nLon = 360, nLat = 181

type(Field) :: testScalar, testVector
type(Field) :: estGrad
type(Field) :: estLap
type(Field) :: exactGrad
type(Field) :: exactLap
type(Field) :: estDoubleDot
type(Field) :: exactDoubleDot

!
!	MPI
!
type(MPISetup) :: interpMPI
integer(kint) :: mpiErrCode

!
!	general
!
real(kreal) :: programStart, programEnd
integer(kint) :: i, j
character(len=MAX_STRING_LENGTH) :: namelistFile, outputDir, outputRoot, matlabFile, vtkFile
character(len=12) :: meshString
real(kreal) :: xij(3), vecA(3)

namelist /mlsqTest/ initNest, faceType, radius, outputDir, outputRoot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROGRAM START 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Initialize computing environment
!
call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD, procRank, mpiErrCode)

call InitLogger(exeLog, procRank)

call ReadNamelistFile(procRank)

programStart = MPI_WTIME()

!
!	Initialize sphere mesh
!
amrLimit = 0
maxNest = initNest
call New(sphere, meshSeed, initNest, maxNest, amrLimit, radius )

call LogStats(sphere, exeLog)

!
!	Initialize interpolation tests
!
allocate(lons(nLon))
allocate(lats(nLat))
allocate(scalarInterp(nLat, nLon))
allocate(scalarGradX(nLat, nLon))
allocate(scalarGradY(nLat, nLon))
allocate(scalarGradZ(nLat, nLon))
allocate(scalarLap(nLat, nLon))
allocate(exactScalar(nLat, nLon))
allocate(exactScalarGradX(nLat, nLon))
allocate(exactScalarGradY(nLat, nLon))
allocate(exactScalarGradZ(nLat, nLon))
allocate(exactScalarLap(nLat,nLon))

do j = 1, nLon
	lons(j) = real(j-1, kreal) * DEG_2_RAD
enddo
do i = 1, nLat
	lats(i) = -0.5_kreal * PI + real(i-1,kreal) * DEG_2_RAD
enddo

call New(testScalar, 1, sphere%particles%N, 'sphHarm', 'n/a')
call New(testVector, 3, sphere%particles%N, 'tangentVector','n/a')
call New(estGrad, 3, sphere%particles%N, 'estGrad', 'n/a')
call New(exactGrad, 3, sphere%particles%N, 'exactGrad', 'n/a')
call New(estLap, 1, sphere%particles%N, 'estLap', 'n/a')
call New(exactLap, 1, sphere%particles%N, 'exactLap', 'n/a')
call New(estDoubleDot, 1, sphere%particles%N, 'estDoubleDot', 'n/a')
call New(exactDoubleDot, 1, sphere%particles%N,'exactDoubleDot','n/a')

do j = 1, nLon
	do i = 1, nLat
		xij = [ radius * cos(lons(j)) * cos(lats(i)), radius * sin(lons(j)) * cos(lats(i)), radius * sin(lats(i)) ]
		exactScalar(i,j) = SphereHarmonic54( xij )
		vecA = HarmGradient( xij )
		exactScalarGradX(i,j) = vecA(1)
		exactScalarGradY(i,j) = vecA(2)
		exactScalarGradZ(i,j) = vecA(3)
		exactScalarLap(i,j) = - 30.0_kreal * SphereHarmonic54( xij )
	enddo
enddo

do i = 1, sphere%particles%n
	xij = PhysCoord(sphere%particles, i)
	call InsertScalarToField( testScalar, SphereHarmonic54(xij) )
	call InsertVectorToField( testVector, VectorFromLinearPotential(xij))
enddo

call MPI_BARRIER(MPI_COMM_WORLD, mpiErrCode)
call LogMessage(exeLog, TRACE_LOGGING_LEVEL, trim(logKey), " initialization complete... Running problem...")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROGRAM RUN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	perform interpolation and derivative approximation
!
call New( scalarMlsq, sphere, testScalar )

call New( vectorMlsq, sphere, testVector )

call New(interpMPI, nLon, numProcs )

do j = interpMPI%indexStart(procRank), interpMPI%indexEnd(procRank)
	do i = 1, nLat
		xij = [ radius * cos(lons(j)) * cos(lats(i)), radius * sin(lons(j)) * cos(lats(i)), radius * sin(lats(i)) ]
		scalarInterp(i,j) = InterpolateScalar( scalarMlsq, sphere, xij )
		vecA = InterpolateScalarGradient( scalarMlsq, sphere, xij )
		scalarGradX(i,j) = vecA(1)
		scalarGradY(i,j) = vecA(2)
		scalarGradZ(i,j) = vecA(3)
		scalarLap(i,j) = InterpolateScalarLaplacian( scalarMlsq, sphere, xij )
	enddo
enddo

call MLSQGradientAtParticles( scalarMlsq, sphere, estGrad )

call MLSQLaplacianAtParticles( scalarMlsq, sphere, estLap )

call MLSQDoubleDotProductAtParticles( vectorMlsq, sphere, estDoubleDot )

do i = 0, numProcs - 1
	call MPI_BCAST( scalarInterp(:, interpMPI%indexStart(i):interpMPI%indexEnd(i)), nLat * interpMPI%messageLength(i), &
		MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	call MPI_BCAST( scalarGradX(:, interpMPI%indexStart(i):interpMPI%indexEnd(i)), nLat * interpMPI%messageLength(i), &
		MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	call MPI_BCAST( scalarGradY(:, interpMPI%indexStart(i):interpMPI%indexEnd(i)), nLat * interpMPI%messageLength(i), &
		MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)		
	call MPI_BCAST( scalarGradZ(:, interpMPI%indexStart(i):interpMPI%indexEnd(i)), nLat * interpMPI%messageLength(i), &
		MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)		
	call MPI_BCAST( scalarLap(:, interpMPI%indexStart(i):interpMPI%indexEnd(i)), nLat * interpMPI%messageLength(i), &
		MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)		
enddo

!
!	write output
!
if ( procRank == 0 ) then
	open(unit=WRITE_UNIT_1, file=vtkFile, status='REPLACE', action='WRITE')
		! vtk header and point locations
		call WriteVTKPoints( sphere%particles, WRITE_UNIT_1, 'SphereMLSQTester')
		! vtk polydata
		call WriteFacesToVTKPolygons( sphere%faces, WRITE_UNIT_1 )
		! point data section header
		call WriteVTKPointDataSectionHeader( WRITE_UNIT_1, sphere%particles%N)
		!  fields @ points
		call WriteVTKLagCoords( sphere%particles, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( testScalar, WRITE_UNIT_1 )
		call WriteFieldToVTKPointData( estGrad, WRITE_UNIT_1 )
		call WriteFieldToVTKPointData( estLap, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( testVector, WRITE_UNIT_1 )
		call WriteFieldToVTKPointData( estDoubleDot, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( exactGrad, WRITE_UNIT_1)
		call WriteFieldToVTKPointData( exactLap, WRITE_UNIT_1 )
		!	fields @ cells
		call WriteFaceAreaToVTKCellData( sphere%faces, sphere%particles, WRITE_UNIT_1)
	close(WRITE_UNIT_1)
	open(unit=WRITE_UNIT_1, file=matlabFile, status='REPLACE', action='WRITE')
		call WriteToMatlab( lons, WRITE_UNIT_1, "lons" )
		call WriteToMatlab( lats, WRITE_UNIT_1, "lats" )
		call WriteToMatlab( scalarInterp, WRITE_UNIT_1, "scalarInterp")
		call WriteToMatlab( exactScalar, WRITE_UNIT_1, "exactScalar")
		call WriteToMatlab( scalarGradX, WRITE_UNIT_1, "scalarGradX")
		call WriteToMatlab( scalarGradY, WRITE_UNIT_1, "scalarGradY")
		call WriteToMatlab( scalarGradZ, WRITE_UNIT_1, "scalarGradZ")
		call WriteToMatlab( exactScalarGradX, WRITE_UNIT_1, "exactScalarGradX")
		call WriteToMatlab( exactScalarGradY, WRITE_UNIT_1, "exactScalarGradY")
		call WriteToMatlab( exactScalarGradZ, WRITE_UNIT_1, "exactScalarGradZ")
		call WriteToMatlab( scalarLap, WRITE_UNIT_1, "scalarLap")
		call WriteToMatlab( exactScalarLap, WRITE_UNIT_1, "exactScalarLap")
	close(WRITE_UNIT_1)
endif

call LogMessage(exeLog,TRACE_LOGGING_LEVEL, "maxScalarInterp err = ", &
	maxval(abs(scalarInterp-exactScalar))/maxval(abs(exactScalar)) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROGRAM END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
programEnd = MPI_WTIME()

if ( procRank == 0 ) then
	write(6,'(A,F12.2,A)') "PROGRAM COMPLETE : elapsed time = ", programEnd - programStart, " seconds."
endif

call Delete(interpMPI)
call Delete(scalarMlsq)
call Delete(vectorMlsq)
call Delete(estGrad)
call Delete(estLap)
call Delete(estDoubleDot)
call Delete(exactDoubleDot)
call Delete(testVector)
call Delete(testScalar)
deallocate(exactScalarLap)
deallocate(exactScalar)
deallocate(exactScalarGradZ)
deallocate(exactScalarGradY)
deallocate(exactScalarGradX)
deallocate(scalarLap)
deallocate(scalarGradZ)
deallocate(scalarGradY)
deallocate(scalarGradX)
deallocate(scalarInterp)
deallocate(lats)
deallocate(lons)
call Delete(exeLog)
call MPI_FINALIZE(mpiErrCode)

contains

pure function SolidBodyRotationVelocity( xyz )
	real(kreal), dimension(3) :: SolidBodyRotationVelocity
	real(kreal), dimension(3), intent(in) :: xyz
	!
	real(kreal), parameter :: OMG = 2.0_kreal * PI
	
	SolidBodyRotationVelocity(1) = - OMG * xyz(2)
	SolidBodyRotationVelocity(2) =   OMG * xyz(1)
	SolidBodyRotationVelocity(3) = 	 0.0_kreal
end function 

pure function VectorFromLinearPotential( xyz )
	real(kreal), dimension(3) :: VectorFromLinearPotential
	real(kreal), dimension(3), intent(in) :: xyz
	VectorFromLinearPotential(1) = 1.0_kreal - xyz(1) * xyz(1)
	VectorFromLinearPotential(2) = - xyz(1) * xyz(2)
	VectorFromLinearPotential(3) = - xyz(1) * xyz(3)
end function

pure function SphereHarmonic54(xyz)
	real(kreal) :: SphereHarmonic54
	real(kreal), intent(in) :: xyz(3)
	SphereHarmonic54 = 3.0_kreal * sqrt(35.0_kreal) * cos( 4.0_kreal * Longitude(xyz)) * sin(Latitude(xyz)) * &
					   (-1.0_kreal + sin(Latitude(xyz)) * sin(Latitude(xyz)))**2
end function

pure function HarmGradient( xyz ) 
	real(kreal) :: HarmGradient(3)
	real(kreal), intent(in) :: xyz(3)
	!
	real(kreal) :: lon, lat, u, v
	
	lat = Latitude(xyz)
	lon = Longitude(xyz)
	
	u = -12.0_kreal * sqrt(35.0_kreal) * sin(4.0_kreal * lon) * (-1.0_kreal + sin(lat)*sin(lat))**2 * tan(lat)
	v =  12.0_kreal * sqrt(35.0_kreal) * cos(4.0_kreal * lon) * sin(lat) * sin(lat)*(-1.0_kreal + sin(lat)*sin(lat)) + &
	3.0_kreal * sqrt(35.0_kreal) * cos(4.0_kreal * lon)*cos(lat)*(-1.0_kreal + sin(lat)*sin(lat))**2
	
	HarmGradient(1) = -u * sin(lon) - v*sin(lat)*cos(lon)
	HarmGradient(2) = u * cos(lon) - v * sin(lat)*sin(lon)
	HarmGradient(3) = v * cos(lat)
end function

subroutine ReadNamelistFile(rank)
	integer(kint), intent(in) :: rank
	!
	integer(kint), parameter :: initBCAST_intSize = 2
	integer(kint), parameter :: initBCAST_realSize = 1
	integer(kint) :: readStat
	integer(kint) :: bcastIntegers(initBCAST_intSize)
	real(kreal) :: bcastReals(initBCAST_realSize)
	
	if ( COMMAND_ARGUMENT_COUNT() /= 1 ) then
		call LogMessage(exeLog,ERROR_LOGGING_LEVEL, trim(logKey), " expected namelist file name as 1st argument.")
		stop
	endif
	
	if ( rank == 0 ) then
		call GET_COMMAND_ARGUMENT(1, namelistFile)
		
		open(unit = READ_UNIT, file=namelistFile, status='OLD', action='READ', iostat = readstat)
			if (readStat /= 0 ) then	
				call LogMessage(exeLog, ERROR_LOGGING_LEVEL, trim(logKey), " cannot read namelist file.")
				stop
			endif
			
			read(READ_UNIT,nml=mlsqTest)
		close(READ_UNIT)			
		
		if ( faceType == 3 ) then
			meshSeed = ICOS_TRI_SPHERE_SEED
			write(meshString,'(A,I1)') 'icosTri', initNest
		elseif ( faceType == 4) then
			meshSeed = CUBED_SPHERE_SEED
			write(meshString,'(A,I1)') 'cubedSphere', initNest
		else
			call LogMessage(exeLog, WARNING_LOGGING_LEVEL, trim(logKey)//" ReadNamelistFile WARNING :", &
				" invalid faceKind -- using triangles.")
			meshSeed = ICOS_TRI_SPHERE_SEED
			write(meshString,'(A,I1)') 'icosTri', initNest
		endif
	
		write(matlabFile,'(A,A,A,A,A,A)') trim(outputDir), '/', trim(outputRoot), '_', trim(meshString), '.m'
		write(vtkFile,'(A,A,A,A,A,A)') trim(outputDir), '/', trim(outputRoot), '_', trim(meshString), '.vtk'
	
		bcastIntegers(1) = meshSeed
		bcastIntegers(2) = initNest
		
		bcastReals(1) = radius
	endif
	
	call MPI_BCAST( bcastIntegers, initBCAST_intSize, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiErrCode)
	call MPI_BCAST( bcastReals, initBCAST_realSize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpiErrCode)
	
	meshSeed = bcastIntegers(1)
	initNest = bcastIntegers(2)
	
	radius = bcastReals(1)
end subroutine

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