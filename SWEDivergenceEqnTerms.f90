program SWEDivergenceEquationTerms

use NumberKindsModule
use LoggerModule
use ParticlesModule
use EdgesModule
use PanelsModule
use SphereMeshModule
use LAPACK95
use STRIPACKInterfaceModule
use SSRFPACKInterfaceModule
use VTKOutputModule

implicit none

include 'mpif.h'

!
! mesh variables
!
type(SphereMesh) :: sphere
integer(kint) :: panelKind, initNest, AMR, nTracer
type(Particles), pointer :: sphereParticles => null()
type(Edges), pointer :: sphereEdges => null()
type(Panels), pointer :: spherePanels => null()
namelist /sphereDefine/ panelKind, initNest, AMR

!
! output variables
!
type(VTKSource) :: vtkOut
character(len=MAX_STRING_LENGTH) :: vtkRoot, vtkFile, outputDir, jobPrefix
character(len=MAX_STRING_LENGTH) :: dataFile, summaryFile
character(len=56) :: amrString
namelist /fileIO/ outputDir, jobPrefix

!
! user input
!
character(len=MAX_STRING_LENGTH) :: namelistInputFile = 'sweDerivativeTests.namelist'

!
! logging
!
type(Logger) :: exeLog
character(len=28) :: logkey
character(len=MAX_STRING_LENGTH) :: logstring

!
! general and mpi variables
!
integer(kint) :: i, j, k, errCode
integer(kint), parameter :: BROADCAST_INT_SIZE = 3, BROADCAST_REAL_SIZE = 1
integer(kint) :: broadcastIntegers(BROADCAST_INT_SIZE)
real(kreal) :: broadcastReals(BROADCAST_REAL_SIZE)
real(kreal) :: wallClock, testClock

!
! program-specific variables
!
real(kreal) :: test1method1Error(4), test1method2Error(4), test1method3Error(4)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	PROGRAM START
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
call MPI_INIT(errCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numProcs,errCode)
call MPI_COMM_RANK(MPI_COMM_WORLD,procRank,errCode)

call InitLogger(exeLog,procRank)

wallClock = MPI_WTIME()

!
!	read user input from namelist file, set starting state
!
call ReadNamelistInputFile(procRank)
call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logkey,' initial broadcast done.')


!
!	initialize sphere
!
nTracer = 6
call New(sphere, panelKind, initNest, AMR, nTracer, SWE_SOLVER)
sphereParticles => sphere%particles
sphereEdges => sphere%edges
spherePanels => sphere%panels

!	TO DO : AMR
write(amrString,'(A,I1)') 'quadNest', initNest

if ( procRank == 0 ) then
	write(vtkRoot, '(5A)') trim(outputDir), '/vtkOut/', trim(jobPrefix), trim(amrString), '_'
	write(vtkFile,'(A,A,I0.4,A)') trim(vtkRoot),'test1', 0, '.vtk'
	call New(vtkOut, sphere, vtkFile, 'quadratic polynomials in r3')
endif

!
!	For each test function, for each method
!
	! store spherical harmonic in relvort variable
	! store exact value of laplacian in potvort variable
	! store computed laplacian in h variable
	! store exact value of gradient in x0 variable
	! store computed gradient in u variable
	! store gradient linf error in tracer
	! store laplacian linf error in tracer
!
!	output results to vtk, store error norms in log messages
!

!
!	test function 1
!
call StartSection(exeLog,'Test function 1')

	!
	!	method 1 : trivariate quadratic polynomials, least squares
	!
	testClock = MPI_WTIME()
	call TrivariateQuadraticApproximations(sphere)

	write(logstring,'(A,F15.4,A)') 'trivariate quadratic elapsed time = ', MPI_WTIME() - testClock, ' seconds.'
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL,'test 1, method 1 : ', trim(logString) )

	call CalculateError(amesh, test1method1(1), test1method1(2), test1method1(3), test1method1(4) )
	if ( procRank == 0 ) then
		call VTKOutput(vtkOut, sphere)
	endif

	!
	!	method 2 : SSRFPACK
	!
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL,'test 1, method 2 : ', trim(logString) )
	testClock = MPI_WTIME()

	write(logstring,'(A,F15.4,A)') 'SSRFPACK elapsed time = ', MPI_WTIME() - testClock, ' seconds.'
	if ( procRank == 0 ) then
		write(vtkFile,'(A,A,I0.4,A)') trim(vtkRoot), 'test1', 1, '.vtk'
		call UpdateFilename(vtkOut, vtkFile)
		call UpdateTitle(vtkOut, 'ssrfpack')
		call VTKOutput(vtkOut, sphere)
	endif
	!
	!	method 3 : singular kernels
	!
	testClock = MPI_WTIME()

	write(logstring,'(A,F15.4,A)') 'singular kernel elapsed time = ', MPI_WTIME() - testClock, ' seconds.'
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL,'test 1, method 3 : ', trim(logString) )

	if ( procRank == 0 ) then
		write(vtkFile,'(A,A,I0.4,A)') trim(vtkRoot), 'test1', 2, '.vtk'
		call UpdateFilename(vtkOut, vtkFile)
		call UpdateTitle(vtkOut, 'singular kernels')
		call VTKOutput(vtkOut, sphere)
	endif
call EndSection(exeLog)

!
!	test function 2
!
!call TestSphericalHarmonic54(sphere)
call TestSphericalHarmonic52(sphere)
call UpdateTitle(vtkOut, 'spherical harmonic 5 2')


!
!	test function 3
!
call TestSphericalHarmonic10_3(sphere)
call UpdateTitle(vtkOut, 'spherical harmonic 10 3')



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	Clear memory and Finalize
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write(logstring,'(A, F8.2,A)') 'elapsed time = ', (MPI_WTIME() - wallClock)/60.0, ' minutes.'
call LogMessage(exelog,TRACE_LOGGING_LEVEL,'PROGRAM COMPLETE : ',trim(logstring))



call Delete(sphere)

call MPI_FINALIZE(mpiErrCode)

contains

subroutine ReadNamelistInputFile(rank)
	integer(kint), intent(in) :: rank
	if ( rank == 0 ) then
		open(unit = READ_UNIT, file=namelistInputFile, status = 'OLD', action = 'READ', iostat = readWriteStat)
		if ( readWriteStat /= 0 ) then
			stop 'cannot read namelist file'
		endif
		read(READ_UNIT,nml=sphereDefine)
		rewind(READ_UNIT)
		read(READ_UNIT,nml=fileIO)
		rewind(READ_UNIT)

		broadcastIntegers(1) = panelKind
		broadcastIntegers(2) = initNest
		broadcastIntegers(3) = AMR
	endif

	call MPI_BCAST(broadcastIntegers, BROADCAST_INT_SIZE, MPI_INTEGER, 0, MPI_COMM_WORLD, errCode)
	panelKind = broadcastIntegers(1)
	initNest = broadcastIntegers(2)
	AMR = broadcastIntegers(3)
end subroutine

subroutine InitLogger(aLog, rank)
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	if (rank == 0 ) then
		call New(aLog,DEBUG_LOGGING_LEVEL)
	else
		call New(aLog,DEBUG_LOGGING_LEVEL)
	endif
	write(logkey,'(A,I0.2,A)') 'EXE_LOG',rank,' : '
end subroutine

function LegendreP10_3(z)
	real(kreal) :: LegendreP103
	real(kreal), intent(in) :: z
	LegendreP103 = sqrt(1.0_kreal-z*z)*(-1.0_kreal+z*z)*(-7.0_kreal*z+105.0_kreal*z**3-357.0_kreal*z**5+323.0_kreal*z**7)
end function

function LegendreP54(z)
! Legendre polynomial P_5^4(z)
	real(kreal) :: Legendre54
	real(kreal), intent(in) :: z
	Legendre54 = z*(-1.0_kreal + z*z)*(-1.0_kreal + z*z)
end function

subroutine TestLinearFunction(aMesh)
	type(SphereMesh), intent(inout) :: aMesh
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j
	aParticles => aMesh%particles
	aPanels => aMesh%panels

	do j = 1, aParticles%N
		aParticles%relVort(j) = aParticles%x(1,j)
		aParticles%x0(:,j) = [1.0_kreal - aParticles%x(1,j)/EARTH_RADIUS/EARTH_RADIUS, &
							  -aparticles%x(2,j)*aParticles%x(1,j)/EARTH_RADIUS/EARTH_RADIUS, &
							  -aParticles%x(3,j)*aparticles%x(1,j)/EARTH_RADIUS/EARTH_RADIUS]
		aParticles%potVort(j) = -2.0_kreal * cos( Longitude(aParticles%x(:,j)) ) * cos( Latitude( aparticles%x(:,j) ) ) / EARTH_RADIUS
	enddo
	do j = 1, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = apanels%x(1,j)
			aPanels%x0(:,j) = [1.0_kreal - aPanels%x(1,j)/EARTH_RADIUS/EARTH_RADIUS, &
							   - aPanels%x(2,j) * aPanels%x(1,j) / EARTH_RADIUS / EARTH_RADIUS
							   - aPanels%x(3,j) * aPanels%x(1,j) /EARTH_RADIUS/EARTH_RADIUS]
			aPanels%potVort(j) = -2.0_kreal * cos( Longitude(aPanels%x(:,j)) ) * cos( Latitude( aPanels%x(:,j) ) ) / EARTH_RADIUS
		else
			aPanels%relVort(j) = 0.0_kreal
			aPanels%potVort(j) = 0.0_kreal
		endif
	enddo
end subroutine

subroutine TestSphericalHarmonic52(amesh)
	type(SphereMesh), intent(inout) :: aMesh
	!
	type(Particles), pointer :: aparticles
	type(Panels), pointer :: apanels
	integer(kint) :: j

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	do j = 1, aParticles%N
		aParticles%relVort(j) = Legendre52( Latitude(aParticles%x(:,j)) ) * cos( 2.0_kreal * Longitude( aParticles%x(:,j) ))

		aParticles%potVort(j) = -30.0_kreal * aParticles%relVort(j)
	enddo
	do j = 1, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = Legendre52( Latitude(aPanels%x(:,j) ) * cos( 2.0_kreal * Longitude( aPanels%x(:,j) ) )
			aPanels%potVort(j) = -30.0_kreal * aPanels%relVort(j)
		else
			aPanels%relVort(j) = 0.0_kreal
			aPanels%potVort(j) = 0.0_kreal
		endif
	enddo
end subroutine

subroutine TestSphericalHarmonic54(amesh)
	type(SphereMesh), intent(inout) :: aMesh
	!
	type(Particles), pointer :: aparticles
	type(Panels), pointer :: apanels
	integer(kint) :: j

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	do j = 1, aParticles%N
		aParticles%relVort(j) = LegendreP54( Latitude(aParticles%x(:,j)) ) * cos( 4.0_kreal * Longitude( aParticles%x(:,j) ))

		aParticles%potVort(j) = -30.0_kreal * aParticles%relVort(j)
	enddo
	do j = 1, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = LegendreP54( Latitude(aPanels%x(:,j) ) * cos( 4.0_kreal * Longitude( aPanels%x(:,j) ) )
			aPanels%potVort(j) = -30.0_kreal * aPanels%relVort(j)
		else
			aPanels%relVort(j) = 0.0_kreal
			aPanels%potVort(j) = 0.0_kreal
		endif
	enddo
end subroutine

subroutine TestSphericalHarmonic10_3(amesh)
	type(SphereMesh), intent(inout) :: aMesh
	!
	type(Particles), pointer :: aparticles
	type(Panels), pointer :: apanels
	integer(kint) :: j

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	do j = 1, aParticles%N
		aParticles%relVort(j) = LegendreP10_3( Latitude(aParticles%x(:,j)) ) * cos( 3.0_kreal * Longitude( aParticles%x(:,j) ))

		aParticles%potVort(j) = -110.0_kreal * aParticles%relVort(j)
	enddo
	do j = 1, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = LegendreP10_3( Latitude(aPanels%x(:,j) ) * cos( 3.0_kreal * Longitude( aPanels%x(:,j) ) )
			aPanels%potVort(j) = -110.0_kreal * aPanels%relVort(j)
		else
			aPanels%relVort(j) = 0.0_kreal
			aPanels%potVort(j) = 0.0_kreal
		endif
	enddo
end subroutine

function LegendreP52(z)
! Legendre polynomial P_5^2(z)
	real(kreal) :: Legendre52
	real(kreal), intent(in) :: z
	Legendre52 = -(-1.0_kreal + z*z)*(-z + 3.0_kreal*z*z*z)
end function

subroutine TrivariateQuadraticApproximations( aMesh )
	type(SphereMesh), intent(inout) :: aMesh
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: i, j, k, n, m, errCode
	real(kreal) :: A(31,10), AT(10,31), ATA(10,10), ATAInvAT(10,31), coeffs(10), scalarData(31), datalocs(3,31)
	integer(kint) :: adjPanels(8), nAdj, edgeList(8), vertList(8), nVerts, nearbyParticles(20), nNear
	logical(klog) :: isNewPoint
	integer(kint), allocatable :: particleChecked(:)

	aParticles => aMesh%particles
	aPanels => aMesh%panels
	allocate(particleChecked(aParticles%N))
	particleChecked = 0

	do j = 1, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			! start fresh
			A = 0.0_kreal
			AT = 0.0_kreal
			ATA = 0.0_kreal
			ATAInvAT = 0.0_kreal
			datalocs = 0.0_kreal
			scalardata = 0.0_kreal

			nearbyParticles = 0
			nNear = 0
			adjPanels = 0
			nAdj = 0

			!
			!	get list of nearby passive particles using adjacent panels
			!
			!	start with panel j
			call CCWEdgesAndParticlesAroundPanel(edgeList, vertList, nVerts, aMesh, j)
			nearbyParticles(1:nverts) = vertList
			nNear = nVerts
			! 	continue with panels adjacent to panel j
			call FindAdjacentPanels(adjPanels, nAdj, aMesh, j)
			do k = 1, nAdj
				call CCWEdgesAndParticlesAroundPanel(edgeList, vertList, nVerts, amesh, adjPanels(k) )
				do i = 1, nVerts
					isNewPoint = .TRUE.
					do m = 1, nNear
						if ( nearbyParticles(m) == vertList(i) ) isNewPoint = .FALSE. ! skip duplicates
					enddo
					if ( isNewPoint ) then
						nearbyParticles( nNear + 1 ) = vertList(i)
						nNear = nNear + 1
					endif
				enddo
			enddo

			!
			! gather data and locations for interpolation matrix
			!
			! start with panel j
			dataLocs(:,1) = aPanels%x(:,j)
			scalarData(1) = aPanels%relVort(j)
			! continue with adjacent panels
			do k = 1, nAdj
				dataLocs(:,k+1) = aPanels%x(:, adjPanels(k))
				scalarData(k+1) = aPanels%relVort(adjPanels(k))
			enddo
			! continue with passive particles
			do k = 1, nNear
				dataLocs(:, 1 + nAdj + k) = aParticles%x(:, nearbyParticles(k))
				scalarData(1 + nAdj + k) = aParticles%relVort(nearbyParticles(k))
			enddo
			!	add two points off the sphere to improve the condition number of the interpolation matrix
			dataLocs(:, 1 + nAdj + nNear + 1) = 2.0_kreal * aPanels%x(:, j)
			scalarData(1 + nAdj + nNear + 1) = aPanels%relVort(j)
			dataLocs(:, 1 + nAdj + nNear + 2) = 0.0_kreal
			scalarData(1 + nAdj + nNear + 2) = apanels%relVort(j)

			n = 1 + nAdj + nNear + 2

			!
			!	data collected. ready for least squares
			!

			! setup normal equations
			A(1:n,1) = 1.0_kreal
			A(1:n,2) = dataLocs(1,1:n)
			A(1:n,3) = dataLocs(2,1:n)
			A(1:n,4) = dataLocs(3,1:n)
			A(1:n,5) = dataLocs(1,1:n) * dataLocs(2,1:n)
			A(1:n,6) = dataLocs(1,1:n) * dataLocs(3,1:n)
			A(1:n,7) = dataLocs(2,1:n) * dataLocs(3,1:n)
			A(1:n,8) = dataLocs(1,1:n) * dataLocs(1,1:n)
			A(1:n,9) = dataLocs(2,1:n) * dataLocs(2,1:n)
			A(1:n,10) = dataLocs(3,1:n) * dataLocs(3,1:n)

			AT = transpose(A)
			ATA = matmul(AT,A)

			!
			!	find Cholesky factorization and invert ATA using LAPACK
			!
			call dpotrf('U', 10, ATA, 10, errCode)
			if ( errCode < 0 ) then
				call LogMessage(exeLog, ERROR_LOGGING_LEVEL, 'potrf ERROR : found illegal value at position ',errCode)
			elseif ( errCode > 0 ) then
				call LogMessage(exeLog, ERROR_LOGGING_LEVEL, 'potrf ERROR : ', 'non symmetric positive definite matrix.')
			endif

			call dpotri( 'U', 10, ATA, 10, errCode)
			if ( errCode < 0 ) then
				call LogMessage(log,ERROR_LOGGING_LEVEL,trim(logkey)//' potri ERROR : found illegal value at position ',errCode)
				call LogMessage(log,ERROR_LOGGING_LEVEL,trim(logKey)//' error found at panel ',j)
			elseif (errCode > 0 ) then
				call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'potri ERROR : found zero on diagonal of Cholesky matrix')
				call LogMessage(log,ERROR_LOGGING_LEVEL,trim(logKey)//' error found at panel ',j)
			endif
			! Unpack LAPACK's symmetric storage
			do m=1,10
				do k=m+1,10
					ata(k,m) = ata(m,k)
				enddo
			enddo

			ATAInvAT = matmul(ATA,AT)

			coeffs = matmul(ATAInvAT(:,1:mm),scalarData(1:mm))

			!
			! compute approximations
			!
			call CCWEdgesAndParticlesAroundPanel(edgeList, vertList, nVerts, aMesh, j)
	! store computed laplacian in h variable
	! store computed gradient in u variable
			aPanels%u(:,j) = QuadraticGradAppx( coeffs, aPanels%x(:,j), aPanels%x(:,j))
			aPanels%h(j) = QuadraticLaplacianAppx(coeffs, aPanels%x(:,j))
			do k = 1, nVerts
				aParticles%u(:,vertList(k)) = aParticles%u(:,vertList(k)) + &
						QuadraticGradAppx( coeffs, aParticles%x(:,vertList(k)), aPanels%x(:,j))
				aParticles%h(j) = aParticles%h(j) + QuadraticLaplacianAppx(coeffs, aParticles%x(:,vertList(k)))
				particleChecked(vertList(k)) = particleChecked(vertList(k)) + 1
			enddo
		endif
	enddo
	do j = 1, aParticles%N
		aParticles%u(:,j) = aParticles%u(:,j)/real(particleChecked(j), kreal)
		aParticles%h(j) = aParticles%h(j)/real(particleChecked(j), kreal)
	enddo
	deallocate(particleChecked)
end subroutine

function QuadraticGradAppx( coeffs, xyz, xyzCent)
	real(kreal) :: QuadraticGradAppx(3)
	real(kreal), intent(in) :: coeffs(10), xyz(3), xyzCent(3)
	!
	real(kreal) ::  proj(3,3), xU(3)

	xU = xyzCent / sum( xyzCent * xyzCent)
	proj = 0.0_kreal
	proj(1,1) = 1.0_kreal - xU(1) * xU(1)
	proj(1,2) = - xU(1) * xU(2)
	proj(1,3) = - xU(1) * xU(3)
	proj(2,1) = - xU(1) * xU(2)
	proj(2,2) = 1.0_kreal - xU(2) * xU(2)
	proj(2,3) = - xU(2) * xU(3)
	proj(3,1) = - xU(1) * xU(3)
	proj(3,2) = - xU(2) * xU(3)
	proj(3,3) = 1.0_kreal - xU(3) * xU(3)

	QuadraticGradAppx = [ coeffs(2) + coeffs(5)*xyz(2) + coeffs(6)*xyz(3) + 2.0_kreal * coeffs(8) * xyz(1), &
						  coeffs(3) + coeffs(5)*xyz(1) + coeffs(7)*xyz(3) + 2.0_kreal * coeffs(9) * xyz(2), &
						  coeffs(4) + coeffs(6)*xyz(1) + coeffs(7)*xyz(2) + 2.0_kreal * coeffs(10)* xyz(3) ]
	QuadraticGradAppx = MatMul( proj, quadraticGradAppx)
end function

function QuadraticLaplacianAppx( coeffs, xyz )
	real(kreal)  :: QuadraticLaplacianAppx
	real(kreal), intent(in) :: coeffs(10), xyz(3)
	!
	QuadraticLaplacianAppx = 2.0_kreal * coeffs(8) * ( 1.0_kreal - xyz(1) * xyz(1) / EARTH_RADIUS / EARTH_RADIUS ) + &
							 2.0_kreal * coeffs(9) * ( 1.0_kreal - xyz(2) * xyz(2) / EARTH_RADIUS / EARTH_RADIUS ) + &
							 2.0_kreal * coeffs(10)* ( 1.0_kreal - xyz(3) * xyz(3) / EARTH_RADIUS / EARTH_RADIUS ) - &
							 2.0_kreal * coeffs(5) * xyz(1) * xyz(2) / EARTH_RADIUS / EARTH_RADIUS - &
							 2.0_kreal * coeffs(6) * xyz(1) * xyz(3) / EARTH_RADIUS / EARTH_RADIUS - &
							 2.0_kreal * coeffs(7) * xyz(2) * xyz(3) / EARTH_RADIUS / EARTH_RADIUS
end function


subroutine CalculateError(amesh, gradLinf, gradL2, lapLinf, lapL2)
	type(SphereMesh), intent(inout) :: amesh
	real(kreal), intent(out) :: gradLinf, gradL2, lapLinf, lapL2
	! store spherical harmonic in relvort variable
	! store exact value of laplacian in potvort variable
	! store computed laplacian in h variable
	! store exact value of gradient in x0 variable
	! store computed gradient in u variable
	! store gradient linf error in tracer
	! store laplacian linf error in tracer
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j
	real(kreal) :: denomGradPinf, denomLapPinf, denomGradAinf, denomGradA2, denomLapAinf, denomLapA2

	aParticles => amesh%particles
	aPanels => amesh%panels

	denomGradPinf = 0.0_kreal
	denomLapPinf = 0.0_kreal
	do j = 1, aParticles%N
		aParticles%tracer(j,1) = sqrt(sum( (aParticles%x0(:,j) - aParticles%u(:,j)) * ( aParticles%x0(:,j) - aParticles%u(:,j) ) ))
		aParticles%tracer(j,2) = abs( aParticles%h(j) - aParticles%potVort(j) )
		if ( sqrt( sum( aParticles%x0(:,j)*aParticles%x0(:,j))) > denomGradPinf ) then
			denomGradPInf = sqrt(sum( aParticles%x0(:,j)*aParticles%x0(:,j)))
		endif
		if ( abs( aParticles%potVort(j)) > denomLapPinf  ) then
			denomLapPinf = abs( aParticles%potVort(j) )
		endif
	enddo

	denomGradAinf = 0.0_kreal
	denomGradA2 = 0.0_kreal
	denomLapAinf = 0.0_kreal
	denomLapA2 = 0.0_kreal
	do j = 1, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			aPanels%tracer(j,1) = sqrt(sum( (aPanels%x0(:,j) - aPanels%u(:,j)) * (aPanels%x0(:,j) - aPanels%u(:,j))))
			aPanels%tracer(j,2) = abs( aPanels%h(j) - aPanels%potVort(j))
			if ( sqrt(sum( aPanels%x0(:,j)*aPanels%x0(:,j))) ) then
				denomGradAinf = sqrt(sum( aPanels%x0(:,j) * aPanels%x0(:,j) ))
			endif
			if ( abs( aPanels%potVort(j) ) > denomLapAinf ) then
				denomLapAinf = abs( aPanels%potVort(j) )
			endif

			denomGradA2 = denomGradA2 + sum(aPanels%x0(:,j) * aPanels%x0(:,j) ) * aPanels%area(j)
			denomLapA2 = denomLapA2 + aPanels%potVort * aPanels%potVort(j) * aPanels%area(j)
		endif
	enddo
	denomGradA2 = sqrt(denomGradA2)
	denomLapA2 = sqrt(denomLapA2)

	gradLinf = max( maxval(aParticles%tracer(1:aParticles%N,1)), maxval(aPanels%tracer(1:aPanels%N,1))) / &
			   min( denomGradPinf, denomGradAinf)
	lapLinf = max( maxval(aParticles%tracer(1:aParticles%N, 2)), maxval(aPanels%tracer(1:aPanels%N,2))) / &
			 min( denomLapPinf, denomLapAInf)

	gradL2 = sqrt( sum( aPanels%tracer(1:aPanels%N,1) * aPanels%tracer(1:aPanels%N,1) * aPanels%area(1:aPanels%N) ) )/denomGradA2
	lapL2 = sqrt( sum( aPanels%tracer(1:aPanels%N,2) * aPanels%tracer(1:aPanels%N,2) * aPanels%area(1:aPanels%N) ) ) / denomLapA2
end subroutine

end program
