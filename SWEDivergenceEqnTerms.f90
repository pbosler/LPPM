program SWEDivergenceEquationTerms

use NumberKindsModule
use LoggerModule
use SphereGeomModule
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
real(kreal) :: test1method1Error(6), test1method2Error(6), test1method3Error(6)

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
nTracer = 3
call New(sphere, panelKind, initNest, AMR, nTracer, SWE_SOLVER)
sphereParticles => sphere%particles
sphereEdges => sphere%edges
spherePanels => sphere%panels

!	TO DO : AMR
write(amrString,'(A,I1)') 'quadNest', initNest

if ( procRank == 0 ) then
	write(vtkRoot, '(5A)') trim(outputDir), '/vtkOut/', trim(jobPrefix), trim(amrString), '_'
	write(vtkFile,'(A,A,I0.4,A)') trim(vtkRoot),'test1_', 0, '.vtk'
	call New(vtkOut, sphere, vtkFile, 'quadratic polynomials in r3')
endif

!
!	For each test function, for each method
!
	! store function values in relvort variable
	! store computed gradient in u variable
	! store computed laplacian in h variable
	! store exact values of laplacian in potvort variable
	! store exact values of gradient in x0 variable
	! store gradient linf error in tracer 1
	! store laplacian linf error in tracer 2 
!
!	output results to vtk, store error norms in log messages
!

!
!	test function 1
!
call StartSection(exeLog,'Test function 1')
call TestLinearFunction(sphere)

	!
	!	method 1 : trivariate quadratic polynomials, least squares
	!
	call StartSection(exeLog, 'trivariate least squares :')
	testClock = MPI_WTIME()
	call TrivariateQuadraticApproximations(sphere)

	write(logstring,'(A,F15.4,A)') 'trivariate quadratic elapsed time = ', MPI_WTIME() - testClock, ' seconds.'
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL,'test 1, method 1 : ', trim(logString) )

	call CalculateError(sphere, test1method1Error(1), test1method1Error(2), test1method1Error(3), test1method1Error(4), &
		test1method1Error(5), test1method1Error(6) )
	call LogMessage(exelog, TRACE_LOGGING_LEVEL, 'scalar interp linf err = ', test1method1Error(5) )
	call LogMessage(exelog, TRACE_LOGGING_LEVEL, 'scalar interp err = ', test1method1Error(6) )
	call LogMessage(exelog, TRACE_LOGGING_LEVEL, 'gradient linf err = ', test1method1Error(1) )
	call LogMessage(exelog, TRACE_LOGGING_LEVEL, 'gradient l2 err = ', test1method1Error(2) )
	call LogMessage(exelog, TRACE_LOGGING_LEVEL, 'Laplacian linf err = ', test1method1Error(3) )
	call LogMessage(exelog, TRACE_LOGGING_LEVEL, 'Laplacian l2 err = ', test1method1Error(4) )
	if ( procRank == 0 ) then
		call VTKOutput(vtkOut, sphere)
	endif
	call EndSection(exeLog)
	
	
	!
	!	method 2 : SSRFPACK
	!
	call TestLinearFunction(sphere)
	call StartSection(exeLog, 'SSRFPACK :')
	testClock = MPI_WTIME()
	call SSRFPACKApproximations(sphere)

	write(logstring,'(A,F15.4,A)') 'SSRFPACK elapsed time = ', MPI_WTIME() - testClock, ' seconds.'
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL,'test 1, method 2 : ', trim(logString) )
	
	call CalculateError(sphere, test1method2Error(1), test1method2Error(2), test1method2Error(3), test1method2Error(4), &
		test1method2Error(5), test1method2Error(6) )
	call LogMessage(exelog, TRACE_LOGGING_LEVEL, 'scalar interp linf err = ', test1method2Error(5) )
	call LogMessage(exelog, TRACE_LOGGING_LEVEL, 'scalar interp err = ', test1method2Error(6) )
	call LogMessage(exelog, TRACE_LOGGING_LEVEL, 'gradient linf err = ', test1method2Error(1) )
	call LogMessage(exelog, TRACE_LOGGING_LEVEL, 'gradient l2 err = ', test1method2Error(2) )
	call LogMessage(exelog, TRACE_LOGGING_LEVEL, 'Laplacian linf err = ', test1method2Error(3) )
	call LogMessage(exelog, TRACE_LOGGING_LEVEL, 'Laplacian l2 err = ', test1method2Error(4) )
	if ( procRank == 0 ) then
		write(vtkFile,'(A,A,I0.4,A)') trim(vtkRoot), 'test1_', 1, '.vtk'
		call UpdateFilename(vtkOut, vtkFile)
		call UpdateTitle(vtkOut, 'ssrfpack')
		call VTKOutput(vtkOut, sphere)
	endif
	call EndSection(exeLog)
	!
	!	method 3 : singular kernels
	!
	testClock = MPI_WTIME()

	write(logstring,'(A,F15.4,A)') 'singular kernel elapsed time = ', MPI_WTIME() - testClock, ' seconds.'
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL,'test 1, method 3 : ', trim(logString) )

	if ( procRank == 0 ) then
		write(vtkFile,'(A,A,I0.4,A)') trim(vtkRoot), 'test1_', 2, '.vtk'
		call UpdateFilename(vtkOut, vtkFile)
		call UpdateTitle(vtkOut, 'singular kernels')
		call VTKOutput(vtkOut, sphere)
	endif
call EndSection(exeLog)

!
!	test function 2
!
!call TestSphericalHarmonic54(sphere)
!call TestSphericalHarmonic52(sphere)
!call UpdateTitle(vtkOut, 'spherical harmonic 5 2')
!
!
!
!	test function 3
!
!call TestSphericalHarmonic10_3(sphere)
!call UpdateTitle(vtkOut, 'spherical harmonic 10 3')



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	Clear memory and Finalize
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write(logstring,'(A, F8.2,A)') 'elapsed time = ', (MPI_WTIME() - wallClock)/60.0, ' minutes.'
call LogMessage(exelog,TRACE_LOGGING_LEVEL,'PROGRAM COMPLETE : ',trim(logstring))



call Delete(sphere)

call MPI_FINALIZE(errCode)

contains

subroutine ReadNamelistInputFile(rank)
	integer(kint), intent(in) :: rank
	integer(kint) :: readWriteStat
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
	real(kreal) :: LegendreP10_3
	real(kreal), intent(in) :: z
	LegendreP10_3 = sqrt(1.0_kreal-z*z)*(-1.0_kreal+z*z)*(-7.0_kreal*z+105.0_kreal*z**3-357.0_kreal*z**5+323.0_kreal*z**7)
end function

function LegendreP54(z)
! Legendre polynomial P_5^4(z)
	real(kreal) :: LegendreP54
	real(kreal), intent(in) :: z
	LegendreP54 = z*(-1.0_kreal + z*z)*(-1.0_kreal + z*z)
end function

function LegendreP52(z)
! Legendre polynomial P_5^2(z)
	real(kreal) :: LegendreP52
	real(kreal), intent(in) :: z
	LegendreP52 = -(-1.0_kreal + z*z)*(-z + 3.0_kreal*z*z*z)
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
		aParticles%x0(:,j) = [1.0_kreal - aParticles%x(1,j) * aParticles%x(1,j) /EARTH_RADIUS/EARTH_RADIUS, &
							  -aparticles%x(2,j)*aParticles%x(1,j)/EARTH_RADIUS/EARTH_RADIUS, &
							  -aParticles%x(3,j)*aparticles%x(1,j)/EARTH_RADIUS/EARTH_RADIUS]
		aParticles%potVort(j) = - 2.0_kreal * aParticles%x(1,j)/EARTH_RADIUS / EARTH_RADIUS
	enddo
	do j = 1, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = apanels%x(1,j)
			aPanels%x0(:,j) = [1.0_kreal - aPanels%x(1,j) * aPanels%x(1,j)  /EARTH_RADIUS/EARTH_RADIUS, &
							   - aPanels%x(2,j) * aPanels%x(1,j) / EARTH_RADIUS / EARTH_RADIUS, &
							   - aPanels%x(3,j) * aPanels%x(1,j) /EARTH_RADIUS/EARTH_RADIUS]
!			aPanels%potVort(j) = -2.0_kreal * cos( Longitude(aPanels%x(:,j)) ) * cos( Latitude( aPanels%x(:,j) ) ) / EARTH_RADIUS
			aPanels%potVort(j) = -2.0_kreal * aPanels%x(1,j) / EARTH_RADIUS / EARTH_RADIUS
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
		aParticles%relVort(j) = LegendreP52( Latitude(aParticles%x(:,j)) ) * cos( 2.0_kreal * Longitude( aParticles%x(:,j) ))

		aParticles%potVort(j) = -30.0_kreal * aParticles%relVort(j)
	enddo
	do j = 1, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = LegendreP52( Latitude(aPanels%x(:,j) )) * cos( 2.0_kreal * Longitude( aPanels%x(:,j) ) )
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
			aPanels%relVort(j) = LegendreP54( Latitude(aPanels%x(:,j) ) )* cos( 4.0_kreal * Longitude( aPanels%x(:,j) ) )
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
			aPanels%relVort(j) = LegendreP10_3( Latitude(aPanels%x(:,j) ) )* cos( 3.0_kreal * Longitude( aPanels%x(:,j) ) )
			aPanels%potVort(j) = -110.0_kreal * aPanels%relVort(j)
		else
			aPanels%relVort(j) = 0.0_kreal
			aPanels%potVort(j) = 0.0_kreal
		endif
	enddo
end subroutine

subroutine SSRFPACKApproximations( aMesh ) 
	type(SphereMesh), intent(inout) :: aMesh
	!
	type(STRIPACKData) :: delTri
	type(SSRFPACKData) :: relVortSource
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j, k, errCode
	real(kreal) :: grad1grad(3), grad2grad(3), grad3grad(3)
	! store computed laplacian in h variable
	! store computed gradient in u variable
	
	aParticles => aMesh%particles
	aPanels => aMesh%panels
	
	! build delaunay triangulation of particles
	call New(delTri, aMesh)
	! use Relative Vorticity as source data for interpolation / approximation
	call New(relVortSource, delTri, .FALSE. )
	call SetSourceRelVort(relVortSource, delTri)
	
	do j = 1, aParticles%N
		aParticles%div(j) = InterpolateScalar(aParticles%x(:,j), relVortSource, delTri)
		aParticles%u(:,j) = relVortSource%grad1(:,aPanels%N_Active + j)
	enddo
	k = 1
	do j = 1, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			aPanels%div(j) = InterpolateScalar(aPanels%x(:,j), relVortSource, delTri)
			aPanels%u(:,delTri%activeMap(k)) = relVortSource%grad1(:,k)
			k = k+1
		endif
	enddo
	
	k = 1
	do j = 1, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			call GRADL(delTri%n, k, delTri%x, delTri%y, delTri%z, relVortSource%grad1(1,:), &
				delTri%list, delTri%lptr, delTri%lend, grad1grad, errCode)
			call GRADL(delTri%n, k, delTri%x, delTri%y, delTri%z, relVortSource%grad1(2,:), &
				delTri%list, delTri%lptr, delTri%lend, grad2grad, errCode)
			call GRADL(delTri%n, k, delTri%x, delTri%y, delTri%z, relVortSource%grad1(3,:), &
				delTri%list, delTri%lptr, delTri%lend, grad3grad, errCode)
			aPanels%h(j) = grad1grad(1) + grad2grad(2) + grad3grad(3)
			k = k + 1
		endif
	enddo
	do j = 1, aParticles%N
		call GRADL(delTri%n, aPanels%N_Active + j, delTri%x, delTri%y, delTri%z, relVortSource%grad1(1,:), &
			delTri%list, delTri%lptr, delTri%lend, grad1grad, errCode)
		call GRADL(delTri%n, aPanels%N_Active + j, delTri%x, delTri%y, delTri%z, relVortSource%grad1(2,:), &
			delTri%list, delTri%lptr, delTri%lend, grad2grad, errCode)
		call GRADL(delTri%n, aPanels%N_Active + j, delTri%x, delTri%y, delTRi%z, relVortSource%grad1(3,:), &
			delTri%list, delTri%lptr, delTri%lend, grad3grad, errCode)
		aParticles%h(j) = grad1grad(1) + grad2grad(2) + grad3grad(3)
	enddo
	
	
	call Delete(delTri)
end subroutine

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
				call LogMessage(exelog,ERROR_LOGGING_LEVEL,trim(logkey)//' potri ERROR : found illegal value at position ',errCode)
				call LogMessage(exelog,ERROR_LOGGING_LEVEL,trim(logKey)//' error found at panel ',j)
			elseif (errCode > 0 ) then
				call LogMessage(exelog,ERROR_LOGGING_LEVEL,logKey,'potri ERROR : found zero on diagonal of Cholesky matrix')
				call LogMessage(exelog,ERROR_LOGGING_LEVEL,trim(logKey)//' error found at panel ',j)
			endif
			! Unpack LAPACK's symmetric storage
			do m=1,10
				do k=m+1,10
					ata(k,m) = ata(m,k)
				enddo
			enddo

			ATAInvAT = matmul(ATA,AT)
			! find coefficients of quadratic approximating polynomial
			coeffs = matmul(ATAInvAT(:,1:n),scalarData(1:n))

			!
			! compute approximations
			!
			! store approximate relVort in div variable
			! store computed laplacian in h variable
			! store computed gradient in u variable

			aPanels%div(j) = QuadraticScalarAppx( coeffs, aPanels%x(:,j))
			aPanels%u(:,j) = QuadraticGradAppx( coeffs, aPanels%x(:,j))
			aPanels%h(j) = QuadraticLaplacianAppx(coeffs, aPanels%x(:,j))
	
			call CCWEdgesAndParticlesAroundPanel(edgeList, vertList, nVerts, aMesh, j)
			do k = 1, nVerts
				aParticles%div(vertList(k)) = aParticles%div(vertList(k)) + QuadraticScalarAppx(coeffs,aParticles%x(:,vertList(k)))
				aParticles%u(:,vertList(k)) = aParticles%u(:,vertList(k)) + &
						QuadraticGradAppx( coeffs, aParticles%x(:,vertList(k)))
				aParticles%h(j) = aParticles%h(j) + QuadraticLaplacianAppx(coeffs, aParticles%x(:,vertList(k)))
				particleChecked(vertList(k)) = particleChecked(vertList(k)) + 1
			enddo
			
		endif
	enddo
	do j = 1, aParticles%N
		aParticles%div(j) = aParticles%div(j) / real(particleChecked(j),kreal)
		aParticles%u(:,j) = aParticles%u(:,j)/real(particleChecked(j), kreal)
		aParticles%h(j) = aParticles%h(j)/real(particleChecked(j), kreal)
	enddo
	deallocate(particleChecked)
end subroutine

function QuadraticScalarAppx( coeffs, xyz)
	real(kreal) :: QuadraticScalarAppx
	real(kreal), intent(in) :: coeffs(10), xyz(3)
	QuadraticScalarAppx = coeffs(1) + coeffs(2) * xyz(1) + coeffs(3) * xyz(2) + coeffs(4) * xyz(3) + &
						  coeffs(5) * xyz(1) * xyz(2) + coeffs(6) * xyz(1) * xyz(3) + coeffs(7) * xyz(2) * xyz(3) + &
						  coeffs(8) * xyz(1) * xyz(1) + coeffs(9) * xyz(2) * xyz(2) + coeffs(10)* xyz(3) * xyz(3)
end function

function QuadraticGradAppx( coeffs, xyz)
	real(kreal) :: QuadraticGradAppx(3)
	real(kreal), intent(in) :: coeffs(10), xyz(3)
	!
	real(kreal) :: cartGrad(3), dotProd
	
	cartGrad = [ coeffs(2) + coeffs(5) * xyz(2) + coeffs(6) * xyz(3) + 2.0_kreal * coeffs(8) * xyz(1) , &
			     coeffs(3) + coeffs(5) * xyz(1) + coeffs(7) * xyz(3) + 2.0_kreal * coeffs(9) * xyz(2) , &
			     coeffs(4) + coeffs(6) * xyz(1) + coeffs(7) * xyz(2) + 2.0_kreal * coeffs(10) * xyz(3) ]
	
	dotProd = sum( xyz * cartGrad)
	
	QuadraticGradAppx = cartGrad - dotProd * xyz / EARTH_RADIUS / EARTH_RADIUS	
end function

function QuadraticLaplacianAppx( coeffs, xyz )
	real(kreal)  :: QuadraticLaplacianAppx
	real(kreal), intent(in) :: coeffs(10), xyz(3)
	!
	real(kreal) :: rsq, rsqmzsq, rsqm2zsq
	
	rsq = EARTH_RADIUS * EARTH_RADIUS
	rSqMZSq = EARTH_RADIUS * EARTH_RADIUS - xyz(3) * xyz(3)
	rsqm2zsq = rsqmzsq - xyz(3) * xyz(3)	
	if ( abs(rsqmzsq) < ZERO_TOL ) then
		QuadraticLaplacianAppx = coeffs(4) * xyz(3) + 2.0_kreal * coeffs(10) * xyz(3) * xyz(3) / rsq - &
								 coeffs(2) * xyz(1) / EARTH_RADIUS - coeffs(3) * xyz(2) / EARTH_RADIUS - coeffs(4) * xyz(3) + &
								 2.0_kreal * coeffs(10)* rsqm2zsq / rsq
	else
		QuadraticLaplacianAppx = -coeffs(2) * xyz(1) / rSqmZsq - coeffs(3) * xyz(2) / rsqmzsq - 4.0_kreal * xyz(1) * xyz(2) / rsqmzsq - &
								 coeffs(6) * xyz(1) * xyz(3) / rsqmzsq - coeffs(7) * xyz(2) * xyz(3) / rsqmzsq - &
								 2.0_kreal * coeffs(8) * ( xyz(1) * xyz(1) - xyz(2) * xyz(2) ) / rsqmzsq + &
								 2.0_kreal * coeffs(9) * ( xyz(1) * xyz(1) - xyz(2) * xyz(2) ) / rsqmzsq - &
								 coeffs(2) * xyz(1) * xyz(3) * xyz(3) / ( rsq * rsqmzsq) - &
								 coeffs(3) * xyz(2) * xyz(3) * xyz(3) / ( rsq * rsqmzsq) + &
								 coeffs(4) * xyz(3) - &
								 2.0_kreal * coeffs(5) * xyz(1) * xyz(2) * xyz(3) * xyz(3) / ( rsq * rsqmzsq) + &
								 coeffs(6) * rsqm2zsq * xyz(1) * xyz(3) / ( rsq * rsqmzsq ) + &
								 coeffs(7) * rsqm2zsq * xyz(2) * xyz(3) / ( rsq * rsqmzsq ) - &
								 2.0_kreal * coeffs(8) * xyz(3) * xyz(3) * xyz(1) * xyz(1) / (rsq * rsqmzsq ) + &
								 2.0_kreal * coeffs(9) * xyz(3) * xyz(3) * xyz(2) * xyz(2) / (rsq * rsqmzsq ) + &
								 2.0_kreal * coeffs(10)* xyz(3) * xyz(3) / rsq - &
								 coeffs(2) * xyz(1) / EARTH_RADIUS - coeffs(3) * xyz(2) / EARTH_RADIUS - coeffs(4) * xyz(3) - &
								 2.0_kreal * coeffs(5) * rsqm2zsq * xyz(1) * xyz(2) / ( rsq * rsqmzsq) - &
								 2.0_kreal * coeffs(8) * rsqm2zsq * xyz(1) * xyz(1) / ( rsq * rsqmzsq) + &
								 2.0_kreal * coeffs(9) * rsqm2zsq * xyz(2) * xyz(2) / ( rsq * rsqmzsq) + &
								 2.0_kreal * coeffs(10)* rsqm2zsq / rsq
	endif
end function


subroutine CalculateError(amesh, gradLinf, gradL2, lapLinf, lapL2, scalarLinf, scalarL2)
	type(SphereMesh), intent(inout) :: amesh
	real(kreal), intent(out) :: gradLinf, gradL2, lapLinf, lapL2, scalarLinf, scalarL2
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
		aParticles%tracer(j,3) = abs( aparticles%div(j) - aParticles%relVort(j))
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
			aPanels%tracer(j,3) = abs( aPanels%div(j) - apanels%relVort(j))
			if ( sqrt(sum( aPanels%x0(:,j)*aPanels%x0(:,j))) > denomGradAinf ) then
				denomGradAinf = sqrt(sum( aPanels%x0(:,j) * aPanels%x0(:,j) ))
			endif
			if ( abs( aPanels%potVort(j) ) > denomLapAinf ) then
				denomLapAinf = abs( aPanels%potVort(j) )
			endif

			denomGradA2 = denomGradA2 + sum(aPanels%x0(:,j) * aPanels%x0(:,j) ) * aPanels%area(j)
			denomLapA2 = denomLapA2 + aPanels%potVort(j) * aPanels%potVort(j) * aPanels%area(j)
		endif
	enddo
	denomGradA2 = sqrt(denomGradA2)
	denomLapA2 = sqrt(denomLapA2)

	if ( max( denomGradPInf, denomGradAinf ) > 2.0_kreal * ZERO_TOL ) then
		aParticles%tracer(:,1) = aParticles%tracer(:,1) / max( denomGradPInf, denomGradAinf )
		aPanels%tracer(:,1) = aPanels%tracer(:,1) / max( denomGradPInf, denomGradAInf)
	endif
	if ( max( denomLapPinf, denomLapAinf) > 2.0_kreal * ZERO_TOL ) then
		aParticles%tracer(:,2) = aParticles%tracer(:,2) / max( denomLapPinf, denomLapAinf)
		aPanels%tracer(:,2) = aPanels%tracer(:,2) / max( denomLapPinf, denomLapAinf)
	endif

	gradLinf = max( maxval(aParticles%tracer(1:aParticles%N,1)), maxval(aPanels%tracer(1:aPanels%N,1))) / &
			   max( denomGradPinf, denomGradAinf)
	lapLinf = max( maxval(aParticles%tracer(1:aParticles%N, 2)), maxval(aPanels%tracer(1:aPanels%N,2))) / &
			 max( denomLapPinf, denomLapAInf)

	gradL2 = sqrt( sum( aPanels%tracer(1:aPanels%N,1) * aPanels%tracer(1:aPanels%N,1) * aPanels%area(1:aPanels%N) ) )/denomGradA2
	lapL2 = sqrt( sum( aPanels%tracer(1:aPanels%N,2) * aPanels%tracer(1:aPanels%N,2) * aPanels%area(1:aPanels%N) ) ) / denomLapA2
	
	scalarLinf = max( maxval(aParticles%tracer(1:aParticles%N,3)), maxval(aPanels%tracer(1:aPanels%N,3)) ) / maxval(aParticles%x(1,:))
	scalarL2 = sqrt(sum( aPanels%tracer(1:aPanels%N,3) * aPanels%tracer(1:aPanels%N,3) * aPanels%area(1:aPanels%N) ))  / &
			   sqrt(sum( aPanels%x(1,:) * aPanels%x(1,:) * aPanels%area(:) ))
end subroutine

end program
