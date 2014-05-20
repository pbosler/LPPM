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

if ( procRank == 0 ) then
	write(vtkRoot, '(5A)') trim(outputDir), '/vtkOut/', trim(jobPrefix), trim(amrString), '_'
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot), 0, '.vtk'
	call New(vtkOut, sphere, vtkFile, 'spherical harmonic 5 2')
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
call TestSphericalHarmonic52(sphere)
	!
	!	method 1 : trivariate quadratic polynomials, least squares
	!
	testClock = MPI_WTIME()



	write(logstring,'(A,F15.4,A)') 'trivariate quadratic elapsed time = ', MPI_WTIME() - testClock, ' seconds.'
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL,'test 1, method 1 : ', trim(logString) )
	!
	!	method 2 : SSRFPACK
	!
	testClock = MPI_WTIME()


	write(logstring,'(A,F15.4,A)') 'SSRFPACK elapsed time = ', MPI_WTIME() - testClock, ' seconds.'
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL,'test 1, method 2 : ', trim(logString) )
	!
	!	method 3 : singular kernels
	!
	testClock = MPI_WTIME()

	write(logstring,'(A,F15.4,A)') 'singular kernel elapsed time = ', MPI_WTIME() - testClock, ' seconds.'
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL,'test 1, method 3 : ', trim(logString) )

	call VTKOutput(vtkOut, sphere)
call EndSection(exeLog)

!
!	test function 2
!
call TestSphericalHarmonic54(sphere)
call UpdateTitle(vtkOut, 'spherical harmonic 5 4')


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


end program
