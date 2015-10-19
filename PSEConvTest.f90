program PSEConvergenceTest

use NumberKindsModule
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

type(Logger) :: exeLog
character(len=MAX_STRING_LENGTH) :: logstring
character(len=15) :: logKey = "planePSE"
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

type(PolyMesh2d) :: planarMesh
integer(kint), parameter :: meshSeed = QUAD_RECT_SEED
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: amrLimit
real(kreal) :: ampFactor 

type(Field) :: scalar
type(Field) :: estGrad, exactGrad, gradError
type(Field) :: estLap, exactLap, lapError
type(Field) :: est2ndPartials, exact2ndPartials, partialErr
real(kreal) :: maxLap, minLap, maxGradMag, minGradMag

type(PSE) :: pseSetup
real(kreal) :: interpLoc(3)

real(kreal), parameter :: b = 3.0_kreal
real(kreal), parameter :: xc = 0.0_kreal, yc = 0.0_kreal
real(kreal), parameter :: maxAbsLap = 36.0_kreal
real(kreal) :: estGradError(9), estLapError(9), interpError(9), meshSize(9)
real(kreal) :: testStart, testEnd, programStart, programEnd

integer(kint) :: i, j
real(kreal) :: xPhys(3), vecA(3), vecB(3)

integer(kint), parameter :: nn = 501
real(kreal), parameter :: dx = 8.0_kreal / real(nn-1,kreal)
real(kreal), parameter :: xmin = -4.0_kreal
real(kreal), parameter :: xmax = 4.0_kreal
real(kreal), parameter :: ymin = xmin
real(kreal), parameter :: ymax = xmax
real(kreal) :: x(nn), y(nn)
real(kreal) :: interpScalar(nn,nn), exactScalar(nn,nn)

character(len=MAX_STRING_LENGTH) :: filename

type(MPISetup) :: particlesMPI, interpMPI
integer(kint) :: mpiErrCode

call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD, procRank, mpiErrCode)

call InitLogger(exeLog, procRank)

programStart = MPI_WTIME()

do i = 1, nn
	x(i) = xmin + dx * (i-1)
	y(i) = ymin + dx * (i-1)
enddo

do initNest = 0, 8
	!call cpu_time(testStart)
	testStart = MPI_WTIME()
	!
	! build the mesh
	!
	maxNest = initNest
	amrLimit = 0
	ampFactor = 3.0_kreal
	
	call New(planarMesh, meshSeed, initNest, maxNest, amrLimit, ampFactor)
	
	call New(particlesMPI, planarMesh%particles%N, numProcs)
	
	call New(interpMPI, nn, numProcs)
	
	!	define the scalar
	call New(scalar, 1, planarMesh%particles%N, "gaussScalar", "n/a")
	call New(estGrad, 2, planarMesh%particles%N, "estGradient", "n/a")
	call New(exactGrad, 2, planarMesh%particles%N, "exactGradient", "n/a")
	call New(gradError, 1, planarMesh%particles%N, "gradError", "n/a")
	call New(estLap, 1, planarMesh%particles%N, "estLap", "n/a")
	call New(exactLap, 1, planarMesh%particles%N, "exactLap", "n/a")
	call New(lapError, 1, planarMesh%particles%N, "lapError", "n/a")
	call New(est2ndPartials, 3, planarMesh%particles%N, "est2ndPartials", "n/a")
	call New(exact2ndPartials, 3, planarMesh%particles%N, "exact2ndPartials", "n/a")
	call New(partialErr, 3, planarMesh%particles%N, "partialErr", "n/a")
	
	do i = 1, planarMesh%particles%N
		xPhys = PhysCoord(planarMesh%particles, i)
		call InsertScalarToField( scalar, Gaussian( xPhys(1:2), b))
		call InsertVectorToField( exactGrad, GaussGrad( xPhys(1:2), b))
		vecA = Gauss2ndDerivs( xPhys(1:2), b)
		call InsertVectorToField( exact2ndPartials, vecA)
		call InsertScalarToField( exactLap, GaussLap( xPhys(1:2), b))
	enddo
	
	call New(pseSetup, planarMesh, 2.0_kreal)
	call PSEPlaneGradientAtParticles( pseSetup, planarMesh, scalar, estGrad, particlesMPI)
	call PSEPlaneSecondPartialsAtParticles( pseSetup, planarMesh, estGrad, est2ndPartials, particlesMPI)
	call PSEPlaneLaplacianAtParticles( pseSetup, planarMesh, scalar, estLap, particlesMPI)
	
	do j = 1, nn
		do i = 1, nn
			exactScalar(i,j) = Gaussian([x(j), y(i)], b)
		enddo
	enddo
	
	interpLoc = 0.0_kreal
	do j = interpMPI%indexStart(procRank), interpMPI%indexEnd(procRank)
		do i = 1, nn
			!exactScalar(i,j) = Gaussian( [x(j),y(i)], b)
			interpLoc(1) = x(j)
			interpLoc(2) = y(i)
			interpScalar(i,j) = PSEPlaneInterpolateScalar(pseSetup, planarMesh, scalar, interpLoc )
		enddo
	enddo
	do j = 0, numProcs-1
		call MPI_BCAST( interpScalar(:, interpMPI%indexStart(j):interpMPI%indexEnd(j)), nn * interpMPI%messageLength(j), &
					    MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, mpiErrCode)
	enddo
	
	do i = 1, planarMesh%particles%N
		vecA = [ estGrad%xComp(i), estGrad%yComp(i), 0.0_kreal]
		vecB = [ exactGrad%xComp(i), exactGrad%yComp(i), 0.0_kreal]
		call InsertScalarToField( gradError, sqrt(sum( (vecB-vecA)*(vecB-vecA))))
		vecA = [ est2ndPartials%xComp(i), est2ndPartials%yComp(i), est2ndPartials%zComp(i)]
		vecB = [ exact2ndPartials%xComp(i), exact2ndPartials%yComp(i), exact2ndPartials%zComp(i)]
		call InsertVectorToField( partialErr, vecB - vecA )
		call InsertScalarToField( lapError, estLap%scalar(i) - exactLap%scalar(i))
	enddo
	
	maxGradMag = MaxMagnitude(exactGrad)
	meshSize(initNest+1) = MaxEdgeLength(planarMesh%edges, planarMesh%particles)
	estGradError(initNest+1) = maxval(gradError%scalar)/maxGradMag
	estLapError(initNest+1) = maxval(abs(lapError%scalar))/maxAbsLap
	interpError(initNest+1) = maxval(abs(interpScalar-exactScalar))
	
	if ( procRank == 0 ) then
		write(logString,'(4A24)') "dx", "gradErr-particles", "lapErr-particles", "interp error"
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, logkey, logString)
		write(logString,'(4F24.10)') meshSize(initNest+1), estGradError(initNest+1), estLapError(initNest+1), interpError(initNest+1)
		call LogMessage(exeLog, TRACE_LOGGING_LEVEL, logkey, logString)
	
		if ( meshSeed == TRI_HEX_SEED ) then
			write(filename, '(A,I1,A)') 'pseTest_triHex', initNest, '.m'
		elseif ( meshSeed == QUAD_RECT_SEED ) then
			write(filename, '(A,I1,A)') 'pseTest_quadRect', initNest, '.m'
		endif

		open(unit=WRITE_UNIT_1, file=filename, status='REPLACE', action='WRITE')
			call WriteParticlesToMatlab( planarMesh%particles, WRITE_UNIT_1)
			call WriteFieldToMatlab( scalar, WRITE_UNIT_1)
			call WriteFieldToMatlab( estGrad, WRITE_UNIT_1)
			call WriteFieldToMatlab( exactGrad, WRITE_UNIT_1)
			call WriteFieldToMatlab( estLap, WRITE_UNIT_1)
			call WriteFieldToMatlab( exactLap, WRITE_UNIT_1)
			call WriteFieldToMatlab( gradError, WRITE_UNIT_1)
			call WriteFieldToMatlab( lapError, WRITE_UNIT_1)
			call WriteFieldToMatlab( est2ndPartials, WRITE_UNIT_1)
			call WriteFieldToMatlab( exact2ndPartials, WRITE_UNIT_1)
			call WriteFieldToMatlab( partialErr, WRITE_UNIT_1)
			write(WRITE_UNIT_1,'(A)',advance='NO') "xi = ["
			do i = 1, nn-1
				write(WRITE_UNIT_1,'(F18.12,A)',advance='NO') x(i), ", "
			enddo
			write(WRITE_UNIT_1,'(F18.12,A)') x(nn), "];"
			write(WRITE_UNIT_1,'(A)') "yi = xi;"
			write(WRITE_UNIT_1,'(A)',advance='NO') "interp = ["
			do i = 1, nn - 1
				do j = 1, nn - 1
					write(WRITE_UNIT_1,'(F18.12,A)',advance='NO') interpScalar(i,j), ", "
				enddo
				write(WRITE_UNIT_1,'(F18.12,A)') interpScalar(i,nn), "; ..."
			enddo
			do j = 1, nn - 1
				write(WRITE_UNIT_1,'(F18.12,A)',advance='NO') interpScalar(nn,j), ", "
			enddo
			write(WRITE_UNIT_1,'(F18.12,A)') interpScalar(nn,nn), "];"
		close(WRITE_UNIT_1)
	endif
	
	!call cpu_time(testEnd)
	testEnd = MPI_WTIME()
	if ( procRank == 0 ) then
		write(logString,'(A,I8,A,F12.2,A)') "nParticles = ", planarMesh%particles%N, ": elapsed time = ", &
			testEnd-testStart, " seconds."
		call LogMessage(exelog, TRACE_LOGGING_LEVEL, logKey, logString)
	endif
	
	!
	! cleanup
	!
	call Delete(particlesMPI)
	call Delete(interpMPI)
	call Delete(pseSetup)
!	call Delete(est2ndPartials)
!	call Delete(exact2ndPartials)
	call Delete(partialErr)
	call Delete(lapError)
	call Delete(exactLap)
	call Delete(estLap)
	call Delete(gradError)
	call Delete(exactGrad)
	call Delete(estGrad)
	call Delete(scalar)
	call Delete(planarMesh)
enddo

if ( procRank == 0 ) then
	write(6,'(4A24)') "dx", "gradErr-particles", "lapErr-particles", "interp error"
	do i = 1, 9
		write(6,'(4F24.10)') meshSize(i), estGradError(i), estLapError(i), interpError(i)
	enddo
	
	programEnd = MPI_WTIME()
	
	write(6,'(A,F12.2,A)') "PROGRAM COMPLETE : elapsed time = ", programEnd - programStart, " seconds."
endif

call Delete(exeLog) 
call MPI_Finalize(mpiErrCode)

contains

function Gaussian( xy, b )
	real(kreal) :: Gaussian
	real(kreal), intent(in) :: xy(2)
	real(kreal), intent(in) :: b
	
	Gaussian = exp( - b * b * ( (xy(1)-xc)*(xy(1)-xc) + (xy(2)-yc)*(xy(2)-yc)))	
end function

function GaussGrad(xy, b)
	real(kreal) :: GaussGrad(2)
	real(kreal), intent(in) :: xy(2)
	real(kreal), intent(in) :: b
	GaussGrad(1) = xy(1)-xc
	GaussGrad(2) = xy(2)-yc
	GaussGrad = -2.0_kreal * GaussGrad * b * b * exp( - b * b * ( (xy(1)-xc)*(xy(1)-xc) + (xy(2)-yc)*(xy(2)-yc) ))
end function

function GaussLap(xy, b)
	real(kreal) :: GaussLap
	real(kreal), intent(in) :: xy(2)
	real(kreal), intent(in) :: b
	GaussLap = ( b*b * ( (xy(1)-xc)*(xy(1)-xc) + (xy(2)-yc)*(xy(2)-yc) ) - 1.0_kreal ) * 4.0_kreal * b * b * &
			exp( - b * b * ( (xy(1)-xc)*(xy(1)-xc) + (xy(2)-yc)*(xy(2)-yc) ))
end function


function Gauss2ndDerivs( xy, b)
	real(kreal) :: Gauss2ndDerivs(3)
	real(kreal), intent(in) :: xy(2)
	real(kreal), intent(in) :: b
	Gauss2ndDerivs(1) = 2.0_kreal * b * b *(2.0_kreal * b*b * (xy(1)-xc)*(xy(1)-xc) - 1.0_kreal)
	Gauss2ndDerivs(2) = 4.0_kreal * b**4 * (xy(1)-xc)*(xy(2)-yc)
	Gauss2ndDerivs(3) = 2.0_kreal * b * b *(2.0_kreal * b*b * (xy(2)-yc)*(xy(2)-yc) - 1.0_kreal)
	Gauss2ndDerivs = Gauss2ndDerivs * exp( - b * b * ( (xy(1)-xc)*(xy(1)-xc) + (xy(2)-yc)*(xy(2)-yc) ))
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
