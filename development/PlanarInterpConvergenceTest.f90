program PlanarInterpConvergence

use NumberKindsModule
use LoggerModule
use PolyMesh2dModule
use ParticlesModule
use EdgesModule
use FacesModule
use STDIntVectorModule
use FieldModule
use TriScatteredInterpModule

implicit none

type(Logger) :: exeLog

type(PolyMesh2d) :: triMesh
integer(kint), parameter :: meshSeed = TRI_HEX_SEED
integer(kint) :: initNest
integer(kint) :: maxNest
integer(kint) :: amrLimit
real(kreal) :: ampFactor

type(Field) :: scalar
type(Field) :: estGrad, exactGrad, gradError
type(Field) :: est2ndPartials, exact2ndPartials, partialsError
type(Field) :: estLap, exactLap, lapError

type(TriScatteredInterp) :: interp

real(kreal), parameter :: b = 3.0_kreal
real(kreal), parameter :: xc = 0.0_kreal, yc = 0.0_kreal

integer(kint) :: i, j

integer(kint), parameter :: nn = 501
real(kreal), parameter :: dx = 8.0_kreal / real( nn - 1, kreal)
real(kreal), parameter :: xmin = -4.0_kreal
real(kreal), parameter :: xmax = 4.0_kreal
real(kreal), parameter :: ymin = xmin
real(kreal), parameter :: ymax = xmax

real(kreal) :: x(nn), y(nn)
real(kreal) :: interpScalar(nn,nn), exactScalar(nn,nn)
real(kreal) :: interpGradX(nn,nn), interpGradY(nn,nn), exactGradX(nn,nn), exactGradY(nn,nn)
real(kreal) :: interpLap(nn,nn), exactLapUnif(nn,nn)
real(kreal) :: xVecA(3), xVecB(3)

real(kreal) :: maxGradMag(9), minGradMag(9), maxLap(9), minLap(9), maxEstGradErr(9), maxEstLapErr(9)
real(kreal) :: maxInterpErr(9), maxGradInterpErr(9), maxLapInterpErr(9)
integer(kint) :: nParticles(9)

character(len=MAX_STRING_LENGTH) :: logstring

!
!----------------
! PROGRAM START
!----------------
!
call New(exeLog, DEBUG_LOGGING_LEVEL)

do initNest = 0, 8
	
	write(logstring,'(A,I3,A)') "test ", initNest+1, ", of 9..."
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL,"Interpolation Convergence : ", logString)

	!
	! build a mesh
	!
	maxNest = initNest
	amrLimit = 0
	ampFactor = 3.0_kreal
	call New(triMesh, meshSeed, initNest, maxNest, amrLimit, ampFactor)
	! define a scalar field on the mesh
	call New(scalar, 1, triMesh%particles%N, "gaussScalar", "n/a")
	call New(estGrad,2, triMesh%particles%N, "estGradient", "n/a")
	call New(exactGrad,2,triMesh%particles%N,"exactGradient","n/a")
	call New(gradError,1,triMesh%particles%N, "gradError", "n/a")
	call New(est2ndPartials,3,triMesh%particles%N, "estPartials","n/a")
	call New(exact2ndPartials,3,triMesh%particles%N,"exactPartials","n/a")
	call New(partialsError,3,triMesh%particles%N, "partialsError","n/a")
	call New(estLap, 1, triMesh%particles%N, "estLaplacian","n/a")
	call New(exactLap,1, triMesh%particles%N,"exactLaplacian", "n/a")
	call New(lapError,1, triMesh%particles%N,"lapError","n/a")
	
	do i = 1, triMesh%particles%N
		call InsertScalarToField( scalar, Gaussian( [triMesh%particles%x(i), triMesh%particles%y(i)], b) )
		call InsertVectorToField(exactGrad, GaussGrad( [triMesh%particles%x(i), triMesh%particles%y(i)], b))
		call InsertVectorToField(exact2ndPartials, Gauss2ndDerivs([triMesh%particles%x(i), triMesh%particles%y(i)], b))
		call InsertScalarToField(exactLap, GaussLap([triMesh%particles%x(i), triMesh%particles%y(i)],b) )
	enddo
	
	!
	! set up interpolation, estimate derivatives
	!
	call New(interp, triMesh)
	call generateFaceUVs(interp, triMesh)
	call estimatePartialDerivatives( interp, triMesh, scalar, estGrad, est2ndPartials)
	call SetPlanarCoefficients(interp, triMesh, scalar, estGrad, est2ndPartials)
	
	!
	! calculate derivative error at particles
	!
	do i = 1, triMesh%particles%N
		xVecA = 0.0_kreal
		xVecB = 0.0_kreal
		xVecA(1:2) = [estGrad%xComp(i), estGrad%yComp(i)]
		xVecB(1:2) = GaussGrad( [triMesh%particles%x(i), triMesh%particles%y(i)], b)
		call InsertScalarToField( gradError, sqrt(sum( xVecA - xVecB)))
		call InsertScalarToField( estLap, est2ndPartials%xComp(i) + est2ndPartials%zComp(i))
		call InsertScalarToField( lapError, estLap%scalar(i) - exactLap%scalar(i))
		xVecA = [est2ndPartials%xComp(i), est2ndPartials%yComp(i), est2ndPartials%zComp(i)]
		xVecB = Gauss2ndDerivs( [triMesh%particles%x(i), triMesh%particles%y(i)], b)
		call InsertVectorToField( partialsError, xVecA - xVecB )
	enddo
	
	nParticles(initNest+1) = triMesh%particles%N
	minGradMag(initNest+1) = MinMagnitude(exactGrad)
	maxGradMag(initNest+1) = MaxMagnitude(exactGrad)
	minLap(initNest+1) = minval(exactLap%scalar)
	maxLap(initNest+1) = maxval(exactLap%scalar)
	maxEstGradErr(initNest+1) = maxval(gradError%scalar)
	maxEstLapErr(initNest+1) = maxval(abs(lapError%scalar))
	
	
	!
	! interpolate to uniform mesh
	!
	do i = 1, nn
		x(i) = xmin + dx * (i-1)
		y(i) = ymin + dx * (i-1)
	enddo
	
	xVecA = 0.0_kreal
	do j = 1, nn
		do i = 1, nn
			exactScalar(i,j) = Gaussian( [x(j), y(i)], b )
			xVecA(1:2) = GaussGrad([x(j),y(i)], b)
			exactGradX(i,j) = xVecA(1)
			exactGradY(i,j) = xVecA(2)
			exactLapUnif(i,j) = GaussLap( [x(j), y(i)], b)
			
			interpScalar(i,j) = InterpolateScalar( interp, triMesh, [x(j),y(i)])
			xVecA(1:2) = InterpolateGradient( interp, triMesh, [x(j), y(i)])
			interpGradX(i,j) = xVecA(1)
			interpGradY(i,j) = xVecA(2)
			interpLap(i,j) = InterpolateLaplacian(interp, triMesh, [x(j), y(i)])
		enddo
	enddo
	
	maxInterpErr(initNest+1) = maxval(abs(interpScalar - exactScalar))
	maxGradInterpErr(initNest+1) = maxval( sqrt( (interpGradX-exactGradX)*(interpGradX-exactGradX) + &
								   (interpGradY - exactGradY)*(interpGradY - exactGradY)))
	maxLapInterpErr(initNest+1) = maxval(abs(interpLap-exactLapUnif))
	
	write(6,'(6A24)') "nparticles", "gradErr-particles", "lapErr-particles", "scalarErr-interp", "gradErr-interp", "lapErr-interp"
	write(6,'(I24,5F24.10)') nParticles(initNest+1), maxEstGradErr(initNest+1)/maxGradMag(initNest+1), &
				  maxEstLapErr(initNest+1)/maxLap(initNest+1), maxInterpErr(initNest+1), &
				  maxGradInterpErr(initNest+1)/maxGradMag(initNest+1), maxLapInterpErr(initNest+1)/maxLap(initNest+1)
	
	!
	! reset for next iteration
	!
	call Delete(interp)
	call Delete(lapError)
	call Delete(exactLap)
	call Delete(estLap)
	call Delete(partialsError)
	call Delete(exact2ndPartials)
	call Delete(est2ndPartials)
	call Delete(gradError)
	call Delete(exactGrad)
	call Delete(estGrad)
	call Delete(scalar)
	call Delete(triMesh)
	
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL, "test complete for initNest = ", initNest)
enddo

write(6,'(6A24)') "nparticles", "gradErr-particles", "lapErr-particles", "scalarErr-interp", "gradErr-interp", "lapErr-interp"
do i = 1,9
	write(6,'(I24,5F24.10)') nParticles(i), maxEstGradErr(i)/maxGradMag(i), maxEstLapErr(i)/maxLap(i), maxInterpErr(i), &
				  maxGradInterpErr(i)/maxGradMag(i), maxLapInterpErr(i)/maxLap(i)
enddo

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
	Gauss2ndDerivs(1) = 2.0_kreal * b * b *(2.0_kreal * b*b *(xy(1)-xc)*(xy(1)-xc) - 1.0_kreal)
	Gauss2ndDerivs(2) = 4.0_kreal * b**4 * (xy(1)-xc)*(xy(2)-yc)
	Gauss2ndDerivs(3) = 2.0_kreal * b * b *(2.0_kreal * b*b * (xy(2)-yc)*(xy(2)-yc) - 1.0_kreal)
	Gauss2ndDerivs = Gauss2ndDerivs * exp( - b * b * ( (xy(1)-xc)*(xy(1)-xc) + (xy(2)-yc)*(xy(2)-yc) ))
end function



end program