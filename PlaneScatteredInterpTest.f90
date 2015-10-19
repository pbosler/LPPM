program PlanarTriScatteredInterpTester

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

real(kreal), parameter :: b = 3.0_kreal

type(Field) :: scalarField, gradientField, secondPartials, gradError, exactGrad, exactPartials, estLaplacian, lapError

type(TriScatteredInterp) :: gaussianInterp

real(kreal), parameter :: xc = 0.25_kreal, yc = 0.01_kreal

integer(kint) :: i, j

integer(kint), parameter :: nn = 101
real(kreal), parameter :: dx = 8.0_kreal/100.0_kreal
real(kreal), parameter :: xmin = -4.0_kreal
real(kreal), parameter :: xmax = 4.0_kreal
real(kreal), parameter :: ymin = xmin
real(kreal), parameter :: ymax = xmax
real(kreal) :: x(nn), y(nn), xVec(3), xVec1(3)
real(kreal) :: interpolatedScalar(nn,nn), interpGradX(nn,nn), interpGradY(nn,nn), interpLap(nn,nn)
real(kreal) :: exactScalar(nn,nn), exactGradX(nn,nn), exactGradY(nn,nn), exactLap(nn,nn)

integer(kint) :: narg
character(len=MAX_STRING_LENGTH) :: filename, arg

!
!----------------
! PROGRAM START
!----------------
!
call New(exeLog, DEBUG_LOGGING_LEVEL)

narg = IARGC()
if (  narg == 0 ) then
	initNest = 0
else 
	call GETARG(1, arg)
	read(arg,*) initNest
endif

maxNest = initNest
amrLimit = 0
ampFactor = 3.0_kreal
call New(triMesh, meshSeed, initNest, maxNest, amrLimit, ampFactor)

call New(scalarField, 1, triMesh%particles%N, "Gaussian", "n/a")
call New(gradientField, 2, triMesh%particles%N, "estGrad", "n/a")
call New(secondPartials, 3, triMesh%particles%N, "est2ndDerivs", "n/a")
call New(gradError, 1, triMesh%particles%N, "estGradError", "length")
call New(exactGrad, 2, triMesh%particles%N, "exactGrad","n/a")
call New(exactPartials,3,triMesh%particles%N,"exact2ndDerivs","n/a")
call New(estLaplacian,1,triMesh%particles%N, "estLap","n/a")
call New(lapError,1,triMesh%particles%N,"lapError", "n/a")

xVec = 0.0_kreal
xVec1 = 0.0_kreal
do i = 1, triMesh%particles%N
	call InsertScalarToField( scalarField, Gaussian([ triMesh%particles%x(i), triMesh%particles%y(i)], b))
	call InsertVectorToField( exactGrad, GaussGrad([triMesh%particles%x(i), triMesh%particles%y(i)], b) )
	call InsertVectorToField( exactPartials, Gauss2ndDerivs([triMesh%particles%x(i), triMesh%particles%y(i)], b) )
enddo
call LogStats(scalarField, exeLog)
call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,"...", "all fields allocated.")

call New(gaussianInterp, triMesh)

call generateFaceUVs(gaussianInterp, triMesh)
call estimatePartialDerivatives(gaussianInterp, triMesh, scalarField, gradientField, secondPartials)
call SetPlanarCoefficients(gaussianInterp, triMesh, scalarField, gradientField, secondPartials)
!call SetPlanarCoefficients(gaussianInterp, triMesh, scalarField, exactGrad, exactPartials)

xVec = 0.0_kreal
xVec1 = 0.0_kreal
do i = 1, triMesh%particles%N
	xVec(1:2) = GaussGrad([triMesh%particles%x(i),triMesh%particles%y(i)],b)
	xVec1(1:2) = [gradientField%xComp(i), gradientField%yComp(i)]
	call InsertScalarToField( gradError, sqrt(sum( (xVec1 - xVec)*(xVec1-xVec))))
	call InsertScalarToField( estLaplacian, secondPartials%xComp(i) + secondPartials%zComp(i))
	call InsertScalarToField( lapError, secondPartials%xComp(i) + secondPartials%zComp(i) - &
					GaussLap( [triMesh%particles%x(i), triMesh%particles%y(i)], b) )
enddo

do i = 1, nn
	x(i) = xmin + dx * (i-1)
	y(i) = ymin + dx * (i-1)
enddo

do j = 1, nn
	do i = 1, nn
		exactScalar(i,j) = Gaussian( [x(j), y(i)], b)
		xVec(1:2) = GaussGrad([x(j), y(i)], b)
		exactGradX(i,j) = xVec(1)
		exactGradY(i,j) = xVec(2)
		exactLap(i,j) = GaussLap([x(j), y(i)], b)
		
		interpolatedScalar(i,j) = InterpolateScalar( gaussianInterp, triMesh, [x(j), y(i)] )
		xVec(1:2) = InterpolateGradient( gaussianInterp, triMesh, [x(j), y(i)])
		interpGradX(i,j) = xVec(1)
		interpGradY(i,j) = xVec(2)
		interpLap(i,j) = InterpolateLaplacian(gaussianInterp, triMesh, [x(j),y(i)])
	enddo
enddo


write(filename,'(A,I1,A)') "triScatteredInterpTest_", initNest, ".m"
open(unit=WRITE_UNIT_1, file=filename, status="REPLACE", action="WRITE")

call WriteParticlesToMatlab(triMesh%particles, WRITE_UNIT_1)
call WriteFieldToMatlab(scalarField, WRITE_UNIT_1)
call WriteFieldToMatlab(gradientField, WRITE_UNIT_1)
call WriteFieldToMatlab(secondPartials, WRITE_UNIT_1)
call WriteFieldToMatlab(gradError, WRITE_UNIT_1)
call WriteFieldToMatlab(estLaplacian, WRITE_UNIT_1)
call WriteFieldToMatlab(lapError, WRITE_UNIT_1)

write(WRITE_UNIT_1,'(A,F8.2,A,F12.8,A,F8.2,A)') "unifx = ", xmin,":", dx,":", xmax, ";"
write(WRITE_UNIT_1,'(A)') "unify = unifx;"
write(WRITE_UNIT_1,'(A)') "[xx, yy] = meshgrid(unifx, unify);"

call WriteMeshToMatlab(triMesh, WRITE_UNIT_1)

write(WRITE_UNIT_1,'(A)',advance='NO') "exactScalar = ["
do i = 1, nn - 1
	do j = 1, nn - 1
		write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') exactScalar(i,j), ","
	enddo
	write(WRITE_UNIT_1,'(F12.7,A)') exactScalar(i,nn), "; ..."
enddo
do j = 1, nn - 1
	write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') exactScalar(nn,j), ","
enddo
write(WRITE_UNIT_1,'(F12.7,A)') exactScalar(nn,nn), "];"

write(WRITE_UNIT_1,'(A)',advance='NO') "exactGradX = ["
do i = 1, nn - 1
	do j = 1, nn - 1
		write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') exactGradX(i,j), ","
	enddo
	write(WRITE_UNIT_1,'(F12.7,A)') exactGradX(i,nn), "; ..."
enddo
do j = 1, nn - 1
	write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') exactGradX(nn,j), ","
enddo
write(WRITE_UNIT_1,'(F12.7,A)') exactGradX(nn,nn), "];"

write(WRITE_UNIT_1,'(A)',advance='NO') "exactGradY = ["
do i = 1, nn - 1
	do j = 1, nn - 1
		write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') exactGradY(i,j), ","
	enddo
	write(WRITE_UNIT_1,'(F12.7,A)') exactGradY(i,nn), "; ..."
enddo
do j = 1, nn - 1
	write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') exactGradY(nn,j), ","
enddo
write(WRITE_UNIT_1,'(F12.7,A)') exactGradY(nn,nn), "];"

write(WRITE_UNIT_1,'(A)',advance='NO') "exactLap = ["
do i = 1, nn - 1
	do j = 1, nn - 1
		write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') exactLap(i,j), ","
	enddo
	write(WRITE_UNIT_1,'(F12.7,A)') exactLap(i,nn), "; ..."
enddo
do j = 1, nn - 1
	write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') exactLap(nn,j), ","
enddo
write(WRITE_UNIT_1,'(F12.7,A)') exactLap(nn,nn), "];"

write(WRITE_UNIT_1,'(A)',advance='NO') "interpolatedScalar = ["
do i = 1, nn - 1
	do j = 1, nn - 1
		write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') interpolatedScalar(i,j), ","
	enddo
	write(WRITE_UNIT_1,'(F12.7,A)') interpolatedScalar(i,nn), "; ..."
enddo
do j = 1, nn - 1
	write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') interpolatedScalar(nn,j), ","
enddo
write(WRITE_UNIT_1,'(F12.7,A)') interpolatedScalar(nn,nn), "];"

write(WRITE_UNIT_1,'(A)',advance='NO') "interpGradX = ["
do i = 1, nn - 1
	do j = 1, nn - 1
		write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') interpGradX(i,j), ","
	enddo
	write(WRITE_UNIT_1,'(F12.7,A)') interpGradX(i,nn), "; ..."
enddo
do j = 1, nn - 1
	write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') interpGradX(nn,j), ","
enddo
write(WRITE_UNIT_1,'(F12.7,A)') interpGradX(nn,nn), "];"

write(WRITE_UNIT_1,'(A)',advance='NO') "interpGradY = ["
do i = 1, nn - 1
	do j = 1, nn - 1
		write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') interpGradY(i,j), ","
	enddo
	write(WRITE_UNIT_1,'(F12.7,A)') interpGradY(i,nn), "; ..."
enddo
do j = 1, nn - 1
	write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') interpGradY(nn,j), ","
enddo
write(WRITE_UNIT_1,'(F12.7,A)') interpGradY(nn,nn), "];"

write(WRITE_UNIT_1,'(A)',advance='NO') "interpLap = ["
do i = 1, nn - 1
	do j = 1, nn - 1
		write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') interpLap(i,j), ","
	enddo
	write(WRITE_UNIT_1,'(F12.7,A)') interpLap(i,nn), "; ..."
enddo
do j = 1, nn - 1
	write(WRITE_UNIT_1, '(F12.7,A)',advance='NO') interpLap(nn,j), ","
enddo
write(WRITE_UNIT_1,'(F12.7,A)') interpLap(nn,nn), "];"

close(WRITE_UNIT_1)

!
!----------------
! PROGRAM END
!----------------
!
call Delete(lapError)
call Delete(estLaplacian)
call Delete(exactGrad)
call Delete(exactPartials)
call Delete(gaussianInterp)
call Delete(gradError)
call Delete(secondPartials)
call Delete(gradientField)
call Delete(scalarField)
call Delete(triMesh)
call Delete(exeLog)

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