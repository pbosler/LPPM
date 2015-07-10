program BIVARConvergenceTest

use NumberKindsModule
use LoggerModule
use PolyMesh2dModule
use ParticlesModule
use EdgesModule
use FacesModule
use STDIntVectorModule
use FieldModule
use BIVARModule

implicit none

type(Logger) :: exeLog
character(len=MAX_STRING_LENGTH) :: logstring

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

! BIVAR
integer(kint) :: nSourcePoints, nTriangles, nBoundarySegments, 
integer(kint), allocatable :: triVerts(:), boundaryEdgesAndTris(:)
integer(kint), allocatable :: intWork1(:), intWork2(:)
real(kreal), allocatable :: realWork(:), realWork2(:)
integer(kint) :: triIndices(nn,nn)
real(kreal), allocatable :: pdd(:)


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
	! use BIVAR to interpolate and estimate derivatives
	!
	!		1. Triangulate the particle set
	allocate(triVerts(6 * triMesh%particles%N - 15))
	allocate(boundaryEdgesAndTris(6*triMesh%particles%N))
	allocate(intWork1(18*triMesh%particles%N))
	allocate(intWork2(triMesh%particles%N))
	allocate(realWork(triMesh%particles%N))
	
	call idtang( triMesh%particles%N, triMesh%particles%x, triMesh%particles%y, nTriangles, triVerts, &
				 nBoundarySegments, boundaryEdgesAndTris, intWork1, intWork2, realWork)
	!		2. Locate output points within the triangulation
	allocate(realWork2(8*triMesh%particles%N))
	do j = 1, nn
		do i = 1, nn
			call idlctn( triMesh%particles%N, triMesh%particles%x, triMesh%particles%y, nTriangles, triVerts, &
				 nBoundarySegments, boundaryEdgesAndTris, x(j), y(i), triIndices(i, j), intWork1, realWork2)
		enddo
	enddo
	!		3. Estimate partial derivatives at each particles
	allocate(pdd(5*triMesh%particles%N))
	call idpdrv( triMesh%particles%N, triMesh%particles%x, triMesh%particles%y, scalar%scalar,&
			     nTriangles, triVerts, pdd, realWork)
	do i = 1, triMesh%particles%N
		call InsertVectorToField(estGrad, [pdd(5*i - 4), pdd(5*i-3)])
		call InsertVectorToField(est2ndPartials, [pdd(5*i-2), pdd(5*i-1), pdd(5*i)]
	enddo
	!	    4. Interpolate the scalar
	do j = 1, nn
		do i = 1, nn
			call idptip( triMesh%particles%N, triMesh%particles%x, triMesH%particles%y, scalar%scalar, &
						  nTriangles, triVerts, nBoundarySegments, boundaryEdgesAndTris, pdd, &
						  triIndices(i,j), x(j), y(i), interpScalar(i,j))
		enddo
	enddo
	
	
	deallocate(pdd)
	deallocate(realWork2)
	deallocate(triVerts)
	deallocate(boundaryEdgesAndTris)
	deallocate(intWork1)
	deallocate(intWork2)
	deallocate(realWork)
	!
	! setup for next iteration
	!
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