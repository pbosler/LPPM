module PSEModule

use NumberKindsModule
use LoggerModule
use PlaneGeomModule
use SphereGeomModule
use ParticlesModule
use EdgesModule
use FacesModule
use PolyMesh2dModule
use FieldModule
use STDIntVectorModule

implicit none
private
public PSE, New, Delete
public PSEGradientAtParticles, PSELaplacianAtParticles, PSESecondPartialsAtParticles, PSEInterpolateScalar

type PSE
	real(kreal) :: eps ! kernel radius
	
	contains
		final :: deletePrivate
end type

interface New
	module procedure :: newPrivate
end interface

interface Delete
	module procedure :: deletePrivate
end interface

contains

subroutine newPrivate(self, aMesh, kernelRadiusMultiple )
	type(PSE), intent(out) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	real(kreal), intent(in), optional :: kernelRadiusMultiple
	!
	real(kreal) :: multiple
	
	if ( present(kernelRadiusMultiple) ) then
		multiple = kernelRadiusMultiple
	else
		multiple = 2.0_kreal
	endif
	
	call setKernelRadius(self, aMesh, multiple)
end subroutine

subroutine deletePrivate(self)
	type(PSE), intent(inout) :: self
end subroutine

subroutine setKernelRadius( self, aMesh, multiple )
	type(PSE), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	real(kreal), intent(in) :: multiple
	self%eps = multiple * MaxEdgeLength(aMesh%edges, aMesh%particles)
end subroutine

subroutine PSEGradientAtParticles( self, aMesh, scalarField, scalarGrad )
	type(PSE), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: scalarField
	type(Field), intent(inout) :: scalarGrad
	!procedure, intent(in) :: gradKernel
	!
	integer(kint) :: i, j
	real(kreal) :: dist, kernelIn, xi(3), xj(3)
	
	call SetFieldToZero(scalarGrad)
	scalarGrad%N = aMesh%particles%N
	
	do i = 1, aMesh%particles%N
		xi = PhysCoord(aMesh%particles, i)
		do j = 1, aMesh%particles%N	
			if ( aMesh%particles%isActive(j) ) then
				xj = PhysCoord(aMesh%particles, j)
				dist = sqrt( sum( (xj - xi)*(xj - xi) ))
				kernelIn = dist / self%eps
				scalarGrad%xComp(i) = scalarGrad%xComp(i) + (scalarField%scalar(j) + scalarField%scalar(i)) * &
					(xi(1)-xj(1)) * bivariateFirstDerivativeKernel8( kernelIn )/(self%eps**3) * aMesh%particles%area(j)
				scalarGrad%yComp(i) = scalarGrad%yComp(i) + (scalarField%scalar(j) + scalarField%scalar(i)) * &
					(xi(2)-xj(2)) * bivariateFirstDerivativeKernel8( kernelIn )/(self%eps**3) * aMesh%particles%area(j)	
			endif
		enddo
	enddo	
	scalarGrad%xComp = scalarGrad%xComp / self%eps
	scalarGrad%yComp = scalarGrad%yComp / self%eps
end subroutine

subroutine PSESecondPartialsAtParticles( self, aMesh, scalarGrad, secondPartials)
	type(PSE), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: scalarGrad
	type(Field), intent(inout) :: secondPartials
	!
	integer(kint) :: i, j
	real(kreal) :: dist, kernelIn, xi(3), xj(3), dxx, dxy, dyx, dyy
	
	call SetFieldToZero(secondPartials)
	secondPartials%N = aMesh%particles%N
	
	do i = 1, aMesh%particles%N
		xi = PhysCoord(aMesh%particles, i)
		dxx = 0.0_kreal
		dxy = 0.0_kreal
		dyx = 0.0_kreal
		dyy = 0.0_kreal
		do j = 1, aMesh%particles%N
			if ( aMesh%particles%isActive(j) ) then
				xj = PhysCoord(aMesh%particles, j)
				dist = sqrt(sum( (xj-xi)*(xj-xi) ))
				kernelIn = dist / self%eps
				dxx = dxx + (scalarGrad%xComp(j) + scalarGrad%xComp(i)) * ( xi(1)-xj(1) ) * &
					bivariateFirstDerivativeKernel8( kernelIn )/(self%eps**3) * aMesh%particles%area(j)
				dxy = dxy + (scalarGrad%xComp(j) + scalarGrad%xComp(i)) * ( xi(2)-xj(2) ) * &
					bivariateFirstDerivativeKernel8( kernelIn )/(self%eps**3) * aMesh%particles%area(j)
				dyx = dyx + (scalarGrad%yComp(j) + scalarGrad%yComp(i)) * ( xi(1)-xj(1) ) * &
					bivariateFirstDerivativeKernel8( kernelIn )/(self%eps**3) * aMesh%particles%area(j)
				dyy = dyy + (scalarGrad%yComp(j) + scalarGrad%yComp(i)) * ( xi(2)-xj(2) ) * &
					bivariateFirstDerivativeKernel8( kernelIn )/(self%eps**3) * aMesh%particles%area(j)
			endif
		enddo
		secondPartials%xComp(i) = dxx / self%eps
		secondPartials%yComp(i) = 0.5_kreal * ( dxy + dyx ) / self%eps
		secondPartials%zComp(i) = dyy / self%eps
	enddo
end subroutine

subroutine PSELaplacianAtParticles( self, aMesh, scalarField, scalarLap )
	type(PSE), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: scalarField
	type(Field), intent(inout) :: scalarLap
	!
	integer(kint) :: i, j
	real(kreal) :: dist, kernelIn, xi(3), xj(3)
	
	call SetFieldToZero(scalarLap)
	scalarLap%N = aMesh%particles%N
	
	do i = 1, aMesh%particles%N
		xi = PhysCoord(aMesh%particles, i)
		do j = 1, aMesh%particles%N	
			if ( aMesh%particles%isActive(j) ) then
				xj = PhysCoord(aMesh%particles, j)
				dist = sqrt( sum( (xj - xi)*(xj - xi) ))
				kernelIn = dist / self%eps
				scalarLap%scalar(i) = scalarLap%scalar(i) + (scalarField%scalar(j) - scalarField%scalar(i)) * &
				  bivariateLaplacianKernel8( kernelIn )/(self%eps**2) * aMesh%particles%area(j)
			endif
		enddo
	enddo		
	scalarLap%scalar = scalarLap%scalar / (self%eps**2)
end subroutine

function PSEInterpolateScalar( self, aMesh, scalarField, xOut, yOut )
	real(kreal) :: PSEInterpolateScalar
	type(PSE), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: scalarField
	real(kreal), intent(in) :: xOut, yOut
	!
	integer(kint) :: j
	real(kreal) :: dist, kernelIn, xi(3), xj(3)
	
	PSEInterpolateScalar = 0.0_kreal
	xi(1) = xOut
	xi(2) = yOut
	xi(3) = 0.0_kreal
	do j = 1, aMesh%particles%N
		if ( aMesh%particles%isActive(j) ) then
			xj = PhysCoord(aMesh%particles, j)
			dist = sqrt( sum( (xj-xi)*(xj-xi)))
			kernelIn = dist / self%eps
			PSEInterpolateScalar = PSEInterpolateScalar + scalarField%scalar(j) * aMesh%particles%area(j) * &
			bivariateDeltaKernel8( kernelIn ) / (self%eps**2)
		endif
	enddo
end function

pure function bivariateFirstDerivativeKernel8( radialDist)
	real(kreal) :: bivariateFirstDerivativeKernel8
	real(kreal), intent(in) :: radialDist
	bivariateFirstDerivativeKernel8 = (-20.0_kreal + 20.0_kreal * radialDist**2 - 5.0_kreal * radialDist**4 + &
									   radialDist**6/3.0_kreal) * exp( - radialDist*radialDist ) / PI
end function

pure function bivariateLaplacianKernel8( radialDist )
	real(kreal) :: bivariateLaplacianKernel8
	real(kreal), intent(in) :: radialDist
	bivariateLaplacianKernel8 = (40.0_kreal - 40.0_kreal * radialDist**2 + 10.0_kreal * radialDist**4 - &
								 2.0_kreal * radialDist**6/3.0_kreal) * exp( -radialDist*radialDist) / PI
end function

pure function bivariateDeltaKernel4( radialDist )
	real(kreal) :: bivariateDeltaKernel4
	real(kreal), intent(in) :: radialDist
!	real(kreal), parameter :: a = 2.0_kreal
!	real(kreal), parameter :: c1 = ( a**2 / PI - 1.0_kreal /( a**2 * PI )) * 1.0_kreal / ( a**2 - 1) 
!	real(kreal), parameter :: c2 = ( 1.0_kreal / (a**2 * PI) - 1.0_kreal / PI) * 1.0_kreal / (a**4 - a**2)
!	 
!	bivariateDeltaKernel4 = c1 * exp( - radialDist * radialDist ) + c2 * exp( - radialDist * radialDist / a / a )
	bivariateDeltaKernel4 = (2.0_kreal - radialDist * radialDist ) * exp( -radialDist * radialDist )/PI
end function

pure function bivariateDeltaKernel8( radialDist ) 
	real(kreal) :: bivariateDeltaKernel8
	real(kreal), intent(in) :: radialDist
	bivariateDeltaKernel8 = (4.0_kreal - 6.0_kreal * radialDist**2 + 2.0_kreal * radialDist**4 - &
			 radialDist**6 / 6.0_kreal) * exp(-radialDist * radialDist) / PI
end function

end module
