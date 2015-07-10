module BIVARInterfaceModule

use NumberKindsModule
use LoggerModule
use ParticlesModule
use FieldModule
use BIVARModule

implicit none

private
public BIVARInterface, New, Delete, Copy
public InterpolateScalarField, InterpolateVectorField
public EstimateDerivatives

!
!----------------
! Types and module constants
!----------------
!
type BIVARInterface
	integer(kint) :: nSourcePoints = 0
	integer(kint) :: nTri = 0
	integer(kint), pointer :: triVerts(:) => null()
	integer(kint) :: nBoundaryTriEdges = 0
	integer(kint), pointer :: boundaryEdgesAndTri(:) => null()
	real(kreal), pointer :: scalarPartials(:) => null()
end type

interface New
	module procedure newPrivate
end interface

interface Copy
	module procedure copyPrivate
end interface

interface Delete
	module procedure deletePrivate
end interface

contains

subroutine newPrivate( self, sourceParticles, sourceField )
	type(BIVARInterface), intent(out) :: self
	type(Particles), intent(in) :: sourceParticles
	type(Field), intent(in) :: sourceField
	
end subroutine

subroutine InterpolateScalarField( aParticles, scalarField, xOut, yOut, interpOut )
	type(Particles), intent(in) :: aParticles
	type(Field), intent(in) :: scalarField
	real(kreal), intent(in) :: xOut(:)
	real(kreal), intent(in) :: yOut(:)
	real(kreal), intent(inout) :: interpOut(:)
	
	
end subroutine

end module
