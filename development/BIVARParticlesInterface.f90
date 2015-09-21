module BIVARInterfaceModule

use NumberKindsModule
use LoggerModule
use ParticlesModule
use FieldModule
use BIVARModule

implicit none

private
public BIVARInterface, New, Delete
public InterpolateScalar, InterpolateVector

!
!----------------
! Types and module constants
!----------------
!
type BIVARInterface
	real(kreal), dimension(:), pointer :: estPartials1 => null()
	real(kreal), dimension(:), pointer :: estPartials2 => null()
end type

type PlaneTri
	integer(kint) :: nTri
	integer(kint), pointer, dimension(:) :: vertIndices => null()
	integer(kint) :: nBorderSegments
	integer(kint), pointer, dimension(:) :: borderIndices => null()
	integer(kint), pointer, dimension(:) :: integerWork => null()
	real(kreal), pointer, dimension(:) :: realWork => null()
	contains
		final :: deleteTri
end type

interface New
	module procedure newPrivateField
	module procedure newPrivateLagParam
	module procedure newTri
end interface

interface Delete
	module procedure deletePrivate
	module procedure deleteTri
end interface

contains

subroutine newTri( delTri, aParticles )
	type(PlaneTri), intent(out) :: delTri
	type(Particles), intent(in) :: aParticles
	!
	integer(kint), allocatable, dimension(:) :: iwp

	allocate(delTri%vertIndices( 6 * aParticles%N - 15) )
	allocate(delTri%borderIndices(6* aParticles%N) )
	allocate(delTri%integerWork( 18 * aParticles%N))
	allocate(delTri%realWork(aParticles%N))
	
	allocate(iwp( aParticles%N))
	
	call IDTANG( aParticles%N, aParticles%x, aParticles%y, delTri%nTri, delTri%vertIndices, delTri%nBorderSegments, &
				 delTri%borderIndices, delTri%integerWork, iwp, delTri%realWork)
	
	deallocate(iwp)
end subroutine

subroutine deleteTri(delTri)
	type(PlaneTri), intent(inout) :: delTri
	if ( associated(delTri%vertIndices)) deallocate(delTri%vertIndices)
	if ( associated(delTri%borderIndices)) deallocate(delTri%borderIndices)
	if ( associated(delTri%integerWork)) deallocate(delTri%integerWork)
	if ( associated(delTri%realWork)) deallocate(delTri%realWork)
end subroutine

subroutine newPrivateField( self, delTri, sourceParticles, sourceField )
	type(BIVARInterface), intent(out) :: self
	type(PlaneTri), intent(in) :: delTri
	type(Particles), intent(in) :: sourceParticles
	type(Field), intent(in) :: sourceField
	!
	real(kreal), allocatable, dimension(:) :: realWork

	allocate(self%estPartials1( 5 * sourceParticles%N ))
	if ( sourceField%nDim == 2 ) allocate(self%estPartials2( 5 * sourceParticles%N))
	
	allocate(realWork(sourceParticles%N))
	
	if ( sourceField%nDim == 1 ) then
		call IDPDRV( sourceParticles%N, sourceParticles%x, sourceParticles%y, sourceField%scalar, &
			delTri%nTri, delTri%vertIndices, self%estPartials1, realWork)
	elseif ( sourceField%nDim == 2 ) then
		call IDPDRV( sourceParticles%N, sourceParticles%x, sourceParticles%y, sourceField%xComp, &
			delTri%nTri, delTri%vertIndices, self%estPartials1, realWork)
		call IDPDRV( sourceParticles%N, sourceParticles%x, sourceParticles%y, sourceField%yComp, &
			delTri%nTri, delTri%vertIndices, self%estPartials1, realWork)	
	endif
	
	deallocate(realWork)
end subroutine

subroutine newPrivateLagParam( self, delTri, sourceParticles )
	type(BIVARInterface), intent(out) :: self
	type(PlaneTri), intent(in) :: delTri
	type(Particles), intent(in) :: sourceParticles
	!
	real(kreal), allocatable, dimension(:) :: realWork

	allocate(self%estPartials1( 5 * sourceParticles%N ))
	allocate(self%estPartials2( 5 * sourceParticles%N ))
	
	allocate(realWork(sourceParticles%N))
	
	call IDPDRV( sourceParticles%N, sourceParticles%x, sourceParticles%y, sourceParticles%x0, &
			delTri%nTri, delTri%vertIndices, self%estPartials1, realWork)
	call IDPDRV( sourceParticles%N, sourceParticles%x, sourceParticles%y, sourceParticles%y0, &
			delTri%nTri, delTri%vertIndices, self%estPartials1, realWork)	

	deallocate(realWork)
end subroutine

subroutine deletePrivate(self)
	type(BIVARInterface), intent(inout) :: self
	if ( associated(self%estPartials1)) deallocate(self%estPartials1)
	if ( associated(self%estPartials2)) deallocate(self%estPartials2)
end subroutine

function InterpolateScalar( self, delTri, sourceParticles, sourceField, xOut, yOut )
	real(kreal) :: InterpolateScalar 
	type(BIVARInterface), intent(in) :: self
	type(PlaneTri), intent(in) :: delTri
	type(Particles), intent(in) :: sourceParticles
	type(Field), intent(in) :: sourceField
	real(kreal), intent(in) :: xOut
	real(kreal), intent(in) :: yOut
	!
	integer(kint) :: inTri
	
	call IDLCTN( sourceParticles%N, sourceParticles%x, sourceParticles%y, delTri%nTri, delTri%vertIndices, &
				 delTri%nBorderSegments, delTri%borderIndices, xOut, yOut, inTri, delTri%integerWork, delTri%realWork)
	
	call IDPTIP( sourceParticles%N, sourceParticles%x, sourceParticles%y, sourceField%scalar, &
		delTri%nTri, delTri%vertIndices, delTri%nBorderSegments, delTri%borderIndices, self%estPartials1, inTri, &
		xOut, yOut, InterpolateScalar)
end function

function InterpolateVector( self, delTri, sourceParticles, sourceField, xOut, yOut )
	real(kreal), dimension(2) :: InterpolateVector 
	type(BIVARInterface), intent(in) :: self
	type(PlaneTri), intent(in) :: delTri
	type(Particles), intent(in) :: sourceParticles
	type(Field), intent(in) :: sourceField
	real(kreal), intent(in) :: xOut
	real(kreal), intent(in) :: yOut
	!
	integer(kint) :: inTri
	
	call IDLCTN( sourceParticles%N, sourceParticles%x, sourceParticles%y, delTri%nTri, delTri%vertIndices, &
				 delTri%nBorderSegments, delTri%borderIndices, xOut, yOut, inTri, delTri%integerWork, delTri%realWork)

	call IDPTIP( sourceParticles%N, sourceParticles%x, sourceParticles%y, sourceField%xComp, &
		delTri%nTri, delTri%vertIndices, delTri%nBorderSegments, delTri%borderIndices, self%estPartials1, inTri, &
		xOut, yOut, InterpolateVector(1))
	call IDPTIP( sourceParticles%N, sourceParticles%x, sourceParticles%y, sourceField%yComp, &
		delTri%nTri, delTri%vertIndices, delTri%nBorderSegments, delTri%borderIndices, self%estPartials2, inTri, &
		xOut, yOut, InterpolateVector(2))
end function		

end module
