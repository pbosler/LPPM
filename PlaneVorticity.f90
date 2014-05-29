module PlaneVorticityModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	Defines the RK4 data structure used by SphereMesh.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
!----------------
! USAGE :  This module provides methods for integrating the barotropic vorticity equation on the sphere.
!----------------
use NumberKindsModule
use LoggerModule
use PlaneGeomModule
use ParticlesModule
use EdgesModule
use PanelsModule
use PlaneMeshModule

implicit none

private
public VorticitySetup
public New, Delete
public InitRankineVortex, SetRankineVortexOnMesh, RankineVortex, RANKINE_N_INT, RANKINE_N_REAL
public NullVorticity

!
!----------------
! Types and module constants
!----------------
!
type VorticitySetup
	integer(kint), pointer :: integers(:) => null() 	! integers required by test case definitions
	real(kreal), pointer :: reals(:) => null()			! reals required by test case definitions
end type

integer(kint), parameter :: RANKINE_N_INT = 0, RANKINE_N_REAL = 4

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'VorticitySetup'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=128) :: logString
character(len=24) :: formatString
!
!----------------
! Interfaces
!----------------
!
!
!----------------
! Interfaces
!----------------
!
interface New
	module procedure NewPrivate
end interface

interface Delete
	module procedure DeletePrivate
end interface

contains
!
!----------------
! Standard methods : Constructor / Destructor
!----------------
!

subroutine NewPrivate(self, nInt, nReal)
	type(VorticitySetup), intent(out) :: self
	integer(kint), intent(in) :: nInt, nReal

	if ( .not. loginit) call InitLogger(log,procRank)

	if (nInt > 0 ) then
		allocate(self%integers(nInt))
		self%integers = 0
	else
		nullify(self%integers)
	endif

	if ( nReal > 0 ) then
		allocate(self%reals(nReal))
		self%reals = 0.0_kreal
	else
		nullify(self%reals)
	endif
end subroutine

subroutine DeletePrivate(self)
	type(VorticitySetup), intent(inout) :: self
	if ( associated(self%integers)) deallocate(self%integers)
	if ( associated(self%reals)) deallocate(self%reals)
	nullify(self%integers)
	nullify(self%reals)
end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine NullVorticity(aMesh, self)
	type(PlaneMesh), intent(inout) :: aMesh
	type(VorticitySetup), intent(in) :: self
end subroutine


subroutine InitRankineVortex(self, xCent, yCent, radius, strength)
	type(VorticitySetup), intent(inout) :: self
	real(kreal), intent(in) :: xCent, yCent, radius, strength

	if ( size(self%reals) /= 4 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'VorticitySetup ERROR : ',' real array size incorrect.')
		return
	endif

	self%reals(1) = xCent
	self%reals(2) = yCent
	self%reals(3) = radius
	self%reals(4) = strength
end subroutine

subroutine SetRankineVortexOnMesh(aMesh, self)
	type(PlaneMesh), intent(inout) :: aMesh
	type(VorticitySetup), intent(in) :: self
	!
	type(Particles), pointer :: aparticles
	type(Panels), pointer :: apanels
	integer(kint) :: j
	real(kreal) :: xyCent(2), radius, strength

	aparticles => aMesh%particles
	apanels => aMesh%panels

	xyCent = [ self%reals(1), self%reals(2)]
	radius = self%reals(3)
	strength = self%reals(4)

	do j=1,aParticles%N
		aParticles%relVort(j) = RankineVortex(aParticles%x0(:,j), xyCent, radius, strength)
	enddo
	do j=1,aPanels%N
		if (aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = 0.0_kreal
		else
			aPanels%relVort(j) = RankineVortex(aPanels%x0(:,j), xyCent, radius, strength)
		endif
	enddo

end subroutine


function RankineVortex(xy, xyCent, radius, strength)
	real(kreal) :: RankineVortex
	real(kreal), intent(in) :: xy(2), xyCent(2), radius, strength
	!
	real(kreal) :: r

	r = sqrt( sum( (xy - xyCent) * (xy - xyCent) ) )

	RankineVortex = 0.0_kreal

	if ( r <= radius) RankineVortex = strength
end function

!
!----------------
! Module methods : module- or type-specific private functions
!----------------
!


subroutine InitLogger(aLog, rank)
	type(Logger), intent(inout) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	if ( rank == 0 ) then
		call New(aLog,logLevel)
	else
		call New(aLog,WARNING_LOGGING_LEVEL)
	endif
	logInit = .TRUE.
end subroutine


end module
