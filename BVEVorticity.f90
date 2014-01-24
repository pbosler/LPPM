module BVESetupModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	Defines vorticity distributions used by SphereMesh objects.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
!----------------
! USAGE :  This module provides methods for integrating the barotropic vorticity equation on the sphere.
!----------------
use NumberKindsModule
use LoggerModule
use SphereGeomModule
use ParticlesModule
use EdgesModule
use PanelsModule
use SphereMeshModule

implicit none

include 'mpif.h'

private
public BVESetup
public New, Delete
public InitSolidBodyRotation, SetSolidBodyRotationOnMesh
public SOLID_BODY_NINT, SOLID_BODY_NREAL
public NullVorticity

!
!----------------
! Types and module constants
!----------------
!
type BVESetup
	integer(kint), pointer :: integers(:) => null()	!	integers required by test case definitions
	real(kreal), pointer :: reals(:) => null() 		! 	real numbers required by test case definitions
end type

integer(kint), parameter :: SOLID_BODY_NINT = 0, SOLID_BODY_NREAL = 1

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'BVESetup'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=128) :: logString
character(len=24) :: formatString
!
!----------------
! Interfaces
!----------------
!
interface New
	module procedure NewPrivate
	module procedure NewPrivateNull
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
	type(BVESetup), intent(out) :: self
	integer(kint), intent(in) :: nINt, nReal

	if (.NOT. logInit) call InitLogger(log,procRank)

	if ( nInt > 0 ) then
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

subroutine NewPrivateNull(self)
	type(BVESetup), intent(out) :: self
	if ( .NOT. logInit) call InitLogger(log,procRank)
	nullify(self%integers)
	nullify(self%reals)
end subroutine

subroutine DeletePrivate(self)
	type(BVESetup), intent(inout) :: self
	if (associated(self%integers))	deallocate(self%integers)
	if (associated(self%reals))	deallocate(self%reals)
	call Delete(log)
	logInit = .FALSE.
end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine InitSolidBodyRotation(solidBody, rotRate)
	type(BVESetup), intent(inout) :: solidBody
	real(kreal), intent(in) :: rotRate

	if ( size(solidBody%reals) /= 1) then
		call LogMessage(log, ERROR_LOGGING_LEVEL,'BVESetup ERROR : ',' real array size incorrect')
		return
	endif

	solidBody%reals(1) = rotRate
end subroutine

subroutine SetSolidBodyRotationOnMesh(aMesh, solidBody)
! sets absolute and relative vorticity for solid body test case
! NOTE : the role of absolute and relative vorticity is reversed for this test case only.
	type(SphereMesh), intent(inout) :: aMesh
	type(BVESetup), intent(in) :: solidBody
	integer(kint) :: j
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	do j=1,aParticles%N
		aParticles%relVort(j) = SolidBodyX(aParticles%x0(:,j), solidBody%reals(1))
	enddo
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = SolidBodyX(aPanels%x0(:,j), solidBody%reals(1))
		else
			aPanels%relVort(j) = 0.0_kreal
		endif
	enddo

end subroutine

subroutine NullVorticity(aMesh,nullVort)
	type(SphereMesh), intent(inout) :: aMesh
	type(BVESetup), intent(in) :: nullVort
end subroutine

!
!----------------
! Module methods : module- or type-specific private functions
!----------------
!

function SolidBodyX(xyz,rotationRate)
! Outputs vorticity associated with solid body rotation with angular velocity = rotation rate
! at position xyz.
	real(kreal) :: SolidBodyX
	real(kreal), intent(in) :: xyz(3), rotationRate
	SolidBodyX = 2.0_kreal*rotationRate*xyz(3)/EARTH_RADIUS
end function


subroutine InitLogger(aLog,rank)
	type(Logger), intent(inout) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

end module
