module SWESetupModule

use NumberKindsModule
use LoggerModule
use SphereGeomModule
use ParticlesModule
use PanelsModule
use SphereMeshModule

implicit none

private
public SWESetup
public New, Delete
public InitNonlinearSteadyStateTC2, SetNonlinearSteadyStateTC2OnMesh, TC2_NREAL, TC2_NINT

!
!----------------
! Types and module constants
!----------------
!
type SWESetup
	integer(kint), pointer :: integers(:) => null()	!	integers required by test case definitions
	real(kreal), pointer :: reals(:) => null() 		! 	real numbers required by test case definitions
end type

integer(kint), parameter :: TC2_NREAL = 3, &
							TC2_NINT = 0 


!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'SWESetup'
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
	type(SWESetup), intent(out) :: self
	integer(kint), intent(in) :: nInt, nReal
	
	if (.NOT. logInit) call InitLogger(log,procRank)
	
	if ( nInt > 0 ) then
		allocate(self%integers(nInt))
		self%integers = 0
	else
		nullify(self%integers)
	endif
	if ( nReal > 0 ) then
		allocate(self%reals(nReal))
		self%real = 0.0_kreal
	else
		nullify(self%reals)
	endif
end subroutine

subroutine DeletePrivate(self)
	type(SWESetup), intent(inout) :: self
	if ( associated( self%integers)) deallocate(self%integers)
	if ( associated(self%reals)) deallocate(self%reals)
end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine InitNonlinearSteadyStateTC2(self, u0, rotAlpha, h0)
	type(SWESetup), intent(inout) :: self
	real(kreal), intent(in) :: u0, rotAlpha, h0
	
	if ( rotAlpha /= 0.0_kreal ) then
		call LogMessage(log, WARNING_LOGGING_LEVEL,logKey,' WARNING : rotated coriolis not implemented in timestepping yet.')
	endif
	
	self%reals(1) = u0
	self%reals(2) = rotAlpha
	self%reals(3) = h0	
end subroutine

subroutine SetNonlinearSteadyStateTC2OnMesh(aMesh, testCase2)
	type(SphereMesh), intent(inout) :: aMesh
	type(SWESetup), intent(in) :: testCase2
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j
	real(kreal) :: u0, rotAlpha, h0
	
	if ( aMesh%problemKind /= SWE_SOLVER ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,' problemKind must be SWE_SOLVER.')
		return
	endif
	
	aParticles => aMesh%particles
	aPanels => aMesh%panels
	
	u0 = testCase2%reals(1)
	rotAlpha = testCase2%reals(2)
	h0 = testCase2%reals(3)
	
	do j=1, aParticles%N
		aParticles%relVort(j) = 2.0_kreal*u0/EARTH_RADIUS * &
			( -sin(rotAlpha)*aParticles%x0(1,j)/EARTH_RADIUS + cos(rotAlpha)*aParticles%x0(3,j)/EARTH_RADIUS )
		aParticles%div(j) = 0.0_kreal
		aParticles%h(j) = GRAV*h0 - (EARTH_RADIUS*OMEGA*u0 + u0*u0/2.0_kreal) * &
			( -sin(rotAlpha)*aParticles%x0(1,j)/EARTH_RADIUS + cos(rotAlpha)*aParticles%x0(3,j)/EARTH_RADIUS)*&
			( -sin(rotAlpha)*aParticles%x0(1,j)/EARTH_RADIUS + cos(rotAlpha)*aParticles%x0(3,j)/EARTH_RADIUS)
		aParticles%potVort(j) = (aParticles%relVort(j) + 2.0_kreal*OMEGA*&
			( -sin(rotAlpha)*aParticles%x0(1,j)/EARTH_RADIUS + cos(rotAlpha)*aParticles%x0(3,j)/EARTH_RADIUS))/&
			aParticles%h(j)
	enddo
	do j=1, aPanels%N
		if (.NOT. aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = 2.0_kreal*u0/EARTH_RADIUS * &
				( -sin(rotAlpha)*aPanels%x0(1,j)/EARTH_RADIUS + cos(rotAlpha)*aPanels%x0(3,j)/EARTH_RADIUS )
			aPanels%div(j) = 0.0_kreal
			aPanels%h(j) = GRAV*h0 - (EARTH_RADIUS*OMEGA*u0 + u0*u0/2.0_kreal) * &
				( -sin(rotAlpha)*aPanels%x0(1,j)/EARTH_RADIUS + cos(rotAlpha)*aPanels%x0(3,j)/EARTH_RADIUS)*&
				( -sin(rotAlpha)*aPanels%x0(1,j)/EARTH_RADIUS + cos(rotAlpha)*aPanels%x0(3,j)/EARTH_RADIUS)
			aPanels%potVort(j) = (aPanels%relVort(j) + 2.0_kreal*OMEGA*&
				( -sin(rotAlpha)*aPanels%x0(1,j)/EARTH_RADIUS + cos(rotAlpha)*aPanels%x0(3,j)/EARTH_RADIUS))/&
				aPanels%h(j)
		else
			aPanels%relVort(j) = 0.0_kreal
			aPanels%div(j) = 0.0_kreal
			aPanels%h(j) = 0.0_kreal
			aPanels%potVort(j) = 0.0_kreal
		endif
	enddo
	
end subroutine

!
!----------------
! Module methods : module- or type-specific private functions
!----------------
!

subroutine InitLogger(aLog,rank)
	type(Logger), intent(inout) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

end module
