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
public InitSingleGaussianVortex, SetSingleGaussianVortexOnMesh
public SINGLE_VORTEX_NINT, SINGLE_VORTEX_NREAL

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
integer(kint), parameter :: SINGLE_VORTEX_NINT = 0, SINGLE_VORTEX_NREAL = 4

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

subroutine InitSingleGaussianVortex(gaussVort, lat0, lon0, beta, maxVal)
	type(BVESetup), intent(inout) :: gaussVort
	real(kreal), intent(in) :: lat0, lon0, beta, maxVal
	gaussVort%reals(1) = lat0
	gaussVort%reals(2) = lon0
	gaussVort%reals(3) = bb
	gaussVort%reals(4) = maxVal
end subroutine

subroutine SetSingleGaussianVortexOnMesh(aMesh,gaussVort)
	type(SphereMesh), intent(inout) :: aMesh
	type(BVESetup), intent(inout) :: gaussVort
	! local variables
	integer(kint) :: j
	real(kreal) :: xyzCent(3)
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	real(kreal) :: allocatable :: gVort(:)
	
	aParticles => aMesh%particles
	aPanels => aMesh%panels
	
	xyzCent = EARTH_RADIUS*[ cos(gaussVort%reals(1))*cos(gaussVort%reals(2)), &
							 cos(gaussVort%reals(1))*sin(gaussVort%reals(2)), &
							 sin(gaussVort%reals(1)) ]
	!
	!	set Gaussian constant if not set already
	!
	if (GAUSS_CONST == 0.0_kreal) then
		allocate(gVort(aPanels%N))
		gVort = 0.0_kreal
		do j=1,aPanels%N
			if ( .NOT. aPanels%hasChildren(j) ) then
				gVort(j) = GeneralGaussian(aPanels%x0(:,j),xyzCent, &
					gaussVort%reals(3), gaussVort%reals(4))
			endif
		enddo
		GAUSS_CONST = sum(gVort*aPanels%area(1:aPanels%N))/(4.0_kreal*PI*EARTH_RADIUS)
		deallocate(gVort)
	endif
	
	do j=1,aParticles%N
		aParticles%absVort(j) = GeneralGaussian(aParticles%x0(:,j), xyzCent, &
				gaussVort%reals(3), gaussVort%reals(4))  &
					- GAUSS_CONST + 2.0_kreal*OMEGA*aParticles%x0(3,j)/EARTH_RADIUS
		aParticles%relVort(j) = aParticles%absVort(j) - 2.0_kreal*OMEGA*aParticles%x(3,j)/EARTH_RADIUS			
	enddo
	
	do j=1,aPanels%N
		if ( aPanels%hasChildren(j) ) then
			aPanels%absVort(j) = 0.0_kreal
			aPanels%relVort(j) = 0.0_kreal
		else
			aPanels%absVort(j) = GeneralGaussian(aPanels%x0(:,j), xyzCent, &
				 gaussVort%reals(3), gaussVort%reals(4))  &
				 - GAUSS_CONST + 2.0_kreal*OMEGA*aPanels%x0(3,j)/EARTH_RADIUS
			aPanels%relVort(j) = aPanels%absVort(j) - 2.0_kreal*OMEGA*aPanels%x(3,j)/EARTH_RADIUS		
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

function GeneralGaussian(xyz,xyzCent, bb, maxVal)
	real(kreal) :: GeneralGaussian
	real(kreal), intent(in) :: xyz(3), xyzCent(3), bb, maxVal
	GeneralGaussian = maxVal * exp( -2.0_kreal*bb*bb*( EARTH_RADIUS*EARTH_RADIUS - sum(xyz*xyzCent)))
end function

subroutine InitLogger(aLog,rank)
	type(Logger), intent(inout) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

end module
