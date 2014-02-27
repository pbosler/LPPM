module TracerSetupModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	Defines tracer distributions used by SphereMesh.
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
public TracerSetup
public New, Delete, NullTracer
public InitCosineBellTracer, SetCosineBellTracerOnMesh
public COS_BELL_NINT, COS_BELL_NREAL
public SetFlowMapLatitudeTracerOnMesh


!
!----------------
! Types and module constants
!----------------
!
type TracerSetup
	integer(kint), pointer :: integers(:) => null()	!	integers required by tracer definitions
	real(kreal), pointer :: reals(:) => null() 		! 	real numbers required by tracer definitions
	integer(kint) :: tracerID						! 	index of starting tracer array in SphereMesh objects
end type

integer(kint), parameter :: COS_BELL_NINT = 0, COS_BELL_NREAL = 4

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'Tracer'
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
	type(TracerSetup), intent(out) :: self
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

subroutine DeletePrivate(self)
	type(TracerSetup), intent(inout) :: self
	if (associated(self%integers))	deallocate(self%integers)
	if (associated(self%reals))	deallocate(self%reals)
	self%tracerID = 0
	call Delete(log)
	logInit = .FALSE.
end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine InitCosineBellTracer(cosBell, lat0, lon0, rr, h0, tracerID)
	type(TracerSetup), intent(inout) :: cosBell
	real(kreal), intent(in) :: lat0, lon0, rr, h0
	integer(kint), intent(in) :: tracerID

	if ( size(cosBell%reals) /= 4) then
		call LogMessage(log, ERROR_LOGGING_LEVEL,'tracerSetup ERROR : ',' real array size incorrect')
		return
	endif

	cosBell%reals(1) = lat0	! latitude of center of cosine bell
	cosBell%reals(2) = lon0 ! longitude of center of cosine bell
	cosBell%reals(3) = rr	! radius of cosine bell
	cosBell%reals(4) = h0	! height of center of cosine bell

	cosBell%tracerID = tracerID
end subroutine

subroutine SetCosineBellTracerOnMesh(aMesh,cosBell)
	type(SphereMesh), intent(inout) :: aMesh
	type(TracerSetup), intent(in) :: cosBell
	real(kreal) :: xyzCent(3)
	integer(kint) :: j
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	xyzCent = [cos(cosBell%reals(1))*cos(cosBell%reals(2)), cos(cosBell%reals(1))*sin(cosBell%reals(2)), sin(cosBell%reals(1)) ]

	do j=1,aParticles%N
		aParticles%tracer(j,cosBell%tracerID) = CosineBellX(aParticles%x0(:,j),xyzCent,cosBell%reals(3),cosBell%reals(4))
	enddo
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			aPanels%tracer(j,cosBell%tracerID) = CosineBellX(aPanels%x0(:,j),xyzCent,cosBell%reals(3), cosBell%reals(4))
		else
			aPanels%tracer(j,cosBell%tracerID) = 0.0_kreal
		endif
	enddo

end subroutine

subroutine SetFlowMapLatitudeTracerOnMesh(aMesh,tracerID)
	type(SphereMesh), intent(inout) :: aMesh
	integer(kint), intent(in) :: tracerID
	! local variables
	integer(kint) :: j
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	do j=1,aParticles%N
		aParticles%tracer(j,tracerID) = Latitude(aParticles%x0(:,j))
	enddo
	do j=1,aPanels%N
		aPanels%traceR(j,tracerID) = Latitude(aPanels%x0(:,j))
	enddo
end subroutine

subroutine NullTracer(aMesh,nullScalar)
	type(SphereMesh), intent(inout) :: aMesh
	type(TracerSetup), intent(in) :: nullScalar
end subroutine


!
!----------------
! Module methods : module- or type-specific private functions
!----------------
!

function CosineBellX(xyz,xyzCent,rr,h0)
	real(kreal) :: CosineBellX
	real(kreal), intent(in) :: xyz(3), xyzCent(3), rr, h0
	real(kreal) :: r
	r = SphereDistance(xyz,xyzCent)
	if ( r < rr ) then
		CosineBellX = (h0/2.0_kreal)*(1.0_kreal + cos(PI*r/rr))
	else
		CosineBellX = 0.0_kreal
	endif
end function


subroutine InitLogger(aLog,rank)
	type(Logger), intent(inout) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

end module
