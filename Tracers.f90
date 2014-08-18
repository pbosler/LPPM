module TracerSetupModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method - Spherical Model
!------------------------------------------------------------------------------
!
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
!> @defgroup Tracers Spherical geometry tracer distributions
!> defines passive tracer distributions for flows in spherical geometry
!
!
! DESCRIPTION:
!> @file
!> defines passive tracer distributions for flows in spherical geometry
!
!------------------------------------------------------------------------------
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
public InitCosineBellTracer, SetCosineBellTracerOnMesh, COS_BELL_NINT, COS_BELL_NREAL
public SetFlowMapLatitudeTracerOnMesh
public InitGaussianHillsTracer, SetGaussianHillsTracerOnMesh, GAUSS_HILLS_N_INT, GAUSS_HILLS_N_REAL
public InitMovingVortsTracer, SetMovingVortsTracerOnMesh
public InitOneGaussianHillTracer, SetOneGaussianHillTracerOnMesh
public TracerVariance


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

integer(kint), parameter :: COS_BELL_NINT = 0, COS_BELL_NREAL = 4, &
							GAUSS_HILLS_N_INT = 1, GAUSS_HILLS_N_REAL = 2

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

	xyzCent = EARTH_RADIUS * [cos(cosBell%reals(1))*cos(cosBell%reals(2)), cos(cosBell%reals(1))*sin(cosBell%reals(2)), sin(cosBell%reals(1)) ]

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

subroutine InitMovingVortsTracer(mvTracer, initLon, initLat, tracerID)
	type(TracerSetup), intent(inout) :: mvTracer
	real(kreal), intent(in) :: initLon, initLat
	integer(kint), intent(in) :: tracerID
	mvTracer%reals(1) = initLon
	mvTracer%reals(2) = initLat
	mvTracer%integers(1) = tracerID
end subroutine

subroutine SetMovingVortsTracerOnMesh( aMesh, mvt)
	type(SphereMesh), intent(inout) :: aMesh
	type(TracerSetup), intent(in) :: mvt
	!
	integer(kint) :: j
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	real(kreal) :: vortStartingLon, vortStartingLat
	integer(kint) :: tracerID
	aParticles => aMesh%particles
	aPanels => aMesh%panels

!	vortStartingLon = mvt%reals(1)
!	vortStartingLat = mvt%reals(1)
	vortStartingLon = 1.5_kreal * PI
	vortStartingLat = 0.0_kreal
	tracerID = mvt%integers(1)

	do j = 1, aParticles%N
		aParticles%tracer(j, tracerID ) = MVTracerAtT0(aParticles%x0(:,j), vortStartingLon, vortStartingLat)
	enddo
	do j = 1, aPanels%N
		if ( aPanels%hasChildren(j) ) then
			aPanels%tracer(j, tracerID ) = 0.0_kreal
		else
			aPanels%tracer(j, tracerID ) = MVTracerAtT0(aPanels%x0(:,j), vortStartingLon, vortStartingLat)
		endif
	enddo
end subroutine

function MVTracerAtT0(xyz, vortStartingLon, vortStartingLat )
	real(kreal) :: MVTracerAtT0
	real(kreal), intent(in) :: xyz(3), vortStartingLon, vortStartingLat
	real(kreal) :: rho, lonPrime, latPrime, lon, lat

	lon = Longitude(xyz)
	lat = Latitude(xyz)

	lonPrime = atan4( cos(lat)*sin( lon - vortStartingLon),  cos(lat)*sin(vortStartingLat)*cos( lon - vortStartingLon) - cos(vortStartingLat)*sin(lat) )
	latPrime = asin( sin(lat)*sin(vortStartingLat) + cos(lat)*cos(vortStartingLat)*cos( lon - vortStartingLon ) )

	rho = 3.0_kreal * cos(latPrime)

	MVTracerAtT0 = 1.0_kreal - tanh( 0.2_kreal * rho * sin( lonPrime ) )
end function

subroutine InitOneGaussianHillTracer(gHill, hMax, beta, tracerID)
	type(TracerSetup), intent(inout) :: gHill
	real(kreal), intent(in) :: hMax, beta
	integer(kint), intent(in) :: tracerID
	!
	if ( size(gHill%reals) /= 2 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL,'tracerSetup ERROR : ',' real array size incorrect')
		return
	endif
	if ( size(gHill%integers) /=1 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL,'tracerSetup ERROR : ',' integer array size incorrect')
		return
	endif

	gHill%reals(1) = hmax
	gHill%reals(2) = beta
	gHill%integers(1) = tracerID
end subroutine

subroutine InitGaussianHillsTracer(gHills, hmax, beta, tracerID)
	type(TracerSetup), intent(inout) :: gHills
	real(kreal), intent(in) :: hmax, beta
	integer(kint), intent(in) :: tracerID
	!
	if ( size(gHills%reals) /= 2 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL,'tracerSetup ERROR : ',' real array size incorrect')
		return
	endif
	if ( size(gHills%integers) /=1 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL,'tracerSetup ERROR : ',' integer array size incorrect')
		return
	endif

	gHills%reals(1) = hmax
	gHills%reals(2) = beta
	gHills%integers(1) = tracerID
end subroutine

subroutine SetGaussianHillsTracerOnMesh(aMesh, gHills)
	type(SphereMesh), intent(inout) :: aMesh
	type(TracerSetup), intent(in) :: gHills
	!
	real(kreal) :: xyzCent1(3), xyzCent2(3)
	integer(kint) :: j
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	xyzCent1 = [ cos(5.0_kreal * PI / 6.0_kreal), sin( 5.0_kreal * PI / 6.0_kreal ), 0.0_kreal ]
	xyzCent2 = [ cos(7.0_kreal * PI / 6.0_kreal), sin( 7.0_kreal * PI / 6.0_kreal ), 0.0_kreal ]

	do j = 1, aParticles%N
		aParticles%tracer(j, gHills%integers(1) ) = GaussianHillsTracer(aParticles%x0(:,j)/EARTH_RADIUS, xyzcent1, xyzcent2, gHills%reals(1), gHills%reals(2))
	enddo
	do j = 1, aPanels%N
		if ( aPanels%hasChildren(j) ) then
			aPanels%tracer(j, gHills%integers(1) ) = 0.0_kreal
		else
			aPanels%tracer(j, gHills%integers(1) ) = GaussianHillsTracer(aPanels%x0(:,j)/EARTH_RADIUS, xyzcent1, xyzcent2, gHills%reals(1), gHills%reals(2))
		endif
	enddo
end subroutine

subroutine SetOneGaussianHillTracerOnMesh(aMesh, gHill)
	type(SphereMesh), intent(inout) :: aMesh
	type(TracerSetup), intent(in) :: gHill
	!
	real(kreal) :: xyzCent(3)
	integer(kint) :: j
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	xyzCent = [ 0.0_kreal, -1.0_kreal, 0.0_kreal ]

	do j = 1, aParticles%N
		aParticles%tracer(j, gHill%integers(1) ) = OneGaussianHillTracer(aParticles%x0(:,j)/EARTH_RADIUS, xyzcent, gHill%reals(1), gHill%reals(2))
	enddo
	do j = 1, aPanels%N
		if ( aPanels%hasChildren(j) ) then
			aPanels%tracer(j, gHill%integers(1) ) = 0.0_kreal
		else
			aPanels%tracer(j, gHill%integers(1) ) = OneGaussianHillTracer(aPanels%x0(:,j)/EARTH_RADIUS, xyzcent, gHill%reals(1), gHill%reals(2))
		endif
	enddo
end subroutine

function OneGaussianHillTracer(xyz, cent, hmax, beta)
		real(kreal) :: OneGaussianHillTracer
		real(kreal), intent(in) :: xyz(3), cent(3), hmax, beta
		OneGaussianHillTracer = exp( -beta * ( sum( (xyz-cent)*(xyz-cent) ) ) )
end function

function GaussianHillsTracer(xyz, cent1, cent2, hmax, beta)
	real(kreal) :: GaussianHillsTracer
	real(kreal), intent(in) :: xyz(3), cent1(3), cent2(3), hmax, beta
	!
	real(kreal) :: h1, h2

	h1 = hmax * exp( -beta * ( sum( (xyz-cent1) * (xyz-cent1) ) ) )
	h2 = hmax * exp( -beta * ( sum( (xyz-cent2) * (xyz-cent2) ) ) )

	GaussianHillsTracer = h1 + h2
end function

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

function TracerVariance(aMesh, tracerID)
	real(kreal) :: TracerVariance
	type(SphereMesh), intent(in) :: aMesh
	integer(kint), intent(in) :: tracerID
	!
	type(Panels), pointer :: aPanels
	aPanels => aMesh%panels
	tracerVariance = sum( aPanels%tracer(1:aPanels%N,tracerID) * aPanels%tracer(1:aPanels%N,tracerID) * aPanels%area(1:aPanels%N) )
end function

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
