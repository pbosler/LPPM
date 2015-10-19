module PlaneTracerModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method
!------------------------------------------------------------------------------
!
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
!> @defgroup PlaneTracer Planar tracer distributions
!> defines passive tracer distributions for incompressible flow in 2d Cartesian geometry
!
!
! DESCRIPTION:
!> @file
!> defines passive tracer distributions for incompressible flow in 2d Cartesian geometry
!
!------------------------------------------------------------------------------
use NumberKindsModule
use LoggerModule
use PlaneGeomModule
use ParticlesModule
use EdgesModule
use PanelsModule
use PlaneMeshModule

implicit none

private
public TracerSetup
public New, Delete
public InitGaussianHill, SetGaussianHillOnMesh, GaussianHill
public InitSlottedCylinder, SetSlottedCylinderOnMesh, SlottedCylinder
public InitThreeTracers, SetThreeTracersOnMesh, ThreeTracers
public NullTracer
public InitDipoleIDTracer, SetDipoleIDTracerOnMesh, DipoleID

!
!----------------
! Types and module constants
!----------------
!
type TracerSetup
	integer(kint), pointer :: integers(:) => null() 	! integers required by test case definitions
	real(kreal), pointer :: reals(:) => null()			! reals required by test case definitions
end type

!integer(kint), parameter :: GAUSS_HILL_N_INT, GAUSS_HILL_N_REAL

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'TracerSetup'
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
	type(TracerSetup), intent(inout) :: self
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
subroutine NullTracer( aMesh, self)
	type(PlaneMesh), intent(inout) :: aMesh
	type(TracerSetup), intent(in) :: self
end subroutine

subroutine InitThreeTracers(self, tracerID)
	type(TracerSetup), intent(inout) :: self
	integer(kint), intent(in) :: tracerID

	call New(self, 1, 0)
	self%integers(1) = tracerID
end subroutine

subroutine SetThreeTracersOnMesh( aMesh, self)
	type(PlaneMesh), intent(inout) :: aMesh
	type(TracerSetup), intent(in) :: self
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	do j = 1, aParticles%N
		aParticles%tracer(j, self%integers(1)) = ThreeTracers(aParticles%x0(:,j))
	enddo
	do j = 1, aPanels%N
		if ( aPanels%hasChildren(j) ) then
			aPanels%tracer(j,self%integers(1)) = 0.0_kreal
		else
			aPanels%tracer(j,self%integers(1)) = ThreeTracers(aPanels%x0(:,j))
		endif
	enddo
end subroutine

function ThreeTracers(xy)
	real(kreal) :: ThreeTracers
	real(kreal), intent(in) :: xy(2)
	!
	real(kreal) :: r, cosArg

	! cosine hill
	r = sqrt( (xy(1) - 0.25_kreal)*(xy(1) - 0.25_kreal) + (xy(2) - 0.5_kreal)*(xy(2)-0.5_kreal))
	if ( r >= 0.15_kreal ) then
		cosArg = 0.15_kreal
	else
		cosArg = r
	endif
	cosArg = cosArg / 0.15_kreal
	ThreeTracers = 0.25_kreal * (1.0_kreal + cos(PI * cosArg))

	! slotted cylinder
	r = sqrt( (xy(1) - 0.5_kreal)*(xy(1)-0.5_kreal) + (xy(2) - 0.75_kreal)*(xy(2) - 0.75_kreal))
	if ( r < 0.15_kreal) ThreeTracers = ThreeTracers + 1.0_kreal
	if ( xy(1) < 0.525_kreal .AND. xy(1) > 0.475_kreal ) then
		if ( xy(2) > 0.6_kreal .AND. xy(2) < 0.85_kreal ) then
			ThreeTracers = ThreeTracers - 1.0_kreal
		endif
	endif

	! cone
	r = 1.0_kreal - sqrt( (xy(1) -0.5_kreal)*(xy(1) - 0.5_kreal) + (xy(2) - 0.25_kreal)*(xy(2) - 0.25_kreal))
	if ( r < 0.0_kreal ) r = 0.0_kreal
	ThreeTracers = ThreeTracers + r
end function

subroutine InitGaussianHill(self, tracerID, xCent, yCent, betaSq, amp)
	type(TracerSetup), intent(inout) :: self
	integer(kint), intent(in) :: tracerID
	real(kreal), intent(in) :: xCent, yCent, betaSq, amp

	call New(self, 1, 4)
	self%integers(1) = tracerID

	self%reals(1) = xCent
	self%reals(2) = yCent
	self%reals(3) = betaSq
	self%reals(4) = amp
end subroutine

subroutine SetGaussianHillOnMesh( aMesh, self )
	type(PlaneMesh), intent(inout) :: aMesh
	type(TracerSetup), intent(in) :: self
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j
	real(kreal) :: xCent, yCent, betaSq, amp

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	xCent = self%reals(1)
	yCent = self%reals(2)
	betaSq = self%reals(3)
	amp = self%reals(4)

	do j = 1, aParticles%N
		aParticles%tracer(j, self%integers(1)) = GaussianHill(aParticles%x0(:,j), xCent, yCent, betaSq, amp)
	enddo
	do j = 1, aPanels%N
		if ( aPanels%hasChildren(j) ) then
			aPanels%tracer(j, self%integers(1)) = 0.0_kreal
		else
			aPanels%tracer(j, self%integers(1)) = GaussianHill(aParticles%x0(:,j), xCent, yCent, betasq, amp)
		endif
	enddo
end subroutine

function GaussianHill( xy, xcent, ycent, betasq, amp)
	real(kreal) :: GaussianHill
	real(kreal), intent(in) :: xy(2), xcent, ycent, betasq, amp
	GaussianHill = amp * exp( - betasq*(xy(1)-xcent)*(xy(1)-xcent) - betasq*(xy(2)-ycent)*(xy(2)-ycent) )
end function

subroutine InitSlottedCylinder(self, tracerID, xcent, ycent, r, slotL, slotR, slotTop)
	type(TracerSetup), intent(inout) :: self
	integer(kint), intent(in) :: tracerID
	real(kreal), intent(in) :: xcent, ycent, r, slotL, slotR, slotTop
	!
	call New(self, 1, 6)
	self%integers(1) = tracerID

	self%reals(1) = xcent
	self%reals(2) = ycent
	self%reals(3) = r
	self%reals(4) = slotL
	self%reals(5) = slotR
	self%reals(6) = slotTop
end subroutine

subroutine SetSlottedCylinderOnMesh(aMesh, self)
	type(PlaneMesh), intent(inout) :: aMesh
	type(TracerSetup), intent(in) :: self
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	do j = 1, aParticles%N
		aParticles%tracer(j, self%integers(1)) = SlottedCylinder(aParticles%x0(:,j), self%reals(1), self%reals(2), &
				          self%reals(3), self%reals(4), self%reals(5), self%reals(6) )
	enddo
	do j = 1, aPanels%N
		if ( aPanels%hasChildren(j) ) then
			aPanels%tracer(j, self%integers(1) )  = 0.0_kreal
		else
			aPanels%tracer(j, self%integers(1) ) = SlottedCylinder( aPanels%x0(:,j), self%reals(1), self%reals(2), &
				          self%reals(3), self%reals(4), self%reals(5), self%reals(6) )
		endif
	enddo
end subroutine

function SlottedCylinder(xy, xcent, ycent, r, slotL, slotR, slotTop)
	real(kreal) :: SlottedCylinder
	real(kreal), intent(in) :: xy(2), xcent, ycent, r, slotL, slotR, slotTop
	!
	real(kreal) :: q, rr
	rr = sqrt( (xy(1) - xcent)*(xy(1) - xcent) + (xy(2) - ycent)*(xy(2) - ycent) )
	q = 0.0_kreal
	if ( rr < r ) q = 1.0_kreal
	if ( (q > 0.0_kreal) .AND. ( xy(2) > slotTop ) ) then
		if ( xy(1) < slotL .AND. xy(1) > slotR	) then
			q = 0.0_kreal
		endif
	endif
	SlottedCylinder = q
end function

function SignedDistance(xy, xc, yc, rad)
	real(kreal) :: SignedDistance
	real(kreal), intent(in) :: xy(2), xc, yc, rad
	!
	SignedDistance = rad - sqrt( (xy(1) - xc) * (xy(1)-xc) + (xy(2)-yc) * (xy(2) - yc) )
end function

subroutine InitDipoleIDTracer(self, xc1, yc1, rad1, xc2, yc2, rad2, tracerID)
	type(TracerSetup), intent(inout) :: self
	real(kreal), intent(in) :: xc1, yc1, rad1, xc2, yc2, rad2
	integer(kint), intent(in) :: tracerID

	call New(self, 1, 6 )

	self%reals(1) = xc1
	self%reals(2) = yc1
	self%reals(3) = rad1
	self%reals(4) = xc2
	self%reals(5) = yc2
	self%reals(6) = rad2
	self%integers(1) = tracerID
end subroutine

function DipoleID(xy, xc1, yc1, rad1, xc2, yc2, rad2)
	real(kreal) :: DipoleID
	real(kreal), intent(in) :: xy(2), xc1, yc1, rad1, xc2, yc2, rad2
	!
	real(kreal) :: d1, d2

	DipoleID = 0.0_kreal

	d1 = SignedDistance(xy, xc1, yc1, rad1)
	d2 = SignedDistance(xy, xc2, yc2, rad2)

	if ( d1 > ZERO_TOL ) DipoleID = DipoleID - 1.0_kreal
	if ( d2 > ZERO_TOL ) DipoleID = DipoleID + 1.0_kreal
end function

subroutine SetDipoleIDTracerOnMesh(aMesh, tsetup)
	type(PlaneMesh), intent(inout) :: aMesh
	type(TracerSetup), intent(in) :: tsetup
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	do j = 1, aParticles%N
		aParticles%tracer(j, tsetup%integers(1) ) = DipoleID( aParticles%x0(:,j), tsetup%reals(1), tsetup%reals(2), tsetup%reals(3), &
														tsetup%reals(4), tsetup%reals(5), tsetup%reals(6) )
	enddo
	do j = 1, aPanels%N
		if ( aPanels%hasChildren(j) ) then
			aPanels%tracer(j, tsetup%integers(1) ) = 0.0_kreal
		else
			aPanels%tracer(j, tsetup%integers(1) ) = DipoleID( aPanels%x0(:,j), tsetup%reals(1), tsetup%reals(2), tsetup%reals(3), &
														tsetup%reals(4), tsetup%reals(5), tsetup%reals(6) )
		endif
	enddo
end subroutine



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
