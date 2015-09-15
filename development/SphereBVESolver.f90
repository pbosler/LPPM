module SphereBVESolverModule

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use ParticlesModule
use EdgesModule
use FacesModule
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use SphereGeomModule
use SphereBVEModule

implicit none

include 'mpif.h'

private
!
!----------------
! Module types and public declarations
!----------------
!
public BVESolver, New, Delete
public Timestep

type BVESolver
	real(kreal), pointer :: xIn(:) => null()
	real(kreal), pointer :: xStage1(:) => null()
	real(kreal), pointer :: xStage2(:) => null()
	real(kreal), pointer :: xStage3(:) => null()
	real(kreal), pointer :: xStage4(:) => null()
	real(kreal), pointer :: yIn(:) => null()
	real(kreal), pointer :: yStage1(:) => null()
	real(kreal), pointer :: yStage2(:) => null()
	real(kreal), pointer :: yStage3(:) => null()
	real(kreal), pointer :: yStage4(:) => null()
	real(kreal), pointer :: zIn(:) => null()
	real(kreal), pointer :: zStage1(:) => null()
	real(kreal), pointer :: zStage2(:) => null()
	real(kreal), pointer :: zStage3(:) => null()
	real(kreal), pointer :: zStage4(:) => null()
	real(kreal), pointer :: relVortIn(:) => null()
	real(kreal), pointer :: relVortStage1(:) => null()
	real(kreal), pointer :: relVortStage2(:) => null()
	real(kreal), pointer :: relVortStage3(:) => null()
	real(kreal), pointer :: relVortStage4(:) => null()
	
	contains
		final :: deletePrivate
end type	

!
!----------------
! Module interfaces
!----------------
!
interface New
	module procedure newPrivate
end interface

interface Delete
	module procedure deletePrivate
end interface

interface Timestep
	module procedure timestepPrivate
end interface

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'BVESolver'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains
!
!----------------
! public methods
!----------------
!
subroutine newPrivate(self, sphereBVE )
	type(BVESolver), intent(out) :: self
	type(BVEMesh), intent(in) :: sphereBVE
	integer(kint) :: nP
	
	if (.NOT. logInit ) call InitLogger(log, procRank)
	
	nP = sphereBVE%mesh%particles%N
	
	allocate(self%xIn(nP))
	allocate(self%xStage1(nP))
	allocate(self%xStage2(nP))
	allocate(self%xStage3(nP))
	allocate(self%xStage4(nP))
	allocate(self%yIn(nP))
	allocate(self%yStage1(nP))
	allocate(self%yStage2(nP))
	allocate(self%yStage3(nP))
	allocate(self%yStage4(nP))
	allocate(self%zIn(nP))
	allocate(self%zStage1(nP))
	allocate(self%zStage2(nP))
	allocate(self%zStage3(nP))
	allocate(self%zStage4(nP))
	allocate(self%relVortIn(nP))
	allocate(self%relVortStage1(nP))
	allocate(self%relVortStage2(nP))
	allocate(self%relVortStage3(nP))
	allocate(self%relVortStage4(nP))
end subroutine

subroutine deletePrivate(self)
	type(BVESolver), intent(inout) :: self
	if ( associated(self%xIn) ) then
		deallocate(self%xIn)
		deallocate(self%xStage1)
		deallocate(self%xStage2)
		deallocate(self%xStage3)
		deallocate(self%xStage4)
		deallocate(self%yIn)
		deallocate(self%yStage1)
		deallocate(self%yStage2)
		deallocate(self%yStage3)
		deallocate(self%yStage4)
		deallocate(self%zIn)
		deallocate(self%zStage1)
		deallocate(self%zStage2)
		deallocate(self%zStage3)
		deallocate(self%zStage4)
		deallocate(self%relVortIn)
		deallocate(self%relVortStage1)
		deallocate(self%relVortStage2)
		deallocate(self%relVortStage3)
		deallocate(self%relVortStage4)
	endif
end subroutine

subroutine timestepPrivate( self, sphereBVE, dt )
	type(BVESolver), intent(inout) :: self
	type(BVEMesh), intent(inout) :: sphereBVE
	real(kreal), intent(in) :: dt
	!
	integer(kint) :: i, nP

	nP = sphereBVE%mesh%particles%N
	!
	!	RK Stage 1
	!
	self%xStage1 = dt * sphereBVE%velocity%xComp(1:nP)
	self%yStage1 = dt * sphereBVE%velocity%yComp(1:nP)
	self%zStage1 = dt * sphereBVE%velocity%zComp(1:nP)
	self%relVortStage1 =  - dt * 2.0_kreal * sphereBVE%rotationRate * sphereBVE%velocity%zComp(1:nP) / sphereBVE%radius
	!
	! RK Stage 2
	!
	self%xIn = sphereBVE%mesh%particles%x(1:nP) + 0.5_kreal * self%xStage1
	self%yIn = sphereBVE%mesh%particles%y(1:nP) + 0.5_kreal * self%yStage1
	self%zIn = sphereBVE%mesh%particles%z(1:nP) + 0.5_kreal * self%zStage1
	self%relVortIn = sphereBVE%relVort%scalar(1:nP) + 0.5_kreal * self%relVortStage1

!	BVESphereVelocity( u, v, w, x, y, z, relVortIn, areaIn, sphereRadius, rotationRate, activeMask, mpiParticles )

	call BVESphereVelocity(self%xStage2, self%yStage2, self%zStage2, self%xIn, self%yIn, self%zIn, &
			self%relVortIn, sphereBVE%mesh%particles%area(1:nP), sphereBVE%radius, sphereBVE%rotationRate, &
			sphereBVE%mesh%particles%isActive(1:nP), sphereBVE%mpiParticles)
	
	self%relVortStage2 = - dt * 2.0_kreal * sphereBVE%rotationRate * self%zStage2 / sphereBVE%radius
	self%xStage2 = dt * self%xStage2
	self%yStage2 = dt * self%yStage2
	self%zStage2 = dt * self%zStage2
	
	!
	! RK Stage 3
	!
	self%xIn = sphereBVE%mesh%particles%x(1:nP) + 0.5_kreal * self%xStage2
	self%yIn = sphereBVE%mesh%particles%y(1:nP) + 0.5_kreal * self%yStage2
	self%zIn = sphereBVE%mesh%particles%z(1:nP) + 0.5_kreal * self%zStage2
	self%relVortIn = sphereBVE%relVort%scalar(1:nP) + 0.5_kreal * self%relVortStage2
	
	call BVESphereVelocity( self%xStage3, self%yStage3, self%zStage3, self%xIn, self%yIn, self%zIn, &
			self%relVortIn, sphereBVE%mesh%particles%area(1:nP), sphereBVE%radius, sphereBVE%rotationRate, &
			sphereBVE%mesh%particles%isActive(1:nP), sphereBVE%mpiParticles)
	
	self%relVortStage3 = - dt * 2.0_kreal * sphereBVE%rotationRate * self%zStage3 / sphereBVE%radius
	self%xStage3 = dt * self%xStage3
	self%yStage3 = dt * self%yStage3
	self%zStage3 = dt * self%zStage3
	
	!
	! RK Stage 4
	!
	self%xIn = sphereBVE%mesh%particles%x(1:nP) + self%xStage3
	self%yIn = sphereBVE%mesh%particles%y(1:nP) + self%yStage3
	self%zIn = sphereBVE%mesh%particles%z(1:nP) + self%zStage3
	self%relVortIn = sphereBVE%relVort%scalar(1:nP) + self%relVortStage3
	
	call BVESphereVelocity( self%xStage4, self%yStage4, self%zStage4, self%xIn, self%yIn, self%zIn, &
			self%relVortIn, sphereBVE%mesh%particles%area(1:nP), sphereBVE%radius, sphereBVE%rotationRate, &
			sphereBVE%mesh%particles%isActive(1:nP), sphereBVE%mpiParticles)
	
	self%relVortStage4 = - dt * 2.0_kreal * sphereBVE%rotationRate * self%zStage4 / sphereBVE%radius
	self%xStage4 = dt * self%xStage4
	self%yStage4 = dt * self%yStage4
	self%zStage4 = dt * self%zStage4
	
	!
	!	update mesh
	!	
	sphereBVE%mesh%particles%x(1:nP) = sphereBVE%mesh%particles%x(1:nP) + self%xStage1/6.0_kreal + &
		self%xStage2/3.0_kreal + self%xStage3/3.0_kreal + self%xStage4/6.0_kreal
	sphereBVE%mesh%particles%y(1:nP) = sphereBVE%mesh%particles%y(1:nP) + self%yStage1/6.0_kreal + &
		self%yStage2/3.0_kreal + self%yStage3/3.0_kreal + self%yStage4/6.0_kreal
	sphereBVE%mesh%particles%z(1:nP) = sphereBVE%mesh%particles%z(1:nP) + self%zStage1/6.0_kreal + &
		self%zStage2/3.0_kreal + self%zStage3/3.0_kreal + self%zStage4/6.0_kreal
	sphereBVE%relVort%scalar(1:nP) = sphereBVE%relVort%scalar(1:nP) + self%relVortStage1/6.0_kreal + &
		self%relVortStage2/3.0_kreal + self%relVortStage3/3.0_kreal + self%relVortStage4/6.0_kreal
	
	call SetVelocityOnMesh(sphereBVE)	
	call SetStreamFunctionsOnMesh(sphereBVE)
end subroutine


!
!----------------
! private methods
!----------------
!
subroutine BVESphereVelocity( u, v, w, x, y, z, relVortIn, areaIn, sphereRadius, rotationRate, activeMask, mpiParticles )
	real(kreal), dimension(:), intent(out) :: u
	real(kreal), dimension(:), intent(out) :: v
	real(kreal), dimension(:), intent(out) :: w
	real(kreal), dimension(:), intent(in) :: x
	real(kreal), dimension(:), intent(in) :: y
	real(kreal), dimension(:), intent(in) :: z
	real(kreal), dimension(:), intent(in) :: relVortIn
	real(kreal), dimension(:), intent(in) :: areaIn
	real(kreal), intent(in) :: sphereRadius
	real(kreal), intent(in) :: rotationRate
	logical(klog), dimension(:), intent(in) :: activeMask
	type(MPISetup), intent(in) :: mpiParticles
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: strength
	
	do i = mpiParticles%indexStart(procRank), mpiParticles%indexEnd(procRank)
		u(i) = 0.0_kreal
		v(i) = 0.0_kreal
		w(i) = 0.0_kreal
		do j = 1, i - 1
			if ( activeMask(j) ) then
				strength = -relVortIn(j)*areaIn(j) / ( 4.0_kreal * PI * sphereRadius * &
					( sphereRadius * sphereRadius - x(i)*x(j) - y(i)*y(j) - z(i)*z(j) ))
				u(i) = u(i) + ( y(i)*z(j) - z(i)*y(j) ) * strength
				v(i) = v(i) + ( z(i)*x(j) - x(i)*z(j) ) * strength
				w(i) = w(i) + ( x(i)*y(j) - y(i)*x(j) ) * strength
			endif
		enddo
		do j = i + 1, size(x)
			if ( activeMask(j) ) then
				strength = -relVortIn(j)*areaIn(j) / ( 4.0_kreal * PI * sphereRadius * &
					( sphereRadius * sphereRadius - x(i)*x(j) - y(i)*y(j) - z(i)*z(j) ))
				u(i) = u(i) + ( y(i)*z(j) - z(i)*y(j) ) * strength
				v(i) = v(i) + ( z(i)*x(j) - x(i)*z(j) ) * strength
				w(i) = w(i) + ( x(i)*y(j) - y(i)*x(j) ) * strength
			endif
		enddo
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST( u(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( v(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( w(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	enddo
end subroutine

subroutine InitLogger(aLog,rank)
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	if ( rank == 0 ) then
		call New(aLog,logLevel)
	else
		call New(aLog,ERROR_LOGGING_LEVEL)
	endif
	logInit = .TRUE.
end subroutine

end module