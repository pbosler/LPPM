module SphereSWESolverModule

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use ParticlesModule
use EdgesModule
use FacesModule
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use PSEDirectSumModule
use SphereGeomModule
use SphereSWEModule

implicit none

include 'mpif.h'

private
!
!----------------
! Module types and public declarations
!----------------
!
public SWESolver, New, Delete
public Timestep

type SWESolver
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
	real(kreal), pointer :: u(:) => null()
	real(kreal), pointer :: v(:) => null()
	real(kreal), pointer :: w(:) => null()
	real(kreal), pointer :: relVortIn(:) => null()
	real(kreal), pointer :: relVortStage1(:) => null()
	real(kreal), pointer :: relVortStage2(:) => null()
	real(kreal), pointer :: relVortStage3(:) => null()
	real(kreal), pointer :: relVortStage4(:) => null()
	real(kreal), pointer :: divIn(:) => null()
	real(kreal), pointer :: divStage1(:) => null()
	real(kreal), pointer :: divStage2(:) => null()
	real(kreal), pointer :: divStage3(:) => null()
	real(kreal), pointer :: divStage4(:) => null()
	real(kreal), pointer :: hIn(:) => null()
	real(kreal), pointer :: hStage1(:) => null()
	real(kreal), pointer :: hStage2(:) => null()
	real(kreal), pointer :: hStage3(:) => null()
	real(kreal), pointer :: hStage4(:) => null()
	real(kreal), pointer :: areaIn(:) => null()
	real(kreal), pointer :: areaStage1(:) => null()
	real(kreal), pointer :: areaStage2(:) => null()
	real(kreal), pointer :: areaStage3(:) => null()
	real(kreal), pointer :: areaStage4(:) => null()
	
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

interface 
	function LapDepthFn( x, y, z )
		real(8) :: LapDepthFn
		real(8), intent(in) :: x, y, z
	end function
end interface 

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'SWESphereK4'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains
!
!----------------
! public methods
!----------------
!
subroutine newPrivate( self, sphereSWE )
	type(SWESolver), intent(out) :: self
	type(SWEMesh), intent(in) :: sphereSWE
	!
	integer(kint) :: nP
	
	if ( .NOT. logInit ) call InitLogger(log, procRank)
	
	nP = sphereSWE%mesh%particles%N
	
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
	allocate(self%u(nP))
	allocate(self%v(nP))
	allocate(self%w(nP))
	allocate(self%relVortIn(nP))
	allocate(self%relVortStage1(nP))
	allocate(self%relVortStage2(nP))
	allocate(self%relVortStage3(nP))
	allocate(self%relVortStage4(nP))
	allocate(self%divIn(nP))
	allocate(self%divStage1(nP))
	allocate(self%divStage2(nP))
	allocate(self%divStage3(nP))
	allocate(self%divStage4(nP))
	allocate(self%hIn(nP))
	allocate(self%hStage1(nP))
	allocate(self%hStage2(nP))
	allocate(self%hStage3(nP))
	allocate(self%hStage4(nP))
	allocate(self%areaIn(nP))
	allocate(self%areaStage1(nP))
	allocate(self%areaStage2(nP))
	allocate(self%areaStage3(nP))
	allocate(self%areaStage4(nP))
end subroutine

subroutine deletePrivate(self)
	type(SWESolver), intent(inout) :: self
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
		deallocate(self%u)
		deallocate(self%v)
		deallocate(self%w)
		deallocate(self%relVortIn)
		deallocate(self%relVortStage1)
		deallocate(self%relVortStage2)
		deallocate(self%relVortStage3)
		deallocate(self%relVortStage4)
		deallocate(self%divIn)
		deallocate(self%divStage1)
		deallocate(self%divStage2)
		deallocate(self%divStage3)
		deallocate(self%divStage4)
		deallocate(self%hIn)
		deallocate(self%hStage1)
		deallocate(self%hStage2)
		deallocate(self%hStage3)
		deallocate(self%hStage4)
		deallocate(self%areaIn)
		deallocate(self%areaStage1)
		deallocate(self%areaStage2)
		deallocate(self%areaStage3)
		deallocate(self%areaStage4)
	endif
end subroutine

subroutine timestepPrivate( self, sphereSWE, dt, pseSetup )
	type(SWESolver), intent(inout) :: self
	type(SWEMesh), intent(inout) :: sphereSWE
	real(kreal), intent(in) :: dt
	type(PSE), intent(in) :: pseSetup
	!
	integer(kint) :: i, nP
	real(kreal) :: omega
	
	nP = sphereSWE%mesh%particles%N

	!
	!	RK Stage1
	!
	call SWESphereRHS( self%relVortStage1, self%divStage1, self%hStage1, self%areaStage1, &
					   sphereSWE%mesh%particles%x(1:nP), sphereSWE%mesh%particles%y(1:nP), sphereSWE%mesh%particles%z(1:nP), &
					   sphereSWE%velocity%xComp(1:nP), sphereSWE%velocity%yComp(1:nP), sphereSWE%velocity%zComp(1:nP), &
					   sphereSWE%relVort%scalar(1:nP), sphereSWE%potVort%scalar(1:nP), &
					   sphereSWE%divergence%scalar(1:nP), sphereSWE%h%scalar(1:nP), &
					   sphereSWE%mesh%particles%area(1:nP), sphereSWE%rotationRate, sphereSWE%g, sphereSWE%radius, &
					   sphereSWE%mesh%particles%isActive(1:nP), pseSetup, sphereSWE%mpiParticles)
	
	self%xStage1 = dt * sphereSWE%velocity%xComp(1:nP)
	self%yStage1 = dt * sphereSWE%velocity%yComp(1:nP)
	self%zStage1 = dt * sphereSWE%velocity%zComp(1:nP)
	self%relVortStage1 = dt * self%relVortStage1
	self%divStage1 = dt * self%divStage1
	self%hStage1 = dt * self%hStage1
	self%areaStage1 = dt * self%areaStage1
	
	!
	!	RK Stage 2
	!
	self%xIn = sphereSWE%mesh%particles%x(1:nP) + 0.5_kreal * self%xStage1
	self%yIn = sphereSWE%mesh%particles%y(1:nP) + 0.5_kreal * self%yStage1
	self%zIn = sphereSWE%mesh%particles%z(1:nP) + 0.5_kreal * self%zStage1
	self%relVortIn = sphereSWE%relVort%scalar(1:nP) + 0.5_kreal * self%relVortStage1
	self%divIn = sphereSWE%divergence%scalar(1:nP) + 0.5_kreal * self%divStage1
	self%hIn = sphereSWE%h%scalar(1:nP) + 0.5_kreal * self%hStage1
	self%areaIn = sphereSWE%mesh%particles%area(1:nP) + 0.5_kreal * self%areaStage1
	
	call SWESphereVelocity( self%u, self%v, self%w, self%xIn, self%yIn, self%zIn, self%relVortIn, self%divIn, self%areaIn, &
			sphereSWE%radius, sphereSWE%mesh%particles%isActive, sphereSWE%mpiParticles )
		
	call SWESphereRHS( self%relVortStage2, self%divStage2, self%hStage2, self%areaStage2, &
					   self%xIn, self%yIn, self%zIn, self%u, self%v, self%w, self%relVortIn, sphereSWE%potVort%scalar(1:nP), &
					   self%divIn, self%hIn, self%areaIn, sphereSWE%rotationRate, sphereSWE%g, sphereSWE%radius, &
					   sphereSWE%mesh%particles%isActive(1:nP), pseSetup, sphereSWE%mpiParticles)
	self%xStage2 = dt * self%u
	self%yStage2 = dt * self%v
	self%zStage2 = dt * self%w
	self%relVortStage2 = dt * self%relVortStage2
	self%divStage2 = dt * self%divStage2
	self%hStage2 = dt * self%hStage2
	self%areaStage2 = dt * self%areaStage2
	
	!
	!	RK Stage 3
	!
	self%xIn = sphereSWE%mesh%particles%x(1:nP) + 0.5_kreal * self%xStage2
	self%yIn = sphereSWE%mesh%particles%y(1:nP) + 0.5_kreal * self%yStage2
	self%zIn = sphereSWE%mesh%particles%z(1:nP) + 0.5_kreal * self%zStage2
	self%relVortIn = sphereSWE%relVort%scalar(1:nP) + 0.5_kreal * self%relVortStage2
	self%divIn = sphereSWE%divergence%scalar(1:nP) + 0.5_kreal * self%divStage2
	self%hIn = sphereSWE%h%scalar(1:nP) + 0.5_kreal * self%hStage2
	self%areaIn = sphereSWE%mesh%particles%area(1:nP) + 0.5_kreal * self%areaStage2
	
	call SWESphereVelocity( self%u, self%v, self%w, self%xIn, self%yIn, self%zIn, self%relVortIn, self%divIn, self%areaIn, &
			sphereSWE%radius, sphereSWE%mesh%particles%isActive, sphereSWE%mpiParticles )
		
	call SWESphereRHS( self%relVortStage3, self%divStage3, self%hStage3, self%areaStage3, &
					   self%xIn, self%yIn, self%zIn, self%u, self%v, self%w, self%relVortIn, sphereSWE%potVort%scalar(1:nP), &
					   self%divIn, self%hIn, self%areaIn, sphereSWE%rotationRate, sphereSWE%g, sphereSWE%radius, &
					   sphereSWE%mesh%particles%isActive(1:nP), pseSetup, sphereSWE%mpiParticles)
	self%xStage3 = dt * self%u
	self%yStage3 = dt * self%v
	self%zStage3 = dt * self%w
	self%relVortStage3 = dt * self%relVortStage3
	self%divStage3 = dt * self%divStage3
	self%hStage3 = dt * self%hStage3
	self%areaStage3 = dt * self%areaStage3
	
	!
	!	RK Stage 4
	!
	self%xIn = sphereSWE%mesh%particles%x(1:nP) + self%xStage3
	self%yIn = sphereSWE%mesh%particles%y(1:nP) + self%yStage3
	self%zIn = sphereSWE%mesh%particles%z(1:nP) + self%zStage3
	self%relVortIn = sphereSWE%relVort%scalar(1:nP) + self%relVortStage3
	self%divIn = sphereSWE%divergence%scalar(1:nP) + self%divStage3
	self%hIn = sphereSWE%h%scalar(1:nP) + self%hStage3
	self%areaIn = sphereSWE%mesh%particles%area(1:nP) + self%areaStage3
	
	call SWESphereVelocity( self%u, self%v, self%w, self%xIn, self%yIn, self%zIn, self%relVortIn, self%divIn, self%areaIn, &
			sphereSWE%radius, sphereSWE%mesh%particles%isActive, sphereSWE%mpiParticles )
		
	call SWESphereRHS( self%relVortstage4, self%divstage4, self%hstage4, self%areastage4, &
					   self%xIn, self%yIn, self%zIn, self%u, self%v, self%w, self%relVortIn, sphereSWE%potVort%scalar(1:nP), &
					   self%divIn, self%hIn, self%areaIn, sphereSWE%rotationRate, sphereSWE%g, sphereSWE%radius, &
					   sphereSWE%mesh%particles%isActive(1:nP), pseSetup, sphereSWE%mpiParticles)

	self%xstage4 = dt * self%u
	self%ystage4 = dt * self%v
	self%zstage4 = dt * self%w
	self%relVortstage4 = dt * self%relVortstage4
	self%divstage4 = dt * self%divstage4
	self%hstage4 = dt * self%hstage4
	self%areastage4 = dt * self%areastage4
	
	!
	!	update mesh state
	!		   
	sphereSWE%mesh%particles%x(1:nP) = sphereSWE%mesh%particles%x(1:nP) + self%xStage1 / 6.0_kreal + self%xStage2 / 3.0_kreal + &
				self%xStage3 / 3.0_kreal + self%xStage4 / 6.0_kreal
	sphereSWE%mesh%particles%y(1:nP) = sphereSWE%mesh%particles%y(1:nP) + self%yStage1 / 6.0_kreal + self%yStage2 / 3.0_kreal + &
				self%yStage3 / 3.0_kreal + self%yStage4 / 6.0_kreal
	sphereSWE%mesh%particles%z(1:nP) = sphereSWE%mesh%particles%z(1:nP) + self%zStage1 / 6.0_kreal + self%zStage2 / 3.0_kreal + &
				self%zStage3 / 3.0_kreal + self%zStage4 / 6.0_kreal			
	sphereSWE%relVort%scalar(1:nP) = sphereSWE%relVort%scalar(1:nP) + self%relVortStage1 / 6.0_kreal + &
			self%relVortStage2 / 3.0_kreal + self%relVortStage3 / 3.0_kreal + self%relVortStage4 / 6.0_kreal 
	sphereSWE%divergence%scalar(1:nP) = sphereSWE%divergence%scalar(1:nP) + self%divStage1 / 6.0_kreal + &
			self%divStage2 / 3.0_kreal + self%divStage3 / 3.0_kreal + self%divStage4 / 6.0_kreal
	sphereSWE%h%scalar(1:nP) = sphereSWE%h%scalar(1:nP) + self%hStage1 / 6.0_kreal + self%hStage2 / 3.0_kreal + &
			self%hStage3 / 3.0_kreal + self%hStage4 / 6.0_kreal
	sphereSWE%mesh%particles%area(1:nP) = sphereSWE%mesh%particles%area(1:nP) + self%areaStage1 / 6.0_kreal + &
			self%areaStage2 / 3.0_kreal + self%areaStage3 / 3.0_kreal + self%areaStage4 / 6.0_kreal
	
	call SWESphereVelocity( sphereSWE%velocity%xComp(1:nP), sphereSWE%velocity%yComp(1:nP), sphereSWE%velocity%zComp(1:nP), &
				sphereSWE%mesh%particles%x(1:nP), sphereSWE%mesh%particles%y(1:nP), sphereSWE%mesh%particles%z(1:nP), &
				sphereSWE%relVort%scalar(1:nP), sphereSWE%divergence%scalar(1:nP), sphereSWE%mesh%particles%area(1:nP), &
				sphereSWE%radius, sphereSWE%mesh%particles%isActive(1:nP), sphereSWE%mpiParticles)
end subroutine


!
!----------------
! private methods
!----------------
!
subroutine SWESphereVelocity( u, v, w, xIn, yIn, zIn, relVortIn, divIn, areaIn, sphereRadius, activeMask, mpiParticles )
	real(kreal), dimension(:), intent(out) :: u
	real(kreal), dimension(:), intent(out) :: v
	real(kreal), dimension(:), intent(out) :: w
	real(kreal), dimension(:), intent(in) :: xIn
	real(kreal), dimension(:), intent(in) :: yIn
	real(kreal), dimension(:), intent(in) :: zIn
	real(kreal), dimension(:), intent(in) :: relVortIn
	real(kreal), dimension(:), intent(in) :: divIn
	real(kreal), dimension(:), intent(in) :: areaIn
	real(kreal), intent(in) :: sphereRadius
	logical(klog), dimension(:), intent(in) :: activeMask
	type(MPISetup), intent(in) :: mpiParticles
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: denom, dotProd, r2
	real(kreal), dimension(3) :: xi, xj
	real(kreal), dimension(3,3) :: proj
	
	r2 = sphereRadius * sphereRadius
	
	do i = mpiParticles%indexStart(procRank), mpiParticles%indexEnd(procRank)
		u(i) = 0.0_kreal
		v(i) = 0.0_kreal
		w(i) = 0.0_kreal
		xi = [xIn(i), yIn(i), zIn(i)]
		do j = 1, i - 1
			if ( activeMask(j) ) then
				xj = [xIn(j), yIn(j), zIn(j) ]
				dotProd = sum( xi * xj )
				denom = r2 - dotProd
				u(i) = u(i) + (-( yIn(i) * zIn(j) - yIn(j) * zIn(i) ) * relVortIn(j) + &
					  (xIn(i) * dotProd - r2 * xIn(j)) * divIn(j) / r2 ) * areaIn(j) / denom
				v(i) = v(i) + (-( zIn(i) * xIn(j) - zIn(j) * xIn(i) ) * relVortIn(j) + &
					  (yIn(i) * dotProd - r2 * yIn(j)) * divIn(j) / r2 ) * areaIn(j) / denom
				w(i) = w(i) + (-( xIn(i) * yIn(j) - xIn(j) * yIn(i) ) * relVortIn(j) + &
					  (zIn(i) * dotProd - r2 * zIn(j)) * divIn(j) / r2 )* areaIn(j) / denom
			endif
		enddo
		do j = i + 1, size(xIn)
			if ( activeMask(j) ) then
				xj = [xIn(j), yIn(j), zIn(j) ]
				dotProd = sum( xi * xj )
				denom = r2 - dotProd
				u(i) = u(i) + (-( yIn(i) * zIn(j) - yIn(j) * zIn(i) ) * relVortIn(j) + &
					  (xIn(i) * dotProd - r2 * xIn(j)) * divIn(j) / r2 ) * areaIn(j) / denom
				v(i) = v(i) + (-( zIn(i) * xIn(j) - zIn(j) * xIn(i) ) * relVortIn(j) + &
					  (yIn(i) * dotProd - r2 * yIn(j)) * divIn(j) / r2 ) * areaIn(j) / denom
				w(i) = w(i) + (-( xIn(i) * yIn(j) - xIn(j) * yIn(i) ) * relVortIn(j) + &
					  (zIn(i) * dotProd - r2 * zIn(j)) * divIn(j) / r2 )* areaIn(j) / denom
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

subroutine SWESphereRHS( dRelVort, dDiv, dH, dArea, &
						 xIn, yIn, zIn, uIn, vIn, wIn, relVortIn, potVortIn, divIn, hIn, areaIn, &
						 rotRate, g, sphereRadius, &
						 activeMask, pseSetup, mpiParticles)
	real(kreal), dimension(:), intent(out) :: dRelVort
	real(kreal), dimension(:), intent(out) :: dDiv
	real(kreal), dimension(:), intent(out) :: dH
	real(kreal), dimension(:), intent(out) :: dArea
	real(kreal), dimension(:), intent(in) :: xIn
	real(kreal), dimension(:), intent(in) :: yIn
	real(kreal), dimension(:), intent(in) :: zIn
	real(kreal), dimension(:), intent(in) :: uIn
	real(kreal), dimension(:), intent(in) :: vIn
	real(kreal), dimension(:), intent(in) :: wIn
	real(kreal), dimension(:), intent(in) :: relVortIn
	real(kreal), dimension(:), intent(in) :: potVortIn
	real(kreal), dimension(:), intent(in) :: divIn
	real(kreal), dimension(:), intent(in) :: hIn
	real(kreal), dimension(:), intent(in) :: areaIn
	real(kreal), intent(in) :: rotRate
	real(kreal), intent(in) :: g
	real(kreal), intent(in) :: sphereRadius
	logical(klog), dimension(:), intent(in) :: activeMask
	type(PSE), intent(in) :: pseSetup
	type(MPISetup), intent(in) :: mpiParticles
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: kIn, kernelD1, kernelLap
	real(kreal) :: ux, uy, uz, vx, vy, vz, wx, wy, wz
	real(kreal) :: lapH, doubleDot
	real(kreal), dimension(3) :: xi, xj, gradU, gradV, gradW
	real(kreal), dimension(3,3) :: projMat
	
	do i = 1, size(xIn)
		dH(i) = - hIn(i) * divIn(i)
		dRelVort(i) = - potVortIn(i) * hIn(i) * divIn(i) - 2.0_kreal * rotRate / sphereRadius * wIn(i)
	enddo
	do i = 1, size(xIn)
		if ( activeMask(i) ) then
			dArea(i) = divIn(i) * areaIn(i)
		endif
	enddo
	
	do i = mpiParticles%indexStart(procRank), mpiParticles%indexEnd(procRank)
		ux = 0.0_kreal
		uy = 0.0_kreal
		uz = 0.0_kreal
		vx = 0.0_kreal
		vy = 0.0_kreal
		vz = 0.0_kreal
		wx = 0.0_kreal
		wy = 0.0_kreal
		wz = 0.0_kreal
		doubleDot = 0.0_kreal
		lapH = 0.0_kreal
		!
		! PSE spatial derivatives
		!
		xi = [ xIn(i), yIn(i), zIn(i) ]
		projMat = SphereProjection( xi )
		do j = 1, size(xIN)
			if ( activeMask(j) ) then
				xj = [ xIn(j), yIn(j), zIn(j) ]
				kIn = SphereDistance( xi, xj ) / pseSetup%eps
				kernelD1 = bivariateFirstDerivativeKernel8( kIn ) / pseSetup%eps / pseSetup%eps
				kernelLap = bivariateLaplacianKernel8( kIn ) / pseSetup%eps / pseSetup%eps
				ux = ux + ( uIn(j) + uIn(i) ) * ( xIn(i) - xIn(j) ) / pseSetup%eps * kernelD1 * areaIn(j)
				uy = uy + ( uIn(j) + uIn(i) ) * ( yIn(i) - yIn(j) ) / pseSetup%eps * kernelD1 * areaIn(j)
				uz = uz + ( uIn(j) + uIn(i) ) * ( zIn(i) - zIn(j) ) / pseSetup%eps * kernelD1 * areaIn(j)
				vx = vx + ( vIn(j) + vIn(i) ) * ( xIn(i) - xIn(j) ) / pseSetup%eps * kernelD1 * areaIn(j)
				vy = vy + ( vIn(j) + vIn(i) ) * ( yIn(i) - yIn(j) ) / pseSetup%eps * kernelD1 * areaIn(j)
				vz = vz + ( vIn(j) + vIn(i) ) * ( zIn(i) - zIn(j) ) / pseSetup%eps * kernelD1 * areaIn(j)
				wx = wx + ( wIn(j) + wIn(i) ) * ( xIn(i) - xIn(j) ) / pseSetup%eps * kernelD1 * areaIn(j)
				wy = wy + ( wIn(j) + wIn(i) ) * ( yIn(i) - yIn(j) ) / pseSetup%eps * kernelD1 * areaIn(j)
				wz = wz + ( wIn(j) + wIn(i) ) * ( zIn(i) - zIn(j) ) / pseSetup%eps * kernelD1 * areaIn(j)
				lapH = lapH + ( hIn(j) - hIn(i) ) * kernelLap / pseSetup%eps / pseSetup%eps * areaIn(j)
			endif
		enddo
		gradU = [ux, uy, uz]
		gradV = [vx, vy, vz]
		gradW = [wx, wy, wz]
		gradU = MATMUL(projMat, gradU)
		gradV = MATMUL(projMat, gradV)
		gradW = MATMUL(projMat, gradW)
		doubleDot = gradU(1)*gradU(1) + gradV(2)*gradV(2) + gradW(3)*gradW(3) + 2.0_kreal * &
			( gradU(2)*gradV(1) + gradU(3)*gradW(1) + gradV(3)*gradW(2) )
		doubleDot = doubleDot / pseSetup%eps / pseSetup%eps	
		lapH = lapH / pseSetup%eps / pseSetup%eps
		dDiv(i) = 2.0_kreal * rotRate * zIn(i) / sphereRadius * relVortIn(i) - doubleDot - g * lapH
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST( dDiv(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	enddo
end subroutine

subroutine InitLogger(aLog,rank)
! Initialize a logger for this module and processor
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
