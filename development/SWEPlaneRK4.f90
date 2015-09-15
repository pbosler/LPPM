module SWEPlaneRK4Module

use NumberKindsModule
use LoggerModule
use ParticlesModule
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use PSEDirectSumModule
use PlanarSWEModule

implicit none

include 'mpif.h'

private
public SWEPlaneRK4, New, Delete
public Timestep
!
!----------------
! types and module variables
!----------------
!
type SWEPlaneRK4
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
	real(kreal), pointer :: u(:) => null()
	real(kreal), pointer :: v(:) => null()
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
! interfaces
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
	function LapDepthFn( x, y )
		real(8) :: LapDepthFn
		real(8), intent(in) :: x, y
	end function
end interface 


!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'SWEPlaneRK4'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains

!
!----------------
! public methods
!----------------
!
subroutine newPrivate(self, planeSWE )
	type(SWEPlaneRK4), intent(out) :: self
	type(SWEMesh), intent(in) :: planeSWE
	
	if ( .NOT. logInit) call InitLogger(log, procRank)
	
	allocate(self%xIn( planeSWE%mesh%particles%N ))
	allocate(self%xStage1( planeSWE%mesh%particles%N))
	allocate(self%xStage2( planeSWE%mesh%particles%N))
	allocate(self%xStage3( planeSWE%mesh%particles%N))
	allocate(self%xStage4( planeSWE%mesh%particles%N))
	allocate(self%yIn( planeSWE%mesh%particles%N ))
	allocate(self%yStage1( planeSWE%mesh%particles%N))
	allocate(self%yStage2( planeSWE%mesh%particles%N))
	allocate(self%yStage3( planeSWE%mesh%particles%N))
	allocate(self%yStage4( planeSWE%mesh%particles%N))
	allocate(self%u( planeSWE%mesh%particles%N))
	allocate(self%v( planeSWE%mesh%particles%N))
	allocate(self%relVortIn( planeSWE%mesh%particles%N ))
	allocate(self%relVortStage1( planeSWE%mesh%particles%N))
	allocate(self%relVortStage2( planeSWE%mesh%particles%N))
	allocate(self%relVortStage3( planeSWE%mesh%particles%N))
	allocate(self%relVortStage4( planeSWE%mesh%particles%N))
	allocate(self%divIn( planeSWE%mesh%particles%N ))
	allocate(self%divStage1( planeSWE%mesh%particles%N))
	allocate(self%divStage2( planeSWE%mesh%particles%N))
	allocate(self%divStage3( planeSWE%mesh%particles%N))
	allocate(self%divStage4( planeSWE%mesh%particles%N))
	allocate(self%hIn( planeSWE%mesh%particles%N ))
	allocate(self%hStage1( planeSWE%mesh%particles%N))
	allocate(self%hStage2( planeSWE%mesh%particles%N))
	allocate(self%hStage3( planeSWE%mesh%particles%N))
	allocate(self%hStage4( planeSWE%mesh%particles%N))
	allocate(self%areaIn( planeSWE%mesh%particles%N ))
	allocate(self%areaStage1( planeSWE%mesh%particles%N))
	allocate(self%areaStage2( planeSWE%mesh%particles%N))
	allocate(self%areaStage3( planeSWE%mesh%particles%N))
	allocate(self%areaStage4( planeSWE%mesh%particles%N))
end subroutine

subroutine deletePrivate( self )
	type(SWEPlaneRK4), intent(inout) :: self
	
	if ( associated( self%xIn))	deallocate(self%xIn)
	if ( associated( self%xStage1)) 	deallocate(self%xStage1)
	if ( associated( self%xStage2))	deallocate(self%xStage2)
	if ( associated( self%xStage3))	deallocate(self%xStage3)
	if ( associated( self%xStage4))	deallocate(self%xStage4)
	if ( associated( self%yIn))	deallocate(self%yIn)
	if ( associated( self%yStage1))	deallocate(self%yStage1)
	if ( associated( self%yStage1))	deallocate(self%yStage2)
	if ( associated( self%yStage2))	deallocate(self%yStage3)
	if ( associated( self%yStage3))	deallocate(self%yStage4)
	if ( associated( self%u))	deallocate(self%u)
	if ( associated( self%v))	deallocate(self%v)
	if ( associated( self%relVortIn))	deallocate(self%relVortIn)
	if ( associated( self%relVortStage1))	deallocate(self%relVortStage1)
	if ( associated( self%relVortStage2))	deallocate(self%relVortStage2)
	if ( associated( self%relVortStage3))	deallocate(self%relVortStage3)
	if ( associated( self%relVortStage4))	deallocate(self%relVortStage4)
	if ( associated( self%divIn))	deallocate(self%divIn)
	if ( associated( self%divStage1))	deallocate(self%divStage1)
	if ( associated( self%divStage2))	deallocate(self%divStage2)
	if ( associated( self%divStage3))	deallocate(self%divStage3)
	if ( associated( self%divStage4))	deallocate(self%divStage4)
	if ( associated( self%hIn))	deallocate(self%hIn)
	if ( associated( self%hStage1))	deallocate(self%hStage1)
	if ( associated( self%hStage2))	deallocate(self%hStage2)
	if ( associated( self%hStage3))	deallocate(self%hStage3)
	if ( associated( self%hStage4))	deallocate(self%hStage4)
	if ( associated( self%areaIn))	deallocate(self%areaIn)
	if ( associated( self%areaStage1)) 	deallocate(self%areaStage1)
	if ( associated( self%areaStage2))	deallocate(self%areaStage2)
	if ( associated( self%areaStage3))	deallocate(self%areaStage3)
	if ( associated( self%areaStage4))	deallocate(self%areaStage4)
end subroutine

subroutine timestepPrivate( self, planarSWE, dt, depthFn, depthLaplacian, pseSetup, mpiParticles )
	type(SWEPlaneRK4), intent(inout) :: self
	type(SWEMesh), intent(inout) :: planarSWE
	real(kreal), intent(in) :: dt
	procedure(LapDepthFn) :: depthFn
	procedure(LapDepthFn) :: depthLaplacian
	type(PSE), intent(in) :: pseSetup
	type(MPISetup), intent(in) :: mpiParticles
	!
	integer(kint) :: i, nP
	
	nP = planarSWE%mesh%particles%N
	!
	!	RK Stage 1 : input = current mesh state
	!
	call SWEPlaneRHS( self%relVortStage1, self%divStage1, self%areaStage1, self%hStage1, &
					  planarSWE%mesh%particles%x(1:nP), planarSWE%mesh%particles%y(1:nP), planarSWE%velocity%xComp(1:nP), &
					  planarSWE%velocity%yComp(1:nP), planarSWE%relvort%scalar(1:nP), planarSWE%divergence%scalar(1:nP), &
					  planarSWE%mesh%particles%area(1:nP), planarSWE%h%scalar(1:nP), depthLaplacian, &
					  planarSWE%f0, planarSWE%beta, planarSWE%g, planarSWE%mesh%particles%isActive(1:nP), &
					  pseSetup, mpiParticles)
	self%xStage1 = dt * planarSWE%velocity%xComp(1:nP)
	self%yStage1 = dt * planarSWE%velocity%yComp(1:nP)
	!self%relVortStage1 = dt * self%relVortStage1
	do i = 1, nP
		self%relVortStage1 = dt * ( planarSWE%h%scalar(i) * planarSWE%potVort%scalar(i) - &
			planarSWE%f0 - planarSWE%beta * planarSWE%mesh%particles%y(i) )
	enddo
	self%divStage1 = dt * self%divStage1
	self%areaStage1 = dt * self%areaStage1
	self%hStage1 = dt * self%hStage1

	!
	!	RK Stage 2 : input = mesh state + 0.5 * stage 1
	!	
	self%xIn = planarSWE%mesh%particles%x(1:nP) + 0.5_kreal * self%xStage1
	self%yIn = planarSWE%mesh%particles%y(1:nP) + 0.5_kreal * self%yStage1
	self%relVortIn = planarSWE%relVort%scalar(1:nP) + 0.5_kreal * self%relvortStage1
	self%divIn = planarSWE%divergence%scalar(1:nP) + 0.5_kreal * self%divStage1
	self%hIn = planarSWE%h%scalar(1:nP) + 0.5_kreal * self%hStage1
	self%areaIn = planarSWE%mesh%particles%area(1:nP) + 0.5_kreal * self%areaStage1
	do i = 1, nP
		if (self%hIn(i) - depthFn(self%xIn(i), self%yIn(i)) <= 0.0_kreal ) then
			self%relVortIn(i) = 0.0_kreal
			self%divIn(i) = 0.0_kreal
			self%u(i) = 0.0_kreal
			self%v(i) = 0.0_kreal
		endif
	enddo
	
	call SWEPlaneVelocity( self%u, self%v, self%xIn, self%yIn, self%relVortIn, self%divIn, self%areaIn, &
						   planarSWE%mesh%particles%isActive, mpiParticles)
	call SWEPlaneRHS( self%relVortStage2, self%divStage2, self%areaStage2, self%hStage2, &
					  self%xIn, self%yIn, self%u, self%v, self%relVortIn, self%divIn, self%areaIn, self%hIn, &
					  depthLaplacian, planarSWE%f0, planarSWE%beta, planarSWE%g, planarSWE%mesh%particles%isActive, &
					  pseSetup, mpiParticles)
	self%xstage2 = dt * self%u
	self%ystage2 = dt * self%v
	!self%relVortstage2 = dt * self%relVortstage2
	do i = 1, nP
		self%relVortStage2(i) = dt * ( self%hIn(i) * planarSWE%potVort%scalar(i) - planarSWE%f0 - &
			planarSWE%beta * self%yIn(i))
	enddo
	self%divstage2 = dt * self%divstage2
	self%areastage2 = dt * self%areastage2
	self%hstage2 = dt * self%hstage2

	!
	!	RK Stage 3 : mesh state + 0.5 * stage 2
	!
	self%xIn = planarSWE%mesh%particles%x(1:nP) + 0.5_kreal * self%xstage2
	self%yIn = planarSWE%mesh%particles%y(1:nP) + 0.5_kreal * self%ystage2
	self%relVortIn = planarSWE%relVort%scalar(1:nP) + 0.5_kreal * self%relvortstage2
	self%divIn = planarSWE%divergence%scalar(1:nP) + 0.5_kreal * self%divstage2
	self%hIn = planarSWE%h%scalar(1:nP) + 0.5_kreal * self%hstage2
	self%areaIn = planarSWE%mesh%particles%area(1:nP) + 0.5_kreal * self%areastage2	
	do i = 1, nP
		if (self%hIn(i) - depthFn(self%xIn(i), self%yIn(i)) <= 0.0_kreal ) then
			self%relVortIn(i) = 0.0_kreal
			self%divIn(i) = 0.0_kreal
			self%u(i) = 0.0_kreal
			self%v(i) = 0.0_kreal
		endif
	enddo
	
	call SWEPlaneVelocity( self%u, self%v, self%xIn, self%yIn, self%relVortIn, self%divIn, self%areaIn, &
						   planarSWE%mesh%particles%isActive, mpiParticles)
	
	call SWEPlaneRHS( self%relVortStage3, self%divStage3, self%areaStage3, self%hStage3, &
					  self%xIn, self%yIn, self%u, self%v, self%relVortIn, self%divIn, self%areaIn, self%hIn, &
					  depthLaplacian, planarSWE%f0, planarSWE%beta, planarSWE%g, planarSWE%mesh%particles%isActive, &
					  pseSetup, mpiParticles)		  
	self%xStage3 = dt * self%u
	self%yStage3 = dt * self%v
	!self%relVortStage3 = dt * self%relVortStage3
	do i = 1, nP
		self%relVortStage3(i) = dt * ( self%hIn(i) * planarSWE%potVort%scalar(i) - planarSWE%f0 - &
			planarSWE%beta * self%yIn(i))
	enddo
	self%divStage3 = dt * self%divStage3
	self%areaStage3 = dt * self%areaStage3
	self%hStage3 = dt * self%hStage3			

	!
	!	RK Stage 4 : input = mesh state + stage 3
	!					  
	self%xIn = planarSWE%mesh%particles%x(1:nP) + self%xStage3
	self%yIn = planarSWE%mesh%particles%y(1:nP) + self%yStage3
	self%relVortIn = planarSWE%relVort%scalar(1:nP) + self%relvortStage3
	self%divIn = planarSWE%divergence%scalar(1:nP) + self%divStage3
	self%hIn = planarSWE%h%scalar(1:nP) + self%hStage3
	self%areaIn = planarSWE%mesh%particles%area(1:nP) + self%areaStage3		
	do i = 1, nP
		if (self%hIn(i) - depthFn(self%xIn(i), self%yIn(i)) <= 0.0_kreal ) then
			self%relVortIn(i) = 0.0_kreal
			self%divIn(i) = 0.0_kreal
			self%u(i) = 0.0_kreal
			self%v(i) = 0.0_kreal
		endif
	enddo
					  
	call SWEPlaneVelocity( self%u, self%v, self%xIn, self%yIn, self%relVortIn, self%divIn, self%areaIn, &
						   planarSWE%mesh%particles%isActive, mpiParticles)
	call SWEPlaneRHS( self%relVortStage4, self%divStage4, self%areaStage4, self%hStage4, &
					  self%xIn, self%yIn, self%u, self%v, self%relVortIn, self%divIn, self%areaIn, self%hIn, &
					  depthLaplacian, planarSWE%f0, planarSWE%beta, planarSWE%g, planarSWE%mesh%particles%isActive,&
					  pseSetup, mpiParticles)
	self%xStage4 = dt * self%u
	self%yStage4 = dt * self%v
	!self%relVortStage4 = dt * self%relVortStage4
	do i = 1, nP
		self%relVortStage4(i) = dt * ( self%hIn(i) * planarSWE%potVort%scalar(i) - planarSWE%f0 - &
			planarSWE%beta * self%yIn(i))
	enddo
	self%divStage4 = dt * self%divStage4
	self%areaStage4 = dt * self%areaStage4
	self%hStage4 = dt * self%hStage4
	
	!
	!	Update mesh state
	!
	planarSWE%mesh%particles%x(1:nP) = planarSWE%mesh%particles%x(1:nP) + self%xStage1 / 6.0_kreal + self%xStage2 / 3.0_kreal + &
					self%xStage3 / 3.0_kreal + self%xStage4 / 6.0_kreal
	planarSWE%mesh%particles%y(1:nP) = planarSWE%mesh%particles%y(1:nP) + self%yStage1 / 6.0_kreal + self%yStage2 / 3.0_kreal + &
					self%yStage3 / 3.0_kreal + self%yStage4 / 6.0_kreal
	planarSWE%relvort%scalar(1:nP) = planarSWE%relvort%scalar(1:nP) + self%relvortStage1 / 6.0_kreal + &
					self%relvortStage2 / 3.0_kreal + self%relvortStage3 / 3.0_kreal + self%relvortStage4 / 6.0_kreal
	planarSWE%divergence%scalar(1:nP) = planarSWE%divergence%scalar(1:nP) + self%divStage1 / 6.0_kreal + &
					self%divStage2 / 3.0_kreal + self%divStage3 / 3.0_kreal + self%divStage4 / 6.0_kreal
	planarSWE%mesh%particles%area(1:nP) = planarSWE%mesh%particles%area(1:nP) + self%areaStage1 / 6.0_kreal + &
					self%areaStage2 / 3.0_kreal + self%areaStage3 / 3.0_kreal + self%areaStage4 / 6.0_kreal
	planarSWE%h%scalar(1:nP) = planarSWE%h%scalar(1:nP) + self%hStage1 / 6.0_kreal + self%hStage2 / 3.0_kreal + &
					self%hStage3 / 3.0_kreal + self%hStage4 / 6.0_kreal 
	
	call SWEPlaneVelocity( planarSWE%velocity%xComp(1:nP), planarSWE%velocity%yComp(1:nP), &
			planarSWE%mesh%particles%x(1:nP), planarSWE%mesh%particles%y(1:nP), planarSWE%relVort%scalar(1:nP), &
			planarSWE%divergence%scalar(1:nP), planarSWE%mesh%particles%area(1:nP), &
			planarSWE%mesh%particles%isActive(1:nP), mpiParticles)
	do i = 1, nP
		planarSWE%depth%scalar(i) = depthFn( planarSWE%mesh%particles%x(i), planarSWE%mesh%particles%y(i))
		if ( planarSWE%h%scalar(i) - planarSWE%depth%scalar(i) <= 0.0_kreal ) then
			planarSWE%relVort%scalar(i) = 0.0_kreal
			planarSWE%divergence%scalar(i) = 0.0_kreal
			planarSWE%velocity%xComp(i) = 0.0_kreal
			planarSWE%velocity%yComp(i) = 0.0_kreal
		endif
	enddo
end subroutine

!
!----------------
! private methods
!----------------
!

subroutine SWEPlaneVelocity( uOut, vOut, xIn, yIn, relVortIn, divIn, areaIn, activeMask, mpiParticles )
	real(kreal), intent(out) :: uOut(:)
	real(kreal), intent(out) :: vOut(:)
	real(kreal), intent(in) :: xIn(:)
	real(kreal), intent(in) :: yIn(:)
	real(kreal), intent(in) :: relVortIn(:)
	real(kreal), intent(in) :: divIn(:)
	real(kreal), intent(in) :: areaIn(:)
	logical(klog), intent(in) :: activeMask(:)
	type(MPISetup), intent(in) :: mpiParticles
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: denom
	
!	uOut(mpiParticles%indexStart(procRank):mpiParticles%indexEnd(procRank)) = 0.0_kreal
!	vOut(mpiParticles%indexStart(procRank):mpiParticles%indexEnd(procRank)) = 0.0_kreal
	uOut = 0.0_kreal
	vOut = 0.0_kreal
	!
	! singular velocity kernels
	!
	do i = mpiParticles%indexStart(procRank), mpiParticles%indexEnd(procRank)
		do j = 1, i - 1
			if ( activeMask(j) ) then
				denom = 2.0_kreal * PI * ( (xIn(j) - xIn(i))**2 + (yIn(j) - yIn(i))**2 )
				uOut(i) = uOut(i) + (-(yIn(i) - yIn(j))*relVortIn(j) + (xIn(i) - xIn(j))*divIn(j) ) * areaIn(j) / denom
				vOut(i) = vOut(i) + ( (xIn(i) - xIn(j))*relVortIn(j) + (yIn(i) - yIn(j))*divIn(j) ) * areaIn(j) / denom
			endif
		enddo
		do j = i + 1, mpiParticles%n
			if ( activeMask(j) ) then
				denom = 2.0_kreal * PI * ( (xIn(j) - xIn(i))**2 + (yIn(j) - yIn(i))**2 )
				uOut(i) = uOut(i) + (-(yIn(i) - yIn(j))*relVortIn(j) + (xIn(i) - xIn(j))*divIn(j) ) * areaIn(j) / denom
				vOut(i) = vOut(i) + ( (xIn(i) - xIn(j))*relVortIn(j) + (yIn(i) - yIn(j))*divIn(j) ) * areaIn(j) / denom
			endif
		enddo
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST( uOut(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( vOut(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)				
	enddo
end subroutine


subroutine SWEPlaneRHS( relVortOut, divOut, areaOut, hOUt, & 
						xIn, yIn, uIn, vIn, relVortIn, divIn, areaIn, hIn, lapDepth, f0, beta, g, &
						activeMask, pseSetup, mpiParticles )
	real(kreal), intent(out) :: relVortOut(:)
	real(kreal), intent(out) :: divOut(:)
	real(kreal), intent(out) :: areaOut(:)
	real(kreal), intent(out) :: hOut(:)
	real(kreal), intent(in) :: xIn(:)
	real(kreal), intent(in) :: yIn(:)
	real(kreal), intent(in) :: uIn(:)
	real(kreal), intent(in) :: vIn(:)
	real(kreal), intent(in) :: relVortIn(:)
	real(kreal), intent(in) :: divIn(:)
	real(kreal), intent(in) :: areaIn(:)
	real(kreal), intent(in) :: hIn(:)
	procedure(LapDepthFn) :: lapDepth
	real(kreal), intent(in) :: f0, beta, g
	logical(klog), intent(in) :: activeMask(:)
	type(PSE), intent(in) :: pseSetup
	type(MPISetup), intent(in) :: mpiParticles
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: kIn, kernelD1, kernelLap
	real(kreal) :: ux, uy, vx, vy, lapD, lapH, doubleDot
	
	!divOut(mpiParticles%indexStart(procRank):mpiParticles%indexEnd(procRank)) = 0.0_kreal
	
	do i = 1, size(xIn)
		hOut(i) = - hIn(i) * divIn(i)
		relVortOut(i) = - ( relVortIn(i) + f0 + beta + yIn(i) ) * divIn(i) - beta * vIn(i)
	enddo
	do i = 1, size(xIn)
		if ( activeMask(i) ) then
			areaOut(i) = divIn(i) * areaIn(i)
		endif
	enddo
		
	do i = mpiParticles%indexStart(procRank), mpiParticles%indexEnd(procRank)
		lapD = lapDepth( xIn(i), yIn(i) )
		ux = 0.0_kreal
		uy = 0.0_kreal
		vx = 0.0_kreal
		vy = 0.0_kreal
		doubleDot = 0.0_kreal
		lapH = 0.0_kreal
		!
		! pse kernels for spatial derivatives
		!
		do j = 1, size(xIn)
			if ( activeMask(j) ) then
				kIn = sqrt( (xIn(i) - xIn(j))**2 + (yIn(i) - yIn(j))**2 ) / pseSetup%eps
				kernelD1 = bivariateFirstDerivativeKernel8( kIn ) / pseSetup%eps**2
				kernelLap = bivariateLaplacianKernel8( kIn ) / pseSetup%eps**2
				ux = ux + ( uIn(j) + uIn(i) ) * ( xIn(i) - xIn(j) ) / pseSetup%eps * kernelD1 * areaIn(j)
				uy = uy + ( uIn(j) + uIn(i) ) * ( yIn(i) - yIn(j) ) / pseSetup%eps * kernelD1 * areaIn(j)
				vx = vx + ( vIn(j) + vIn(i) ) * ( xIn(i) - xIn(j) ) / pseSetup%eps * kernelD1 * areaIn(j)
				vy = vy + ( vIn(j) + vIn(i) ) * ( yIn(i) - yIn(j) ) / pseSetup%eps * kernelD1 * areaIn(j)
				lapH = lapH + ( hIn(j) - hIn(i) ) * kernelLap * areaIn(j)
			endif
		enddo		
		doubleDot = (ux*ux + 2.0_kreal*uy*vx + vy*vy) / pseSetup%eps**2
		lapH = lapH / pseSetup%eps**2
		divOut(i) = (f0 + beta * yIn(i) ) * relVortIn(i) - doubleDot - g * lapD - g * lapH
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST( divOut(mpiParticles%indexStart(i):mpiParticles%indexEnd(i)), mpiParticles%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode )								
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