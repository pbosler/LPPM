module SWEDirectSumModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	Defines the RK4 data structure and methods used by SphereMesh.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
!----------------
! USAGE :  This module provides methods for integrating the shallow water equations on the sphere.
!----------------
use NumberKindsModule
use LoggerModule
use SphereGeomModule
use ParticlesModule
use EdgesModule
use PanelsModule
use SphereMeshModule
use STRIPACKInterfaceModule
use SSRFPACKInterfaceModule

implicit none

include 'mpif.h'

private
public SWERK4Data
public New, Delete
public SWERK4Timestep
public VELOCITY_SMOOTH, SetVelocitySmoothingParameter
!
!----------------
! Types and module constants
!----------------
!
real(kreal), protected, save :: VELOCITY_SMOOTH = 0.01_kreal

type SWERK4Data
	!
	!	MPI / load balancing variables
	!
	integer(kint), pointer :: particlesIndexStart(:) => null(), &
							  particlesIndexEnd(:) => null(), &
							  particlesMessageSize(:) => null(), &
							  activePanelsIndexStart(:) => null(), &
							  activePanelsIndexEnd(:) => null(), &
							  activePanelsMessageSize(:) => null(), &
							  passivePanelsIndexStart(:) => null(), &
							  passivePanelsIndexEnd(:) => null(), &
							  passivePanelsMessageSize(:) => null()
	logical(klog) :: rk4IsReady, &
					 mpiIsReady

	!
	! parameters
	!
	real(kreal) :: smooth

	!
	! mesh & physical variables
	!
	type(Panels), pointer :: activePanels => null(), &
							 passivePanels => null()
	integer(kint), pointer :: activeMap(:) => null(), &
							  passiveMap(:) => null()
	type(STRIPACKData), pointer :: delTri => null()
	type(SSRFPACKData), pointer :: velocitySource1 => null(), &
								   thicknessSource1 => null(), &
								   thicknessSource2 => null()
	!
	! Runge-Kutta variables
	!
	real(kreal), pointer :: particlesInput(:,:) => null(), &
							particlesStage1(:,:) => null(), &
							particlesStage2(:,:) => null(), &
							particlesStage3(:,:) => null(), &
							particlesStage4(:,:) => null(), &
							newParticlesX(:,:) => null(), &
							&
							particlesVelocity(:,:) => null(), &
							particlesDoubleDotU(:) => null(), &
							particlesLapH(:) => null(), &
							&
							particlesRelVortInput(:) => null(), &
							particlesRelVortStage1(:) => null(), &
							particlesRelVortStage2(:) => null(), &
							particlesRelVortStage3(:) => null(), &
							particlesRelVortStage4(:) => null(), &
							newparticlesRelVort(:) => null(), &
							&
							particlesDivInput(:) => null(), &
							particlesDivStage1(:) => null(), &
							particlesDivStage2(:) => null(), &
							particlesDivStage3(:) => null(), &
							particlesDivStage4(:) => null(), &
							newparticlesDiv(:) => null(), &
							&
							particlesHInput(:) => null(), &
							particlesHStage1(:) => null(), &
							particlesHStage2(:) => null(), &
							particlesHStage3(:) => null(), &
							particlesHStage4(:) => null(), &
							newParticlesH(:) => null(), &
							&					
							activePanelsInput(:,:) => null(), &
							activePanelsStage1(:,:) => null(), &
							activePanelsStage2(:,:) => null(), &
							activePanelsStage3(:,:) => null(), &
							activePanelsStage4(:,:) => null(), &
							newActivePanelsX(:,:) => null(), &
							&
							activePanelsVelocity(:,:) => null(), &
							activePanelsDoubleDotU(:) => null(), &
							activePanelsLapH(:) => null(), &
							&
							activeRelVortInput(:) => null(), &
							activeRelVortStage1(:) => null(), &
							activeRelVortStage2(:) => null(), &
							activeRelVortStage3(:) => null(), &
							activeRelVortStage4(:) => null(), &
							newActiveRelVort(:) => null(), &
							&
							activeDivInput(:) => null(), &
							activeDivStage1(:) => null(), &
							activeDivStage2(:) => null(), &
							activeDivStage3(:) => null(), &
							activeDivStage4(:) => null(), &
							newActiveDiv(:) => null(), &
							&
							activeHInput(:) => null(), &
							activeHStage1(:) => null(), &
							activeHStage2(:) => null(), &
							activeHStage3(:) => null(), &
							activeHStage4(:) => null(), &
							newActiveH(:) => null(), &
							&
							activeAreaInput(:) => null(), &
							activeAreaStage1(:) => null(), &
							activeAreaStage2(:) => null(), &
							activeAreaStage3(:) => null(), &
							activeAreaStage4(:) => null(), &
							newActiveArea(:) => null(), &
							&
							passivePanelsInput(:,:) => null(), &
							passivePanelsStage1(:,:) => null(), &
							passivePanelsStage2(:,:) => null(), &
							passivePanelsStage3(:,:) => null(), &
							passivePanelsStage4(:,:) => null(), &
							newPassivePanelsX(:,:) => null()
end type

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'SWE'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=128) :: logstring
character(len=24) :: formatstring
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
subroutine NewPrivate(self, aMesh, nProcs)
	type(SWERK4Data), intent(out) :: self
	type(SphereMesh), intent(in) :: aMesh
	integer(kint), intent(in) :: nProcs
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	!
	integer(kint) :: nParticles, nActive, nPassive, nTracer, panelKind, problemKind

	if ( .NOT. logInit ) call InitLogger(log,procRank)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'creating new SWERK4Data.')
	!
	! error checking
	!
	if ( aMesh%problemKind /= SWE_SOLVER) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'ERROR : SWERK4 requires problemKind = SWE_SOLVER')
		return
	endif

	!
	! initialize to 0 / null
	!
	self%rk4IsReady = .FALSE.
	self%mpiIsReady = .FALSE.

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	nParticles = aParticles%N
	nActive = aPanels%N_Active
	nPassive = aPanels%N - nActive

	nTracer = aMesh%nTracer

	panelKind = aMesh%panelKind

	problemKind = SWE_SOLVER

	!
	! allocate MPI variables
	!
	allocate(self%particlesIndexStart(0:nProcs-1))
	allocate(self%particlesIndexEnd(0:nProcs-1))
	allocate(self%particlesMessageSize(0:nProcs-1))
	self%particlesIndexStart = -1
	self%particlesIndexEnd = -1
	self%particlesMessageSize = -1
	allocate(self%activePanelsIndexStart(0:nProcs-1))
	allocate(self%activePanelsIndexEnd(0:nProcs-1))
	allocate(self%activePanelsMessageSize(0:nProcs-1))
	self%activePanelsIndexStart = -1
	self%activePanelsIndexEnd = -1
	self%activePanelsMessageSize = -1
	allocate(self%passivePanelsIndexStart(0:nProcs-1))
	allocate(self%passivePanelsIndexEnd(0:nProcs-1))
	allocate(self%passivePanelsMessageSize(0:nProcs-1))
	self%passivePanelsIndexStart = -1
	self%passivePanelsIndexEnd = -1
	self%passivePanelsMessageSize = -1
	!
	! set physical parameters
	!
	self%smooth = VELOCITY_SMOOTH

	!
	! separate activePanels from passivePanels
	!
	allocate(self%activePanels)
	call New(self%activePanels,nActive,panelKind,nTracer,problemKind)
	self%activePanels%N = nActive
	self%activePanels%N_Active = nActive
	allocate(self%activeMap(nActive))
	self%activeMap = 0

	allocate(self%passivePanels)
	call New(self%passivePanels,nPassive,panelKind,nTracer,problemKind)
	self%passivePanels%N = nPassive
	allocate(self%passiveMap(nPassive))
	self%passiveMap = 0

	call GatherPanels(aPanels,self%activePanels,self%activeMap, self%passivePanels,self%passiveMap)

	!
	! find Delaunay triangulation of particles, allocate interpolation sources for derivatives
	!
	allocate(self%delTri)
	call New(self%delTri,aMesh)
	call DelaunayTriangulation(self%delTri)

	allocate(self%velocitySource1)
	call New(self%velocitySource1,self%delTri,.TRUE.)

	allocate(self%thicknessSource1)
	call New(self%thicknessSource1,self%delTri,.FALSE.)
	allocate(self%thicknessSource2)
	call New(self%thicknessSource2,self%delTri,.FALSE.)

	!
	! allocate RK4 variables
	!
	allocate(self%particlesInput(3,nParticles))
	allocate(self%particlesStage1(3,nParticles))
	allocate(self%particlesStage2(3,nParticles))
	allocate(self%particlesStage3(3,nParticles))
	allocate(self%particlesStage4(3,nParticles))
	allocate(self%newParticlesX(3,nParticles))

	allocate(self%activePanelsInput(3,nActive))
	allocate(self%activePanelsStage1(3,nActive))
	allocate(self%activePanelsStage2(3,nActive))
	allocate(self%activePanelsStage3(3,nActive))
	allocate(self%activePanelsStage4(3,nActive))
	allocate(self%newActivePanelsX(3,nActive))

	allocate(self%passivePanelsInput(3,nPassive))
	allocate(self%passivePanelsStage1(3,nPassive))
	allocate(self%passivePanelsStage2(3,nPassive))
	allocate(self%passivePanelsStage3(3,nPassive))
	allocate(self%passivePanelsStage4(3,nPassive))
	allocate(self%newPassivePanelsX(3,nPassive))

	allocate(self%activeVelocityInput(3,nActive))
	allocate(self%activeVelocityStage1(3,nActive))
	allocate(self%activeVelocityStage2(3,nActive))
	allocate(self%activeVelocityStage3(3,nActive))
	allocate(self%activeVelocityStage4(3,nActive))
	allocate(self%newActiveVelocity(3,nActive))

	allocate(self%activeRelVortInput(nActive))
	allocate(self%activeRelVortStage1(nActive))
	allocate(self%activeRelVortStage2(nActive))
	allocate(self%activeRelVortStage3(nActive))
	allocate(self%activeRelVortStage4(nActive))
	allocate(self%newActiveRelVort(nActive))

	allocate(self%activeDivInput(nActive))
	allocate(self%activeDivStage1(nActive))
	allocate(self%activeDivStage2(nActive))
	allocate(self%activeDivStage3(nActive))
	allocate(self%activeDivStage4(nActive))
	allocate(self%newActiveDiv(nActive))

	allocate(self%activeHInput(nActive))
	allocate(self%activeHStage1(nActive))
	allocate(self%activeHStage2(nActive))
	allocate(self%activeHStage3(nActive))
	allocate(self%activeHStage4(nActive))
	allocate(self%newActiveH(nActive))

	allocate(self%activeAreaInput(nActive))
	allocate(self%activeAreaStage1(nActive))
	allocate(self%activeAreaStage2(nActive))
	allocate(self%activeAreaStage3(nActive))
	allocate(self%activeAreaStage4(nActive))
	allocate(self%newActiveArea(nActive))

	self%rk4IsReady = .TRUE.
	call ZeroRK4(self)

	call LoadBalance(self,nParticles,nActive,nPassive,nProcs)

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'RK4 data ready.')
end subroutine

subroutine DeletePrivate(self)
	type(SWERK4Data), intent(inout) :: self
	
	self%rk4IsReady = .FALSE.
	self%mpiIsReady = .FALSE.
	
	call Delete(self%activePanels)
	call Delete(self%passivePanels)
	deallocate(self%activeMap)
	deallocate(self%passiveMap)
	deallocate(self%activePanels)
	deallocate(self%passivePanels)
	
	call Delete(self%delTri)
	call Delete(self%velocitySource1)
	call Delete(self%thicknessSource1)
	call Delete(self%thicknessSource2)
	deallocate(self%delTri)
	deallocate(self%velocitySource1)
	deallocate(self%thicknessSource1)
	deallocate(self%thicknessSource2)
	
	deallocate(self%particlesIndexStart)
	deallocate(self%particlesIndexEnd)
	deallocate(self%particlesMessageSize)
	deallocate(self%activePanelsIndexStart)
	deallocate(self%activePanelsIndexEnd)
	deallocate(self%activePanelsMessageSize)
	deallocate(self%passivePanelsIndexStart)
	deallocate(self%passivePanelsIndexEnd)
	deallocate(self%passivePanelsMessageSize)
	
	deallocate(self%particlesInput)
	deallocate(self%particlesStage1)
	deallocate(self%particlesStage2)
	deallocate(self%particlesStage3)
	deallocate(self%particlesStage4)
	deallocate(self%newParticlesX)

	deallocate(self%activePanelsInput)
	deallocate(self%activePanelsStage1)
	deallocate(self%activePanelsStage2)
	deallocate(self%activePanelsStage3)
	deallocate(self%activePanelsStage4)
	deallocate(self%newActivePanelsX)

	deallocate(self%passivePanelsInput)
	deallocate(self%passivePanelsStage1)
	deallocate(self%passivePanelsStage2)
	deallocate(self%passivePanelsStage3)
	deallocate(self%passivePanelsStage4)
	deallocate(self%newPassivePanelsX)

	deallocate(self%activeVelocityInput)
	deallocate(self%activeVelocityStage1)
	deallocate(self%activeVelocityStage2)
	deallocate(self%activeVelocityStage3)
	deallocate(self%activeVelocityStage4)
	deallocate(self%newActiveVelocity)

	deallocate(self%activeRelVortInput)
	deallocate(self%activeRelVortStage1)
	deallocate(self%activeRelVortStage2)
	deallocate(self%activeRelVortStage3)
	deallocate(self%activeRelVortStage4)
	deallocate(self%newActiveRelVort)

	deallocate(self%activeDivInput)
	deallocate(self%activeDivStage1)
	deallocate(self%activeDivStage2)
	deallocate(self%activeDivStage3)
	deallocate(self%activeDivStage4)
	deallocate(self%newActiveDiv)

	deallocate(self%activeHInput)
	deallocate(self%activeHStage1)
	deallocate(self%activeHStage2)
	deallocate(self%activeHStage3)
	deallocate(self%activeHStage4)
	deallocate(self%newActiveH)

	deallocate(self%activeAreaInput)
	deallocate(self%activeAreaStage1)
	deallocate(self%activeAreaStage2)
	deallocate(self%activeAreaStage3)
	deallocate(self%activeAreaStage4)
	deallocate(self%newActiveArea)
end subroutine
!
!----------------
! Public member functions
!----------------
!
subroutine SWERK4Timestep(self, aMesh, dt, procRank, nProcs)
	type(SWERK4Data), intent(inout) :: self
	type(SphereMesh), intent(inout) :: aMesh
	real(kreal), intent(in) :: dt
	integer(kint), intent(in) :: procRank, nProcs
	!
	type(Particles), pointer :: aParticles
	integer(kint) :: j, errCode
	
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'entering SWERK4Timestep.')
	
	if ( self%rk4isReady .AND. self%mpiIsReady) then
		call ZeroRK4(self)
	else
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'SWERK4 ERROR : data not ready.')
		return
	endif
	
	aParticles => aMesh%particles
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 1       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 1
	self%particlesInput = aParticles%x(:,1:aParticles%N)
	self%particlesRelVortInput = aParticles%relVort(1:aParticles%N)
	self%particlesDivInput = aParticles%div(1:aParticles%N)
	self%particlesHInput = aParticles%h(1:aParticles%N)
	
	self%activePanelsInput = self%activePanels%x
	self%activeRelVortInput = self%activePanels%relVort
	self%activeDivInput = self%activePanels%div
	self%activeHInput = self%activePanels%h
	self%activeAreaInput = self%activePanels%area
	
	self%passivePanelsInput = self%passivePanels%x
	!
	! PARALLEL : compute particle velocities
	!
	call SWEVelocityActive(self%activePanelsVelocity, self%activePanelsInput, &
				self%activeRelVortInput, self%activeDivInput, self%activeAreaInput, &
				self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank))
	call SWEVelocityPassive(self%particlesVelocity, self%particlesInput, self%activePanelsInput, &
				self%activeRelVortInput, self%activeDivInput, self%activeAreaInput, &
				self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))
	call SWEVelocitySmooth(self%passivePanelsStage1, self%passivePanelsInput, self%activePanelsInput, &
				self%activeRelVortInput, self%activeDivInput, self%activeAreaInput, &
				self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank))
	!
	! broadcast velocities (complete set needed for divergence equation)
	!			
	do j=0, nProcs-1
		call MPI_BCAST(self%particlesVelocity(:,self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
					   3*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsVelocity(:,self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
					   3*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)	
		call MPI_BCAST(self%passivePanelsStage1(:,self%passivePanelsIndexStart(j):passivePanelsIndexEnd(j)), &
					   3*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)			   
	enddo
	!
	! END PARALLEL
	!
	call SetSourceVelocity(self%velocitySource, self%delTri, self%particlesVelocity, self%activePanelsVelocity)
	call SetSourceH(self%thicknessSource, self%delTri, self%particlesHInput, self%activePanelsHInput)
	
	!
	! PARALLEL : compute divergence equation forcing terms
	!
	call ComputeDoubleDotU(self%particlesDoubleDotU, self%delTri, self%velocitySource, &
		self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))
	call ComputeLaplacianH(self%particlesLapH, self%delTri, self%hSource1, self%hSource2, &	
		self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))

	call ComputeDoubleDotU(self%activePanelsDoubleDotU, self%delTri, self%velocitySource, &
		self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank))	
	call ComputeLaplacianH(self%activePanelsLapH, self%delTri, self%hSource1, self%hSource2, &
		self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))
	!
	! compute RHS for vorticity, divergence, h, and area
	!
	do j=self%particlesIndexStart(procRank),self%particlesIndexEnd(procRank)
		self%particlesRelVortStage1(j) = -(self%particlesRelVortInput(j) + 2.0_kreal*OMEGA*self%particlesInput(3,j)/EARTH_RADIUS)*&
					self%particlesDivInput(j) - 2.0_kreal*OMEGA*self%particlesVelocityStage1(j)/EARTH_RADIUS
		self%particlesDivStage1(j) = -self%particlesDoubleDotU(j) + &
			2.0_kreal*OMEGA*self%particlesInput(3,j)/EARTH_RADIUS*self%particlesRelVortInput(j) - &
			GRAV*self%particlesLapH(j)
		self%particlesHStage1(j) = -self%particlesDivInput(j)*self%particlesHInput(j)
	enddo
	
	do j=self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank)
		self%activeRelVortStage1(j) = -(self%activeRelVortInput(j) + 2.0_kreal*OMEGA*self%activePanelsInput(3,j)/EARTH_RADIUS)*&
			self%activeDivInput(j) - 2.0_kreal*OMEGA*self%activeVelocityStage1(3,j)/EARTH_RADIUS
		self%activeDivStage1(j) = - self%activePanelsDoubleDotU(j) + &
			2.0_kreal*OMEGA*self%activePanelsInput(3,j)/EARTH_RADIUS*self%activeRelVortInput(j) - &
			GRAV*self%activePanelsLapH(j)
		self%activeHStage1(j) = -self%activeDivInput(j)*self%activeHInput(j)
		self%activeAreaStage1(j) = self%activeDivInput(j)*self%activeAreaInput(j)	
	enddo
	
	do j=1, nProcs-1
		!
		! broadcast stage1 particles' flow quantities 
		!
		call MPI_BCAST(self%particlesRelVortStage1(self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
						self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%particlesDivStage1(self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
						self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%particlesHStage1(self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
						self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)							   
		
		!
		! broadcast stage1 active panels' flow quantities
		!
		call MPI_BCAST(self%activeRelVortStage1(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activeDivStage1(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activeHStage1(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activeAreaStage1(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)												
	enddo
	!
	! END PARALLEL
	!
	
	!
	! STAGE 1 ONLY : Store velocity and kinetic energy
	!
	do j=1,aParticles%N
		aParticles%ke(j) = sum(self%particlesVelocity(:,j)*self%particlesVelocity(:,j))
		aParticles%u(:,j) = self%particlesVelocity(:,j)
	enddo
	do j=1,self%activePanels%N
		self%activePanels%ke(j) = sum(self%activePanelsVelocity(:,j)*self%activePanelsVelocity(:,j))
		self%activePanels%u(:,j) = self%activePanelsVelocity(:,j)
	enddo

	self%particlesStage1 = dt*self%particlesVelocity
	self%particlesRelVortStage1 = dt*self%particlesRelVortStage1
	self%particlesDivStage1 = dt*self%particlesDivStage1
	self%particlesHStage1 = dt*self%particlesHStage1
	
	self%activePanelsStage1 = dt*self%activePanelsVelocity
	self%activePanelsRelVortStage1 = dt*self%activePanelsRelVortStage1
	self%activePanelsDivStage1 = dt*self%activePanels%
		
	self%passivePanelsStage1 = dt*self%passivePanelsStage1

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 2       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 2
	self%particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*self%particlesStage1
	self%particlesRelVortInput = aParticles%relVort(1:aParticles%N) + 0.5_kreal*self%particlesRelVortStage1
	self%particlesDivInput = aParticles%div(1:aParticles%N) + 0.5_kreal*self%particlesDivStage1
	self%particlesHInput = aParticles%h(1:aParticles%N) + 0.5_kreal*self%particlesHStage1
	
	
	
	!
	! PARALLEL : compute particle velocities
	!
	
	!
	! broadcast velocities (needed for divergence equation)
	!
	
	!
	! END PARALLEL
	!
	
	! SET SSRFPACK SOURCE TO NEW VELOCITY (not mesh velocity)
	! SET SSRFPACK H SOURCE TO HINPUT (not mesh h)
	
	!
	! PARALLEL : compute divergence equation forcing terms
	!
	
	!
	! compute RHS for vorticity, divergence, h, and area
	!

	!
	! broadcast stage2 vorticity, divergence, h, area
	!

	!
	! END PARALLEL
	!

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 3       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 3
	
	!
	! PARALLEL : compute particle velocities
	!
	
	!
	! broadcast velocities (needed for divergence equation)
	!
	
	!
	! END PARALLEL
	!
	
	! SET SSRFPACK SOURCE TO NEW VELOCITY (not mesh velocity)
	! SET SSRFPACK H SOURCE TO HINPUT (not mesh h)
	
	!
	! PARALLEL : compute divergence equation forcing terms
	!
	
	!
	! compute RHS for vorticity, divergence, h, and area
	!

	!
	! broadcast stage3 vorticity, divergence, h, area
	!

	!
	! END PARALLEL
	!	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 4      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 4
	
	!
	! PARALLEL : compute particle velocities
	!
	
	!
	! broadcast velocities (needed for divergence equation)
	!
	
	!
	! END PARALLEL
	!
	
	! SET SSRFPACK SOURCE TO NEW VELOCITY (not mesh velocity)
	! SET SSRFPACK H SOURCE TO HINPUT (not mesh h)
	
	!
	! PARALLEL : compute divergence equation forcing terms
	!
	
	!
	! compute RHS for vorticity, divergence, h, and area
	!

	!
	! broadcast stage4 vorticity, divergence, h, area
	!

	!
	! END PARALLEL
	!
	
end subroutine

subroutine SetVelocitySmoothingParameter(newSmooth)
	real(kreal), intent(in) :: newSmooth
	VELOCITY_SMOOTH = newSmooth
end subroutine

!
!----------------
! Module methods : module- or type-specific private functions
!----------------
!
subroutine ZeroRK4(self)
	type(SWERK4Data), intent(inout) :: self
	self%particlesInput = 0.0_kreal
	self%particlesStage1 = 0.0_kreal
	self%particlesStage2 = 0.0_kreal
	self%particlesStage3 = 0.0_kreal
	self%particlesStage4 = 0.0_kreal
	self%newParticlesX = 0.0_kreal

	self%activePanelsInput = 0.0_kreal
	self%activePanelsStage1 = 0.0_kreal
	self%activePanelsStage2 = 0.0_kreal
	self%activePanelsStage3 = 0.0_kreal
	self%activePanelsStage4 = 0.0_kreal
	self%newActivePanelsX = 0.0_kreal

	self%passivePanelsInput = 0.0_kreal
	self%passivePanelsStage1 = 0.0_kreal
	self%passivePanelsStage2 = 0.0_kreal
	self%passivePanelsStage3 = 0.0_kreal
	self%passivePanelsStage4 = 0.0_kreal
	self%newPassivePanelsX = 0.0_kreal

	self%activeVelocityInput = 0.0_kreal
	self%activeVelocityStage1 = 0.0_kreal
	self%activeVelocityStage2 = 0.0_kreal
	self%activeVelocityStage3 = 0.0_kreal
	self%activeVelocityStage4 = 0.0_kreal
	self%newActiveVelocity = 0.0_kreal

	self%activeRelVortInput = 0.0_kreal
	self%activeRelVortStage1 = 0.0_kreal
	self%activeRelVortStage2 = 0.0_kreal
	self%activeRelVortStage3 = 0.0_kreal
	self%activeRelVortStage4 = 0.0_kreal
	self%newActiveRelVort = 0.0_kreal

	self%activeDivInput = 0.0_kreal
	self%activeDivStage1 = 0.0_kreal
	self%activeDivStage2 = 0.0_kreal
	self%activeDivStage3 = 0.0_kreal
	self%activeDivStage4 = 0.0_kreal
	self%newActiveDiv = 0.0_kreal

	self%activeHInput = 0.0_kreal
	self%activeHStage1 = 0.0_kreal
	self%activeHStage2 = 0.0_kreal
	self%activeHStage3 = 0.0_kreal
	self%activeHStage4 = 0.0_kreal
	self%newActiveH = 0.0_kreal

	self%activeAreaInput = 0.0_kreal
	self%activeAreaStage1 = 0.0_kreal
	self%activeAreaStage2 = 0.0_kreal
	self%activeAreaStage3 = 0.0_kreal
	self%activeAreaStage4 = 0.0_kreal
	self%newActiveArea = 0.0_kreal
end subroutine

subroutine ComputeDoubleDotU(ddU, delTri, velocitySource1, indexStart, indexEnd)
	real(kreal), intent(out) :: ddU(:)
	type(STRIPACKData), intent(in) :: delTri
	type(SSRFPACKData), intent(in) :: velocitySource1
	integer(kint), intent(in) :: indexStart, indexEnd
	!
	integer(kint) :: j, errCode
	
	do j=indexStart, indexEnd
		ddU(j) = velocitySource1%grad1(1,j)*velocitySource1%grad1(1,j) + &
				 velocitySource1%grad2(2,j)*velocitySource1%grad2(2,j) + &
				 velocitySource1%grad3(3,j)*velocitySource1%grad3(3,j) + &
				 2.0_kreal*velocitySource1%grad1(2,j)*velocitySource1%grad2(1,j) + &
				 2.0_kreal*velocitySource1%grad1(3,j)*velocitySource1%grad3(1,j) + &
				 2.0_kreal*velocitySource1%grad2(3,j)*velocitySource1%grad3(2,j)
	enddo
end subroutine

subroutine ComputeLaplacianH(lapH, delTri, hSource1, hsource2, indexStart, indexEnd)
	real(kreal), intent(out) :: lapH(:)
	type(STRIPACKData), intent(in) :: delTri
	type(SSRFPACKData), intent(inout) :: hSource1, hSource2
	integer(kint), intent(in) :: indexStart, indexEnd
	!
	integer(kint) :: j, errCode
	
	do j=1, delTri%n
		hSource2%data1(j) = hSource1%grad1(1,j)
		hSource2%data2(j) = hSource1%grad1(2,j)
		hSource2%data3(j) = hSource1%grad1(3,j)
	enddo
	
	do j = indexStart, indexEnd
		call GRADL(delTri%n, j, delTri%x, delTri%y, delTri%z, hSource2%data1, &
				   delTri%list, delTri%lptr, delTri%lend, hSource2%grad1(:,j), errCode)
		call GRADL(delTri%n, j, delTri%x, delTri%y, delTri%z, hSource2%data2, &
				   delTri%list, delTri%lptr, delTri%lend, hSource2%grad2(:,j), errCode)
		call GRADL(delTri%n, j, delTri%x, delTri%y, delTri%z, hSource2%data3, &
				   delTri%list, delTri%lptr, delTri%lend, hSource2%grad3(:,j), errCode)
		lapH(j) = hSource2%grad1(1,j) + hSource2%grad2(2,j) + hSource2%grad3(3,j)				   		   		   
	enddo
end subroutine

subroutine SWEVelocityActive(u, x, relVort, div, area, indexStart, indexEnd)
	real(kreal), intent(out) :: u(:,:)
	real(kreal), intent(in) :: x(:,:)
	real(kreal), intent(in) :: relVort(:), div(:), area(:)
	integer(kint), intent(in) :: indexStart, indexEnd
	!
	integer(kint) :: j, k, nn
	real(kreal) :: denom
	
	u = 0.0_kreal
	nn = size(x,2)
	do j=indexStart, indexEnd
		do k=1,j-1
			denom = EARTH_RADIUS*EARTH_RADIUS- sum(x(:,j)*x(:,k))
			u(1,j) = u(1,j) + ( x(2,j)*x(3,k) - x(3,j)*x(2,k) )*relVort(k)*area(k)/denom + &
							  ( x(1,j)*x(1,k)/EARTH_RADIUS/EARTH_RADIUS - x(1,k) )*div(k)*area(k)/denom
			u(2,j) = u(2,j) + ( x(3,j)*x(1,k) - x(1,j)*x(3,k) )*relVort(k)*area(k)/denom + &
							  ( x(2,j)*x(2,k)/EARTH_RADIUS/EARTH_RADIUS - x(2,k) )*div(k)*area(k)/denom
			u(3,j) = u(3,j) + ( x(1,j)*x(2,k) - x(2,j)*x(1,k) )*relVort(k)*area(k)/denom + 
					 		  ( x(3,j)*x(3,k)/EARTH_RADIUS/EARTH_RADIUS - x(3,k) )*div(k)*area(k)/denom		  				  
		enddo
		do k=j+1,nn
			denom = EARTH_RADIUS*EARTH_RADIUS- sum(x(:,j)*x(:,k))
			u(1,j) = u(1,j) + ( x(2,j)*x(3,k) - x(3,j)*x(2,k) )*relVort(k)*area(k)/denom + &
							  ( x(1,j)*x(1,k)/EARTH_RADIUS/EARTH_RADIUS - x(1,k) )*div(k)*area(k)/denom
			u(2,j) = u(2,j) + ( x(3,j)*x(1,k) - x(1,j)*x(3,k) )*relVort(k)*area(k)/denom + &
							  ( x(2,j)*x(2,k)/EARTH_RADIUS/EARTH_RADIUS - x(2,k) )*div(k)*area(k)/denom
			u(3,j) = u(3,j) + ( x(1,j)*x(2,k) - x(2,j)*x(1,k) )*relVort(k)*area(k)/denom + 
					 		  ( x(3,j)*x(3,k)/EARTH_RADIUS/EARTH_RADIUS - x(3,k) )*div(k)*area(k)/denom				
		enddo
	enddo
	u(:,indexStart:indexEnd) = u(:,indexStart:indexEnd)/(-4.0_kreal*PI*EARTH_RADIUS)
end subroutine

subroutine SWEVelocityPassive(u, x, xActive, relVort, div, area, indexStart, indexEnd)
	real(kreal), intent(out) :: u(:,:)
	real(kreal), intent(in) :: x(:,:), xActive(:,:)
	real(kreal), intent(in) :: relVort(:), div(:), area(:)
	integer(kint), intent(in) :: indexStart, indexEnd
	!
	integer(kint) :: j, k, nn
	real(kreal) :: denom
	
	u = 0.0_kreal
	nn = size(xActive,2)
	do j=indexStart, indexEnd
		do k=1,nn
			denom = EARTH_RADIUS*EARTH_RADIUS- sum(x(:,j)*x(:,k))
			u(1,j) = u(1,j) + ( x(2,j)*x(3,k) - x(3,j)*x(2,k) )*relVort(k)*area(k)/denom + &
							  ( x(1,j)*x(1,k)/EARTH_RADIUS/EARTH_RADIUS - x(1,k) )*div(k)*area(k)/denom
			u(2,j) = u(2,j) + ( x(3,j)*x(1,k) - x(1,j)*x(3,k) )*relVort(k)*area(k)/denom + &
							  ( x(2,j)*x(2,k)/EARTH_RADIUS/EARTH_RADIUS - x(2,k) )*div(k)*area(k)/denom
			u(3,j) = u(3,j) + ( x(1,j)*x(2,k) - x(2,j)*x(1,k) )*relVort(k)*area(k)/denom + 
					 		  ( x(3,j)*x(3,k)/EARTH_RADIUS/EARTH_RADIUS - x(3,k) )*div(k)*area(k)/denom		  				  
		enddo
	enddo
	u(:,indexStart:indexEnd) = u(:,indexStart:indexEnd)/(-4.0_kreal*PI*EARTH_RADIUS)
end subroutine

subroutine SWEVelocitySmooth(u, x, xActive, relVort, div, area, indexStart, indexEnd)
	real(kreal), intent(out) :: u(:,:)
	real(kreal), intent(in) :: x(:,:), xActive(:,:)
	real(kreal), intent(in) :: relVort(:), div(:), area(:)
	integer(kint), intent(in) :: indexStart, indexEnd
	!
	integer(kint) :: j, k, nn
	real(kreal) :: denom
	
	u = 0.0_kreal
	nn = size(xActive,2)
	do j=indexStart, indexEnd
		do k=1,nn
			denom = EARTH_RADIUS*EARTH_RADIUS- sum(x(:,j)*x(:,k)) + VELOCITY_SMOOTH*VELOCITY_SMOOTH
			u(1,j) = u(1,j) + ( x(2,j)*x(3,k) - x(3,j)*x(2,k) )*relVort(k)*area(k)/denom + &
							  ( x(1,j)*x(1,k)/EARTH_RADIUS/EARTH_RADIUS - x(1,k) )*div(k)*area(k)/denom
			u(2,j) = u(2,j) + ( x(3,j)*x(1,k) - x(1,j)*x(3,k) )*relVort(k)*area(k)/denom + &
							  ( x(2,j)*x(2,k)/EARTH_RADIUS/EARTH_RADIUS - x(2,k) )*div(k)*area(k)/denom
			u(3,j) = u(3,j) + ( x(1,j)*x(2,k) - x(2,j)*x(1,k) )*relVort(k)*area(k)/denom + 
					 		  ( x(3,j)*x(3,k)/EARTH_RADIUS/EARTH_RADIUS - x(3,k) )*div(k)*area(k)/denom		  				  
		enddo
	enddo
	u(:,indexStart:indexEnd) = u(:,indexStart:indexEnd)/(-4.0_kreal*PI*EARTH_RADIUS)
end subroutine

subroutine LoadBalance(self,nParticles,nActive,nPassive,nProcs)
	type(SWERK4Data), intent(inout) :: self
	integer(kint), intent(in) :: nParticles, nActive, nPassive, nProcs
	!
	integer(kint) :: j, chunkSize

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'entering LoadBalance.')

	chunkSize = nParticles/nProcs
	do j=0,nProcs-1
		self%particlesIndexStart(j) = j*chunkSize + 1
		self%particlesIndexEnd(j) = (j+1)*chunkSize
	enddo
	self%particlesIndexEnd(nProcs-1) = nParticles
	self%particlesMessageSize = self%particlesIndexEnd - self%particlesIndexStart + 1

	chunkSize = nActive / nProcs
	do j=0,nProcs-1
		self%activePanelsIndexStart(j) = j*chunkSize + 1
		self%activePanelsIndexEnd(j) = (j+1)*chunkSize
	enddo
	self%activePanelsIndexEnd(nProcs-1) = nActive
	self%activePanelsMessageSize = self%activePanelsIndexEnd - self%activePanelsIndexStart + 1

	chunkSize = nPassive / nProcs
	do j=0,nProcs-1
		self%passivePanelsIndexStart(j) = j*chunkSize + 1
		self%passivePanelsIndexEnd(j) = (j+1)*chunkSize
	enddo
	self%passivePanelsIndexEnd(nProcs-1) = nPassive
	self%passivePanelsMessageSize = self%passivePanelsIndexEnd - self%passivePanelsIndexStart + 1

	self%mpiIsReady = .TRUE.
end subroutine


subroutine InitLogger(aLog,rank)
	type(Logger), intent(inout) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine


end module
