module BVEDirectSumModule
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
public BVERK4Data
public New, Delete
public BVERK4Timestep

!
!----------------
! Types and module constants
!----------------
!
type BVERK4Data
	! MPI load balancing
	integer(kint), pointer :: particlesIndexStart(:) => null(), &
							  particlesIndexEnd(:) =>null(), &
							  particlesMessageSize(:) => null(), &
							  activePanelsIndexStart(:) => null(), &
							  activePanelsIndexEnd(:) => null(), &
							  activePanelsMessageSize(:) => null(), &
							  passivePanelsIndexStart(:) => null() ,&
							  passivePanelsIndexEnd(:) => null(), &
							  passivePanelsMessageSize(:) => null()

	logical(klog) :: rk4isReady, &
				     mpiIsReady
	! parameters
	real(kreal) :: smooth

	! mesh variables
	type(Panels), pointer :: activePanels => null(),&
							 passivePanels => null()
	integer(kint), pointer :: activeMap(:) => null(), &
							  passiveMap(:) => null()

	! Runge-Kutta variables
	real(kreal), pointer :: particlesInput(:,:) => null(), &
							particlesStage1(:,:) => null(), &
							particlesStage2(:,:) => null(), &
							particlesStage3(:,:) => null(), &
							particlesStage4(:,:) => null(), &
							newParticlesX(:,:) => null(), &
							activePanelsInput(:,:) => null(), &
							activePanelsStage1(:,:) => null() , &
							activePanelsStage2(:,:) => null() , &
							activePanelsStage3(:,:) => null() , &
							activePanelsStage4(:,:) => null() , &
							activeVortInput(:) => null(), &
							activeVortStage1(:) => null(), &
							activeVortStage2(:) => null(), &
							activeVortStage3(:) => null(), &
							activeVortStage4(:) => null(), &
							newActiveVort(:) => null(),&
							area(:) => null(), &
							newActivePanelsX(:,:) => null() , &
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
character(len=28), save :: logKey = 'BVE'
integer(kint), parameter :: logLevel = TRACE_LOGGING_LEVEL
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
subroutine NewPrivate(self,aMesh,nProcs)
! Creates memory for an RK4 object.
!	Calls LoadBalance and ZeroRK to return an object ready for use.
	type(BVERK4Data), intent(out) :: self
	type(SphereMesh), intent(in) :: aMesh
	integer(kint), intent(in) :: nProcs
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: nParticles, nActive, nPassive, nTracer, panelKind, problemKind

	if ( .NOT. loginit) call InitLogger(log,procRank)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Creating New BVERK4Data.')

	self%rk4isReady = .FALSE.
	self%rk4isReady = .FALSE.

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	nParticles = aParticles%N
	nActive = aPanels%N_Active
	nPassive = aPanels%N - nActive

	nTracer = aMesh%nTracer

	panelKind = aMesh%panelKind

	problemKind = BVE_SOLVER

	!
	!	Allocate MPI variables
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
	!	Set physical parameters
	!
	self%smooth = VELOCITY_SMOOTH

	!
	!	Separate active panels from passive panels
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

	call GatherPanels(aPanels,self%activePanels,self%activeMap,self%passivePanels,self%passiveMap)

	!
	!	Allocate RK4 variables
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

	allocate(self%activeVortInput(nActive))
	allocate(self%activeVortStage1(nActive))
	allocate(self%activeVortStage2(nActive))
	allocate(self%activeVortStage3(nActive))
	allocate(self%activeVortStage4(nActive))
	allocate(self%newActiveVort(nActive))

	allocate(self%area(nActive))

	allocate(self%passivePanelsInput(3,nPassive))
	allocate(self%passivePanelsStage1(3,nPassive))
	allocate(self%passivePanelsStage2(3,nPassive))
	allocate(self%passivePanelsStage3(3,nPassive))
	allocate(self%passivePanelsStage4(3,nPassive))
	allocate(self%newPassivePanelsX(3,nPassive))

	self%rk4isReady = .TRUE.
	call ZeroRK4(self)

	call LoadBalance(self,nParticles,nActive,nPassive,nProcs)

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,' RK4 data ready.')
end subroutine

subroutine DeletePrivate(self)
!	Deletes / frees memory associated with an RK4 object
!
	type(BVERK4Data), intent(inout) :: self

	self%mpiIsReady = .FALSE.
	self%rk4isReady = .FALSE.

	deallocate(self%activeMap)
	deallocate(self%passiveMap)
	call Delete(self%activePanels)
	call Delete(self%passivePanels)
	deallocate(self%activePanels)
	deallocate(self%passivePanels)

	! MPI variables
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

	deallocate(self%activeVortInput)
	deallocate(self%activeVortStage1)
	deallocate(self%activeVortStage2)
	deallocate(self%activeVortStage3)
	deallocate(self%activeVortStage4)
	deallocate(self%newActiveVort)

	deallocate(self%passivePanelsInput)
	deallocate(self%passivePanelsStage1)
	deallocate(self%passivePanelsStage2)
	deallocate(self%passivePanelsStage3)
	deallocate(self%passivePanelsStage4)
	deallocate(self%newPassivePanelsX)

end subroutine

!
!----------------
! Public member functions
!----------------
!
subroutine BVERK4Timestep(self,aMesh,dt,procRank, nProcs)
!	Advances aMesh one time step forward using 4th order Runge-Kutta
!
!	Biot-Savart law is computed in parallel using replicated data on each processor.
!	Each process computes its portion of the integral, then broadcasts its results to MPI_WORLD_COMM.
	type(BVERK4Data), intent(inout) :: self
	type(SphereMesh), intent(inout) :: aMesh
	real(kreal), intent(in) :: dt
	integer(kint), intent(in) :: procRank, nProcs
	type(Particles), pointer :: aParticles
	integer(kint) :: j, mpiErrCode

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering RK4.')

	if ( self%rk4IsReady .AND. self%mpiIsReady ) then
		call ZeroRK4(self)
	else
		call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'BVERK4 ERROR : data not ready.')
		return
	endif

	aParticles => aMesh%particles

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 1       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 1
	self%particlesInput = aParticles%x(:,1:aParticles%N)
	self%activePanelsInput = self%activePanels%x
	self%activeVortInput = self%activePanels%relVort
	self%area = self%activePanels%area
	self%passivePanelsInput = self%passivePanels%x

	call BVEPassiveRHS(self%particlesStage1, self%particlesInput,&
					   self%activePanelsInput, self%activeVortInput,self%area, &
					   self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))
	call BVESmoothVelocity(self%passivePanelsStage1, self%passivePanelsInput, &
						   self%activePanelsInput, self%activeVortInput, self%area, &
						   self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank))
	call BVEActiveRHS(self%activePanelsStage1, self%activeVortStage1, &
					 self%activePanelsInput, self%activeVortInput, self%area, &
					 self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank))

	do j=0,nProcs-1
		call MPI_BCAST(self%particlesStage1(:,self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
					   3*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%passivePanelsStage1(:,self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
					   3*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%activePanelsStage1(:,self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						3*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%activeVortStage1(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
	enddo

	! STAGE 1 ONLY : Store velocity and kinetic energy
	do j=1,aParticles%N
		aParticles%ke(j) = sum(self%particlesStage1(:,j)*self%particlesStage1(:,j))
		aParticles%u(:,j) = self%particlesStage1(:,j)
	enddo
	do j=1,self%activePanels%N
		self%activePanels%ke(j) = sum(self%activePanelsStage1(:,j)*self%activePanelsStage1(:,j))
		self%activePanels%u(:,j) = self%activePanelsStage1(:,j)
	enddo
	aMesh%totalKE = 0.5*sum(self%activePanels%ke*self%area)


	self%particlesStage1 = dt*self%particlesStage1
	self%passivePanelsStage1 = dt*self%passivePanelsStage1
	self%activePanelsStage1 = dt*self%activePanelsStage1
	self%activeVortStage1 = dt*self%activeVortStage1

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 2       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 2
	self%particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*self%particlesStage1
	self%activePanelsInput = self%activePanels%x + 0.5_kreal*self%activePanelsStage1
	self%activeVortInput = self%activePanels%relVort + 0.5_kreal*self%activeVortStage1
	self%area = self%activePanels%area
	self%passivePanelsInput = self%passivePanels%x + 0.5_kreal*self%passivePanelsStage1

	call BVEPassiveRHS(self%particlesStage2, self%particlesInput,&
					   self%activePanelsInput, self%activeVortInput,self%area, &
					   self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))
	call BVESmoothVelocity(self%passivePanelsStage2, self%passivePanelsInput, &
						   self%activePanelsInput, self%activeVortInput, self%area, &
						   self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank))
	call BVEActiveRHS(self%activePanelsStage2, self%activeVortStage2, &
					 self%activePanelsInput, self%activeVortInput, self%area, &
					 self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank))

	do j=0,nProcs-1
		call MPI_BCAST(self%particlesStage2(:,self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
					   3*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%passivePanelsStage2(:,self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
					   3*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%activePanelsStage2(:,self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						3*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%activeVortStage2(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
	enddo

	self%particlesStage2 = dt*self%particlesStage2
	self%passivePanelsStage2 = dt*self%passivePanelsStage2
	self%activePanelsStage2 = dt*self%activePanelsStage2
	self%activeVortStage2 = dt*self%activeVortStage2

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 3       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 3
	self%particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*self%particlesStage2
	self%activePanelsInput = self%activePanels%x + 0.5_kreal*self%activePanelsStage2
	self%activeVortInput = self%activePanels%relVort + 0.5_kreal*self%activeVortStage2
	self%area = self%activePanels%area
	self%passivePanelsInput = self%passivePanels%x + 0.5_kreal*self%passivePanelsStage2

	call BVEPassiveRHS(self%particlesStage3, self%particlesInput,&
					   self%activePanelsInput, self%activeVortInput,self%area, &
					   self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))
	call BVESmoothVelocity(self%passivePanelsStage3, self%passivePanelsInput, &
						   self%activePanelsInput, self%activeVortInput, self%area, &
						   self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank))
	call BVEActiveRHS(self%activePanelsStage3, self%activeVortStage3, &
					 self%activePanelsInput, self%activeVortInput, self%area, &
					 self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank))

	do j=0,nProcs-1
		call MPI_BCAST(self%particlesStage3(:,self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
					   3*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%passivePanelsStage3(:,self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
					   3*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%activePanelsStage3(:,self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						3*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%activeVortStage3(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
	enddo

	self%particlesStage3 = dt*self%particlesStage3
	self%passivePanelsStage3 = dt*self%passivePanelsStage3
	self%activePanelsStage3 = dt*self%activePanelsStage3
	self%activeVortStage3 = dt*self%activeVortStage3


	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 4       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 4
	self%particlesInput = aParticles%x(:,1:aParticles%N) + self%particlesStage3
	self%activePanelsInput = self%activePanels%x + self%activePanelsStage3
	self%activeVortInput = self%activePanels%relVort + self%activeVortStage3
	self%area = self%activePanels%area
	self%passivePanelsInput = self%passivePanels%x + self%passivePanelsStage3

	call BVEPassiveRHS(self%particlesStage4, self%particlesInput,&
					   self%activePanelsInput, self%activeVortInput,self%area, &
					   self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))
	call BVESmoothVelocity(self%passivePanelsStage4, self%passivePanelsInput, &
						   self%activePanelsInput, self%activeVortInput, self%area, &
						   self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank))
	call BVEActiveRHS(self%activePanelsStage4, self%activeVortStage4, &
					 self%activePanelsInput, self%activeVortInput, self%area, &
					 self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank))

	do j=0,nProcs-1
		call MPI_BCAST(self%particlesStage4(:,self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
					   3*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%passivePanelsStage4(:,self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
					   3*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%activePanelsStage4(:,self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						3*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST(self%activeVortStage4(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
	enddo

	self%particlesStage4 = dt*self%particlesStage4
	self%passivePanelsStage4 = dt*self%passivePanelsStage4
	self%activePanelsStage4 = dt*self%activePanelsStage4
	self%activeVortStage4 = dt*self%activeVortStage4

	!!!!!!!!!!!!!!!!!
	! 	RK update   !
	!!!!!!!!!!!!!!!!!

	self%newParticlesX = aParticles%x(:,1:aParticles%N) + self%particlesStage1/6.0_kreal + &
			self%particlesStage2/3.0_kreal + self%particlesStage3/3.0_kreal + self%particlesStage4/6.0_kreal
	self%newActivePanelsX = self%activePanels%x + self%activePanelsStage1/6.0_kreal + &
			self%activePanelsStage2/3.0_kreal + self%activePanelsStage3/3.0_kreal + self%activePanelsStage4/6.0_kreal
	self%newPassivePanelsX = self%passivePanels%x + self%passivePanelsStage1/6.0_kreal + &
			self%passivePanelsStage2/3.0_kreal + self%passivePanelsStage3/3.0_kreal + self%passivePanelsStage4/6.0_kreal
	self%newActiveVort = self%activePanels%relVort + self%activeVortStage1/6.0_kreal + &
			self%activeVortStage2/3.0_kreal + self%activeVortStage3/3.0_kreal + self%activeVortStage4/6.0_kreal

	aParticles%x(:,1:aParticles%N) = self%newParticlesX
	aParticles%relVort(1:aParticles%N) = aParticles%absVort(1:aParticles%N) - &
		2.0_kreal*OMEGA/EARTH_RADIUS*self%newParticlesX(3,1:aParticles%N)
	self%activePanels%x = self%newActivePanelsX
	self%activePanels%relVort = self%newActiveVort
	self%passivePanels%x = self%newPassivePanelsX

	call ScatterPanels(aMesh%panels,self%activePanels,self%activeMap,self%passivePanels,self%passiveMap)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'... timestep complete.')
end subroutine

!
!----------------
! Module methods : module- or type-specific private functions
!----------------
!
subroutine ZeroRK4(self)
!	Sets RK4 variables to zero
!
	type(BVERK4Data), intent(inout) :: self
	if ( .not. self%rk4isReady ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'ZeroRK4 ERROR : ','memory not ready.')
		return
	endif

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering ZeroRK4.')

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

	self%activeVortInput = 0.0_kreal
	self%activeVortStage1 = 0.0_kreal
	self%activeVortStage2 = 0.0_kreal
	self%activeVortStage3 = 0.0_kreal
	self%activeVortStage4 = 0.0_kreal
	self%newActiveVort = 0.0_kreal

	self%area = 0.0_kreal

	self%passivePanelsInput = 0.0_kreal
	self%passivePanelsStage1 = 0.0_kreal
	self%passivePanelsStage2 = 0.0_kreal
	self%passivePanelsStage3 = 0.0_kreal
	self%passivePanelsStage4 = 0.0_kreal
	self%newPassivePanelsX = 0.0_kreal

end subroutine

subroutine LoadBalance(self,nParticles,nActive,nPassive,nProcs)
!	Calculates the portion of each work array to be done by each process.
!
	type(BVERK4Data), intent(inout) :: self
	integer(kint), intent(in) :: nParticles,nActive, nPassive, nProcs
	integer(kint) :: j, chunkSize

	call LogMessage(log,DEBUG_LOGGING_LEVEL,' Entering LoadBalance.')

	chunkSize = nParticles / nProcs
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

subroutine BVEActiveRHS(dX, dVort, pointVorts, vort, area, indexStart, indexEnd)
!	Computes the Biot-Savart law for active particles.
!	Skips singular panel.
!
	real(kreal), intent(out) :: dX(:,:)
	real(kreal), intent(out) :: dVort(:)
	real(kreal), intent(in) :: pointVorts(:,:)
	real(kreal), intent(in) :: vort(:)
	real(kreal), intent(in) :: area(:)
	integer(kint), intent(in) :: indexStart, indexEnd
	integer(kint) :: nn, j, k
	real(kreal) :: denom

	! Error checking
	nn = size(dX,2)
	if ( nn /= size(dVort) .OR. nn /= size(pointVorts,2) ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'BVEActiveRHS ERROR : ','size mismatch 1.')
		return
	endif
	if ( nn /= size(vort) .OR. nn /= size(area) ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'BVEActiveRHS ERROR : ','size mismatch 2.')
		return
	endif

	dX = 0.0_kreal
	dVort = 0.0_kreal
	do j=indexStart, indexEnd
		do k=1,j-1
				denom = EARTH_RADIUS*EARTH_RADIUS - sum(pointVorts(:,k)*pointVorts(:,j))
				dX(1,j) = dX(1,j) + (pointVorts(2,j)*pointVorts(3,k) - pointVorts(3,j)*pointVorts(2,k))*vort(k)*area(k)/denom
				dX(2,j) = dX(2,j) + (pointVorts(3,j)*pointVorts(1,k) - pointVorts(1,j)*pointVorts(3,k))*vort(k)*area(k)/denom
				dX(3,j) = dX(3,j) + (pointVorts(1,j)*pointVorts(2,k) - pointVorts(2,j)*pointVorts(1,k))*vort(k)*area(k)/denom
		enddo
		do k=j+1,nn
			denom = EARTH_RADIUS*EARTH_RADIUS - sum(pointVorts(:,k)*pointVorts(:,j))
			dX(1,j) = dX(1,j) + (pointVorts(2,j)*pointVorts(3,k) - pointVorts(3,j)*pointVorts(2,k))*vort(k)*area(k)/denom
			dX(2,j) = dX(2,j) + (pointVorts(3,j)*pointVorts(1,k) - pointVorts(1,j)*pointVorts(3,k))*vort(k)*area(k)/denom
			dX(3,j) = dX(3,j) + (pointVorts(1,j)*pointVorts(2,k) - pointVorts(2,j)*pointVorts(1,k))*vort(k)*area(k)/denom
		enddo
	enddo
	dX = dX/(-4.0_kreal*PI*EARTH_RADIUS)
	dVort = -2.0_kreal*Omega*dX(3,:)/EARTH_RADIUS
end subroutine


subroutine BVEPassiveRHS(dX, X,pointVorts,vort,area,indexStart,indexEnd)
!  Performs the Biot-Savart integral summation in the rotating frame using midpoint rule
!  quadrature for a set of passive particles (panel vertices).
	real(kreal), intent(out) :: dX(:,:)
!	real(kreal), intent(out) :: dVort(:)
	real(kreal), intent(in) :: X(:,:)
	real(kreal), intent(in) :: pointVorts(:,:)
	real(kreal), intent(in) :: vort(:)
	real(kreal), intent(in) :: area(:)
	integer(kint), intent(in) :: indexStart, indexEnd
	integer(kint) :: mm, nn, j, k
	real(kreal) :: denom
	! Check for size mismatch errors
	mm = size(x,2)
	nn = size(pointVorts,2)
	if ( mm /= size(dX,2) ) then
		print *,"BVEPassiveRHS ERROR : size mismatch 1."
		return
	endif
	if ( nn /= size(vort) ) then
		print *,"BVEPassiveRHS ERROR : size mismatch 2."
		return
	endif
	if ( nn /= size(area) ) then
		print *,"BVEPassiveRHS ERROR : size mismatch 3."
		return
	endif
	dX = 0.0_kreal
!	dVort = 0.0_kreal
	do j=indexStart,indexEnd
		do k=1,nn
			denom = EARTH_RADIUS*EARTH_RADIUS - sum(pointVorts(:,k)*X(:,j))
			dX(1,j) = dX(1,j) + (X(2,j)*pointVorts(3,k) - X(3,j)*pointVorts(2,k))*vort(k)*area(k)/denom
			dX(2,j) = dX(2,j) + (X(3,j)*pointVorts(1,k) - X(1,j)*pointVorts(3,k))*vort(k)*area(k)/denom
			dX(3,j) = dX(3,j) + (X(1,j)*pointVorts(2,k) - X(2,j)*pointVorts(1,k))*vort(k)*area(k)/denom
		enddo
	enddo
	dX = dX/(-4.0_kreal*PI*EARTH_RADIUS)
!	dVort = -2.0_kreal*OMEGA*dX(3,:)
end subroutine


subroutine BVESmoothVelocity(dX,X,pointVorts,vort,Area,indexStart,indexEnd)
!  Performs the regularized Biot-Savart integral summation in the rotating frame using midpoint rule
!  quadrature for a given smoothing parameter.
	real(kreal), intent(out) :: dX(:,:)
	real(kreal), intent(in) :: X(:,:)
	real(kreal), intent(in) :: pointVorts(:,:)
	real(kreal), intent(in) :: vort(:)
	real(kreal), intent(in) :: area(:)
	integer(kint), intent(in) :: indexStart, indexEnd
	integer(kint) :: j,k, nn, mm
	real(kreal) :: denom
	! Check for size mismatch errors
	mm = size(x,2)
	if ( mm /= size(dx,2) ) then
		print *,"BVEInertialVelocity ERROR : size mismatch 1."
		return
	endif
	nn = size(pointVorts,2)
	if ( nn /= size(vort) ) then
		print *,"BVEInertialVelocity ERROR : size mismatch 2."
		return
	endif
	if ( nn /= size(area) ) then
		print *,"BVEInertialVelocity ERROR : size mismatch 3."
		return
	endif
	dX = 0.0_kreal
	do j=indexStart,indexEnd
		do k=1,nn
			denom = EARTH_RADIUS*EARTH_RADIUS - sum(pointVorts(:,k)*x(:,j)) + VELOCITY_SMOOTH*VELOCITY_SMOOTH
			dX(1,j) = dX(1,j) + (X(2,j)*pointVorts(3,k) - X(3,j)*pointVorts(2,k))*vort(k)*area(k)/denom
			dX(2,j) = dX(2,j) + (X(3,j)*pointVorts(1,k) - X(1,j)*pointVorts(3,k))*vort(k)*area(k)/denom
			dX(3,j) = dX(3,j) + (X(1,j)*pointVorts(2,k) - X(2,j)*pointVorts(1,k))*vort(k)*area(k)/denom
		enddo
	enddo
	dX = dX/(-4.0_kreal*PI*EARTH_RADIUS)
end subroutine


subroutine InitLogger(aLog,rank)
	type(Logger), intent(inout) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

end module
