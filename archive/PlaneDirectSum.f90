module PlaneDirectSumModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method
!------------------------------------------------------------------------------
!
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
!> @defgroup PlaneDirectSum Planar Direct Sum RK4
!> Provides RK4 time stepping for the 2d incompressible equations using direct summation.
!
!
! DESCRIPTION:
!> @file
!> Provides RK4 time stepping for the 2d incompressible equations using direct summation.
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

include 'mpif.h'

private
public PlaneRK4DirectSum
public New, Delete
public RK4TimestepNoRotation, RK4TimestepIncompressibleAdvection
public VELOCITY_SMOOTH
public LeVeque93

!
!----------------
! Types and module constants
!----------------
!
real(kreal), protected, save :: VELOCITY_SMOOTH = 0.01_kreal

type PlaneRK4DirectSum
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
character(len=28), save :: logKey = 'PlaneRK4'
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

interface
	function AdvectionVelocity(xy, t)
		real(8) :: AdvectionVelocity(2)
		real(8), intent(in) :: xy(2), t
	end function
end interface


contains
!
!----------------
! Standard methods : Constructor / Destructor
!----------------
!
subroutine NewPrivate(self, aMesh, nProcs)
	type(PlaneRK4DirectSum), intent(out) :: self
	type(PlaneMesh), intent(in) :: aMesh
	integer(kint), intent(in) :: nProcs
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: nParticles, nActive, nPassive, nTracer, panelKind, problemKind

	if ( .NOT. loginit) call InitLogger(log, procRank)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Creating new PlaneRK4DirectSum.')

	self%rk4isReady = .FALSE.

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	nParticles = aParticles%N
	nActive = aPanels%N_Active
	nPassive = aPanels%N - aPanels%N_Active

	nTracer = aMesh%nTracer

	panelkind = QUAD_PANEL
	problemKind = PLANE_SOLVER

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
	! separate active panels from passive panels
	!
	allocate(self%activePanels)
	call New(self%activePanels, nActive, panelKind, nTracer, problemKind)
	allocate(self%activeMap(nActive))
	self%activePanels%N = nActive
	self%activePanels%N_active = nActive
	self%activeMap = 0

	allocate(self%passivePanels)
	call New(self%passivePanels, nPassive, panelKind, nTracer, problemKind)
	allocate(self%passiveMap(nPassive))
	self%passivePanels%N = nPassive
	self%passiveMap = 0

	call GatherPanels(aPanels, self%activePanels, self%activeMap, self%passivePanels, self%passiveMap)

	!
	!	Allocate RK4 variables
	!
	allocate(self%particlesInput(2,nParticles))
	allocate(self%particlesStage1(2,nParticles))
	allocate(self%particlesStage2(2,nParticles))
	allocate(self%particlesStage3(2,nParticles))
	allocate(self%particlesStage4(2,nParticles))
	allocate(self%newParticlesX(2,nParticles))

	allocate(self%activePanelsInput(2,nActive))
	allocate(self%activePanelsStage1(2,nActive))
	allocate(self%activePanelsStage2(2,nActive))
	allocate(self%activePanelsStage3(2,nActive))
	allocate(self%activePanelsStage4(2,nActive))
	allocate(self%newActivePanelsX(2,nActive))

	allocate(self%activeVortInput(nActive))
	allocate(self%activeVortStage1(nActive))
	allocate(self%activeVortStage2(nActive))
	allocate(self%activeVortStage3(nActive))
	allocate(self%activeVortStage4(nActive))
	allocate(self%newActiveVort(nActive))

	allocate(self%area(nActive))

	allocate(self%passivePanelsInput(2,nPassive))
	allocate(self%passivePanelsStage1(2,nPassive))
	allocate(self%passivePanelsStage2(2,nPassive))
	allocate(self%passivePanelsStage3(2,nPassive))
	allocate(self%passivePanelsStage4(2,nPassive))
	allocate(self%newPassivePanelsX(2,nPassive))

	self%rk4isReady = .TRUE.
	call ZeroRK4(self)

	call LoadBalance(self,nParticles,nActive,nPassive,nProcs)

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,' RK4 data ready.')
end subroutine


subroutine DeletePrivate(self)
!	Deletes / frees memory associated with an RK4 object
!
	type(PlaneRK4DirectSum), intent(inout) :: self

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
subroutine RK4TimestepNoRotation(self, aMesh, dt, procRank, nProcs)
	type(PlaneRK4DirectSum), intent(inout) :: self
	type(PlaneMesh), intent(inout) :: aMesh
	real(kreal), intent(in) :: dt
	integer(kint), intent(in) :: procRank, nProcs
	!
	type(Particles), pointer :: aParticles
	integer(kint) :: j, mpiErrCode

	call LogMessage(log, DEBUG_LOGGING_LEVEL, logkey,' entering rk4 timestep.')

	if ( self%rk4isReady .AND. self%mpiIsReady ) then
		call ZeroRK4(self)
	else
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, ' RK4Timestep ERROR : data not ready.')
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

	call PlanePassiveVel_noRotation( self%particlesStage1, self%particlesInput, &
		 	self%activePanelsInput, self%activeVortInput, self%area, &
		 	self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))
	call PlaneActiveVel_noRotation( self%activePanelsStage1, self%activePanelsInput, self%activeVortInput, self%area, &
		 	self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank))
	call PlaneSmoothVel_noRotation( self%passivePanelsStage1, self%passivePanelsINput, &
			self%activePanelsINput, self%activeVortInput, self%area, &
			self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank))
	do j = 0, nProcs - 1
		call MPI_BCAST( self%particlesStage1(:, self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
					    2*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					    MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( self%activePanelsStage1(:, self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						2*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( self%passivePanelsStage1(:, self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
						2*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
	enddo

	!
	!	STAGE 1 ONLY : store velocity and kinetic energy
	!
	do j = 1, aParticles%N
		aParticles%ke(j) = sum( self%particlesStage1(:,j) * self%particlesStage1(:,j))
		aParticles%u(:,j) = self%particlesStage1(:,j)
	enddo
	do j = 1, self%activePanels%N
		self%activePanels%ke(j) = sum( self%activePanelsStage1(:,j) * self%activePanelsStage1(:,j) )
		self%activePanels%u(:,j) = self%activePanelsStage1(:,j)
	enddo

	self%particlesStage1 = dt * self%particlesStage1
	self%activePanelsStage1 = dt * self%activePanelsStage1
	self%passivePanelsStage1 = dt * self%passivePanelsStage1


	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 2       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 2
	self%particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*self%particlesStage1
	self%activePanelsInput = self%activePanels%x + 0.5_kreal*self%activePanelsStage1
	self%passivePanelsInput = self%passivePanels%x + 0.5_kreal*self%passivePanelsStage1

	call PlanePassiveVel_noRotation( self%particlesStage2, self%particlesInput, &
		 	self%activePanelsInput, self%activeVortInput, self%area, &
		 	self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))
	call PlaneActiveVel_noRotation( self%activePanelsStage2, self%activePanelsInput, self%activeVortInput, self%area, &
		 	self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank))
	call PlaneSmoothVel_noRotation( self%passivePanelsStage2, self%passivePanelsINput, &
			self%activePanelsINput, self%activeVortInput, self%area, &
			self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank))
	do j = 0, nProcs - 1
		call MPI_BCAST( self%particlesStage2(:, self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
					    2*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					    MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( self%activePanelsStage2(:, self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						2*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( self%passivePanelsStage2(:, self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
						2*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
	enddo
	self%particlesStage2 = dt * self%particlesStage2
	self%activePanelsStage2 = dt * self%activePanelsStage2
	self%passivePanelsStage2 = dt * self%passivePanelsStage2

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 3       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 3
	self%particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*self%particlesStage2
	self%activePanelsInput = self%activePanels%x + 0.5_kreal*self%activePanelsStage2
	self%passivePanelsInput = self%passivePanels%x + 0.5_kreal*self%passivePanelsStage2
	call PlanePassiveVel_noRotation( self%particlesStage3, self%particlesInput, &
		 	self%activePanelsInput, self%activeVortInput, self%area, &
		 	self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))
	call PlaneActiveVel_noRotation( self%activePanelsStage3, self%activePanelsInput, self%activeVortInput, self%area, &
		 	self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank))
	call PlaneSmoothVel_noRotation( self%passivePanelsStage3, self%passivePanelsINput, &
			self%activePanelsINput, self%activeVortInput, self%area, &
			self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank))
	do j = 0, nProcs - 1
		call MPI_BCAST( self%particlesStage3(:, self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
					    2*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					    MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( self%activePanelsStage3(:, self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						2*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( self%passivePanelsStage3(:, self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
						2*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
	enddo
	self%particlesStage3 = dt * self%particlesStage3
	self%activePanelsStage3 = dt * self%activePanelsStage3
	self%passivePanelsStage3 = dt * self%passivePanelsStage3

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 4       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 4
	self%particlesInput = aParticles%x(:,1:aParticles%N) + self%particlesStage3
	self%activePanelsInput = self%activePanels%x + self%activePanelsStage3
	self%passivePanelsInput = self%passivePanels%x + self%passivePanelsStage3
	call PlanePassiveVel_noRotation( self%particlesStage4, self%particlesInput, &
		 	self%activePanelsInput, self%activeVortInput, self%area, &
		 	self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))
	call PlaneActiveVel_noRotation( self%activePanelsStage4, self%activePanelsInput, self%activeVortInput, self%area, &
		 	self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank))
	call PlaneSmoothVel_noRotation( self%passivePanelsStage4, self%passivePanelsINput, &
			self%activePanelsINput, self%activeVortInput, self%area, &
			self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank))
	do j = 0, nProcs - 1
		call MPI_BCAST( self%particlesStage4(:, self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
					    2*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					    MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( self%activePanelsStage4(:, self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						2*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
		call MPI_BCAST( self%passivePanelsStage4(:, self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
						2*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
						MPI_COMM_WORLD, mpiErrCode)
	enddo
	self%particlesStage4 = dt * self%particlesStage4
	self%activePanelsStage4 = dt * self%activePanelsStage4
	self%passivePanelsStage4 = dt * self%passivePanelsStage4

	!!!!!!!!!!!!!!!!!
	! 	RK update   !
	!!!!!!!!!!!!!!!!!

	self%newParticlesX = aParticles%x(:,1:aParticles%N) + self%particlesStage1/6.0_kreal + &
			self%particlesStage2/3.0_kreal + self%particlesStage3/3.0_kreal + self%particlesStage4/6.0_kreal
	self%newActivePanelsX = self%activePanels%x + self%activePanelsStage1/6.0_kreal + &
			self%activePanelsStage2/3.0_kreal + self%activePanelsStage3/3.0_kreal + self%activePanelsStage4/6.0_kreal
	self%newPassivePanelsX = self%passivePanels%x + self%passivePanelsStage1/6.0_kreal + &
			self%passivePanelsStage2/3.0_kreal + self%passivePanelsStage3/3.0_kreal + self%passivePanelsStage4/6.0_kreal

	aParticles%x(:,1:aParticles%N) = self%newParticlesX
	self%activePanels%x = self%newActivePanelsX
	self%passivePanels%x = self%newPassivePanelsX

	call ScatterPanels( aMesh%panels, self%activePanels, self%activeMap, self%passivePanels, self%passiveMap)

	call LogMessage(log, DEBUG_LOGGING_LEVEL, logkey, '... timestep complete.')
end subroutine

subroutine RK4TimestepIncompressibleAdvection(self, aMesh, dt, t, procRank, nProcs, velocityFunction)
	type(PlaneRK4DirectSum), intent(inout) :: self
	type(PlaneMesh), intent(inout) :: aMesh
	real(kreal), intent(in) :: dt, t
	integer(kint), intent(in) :: procRank, nProcs
	procedure(AdvectionVelocity) :: velocityFunction
	!
	type(Particles), pointer :: aParticles
	integer(kint) :: j, errCode

	call LogMessage(log, DEBUG_LOGGING_LEVEL, logkey, 'entering advection rk4 timestep.')

	if ( self%rk4isReady .AND. self%mpiIsReady ) then
		call ZeroRK4(self)
	else
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, ' RK4Timestep ERROR : data not ready.')
		return
	endif

	aParticles => aMesh%particles

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 1       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 1
	self%particlesInput = aParticles%x(:,1:aParticles%N)
	self%activePanelsInput = self%activePanels%x
	self%passivePanelsInput = self%passivePanels%x

	do j = self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank)
		self%particlesStage1(:,j) = VelocityFunction(self%particlesInput(:,j), t)
	enddo
	do j = self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank)
		self%activePanelsStage1(:,j) = VelocityFunction(self%activePanelsInput(:,j), t)
	enddo
	do j = self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank)
		self%passivePanelsStage1(:,j) = VelocityFunction(self%passivePanelsInput(:,j), t)
	enddo

	do j = 0, nProcs - 1
		call MPI_BCAST(self%particlesStage1(:,self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
				2*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsStage1(:,self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
				2*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%passivePanelsStage1(:,self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
				2*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
	enddo

	!
	!	STAGE 1 ONLY : store velocity and kinetic energy
	!
	do j = 1, aParticles%N
		aParticles%ke(j) = sum( self%particlesStage1(:,j) * self%particlesStage1(:,j))
		aParticles%u(:,j) = self%particlesStage1(:,j)
	enddo
	do j = 1, self%activePanels%N
		self%activePanels%ke(j) = sum( self%activePanelsStage1(:,j) * self%activePanelsStage1(:,j) )
		self%activePanels%u(:,j) = self%activePanelsStage1(:,j)
	enddo

	self%particlesStage1 = dt * self%particlesStage1
	self%activePanelsStage1 = dt * self%activePanelsStage1
	self%passivePanelsStage1 = dt * self%passivePanelsStage1


	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 2       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 2
	self%particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*self%particlesStage1
	self%activePanelsInput = self%activePanels%x + 0.5_kreal*self%activePanelsStage1
	self%passivePanelsInput = self%passivePanels%x + 0.5_kreal*self%passivePanelsStage1

	do j = self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank)
		self%particlesStage2(:,j) = VelocityFunction(self%particlesInput(:,j), t)
	enddo
	do j = self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank)
		self%activePanelsStage2(:,j) = VelocityFunction(self%activePanelsInput(:,j), t)
	enddo
	do j = self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank)
		self%passivePanelsStage2(:,j) = VelocityFunction(self%passivePanelsInput(:,j), t)
	enddo

	do j = 0, nProcs - 1
		call MPI_BCAST(self%particlesStage2(:,self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
				2*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsStage2(:,self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
				2*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%passivePanelsStage2(:,self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
				2*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
	enddo

	self%particlesStage2 = dt * self%particlesStage2
	self%activePanelsStage2 = dt * self%activePanelsStage2
	self%passivePanelsStage2 = dt * self%passivePanelsStage2

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 3       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 3
	self%particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*self%particlesStage2
	self%activePanelsInput = self%activePanels%x + 0.5_kreal*self%activePanelsStage2
	self%passivePanelsInput = self%passivePanels%x + 0.5_kreal*self%passivePanelsStage2

	do j = self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank)
		self%particlesStage3(:,j) = VelocityFunction(self%particlesInput(:,j), t)
	enddo
	do j = self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank)
		self%activePanelsStage3(:,j) = VelocityFunction(self%activePanelsInput(:,j), t)
	enddo
	do j = self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank)
		self%passivePanelsStage3(:,j) = VelocityFunction(self%passivePanelsInput(:,j), t)
	enddo

	do j = 0, nProcs - 1
		call MPI_BCAST(self%particlesStage3(:,self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
				2*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsStage3(:,self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
				2*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%passivePanelsStage3(:,self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
				2*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
	enddo

	self%particlesStage3 = dt * self%particlesStage3
	self%activePanelsStage3 = dt * self%activePanelsStage3
	self%passivePanelsStage3 = dt * self%passivePanelsStage3

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 4       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 4
	self%particlesInput = aParticles%x(:,1:aParticles%N) + self%particlesStage3
	self%activePanelsInput = self%activePanels%x + self%activePanelsStage3
	self%passivePanelsInput = self%passivePanels%x + self%passivePanelsStage3

	do j = self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank)
		self%particlesStage4(:,j) = VelocityFunction(self%particlesInput(:,j), t)
	enddo
	do j = self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank)
		self%activePanelsStage4(:,j) = VelocityFunction(self%activePanelsInput(:,j), t)
	enddo
	do j = self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank)
		self%passivePanelsStage4(:,j) = VelocityFunction(self%passivePanelsInput(:,j), t)
	enddo

	do j = 0, nProcs - 1
		call MPI_BCAST(self%particlesStage4(:,self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
				2*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsStage4(:,self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
				2*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%passivePanelsStage4(:,self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
				2*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
	enddo
	self%particlesStage4 = dt * self%particlesStage4
	self%activePanelsStage4 = dt * self%activePanelsStage4
	self%passivePanelsStage4 = dt * self%passivePanelsStage4
	!!!!!!!!!!!!!!!!!
	! 	RK update   !
	!!!!!!!!!!!!!!!!!

	self%newParticlesX = aParticles%x(:,1:aParticles%N) + self%particlesStage1/6.0_kreal + &
			self%particlesStage2/3.0_kreal + self%particlesStage3/3.0_kreal + self%particlesStage4/6.0_kreal
	self%newActivePanelsX = self%activePanels%x + self%activePanelsStage1/6.0_kreal + &
			self%activePanelsStage2/3.0_kreal + self%activePanelsStage3/3.0_kreal + self%activePanelsStage4/6.0_kreal
	self%newPassivePanelsX = self%passivePanels%x + self%passivePanelsStage1/6.0_kreal + &
			self%passivePanelsStage2/3.0_kreal + self%passivePanelsStage3/3.0_kreal + self%passivePanelsStage4/6.0_kreal

	aParticles%x(:,1:aParticles%N) = self%newParticlesX
	self%activePanels%x = self%newActivePanelsX
	self%passivePanels%x = self%newPassivePanelsX

	call ScatterPanels( aMesh%panels, self%activePanels, self%activeMap, self%passivePanels, self%passiveMap)

	call LogMessage(log, DEBUG_LOGGING_LEVEL, logkey, '... timestep complete.')
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
	type(PlaneRK4DirectSum), intent(inout) :: self
	!
	if ( .not. self%rk4isReady ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'ZeroRK4 ERROR : ','memory not ready.')
		return
	endif

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
	type(PlaneRK4DirectSum), intent(inout) :: self
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

subroutine PlaneActiveVel_noRotation(dX, x, vort, area, indexStart, indexEnd)
	real(kreal), intent(out) :: dX(:,:)
	real(kreal), intent(in) :: x(:,:), vort(:), area(:)
	integer(kint), intent(in) :: indexStart, indexEnd
	!
	integer(kint) :: j, k, n
	real(kreal) :: denom

	n = size(vort)

	dX = 0.0_kreal
	do j = indexStart, indexEnd
		do k = 1, j-1
			denom = 2.0_kreal * PI * ( sum( (x(:,j) - x(:,k)) * (x(:,j) - x(:,k))))
			dX(1,j) = dx(1,j) - ( x(2,j) - x(2,k) ) * vort(k) * area(k) / denom
			dX(2,j) = dx(2,j) + ( x(1,j) - x(1,k) ) * vort(k) * area(k) / denom
		enddo
		do k = j+1, n
			denom = 2.0_kreal * PI * ( sum( (x(:,j) - x(:,k)) * (x(:,j) - x(:,k))))
			dX(1,j) = dx(1,j) - ( x(2,j) - x(2,k) ) * vort(k) * area(k) / denom
			dX(2,j) = dx(2,j) + ( x(1,j) - x(1,k) ) * vort(k) * area(k) / denom
		enddo
	enddo
end subroutine

subroutine PlanePassiveVel_noRotation(dxP, xP, xA, vort, area, indexStart, indexEnd)
	real(kreal), intent(out) :: dxP(:,:)
	real(kreal), intent(in) :: xP(:,:), xA(:,:), vort(:), area(:)
	integer(kint), intent(in) :: indexStart, indexEnd
	!
	integer(kint) :: j, k, n
	real(kreal) :: denom

	n = size(vort)

	dxP = 0.0_kreal
	do j = indexStart, indexEnd
		do k = 1, n
			denom = 2.0_kreal * PI * ( sum( (xP(:,j) - xA(:,k))*(xP(:,j) - xA(:,k))))
			dxp(1,j) = dxp(1,j) - ( xp(2,j) - xa(2,k) ) * vort(k) * area(k) / denom
			dxp(2,j) = dxp(2,j) + ( xp(1,j) - xa(1,k) ) * vort(k) * area(k) / denom
		enddo
	enddo
end subroutine

subroutine PlaneSmoothVel_noRotation(dxP, xP, xA, vort, area, indexStart, indexEnd)
	real(kreal), intent(out) :: dxP(:,:)
	real(kreal), intent(in)  :: xP(:,:), xa(:,:), vort(:), area(:)
	integer(kint), intent(in) :: indexStart, indexEnd
	!
	integer(kint) :: j, k, n
	real(kreal) :: denom

	n = size(vort)

	do j = indexStart, indexEnd
		do k = 1, n
			denom = 2.0_kreal * PI * ( sum( (xP(:,j) - xA(:,k)) * (xP(:,j) - xA(:,k))) + VELOCITY_SMOOTH * VELOCITY_SMOOTH)
			dxP(1,j) = dxP(1,j) - ( xp(2,j) - xA(2,k) ) * vort(k) * area(k) / denom
			dxP(2,j) = dxP(2,j) + ( xP(1,j) - xa(1,k) ) * vort(k) * area(k) / denom
		enddo
	enddo
end subroutine

function LeVeque93(xy, t)
	real(kreal) :: LeVeque93(2)
	real(kreal), intent(in) :: xy(2), t
	LeVeque93(1) = - xy(2) + 0.5_kreal
	LeVeque93(2) =   xy(1) - 0.5_kreal
	LeVeque93 = LeVeque93 / (2.0_kreal * PI)
end function

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
