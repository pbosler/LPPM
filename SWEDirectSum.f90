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
public LogStats
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
	type(SSRFPACKData), pointer :: velocitySource => null(), &
								   thicknessSource1 => null(), &
								   thicknessSource2 => null()
	!
	! Runge-Kutta variables
	!
	real(kreal), pointer :: particlesInput(:,:) => null(), &		! particle positions
							particlesStage1(:,:) => null(), &		! particle positions
							particlesStage2(:,:) => null(), &		! particle positions
							particlesStage3(:,:) => null(), &		! particle positions
							particlesStage4(:,:) => null(), &		! particle positions
							newParticlesX(:,:) => null(), &			! particle positions
							particlesVelocity(:,:) => null(), &		! particle velocity working array
							particlesDoubleDotU(:) => null(), &		! particle velocity double dot product working array
							particlesLapH(:) => null(), &			! particle laplacian h working array
							particlesRelVortInput(:) => null(), &	! particle relative vorticity
							particlesRelVortStage1(:) => null(), &	! particle relative vorticity
							particlesRelVortStage2(:) => null(), &	! particle relative vorticity
							particlesRelVortStage3(:) => null(), &	! particle relative vorticity
							particlesRelVortStage4(:) => null(), &	! particle relative vorticity
							newparticlesRelVort(:) => null(), &		! particle relative vorticity
							particlesDivInput(:) => null(), &		! particle divergence
							particlesDivStage1(:) => null(), &		! particle divergence
							particlesDivStage2(:) => null(), &		! particle divergence
							particlesDivStage3(:) => null(), &		! particle divergence
							particlesDivStage4(:) => null(), &		! particle divergence
							newparticlesDiv(:) => null(), &			! particle divergence
							particlesHInput(:) => null(), &			! particle thickness
							particlesHStage1(:) => null(), &		! particle thickness
							particlesHStage2(:) => null(), &		! particle thickness
							particlesHStage3(:) => null(), &		! particle thickness
							particlesHStage4(:) => null(), &		! particle thickness
							newParticlesH(:) => null()				! particle thickness

	real(kreal), pointer ::	activePanelsInput(:,:) => null(), &		! panel positions
							activePanelsStage1(:,:) => null(), &	! panel positions
							activePanelsStage2(:,:) => null(), &	! panel positions
							activePanelsStage3(:,:) => null(), &	! panel positions
							activePanelsStage4(:,:) => null(), &	! panel positions
							newActivePanelsX(:,:) => null(), &		! panel positions
							activePanelsVelocity(:,:) => null(), &	! panel velocity working array
							activePanelsDoubleDotU(:) => null(), &	! panel velocity double dot product working array
							activePanelsLapH(:) => null(), &		! panel laplacian h working array
							activePanelsRelVortInput(:) => null(), &! panel relative vorticity
							activePanelsRelVortStage1(:) => null(),&! panel relative vorticity
							activePanelsRelVortStage2(:) => null(),&! panel relative vorticity
							activePanelsRelVortStage3(:) => null(),&! panel relative vorticity
							activePanelsRelVortStage4(:) => null(),&! panel relative vorticity
							newActivePanelsRelVort(:) => null(), &	! panel relative vorticity
							activePanelsDivInput(:) => null(), &	! panel divergence
							activePanelsDivStage1(:) => null(), &	! panel divergence
							activePanelsDivStage2(:) => null(), &	! panel divergence
							activePanelsDivStage3(:) => null(), &	! panel divergence
							activePanelsDivStage4(:) => null(), &	! panel divergence
							newActivePanelsDiv(:) => null(), &		! panel divergence
							activePanelsHInput(:) => null(), &		! panel thickness
							activePanelsHStage1(:) => null(), &		! panel thickness
							activePanelsHStage2(:) => null(), &		! panel thickness
							activePanelsHStage3(:) => null(), &		! panel thickness
							activePanelsHStage4(:) => null(), &		! panel thickness
							newActivePanelsH(:) => null(), &		! panel thickness
							activePanelsAreaInput(:) => null(), &	! panel area
							activePanelsAreaStage1(:) => null(), &	! panel area
							activePanelsAreaStage2(:) => null(), &	! panel area
							activePanelsAreaStage3(:) => null(), &	! panel area
							activePanelsAreaStage4(:) => null(), &	! panel area
							newActivePanelsArea(:) => null()		! panel area

	real(kreal), pointer ::	passivePanelsInput(:,:) => null(), &	! passive panel positions
							passivePanelsStage1(:,:) => null(), &	! passive panel positions
							passivePanelsStage2(:,:) => null(), &	! passive panel positions
							passivePanelsStage3(:,:) => null(), &	! passive panel positions
							passivePanelsStage4(:,:) => null(), &	! passive panel positions
							newPassivePanelsX(:,:) => null()		! passive panel positions
end type

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'SWE'
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

interface LogStats
	module procedure LogStatsSWERK4
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
	call New(self%delTri,aMesh%particles%x(:,1:aMesh%particles%N),amesh%particles%N,self%activePanels%x,self%activePanels%N_Active)

	allocate(self%velocitySource)
	call New(self%velocitySource,self%delTri,.TRUE.)

	allocate(self%thicknessSource1)
	call New(self%thicknessSource1,self%delTri,.FALSE.)
	allocate(self%thicknessSource2)
	call New(self%thicknessSource2,self%delTri,.TRUE.)

	!
	! allocate RK4 variables
	!
	allocate(self%particlesInput(3,nParticles))
	allocate(self%particlesStage1(3,nParticles))
	allocate(self%particlesStage2(3,nParticles))
	allocate(self%particlesStage3(3,nParticles))
	allocate(self%particlesStage4(3,nParticles))
	allocate(self%newParticlesX(3,nParticles))

	allocate(self%particlesRelVortInput(nParticles))
	allocate(self%particlesRelVortStage1(nParticles))
	allocate(self%particlesRelVortStage2(nParticles))
	allocate(self%particlesRelVortStage3(nParticles))
	allocate(self%particlesRelVortStage4(nParticles))
	allocate(self%newParticlesRelVort(nParticles))

	allocate(self%particlesDivInput(nParticles))
	allocate(self%particlesDivStage1(nParticles))
	allocate(self%particlesDivStage2(nParticles))
	allocate(self%particlesDivStage3(nParticles))
	allocate(self%particlesDivStage4(nParticles))
	allocate(self%newParticlesDiv(nParticles))

	allocate(self%particlesHInput(nParticles))
	allocate(self%particlesHStage1(nParticles))
	allocate(self%particlesHStage2(nParticles))
	allocate(self%particlesHStage3(nParticles))
	allocate(self%particlesHStage4(nParticles))
	allocate(self%newParticlesH(nParticles))

	allocate(self%particlesVelocity(3,nParticles))
	allocate(self%particlesDoubleDotU(nParticles))
	allocate(self%particlesLapH(nParticles))

	allocate(self%activePanelsInput(3,nActive))
	allocate(self%activePanelsStage1(3,nActive))
	allocate(self%activePanelsStage2(3,nActive))
	allocate(self%activePanelsStage3(3,nActive))
	allocate(self%activePanelsStage4(3,nActive))
	allocate(self%newActivePanelsX(3,nActive))

	allocate(self%activePanelsVelocity(3,nActive))
	allocate(self%activePanelsDoubleDotU(nActive))
	allocate(self%activePanelsLapH(nActive))

	allocate(self%activePanelsRelVortInput(nActive))
	allocate(self%activePanelsRelVortStage1(nActive))
	allocate(self%activePanelsRelVortStage2(nActive))
	allocate(self%activePanelsRelVortStage3(nActive))
	allocate(self%activePanelsRelVortStage4(nActive))
	allocate(self%newactivePanelsRelVort(nActive))

	allocate(self%activePanelsDivInput(nActive))
	allocate(self%activePanelsDivStage1(nActive))
	allocate(self%activePanelsDivStage2(nActive))
	allocate(self%activePanelsDivStage3(nActive))
	allocate(self%activePanelsDivStage4(nActive))
	allocate(self%newactivePanelsDiv(nActive))

	allocate(self%activePanelsHInput(nActive))
	allocate(self%activePanelsHStage1(nActive))
	allocate(self%activePanelsHStage2(nActive))
	allocate(self%activePanelsHStage3(nActive))
	allocate(self%activePanelsHStage4(nActive))
	allocate(self%newactivePanelsH(nActive))

	allocate(self%activePanelsAreaInput(nActive))
	allocate(self%activePanelsAreaStage1(nActive))
	allocate(self%activePanelsAreaStage2(nActive))
	allocate(self%activePanelsAreaStage3(nActive))
	allocate(self%activePanelsAreaStage4(nActive))
	allocate(self%newactivePanelsArea(nActive))

	allocate(self%passivePanelsInput(3,nPassive))
	allocate(self%passivePanelsStage1(3,nPassive))
	allocate(self%passivePanelsStage2(3,nPassive))
	allocate(self%passivePanelsStage3(3,nPassive))
	allocate(self%passivePanelsStage4(3,nPassive))
	allocate(self%newPassivePanelsX(3,nPassive))

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
	call Delete(self%velocitySource)
	call Delete(self%thicknessSource1)
	call Delete(self%thicknessSource2)
	deallocate(self%delTri)
	deallocate(self%velocitySource)
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

	deallocate(self%particlesVelocity)
	deallocate(self%particlesDoubleDotU)
	deallocate(self%particlesLapH)

	deallocate(self%particlesRelVortInput)
	deallocate(self%particlesRelVortStage1)
	deallocate(self%particlesRelVortStage2)
	deallocate(self%particlesRelVortStage3)
	deallocate(self%particlesRelVortStage4)
	deallocate(self%newParticlesRelVort)

	deallocate(self%particlesDivInput)
	deallocate(self%particlesDivStage1)
	deallocate(self%particlesDivStage2)
	deallocate(self%particlesDivStage3)
	deallocate(self%particlesDivStage4)
	deallocate(self%newParticlesDiv)

	deallocate(self%particlesHInput)
	deallocate(self%particlesHStage1)
	deallocate(self%particlesHStage2)
	deallocate(self%particlesHStage3)
	deallocate(self%particlesHStage4)
	deallocate(self%newParticlesH)

	deallocate(self%activePanelsInput)
	deallocate(self%activePanelsStage1)
	deallocate(self%activePanelsStage2)
	deallocate(self%activePanelsStage3)
	deallocate(self%activePanelsStage4)
	deallocate(self%newActivePanelsX)

	deallocate(self%activePanelsVelocity)
	deallocate(self%activePanelsDoubleDotU)
	deallocate(self%activePanelsLapH)

	deallocate(self%activePanelsRelVortInput)
	deallocate(self%activePanelsRelVortStage1)
	deallocate(self%activePanelsRelVortStage2)
	deallocate(self%activePanelsRelVortStage3)
	deallocate(self%activePanelsRelVortStage4)
	deallocate(self%newActivePanelsRelVort)

	deallocate(self%activePanelsDivInput)
	deallocate(self%activePanelsDivStage1)
	deallocate(self%activePanelsDivStage2)
	deallocate(self%activePanelsDivStage3)
	deallocate(self%activePanelsDivStage4)
	deallocate(self%newActivePanelsDiv)

	deallocate(self%activePanelsHInput)
	deallocate(self%activePanelsHStage1)
	deallocate(self%activePanelsHStage2)
	deallocate(self%activePanelsHStage3)
	deallocate(self%activePanelsHStage4)
	deallocate(self%newActivePanelsH)

	deallocate(self%activePanelsAreaInput)
	deallocate(self%activePanelsAreaStage1)
	deallocate(self%activePanelsAreaStage2)
	deallocate(self%activePanelsAreaStage3)
	deallocate(self%activePanelsAreaStage4)
	deallocate(self%newActivePanelsArea)

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
subroutine LogStatsSWERK4(self, aLog, message)
	type(SWERK4Data), intent(in) :: self
	type(Logger), intent(inout) :: aLog
	character(len=*), intent(in), optional :: message
	!
	character(len=64) :: key
	if ( present(message) ) then
		call StartSection(aLog, 'SWERK4 Stats :', message)
	else
		call StartSection(aLog,'SWERK4 Stats :')
	endif
	key = 'size(activeInput) = '
	call LogMessage(alog,TRACE_LOGGING_LEVEL,key,size(self%activePanelsInput,2))
	key = 'activePanelsIndexStart(procRank) = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,self%activePanelsIndexStart(procRank))
	key = 'activePanelsIndexEnd(procRank) = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,self%activePanelsIndexEnd(procRank))
	key = 'size(particlesInput) = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,size(self%particlesInput,2))
	key = 'particlesIndexStart(procRank) = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,self%particlesIndexStart(procRank))
	key = 'particlesIndexEnd(procRank) = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,self%particlesIndexEnd(procRank))
	key = 'size(passivePanelsInput) = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,size(self%passivePanelsInput,2))
	key = 'passivePanelsIndexStart(procRank) = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,self%passivePanelsIndexStart(procRank))
	key = 'passivePanelsIndexEnd(procRank) = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,self%passivePanelsIndexEnd(procRank))


	call EndSection(aLog)
end subroutine

subroutine SWERK4Timestep(self, aMesh, dt, procRank, nProcs)
	type(SWERK4Data), intent(inout) :: self
	type(SphereMesh), intent(inout) :: aMesh
	real(kreal), intent(in) :: dt
	integer(kint), intent(in) :: procRank, nProcs
	!
	type(Particles), pointer :: aParticles
	integer(kint) :: j, errCode

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'entering SWERK4Timestep.')
	!call LogStats(self,log)

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
	self%activePanelsRelVortInput = self%activePanels%relVort
	self%activePanelsDivInput = self%activePanels%div
	self%activePanelsHInput = self%activePanels%h
	self%activePanelsAreaInput = self%activePanels%area

	self%passivePanelsInput = self%passivePanels%x

	!call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'stage 1 ready.')

	!
	! PARALLEL : compute particle velocities
	!
	call SWEVelocityActive(self%activePanelsVelocity, self%activePanelsInput, &
				self%activePanelsRelVortInput, self%activePanelsDivInput, self%activePanelsAreaInput, &
				self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank))
	call SWEVelocityPassive(self%particlesVelocity, self%particlesInput, self%activePanelsInput, &
				self%activePanelsRelVortInput, self%activePanelsDivInput, self%activePanelsAreaInput, &
				self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))
	call SWEVelocitySmooth(self%passivePanelsStage1, self%passivePanelsInput, self%activePanelsInput, &
				self%activePanelsRelVortInput, self%activePanelsDivInput, self%activePanelsAreaInput, &
				self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank))

	!
	! broadcast velocities (complete set needed for divergence equation)
	!
	do j=0, nProcs-1
		call MPI_BCAST(self%particlesVelocity(:,self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
					   3*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsVelocity(:,self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
					   3*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%passivePanelsStage1(:,self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
					   3*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
	enddo

	!call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'stage 1 velocity done.')
	!
	! END PARALLEL
	!
	call UpdateNodePositions(self%delTri, self%particlesInput, aMesh%particles%N, self%activePanelsInput, self%activePanels%N_Active)
	call SetSourceVelocity(self%velocitySource, self%delTri, self%particlesVelocity, aMesh%particles%N, &
		self%activePanelsVelocity, self%activePanels%N_Active)
	call SetSourceH(self%thicknessSource1, self%delTri, self%particlesHInput, aMesh%particles%N,&
		 self%activePanelsHInput, self%activePanels%N_Active)
	!call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'stage1 derivative sources ready.')

	call ComputeDoubleDotU(self%particlesDoubleDotU, aParticles%N, self%activePanelsDoubleDotU, self%activePanels%N, self%velocitySource)
	call ComputeLaplacianH(self%particlesLapH, aParticles%N, self%activePanelsLapH, self%activePanels%N, self%delTri, self%thicknessSource1)
	!
	! PARALLEL : compute RHS for vorticity, divergence, h, and area
	!
	do j=self%particlesIndexStart(procRank),self%particlesIndexEnd(procRank)
		self%particlesRelVortStage1(j) = -(self%particlesRelVortInput(j) + 2.0_kreal*OMEGA*self%particlesInput(3,j)/EARTH_RADIUS)*&
					self%particlesDivInput(j) - 2.0_kreal*OMEGA*self%particlesVelocity(3,j)/EARTH_RADIUS
		self%particlesDivStage1(j) = -self%particlesDoubleDotU(j) + &
			(2.0_kreal*OMEGA*self%particlesInput(3,j)/EARTH_RADIUS)*self%particlesRelVortInput(j) - &
			GRAV*self%particlesLapH(j)
		self%particlesHStage1(j) = -self%particlesDivInput(j)*self%particlesHInput(j)
	enddo

	do j=self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank)
		self%activePanelsRelVortStage1(j) = -(self%activePanelsRelVortInput(j) + 2.0_kreal*OMEGA*self%activePanelsInput(3,j)/EARTH_RADIUS)*&
			self%activePanelsDivInput(j) - 2.0_kreal*OMEGA*self%activePanelsVelocity(3,j)/EARTH_RADIUS
		self%activePanelsDivStage1(j) = - self%activePanelsDoubleDotU(j) + &
			(2.0_kreal*OMEGA*self%activePanelsInput(3,j)/EARTH_RADIUS)*self%activePanelsRelVortInput(j) - &
			GRAV*self%activePanelsLapH(j)
		self%activePanelsHStage1(j) = -self%activePanelsDivInput(j)*self%activePanelsHInput(j)
		self%activePanelsAreaStage1(j) = self%activePanelsDivInput(j)*self%activePanelsAreaInput(j)
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
		call MPI_BCAST(self%activePanelsRelVortStage1(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsDivStage1(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsHStage1(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsAreaStage1(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)

		call MPI_BCAST(self%activePanelsDoubleDotU(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsLapH(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
	enddo
	!call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'stage 1 data done.')
	!
	! END PARALLEL
	!

	!
	! STAGE 1 ONLY : Store velocity and kinetic energy
	!
	do j=1,aParticles%N
		aParticles%ke(j) = sum(self%particlesVelocity(:,j)*self%particlesVelocity(:,j))
		aParticles%u(:,j) = self%particlesVelocity(:,j)
		aParticles%tracer(j,3) = self%particlesDoubleDotU(j)
		aParticles%tracer(j,4) = self%particlesLapH(j)
	enddo
	do j=1,self%activePanels%N
		self%activePanels%ke(j) = sum(self%activePanelsVelocity(:,j)*self%activePanelsVelocity(:,j))
		self%activePanels%u(:,j) = self%activePanelsVelocity(:,j)
		self%activePanels%tracer(j,3) = self%activePanelsDoubleDotU(j)
		self%activePanels%tracer(j,4) = self%activePanelsLapH(j)
	enddo

!	call LogMessage(log,DEBUG_LOGGING_LEVEL,'particles max(abs(drelVort)) = ',maxval(abs(self%particlesRelVortStage1)))
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,'panels max(abs(drelVort)) = ',maxval(abs(self%activePanelsRelVortStage1)))
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,'particles max(abs(dDiv)) = ',maxval(abs(self%particlesDivStage1)))
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,'panels max(abs(dDiv)) = ',maxval(abs(self%activePanelsDivStage1)))
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,'particles max(abs(dH)) = ',maxval(abs(self%particlesHStage1)))
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,'panels max(abs(dH)) = ',maxval(abs(self%activePanelsHStage1)))
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,'panels max(abs(dArea)) = ',maxval(abs(self%activePanelsAreaStage1)))

	self%particlesStage1 = dt*self%particlesVelocity
	self%particlesRelVortStage1 = dt*self%particlesRelVortStage1
	self%particlesDivStage1 = dt*self%particlesDivStage1
	self%particlesHStage1 = dt*self%particlesHStage1

	self%activePanelsStage1 = dt*self%activePanelsVelocity
	self%activePanelsRelVortStage1 = dt*self%activePanelsRelVortStage1
	self%activePanelsDivStage1 = dt*self%activePanelsDivStage1
	self%activePanelsHStage1 = dt*self%activePanelsHStage1
	self%activePanelsAreaStage1 = dt*self%activePanelsAreaStage1

	self%passivePanelsStage1 = dt*self%passivePanelsStage1

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 2       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 2
	self%particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*self%particlesStage1
	self%particlesRelVortInput = aParticles%relVort(1:aParticles%N) + 0.5_kreal*self%particlesRelVortStage1
	self%particlesDivInput = aParticles%div(1:aParticles%N) + 0.5_kreal*self%particlesDivStage1
	self%particlesHInput = aParticles%h(1:aParticles%N) + 0.5_kreal*self%particlesHStage1

	self%activePanelsInput = self%activePanels%x + 0.5_kreal*self%activePanelsStage1
	self%activePanelsRelVortInput = self%activePanels%relVort + 0.5_kreal*self%activePanelsRelVortStage1
	self%activePanelsDivInput = self%activePanels%div + 0.5_kreal*self%activePanelsDivStage1
	self%activePanelsHInput = self%activePanels%h + 0.5_kreal*self%activePanelsHStage1
	self%activePanelsAreaInput = self%activePanels%area + 0.5_kreal*self%activePanelsAreaStage1

	self%passivePanelsInput = self%passivePanels%x + 0.5_kreal*self%passivePanelsStage1

!	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'stage 2 ready.')
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,'stage 2 surf area err = ',abs(4.0_kreal*EARTH_RADIUS*EARTH_RADIUS - sum(self%activePanelsAreaInput))/4.0_kreal/EARTH_RADIUS/EARTH_RADIUS)
	!
	! PARALLEL : compute particle velocities
	!
	call SWEVelocityActive(self%activePanelsVelocity, self%activePanelsInput, &
				self%activePanelsRelVortInput, self%activePanelsDivInput, self%activePanelsAreaInput, &
				self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank))
	call SWEVelocityPassive(self%particlesVelocity, self%particlesInput, self%activePanelsInput, &
				self%activePanelsRelVortInput, self%activePanelsDivInput, self%activePanelsAreaInput, &
				self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))
	call SWEVelocitySmooth(self%passivePanelsStage2, self%passivePanelsInput, self%activePanelsInput, &
				self%activePanelsRelVortInput, self%activePanelsDivInput, self%activePanelsAreaInput, &
				self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank))
	!
	! broadcast velocities (needed for divergence equation)
	!
	do j=0, nProcs-1
		call MPI_BCAST(self%particlesVelocity(:,self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
					   3*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsVelocity(:,self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
					   3*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%passivePanelsStage2(:,self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
					   3*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
	enddo
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'stage 2 velocity done.')
	!
	! END PARALLEL
	!
	call UpdateNodePositions(self%delTri, self%particlesInput, aMesh%particles%N, self%activePanelsInput, self%activePanels%N_Active)
	call SetSourceVelocity(self%velocitySource, self%delTri, self%particlesVelocity, aMesh%particles%N, &
		self%activePanelsVelocity, self%activePanels%N_Active)
	call SetSourceH(self%thicknessSource1, self%delTri, self%particlesHInput, aMesh%particles%N,&
		 self%activePanelsHInput, self%activePanels%N_Active)


	!
	! PARALLEL : compute divergence equation forcing terms
	!
	call ComputeDoubleDotU(self%particlesDoubleDotU, aParticles%N, self%activePanelsDoubleDotU, self%activePanels%N, self%velocitySource)
	call ComputeLaplacianH(self%particlesLapH, aParticles%N, self%activePanelsLapH, self%activePanels%N, self%delTri, self%thicknessSource1)

	!
	! compute RHS for vorticity, divergence, h, and area
	!
	do j=self%particlesIndexStart(procRank),self%particlesIndexEnd(procRank)
		self%particlesRelVortStage2(j) = -(self%particlesRelVortInput(j) + 2.0_kreal*OMEGA*self%particlesInput(3,j)/EARTH_RADIUS)*&
					self%particlesDivInput(j) - 2.0_kreal*OMEGA*self%particlesVelocity(3,j)/EARTH_RADIUS
		self%particlesDivStage2(j) = -self%particlesDoubleDotU(j) + &
			2.0_kreal*OMEGA*self%particlesInput(3,j)/EARTH_RADIUS*self%particlesRelVortInput(j) - &
			GRAV*self%particlesLapH(j)
		self%particlesHStage2(j) = -self%particlesDivInput(j)*self%particlesHInput(j)
	enddo

	do j=self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank)
		self%activePanelsRelVortStage2(j) = -(self%activePanelsRelVortInput(j) + 2.0_kreal*OMEGA*self%activePanelsInput(3,j)/EARTH_RADIUS)*&
			self%activePanelsDivInput(j) - 2.0_kreal*OMEGA*self%activePanelsVelocity(3,j)/EARTH_RADIUS
		self%activePanelsDivStage2(j) = - self%activePanelsDoubleDotU(j) + &
			2.0_kreal*OMEGA*self%activePanelsInput(3,j)/EARTH_RADIUS*self%activePanelsRelVortInput(j) - &
			GRAV*self%activePanelsLapH(j)
		self%activePanelsHStage2(j) = -self%activePanelsDivInput(j)*self%activePanelsHInput(j)
		self%activePanelsAreaStage2(j) = self%activePanelsDivInput(j)*self%activePanelsAreaInput(j)
	enddo
	!
	! broadcast stage2 vorticity, divergence, h, area
	!
	do j=1, nProcs-1
		call MPI_BCAST(self%particlesRelVortStage2(self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
						self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%particlesDivStage2(self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
						self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%particlesHStage2(self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
						self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsRelVortStage2(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsDivStage2(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsHStage2(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsAreaStage2(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
	enddo
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'stage 2 data done.')
	!
	! END PARALLEL
	!

	self%particlesStage2 = dt*self%particlesVelocity
	self%particlesRelVortStage2 = dt*self%particlesRelVortStage2
	self%particlesDivStage2 = dt*self%particlesDivStage2
	self%particlesHStage2 = dt*self%particlesHStage2

	self%activePanelsStage2 = dt*self%activePanelsVelocity
	self%activePanelsRelVortStage2 = dt*self%activePanelsRelVortStage2
	self%activePanelsDivStage2 = dt*self%activePanelsDivStage2
	self%activePanelsHStage2 = dt*self%activePanelsHStage2
	self%activePanelsAreaStage2 = dt*self%activePanelsAreaStage2

	self%passivePanelsStage2 = dt*self%passivePanelsStage2

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 3       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 3
	self%particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*self%particlesStage2
	self%particlesRelVortInput = aParticles%relVort(1:aParticles%N) + 0.5_kreal*self%particlesRelVortStage2
	self%particlesDivInput = aParticles%div(1:aParticles%N) + 0.5_kreal*self%particlesDivStage2
	self%particlesHInput = aParticles%h(1:aParticles%N) + 0.5_kreal*self%particlesHStage2

	self%activePanelsInput = self%activePanels%x + 0.5_kreal*self%activePanelsStage2
	self%activePanelsRelVortINput = self%activePanels%relVort + 0.5_kreal*self%activePanelsRelVortStage2
	self%activePanelsDivINput = self%activePanels%div + 0.5_kreal*self%activePanelsDivStage2
	self%activePanelsHInput = self%activePanels%h + 0.5_kreal*self%activePanelsHStage2
	self%activePanelsAreaINput = self%activePanels%area + 0.5_kreal*self%activePanelsAreaStage2

	self%passivePanelsInput = self%passivePanels%x + 0.5_kreal*self%passivePanelsStage2
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'stage 3 ready.')
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,'stage 3 surf area err = ',abs(4.0_kreal*EARTH_RADIUS*EARTH_RADIUS - sum(self%activePanelsAreaInput))/4.0_kreal/EARTH_RADIUS/EARTH_RADIUS)
	!
	! PARALLEL : compute particle velocities
	!
	call SWEVelocityActive(self%activePanelsVelocity, self%activePanelsInput, &
				self%activePanelsRelVortInput, self%activePanelsDivInput, self%activePanelsAreaInput, &
				self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank))
	call SWEVelocityPassive(self%particlesVelocity, self%particlesInput, self%activePanelsInput, &
				self%activePanelsRelVortInput, self%activePanelsDivInput, self%activePanelsAreaInput, &
				self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))
	call SWEVelocitySmooth(self%passivePanelsStage3, self%passivePanelsInput, self%activePanelsInput, &
				self%activePanelsRelVortInput, self%activePanelsDivInput, self%activePanelsAreaInput, &
				self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank))
	!
	! broadcast velocities (needed for divergence equation)
	!
	do j=0, nProcs-1
		call MPI_BCAST(self%particlesVelocity(:,self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
					   3*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsVelocity(:,self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
					   3*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%passivePanelsStage3(:,self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
					   3*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
	enddo
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'stage 3 velocity done.')
	!
	! END PARALLEL
	!
	call UpdateNodePositions(self%delTri, self%particlesInput, aMesh%particles%N, self%activePanelsInput, self%activePanels%N_Active)
	call SetSourceVelocity(self%velocitySource, self%delTri, self%particlesVelocity, aMesh%particles%N, &
		self%activePanelsVelocity, self%activePanels%N_Active)
	call SetSourceH(self%thicknessSource1, self%delTri, self%particlesHInput, aMesh%particles%N,&
		 self%activePanelsHInput, self%activePanels%N_Active)


	!
	! PARALLEL : compute divergence equation forcing terms
	!
	call ComputeDoubleDotU(self%particlesDoubleDotU, aParticles%N, self%activePanelsDoubleDotU, self%activePanels%N, self%velocitySource)
	call ComputeLaplacianH(self%particlesLapH, aParticles%N, self%activePanelsLapH, self%activePanels%N, self%delTri, self%thicknessSource1)

	!
	! compute RHS for vorticity, divergence, h, and area
	!
	do j=self%particlesIndexStart(procRank),self%particlesIndexEnd(procRank)
		self%particlesRelVortStage3(j) = -(self%particlesRelVortInput(j) + 2.0_kreal*OMEGA*self%particlesInput(3,j)/EARTH_RADIUS)*&
					self%particlesDivInput(j) - 2.0_kreal*OMEGA*self%particlesVelocity(3,j)/EARTH_RADIUS
		self%particlesDivStage3(j) = -self%particlesDoubleDotU(j) + &
			2.0_kreal*OMEGA*self%particlesInput(3,j)/EARTH_RADIUS*self%particlesRelVortInput(j) - &
			GRAV*self%particlesLapH(j)
		self%particlesHStage3(j) = -self%particlesDivInput(j)*self%particlesHInput(j)
	enddo

	do j=self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank)
		self%activePanelsRelVortStage3(j) = -(self%activePanelsRelVortInput(j) + 2.0_kreal*OMEGA*self%activePanelsInput(3,j)/EARTH_RADIUS)*&
			self%activePanelsDivInput(j) - 2.0_kreal*OMEGA*self%activePanelsVelocity(3,j)/EARTH_RADIUS
		self%activePanelsDivStage3(j) = - self%activePanelsDoubleDotU(j) + &
			2.0_kreal*OMEGA*self%activePanelsInput(3,j)/EARTH_RADIUS*self%activePanelsRelVortInput(j) - &
			GRAV*self%activePanelsLapH(j)
		self%activePanelsHStage3(j) = -self%activePanelsDivInput(j)*self%activePanelsHInput(j)
		self%activePanelsAreaStage3(j) = self%activePanelsDivInput(j)*self%activePanelsAreaInput(j)
	enddo
	!
	! broadcast stage3 vorticity, divergence, h, area
	!
	do j=1, nProcs-1
		call MPI_BCAST(self%particlesRelVortStage3(self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
						self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%particlesDivStage3(self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
						self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%particlesHStage3(self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
						self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsRelVortStage3(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsDivStage3(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsHStage3(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsAreaStage3(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
	enddo
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'stage 3 data done.')
	!
	! END PARALLEL
	!

	self%particlesStage3 = dt*self%particlesVelocity
	self%particlesRelVortStage3 = dt*self%particlesRelVortStage3
	self%particlesDivStage3 = dt*self%particlesDivStage3
	self%particlesHStage3 = dt*self%particlesHStage3

	self%activePanelsStage3 = dt*self%activePanelsVelocity
	self%activePanelsRelVortStage3 = dt*self%activePanelsRelVortStage3
	self%activePanelsDivStage3 = dt*self%activePanelsDivStage3
	self%activePanelsHStage3 = dt*self%activePanelsHStage3
	self%activePanelsAreaStage3 = dt*self%activePanelsAreaStage3

	self%passivePanelsStage3 = dt*self%passivePanelsStage3


	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 4      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 4
	self%particlesInput = aParticles%x(:,1:aParticles%N) + self%particlesStage3
	self%particlesRelVortInput = aParticles%relVort(1:aParticles%N) + self%particlesRelVortStage3
	self%particlesDivInput = aParticles%div(1:aParticles%N) + self%particlesDivStage3
	self%particlesHInput = aParticles%h(1:aParticles%N) + self%particlesHStage3

	self%activePanelsInput = self%activePanels%x + self%activePanelsStage3
	self%activePanelsRelVortInput = self%activePanels%relVort + self%activePanelsRelVortStage3
	self%activePanelsDivInput = self%activePanels%div + self%activePanelsDivStage3
	self%activePanelsHInput = self%activePanels%h + self%activePanelsHStage3
	self%activePanelsAreaInput = self%activePanels%area + self%activePanelsAreaStage3

	self%passivePanelsInput = self%passivePanels%x + self%passivePanelsStage3
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'stage 4 ready.')
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,'stage 4 surf area err = ',abs(4.0_kreal*EARTH_RADIUS*EARTH_RADIUS - sum(self%activePanelsAreaInput))/4.0_kreal/EARTH_RADIUS/EARTH_RADIUS)
	!
	! PARALLEL : compute particle velocities
	!
	call SWEVelocityActive(self%activePanelsVelocity, self%activePanelsInput, &
				self%activePanelsRelVortInput, self%activePanelsDivInput, self%activePanelsAreaInput, &
				self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank))
	call SWEVelocityPassive(self%particlesVelocity, self%particlesInput, self%activePanelsInput, &
				self%activePanelsRelVortInput, self%activePanelsDivInput, self%activePanelsAreaInput, &
				self%particlesIndexStart(procRank), self%particlesIndexEnd(procRank))
	call SWEVelocitySmooth(self%passivePanelsStage4, self%passivePanelsInput, self%activePanelsInput, &
				self%activePanelsRelVortInput, self%activePanelsDivInput, self%activePanelsAreaInput, &
				self%passivePanelsIndexStart(procRank), self%passivePanelsIndexEnd(procRank))
	!
	! broadcast velocities (needed for divergence equation)
	!
	do j=0, nProcs-1
		call MPI_BCAST(self%particlesVelocity(:,self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
					   3*self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsVelocity(:,self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
					   3*self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%passivePanelsStage4(:,self%passivePanelsIndexStart(j):self%passivePanelsIndexEnd(j)), &
					   3*self%passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
	enddo
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'stage 4 velocity done.')
	!
	! END PARALLEL
	!
	call UpdateNodePositions(self%delTri, self%particlesInput, aMesh%particles%N, self%activePanelsInput, self%activePanels%N_Active)
	call SetSourceVelocity(self%velocitySource, self%delTri, self%particlesVelocity, aMesh%particles%N, &
		self%activePanelsVelocity, self%activePanels%N_Active)
	call SetSourceH(self%thicknessSource1, self%delTri, self%particlesHInput, aMesh%particles%N,&
		 self%activePanelsHInput, self%activePanels%N_Active)


	!
	! PARALLEL : compute divergence equation forcing terms
	!
	call ComputeDoubleDotU(self%particlesDoubleDotU, aParticles%N, self%activePanelsDoubleDotU, self%activePanels%N, self%velocitySource)
	call ComputeLaplacianH(self%particlesLapH, aParticles%N, self%activePanelsLapH, self%activePanels%N, self%delTri, self%thicknessSource1)
	!
	! compute RHS for vorticity, divergence, h, and area
	!
	do j=self%particlesIndexStart(procRank),self%particlesIndexEnd(procRank)
		self%particlesRelVortStage4(j) = -(self%particlesRelVortInput(j) + 2.0_kreal*OMEGA*self%particlesInput(3,j)/EARTH_RADIUS)*&
					self%particlesDivInput(j) - 2.0_kreal*OMEGA*self%particlesVelocity(3,j)/EARTH_RADIUS
		self%particlesDivStage4(j) = -self%particlesDoubleDotU(j) + &
			2.0_kreal*OMEGA*self%particlesInput(3,j)/EARTH_RADIUS*self%particlesRelVortInput(j) - &
			GRAV*self%particlesLapH(j)
		self%particlesHStage4(j) = -self%particlesDivInput(j)*self%particlesHInput(j)
	enddo

	do j=self%activePanelsIndexStart(procRank), self%activePanelsIndexEnd(procRank)
		self%activePanelsRelVortStage4(j) = -(self%activePanelsRelVortInput(j) + 2.0_kreal*OMEGA*self%activePanelsInput(3,j)/EARTH_RADIUS)*&
			self%activePanelsDivInput(j) - 2.0_kreal*OMEGA*self%activePanelsVelocity(3,j)/EARTH_RADIUS
		self%activePanelsDivStage4(j) = - self%activePanelsDoubleDotU(j) + &
			2.0_kreal*OMEGA*self%activePanelsInput(3,j)/EARTH_RADIUS*self%activePanelsRelVortInput(j) - &
			GRAV*self%activePanelsLapH(j)
		self%activePanelsHStage4(j) = -self%activePanelsDivInput(j)*self%activePanelsHInput(j)
		self%activePanelsAreaStage4(j) = self%activePanelsDivInput(j)*self%activePanelsAreaInput(j)
	enddo
	!
	! broadcast stage4 vorticity, divergence, h, area
	!
	do j=1, nProcs-1
		call MPI_BCAST(self%particlesRelVortStage4(self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
						self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%particlesDivStage4(self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
						self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%particlesHStage4(self%particlesIndexStart(j):self%particlesIndexEnd(j)), &
						self%particlesMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsRelVortStage4(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsDivStage4(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsHStage4(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
		call MPI_BCAST(self%activePanelsAreaStage4(self%activePanelsIndexStart(j):self%activePanelsIndexEnd(j)), &
						self%activePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, MPI_COMM_WORLD, errCode)
	enddo
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'stage 4 data done.')
	!
	! END PARALLEL
	!

	self%particlesStage4 = dt*self%particlesVelocity
	self%particlesRelVortStage4 = dt*self%particlesRelVortStage4
	self%particlesDivStage4 = dt*self%particlesDivStage4
	self%particlesHStage4 = dt*self%particlesHStage4

	self%activePanelsStage4 = dt*self%activePanelsVelocity
	self%activePanelsRelVortStage4 = dt*self%activePanelsRelVortStage4
	self%activePanelsDivStage4 = dt*self%activePanelsDivStage4
	self%activePanelsHStage4 = dt*self%activePanelsHStage4
	self%activePanelsAreaStage4 = dt*self%activePanelsAreaStage4

	self%passivePanelsStage4 = dt*self%passivePanelsStage4

	!!!!!!!!!!!!!!!!!
	! 	RK update   !
	!!!!!!!!!!!!!!!!!

	self%newParticlesX = aParticles%x(:,1:aParticles%N) + self%particlesStage1/6.0_kreal + self%particlesStage2/3.0_kreal + &
		self%particlesStage3/3.0_kreal + self%particlesStage4/6.0_kreal
	self%newParticlesRelVort = aParticles%relVort(1:aParticles%N) + self%particlesRelVOrtStage1/6.0_kreal + self%particlesRelVortStage2/3.0_kreal + &
		self%particlesRelVortStage3/3.0_kreal + self%particlesRelVortStage4/6.0_kreal
	self%newParticlesDiv = aParticles%div(1:aParticles%N) + self%particlesDivStage1/6.0_kreal + self%particlesDivStage2/3.0_kreal + &
		self%particlesDivStage3/3.0_kreal + self%particlesDivStage4/6.0_kreal
	self%newParticlesH = aParticles%h(1:aParticles%N) + self%particlesHStage1/6.0_kreal + self%particlesHStage2/6.0_kreal + &
		self%particlesHStage3/3.0_kreal + self%particlesHStage4/6.0_kreal

	self%newActivePanelsX = self%activePanels%x + self%activePanelsStage1/6.0_kreal + self%activePanelsStage2/3.0_kreal + &
		self%activePanelsStage3/3.0_kreal + self%activePanelsStage4/6.0_kreal
	self%newActivePanelsRelVort = self%activePanels%relVort + self%activePanelsRelVortStage1/6.0_kreal + self%activePanelsRelVortStage2/3.0_kreal + &
		self%activePanelsRelVortStage3/3.0_kreal + self%activePanelsRelVortStage4/6.0_kreal
	self%newActivePanelsDiv = self%activePanels%div + self%activePanelsDivStage1/6.0_kreal + self%activePanelsDivStage2/3.0_kreal + &
		self%activePanelsDivStage3/3.0_kreal + self%activePanelsDivStage4/6.0_kreal
	self%newActivePanelsH = self%activePanels%h + self%activePanelsHStage1/6.0_kreal + self%activePanelsHStage2/3.0_kreal + &
		self%activePanelsHStage3/3.0_kreal + self%activePanelsHStage4/6.0_kreal
	self%newActivePanelsArea = self%activePanels%area + self%activePanelsAreaStage1/6.0_kreal + self%activePanelsAreaStage2/3.0_kreal + &
		self%activePanelsAreaStage3/3.0_kreal + self%activePanelsAreaStage4/6.0_kreal

	self%newPassivePanelsX = self%passivePanels%x + self%passivePanelsStage1/6.0_kreal + self%passivePanelsStage2/3.0_kreal + &
		self%passivePanelsStage3/3.0_kreal + self%passivePanelsStage4/6.0_kreal

	aParticles%x(:,1:aParticles%N) = self%newParticlesX
	aParticles%relVort(1:aParticles%N) = self%newParticlesRelVort
	aParticles%div(1:aParticles%N) = self%newParticlesDiv
	aParticles%h(1:aParticles%N) = self%newParticlesH

	self%activePanels%x = self%newActivePanelsX
	self%activePanels%relVort = self%newActivePanelsRelVort
	self%activePanels%div = self%newActivePanelsDiv
	self%activePanels%h = self%newActivePanelsH
	self%activePanels%area = self%newActivePanelsArea

	self%passivePanels%x = self%newPassivePanelsX

	call ScatterPanels(aMesh%panels, self%activePanels, self%activeMap, self%passivePanels, self%passiveMap)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'... timestep complete.')
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

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'entering ZeroRK4.')

	self%particlesInput = 0.0_kreal
	self%particlesStage1 = 0.0_kreal
	self%particlesStage2 = 0.0_kreal
	self%particlesStage3 = 0.0_kreal
	self%particlesStage4 = 0.0_kreal
	self%newParticlesX = 0.0_kreal

	self%particlesVelocity = 0.0_kreal
	self%particlesDoubleDotU = 0.0_kreal
	self%particlesLapH = 0.0_kreal

	self%particlesRelVortInput = 0.0_kreal
	self%particlesRelVortStage1 = 0.0_kreal
	self%particlesRelVortStage2 = 0.0_kreal
	self%particlesRelVortStage3 = 0.0_kreal
	self%particlesRelVortStage4 = 0.0_kreal
	self%newParticlesRelVort = 0.0_kreal

	self%particlesDivInput = 0.0_kreal
	self%particlesDivStage1 = 0.0_kreal
	self%particlesDivStage2 = 0.0_kreal
	self%particlesDivStage3 = 0.0_kreal
	self%particlesDivStage4 = 0.0_kreal
	self%newParticlesDiv = 0.0_kreal

	self%particlesHInput = 0.0_kreal
	self%particlesHStage1 = 0.0_kreal
	self%particlesHStage2 = 0.0_kreal
	self%particlesHStage3 = 0.0_kreal
	self%particlesHStage4 = 0.0_kreal
	self%newParticlesH = 0.0_kreal

	self%activePanelsInput = 0.0_kreal
	self%activePanelsStage1 = 0.0_kreal
	self%activePanelsStage2 = 0.0_kreal
	self%activePanelsStage3 = 0.0_kreal
	self%activePanelsStage4 = 0.0_kreal
	self%newActivePanelsX = 0.0_kreal

	self%activePanelsVelocity = 0.0_kreal
	self%activePanelsDoubleDotU = 0.0_kreal
	self%activePanelsLapH = 0.0_kreal

	self%activePanelsRelVortInput = 0.0_kreal
	self%activePanelsRelVortStage1 = 0.0_kreal
	self%activePanelsRelVortStage2 = 0.0_kreal
	self%activePanelsRelVortStage3 = 0.0_kreal
	self%activePanelsRelVortStage4 = 0.0_kreal
	self%newActivePanelsRelVort = 0.0_kreal

	self%activePanelsDivInput = 0.0_kreal
	self%activePanelsDivStage1 = 0.0_kreal
	self%activePanelsDivStage2 = 0.0_kreal
	self%activePanelsDivStage3 = 0.0_kreal
	self%activePanelsDivStage4 = 0.0_kreal
	self%newActivePanelsDiv = 0.0_kreal

	self%activePanelsHInput = 0.0_kreal
	self%activePanelsHStage1 = 0.0_kreal
	self%activePanelsHStage2 = 0.0_kreal
	self%activePanelsHStage3 = 0.0_kreal
	self%activePanelsHStage4 = 0.0_kreal
	self%newActivePanelsH = 0.0_kreal

	self%activePanelsAreaInput = 0.0_kreal
	self%activePanelsAreaStage1 = 0.0_kreal
	self%activePanelsAreaStage2 = 0.0_kreal
	self%activePanelsAreaStage3 = 0.0_kreal
	self%activePanelsAreaStage4 = 0.0_kreal
	self%newActivePanelsArea = 0.0_kreal

	self%passivePanelsInput = 0.0_kreal
	self%passivePanelsStage1 = 0.0_kreal
	self%passivePanelsStage2 = 0.0_kreal
	self%passivePanelsStage3 = 0.0_kreal
	self%passivePanelsStage4 = 0.0_kreal
	self%newPassivePanelsX = 0.0_kreal
end subroutine

!subroutine ComputeDoubleDotU(ddU, delTri, velocitySource1, indexStart, indexEnd)
!	real(kreal), intent(out) :: ddU(:)
!	type(STRIPACKData), intent(in) :: delTri
!	type(SSRFPACKData), intent(in) :: velocitySource1
!	integer(kint), intent(in) :: indexStart, indexEnd
!	!
!	integer(kint) :: j, errCode
!
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'entering ComputeDoubleDotU.')
!
!	do j=indexStart, indexEnd
!		ddU(j) = velocitySource1%grad1(1,j)*velocitySource1%grad1(1,j) + &
!				 velocitySource1%grad2(2,j)*velocitySource1%grad2(2,j) + &
!				 velocitySource1%grad3(3,j)*velocitySource1%grad3(3,j) + &
!				 2.0_kreal*velocitySource1%grad1(2,j)*velocitySource1%grad2(1,j) + &
!				 2.0_kreal*velocitySource1%grad1(3,j)*velocitySource1%grad3(1,j) + &
!				 2.0_kreal*velocitySource1%grad2(3,j)*velocitySource1%grad3(2,j)
!	enddo
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,'max(ddotU) = ',maxval(ddu))
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,'min(ddotU) = ',minval(ddu))
!end subroutine

subroutine ComputeDoubleDotU(particlesDDU, nParticles, activePanelsDDU, nActive, velocitySource)
	real(kreal), intent(out) :: particlesDDU(:), activePanelsDDU(:)
	integer(kint), intent(in) :: nParticles, nActive
	type(SSRFPACKData), intent(in) :: velocitySource
	!
	integer(kint) :: j

	do j=1,nActive
		activePanelsDDU(j) =  velocitySource%grad1(1,j)*velocitySource%grad1(1,j) + &
							  velocitySource%grad2(2,j)*velocitySource%grad2(2,j) + &
							  velocitySource%grad3(3,j)*velocitySource%grad3(3,j) + &
							  2.0_kreal*velocitySource%grad1(2,j)*velocitySource%grad2(1,j) + &
							  2.0_kreal*velocitySource%grad1(3,j)*velocitySource%grad3(1,j) + &
							  2.0_kreal*velocitySource%grad2(3,j)*velocitySource%grad3(2,j)
	enddo
	activePanelsDDU = activePanelsDDU/EARTH_RADIUS/EARTH_RADIUS
	do j=1,nParticles
		particlesDDU(j) = velocitySource%grad1(1,nActive+j)*velocitySource%grad1(1,nActive+j) + &
					  	  velocitySource%grad2(2,nActive+j)*velocitySource%grad2(2,nActive+j) + &
					  	  velocitySource%grad3(3,nActive+j)*velocitySource%grad3(3,nActive+j) + &
					  	  2.0_kreal*velocitySource%grad1(2,nActive+j)*velocitySource%grad2(1,nActive+j) + &
					  	  2.0_kreal*velocitySource%grad1(3,nActive+j)*velocitySource%grad3(1,nActive+j) + &
					  	  2.0_kreal*velocitySource%grad2(3,nActive+j)*velocitySource%grad3(2,nActive+j)
	enddo
	particlesDDU = particlesDDU/EARTH_RADIUS/EARTH_RADIUS
end subroutine

subroutine ComputeLaplacianH(particlesLapH, nParticles, activePanelsLapH, nActive, delTri, hSource)
	real(kreal), intent(out) :: particlesLapH(:), activePanelsLapH(:)
	type(STRIPACKData), intent(in) :: delTri
	type(SSRFPACKData), intent(in) :: hSource
	integer(kint), intent(in) :: nparticles, nactive
	!
	real(kreal) :: grad1grad(3), grad2grad(3), grad3grad(3)
	integer(kint) :: j, errCode
	do j=1,nactive
		call GRADL(delTri%n, j, delTri%x, delTri%y, delTri%z, hSource%grad1(1,:), &
			delTri%list, delTri%lptr, delTri%lend, grad1grad, errCode)
		call GRADL(delTri%n, j, delTri%x, delTri%y, delTri%z, hSource%grad1(2,:), &
			delTri%list, delTri%lptr, delTri%lend, grad2grad, errCode)
		call GRADL(delTri%n, j, delTri%x, delTri%y, delTri%z, hSource%grad1(3,:), &
			delTri%list, delTri%lptr, delTri%lend, grad3grad, errCode)
		activePanelsLapH(j) = grad1grad(1) + grad2grad(2) + grad3grad(3)
	enddo
	activePanelsLapH = activePanelsLapH/EARTH_RADIUS/EARTH_RADIUS
	do j=1,nParticles
		call GRADL(delTri%n, nActive + j, delTri%x, delTri%y, delTri%z, hSource%grad1(1,:), &
			delTri%list, delTri%lptr, delTri%lend, grad1grad, errCode)
		call GRADL(delTri%n, nActive + j, delTri%x, delTri%y, delTri%z, hSource%grad1(2,:), &
			delTri%list, delTri%lptr, delTri%lend, grad2grad, errCode)
		call GRADL(delTri%n, nActive + j, delTri%x, delTri%y, delTri%z, hSource%grad1(3,:), &
			delTri%list, delTri%lptr, delTri%lend, grad3grad, errCode)
		particlesLapH(j) = grad1grad(1) + grad2grad(2) + grad3grad(3)
	enddo
	particlesLapH = particlesLapH/EARTH_RADIUS/EARTH_RADIUS
end subroutine

subroutine SWEVelocityActive(u, x, relVort, div, area, indexStart, indexEnd)
	real(kreal), intent(out) :: u(:,:)
	real(kreal), intent(in) :: x(:,:)
	real(kreal), intent(in) :: relVort(:), div(:), area(:)
	integer(kint), intent(in) :: indexStart, indexEnd
	!
	integer(kint) :: j, k, nn
	real(kreal) :: denom

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'entering SWEVelocityActive.')

	u = 0.0_kreal
	nn = size(x,2)
	do j=indexStart, indexEnd
		do k=1,j-1
			denom = EARTH_RADIUS*EARTH_RADIUS- sum(x(:,j)*x(:,k))
			u(1,j) = u(1,j) + ( x(2,j)*x(3,k) - x(3,j)*x(2,k) )*relVort(k)*area(k)/denom + &
							  ( x(1,j)*x(1,k)/EARTH_RADIUS/EARTH_RADIUS - x(1,k) )*div(k)*area(k)/(denom)
			u(2,j) = u(2,j) + ( x(3,j)*x(1,k) - x(1,j)*x(3,k) )*relVort(k)*area(k)/denom + &
							  ( x(2,j)*x(2,k)/EARTH_RADIUS/EARTH_RADIUS - x(2,k) )*div(k)*area(k)/(denom)
			u(3,j) = u(3,j) + ( x(1,j)*x(2,k) - x(2,j)*x(1,k) )*relVort(k)*area(k)/denom + &
					 		  ( x(3,j)*x(3,k)/EARTH_RADIUS/EARTH_RADIUS - x(3,k) )*div(k)*area(k)/(denom)
		enddo
		do k=j+1,nn
			denom = EARTH_RADIUS*EARTH_RADIUS- sum(x(:,j)*x(:,k))
			u(1,j) = u(1,j) + ( x(2,j)*x(3,k) - x(3,j)*x(2,k) )*relVort(k)*area(k)/denom + &
							  ( x(1,j)*x(1,k)/EARTH_RADIUS/EARTH_RADIUS - x(1,k) )*div(k)*area(k)/(denom)
			u(2,j) = u(2,j) + ( x(3,j)*x(1,k) - x(1,j)*x(3,k) )*relVort(k)*area(k)/denom + &
							  ( x(2,j)*x(2,k)/EARTH_RADIUS/EARTH_RADIUS - x(2,k) )*div(k)*area(k)/(denom)
			u(3,j) = u(3,j) + ( x(1,j)*x(2,k) - x(2,j)*x(1,k) )*relVort(k)*area(k)/denom + &
					 		  ( x(3,j)*x(3,k)/EARTH_RADIUS/EARTH_RADIUS - x(3,k) )*div(k)*area(k)/(denom)
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

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'entering SWEVelocityPassive.')
	u = 0.0_kreal
	nn = size(xActive,2)
	do j=indexStart, indexEnd
		do k=1,nn
			denom = EARTH_RADIUS*EARTH_RADIUS- sum(x(:,j)*xActive(:,k))
			u(1,j) = u(1,j) + ( x(2,j)*xActive(3,k) - x(3,j)*xActive(2,k) )*relVort(k)*area(k)/denom + &
							  ( x(1,j)*xActive(1,k)/EARTH_RADIUS/EARTH_RADIUS - xActive(1,k) )*div(k)*area(k)/(denom)
			u(2,j) = u(2,j) + ( x(3,j)*xActive(1,k) - x(1,j)*xActive(3,k) )*relVort(k)*area(k)/denom + &
							  ( x(2,j)*xActive(2,k)/EARTH_RADIUS/EARTH_RADIUS - xActive(2,k) )*div(k)*area(k)/(denom)
			u(3,j) = u(3,j) + ( x(1,j)*xActive(2,k) - x(2,j)*xActive(1,k) )*relVort(k)*area(k)/denom + &
					 		  ( x(3,j)*xActive(3,k)/EARTH_RADIUS/EARTH_RADIUS - xActive(3,k) )*div(k)*area(k)/(denom )
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
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'entering SWEVelocitySmooth.')
	u = 0.0_kreal
	nn = size(xActive,2)
	do j=indexStart, indexEnd
		do k=1,nn
			denom = EARTH_RADIUS*EARTH_RADIUS- sum(x(:,j)*xActive(:,k)) + VELOCITY_SMOOTH*VELOCITY_SMOOTH
			u(1,j) = u(1,j) + ( x(2,j)*xActive(3,k) - x(3,j)*xActive(2,k) )*relVort(k)*area(k)/denom + &
							  ( x(1,j)*xActive(1,k)/EARTH_RADIUS/EARTH_RADIUS - xActive(1,k) )*div(k)*area(k)/denom
			u(2,j) = u(2,j) + ( x(3,j)*xActive(1,k) - x(1,j)*xActive(3,k) )*relVort(k)*area(k)/denom + &
							  ( x(2,j)*xActive(2,k)/EARTH_RADIUS/EARTH_RADIUS - xActive(2,k) )*div(k)*area(k)/denom
			u(3,j) = u(3,j) + ( x(1,j)*xActive(2,k) - x(2,j)*xActive(1,k) )*relVort(k)*area(k)/denom + &
					 		  ( x(3,j)*xActive(3,k)/EARTH_RADIUS/EARTH_RADIUS - xActive(3,k) )*div(k)*area(k)/denom
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
