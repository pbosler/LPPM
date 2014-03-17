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
	type(SSRFPACKData), pointer :: velocitySource => null(), &
								   thicknessSource => null()

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
							activePanelsInput(:,:) => null(), &
							activePanelsStage1(:,:) => null(), &
							activePanelsStage2(:,:) => null(), &
							activePanelsStage3(:,:) => null(), &
							activePanelsStage4(:,:) => null(), &
							newActivePanelsX(:,:) => null(), &
							&
							passivePanelsInput(:,:) => null(), &
							passivePanelsStage1(:,:) => null(), &
							passivePanelsStage2(:,:) => null(), &
							passivePanelsStage3(:,:) => null(), &
							passivePanelsStage4(:,:) => null(), &
							newPassivePanelsX(:,:) => null(), &
							&
							activeVelocityInput(:,:) => null(), &
							activeVelocityStage1(:,:) => null(), &
							activeVelocityStage2(:,:) => null(), &
							activeVelocityStage3(:,:) => null(), &
							activeVelocityStage4(:,:) => null(), &
							newActiveVelocity(:,:) => null(), &
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
							newActiveArea(:) => null()
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

	allocate(self%velocitySource)
	call New(self%velocitySource,self%delTri,.TRUE.)

	allocate(self%thicknessSource)
	call New(self%thicknessSource,self%delTri,.FALSE.)

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
!
!----------------
! Public member functions
!----------------
!


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
