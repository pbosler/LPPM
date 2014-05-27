module PlaneDirectSumModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	Defines the RK4 data structure and methods used by PlaneMesh.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
!----------------
! USAGE :  This module provides methods for integrating the incompressible Euler equations in the plane.
!----------------
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
public RK4Timestep
public VELOCITY_SMOOTH

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
