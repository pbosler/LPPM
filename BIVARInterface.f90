module BIVARInterfaceModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
!----------------
! USAGE :  This module provides methods for defining a mesh in the plane.
!----------------
use NumberKindsModule
use LoggerModule
use ParticlesModule
use PanelsModule
use PlaneMeshModule
use BIVARModule

implicit none

private
public BIVARSetup
public New, Delete
public AssignTracer, AssignVorticity

!
!----------------
! Types and module constants
!----------------
!
type BIVARSetup
	real(kreal), pointer :: x(:) => null()
	real(kreal), pointer :: y(:) => null()
	real(kreal), pointer :: x0(:) => null()
	real(kreal), pointer :: y0(:) => null()
	real(kreal), pointer :: vort(:) => null()
	real(kreal), pointer :: tracer(:,:) => null()
	real(kreal), pointer :: realWork(:) => null()
	
	integer(kint), pointer :: activeMap(:) => null()
	type(Particles), pointer :: particles => null()
	type(Panels), pointer :: activePanels => null()
	integer(kint) :: n	
	integer(kint) :: nParticles
	integer(kint) :: nActive
end type

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'BIVAR : '
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

interface ResetLagrangianParameter
	module procedure ResetAlphaPrivate
end interface

contains
!
!----------------
! Standard methods : Constructor / Destructor
!----------------
!
subroutine NewPrivate(self, aMesh)
	type(BIVARSetup), intent(out) :: self
	type(PlaneMesh), intent(in) :: aMesh
	!
	integer(kint) :: j, k, nActive, nPassive, nTracer, n
	type(Panels), pointer :: aPanels, passivePanels
	integer(kint), allocatable :: passiveMap(:)
	
	if (.NOT. logInit) call InitLogger(log, procRank)
	
	aPanels => aMesh%panels
	nActive = aPanels%N_Active
	nPassive = aPanels%N - aPanels%N_Active
	nTracer = aMesh%nTracer
	
	self%particles => aMesh%particles
	
	self%n = nActive + self%particles%N
	n = self%n
	self%nParticles = self%particles%N
	self%nActive = nActive
	
	allocate(self%activePanels)
	call New(self%activePanels, nActive, QUAD_PANEL, nTracer, PLANE_SOLVER)
	self%activePanels%n = nActive
	self%activePanels%N_Active = nActive
	allocate(self%activeMap(nActive))
	self%activeMap = 0
	
	allocate(passivePanels)
	call New(passivePanels, nPassive, QUAD_PANEL, nTracer, PLANE_SOLVER)
	passivePanels%N = nPassive
	allocate(passiveMap(nPassive))
	passiveMap = 0
	
	call GatherPanels(aPanels, self%activePanels, self%activeMap, passivePanels, passiveMap)
	
	allocate(self%x(n))
	self%x = 0.0_kreal
	allocate(self%y(n))
	self%y = 0.0_kreal
	allocate(self%x0(n))
	self%x0 = 0.0_kreal
	allocate(self%y0(n))
	self%y0 = 0.0_kreal
	allocate(self%vort(n))
	self%vort = 0.0_kreal
	if ( ntracer > 0 ) then
		allocate(self%tracer(n,nTracer))
		self%tracer = 0.0_kreal
	endif
	allocate(self%realWork( 8 * n ) )
	self%realWork = 0.0_kreal
	
	!
	! populate interpolation source data
	!
	do j = 1, self%particles%N
		self%x(j) = self%particles%x(1,j)
		self%y(j) = self%particles%x(2,j)
		self%x0(j) = self%particles%x0(1,j)
		self%y0(j) = self%particles%x0(2,j)
		self%vort(j) = self%particles%relVort(j)
	enddo
	do k = 1, nTracer
		do j = 1, self%particles%N
			self%tracer(j, k) = self%particles%tracer(j,k)
		enddo
	enddo
	
	do j = 1, self%activePanels%N
		self%x(self%nParticles + j) = self%activePanels%x(1,j)
		self%y(self%nParticles + j) = self%activePanels%x(2,j)
		self%x0(self%nParticles + j) = self%activePanels%x0(1,j)
		self%y0(self%nParticles + j) = self%activePanels%x0(2,j)
		self%vort(self%nParticles+j) = self%activePanels%relVort(j)
	enddo
	do k = 1, nTracer
		do j = 1, self%nActive
			self%tracer(self%nParticles+j, k) = self%activePanels%tracer(j,k)
		enddo
	enddo
	
	call Delete(passivePanels)
	deallocate(passivePanels)
	deallocate(passiveMap)
end subroutine


subroutine DeletePrivate(self)
	type(BIVARSetup), intent(inout) :: self
	
	call Delete(self%activePanels)
	deallocate(self%activeMap)
	deallocate(self%activePanels)
	
	self%nParticles = 0
	nullify(self%particles)
	
	deallocate(self%x)
	deallocate(self%y) 
	deallocate(self%x0)
	deallocate(self%y0)
	deallocate(self%vort)
	if ( associated(self%tracer)) deallocate(self%tracer)
end subroutine	

!
!----------------
! Public functions
!----------------
!
subroutine AssignTracer(tracerOut, alphaIn, self, tracerID)
	real(kreal), intent(out) :: tracerOut(:)
	real(kreal), intent(in) :: alphaIn(:,:)
	type(BIVARSetup), intent(in) :: self
	integer(kint), intent(in) :: tracerID
	!
	integer(kint) :: nOut, md
	integer(kint), allocatable :: integerWork(:)
	
	
	nOut = size(alphaIn,2)
	if ( nOut /= size(tracerOut) ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL,trim(logkey)//' AssignTracer ERROR : ', ' size mismatch.')
		return
	endif
	md = 1
	allocate(integerWork( 31*self%n + nOut))
	integerWork = 0
	call IDBVIP( md, self%n, self%x, self%y, self%tracer(:,tracerId), &
				 nOut, alphaIn(1,:), alphaIn(2,:), tracerOut, integerWork, self%realWork)
	deallocate(integerWork)
end subroutine

subroutine AssignVorticity(vorticityOut, alphaIn, self)
	real(kreal), intent(out) :: vorticityOut(:)
	real(kreal), intent(in) :: alphaIn(:,:)
	type(BIVARSetup), intent(in) :: self
	!
	integer(kint) :: nOut, md
	integer(kint), allocatable :: integerWork(:)
	
	nOut = size(alphaIn,2)
	if ( nOut /= size(vorticityOut) ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL,trim(logkey)//' AssignVorticity ERROR : ', ' size mismatch.')
		return
	endif
	md = 1
	allocate(integerWork( 31*self%n + nOut))
	integerWork = 0
	call IDBVIP(md, self%n, self%x, self%y, self%vort, &
				nOut, alphaIn(1,:), alphaIn(2,:), vorticityOut, integerWork, self%realWork)
	deallocate(integerWork)
end subroutine

subroutine ResetAlphaPrivate(self)
	type(BIVARSetup), intent(inout) :: Self
	!
	integer(kint) :: j
	do j = 1, self%n
		self%x0(j) = self%x(j)
		self%y0(j) = self%y(j)
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
