module PanelsModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method - Spherical Model
!------------------------------------------------------------------------------
!
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
!> @defgroup Panels Panels module
!> Provides the primitive Panels data structure that defines the passive particles of LPPM meshes.
!
!
! DESCRIPTION:
!> @file
!> Provides the primitive Panels data structure that defines the passive particles of LPPM meshes.
!
!------------------------------------------------------------------------------
use NumberKindsModule
use LoggerModule

implicit none
private
public Panels
public New, Delete, Copy
public GetNTracer, GetPanelKind
public GatherPanels, ScatterPanels
public PanelMax
public LogStats

!
!----------------
! Types and module constants
!----------------
!
type Panels
	! Grid variables
	real(kreal), pointer :: x(:,:) => null()	! physical position
	real(kreal), pointer :: x0(:,:) => null()	! Lagrangian coordinate
	real(kreal), pointer :: area(:) => null()	! panel area
	integer(kint), pointer :: edges(:,:) => null() ! edges around panel (pointers to Edges)
	integer(kint), pointer :: vertices(:,:) => null() ! particles at panel vertices (pointers to Particles)
	integer(kint), pointer :: children(:,:) => null() ! child panels (pointers to Panels)
	integer(kint), pointer :: nest(:) => null()	! nest level of panels
	logical(klog), pointer :: hasChildren(:) => null()	! true if panel has been divided
	integer(kint) :: N	! number of panels in memory
	integer(kint) :: N_Active ! number of panels in computation
	integer(kint) :: N_Max	! max number of panels allowed in memory
	! Data Variables
	real(kreal), pointer :: tracer(:,:) => null()	! passive tracers
	real(kreal), pointer :: absVort(:) => null()	! absolute vorticity
	real(kreal), pointer :: relVort(:) => null()	! relative vorticity
	real(kreal), pointer :: stream(:) => null() 	! stream function
	real(kreal), pointer :: potVort(:) => null()	! potential vorticity
	real(kreal), pointer :: h(:) => null()			! fluid thickness
	real(kreal), pointer :: div(:) => null()		! divergence
	real(kreal), pointer :: u(:,:) => null()		! fluid velocity
	real(kreal), pointer :: ke(:) => null()			! local kinetic pe
end type

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .False.
type(Logger) :: log
character(len=28), save :: logKey = 'Panels'
integer(kint), parameter :: logLevel = TRACE_LOGGING_LEVEL
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

interface Copy
	module procedure CopyPanelByIndex
	module procedure CopyPanels
end interface

interface GetNTracer
	module procedure GetNTracerPanels
end interface

interface LogStats
	module procedure LogPanelStats
end interface

contains
!
!----------------
! Standard methods : Constructor / Destructor
!----------------
!
subroutine NewPrivate(self,nMax,panelKind,nTracer,problemKind)
! allocates memory for a panels object, initializes state to zero / null.
	! calling parameters
	type(Panels), intent(out) :: self
	integer(kint), intent(in) :: nMax, &
								 panelKind, &
								 nTracer, &
								 problemKind
	if ( .NOT. logInit) then
		call InitLogger(log,procRank)
	endif

	! Error checking
	if ( nMax <= 0 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey," New ERROR: invalid nMax.")
		return
	endif
	if ( panelKind /= TRI_PANEL .AND. panelKind /= QUAD_PANEL ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey," New ERROR: invalid panelKind.")
		return
	endif
	if ( problemKind < ADVECTION_SOLVER .AND. problemKind > SWE_SOLVER ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey," New ERROR: invalid problemKind.")
		return
	endif
	!
	! Allocate data structure, set to zero/null
	!
	if ( problemKind == PLANE_SOLVER ) then
		allocate(self%x(2,nMax))
		allocate(self%x0(2,nMax))
		allocate(self%u(2,nMax))
	else
		allocate(self%x(3,nMax))
		allocate(self%x0(3,nMax))
		allocate(self%u(3,nMax))
	endif
	self%x = 0.0_kreal
	self%x0 = 0.0_kreal
	self%u = 0.0_kreal

	allocate(self%area(nMax))
	self%area = 0.0_kreal

	if ( panelKind == TRI_PANEL ) then
		allocate(self%edges(3,nMax))
		allocate(self%vertices(3,nMax))
	elseif ( panelKind == QUAD_PANEL) then
		allocate(self%edges(4,nMax))
		allocate(self%vertices(4,nMax))
	endif
	self%edges = 0
	self%vertices = 0
	allocate(self%nest(nMax))
	self%nest = -1
	allocate(self%children(4,nMax))
	self%children = 0
	allocate(self%hasChildren(nMax))
	self%hasChildren = .FALSE.
	if ( nTracer > 0 ) then
		allocate(self%tracer(nMax,nTracer))
		self%tracer = 0.0_kreal
	else
		nullify(self%tracer)
	endif
	if ( problemKind == ADVECTION_SOLVER) then
		! nullify bve / swe variables
		nullify(self%absVort)
		nullify(self%relVort)
		nullify(self%potVort)
		nullify(self%h)
		nullify(self%div)
		nullify(self%ke)
		nullify(self%stream)
	elseif (problemKind == BVE_SOLVER .OR. problemKind == PLANE_SOLVER ) then
		! allocate bve variables
		allocate(self%relVort(nMax))
		self%relVort = 0.0_kreal
		allocate(self%absVort(nMax))
		self%absVort = 0.0_kreal
		allocate(self%ke(nMax))
		self%ke = 0.0_kreal
		allocate(self%stream(nMax))
		self%stream = 0.0_kreal
		! nullify swe variables
		nullify(self%h)
		nullify(self%div)
		nullify(self%potVort)
	elseif (problemKind == SWE_SOLVER) then
		! allocate swe variables
		allocate(self%relVort(nMax))
		self%relVort = 0.0_kreal
		allocate(self%potVort(nMax))
		self%potVort = 0.0_kreal
		allocate(self%h(nMax))
		self%h = 0.0_kreal
		allocate(self%div(nMax))
		self%div = 0.0_kreal
		allocate(self%ke(nMax))
		self%ke = 0.0_kreal
		! nullify bve variables
		nullify(self%absVort)
		nullify(self%stream)
	endif
	self%N = 0
	self%N_Active = 0
	self%N_Max = nMax
end subroutine

subroutine DeletePrivate(self)
! frees memory allocated for a panels object
	type(Panels), intent(inout) :: self
	deallocate(self%x)
	deallocate(self%x0)
	deallocate(self%area)
	deallocate(self%edges)
	deallocate(self%vertices)
	deallocate(self%children)
	deallocate(self%hasChildren)
	deallocate(self%nest)
	deallocate(self%u)
	if ( associated(self%tracer)) deallocate(self%tracer)
	if ( associated(self%absVort)) deallocate(self%absVort)
	if ( associated(self%relVort)) deallocate(self%relVort)
	if ( associated(self%stream)) deallocate(self%stream)
	if ( associated(self%potVort)) deallocate(self%potVort)
	if ( associated(self%h)) deallocate(self%h)
	if ( associated(self%div)) deallocate(self%div)
	if ( associated(self%ke)) deallocate(self%ke)
	self%N = 0
	self%N_Active = 0
	self%N_Max = 0
end subroutine

subroutine CopyPanels(newPanels,oldPanels)
! copies an entire panels object
	type(Panels), intent(inout) :: newPanels
	type(Panels), intent(in) :: oldPanels
	integer(kint) :: j

	! Error checking
	if ( newPanels%N_Max < oldPanels%N) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'Copy ERROR : not enough memory.')
		return
	endif
	if ( associated(oldPanels%tracer)) then
		if ( .NOT. associated(newPanels%tracer)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'Copy ERROR : cannot assign tracer.')
			return
		endif
	endif
	if ( associated(oldPanels%absVort)) then
		if ( .NOT. associated(newPanels%absVort)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'Copy ERROR : cannot assign absVort.')
			return
		endif
	endif
	if ( associated(oldPanels%relVort)) then
		if ( .NOT. associated(oldPanels%relVort)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'Copy ERROR : cannot assign relVort.')
			return
		endif
	endif
	if ( associated(oldPanels%stream)) then
		if ( .NOT. associated(oldPanels%stream)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'Copy ERROR : cannot assign stream.')
			return
		endif
	endif
	if ( size(newPanels%edges,1) /= size(oldPanels%edges,1)) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'Copy ERROR : mismatched panelKind.')
		return
	endif
	if ( associated(oldPanels%potVort) ) then
		if (.NOT. associated(newPanels%potVort)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticles ERROR : cannot assign potVort.')
			return
		endif
	endif
	if ( associated(oldPanels%h) ) then
		if (.NOT. associated(newPanels%h)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticles ERROR : cannot assign h.')
			return
		endif
	endif
	if ( associated(oldPanels%div) ) then
		if (.NOT. associated(newPanels%div)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticles ERROR : cannot assign div.')
			return
		endif
	endif
!	if ( associated(oldPanels%pe) ) then
!		if (.NOT. associated(newPanels%pe)) then
!			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticles ERROR : cannot assign pe.')
!			return
!		endif
!	endif
	if ( associated(oldPanels%ke) ) then
		if (.NOT. associated(newPanels%ke)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticles ERROR : cannot assign ke.')
			return
		endif
	endif


	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering CopyParticles.')

	do j=1,oldPanels%N
		newPanels%x(:,j) = oldPanels%x(:,j)
		newPanels%x0(:,j) = oldPanels%x0(:,j)
		newPanels%area(j) = oldPanels%area(j)
		newPanels%edges(:,j) = oldPanels%edges(:,j)
		newPanels%vertices(:,j) = oldPanels%vertices(:,j)
		newPanels%nest(j) = oldPanels%nest(j)
		newPanels%children(:,j) = oldPanels%children(:,j)
		newPanels%hasChildren(j) = oldPanels%hasChildren(j)
		newPanels%u(:,j) = oldPanels%u(:,j)
	enddo
	do j = oldPanels%N + 1, newPanels%N_Max
		newPanels%x(:,j) = 0.0_kreal
		newPanels%x0(:,j) = 0.0_kreal
		newPanels%area(j) = 0.0_kreal
		newPanels%edges(:,j) = 0
		newPanels%vertices(:,j) = 0
		newPanels%nest(j) = -1
		newPanels%children(:,j) = 0
		newPanels%hasChildren(j) = .FALSE.
		newPanels%u(:,j) = 0.0_kreal
	enddo
	if ( associated(oldPanels%tracer)) then
		do j=1,oldPanels%N
			newPanels%tracer(j,:) = oldPanels%tracer(j,:)
		enddo
		do j = oldPanels%N + 1, newPanels%N_Max
			newPanels%tracer(j,:) = 0.0_kreal
		enddo
	endif
	if ( associated(oldPanels%relVort)) then
		do j=1,oldPanels%N
			newPanels%relVort(j) = oldPanels%relVort(j)
		enddo
		do j = oldPanels%N + 1, newPanels%N_Max
			newPanels%relVort(j) = 0.0_kreal
		enddo
	endif
	if ( associated(oldPanels%stream)) then
		do j=1,oldPanels%N
			newPanels%stream(j) = oldPanels%stream(j)
		enddo
		do j = oldPanels%N + 1, newPanels%N_Max
			newPanels%stream(j) = 0.0_kreal
		enddo
	endif
	if ( associated(oldPanels%absVort)) then
		do j=1,oldPanels%N
			newPanels%absVort(j) = oldPanels%absVort(j)
		enddo
		do j = oldPanels%N + 1, newPanels%N_Max
			newPanels%absVort(j) = 0.0_kreal
		enddo
	endif
	if ( associated(oldPanels%potVort)) then
		do j=1,oldPanels%N
			newPanels%potVort(j) = oldPanels%potVort(j)
		enddo
		do j = oldPanels%N + 1, newPanels%N_Max
			newPanels%potVort(j) = 0.0_kreal
		enddo
	endif
	if ( associated(oldPanels%h)) then
		do j=1,oldPanels%N
			newPanels%h(j) = oldPanels%h(j)
		enddo
		do j = oldPanels%N + 1, newPanels%N_Max
			newPanels%h(j) = 0.0_kreal
		enddo
	endif
	if ( associated(oldPanels%div)) then
		do j=1,oldPanels%N
			newPanels%div(j) = oldPanels%div(j)
		enddo
		do j = oldPanels%N + 1, newPanels%N_Max
			newPanels%div(j) = 0.0_kreal
		enddo
	endif
	if ( associated(oldPanels%ke)) then
		do j=1,oldPanels%N
			newPanels%ke(j) = oldPanels%ke(j)
		enddo
		do j = oldPanels%N + 1, newPanels%N_Max
			newPanels%ke(j) = 0.0_kreal
		enddo
	endif
!	if ( associated(oldPanels%pe)) then
!		do j=1,oldPanels%N
!			newPanels%pe(j) = oldPanels%pe(j)
!		enddo
!	endif
	newPanels%N = oldPanels%N
	newPanels%N_Active = oldPanels%N_Active
end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine CopyPanelByIndex(newPanels,newIndex,oldPanels,oldIndex)
	type(Panels), intent(inout) :: newPanels
	type(Panels), intent(in) :: oldPanels
	integer(kint), intent(in) :: newIndex, oldIndex

	! Error checking
	if ( newIndex > newPanels%N_Max .OR. oldIndex > oldPanels%N ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyPanelByIndex ERROR : out of bounds.')
		return
	endif

	newPanels%x(:,newIndex) = oldPanels%x(:,oldIndex)
	newPanels%x0(:,newIndex) = oldPanels%x0(:,oldIndex)
	newPanels%area(newIndex) = oldPanels%area(oldIndex)
	newPanels%edges(:,newIndex) = oldPanels%edges(:,oldIndex)
	newPanels%vertices(:,newIndex) = oldPanels%vertices(:,oldIndex)
	newPanels%nest(newIndex) = oldPanels%nest(oldIndex)
	newPanels%children(:,newINdex) = oldPanels%children(:,oldIndex)
	newPanels%hasChildren(newIndex) = oldPanels%hasChildren(oldIndex)
	newPanels%u(:,newIndex) = oldPanels%u(:,oldINdex)
	if ( associated(oldPanels%tracer) ) then
		if ( associated(newPanels%tracer) ) then
			newPanels%tracer(newIndex,:) = oldPanels%tracer(oldIndex,:)
		else
			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyPanelByIndex WARNING : cannot assign tracer.')
		endif
	endif
	if ( associated(oldPanels%relVort)) then
		if ( associated(newPanels%relVort)) then
			newPanels%relVort(newIndex) = oldPanels%relVort(oldIndex)
		else
			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyPanelByIndex WARNING: cannot assign relVort.')
		endif
	endif
	if ( associated(oldPanels%stream)) then
		if ( associated(newPanels%stream)) then
			newPanels%stream(newIndex) = oldPanels%stream(oldIndex)
		else
			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyPanelByIndex WARNING: cannot assign stream.')
		endif
	endif
	if ( associated(oldPanels%absVort)) then
		if ( associated(newPanels%absVort)) then
			newPanels%absVort(newIndex) = oldPanels%absVort(oldIndex)
		else
			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyPanelByIndex WARNING: cannot assign absVort.')
		endif
	endif
	if ( associated(oldPanels%potVort)) then
		if ( associated(newPanels%potVort)) then
			newPanels%potVort(newIndex) = oldPanels%potVort(oldIndex)
		else
			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyPanelByIndex WARNING: cannot assign potVort.')
		endif
	endif
	if ( associated(oldPanels%h)) then
		if ( associated(newPanels%h)) then
			newPanels%h(newIndex) = oldPanels%h(oldIndex)
		else
			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyPanelByIndex WARNING: cannot assign h.')
		endif
	endif
	if ( associated(oldPanels%div)) then
		if ( associated(newPanels%div)) then
			newPanels%div(newIndex) = oldPanels%div(oldIndex)
		else
			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyPanelByIndex WARNING: cannot assign div.')
		endif
	endif
	if ( associated(oldPanels%ke)) then
		if ( associated(newPanels%ke)) then
			newPanels%ke(newIndex) = oldPanels%ke(oldIndex)
		else
			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyPanelByIndex WARNING: cannot assign ke.')
		endif
	endif
!	if ( associated(oldPanels%pe)) then
!		if ( associated(newPanels%pe)) then
!			newPanels%pe(newIndex) = oldPanels%pe(oldIndex)
!		else
!			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyPanelByIndex WARNING: cannot assign pe.')
!		endif
!	endif
end subroutine

function GetPanelKind(self)
! returns the integer identifier of panel kind (TRI_PANEL = 3, QUAD_PANEL = 4)
	type(Panels), intent(in) :: self
	integer(kint) :: GetPanelKind
	GetPanelKind = size(self%vertices,1)
end function

function GetNTracerPanels(self)
! returns the number of tracers carried by a panel set
	type(Panels), intent(in) :: self
	integer(kint) :: GetNTracerPanels
	if ( associated(self%tracer)) then
		GetNTracerPanels = size(self%tracer,2)
	else
		GetNTracerPanels = 0
	endif
end function


function GetProblemKind(self)
! returns the problem kind of a panel set
	type(Panels), intent(in) :: self
	integer(kint) :: GetProblemKind
	if ( associated(self%relVort) .AND. associated(self%absVort)) then
		GetProblemKind = BVE_SOLVER
	elseif ( associated(self%potVort) .AND. associated(self%h) ) then
		GetProblemKind = SWE_SOLVER
	else
		GetProblemKind = ADVECTION_SOLVER
	endif
end function

subroutine GatherPanels(self,activePanels,activeMap,passivePanels,passiveMap)
! separates low-level panels from subdivided panels.  stores indices of separated panels in activeMap and passiveMap arrays.
	! calling paramters
	type(Panels), intent(in) :: self
	type(Panels), intent(inout) :: activePanels, &
								 passivePanels
	integer(kint), intent(inout) :: activeMap(:), &
								  passiveMap(:)
	! local variables
	integer(kint) :: j, activeK, passiveK, nActive, nPassive

	nActive = self%N_Active
	nPassive = self%N - self%N_Active

	! Error checking
	if ( size(activeMap) /= nActive ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'GatherPanels ERROR : activeMap size mismatch.')
		return
	endif
	if ( size(passiveMap) /= nPassive) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'GatherPanels ERROR : passiveMap size mismatch.')
		return
	endif
	if ( activePanels%N_Max < nActive ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'GatherPanels ERROR : not enough memory for active panels.')
		call LogMessage(log,ERROR_LOGGING_LEVEL,'activePanels%N_Max = ',activePanels%N_Max)
		call LogMessage(log,ERROR_LOGGING_LEVEL,'nActive (whole mesh) = ',nActive)
		return
	endif
	if ( passivePanels%N_Max < nPassive) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'GatherPanels ERROR : not enough memory for passive panels.')
		return
	endif

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'Entering GatherPanels...')

	activeK = 1
	passiveK = 1
	do j=1,self%N
		if ( self%hasChildren(j) ) then ! panel j is passive
			call CopyPanelByIndex(passivePanels,passiveK,self,j)
			passiveMap(passiveK) = j
			passiveK = passiveK + 1
		else ! panel j is active
			call CopyPanelByIndex(activePanels,activeK,self,j)
			activeMap(activeK) = j
			activeK = activeK + 1
		endif
	enddo

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'GatherPanels complete.')
end subroutine


subroutine ScatterPanels(self,activePanels,activeMap,passivePanels,passiveMap)
! copies separated active / passive panels into a whole panels object. (inverts the GatherPanels subroutine)
	! calling parameters
	type(Panels), intent(inout) :: self
	type(Panels), intent(in) :: activePanels, &
								passivePanels
	integer(kint), intent(in) :: activeMap(:), &
								 passiveMap(:)
	! local variables
	integer(kint) :: j

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering ScatterPanels...')

	do j=1,activePanels%N
		call CopyPanelByIndex(self,activeMap(j),activePanels,j)
	enddo
	do j=1,passivePanels%N
		call CopyPanelByIndex(self,passiveMap(j),passivePanels,j)
	enddo

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'ScatterPanels complete.')
end subroutine


function PanelMax(panelKind,maxNest)
! computes the maximum number of panels for preallocation
	integer(kint) :: PanelMax
	integer(kint), intent(in) :: panelKind, &
								 maxNest
	integer(kint) :: j
	! Error checking
	if ( maxNest < 0 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'PanelMax ERROR : invalid maxNest.')
		return
	endif
	PanelMax = 0
	if ( panelKind == TRI_PANEL ) then
		do j=0,maxNest
			PanelMax = PanelMax + 20*4**j
		enddo
	elseif ( panelKind == QUAD_PANEL ) then
		do j=0,maxNest
			PanelMax = PanelMax + 6*4**j
		enddo
	else
		PanelMax = -1
	endif
end function

!
!----------------
! Module methods : type-specific functions
!----------------
!

subroutine LogPanelStats(self,aLog,message)
! outputs data associated with a panels object to a Logger object
	type(Panels), intent(in) :: self
	type(Logger), intent(inout) :: aLog
	character(len=*), intent(in), optional :: message
	! local variables
	character(len=32) :: key
	integer(kint) :: k, nTracer, j
	real(kreal) :: maxU, minU, uj

	if (present(message) ) then
		if ( GetPanelKind(self) == TRI_PANEL) then
			call StartSection(aLog,'TriPanel Stats : ',message)
		elseif ( GetPanelKind(self) == QUAD_PANEL ) then
			call StartSection(aLog,'QuadPanel Stats : ',message)
		endif
	else
		if ( GetPanelKind(self) == TRI_PANEL) then
			call StartSection(aLog,'TriPanel Stats : ')
		elseif ( GetPanelKind(self) == QUAD_PANEL ) then
			call StartSection(aLog,'QuadPanel Stats : ')
		endif
	endif
	key = 'N_Active = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,self%N_Active)
	key = 'N = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,self%N)
	key = 'N_Max = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,self%N_Max)
	key = 'MaxNest = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,maxval(self%nest(1:self%N)))
	key = 'Surface Area error <m^2> = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,sum(self%area(1:self%N)))
	maxU = 0.0_kreal
	minU = 0.0_kreal
	do j=1,self%N
		 uj = sqrt(sum( self%u(:,j)*self%u(:,j)))
		 if ( uj > maxU ) maxU = uj
		 if ( uj < minU ) minU = uj
	enddo
	key = 'max. velocity <m s^(-1)> = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,maxU)
	key = 'min. velocity <m s^(-1)> = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,minU)
	if ( associated(self%absVort)) then
		key = 'Max absVort <s^(-1)> = '
		call LogMessage(alog,TRACE_LOGGING_LEVEL,key,maxVal(self%absVort(1:self%N)))
		key = 'Min absVort <s^(-1)> = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,minVal(self%absVort(1:self%N)))
		key = 'absVort integral = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,sum(self%absVort(1:self%N)*self%area(1:self%N))&
				/(4.0_kreal*PI*EARTH_RADIUS*EARTH_RADIUS))
	endif
	if ( associated(self%relVort)) then
		key = 'Max relVort <s^(-1)> = '
		call LogMessage(alog,TRACE_LOGGING_LEVEL,key,maxVal(self%relVort(1:self%N)))
		key = 'Min relVort <s^(-1)> = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,minVal(self%relVort(1:self%N)))
		key = 'normalized relVort integral = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,sum(self%relVort(1:self%N)*self%area(1:self%N))&
			/(4.0_kreal*PI*EARTH_RADIUS*EARTH_RADIUS))
	endif
	if ( associated(self%potVort)) then
		key = 'Max potVort < m^(-1) s^(-1)> = '
		call LogMessage(alog,TRACE_LOGGING_LEVEL,key,maxVal(self%potVort(1:self%N)))
		key = 'Min potVort = <m^(-1) s^(-1)> '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,minVal(self%potVort(1:self%N)))
		key = 'normalized potVort integral = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,sum(self%potVort(1:self%N)*self%area(1:self%N)) &
				/(4.0_kreal*PI*EARTH_RADIUS*EARTH_RADIUS))
	endif
	if ( associated(self%h)) then
		key = 'Max h <m> = '
		call LogMessage(alog,TRACE_LOGGING_LEVEL,key,maxVal(self%h(1:self%N)))
		key = 'Min h <m> = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,minVal(self%h(1:self%N)))
		key = 'normalized h integral = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,sum(self%h(1:self%N)*self%area(1:self%N)) &
			/(4.0_kreal*PI*EARTH_RADIUS*EARTH_RADIUS))
	endif
	if ( associated(self%div)) then
		key = 'Max div = '
		call LogMessage(alog,TRACE_LOGGING_LEVEL,key,maxVal(self%div(1:self%N)))
		key = 'Min div = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,minVal(self%div(1:self%N)))
		key = 'normalized div integral = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,sum(self%div(1:self%N)*self%area(1:self%N)) & 
			/(4.0_kreal*PI*EARTH_RADIUS*EARTH_RADIUS))
	endif
	if ( associated(self%tracer)) then
		nTracer = GetNTracer(self)
		do k=1,nTracer
			write(key,'(A,I2,A)') 'Max tracer',k,' = '
			call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,maxVal(self%tracer(1:self%N,k)))
			write(key,'(A,I2,A)') 'Min tracer',k,' = '
			call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,minVal(self%tracer(1:self%N,k)))
			write(key,'(A,I2,A)') 'normalized tracer',k,' integral = '
			call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,sum(self%tracer(1:self%N,k)*self%area(1:self%N)) &
				/(4.0_kreal*PI*EARTH_RADIUS*EARTH_RADIUS))
		enddo
	endif

	call EndSection(aLog)
end subroutine


subroutine InitLogger(aLog,rank)
! initializes a logger for this module and processor
	type(Logger), intent(inout) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

end module
