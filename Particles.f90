module ParticlesModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	Defines the passive particle data structure used by icosahedral triangle and cubed
!	sphere Lagrangian meshes of the sphere.
!
!	No Dual : no information about the dual mesh is provided by this implementation.
!
!	No lists : Dynamic memory allocation is avoided by assuming constant polygonal structures
!		 (i.e., quadrilaterals or triangles only), and not maintaining dual mesh methods.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
use NumberKindsModule
use LoggerModule

implicit none
private
public Particles
public New, Delete, Copy
public GetNTracer
public ParticleMax
public LogStats
!
!----------------
! Types and module constants
!----------------
!
type Particles
	! Grid variables
	real(kreal), pointer :: x(:,:)	=> null() ! physical coordinate
	real(kreal), pointer :: x0(:,:) => null() ! Lagrangian coordinate
	integer(kint) :: N				! N particles in computation
	integer(kint) :: N_Max			! Max particles allowed in memory
	! Data variables
	real(kreal), pointer :: tracer(:,:)	=> null() ! passive tracers
	real(kreal), pointer :: absVort(:)	=> null() ! absolute vorticity
	real(kreal), pointer :: relVort(:)	=> null() ! relative vorticity
	real(kreal), pointer :: potVort(:) => null()  ! potential vorticity
	real(kreal), pointer :: h(:) => null()  ! fluid thickness
	real(kreal), pointer :: div(:) => null() ! divergence
	real(kreal), pointer :: u(:,:) => null() ! fluid velocity
	real(kreal), pointer :: ke(:) => null()  ! local kinetic pe
	!real(kreal), pointer :: pe(:) => null() ! total pe of each particle
end type

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'Particles'
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
	module procedure CopyParticleByIndex
	module procedure CopyParticles
end interface

interface GetNTracer
	module procedure GetNTracerParticles
end interface

interface LogStats
	module procedure LogStatsParticles
end interface

contains
!
!----------------
! Standard methods : Constructor / Destructor, Copy
!----------------
!
subroutine NewPrivate(self,nMax,panelKind,nTracer,problemKind)
! Allocates memory for a particles data structure. Sets initial state to zero / null.
	! Calling parameters
	type(Particles), intent(out) :: self
	integer(kint), intent(in) :: nMax, &
								 panelKind, &
								 nTracer, &
								 problemKind
	if (.NOT. logInit ) then
		call InitLogger(log,procRank)
	endif

	! Error checking
	if ( nMax <= 0 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey," invalid nMax.")
		return
	endif
	if ( panelKind /= TRI_PANEL .AND. panelKind /= QUAD_PANEL ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey," invalid panelKind.")
		return
	endif
	if ( problemKind < ADVECTION_SOLVER .AND. problemKind > SWE_SOLVER) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey," invalid panelKind.")
		return
	endif
	!
	! Allocate data structure, set to zero/null
	!
	allocate(self%x(3,nMax))
	self%x = 0.0_kreal
	allocate(self%x0(3,nMax))
	self%x0 = 0.0_kreal
	allocate(self%u(3,nMax))
	self%u = 0.0_kreal

	if ( nTracer > 0 ) then
		allocate(self%tracer(nMax,nTracer))
		self%tracer = 0.0_kreal
	else
		nullify(self%tracer)
	endif
	if ( problemKind == ADVECTION_SOLVER ) then
		! nullify bve/swe variables
		nullify(self%absVort)
		nullify(self%relVort)
		nullify(self%potVort)
		nullify(self%h)
		nullify(self%div)
		nullify(self%ke)
	elseif (problemKind == BVE_SOLVER) then
		! allocate bve variables
		allocate(self%absVort(nMax))
		self%absVort = 0.0_kreal
		allocate(self%relVort(nMax))
		self%relVort = 0.0_kreal
		allocate(self%ke(nMax))
		self%ke = 0.0_kreal
		! nullify swe variables
		nullify(self%potVort)
		nullify(self%h)
		nullify(self%div)
	elseif (problemKind == SWE_SOLVER ) then
		! nullify bve variables
		nullify(self%absVort)
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
	endif
	self%N = 0
	self%N_Max = nMax
end subroutine


subroutine DeletePrivate(self)
! Frees memory associated with an instance of a particles data structure.
	type(Particles), intent(inout) :: self
	deallocate(self%x)
	deallocate(self%x0)
	deallocate(self%u)
	if ( associated(self%tracer)) deallocate(self%tracer)
	if ( associated(self%absVort)) deallocate(self%absVort)
	if ( associated(self%relVort)) deallocate(self%relVort)
	if ( associated(self%potVort)) deallocate(self%potVort)
	if ( associated(self%h)) deallocate(self%h)
	if ( associated(self%div)) deallocate(self%div)
	if ( associated(self%ke)) deallocate(self%ke)
	!if ( associated(self%pe)) deallocate(self%pe)
	self%N = 0
	self%N_Max = 0
end subroutine


subroutine CopyParticles(newParticles,oldParticles)
! Copies an entire particles data structure by replicating data.
	type(Particles), intent(inout) :: newParticles
	type(Particles), intent(in) :: oldParticles
	integer(kint) :: j

	! Error checking
	if ( newParticles%N_Max < oldParticles%N) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticles ERROR : not enough memory.')
		return
	endif
	if ( associated(oldParticles%tracer)) then
		if ( .NOT. associated(newParticles%tracer)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticles ERROR : cannot assign tracer.')
			return
		endif
	endif
	if (associated(oldParticles%absVort)) then
		if ( .NOT. associated(newParticles%absVort)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticles ERROR : cannot assign absVort.')
			return
		endif
	endif
	if ( associated(oldParticles%relVort) ) then
		if (.NOT. associated(newParticles%relVort)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticles ERROR : cannot assign relVort.')
			return
		endif
	endif
	if ( associated(oldParticles%potVort) ) then
		if (.NOT. associated(newParticles%potVort)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticles ERROR : cannot assign potVort.')
			return
		endif
	endif
	if ( associated(oldParticles%h) ) then
		if (.NOT. associated(newParticles%h)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticles ERROR : cannot assign h.')
			return
		endif
	endif
	if ( associated(oldParticles%div) ) then
		if (.NOT. associated(newParticles%div)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticles ERROR : cannot assign div.')
			return
		endif
	endif
	if ( associated(oldParticles%u) ) then
		if (.NOT. associated(newParticles%u)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticles ERROR : cannot assign u.')
			return
		endif
	endif
	if ( associated(oldParticles%ke) ) then
		if (.NOT. associated(newParticles%ke)) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticles ERROR : cannot assign ke.')
			return
		endif
	endif
!	if ( associated(oldParticles%pe) ) then
!		if ( .NOT. associated(newParticles%pe) ) then
!			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticles ERROR : cannot assign pe.')
!			return
!		endif
!	endif
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering CopyParticles.')

	do j=1,oldParticles%N
		newParticles%x(:,j) = oldParticles%x(:,j)
		newParticles%x0(:,j) = oldParticles%x0(:,j)
		newParticles%u(:,j) = oldParticles%u(:,j)
	enddo
	if ( associated(oldParticles%tracer)) then
		do j=1,oldParticles%N
			newParticles%tracer(j,:) = oldParticles%tracer(j,:)
		enddo
	endif
	if ( associated(oldParticles%relVort)) then
		do j=1,oldParticles%N
			newParticles%relVort(j) = oldParticles%relVort(j)
		enddo
	endif
	if ( associated(oldParticles%absVort) ) then
		do j=1,oldParticles%N
			newParticles%absVort(j) = oldParticles%absVort(j)
		enddo
	endif
	if ( associated(oldParticles%potVort) ) then
		do j=1,oldParticles%N
			newParticles%potVort(j) = oldParticles%potVort(j)
		enddo
	endif
	if ( associated(oldParticles%h) ) then
		do j=1,oldParticles%N
			newParticles%h(j) = oldParticles%h(j)
		enddo
	endif
	if ( associated(oldParticles%div) ) then
		do j=1,oldParticles%N
			newParticles%div(j) = oldParticles%div(j)
		enddo
	endif
	if ( associated(oldParticles%ke) ) then
		do j=1,oldParticles%N
			newParticles%ke(j) = oldParticles%ke(j)
		enddo
	endif
!	if ( associated(oldParticles%pe)) then
!		do j=1,oldParticles%N
!			newParticles%pe(j) = oldParticles%pe(j)
!		enddo
!	endif
	newParticles%N = oldParticles%N
end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine CopyParticleByIndex(newParticles,newIndex,oldParticles,oldIndex)
! Copies one particle from one data structure to another.
	type(Particles), intent(inout) :: newParticles
	type(Particles), intent(in) :: oldParticles
	integer(kint), intent(in) :: newIndex, &
								 oldIndex
	! Error checking
	if ( newIndex > newParticles%N_Max .OR. oldIndex > oldParticles%N ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyParticleByIndex ERROR: invalid index')
		return
	endif

	newParticles%x(:,newIndex) = oldParticles%x(:,oldIndex)
	newParticles%x0(:,newIndex) = oldParticles%x0(:,oldIndex)
	newParticles%u(:,newIndex) = oldParticles%u(:,oldIndex)
	if ( associated(oldParticles%tracer) ) then
		if ( associated(newParticles%tracer) ) then
			newParticles%tracer(newIndex,:) = oldParticles%tracer(oldIndex,:)
		else
			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyParticleByIndex WARNING: Cannot assign tracer.')
		endif
	endif
	if ( associated(oldParticles%relVort)) then
		if ( associated(newParticles%relVort)) then
			newParticles%relVort(newIndex) = oldParticles%relVort(oldIndex)
		else
			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyParticleByIndex WARNING: Cannot assign relVort.')
		endif
	endif
	if ( associated(oldParticles%absVort)) then
		if ( associated(newParticles%absVort)) then
			newParticles%absVort(newIndex) = oldParticles%absVort(oldIndex)
		else
			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyParticleByIndex WARNING: Cannot assign absVort.')
		endif
	endif
	if ( associated(oldParticles%potVort)) then
		if ( associated(newParticles%potVort)) then
			newParticles%potVort(newIndex) = oldParticles%potVort(oldIndex)
		else
			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyParticleByIndex WARNING: Cannot assign potVort.')
		endif
	endif
	if ( associated(oldParticles%h)) then
		if ( associated(newParticles%h)) then
			newParticles%h(newIndex) = oldParticles%h(oldIndex)
		else
			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyParticleByIndex WARNING: Cannot assign h.')
		endif
	endif
	if ( associated(oldParticles%div)) then
		if ( associated(newParticles%div)) then
			newParticles%div(newIndex) = oldParticles%div(oldIndex)
		else
			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyParticleByIndex WARNING: Cannot assign div.')
		endif
	endif
	if ( associated(oldParticles%ke)) then
		if ( associated(newParticles%ke)) then
			newParticles%ke(newIndex) = oldParticles%ke(oldIndex)
		else
			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyParticleByIndex WARNING: Cannot assign ke.')
		endif
	endif
!	if ( associated(oldParticles%pe)) then
!		if ( associated(newParticles%pe)) then
!			newParticles%pe(newIndex) = oldParticles%pe(oldIndex)
!		else
!			call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'CopyParticleByIndex WARNING: Cannot assign pe.')
!		endif
!	endif
end subroutine


function GetNTracerParticles(self)
	type(Particles), intent(in) :: self
	integer(kint) :: GetNTracerParticles
	if ( associated(self%tracer)) then
		GetNTracerParticles = size(self%tracer,2)
	else
		GetNTracerParticles = 0
	endif
end function


function ParticleMax(panelKind,maxNest)
! Returns the maximum number of particles needed for memory allocation.
	integer(kint) :: particleMax
	integer(kint), intent(in) :: panelKind, maxNest
	if ( panelKind == TRI_PANEL ) then
		particleMax = 2 + 10*4**maxNest
	elseif (panelKind == QUAD_PANEL ) then
		particleMax = 2 + 6*4**maxNest
	endif
end function

!
!----------------
! Module methods : type-specific functions
!----------------
!
subroutine LogStatsParticles(self,aLog,message)
! send data about particles object to log (usually console) for output.
	! Calling parameters
	type(Particles), intent(in) :: self
	type(Logger), intent(inout) :: aLog
	character(len=*), intent(in), optional :: message
	! local variables
	character(len=24) :: key
	integer(kint) :: k, nTracer, j
	real(kreal) :: maxU, minU, uj

	if (present(message)) then
		call StartSection(aLog,'Particles Stats : ',message)
	else
		call StartSection(aLog,'Particles Stats : ')
	endif

	key = 'N = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,self%N)
	key = 'N_Max = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,self%N_Max)
	key = 'max. velocity = '
	maxU = 0.0_kreal
	minU = 0.0_kreal
	do j=1,self%N
		 uj = sqrt(sum( self%u(:,j)*self%u(:,j)))
		 if ( uj > maxU ) maxU = uj
		 if ( uj < minU ) minU = uj
	enddo
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,maxU)
	key = 'min. velocity = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,minU)
	if ( associated(self%absVort) ) then
		key = 'Max absVort = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,maxVal(self%absVort(1:self%N)))
		key = 'Min absVort = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,minval(self%absVort(1:self%N)))
	endif
	if ( associated(self%relVort) ) then
		key = 'Max relVort = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,maxVal(self%relVort(1:self%N)))
		key = 'Min relVort = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,minval(self%relVort(1:self%N)))
	endif
	if ( associated(self%potVort) ) then
		key = 'Max potVort = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,maxVal(self%potVort(1:self%N)))
		key = 'Min potVort = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,minval(self%potVort(1:self%N)))
	endif
	if ( associated(self%h) ) then
		key = 'Max h = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,maxVal(self%h(1:self%N)))
		key = 'Min h = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,minval(self%h(1:self%N)))
	endif
	if ( associated(self%div) ) then
		key = 'Max div = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,maxVal(self%div(1:self%N)))
		key = 'Min div = '
		call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,minval(self%div(1:self%N)))
	endif
	if ( associated(self%tracer)) then
		nTracer = GetNTracer(self)
		do k=1,nTracer
			write(key,'(A,I2,A)') 'Max tracer',k,' = '
			call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,maxVal(self%tracer(1:self%N,k)))
			write(key,'(A,I2,A)') 'Min tracer',k,' = '
			call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,minval(self%tracer(1:self%N,k)))
		enddo
	endif
	call EndSection(aLog)
end subroutine

subroutine InitLogger(aLog,rank)
! Initialize a logger for this module and processor
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I3,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

end module
