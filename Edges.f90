module EdgesModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	Defines the edges data structure used by icosahedral triangle and cubed
!	sphere Lagrangian meshes of the sphere.
!
!	no wings : simple edge data structure (as oppose to winged edges used by Voronoi meshes) relies on
!		constant polygonal structure (triangles or quadrilaterals) across the whole mesh.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
use NumberKindsModule
use LoggerModule


implicit none
private
public Edges
public New, Delete, Copy
public LogStats
public EdgeMax
!
!----------------
! Types and module constants
!----------------
!
type Edges
	! Grid variables
	integer(kint), pointer :: verts(:,:)	! vert(1,j) and vert(2,j) are the start and end vertices of edge j
	integer(kint), pointer :: leftPanel(:)	! leftPanel(j) is index of the panel to the left of edge j
	integer(kint), pointer :: rightPanel(:) ! rightPanel(j) is index of panel to the right of edge j
	integer(kint), pointer :: children(:,:)	! if edge is divided, children(1,j) and children(2,j) are indices of edge children
	logical(klog), pointer :: hasChildren(:)! hasChildren(j) is .TRUE. if edge(j) has been divided
	integer(kint) :: N						! current number of edges in memory
	integer(kint) :: N_Max					! maximum number of edges allowed in memory
	! Physical variables
	real(kreal), pointer :: length(:)		! length(j) is the length of edge j in physical space
	real(kreal), pointer :: length0(:)		! length0(j) is the length of edge j in Lagrangian parameter space
end type

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'Edges'
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
	module procedure CopyEdgeByIndex
	module procedure CopyEdges
end interface

interface LogStats
	module procedure LogEdgesStats
end interface

contains
!
!----------------
! Standard methods : Constructor / Destructor, Copy
!----------------
!
subroutine NewPrivate(self,nMax)
! Allocates memory for an instance of the edges data structure. Sets initial state to zero/null.
	type(Edges), intent(out) :: self
	integer(kint), intent(in) :: nMax

	if (.NOT. logInit) then
		call InitLogger(log,procRank)
	endif

	! Error checking
	if ( nMax <= 0 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,' invalid nMax.')
		return
	endif

	allocate(self%verts(2,nMax))
	self%verts = 0
	allocate(self%leftPanel(nMax))
	self%leftPanel = 0
	allocate(self%rightPanel(nMax))
	self%rightPanel = 0
	allocate(self%children(2,nMax))
	self%children = 0
	allocate(self%hasChildren(nMax))
	self%hasChildren = .False.
	allocate(self%length(nMax))
	self%length = 0.0_kreal
	allocate(self%length0(nMax))
	self%length0 = 0.0_kreal
	self%N = 0
	self%N_Max = nMax
end subroutine


subroutine DeletePrivate(self)
! Frees memory associated with an instance of edges data structure.
	type(Edges), intent(inout) :: self
	deallocate(self%verts)
	deallocate(self%leftPanel)
	deallocate(self%rightPanel)
	deallocate(self%children)
	deallocate(self%hasChildren)
	deallocate(self%length)
	deallocate(self%length0)
	self%N = 0
	self%N_Max = 0
end subroutine

subroutine CopyEdges(newEdges,oldEdges)
! copies an entire edges data structure
	type(Edges), intent(inout) :: newEdges
	type(Edges), intent(in) :: oldEdges
	integer(kint) :: j
	! Error checking
	if ( newEdges%N_Max < oldEdges%N ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyEdges ERROR : not enough memory')
		return
	endif

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering CopyParticles.')

	do j=1,oldEdges%N
		newEdges%verts(:,j) = oldEdges%verts(:,j)
		newEdges%leftPanel(j) = oldEdges%leftPanel(j)
		newEdges%rightPanel(j) = oldEdges%rightPanel(j)
		newEdges%children(:,j) = oldEdges%children(:,j)
		newEdges%hasChildren(j) = oldEdges%hasChildren(j)
		newEdges%length(j) = oldEdges%length(j)
		newEdges%length0(j) = oldEdges%length0(j)
	enddo
	newEdges%N = oldEdges%N
end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine CopyEdgeByIndex(newEdges,newIndex,oldEdges,oldIndex)
! copies one edge's information from one edges data structure to another.
	type(Edges), intent(inout) :: newEdges
	type(Edges), intent(in) :: oldEdges
	integer(kint), intent(in) :: newIndex, &
							 	 oldIndex
	! Error checking
	if ( oldIndex > oldEdges%N .OR. newIndex > newEdges%N_Max) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'CopyEdgeByIndex ERROR : out of bound')
		return
	endif


	newEdges%verts(:,newIndex) = oldEdges%verts(:,oldIndex)
	newEdges%leftPanel(newIndex) = oldEdges%leftPanel(oldIndex)
	newEdges%rightPanel(newIndex) = oldEdges%rightPanel(oldIndex)
	newEdges%children(:,newIndex) = oldEdges%children(:,oldIndex)
	newEdges%hasChildren(newIndex) = oldEdges%hasChildren(oldIndex)
	newEdges%length(newIndex) = oldEdges%length(oldIndex)
	newEdges%length0(newIndex) = oldEdges%length0(oldIndex)
end subroutine

function EdgeMax(panelKind,maxNest)
! Returns the maximum number of edges needed for memory allocation.
	integer(kint) :: edgeMax
	integer(kint), intent(in) :: panelKind, maxNest
	integer(kint) :: j
	edgeMax = 0
	if ( panelKind == TRI_PANEL ) then
		do j=0,maxNest
			edgeMax = edgeMax + 30*4**j
		enddo
	elseif (panelKind == QUAD_PANEL ) then
		do j=0,maxNest
			edgeMax = edgeMax + 12*4**j
		enddo
	endif
end function

!
!----------------
! Module methods : type-specific functions
!----------------
!
subroutine LogEdgesStats(self,aLog,message)
! outputs data associated with an Edges object to a Logger
	type(Edges), intent(in) :: self
	type(Logger), intent(inout) :: aLog
	character(len=*), intent(in), optional :: message
	character(len=24) :: key
	real(kreal) :: lmax, lmin
	integer(kint) :: j

	if ( present(message)) then
		call StartSection(aLog,'Edges Stats : ',message)
	else
		call StartSection(aLog,'Edges Stats : ')
	endif



	key = 'N = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,self%N)
	key = 'N_Max = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,self%N_Max)

	lmax = maxval(self%length(1:self%N))
	lmin = lmax
	do j=1,self%N
		if ( .NOT. self%hasChildren(j) .AND. self%length(j) < lmin ) lmin = self%length(j)
	enddo
	key = 'Max Edge length = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,lmax)
	key = 'Min Edge Length = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,lmin)

	lmax = maxval(self%length0(1:self%N))
	lmin = lmax
	do j=1,self%N
		if ( .NOT. self%hasChildren(j) .AND. self%length(j) < lmin) lmin = self%length0(j)
	enddo

	key = 'Max init length = '
	call LogMessage(alog,TRACE_LOGGING_LEVEL,key, lmax)
	key = 'Min init length = '
	call LogMessage(aLog,TRACE_LOGGING_LEVEL,key,lmin)
	call EndSection(aLog)
end subroutine

subroutine InitLogger(aLog,rank)
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I3,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

end module
