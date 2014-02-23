module STRIPACKInterfaceModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	Provides an interface for SphereMesh objects and SphereMeshOutput objects to
!	use the cubic Hermite interpolation procedures provided by ssrfpack.f.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
!	Renka, R. "Algorithm 772 : STRIPACK : Delaunay triangulation and Voronoi diagram on the surface of a sphere."
!		ACM Transactions on Mathematical Software, 23 : 416-434. 1997.
!
use NumberKindsModule
use LoggerModule
use ParticlesModule
use EdgesModule
use PanelsModule
use SphereMeshModule

implicit none
private
public STRIPACKData
public New, Delete
public DelaunayTriangulation

!
!----------------
! Types and module constants
!----------------
!

type STRIPACKData
	real(kreal), pointer :: x(:) => null(), &	! x coordinates of Delaunay triangulation nodes
						    y(:) => null(), &	! y coordinates of Delaunay triangulation nodes
						    z(:) => null()		! z coordinates of Delaunay triangulation nodes

	integer(kint), pointer :: list(:) => null(), &	! pointers to Delaunay nodes according to STRIPACK data structures
							  lptr(:) => null(), &	! pointers to Delaunay nodes according to STRIPACK data structures
							  lend(:) => null(), &	! pointers to Delaunay nodes according to STRIPACK data structures
							  activeMap(:) => null() ! indices of active panels in sphereMesh, set by GatherPanels

	type(Particles), pointer :: particles => null()	! pointers to SphereMesh data, input to STRIPACK
	type(Panels), pointer :: activePanels => null()	! pointers to SphereMesh data, input to STRIPACK
	integer(kint) :: n								! number of nodes from SphereMesh data to use with STRIPACK
	logical(klog) :: memoryReady = .FALSE.
end type


!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28) :: logKey='STRIPACK'
integer(kint), parameter :: logLevel = TRACE_LOGGING_LEVEL
character(len=28) :: formatString
character(len=128) :: logString
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
subroutine NewPrivate(self,aMesh)
!	Allocates memory for a STRIPACK Delaunay triangulation of a SphereMesh object.
	type(STRIPACKData), intent(out) :: self
	type(SphereMesh), intent(in) :: aMesh
	type(Panels), pointer :: aPanels
	type(Panels), pointer :: passivePanels
	integer(kint), allocatable :: passiveMap(:)
	integer(kint) :: nActive, nPassive, panelKind, nTracer, problemKind, j, n
	real(kreal) :: norm

	if ( .NOT. logInit) call InitLogger(log,procRank)

	self%particles => aMesh%particles
	aPanels => aMesh%panels

	n = aPanels%N_Active + self%particles%N
	self%n = n

	nActive = aPanels%N_Active
	nPassive = aPanels%N - aPanels%N_active
	nTracer = GetNTracer(aPanels)
	panelKind = GetPanelKind(aPanels)
	problemKind = GetProblemKind(aPanels)

	! Allocate Delaunay triangulation data arrays
	allocate(self%x(n))
	self%x = 0.0_kreal
	allocate(self%y(n))
	self%y = 0.0_kreal
	allocate(self%z(n))
	self%z = 0.0_kreal
	allocate(self%lend(n))
	self%lend = 0
	allocate(self%list(6*n-12))
	self%list = 0
	allocate(self%lptr(6*n-12))
	self%lptr = 0

	! Associate pointers with source data
	allocate(self%activePanels)
	call New(self%activePanels,nActive,panelKind,nTracer,problemKind)
	self%activePanels%N = nActive
	self%activePanels%N_Active = nActive
	allocate(self%activeMap(aPanels%N_Active))
	self%activeMap = 0

	allocate(passivePanels)
	call New(passivePanels,nPassive,panelKind,nTracer,problemKind)
	passivePanels%N = nPassive
	allocate(passiveMap(aPanels%N-aPanels%N_Active))
	passiveMap = 0

	call GatherPanels(aPanels,self%activePanels,self%activeMap,passivePanels,passiveMap)

	!  Get STRIPACK Delauanay source data from SphereMesh
	!	Renormalize input to unit sphere for STRIPACK
	do j=1,nActive
		norm = sqrt(sum(self%activePanels%x(:,j)*self%activePanels%x(:,j)))
		self%x(j) = self%activePanels%x(1,j)/norm
		self%y(j) = self%activePanels%x(2,j)/norm
		self%z(j) = self%activePanels%x(3,j)/norm
	enddo
	do j=1,self%particles%N
		norm = sqrt(sum(self%particles%x(:,j)*self%particles%x(:,j)))
		self%x(nActive+j) = self%particles%x(1,j)/norm
		self%y(nActive+j) = self%particles%x(2,j)/norm
		self%z(nActive+j) = self%particles%x(3,j)/norm
	enddo

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,' STRIPACK memory allocated.')

	self%memoryReady = .TRUE.

	call DelaunayTriangulation(self)

	call Delete(passivePanels)
	deallocate(passiveMap)
	nullify(aPanels)
end subroutine

subroutine DeletePrivate(self)
!	Free memory associated with an instance of STRIPACKData
	type(STRIPACKData), intent(inout) :: self
	deallocate(self%x)
	deallocate(self%y)
	deallocate(self%z)
	deallocate(self%list)
	deallocate(self%lptr)
	deallocate(self%lend)
	deallocate(self%activeMap)
	call Delete(self%activePanels)
	nullify(self%activePanels)
	nullify(self%particles)
	self%memoryReady = .FALSE.
end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine DelaunayTriangulation(self)
!	Build the Delaunay triangulation of the SphereMesh data associated with self.
!	NOTE : STRIPACKData object must be preallocated and assigned to a SphereMesh object (using "New")
!		   prior to calling this subroutine.
!
	type(STRIPACKData), intent(inout) :: self
	real(kreal), allocatable :: dist(:)
	integer(kint), allocatable :: near(:), next(:)
	integer(kint) :: lnew, errCode

	if ( .NOT. self%memoryReady ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,' memory not allocated.')
		return
	endif

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,' Entering DelaunayTriangulation.')

	! Allocate STRIPACK TRMESH working arrays
	allocate(near(self%n))
	near = 0
	allocate(next(self%n))
	next = 0
	allocate(dist(self%n))
	dist = 0.0_kreal

	! Build the Delaunay triangulation
	call TRMESH(self%N,&
				self%x,self%y,self%z,&
				self%list,self%lptr,self%lend,&
				lnew,near,next,dist,errCode)

	! Check for errors
	if ( errCode == -1 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'TRMESH ERROR :',' found n < 3 points.')
	elseif (errCode == -2 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'TRMESH ERROR :',' found first three nodes to be colinear.')
	elseif ( errCode > 0 ) then
		write(logString,'(A,I8,A)') ' node ',errCode,' is a duplicate.'
		call LogMessage(log,ERROR_LOGGING_LEVEL,'TRMESH ERROR :',logString)
	endif

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,' TRMESH complete.')

	! delete working arrays
	deallocate(near)
	deallocate(next)
	deallocate(dist)
end subroutine

!
!----------------
! Module methods : type-specific functions
!----------------
!
subroutine InitLogger(aLog,rank)
	type(Logger), intent(inout) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine


end module
