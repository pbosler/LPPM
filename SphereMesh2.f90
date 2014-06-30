module SphereMeshModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	Defines the mesh data structure used by icosahedral triangle and cubed
!	sphere Lagrangian meshes of the sphere.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
!----------------
! USAGE :  This module provides the basic data structures for particle-panel methods on the
! unit sphere.  Memory is preallocated statically for speed, assuming the column-priority indexing
! typical to Fortran.  This means that information associated with panel j will be contained in
! column j of the panels' structure's various arrays.
!		For example : the physical position of panel j is stored at panels%x(:,j), and the area
!				associated with panel j is stored at panels%area(j).  Vertices that make up the corners
!				of panel j are stored at panels%vertices(:,j), etc...  The same is true for particles and
!				edges, so that the indices to the particles that make up edge k are stored at
!				edges%verts(:,k).
!		Memory is handled with the "New" and "Delete" interfaces and procedures; grids are initialized
! 				with subroutines labeled "Init*".
!----------------
use NumberKindsModule
use LoggerModule
use IntegerListModule
use SphereGeomModule
!use lapack95
use ParticlesModule
use EdgesModule
use PanelsModule


implicit none

private
public SphereMesh
public New, Delete, Copy
public LogStats
public CCWEdgesAndParticlesAroundPanel, FindAdjacentPanels
public Renormalize, ResetSphereArea
public TotalMass, TotalAbsVort, TotalRelVort, TotalEnstrophy
public CountSubTriangles
public LocatePoint
public DividePanel
public GetRossbyNumber
public ResetLagrangianParameter
public MaximumCirculation, MaximumVorticityVariation
public MaximumLagrangianParameterVariation, AverageLagrangianParameterVariation
public MaximumTracerMass, MaximumTracerVariation
!
!----------------
! Types and module constants
!----------------
!
type SphereMesh
	type(Particles), pointer :: particles => null()
	type(Edges), pointer :: edges => null()
	type(Panels), pointer :: panels => null()
	integer(kint) :: panelKind
	integer(kint) :: initNest
	integer(kint) :: AMR
	integer(kint) :: nTracer
	integer(kint) :: problemKind
	real(kreal) :: totalE
	real(kreal) :: totalEnstrophy
	real(kreal), pointer :: totalMass(:)
end type
!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'SphereMesh'
integer(kint), parameter :: logLevel = TRACE_LOGGING_LEVEL
character(len=128) :: logString
character(len=24) :: formatString
!
!----------------
! Interfaces
!----------------
!
interface New
	module procedure NewMesh
end interface

interface Delete
	module procedure DeletePrivate
end interface

interface Copy
	module procedure CopyPrivate
end interface

interface LogStats
	module procedure LogStatsSphereMesh
end interface

interface GetNTracer
	module procedure GetNTracerMesh
end interface

interface ResetLagrangianParameter
	module procedure ResetLagrangianParameterPrivate
end interface

contains
!
!----------------
! Standard methods : Constructor / Destructor
!----------------
!


subroutine NewMesh(self,panelKind,initNest,AMR,nTracer,problemKind)
! Allocates memory for a SphereMesh object and initializes a discretization of the sphere with
! the desired panelKind.
!	Tracers and Vorticity data must be set separately.
!	Returns a uniform mesh only; AMR routines must be called separately using RefineRemeshModule.
	type(SphereMesh), intent(out) :: self
	integer(kint), intent(in) :: panelKind, &
								 initNest, &
								 AMR, &
								 nTracer, &
								 problemKind
	integer(kint) :: nPanels, nParticles, nEdges

	if ( .NOT. logInit ) then
		call InitLogger(log,procRank)
	endif
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering NewMesh...')

	! Allocate the pointers
	allocate(self%particles)
	allocate(self%edges)
	allocate(self%panels)

	self%panelKind = panelKind
	self%initNest = initNest
	self%AMR = AMR
	self%nTracer = nTracer
	self%problemKind = problemKind
	allocate(self%totalMass(nTracer))
	self%totalMass = 0.0_kreal

	! compute the array sizes
	nPanels = PanelMax(panelKind,initNest + AMR)
	nParticles = ParticleMax(panelKind,initNest + AMR)
	nEdges = EdgeMax(panelKind,initNest+AMR)

	call New(self%particles,nParticles,panelKind,nTracer,problemKind)
	call New(self%edges,nEdges)
	call New(self%panels,nPanels,panelKind,nTracer,problemKind)

	if ( panelKind == TRI_PANEL) then
		call InitIcosTriMesh(self,initNest)
	elseif ( panelKind == QUAD_PANEL) then
		call InitCubedSphere(self,initNest)
	endif

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'NewMesh complete.')
end subroutine

subroutine DeletePrivate(self)
! 	Deletes and frees memory associated with a SphereMesh object.
	type(SphereMesh), intent(inout) :: self
	self%panelKind = 0
	self%initNest = -1
	self%AMR = 0
	self%nTracer = -1
	self%problemKind = 0
	call Delete(self%particles)
	call Delete(self%edges)
	call Delete(self%panels)
	deallocate(self%particles)
	deallocate(self%edges)
	deallocate(self%panels)
end subroutine


subroutine CopyPrivate(newMesh,oldMesh)
! Copies by replicating memory from one SphereMesh object into another.
! Overwrites data currently in destination object.
	type(SphereMesh), intent(inout) :: newMesh
	type(SphereMesh), intent(in) :: oldMesh

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering CopyMesh.')

	call Copy(newMesh%particles,oldMesh%particles)
	call Copy(newMesh%edges,oldMesh%edges)
	call Copy(newMesh%panels,oldMesh%panels)
	newMesh%panelKind = oldMesh%panelKind
	newMesh%initNest = oldMesh%initNest
	newMesh%AMR = oldMesh%AMR
	newMesh%nTracer = oldMesh%nTracer
	newMesh%problemKind = oldMesh%problemKind
end subroutine

!
!----------------
! Public functions
!----------------
!

subroutine LogStatsSphereMesh(self,alog,message)
! Prints SphereMesh statistical data to log.
	type(SphereMesh), intent(in) :: self
	type(Logger), intent(inout) :: aLog
	character(len=*), intent(in), optional :: message

	if ( present(message) ) then
		call LogStats(self%particles,aLog,message)
		call LogStats(self%edges,aLog,message)
		call LogStats(self%panels,aLog,message)
	else
		call LogStats(self%particles,alog)
		call LogStats(self%edges,alog)
		call LogStats(self%panels,alog)
	endif
end subroutine


function GetNTracerMesh(self)
	type(SphereMesh), intent(in) :: self
	integer(kint) :: GetNTracerMesh
	GetNTracerMesh = self%nTracer
end function


subroutine Renormalize(self)
! 	Reprojects SphereMesh physical coordinates back to spherical surface.
	type(SphereMesh), intent(inout) :: self
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j
	aParticles => self%particles
	aPanels => self%panels
	do j=1,aParticles%N
		aParticles%x(:,j) = aParticles%x(:,j)/sqrt(sum(aParticles%x(:,j)*aParticles%x(:,j)))*EARTH_RADIUS
	enddo
	do j=1,aPanels%N
		aPanels%x(:,j) = aPanels%x(:,j)/sqrt(sum(aPanels%x(:,j)*aPanels%x(:,j)))*EARTH_RADIUS
	enddo
end subroutine


function PanelArea(self,panelIndex)
! 	Interface for calling the appropriate private panel area function.
	type(SphereMesh), intent(inout) :: self
	integer(kint), intent(in) :: panelIndex
	real(kreal) :: PanelArea
	integer(kint) :: panelKind
	panelKind = GetPanelKind(self%panels)
	if ( panelKind == QUAD_PANEL ) then
		PanelArea =  QuadPanelArea(self,panelIndex)
	endif
end function


subroutine ResetSphereArea(self)
! Resets SphereMesh panel areas for use in flows with nonzero divergence.
	type(SphereMesh), intent(inout) :: self
	integer(kint) :: j
	do j=1,self%panels%N
		if ( .NOT. self%panels%hasChildren(j) ) then
			self%panels%area(j) = PanelArea(self,j)
		else
			self%panels%area(j) = 0.0_kreal
		endif
	enddo
end subroutine

subroutine CCWEdgesAndParticlesAroundPanel(edgeList,vertList,nVerts,self,panelIndex)
! 	Returns a counter-clockwise ordered list of edges and a counter-clockwise ordered list
! 	of vertices around a panel.  For use with AMR, this subroutine does not assume a panel
!	has only 4 vertices-- panels may have as many as 8.
	integer(kint), intent(inout) :: edgeList(:), vertList(:), nVerts
	type(SphereMesh), intent(in) :: self
	integer(kint), intent(in) :: panelIndex
	type(Particles), pointer :: aParticles
	type(Edges), pointer :: anEdges
	type(Panels), pointer :: aPanels
	integer(kint) :: panelKind, k, edgeK, vertK

	aParticles=>self%particles
	anEdges=>self%edges
	aPanels=>self%panels

	panelKind = GetPanelKind(aPanels)
	if ( panelKind == QUAD_PANEL ) then
		if ( size(edgeList) < 8 .OR. size(vertList) < 8 ) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,'CCWEdgesAndParticlesAroundPanel ERROR : ',' output array too small.')
			return
		endif
	endif

	edgelist = 0
	vertlist = 0
	edgeK = 0
	vertK = 0

	do k=1,panelKind
		if ( anEdges%hasChildren(aPanels%edges(k,panelIndex)) ) then
			if ( anEdges%leftPanel(aPanels%edges(k,panelIndex)) == panelIndex) then
				edgeList(edgeK+1) = anEdges%children(1,aPanels%edges(k,panelIndex))
				vertList(vertK+1) = anEdges%verts(1,anEdges%children(1,aPanels%edges(k,panelIndex)))
				edgeList(edgeK+2) = anEdges%children(2,aPanels%edges(k,panelIndex))
				vertList(vertK+2) = anEdges%verts(2,anEdges%children(1,aPanels%edges(k,panelIndex)))
			else
				edgeList(edgeK+1) = anEdges%children(2,aPanels%edges(k,panelIndex))
				vertList(vertK+1) = anEdges%verts(2,anEdges%children(2,aPanels%edges(k,panelIndex)))
				edgeList(edgeK+2) = anEdges%children(1,aPanels%edges(k,panelIndex))
				vertList(vertK+2) = anEdges%verts(1,anEdges%children(2,aPanels%edges(k,panelIndex)))
			endif
			edgeK = edgeK + 2
			vertK = vertK+2
		else
			edgeList(edgeK+1) = aPanels%edges(k,panelIndex)
			if ( anEdges%leftPanel(aPanels%edges(k,panelIndex)) == panelIndex) then
				vertList(vertK+1) = anEdges%verts(1,aPanels%edges(k,panelIndex))
			else
				vertList(vertK+1) = anEdges%verts(2,aPanels%edges(k,panelIndex))
			endif
			edgeK = edgeK+1
			vertK = vertK+1
		endif
	enddo
	nVerts = vertK
end subroutine

function TotalMass(self,tracerNumber)
!	Computes the mass integral of passive tracers over the sphere.
	type(SphereMesh), intent(in) :: self
	integer(kint), intent(in) :: tracerNumber
	real(kreal) :: TotalMass
	type(Panels), pointer :: aPanels
	aPanels=>self%panels

	! Error checking
	if ( tracerNumber > GetNTracer(aPanels) ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'TotalMass ERROR : ','invalid tracerNumber.')
		return
	endif
	TotalMass = sum(aPanels%tracer(:,tracerNumber)*aPanels%area)
end function


function TotalAbsVort(self)
!	Calculates the integral of absolute vorticity over the sphere.
!	For diagnostic purposes (this quantity should be zero).
	type(SphereMesh), intent(in) :: self
	real(kreal) :: TotalAbsVort
	type(Panels), pointer :: aPanels

	aPanels=>self%panels

	! Error checking
	if ( .NOT. associated(aPanels%absVort) ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'TotalAbsVort ERROR : ','absVort not found.')
		return
	endif
	TotalAbsVort = sum(aPanels%absVort*aPanels%area)
end function

function TotalEnergy(self)
!	Calculates the integral of absolute vorticity over the sphere.
!	For diagnostic purposes (this quantity should be zero).
	type(SphereMesh), intent(in) :: self
	real(kreal) :: TotalEnergy
	type(Panels), pointer :: aPanels

	aPanels=>self%panels

	! Error checking
	if ( .NOT. associated(aPanels%ke) ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'TotalAbsVort ERROR : ','absVort not found.')
		return
	endif
	TotalEnergy = 0.5_kreal*sum(aPanels%ke*aPanels%area)
	if ( associated(aPanels%h) ) then
		TotalEnergy = TotalEnergy + GRAV*sum(aPanels%h*aPanels%area)
	endif
end function

function GetRossbyNumber(self)
	real(kreal) :: GetRossbyNumber
	type(SphereMesh), intent(in) :: self
	GetRossbyNumber = maxval(abs(self%panels%relVort))/(2.0_kreal*OMEGA)
end function

function TotalRelVort(self)
!	Calculates the integral of relative vorticity over the sphere.
!	For diagnostic purposes (this quantity should be zero).
	type(SphereMesh), intent(in) :: self
	real(kreal) :: TotalRelVort
	type(Panels), pointer :: aPanels
	aPanels=>self%panels

	! Error checking
	if ( .NOT. associated(aPanels%relVort) ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'TotalRelVort ERROR : ','relVort not found.')
		return
	endif
	TotalRelVort = sum(aPanels%relVort*aPanels%area)
end function


function TotalEnstrophy(self)
!	Computes the total enstrophy as a surface integral over the sphere.
	type(SphereMesh), intent(in) :: self
	real(kreal) :: TotalEnstrophy
	type(Panels), pointer :: aPanels

	aPanels=>self%panels

	! Error checking
	if ( .NOT. associated(aPanels%relVort) ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'TotalEnstrophy ERROR : ',' relVort not found.')
		return
	endif

	totalEnstrophy = 0.5*sum(aPanels%relVort*aPanels%relVort*aPanels%area)
end function


function LocatePoint(self,xyz)
!  	Primary function to locate what panel of a SphereMesh object contains the point xyz.
! 	Calls private subroutines TreeSearch and WalkSearch.
	integer(kint) :: LocatePoint
	type(SphereMesh), intent(in) :: self
	real(kreal), intent(in) :: xyz(3)
	integer(kint) :: rootPanel, startWalk
	! Part 1 : tree search
	rootPanel = FindClosestRootPanel(self,xyz)
	call LocatePointTreeSearch(startWalk,self,xyz,rootPanel)
	! Part 2 : walk search
	call LocatePointWalkSearch(LocatePoint,self,xyz,startWalk)
end function


function CountSubtriangles(self)
!	Returns the total number of panel subtriangles in a SphereMesh object.
	integer(kint) :: CountSubTriangles
	type(SphereMesh), intent(in) :: self
	integer(kint) :: edgeList(8), vertList(8), nVerts, j
	CountSubTriangles = 0
	do j=1,self%panels%N
		if ( .NOT. self%panels%hasChildren(j) ) then
			call CCWEdgesAndParticlesAroundPanel(edgeList,vertList,nVerts,self,j)
			CountSubTriangles = CountSubTriangles + nVerts
		endif
	enddo
end function

subroutine ResetLagrangianParameterPrivate(self)
  type(SphereMesh), intent(in) :: self
  integer(kint) :: j
  type(Particles), pointer :: aParticles
  type(Panels), pointer :: aPanels

  aParticles => self%particles
  aPanels => self%panels

  do j=1,aParticles%N
     aParticles%x0(:,j) = aParticles%x(:,j)
  enddo
  do j=1,aPanels%N
     aPanels%x0(:,j) = aPanels%x(:,j)
  enddo

end subroutine

function MaximumCirculation(self)
	real(kreal) :: MaximumCirculation
	type(SphereMesh), intent(in) :: self
	!
	type(Panels), pointer :: aPanels
	aPanels => self%panels
	MaximumCirculation = maxval( abs(aPanels%relvort(1:aPanels%N)) * aPanels%area(1:aPanels%N) )
end function

function MaximumTracerMass(self, tracerID)
	real(kreal) :: MaximumTracerMass
	type(SphereMesh), intent(in) :: self
	integer(kint), intent(in) :: tracerID
	!
	MaximumTracerMass = maxval( self%panels%tracer(1:self%panels%N, tracerID) * self%panels%area(1:self%panels%N) )
end function

function MaximumTracerVariation(self, tracerID)
	real(kreal) :: MaximumTracerVariation
	type(SphereMesh), intent(in) :: self
	integer(kint), intent(in) :: tracerID
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	real(kreal) :: tmax, tmin

	aParticles => self%particles
	aPanels => self%panels

	tmax = max( maxval(aParticles%tracer(1:aParticles%N, tracerID)), maxval( aPanels%tracer(1:aPanels%N, tracerID)) )
	tmin = min( minval(aParticles%tracer(1:aParticles%N, tracerID)), minval( aPanels%tracer(1:aPanels%N, tracerID)) )

	MaximumTracerVariation = tmax - tmin
end function


function MaximumVorticityVariation(self)
	real(kreal) :: MaximumVorticityVariation
	type(SphereMesh), intent(in) :: self
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	real(kreal) :: maxVort, minVort

	aParticles => self%particles
	aPanels => self%panels

	maxVort = max( maxval(aParticles%relVort(1:aParticles%N)), maxval(aPanels%relVort(1:aPanels%N)) )
	minVort = min( minval(aParticles%relVort(1:aParticles%N)), minval(aPanels%relVort(1:aPanels%N)) )

	MaximumVorticityVariation = maxVort - minVort
end function

function AverageLagrangianParameterVariation(self)
	real(kreal) :: AverageLagrangianParameterVariation
	type(SphereMesh), intent(in) :: self
 	AverageLagrangianParameterVariation = sqrt( 4.0_kreal * PI * EARTH_RADIUS * EARTH_RADIUS / self%panels%N_Active )
end function

function MaximumLagrangianParameterVariation(self)
	real(kreal) :: MaximumLagrangianParameterVariation
	type(SphereMesh), intent(in) :: self
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j, k, edgeList(8), vertList(8), nVerts
	real(kreal) :: maxx0(3), minx0(3)

	MaximumLagrangianParameterVariation = 0.0_kreal

	aParticles => self%particles
	aPanels => self%panels

	do j = 1, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			maxx0 = aPanels%x0(:,j)
			minx0 = aPanels%x0(:,j)
			call CCWEdgesAndParticlesAroundPanel(edgeList, vertList, nVerts, self, j)
			do k = 1, nVerts
				if ( aParticles%x0(1,vertLIst(k)) > maxx0(1) ) maxx0(1) = aParticles%x0(1,vertlist(k))
				if ( aParticles%x0(1,vertList(k)) < minx0(1) ) minx0(1) = aParticles%x0(1,vertList(k))
				if ( aParticles%x0(2,vertLIst(k)) > maxx0(2) ) maxx0(2) = aParticles%x0(2,vertlist(k))
				if ( aParticles%x0(2,vertList(k)) < minx0(2) ) minx0(2) = aParticles%x0(2,vertList(k))
				if ( aParticles%x0(3,vertLIst(k)) > maxx0(3) ) maxx0(3) = aParticles%x0(3,vertlist(k))
				if ( aParticles%x0(3,vertList(k)) < minx0(3) ) minx0(3) = aParticles%x0(3,vertList(k))
			enddo
			if ( sum( maxx0 - minx0 ) > MaximumLagrangianParameterVariation ) &
				MaximumLagrangianParameterVariation = sum(maxx0 - minx0)
		endif
	enddo

end function

!
!----------------
! Module methods : type-specific functions
!----------------
!
subroutine InitCubedSphere(self,initNest)
!	Initializes a SphereMesh object with a cubed sphere discretization to a uniform resolution.
!	NOTE : memory must be preallocated before calling this subroutine.
	! Calling Parameters
	type(SphereMesh), intent(inout) :: self
	integer(kint), intent(in) :: initNest
	! Local variables
	real(kreal), parameter :: a = 1.0_kreal/dsqrt(3.0_kreal)
	integer(kint) :: j, k, startIndex, nOldPanels
	type(Particles), pointer :: aParticles
	type(Edges), pointer :: anEdges
	type(Panels), pointer :: aPanels

	! Error checking
	if ( initNest < 0 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'InitCubedSphere ERROR : invalid initNest.')
		return
	endif
	if ( (.NOT. associated(self%particles) .OR. (.NOT. associated(self%edges))) .OR. .NOT. associated(self%edges) ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'InitCubedSphere ERROR : memory not allcoated.')
		return
	endif

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'Initializing Cubed Sphere...')

	aParticles => self%particles
	anEdges => self%edges
	aPanels => self%panels

	! Setup the root cube
	!	Vertices = passive particles
	aParticles%N = 8
	!		face 1 = Atlantic, vertices have positve x coordinates
	aParticles%x(:,1) = [ a,-a, a]
	aParticles%x(:,2) = [ a,-a,-a]
	aParticles%x(:,3) = [ a, a,-a]
	aParticles%x(:,4) = [ a, a, a]
	!		face 2 = Asia, vertices have positive y coordinates
	aParticles%x(:,5) = [-a, a,-a]
	aParticles%x(:,6) = [-a, a, a]
	! 		face 3 = Pacific, vertices have negative x coordinates
	aParticles%x(:,7) = [-a,-a,-a]
	aParticles%x(:,8) = [-a,-a, a]
	!		face 4 = Americas, vertices have negative y coordinates
	!		face 5 = north pole, vertices have positive z coordinates
	!		face 6 = south pole, vertices have negative z coordinates
	aParticles%x(:,1:8) = aParticles%x(:,1:8)*EARTH_RADIUS
	!	Edges connect vertices
	anEdges%N = 12
	anEdges%verts(:,1) = [1,2]
	anEdges%verts(:,2) = [2,3]
	anEdges%verts(:,3) = [3,4]
	anEdges%verts(:,4) = [4,1]
	anEdges%verts(:,5) = [3,5]
	anEdges%verts(:,6) = [5,6]
	anEdges%verts(:,7) = [6,4]
	anEdges%verts(:,8) = [5,7]
	anEdges%verts(:,9) = [7,8]
	anEdges%verts(:,10)= [8,6]
	anEdges%verts(:,11)= [7,2]
	anEdges%verts(:,12)= [1,8]
	!	Edges separate faces
	anEdges%leftPanel(1) = 1
	anEdges%leftPanel(2) = 1
	anEdges%leftPanel(3) = 1
	anEdges%leftPanel(4) = 1
	anEdges%leftPanel(5) = 2
	anEdges%leftPanel(6) = 2
	anEdges%leftPanel(7) = 2
	anEdges%leftPanel(8) = 3
	anEdges%leftPanel(9) = 3
	anEdges%leftPanel(10) = 3
	anEdges%leftPanel(11) = 4
	anEdges%leftPanel(12) = 4
	anEdges%rightPanel(1) = 4
	anEdges%rightPanel(2) = 6
	anEdges%rightPanel(3) = 2
	anEdges%rightPanel(4) = 5
	anEdges%rightPanel(5) = 6
	anEdges%rightPanel(6) = 3
	anEdges%rightPanel(7) = 5
	anEdges%rightPanel(8) = 6
	anEdges%rightPanel(9) = 4
	anEdges%rightPanel(10) = 5
	anEdges%rightPanel(11) = 6
	anEdges%rightPanel(12) = 5
	! 	Edges have arc length
	!anEdges%length(1:12) = 1.230959417340776_kreal*EARTH_RADIUS

	! 	Faces of cube = Panels
	aPanels%N = 6
	aPanels%nest(1:6) = 0
	aPanels%N_Active = 6
	!	Panels connect edges
	aPanels%edges(:,1) = [1,2,3,4]
	aPanels%edges(:,2) = [3,5,6,7]
	aPanels%edges(:,3) = [6,8,9,10]
	aPanels%edges(:,4) = [9,11,1,12]
	aPanels%edges(:,5) = [12,4,7,10]
	aPanels%edges(:,6) = [11,8,5,2]
	! 	Panels have active particles at their centriods
	aPanels%x(:,1) = [1.0_kreal,0.0_kreal,0.0_kreal]
	aPanels%x(:,2) = [0.0_kreal,1.0_kreal,0.0_kreal]
	aPanels%x(:,3) = [-1.0_kreal,0.0_kreal,0.0_kreal]
	aPanels%x(:,4) = [0.0_kreal,-1.0_kreal,0.0_kreal]
	aPanels%x(:,5) = [0.0_Kreal,0.0_kreal,1.0_kreal]
	aPanels%x(:,6) = [0.0_kreal,0.0_kreal,-1.0_kreal]
	aPanels%x(:,1:6) = aPanels%x(:,1:6)*EARTH_RADIUS
	!	Panels connect to particles
	aPanels%vertices(:,1) = [1,2,3,4]
	aPanels%vertices(:,2) = [4,3,5,6]
	aPanels%vertices(:,3) = [6,5,7,8]
	aPanels%vertices(:,4) = [8,7,2,1]
	aPanels%vertices(:,5) = [8,1,4,6]
	aPanels%vertices(:,6) = [2,7,5,3]
	! 	Panels have area
	aPanels%area(1:6) = 4.0_kreal*PI*EARTH_RADIUS*EARTH_RADIUS/6.0_kreal

	! Root cube is set up
	!	Finish initialization by recording initial state
	aParticles%x0(:,1:8) = aParticles%x(:,1:8)
!	anEdges%length0(1:12) = anEdges%length(1:12)
	aPanels%x0(:,1:6) = aPanels%x(:,1:6)

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'... root cube ready; dividing panels.')

	! Divide root cube until desired initial nest level
	if ( initNest > 0 ) then
		startIndex = 1
		do k=1,initNest
			nOldPanels = aPanels%N
			do j=startIndex,nOldPanels
				if (.NOT. aPanels%hasChildren(j) ) then
					call DivideQuadPanel(self,j)
				endif
			enddo
			startIndex = nOldPanels
		enddo
	endif

!	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'... mesh ready; initializing dual mesh.')

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'... uniform quadrilateral mesh ready.')
end subroutine

subroutine InitIcosTriMesh(self, initNest)
	type(SphereMesh), intent(out) :: self
	integer(kint), intent(in) :: initNest
	!
	integer(kint) :: j, k
	type(Particles), pointer :: aparticles
	type(Edges), pointer :: anEdges
	type(Panels), pointer :: aPanels
	
	if ( initNest < 0 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, 'InitIcosTri ERROR : invalid initNest.')
		return
	endif
	
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey," initializing icosTri grid.")
	
	aParticles => self%particles
	anEdges => self%edges
	aPanels => self%panels
	
	!
	! set up root icosahedron
	!
	!	Particles
	aParticles%N = 12
	aParticles%x(:,1) = [0.0_kreal, 0.0_kreal, 1.0_kreal]
	aParticles%x(:,2) = [0.723606797749978969640917366873_kreal,&
						 0.525731112119133606025669084848_kreal,&
						 0.447213595499957939281834733746_kreal]
	aParticles%x(:,3) =[-0.276393202250021030359082633126_kreal,&
						 0.850650808352039932181540497063_kreal,&
						 0.447213595499957939281834733746_kreal]
	aParticles%x(:,4) =[-0.894427190999915878563669467492_kreal,&
						 0.0_kreal,&
						 0.447213595499957939281834733746_kreal]
	aParticles%x(:,5) =[-0.276393202250021030359082633127_kreal,&
						-0.850650808352039932181540497063_kreal,&
						 0.447213595499957939281834733746_kreal]
	aParticles%x(:,6) = [0.723606797749978969640917366873_kreal,&
						-0.525731112119133606025669084848_kreal,&
						 0.447213595499957939281834733746_kreal]
	aParticles%x(:,7) = [0.894427190999915878563669467492_kreal,&
						 0.0_kreal,&
						-0.447213595499957939281834733746_kreal]
	aParticles%x(:,8) = [0.276393202250021030359082633127_kreal,&
						 0.850650808352039932181540497063_kreal,&
						-0.447213595499957939281834733746_kreal]
	aParticles%x(:,9) =[-0.723606797749978969640917366873_kreal,&
						 0.525731112119133606025669084848_kreal,&
						-0.447213595499957939281834733746_kreal]
	aParticles%x(:,10)=[-0.723606797749978969640917366873_kreal,&
						-0.525731112119133606025669084848_kreal,&
						-0.447213595499957939281834733746_kreal]
	aParticles%x(:,11)= [0.276393202250021030359082633127_kreal,&
						-0.850650808352039932181540497063_kreal,&
						-0.447213595499957939281834733746_kreal]
	aParticles%x(:,12)= [0.0_kreal,0.0_kreal,-1.0_kreal]
	do j = 1, 12
		aParticles%x(:,j) = aParticles%x(:,j) / sqrt( sum( aParticles%x(:,j) * aParticles%x(:,j) ) ) * EARTH_RADIUS
		aParticles%x0(:,j) = aParticles%x(:,j)
	enddo
		
	!	Edges
	anEdges%N = 30
	anEdges%verts(:,1) = [1,2]
	anEdges%verts(:,2) = [2,3]
	anEdges%verts(:,3) = [3,1]
	anEdges%verts(:,4) = [3,4]
	anEdges%verts(:,5) = [4,1]
	anEdges%verts(:,6) = [4,5]
	anEdges%verts(:,7) = [5,1]
	anEdges%verts(:,8) = [5,6]
	anEdges%verts(:,9) = [6,1]
	anEdges%verts(:,10)= [6,2]
	anEdges%verts(:,11)= [2,8]
	anEdges%verts(:,12)= [8,3]
	anEdges%verts(:,13)= [8,9]
	anEdges%verts(:,14)= [9,3]
	anEdges%verts(:,15)= [9,4]
	anEdges%verts(:,16)= [9,10]
	anEdges%verts(:,17)= [10,4]
	anEdges%verts(:,18)= [10,5]
	anEdges%verts(:,19)= [10,11]
	anEdges%verts(:,20)= [11,5]
	anEdges%verts(:,21)= [11,6]
	anEdges%verts(:,22)= [11,7]
	anEdges%verts(:,23)= [7,6]
	anEdges%verts(:,24)= [7,2]
	anEdges%verts(:,25)= [7,8]
	anEdges%verts(:,26)= [8,12]
	anEdges%verts(:,27)= [12,9]
	anEdges%verts(:,28)= [12,10]
	anEdges%verts(:,29)= [12,11]
	anEdges%verts(:,30)= [12,7]
	anEdges%leftPanel(1) = 1
	anEdges%leftPanel(2) = 1
	anEdges%leftPanel(3) = 1
	anEdges%leftPanel(4) = 2
	anEdges%leftPanel(5) = 2
	anEdges%leftPanel(6) = 3
	anEdges%leftPanel(7) = 3
	anEdges%leftPanel(8) = 4
	anEdges%leftPanel(9) = 4
	anEdges%leftPanel(10)= 5
	anEdges%leftPanel(11)= 6
	anEdges%leftPanel(12)= 6
	anEdges%leftPanel(13)= 7
	anEdges%leftPanel(14)= 7
	anEdges%leftPanel(15)= 8
	anEdges%leftPanel(16)= 9
	anEdges%leftPanel(17)= 9
	anEdges%leftPanel(18)= 10
	anEdges%leftPanel(19)= 11
	anEdges%leftPanel(20)= 11
	anEdges%leftPanel(21)= 12
	anEdges%leftPanel(22)= 13
	anEdges%leftPanel(23)= 13
	anEdges%leftPanel(24)= 14
	anEdges%leftPanel(25)= 15
	anEdges%leftPanel(26)= 16
	anEdges%leftPanel(27)= 16
	anEdges%leftPanel(28)= 17
	anEdges%leftPanel(29)= 18
	anEdges%leftPanel(30)= 19
	anEdges%rightPanel(1) = 5
	anEdges%rightPanel(2) = 6
	anEdges%rightPanel(3) = 2
	anEdges%rightPanel(4) = 8
	anEdges%rightPanel(5) = 3
	anEdges%rightPanel(6) = 10
	anEdges%rightPanel(7) = 4
	anEdges%rightPanel(8) = 12
	anEdges%rightPanel(9) = 5
	anEdges%rightPanel(10)= 14
	anEdges%rightPanel(11)= 15
	anEdges%rightPanel(12)= 7
	anEdges%rightPanel(13)= 16
	anEdges%rightPanel(14)= 8
	anEdges%rightPanel(15)= 9
	anEdges%rightPanel(16)= 17
	anEdges%rightPanel(17)= 10
	anEdges%rightPanel(18)= 11
	anEdges%rightPanel(19)= 18
	anEdges%rightPanel(20)= 12
	anEdges%rightPanel(21)= 13
	anEdges%rightPanel(22)= 19
	anEdges%rightPanel(23)= 14
	anEdges%rightPanel(24)= 15
	anEdges%rightPanel(25)= 20
	anEdges%rightPanel(26)= 20
	anEdges%rightPanel(27)= 17
	anEdges%rightPanel(28)= 18
	anEdges%rightPanel(29)= 19
	anEdges%rightPanel(30)= 20
	
	! Panels
	aPanels%N = 20
	aPanels%N_Active = 20
	aPanels%edges(:,1) = [1,2,3]
	aPanels%edges(:,2) = [3,4,5]
	aPanels%edges(:,3) = [5,6,7]
	aPanels%edges(:,4) = [7,8,9]
	aPanels%edges(:,5) = [9,10,1]
	aPanels%edges(:,6) = [12,2,11]
	aPanels%edges(:,7) = [12,13,14]
	aPanels%edges(:,8) = [15,4,14]
	aPanels%edges(:,9) = [15,16,17]
	aPanels%edges(:,10)= [18,6,17]
	aPanels%edges(:,11)= [18,19,20]
	aPanels%edges(:,12)= [21,8,20]
	aPanels%edges(:,13)= [21,22,23]
	aPanels%edges(:,14)= [24,10,23]
	aPanels%edges(:,15)= [24,25,11]
	aPanels%edges(:,16)= [27,13,26]
	aPanels%edges(:,17)= [28,16,27]
	aPanels%edges(:,18)= [29,19,28]
	aPanels%edges(:,19)= [30,22,29]
	aPanels%edges(:,20)= [26,25,30]
	aPanels%vertices(:,1) = [1,2,3]
	aPanels%vertices(:,2) = [1,3,4]
	aPanels%vertices(:,3) = [1,4,5]
	aPanels%vertices(:,4) = [1,5,6]
	aPanels%vertices(:,5) = [1,6,2]
	aPanels%vertices(:,6) = [8,3,2]
	aPanels%vertices(:,7) = [3,8,9]
	aPanels%vertices(:,8) = [9,4,3]
	aPanels%vertices(:,9) = [4,9,10]
	aPanels%vertices(:,10)= [10,5,4]
	aPanels%vertices(:,11)= [5,10,11]
	aPanels%vertices(:,12)= [11,6,5]
	aPanels%vertices(:,13)= [6,11,7]
	aPanels%vertices(:,14)= [7,2,6]
	aPanels%vertices(:,15)= [2,7,8]
	aPanels%vertices(:,16)= [12,9,8]
	aPanels%vertices(:,17)= [12,10,9]
	aPanels%vertices(:,18)= [12,11,10]
	aPanels%vertices(:,19)= [12,7,11]
	aPanels%vertices(:,20)= [12,8,7]
	
	do j = 1, 20
		aPanels%x(:,j) = SphereTriCenter( aParticles%x(:, aPanels%vertices(j,1) ), &
										  aParticles%x(:, aPanels%vertices(j,2) ), &
										  aParticles%x(:, aPanels%vertices(j,3) ) )
		aPanels%x0(:,j) = aPanels%x(:,j)
		do k = 1, 3
			aPanels%area(j) = aPanels%area(j) + SphereTriArea( aParticles%x(:,aPanels%vertices(j, k)), &
															   aPanels%x(:, j), &
															   aParticles%x(:,aPanels%vertices(j, mod(k+1,3) + 1)) ) 
		enddo
	enddo
	aPanels%nest(1:20) = 0
	
	! Divide root cube until desired initial nest level
	if ( initNest > 0 ) then
		startIndex = 1
		do k=1,initNest
			nOldPanels = aPanels%N
			do j=startIndex,nOldPanels
				if (.NOT. aPanels%hasChildren(j) ) then
					call DivideQuadPanel(self,j)
				endif
			enddo
			startIndex = nOldPanels
		enddo
	endif
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'... uniform quadrilateral mesh ready.')
end subroutine

subroutine DividePanel(self,panelIndex)
! Interface to call the appropriate DividePanel subroutine.
	type(SphereMesh), intent(inout) :: self
	integer(kint), intent(in) :: panelIndex

	if (self%panelKind == QUAD_PANEL ) then
		call DivideQuadPanel(self,panelIndex)
	elseif ( self%panelKind == TRI_PANEL ) then
		call DivideTriPanel(self, panelIndex)	
	else
		call LogMessage(log,ERROR_LOGGING_LEVEL,'DividePanel ERROR : ','invalid panelKind.')
	endif

end subroutine


subroutine DivideQuadPanel(self,panelIndex)
!	Divides a quadrilateral panel from a cubed sphere mesh into 4 subpanels.
!	Vorticity and Tracer data must be assigned separately.
!	New particles and panels are added to the end of the SphereMesh data arrays.
	! Calling parameters
	type(SphereMesh), intent(inout) :: self
	integer(kint), intent(in) :: panelIndex
	! Local variables
	type(Particles), pointer :: aParticles
	type(Edges), pointer :: anEdges
	type(Panels), pointer :: aPanels
	integer(kint) :: vertexIndices(4), edgeIndices(4), nestLevel
	integer(kint) :: nParticles, nPanels, nEdges
	integer(kint) :: j, childEdges(2)
	logical(klog) :: edgeOrientation(4), alreadyDivided(4)


	aParticles => self%particles
	anEdges => self%edges
	aPanels => self%Panels

	! Error checking
	if ( aPanels%N_Max - aPanels%N < 4) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'DivideQuadPanel ERROR : not enough memory.')
		return
	endif
	if ( panelIndex <= 0 .OR. panelIndex > aPanels%N ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'DivideQuadPanel ERROR : invalid panelIndex.')
		return
	endif
	if ( aPanels%hasChildren(panelIndex)) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'DivideQuadPanel ERROR : panel has already been divided.')
		return
	endif

	! Current mesh state
	vertexIndices = aPanels%vertices(:,panelIndex)
	edgeIndices = aPanels%edges(:,panelIndex)
	nestLevel = aPanels%nest(panelIndex)
	nParticles = aParticles%N
	nPanels = aPanels%N
	nEdges = anEdges%N

	edgeOrientation = .FALSE.
	alreadyDivided = .FALSE.
	do j=1,4
		if ( anEdges%hasChildren(edgeIndices(j)) )  alreadyDivided(j) = .TRUE.
		if ( anEdges%leftPanel(edgeIndices(j)) == panelIndex ) edgeOrientation(j) = .TRUE.
	enddo

	! Panel will be divided from outside to inside, starting at boundaries of parent panel
	! 	New subpanels connect to parent vertices
	aPanels%vertices(1,nPanels+1) = vertexIndices(1)
	aPanels%vertices(2,nPanels+2) = vertexIndices(2)
	aPanels%vertices(3,nPanels+3) = vertexIndices(3)
	aPanels%vertices(4,nPanels+4) = vertexIndices(4)
	!
	!	Parent edge 1, subpanels nPanels + 1 and nPanels + 2
	!
	if ( alreadyDivided(1) ) then
		childEdges = anEdges%children(:,edgeIndices(1))
		if ( edgeOrientation(1) ) then
			aPanels%edges(1,nPanels+1) = childEdges(1)
			aPanels%vertices(2,nPanels+1) = anEdges%verts(2,childEdges(1))
			anEdges%leftPanel(childEdges(1)) = nPanels+1

			aPanels%edges(1,nPanels+2) = childEdges(2)
			aPanels%vertices(1,nPanels+2) = anEdges%verts(2,childEdges(1))
			anEdges%leftPanel(childEdges(2)) = nPanels+2
		else
			aPanels%edges(1,nPanels+1) = childEdges(2)
			aPanels%vertices(2,nPanels+1) = anEdges%verts(2,childEdges(1))
			anEdges%rightPanel(childEdges(2)) = nPanels+1

			aPanels%edges(1,nPanels+2) = childEdges(1)
			aPanels%vertices(1,nPanels+2) = anEdges%verts(2,childEdges(1))
			anEdges%rightPanel(childEdges(1)) = nPanels+2
		endif
	else ! divide edge
		anEdges%hasChildren(edgeIndices(1)) = .TRUE.
		anEdges%children(:,edgeIndices(1)) = [nEdges+1,nEdges+2]

		! create new particle
		aParticles%x(:,nParticles+1) = SphereMidpoint(aParticles%x(:,vertexIndices(1)),&
													  aParticles%x(:,vertexIndices(2)))
		aParticles%x0(:,nParticles+1) = SphereMidpoint(aParticles%x0(:,vertexIndices(1)),&
											  aParticles%x0(:,vertexIndices(2)))
		aPanels%vertices(1,nPanels+2) = nParticles+1
		aPanels%vertices(2,nPanels+1) = nParticles+1
		! create new edges
		if ( edgeOrientation(1) ) then
			anEdges%verts(:,nEdges+1) = [vertexIndices(1),nParticles+1]
			anEdges%leftPanel(nEdges+1) = nPanels+1
			anEdges%rightPanel(nEdges+1) = anEdges%rightPanel(edgeIndices(1))

			aPanels%edges(1,nPanels+1) = nEdges+1

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(2)]
			anEdges%leftPanel(nEdges+2) = nPanels+2
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeIndices(1))

			aPanels%edges(1,nPanels+2) = nEdges+2
		else
			anEdges%verts(:,nEdges+1) = [vertexIndices(2),nParticles+1]
			anEdges%leftPanel(nEdges+1) = anEdges%leftPanel(edgeIndices(1))
			anEdges%rightPanel(nEdges+1) = nPanels+2

			aPanels%edges(1,nPanels+1) = nEdges+2

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(1)]
			anEdges%leftPanel(nEdges+2) = anEdges%leftPanel(edgeIndices(1))
			anEdges%rightPanel(nEdges+2) = nPanels+1

			aPanels%edges(1,nPanels+2) = nEdges+1
		endif

		! update insertion points
		nParticles = nParticles + 1
		nEdges = nEdges + 2
	endif
	!
	!	Parent edge 2
	!
	if ( alreadyDivided(2) ) then
		! connect to edge children
		childEdges = anEdges%children(:,edgeIndices(2))
		aPanels%vertices(3,nPanels+2) = anEdges%verts(2,childEdges(1))
		aPanels%vertices(2,nPanels+3) = anEdges%verts(2,childEdges(1))
		if ( edgeOrientation(2) ) then
			aPanels%edges(2,nPanels+2) = childEdges(1)
			anEdges%leftPanel(childEdges(1)) = nPanels+2

			aPanels%edges(2,nPanels+3) = childEdges(2)
			anEdges%leftPanel(childEdges(2)) = nPanels+3
		else
			aPanels%edges(2,nPanels+2) = childEdges(2)
			anEdges%rightPanel(childEdges(2)) = nPanels+2

			aPanels%edges(2,nPanels+3) = childEdges(1)
			anEdges%rightPanel(childEdges(1)) = nPanels+3
		endif
	else ! divide edge 2
		anEdges%hasChildren(edgeIndices(2)) = .TRUE.
		anEdges%children(:,edgeIndices(2)) = [nEdges+1,nEdges+2]

		! create new particle
		aParticles%x(:,nParticles+1) = SphereMidpoint(aParticles%x(:,vertexIndices(2)),&
													  aParticles%x(:,vertexIndices(3)))
		aParticles%x0(:,nParticles+1) = SphereMidpoint(aParticles%x0(:,vertexIndices(2)),&
													  aParticles%x0(:,vertexIndices(3)))
		aPanels%vertices(3,nPanels+2) = nParticles+1
		aPanels%vertices(2,nPanels+3) = nparticles+1
		! create new edges
		if ( edgeOrientation(2) ) then
			anEdges%verts(:,nEdges+1) = [vertexIndices(2),nParticles+1]
			anEdges%leftPanel(nEdges+1) = nPanels+2
			anEdges%rightPanel(nEdges+1) = anEdges%rightPanel(edgeIndices(2))
			aPanels%edges(2,nPanels+2) = nEdges+1

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(3)]
			anEdges%leftPanel(nEdges+2) = nPanels+3
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeIndices(2))
			aPanels%edges(2,nPanels+3) = nEdges+2
		else
			anEdges%verts(:,nEdges+1) = [vertexIndices(3),nParticles+1]
			anEdges%leftPanel(nEdges+1) = anEdges%leftPanel(edgeIndices(2))
			anEdges%rightPanel(nEdges+1) = nPanels+3
			aPanels%edges(2,npanels+3) = nEdges+1

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(2)]
			anEdges%leftPanel(nEdges+2) = anEdges%leftPanel(edgeIndices(2))
			anEdges%rightPanel(nEdges+2) = nPanels+2
			aPanels%edges(2,nPanels+2) = nEdges+2
		endif

		! update insertion points
		nParticles = nParticles + 1
		nEdges = nEdges + 2
	endif
	!
	!	Parent edge 3
	!
	if ( alreadyDivided(3) ) then
		childEdges = anEdges%children(:,edgeIndices(3))
		aPanels%vertices(4,nPanels+3) = anEdges%verts(2,childEdges(1))
		aPanels%vertices(3,nPanels+4) = anEdges%verts(2,childEdges(1))

		if ( edgeOrientation(3) ) then
			aPanels%edges(3,nPanels+3) = childEdges(1)
			anEdges%leftPanel(childEdges(1)) = nPanels+3

			aPanels%edges(3,nPanels+4) = childEdges(2)
			anEdges%leftPanel(childedges(2)) = nPanels+4
		else
			aPanels%edges(3,nPanels+3) = childEdges(2)
			anEdges%rightPanel(childEdges(2)) = nPanels+3

			aPanels%edges(3,nPanels+4) = childEdges(1)
			anEdges%rightPanel(childEdges(1)) = nPanels+4
		endif
	else ! divide edge 3
		anEdges%hasChildren(edgeIndices(3)) = .TRUE.
		anEdges%children(:,edgeIndices(3)) = [nEdges+1,nEdges+2]

		aParticles%x(:,nParticles+1) = SphereMidpoint(aParticles%x(:,vertexIndices(3)),&
													  aParticles%x(:,vertexIndices(4)))
		aParticles%x0(:,nParticles+1) = SphereMidpoint(aParticles%x0(:,vertexIndices(3)),&
													  aParticles%x0(:,vertexIndices(4)))
		aPanels%vertices(4,nPanels+3) = nParticles+1
		aPanels%vertices(3,nPanels+4) = nParticles+1
		if ( edgeOrientation(3) ) then
			anEdges%verts(:,nEdges+1) = [vertexIndices(3),nParticles+1]
			anEdges%leftPanel(nEdges+1) = nPanels+3
			anEdges%rightPanel(nEdges+1) = anEdges%rightPanel(edgeIndices(3))
			aPanels%edges(3,nPanels+3) = nEdges+1

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(4)]
			anEdges%leftPanel(nEdges+2) = nPanels+4
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeIndices(3))
			aPanels%edges(3,nPanels+4) = nEdges+2

		else
			anEdges%verts(:,nEdges+1) = [vertexIndices(4),nParticles+1]
			anEdges%leftPanel(nEdges+1) = anEdges%leftPanel(edgeIndices(3))
			anEdges%rightPanel(nEdges+1) = nPanels+4
			aPanels%edges(3,nPanels+4) = nEdges+1

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(3)]
			anEdges%leftPanel(nEdges+2) = anEdges%leftPanel(edgeIndices(3))
			anEdges%rightPanel(nEdges+2) = nPanels+3
			aPanels%edges(3,nPanels+3) = nEdges+2

		endif

		! update insertion points
		nParticles = nParticles + 1
		nEdges = nEdges + 2
	endif
	!
	!	Parent edge 4
	!
	if ( alreadyDivided(4) ) then
		childEdges = anEdges%children(:,edgeIndices(4))
		aPanels%vertices(1,nPanels+4) = anEdges%verts(2,childEdges(1))
		aPanels%vertices(4,nPanels+1) = anEdges%verts(2,childEdges(1))
		if ( edgeOrientation(4) ) then
			aPanels%edges(4,nPanels+4) = childEdges(1)
			anEdges%leftPanel(childEdges(1)) = nPanels+4

			aPanels%edges(4,nPanels+1) = childEdges(2)
			anEdges%leftPanel(childEdges(2)) = nPanels+1
		else
			aPanels%edges(4,nPanels+4) = childEdges(2)
			anEdges%rightPanel(childEdges(2)) = nPanels+4

			aPanels%edges(4,nPanels+1) = childEdges(1)
			anEdges%rightPanel(childEdges(1)) = nPanels+1
		endif
	else ! divide edge 4
		anEdges%hasChildren(edgeIndices(4)) = .TRUE.
		anEdges%children(:,edgeIndices(4)) = [nEdges+1,nEdges+2]

		aParticles%x(:,nParticles+1) = SphereMidpoint(aParticles%x(:,vertexIndices(4)),&
													  aParticles%x(:,vertexIndices(1)))
		aParticles%x0(:,nParticles+1) = SphereMidpoint(aParticles%x0(:,vertexIndices(4)),&
													  aParticles%x0(:,vertexIndices(1)))
		aPanels%vertices(1,nPanels+4) = nParticles+1
		aPanels%vertices(4,nPanels+1) = nParticles+1
		if ( edgeOrientation(4) ) then
			anEdges%verts(:,nEdges+1) = [vertexIndices(4),nParticles+1]
			anEdges%leftPanel(nEdges+1) = nPanels+4
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeIndices(4))
			aPanels%edges(4,nPanels+4) = nEdges+1

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(1)]
			anEdges%leftPanel(nEdges+2) = npanels+1
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeIndices(4))
			aPanels%edges(4,nPanels+1) = nEdges+2

		else
			anEdges%verts(:,nEdges+1) = [vertexIndices(1),nParticles+1]
			anEdges%leftPanel(nEdges+1) = anEdges%leftPanel(edgeIndices(4))
			anEdges%rightPanel(nEdges+1) = nPanels+1
			aPanels%edges(4,nPanels+1) = nEdges+1

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(4)]
			anEdges%leftPanel(nEdges+2) = anEdges%leftPanel(edgeIndices(4))
			anEdges%rightPanel(nEdges+2) = nPanels+4
			aPanels%edges(4,nPanels+4) = nEdges+2
		endif

		! update insertion points
		nParticles = nParticles + 1
		nEdges = nEdges + 2
	endif
	!
	!	Center particle
	!
	aParticles%x(:,nParticles+1) = aPanels%x(:,panelIndex)
	aParticles%x0(:,nParticles+1) = aPanels%x0(:,panelIndex)


	aPanels%vertices(3,nPanels+1) = nParticles+1
	aPanels%vertices(4,nPanels+2) = nParticles+1
	aPanels%vertices(1,nPanels+3) = nParticles+1
	aPanels%vertices(2,nPanels+4) = nParticles+1
	!
	!	Interior edges
	!
	anEdges%verts(:,nEdges+1) = [aPanels%vertices(2,nPanels+1),nParticles+1]
	anEdges%leftPanel(nEdges+1) = nPanels+1
	anEdges%rightPanel(nEdges+1) = nPanels+2

	anEdges%verts(:,nEdges+2) = [nParticles+1,aPanels%vertices(3,nPanels+4)]
	anEdges%leftPanel(nEdges+2) = nPanels+4
	anEdges%rightPanel(nEdges+2) = nPanels+3

	anEdges%verts(:,nEdges+3) = [aPanels%vertices(3,nPanels+2),nParticles+1]
	anEdges%leftPanel(nEdges+3) = nPanels+2
	anEdges%rightPanel(nEdges+3) = nPanels+3
	anEdges%verts(:,nEdges+4) = [nParticles+1,apanels%vertices(4,nPanels+1)]
	anEdges%leftPanel(nEdges+4) = nPanels+1
	anEdges%rightPanel(nEdges+4) = nPanels+4

	aPanels%edges(2,nPanels+1) = nEdges+1
	aPanels%edges(4,nPanels+2) = nEdges+1
	aPanels%edges(2,nPanels+4) = nEdges+2
	aPanels%edges(4,nPanels+3) = nEdges+2
	aPanels%edges(3,nPanels+2) = nEdges+3
	aPanels%edges(1,nPanels+3) = nEdges+3
	aPanels%edges(3,nPanels+1) = nEdges+4
	aPanels%edges(1,nPanels+4) = nEdges+4

	nParticles = nParticles+1
	nEdges = nEdges+4

	!
	!	Update data structures
	!
	aPanels%area(panelIndex) = 0.0_kreal
	aPanels%hasChildren(panelIndex) = .TRUE.
	aPanels%children(:,panelIndex) = [1,2,3,4] + nPanels
	aPanels%N = nPanels+4
	aPanels%N_Active = aPanels%N_Active + 3

	anEdges%N = nEdges

	aParticles%N = nParticles

	!
	!	New subpanels : centers, areas, and nest level
	!
	do j=1,4
		aPanels%x(:,nPanels+j) = SphereQuadCenter(aParticles%x(:,aPanels%vertices(1,nPanels+j)),&
												  aParticles%x(:,aPanels%vertices(2,nPanels+j)),&
												  aParticles%x(:,aPanels%vertices(3,nPanels+j)),&
												  aParticles%x(:,aPanels%vertices(4,nPanels+j)))
		aPanels%x0(:,nPanels+j) = SphereQuadCenter(aParticles%x0(:,aPanels%vertices(1,nPanels+j)),&
												  aParticles%x0(:,aPanels%vertices(2,nPanels+j)),&
												  aParticles%x0(:,aPanels%vertices(3,nPanels+j)),&
												  aParticles%x0(:,aPanels%vertices(4,nPanels+j)))
		aPanels%area(nPanels+j) = QuadPanelArea(self,nPanels+j)
		aPanels%nest(nPanels+j) = nestLevel+1
	enddo

end subroutine

subroutine DivideTriPanel(self, panelIndex) 
	type(SphereMesh), intent(inout) :: self
	integer(kint), intent(in) :: panelIndex
	!
	type(Particles), pointer :: aParticles
	type(Edges), pointer :: anEdges
	type(Panels), pointer :: aPanels
	logical(klog) :: edgeOrientation(3), alreadyDivided(3)
	
	!
	if ( aPanels%N_Max - aPanels%N < 4 ) then
		call LogMessage(log, WARNING_LOGGING_LEVEL, logkey, 'DivideTriPanel WARNING : ', 'not enough memory.')
		return
	endif
	if ( panelIndex <=0 .OR. panelIndex > aPanels%N ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, 'DivideTriPanel ERROR : ', ' invalid panel index.')
		return
	endif
	if ( aPanels%hasChildren(panelIndex) ) then
		call LogMessage(log, WARNING_LOGGING_LEVEL, logKey, 'DivideTriPanel WARNING : ', ' panel already has children.')
		return
	endif
	
	do j = 1, 3
		if ( anEdges%hasChildren( aPanels%edges(j, panelIndex) ) ) alreadyDivided(j) = .TRUE.
		if ( anEdges%leftPanel( aPanels%edges(j, panelIndex) ) == panelIndex ) edgeOrientation(j) = .TRUE.
	enddo
	
end subroutine 


function QuadPanelArea(self,panelIndex)
! Determines the spherical area contained within a quadrilateral panel from a cubed sphere mesh.
	! calling parameters
	real(kreal) :: QuadPanelArea
	type(SphereMesh), intent(in) :: self
	integer(kint), intent(in) :: panelIndex
	! local variables
	type(Particles), pointer :: aParticles
	type(Edges), pointer :: anEdges
	type(Panels), pointer :: aPanels
	integer(kint) :: leafEdges(8), ccwVerts(8), j, nLeafEdges
	real(kreal) :: centerX(3)

	aParticles => self%particles
	anEdges => self%edges
	aPanels => self%panels

	! Error checking
	if ( panelIndex > aPanels%N) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'QuadPanelArea ERROR : out of bounds.')
		return
	endif
	centerX = aPanels%x(:,panelIndex)

	if ( aPanels%hasChildren(panelIndex)) then
		write(logString,'(A,I4,A)') 'panel ',panelIndex,' has children; area = 0.'
		call LogMessage(log,ERROR_LOGGING_LEVEL,'AREA ERROR : ',trim(logString))
	endif

	! Find leaf edges for each of the parent panel edges (in case adjacent panels are not at same resolution)
	call CCWEdgesAndParticlesAroundPanel(leafEdges,ccwVerts,nleafEdges,self,panelIndex)

	QuadPanelArea = 0.0_kreal
	do j=1,nLeafEdges-1
		QuadPanelArea = QuadPanelArea + SphereTriArea(aParticles%x(:,ccwVerts(j)),centerX,aParticles%x(:,ccwVerts(j+1)))
	enddo
	QuadPanelArea = QuadPanelArea + SphereTriArea(aParticles%x(:,ccwVerts(nLeafEdges)),centerX,aParticles%x(:,ccwVerts(1)))
!	DEBUG
!	if (QuadPanelArea == 0.0_kreal .OR. isnan(QuadPanelArea) ) then
!		call StartSection(log,'PanelEdges')
!			call LogMessage(log,DEBUG_LOGGING_LEVEL,'Panel = ',panelIndex)
!			write(formatString,'(A,I4,A)') '(',size(leafEdges),'I6)'
!			write(logString,formatString) leafEdges
!			call LogMessage(log,DEBUG_LOGGING_LEVEL,'leafEdges = ',trim(logString))
!		call EndSection(log)
!	endif
end function

!subroutine SetEdgeLengths(self)
!	type(SphereMesh), intent(inout) :: self
!	integer(kint) :: j
!	type(Edges), pointer :: anEdges
!	type(Particles), pointer :: aParticles
!
!	anEdges => self%edges
!	aParticles => self%particles
!
!	do j=1,self%edges%N
!		if ( anEdges%hasChildren(j) ) then
!			anEdges%length(j) = 0.0_kreal
!			anEdges%length0(j) = 0.0_kreal
!		else
!			anEdges%length(j) = SphereDistance(aparticles%x(:,anEdges%verts(1,j)),aParticles%x(:,anEdges%verts(2,j)))
!			anEdges%length0(j) = SphereDistance(aparticles%x0(:,anEdges%verts(1,j)),aParticles%x0(:,anEdges%verts(2,j)))
!		endif
!	enddo
!end subroutine


function FindClosestRootPanel(self,xyz)
!	Finds the nearest panel of the original cube's 6 faces closest to the point xyz.
!	Used to initialize a tree search.
	integer(kint) :: FindClosestRootPanel
	type(SphereMesh), intent(in) :: self
	real(kreal), intent(in) :: xyz(3)
	integer(kint) :: panelKind, j
	type(Panels), pointer :: aPanels
	real(kreal) :: currentMin, testDist
	aPanels=>self%panels
	panelKind = GetPanelKind(aPanels)
	FindClosestRootPanel = 0
	currentMin = PI
	if ( panelKind == QUAD_PANEL) then
		do j=1,6
			testDist = SphereDistance(aPanels%x(:,j),xyz)
			if ( testDist < currentMin ) then
				currentMin = testDist
				FindClosestRootPanel = j
			endif
		enddo
	endif
end function


recursive subroutine LocatePointTreeSearch(inPanel,self,xyz,rootPanel)
!	Performs a tree search via panels children (its subpanels) to find the closest centroid of a
!	child of rootPanel to the point xyz.
	integer(kint), intent(out) :: inPanel
	type(SphereMesh), intent(in) :: self
	real(kreal), intent(in) :: xyz(3)
	integer(kint), intent(in) :: rootPanel
	integer(kint) :: j, currentPanel
	real(kreal) :: currentMin, testDist, centroid(3)
	type(Panels), pointer :: aPanels

	aPanels => self%panels

	if ( aPanels%hasChildren(rootPanel) ) then
		currentMin = PI
		do j=1,4
			centroid = PanelCentroid(self,aPanels%children(j,rootPanel))
			testDist = SphereDistance(centroid,xyz)
			if ( testDist < currentMin ) then
				currentMin = testDist
				currentPanel = aPanels%children(j,rootPanel)
			endif
		enddo
		call LocatePointTreeSearch(inPanel,self,xyz,currentPanel)
	else
		inPanel = rootPanel
	endif
end subroutine


recursive subroutine LocatePointWalkSearch(inPanel,self,xyz,startPanel)
!	Performs a walk search to find the closest panel centroid of a SphereMesh to the point xyz.
!	Assuming panels are convex, this panel contains the point xyz.
	integer(kint), intent(out) :: inPanel
	type(SphereMesh), intent(in) :: self
	real(kreal), intent(in) :: xyz(3)
	integer(kint), intent(in) :: startPanel
	real(kreal) :: currentMin, testDist, centroid(3)
	integer(kint) :: neighborPanels(8), nNeighbors, currentPanel, j


	! Error checking
	if ( startPanel < 0 .OR. startPanel > self%panels%N) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'WalkSearch ERROR : ','invalid initial guess.')
		return
	endif
	centroid = PanelCentroid(self,startPanel)
	currentMin = SphereDistance(centroid,xyz)
	currentPanel = startPanel
	call FindAdjacentPanels(neighborPanels,nNeighbors,self,startpanel)

	do j=1,nNeighbors
		centroid = PanelCentroid(self,neighborPanels(j))
		testDist = SphereDistance(centroid,xyz)
		if ( testDist < currentMin ) then
			currentMin = testDist
			currentPanel = neighborPanels(j)
		endif
	enddo
	if ( currentPanel == startPanel) then
		inPanel = currentPanel
		return
	else
		call LocatePointWalkSearch(inPanel,self,xyz,currentPanel)
	endif

end subroutine


subroutine FindAdjacentPanels(adjPanels,nNeighbors,self,panelIndex)
!	Returns a counter-clockwise ordered list of a panel's neighbors.
	integer(kint), intent(inout) :: adjPanels(:), nNeighbors
	type(Spheremesh), intent(in) :: self
	integer(kint), intent(in) :: panelIndex
	integer(kint) :: edgeList(8), vertList(8), nVerts, k
	type(Edges), pointer :: anEdges

	! Error checking
	if ( GetPanelKind(self%panels) == QUAD_PANEL ) then
		if ( size(adjPanels) < 8 ) then
			call LogMessage(log,WARNING_LOGGING_LEVEL,'FindAdjacentPanels WARNING : ',' not enough memory.')
		endif
	endif

	anEdges => self%edges
	call CCWEdgesAndParticlesAroundPanel(edgeList,vertList,nVerts,self,panelIndex)
	nNeighbors = nVerts
	do k=1,nVerts
		if ( anEdges%leftPanel(edgeList(k)) == panelIndex) then
			adjPanels(k) = anEdges%rightPanel(edgeList(k))
		else
			adjPanels(k) = anEdges%leftPanel(edgeList(k))
		endif
	enddo

end subroutine

function PanelCentroid(self,panelIndex)
!	Determines the centroid of a quadrilateral panel with only 4 vertices.
!	NOTE : This function assumes quad panels have only 4 vertices, tri panels have only 3 vertices
!	Its purpose is to facilitate the TreeSearch algorithm, where inexact results are ok, and the
!	WalkSearch algorithm, where every panel is a low-level panel with only 3 or 4 vertices.
!
!	TO DO :
!	-- does not account for AMR --
!
	real(kreal) :: PanelCentroid(3)
	type(SphereMesh), intent(in) :: self
	integer(kint), intent(in) :: panelIndex
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: panelKind

	aParticles => self%particles
	aPanels => self%panels

	panelKind = GetPanelKind(aPanels)

	if ( panelKind == QUAD_PANEL) then
		PanelCentroid = SphereQuadCenter(aParticles%x(:,aPanels%vertices(1,panelIndex)), &
										 aParticles%x(:,aPanels%vertices(2,panelIndex)), &
										 aParticles%x(:,aPanels%vertices(3,panelIndex)), &
										 aParticles%x(:,aPanels%vertices(4,panelIndex)))
	endif
end function



subroutine InitLogger(aLog,rank)
	type(Logger), intent(inout) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine



end module
