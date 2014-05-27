module PlaneMeshModule
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
use PlaneGeomModule
use ParticlesModule
use EdgesModule
use PanelsModule

implicit none

include 'mpif.h'

private
public PlaneMesh
public New, Delete, Copy
public InitializeRectangle
public CCWEdgesAndParticlesAroundPanel, FindAdjacentPanels
public LocatePoint
public LogStats
public TotalArea

!
!----------------
! Types and module constants
!----------------
!
type PlaneMesh
	type(Particles), pointer :: particles => null()
	type(Edges), pointer :: edges => null()
	type(Panels), pointer :: panels => null()
	integer(kint) :: initNest
	integer(kint) :: AMR
	integer(kint) :: nTracer
	integer(kint) :: nRootPanels
end type

integer(kint), parameter :: panelKind = QUAD_PANEL, problemKind = PLANE_SOLVER
!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'PlaneMesh'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=128) :: logString
character(len=24) :: formatString
!
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
	module procedure CopyPrivate
end interface

interface LogStats
	module procedure LogStatsPrivate
end interface


contains
!
!----------------
! Standard methods : Constructor / Destructor
!----------------
!
subroutine NewPrivate(self, initNest, AMR, nTracer)
	type(PlaneMesh), intent(out) :: self
	integer(kint), intent(in) :: initNest, AMR, nTracer
	!
	integer(kint) :: nPanels, nParticles, nEdges

	if ( .NOT. logInit ) call InitLogger(log, procRank)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'allocating memory for new plane mesh...')
	!
	! Allocate the pointers
	!
	allocate(self%particles)
	allocate(self%edges)
	allocate(self%panels)

	self%initNest = initNest
	self%AMR = AMR
	self%nTracer = nTracer

	!
	! TO DO : compute the array sizes
	!
	nPanels = PlanePanelMax(initNest + AMR)
	nParticles = PlaneParticleMax(initNest + AMR)
	nEdges = PlaneEdgeMax(initNest+AMR)

!	nPanels = 10000
!	nEdges = 10000
!	nParticles = 10000

	call New(self%particles,nParticles,panelKind,nTracer,problemKind)
	call New(self%edges, nEdges)
	call New(self%panels, nPanels, panelKind, nTracer, problemKind)

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'new mesh ready for initialization.')
end subroutine

subroutine DeletePrivate(self)
	type(PlaneMesh), intent(inout) :: self

	call Delete(self%panels)
	call Delete(self%edges)
	call Delete(self%particles)

	deallocate(self%particles)
	deallocate(self%edges)
	deallocate(self%panels)

	self%initNest = 0
	self%AMR = 0
	self%nTracer = 0
end subroutine

subroutine CopyPrivate(newMesh, oldMesh)
	type(PlaneMesh), intent(inout) :: newMesh
	type(PlaneMesh), intent(in) :: oldMesh

	newMesh%initNest = oldMesh%initNest
	newMesh%AMR = oldMesh%AMR
	newMesh%nTracer = oldMesh%nTracer
	newMesh%nRootPanels = oldMesh%nRootPanels

	call Copy(newMesh%particles, oldMesh%particles)
	call Copy(newMesh%edges, oldMesh%edges)
	call Copy(newMesh%panels, oldMesh%panels)

end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine LogStatsPrivate(self, aLog, msg)
	type(PlaneMesh), intent(in) :: self
	type(Logger), intent(inout) :: aLog
	character(len=*), intent(in), optional :: msg
	if ( present(msg) ) then
		call LogStats(self%particles,aLog,msg)
		call LogStats(self%edges,aLog,msg)
		call LogStats(self%panels,aLog,msg)
	else
		call LogStats(self%particles,aLog)
		call LogStats(self%edges,aLog)
		call LogStats(self%panels,aLog)
	endif
end subroutine

subroutine InitializeRectangle(self, xmin, xmax, ymin, ymax, boundaryType)
	type(PlaneMesh), intent(inout) :: self
	real(kreal), intent(in) :: xmin, xmax, ymin, ymax
	integer(kint), intent(in) :: boundaryType
	!
	type(Particles), pointer :: aParticles
	type(Edges), pointer :: anEdges
	type(Panels), pointer :: aPanels
	integer(kint) :: j, k, startIndex, nOldPanels
	real(kreal) :: area0

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,' initializing rectangular mesh.')

	if ( boundaryType /= FREE_BOUNDARIES .AND. boundaryType /= PERIODIC_BOUNDARIES ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL,logkey,'ERROR : invalid boundary condition.')
		return
	endif

	aParticles => self%particles
	anEdges => self%edges
	aPanels => self%panels

	!
	! initialize mesh with 4 quadrilateral panels
	!
	aParticles%N = 9
	aParticles%x(:,1) = [ xmin, ymax]
	aParticles%x(:,2) = [ xmin, 0.5_kreal*(ymin + ymax)]
	aParticles%x(:,3) = [ xmin, ymin]
	aParticles%x(:,4) = [ 0.5_kreal*(xmin + xmax), ymin]
	aParticles%x(:,5) = [ xmax, ymin]
	aParticles%x(:,6) = [ xmax, 0.5_kreal*(ymin + ymax)]
	aParticles%x(:,7) = [ xmax, ymax]
	aParticles%x(:,8) = [ 0.5_kreal*(xmin + xmax), ymax]
	aParticles%x(:,9) = [ 0.5_kreal*(xmin + xmax), 0.5_kreal*(ymin + ymax)]
	do j=1,9
		aParticles%x0(:,j) = aParticles%x(:,j)
	enddo

	anEdges%N = 12
	do j=1,7
		anEdges%verts(:,j) = [j,j+1]
	enddo
	anEdges%verts(:,8) = [8,1]
	anEdges%verts(:,9) = [2,9]
	anEdges%verts(:,10) = [9,6]
	anEdges%verts(:,11) = [4,9]
	anEdges%verts(:,12) = [9,8]
	anEdges%leftPanel(1) = 1
	anEdges%leftPanel(2) = 2
	anEdges%leftPanel(3) = 2
	anEdges%leftPanel(4) = 3
	anEdges%leftPanel(5) = 3
	anEdges%leftPanel(6) = 4
	anEdges%leftPanel(7) = 4
	anEdges%leftPanel(8) = 1
	anEdges%leftPanel(9) = 1
	anEdges%leftPanel(10) = 4
	anEdges%leftPanel(11) = 2
	anEdges%leftPanel(12) = 1
	if ( boundaryType == FREE_BOUNDARIES ) then
		anEdges%rightPanel(1) = 0
		anEdges%rightPanel(2) = 0
		anEdges%rightPanel(3) = 0
		anEdges%rightPanel(4) = 0
		anEdges%rightPanel(5) = 0
		anEdges%rightPanel(6) = 0
		anEdges%rightPanel(7) = 0
		anEdges%rightPanel(8) = 0
		anEdges%rightPanel(9) = 2
		anEdges%rightPanel(10) = 3
		anEdges%rightPanel(11) = 3
		anEdges%rightPanel(12) = 4
	elseif (boundaryType == PERIODIC_BOUNDARIES) then
		anEdges%rightPanel(1) = 4
		anEdges%rightPanel(2) = 3
		anEdges%rightPanel(3) = 1
		anEdges%rightPanel(4) = 4
		anEdges%rightPanel(5) = 2
		anEdges%rightPanel(6) = 1
		anEdges%rightPanel(7) = 3
		anEdges%rightPanel(8) = 2
		anEdges%rightPanel(9) = 2
		anEdges%rightPanel(10) = 3
		anEdges%rightPanel(11) = 3
		anEdges%rightPanel(12) = 4
	endif

	self%nRootPanels = 4
	aPanels%N = 4
	aPanels%N_Active = 4
	aPanels%nest(1:4) = 0
	aPanels%edges(:,1) = [1,9,12,8]
	aPanels%edges(:,2) = [2,3,11,9]
	aPanels%edges(:,3) = [11,4,5,10]
	aPanels%edges(:,4) = [12,10,6,7]
	aPanels%vertices(:,1) = [1,2,9,8]
	aPanels%vertices(:,2) = [2,3,4,9]
	aPanels%vertices(:,3) = [9,4,5,6]
	aPanels%vertices(:,4) = [8,9,6,7]
	area0 = (xmax - xmin) * (ymax - ymin)
	do j=1,4
		aPanels%x(:,j) = QuadCentroid( aParticles%x(:,aPanels%vertices(1,j)), aParticles%x(:,aPanels%vertices(2,j)), &
									      aParticles%x(:,aPanels%vertices(3,j)), aParticles%x(:,aPanels%vertices(4,j)) )
		aPanels%x0(:,j) = aPanels%x(:,j)
		aPanels%area(j) = area0/4.0_kreal
	enddo

	!
	! recursively divide panels to desired nest level
	!
	if ( self%initNest > 0 ) then
		startIndex = 1
		call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,' dividing rectangular mesh ...')
		do k=1, self%initNest
			nOldPanels = aPanels%N
			do j=startIndex, nOldPanels
				call DividePanel(self, j)
			enddo
			startIndex = nOldPanels+1
		enddo
	endif
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'... rectangular quadrilateral mesh ready.')
end subroutine

subroutine CCWEdgesAndParticlesAroundPanel(edgeList, vertlist, nVerts, self, panelIndex)
! 	Returns a counter-clockwise ordered list of edges and a counter-clockwise ordered list
! 	of vertices around a panel.  For use with AMR, this subroutine does not assume a panel
!	has only 4 vertices-- panels may have as many as 8.
	integer(kint), intent(inout) :: edgeList(8), vertList(8), nVerts
	type(PlaneMesh), intent(in) :: self
	integer(kint), intent(in) :: panelIndex
	!
	type(Particles), pointer :: aParticles
	type(Edges), pointer :: anEdges
	type(Panels), pointer :: aPanels
	integer(kint) :: k, edgeK, vertK

	aParticles=>self%particles
	anEdges=>self%edges
	aPanels=>self%panels

	!DEBUG
	!call StartSection(log,'entering CCWEdgesAndParticlesAroundPanel')

	if ( size(edgeList) < 8 .OR. size(vertList) < 8 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'CCWEdgesAndParticlesAroundPanel ERROR : ',' output array too small.')
		return
	endif
	if ( aPanels%hasChildren(panelIndex) ) then
		call LogMessage(log,WARNING_LOGGING_LEVEL,logkey//'CCWEdgesAndParticlesAroundPanel ERROR : ','panel is divided.')
		return
	endif
	edgelist = 0
	vertlist = 0
	edgek = 0
	vertk = 0
	do k=1,4
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

	!DEBUG
	!call EndSection(log)
end subroutine


subroutine FindAdjacentPanels(adjPanels, nAdj, self, panelIndex)
	integer(kint), intent(inout) :: adjPanels(:), nAdj
	type(PlaneMesh), intent(in) :: self
	integer(kint), intent(in) :: panelindex
	!
	integer(kint) :: edgeList(8), vertList(8), nVerts, j
	type(Edges), pointer :: anEdges

	call CCWEdgesAndParticlesAroundPanel(edgelist,vertlist, nverts, self, panelindex)

	anEdges => self%edges
	nAdj = nVerts
	do j=1, nVerts
		if ( anEdges%leftPanel(edgelist(j)) == panelIndex) then
			adjpanels(j) = anEdges%rightPanel(edgelist(j))
		elseif ( anEdges%rightPanel(edgeList(j)) == panelIndex) then
			adjpanels(j) = anEdges%leftPanel(edgelist(j))
		else
			write(logstring,'(A, I8, A, I8)') ' found unconnected edge : ', j, ' in panel ', panelIndex
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey//'AdjacentPanel ERROR :', logString)
			return
		endif
	enddo
end subroutine

function LocatePoint(self, xy)
! returns the index of the panel that contains the point xy
	integer(kint) :: LocatePoint
	type(PlaneMesh), intent(in) :: self
	real(kreal), intent(in) :: xy(2)
	!
	integer(kint) :: rootPanel, startWalk

	rootPanel = NearestRootPanel(self,xy)
	call LocatePointTreeSearch(startWalk, self, xy, rootPanel)

	call LocatePointWalkSearch(LocatePoint, self, xy, startWalk)
end function

function TotalArea(self)
	real(kreal) :: TotalArea
	type(PlaneMesh), intent(in) :: self
	!
	type(Panels), pointer :: apanels
	apanels => self%panels
	TotalArea = sum(apanels%area)
end function

!
!----------------
! Module methods : module- or type-specific private functions
!----------------
!

subroutine DividePanel(self, panelIndex)
	type(PlaneMesh), intent(inout) :: self
	integer(kint), intent(in) :: panelIndex
	!
	type(Particles), pointer :: aParticles
	type(Edges), pointer :: anEdges
	type(Panels), pointer :: aPanels
	integer(kint) :: vertIndices(4), edgeIndices(4), nestlevel
	integer(kint) :: nParticles, nPanels, nEdges
	integer(kint) :: j, k, childEdges(2)
	logical(klog) :: edgeOrientation(4), alreadyDivided(4)

	aParticles => self%particles
	anEdges => self%edges
	aPanels => self%Panels

	!
	! Error checking
	!
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

	!
	! get current mesh state
	!
	vertIndices = aPanels%vertices(:,panelIndex)
	edgeIndices = aPanels%edges(:,panelIndex)
	nestLevel = aPanels%nest(panelIndex)
	nParticles = aParticles%N
	nPanels = aPanels%N
	nEdges = anEdges%N
!DEBUG
!	write(logstring,'(A,I4)') 'dividing panel ',panelIndex
!	call StartSection(log,logstring)
!	write(logstring,'(A,4I8)') 'vertices = ', vertIndices
!	call LogMessage(log, DEBUG_LOGGING_LEVEL, '  ', logstring)
!	write(logstring,'(A,4I8)') 'edges    = ', edgeIndices
!	call LogMessage(log,DEBUG_LOGGING_LEVEL,'  ',logstring)



	edgeOrientation = .FALSE.
	alreadyDivided = .FALSE.
	do j=1,4
		if ( anEdges%hasChildren(edgeIndices(j)) )  alreadyDivided(j) = .TRUE.
		if ( anEdges%leftPanel(edgeIndices(j)) == panelIndex ) edgeOrientation(j) = .TRUE.
	enddo

	!
	! connect parent vertices to child panels
	!
	aPanels%vertices(1,nPanels+1) = vertIndices(1)
	aPanels%vertices(2,nPanels+2) = vertIndices(2)
	aPanels%vertices(3,npanels+3) = vertIndices(3)
	aPanels%vertices(4,nPanels+4) = vertIndices(4)
	!
	! Panel will be divided from outside to inside, starting at boundaries of parent panel
	!
	if ( alreadyDivided(1) ) then
		!
		! connect to existing child edges
		!
		childEdges = anEdges%children(:,edgeIndices(1))
		if ( edgeOrientation(1) ) then
			aPanels%edges(1,nPanels+1) = childEdges(1)
			aPanels%edges(1,nPanels+2) = childEdges(2)
			anEdges%leftPanel(childEdges(1)) = nPanels+1
			anEdges%leftPanel(childEdges(2)) = nPanels+2
		else
			aPanels%edges(1,nPanels+1) = childEdges(2)
			aPanels%edges(1,nPanels+2) = childEdges(1)
			anEdges%rightPanel(childEdges(1)) = nPanels+2
			anEdges%rightPanel(childEdges(2)) = nPanels+1
		endif
		aPanels%vertices(2,nPanels+1) = anEdges%verts(2,childEdges(1))
		aPanels%vertices(1,nPanels+2) = anEdges%verts(2,childEdges(1))
	else
		!
		! divide parent edge 1
		!
		anEdges%hasChildren(edgeIndices(1)) = .TRUE.
		anEdges%children(:,edgeIndices(1)) = [nEdges+1,nEdges+2]
		if ( edgeOrientation(1) ) then
			aPanels%edges(1,nPanels+1) = nEdges+1
			aPanels%edges(1,nPanels+2) = nEdges+2

			anEdges%verts(:,nEdges+1) = [vertIndices(1),nParticles+1]
			anEdges%leftPanel(nEdges+1) = nPanels+1
			anEdges%rightPanel(nEdges+1) = anEdges%rightPanel(edgeIndices(1))

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertIndices(2)]
			anEdges%leftPanel(nEdges+2) = nPanels+2
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeIndices(1))
		else
			aPanels%edges(1,nPanels+1) = nEdges+2
			aPanels%edges(1,nPanels+2) = nEdges+1

			anEdges%verts(:,nEdges+1) = [vertIndices(2),nParticles+1]
			anEdges%leftPanel(nEdges+1) = anEdges%leftPanel(edgeIndices(1))
			anEdges%rightPanel(nEdges+1) = nPanels+2

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertIndices(1)]
			anEdges%leftPanel(nEdges+2) = anEdges%leftPanel(edgeIndices(1))
			anEdges%rightPanel(nEdges+2) = nPanels+1
		endif

		aParticles%x(:,nParticles+1) = Midpoint(aParticles%x(:,vertIndices(1)), aParticles%x(:,vertIndices(2)))
		aParticles%x0(:,nParticles+1) = Midpoint(aParticles%x0(:,vertIndices(1)), aParticles%x0(:,vertIndices(2)))

		aPanels%vertices(2,nPanels+1) = nParticles+1
		aPanels%vertices(1,nPanels+2) = nParticles+1

		nParticles = nParticles + 1
		nEdges = nEdges + 2
	endif
	if ( alreadyDivided(2) ) then
		!
		! connect to existing child edges
		!
		childEdges = anEdges%children(:,edgeIndices(2))
		if ( edgeOrientation(2) ) then
			aPanels%edges(2,nPanels+2) = childEdges(1)
			aPanels%edges(2,nPanels+3) = childEdges(2)
			anEdges%leftPanel(childEdges(1)) = nPanels+2
			anEdges%leftPanel(childEdges(2)) = nPanels+3
		else
			aPanels%edges(2,nPanels+2) = childEdges(2)
			aPanels%edges(2,nPanels+3) = childEdges(1)
			anEdges%rightPanel(childEdges(1)) = nPanels+3
			anEdges%rightPanel(childEdges(2)) = nPanels+2
		endif
		aPanels%vertices(3,nPanels+2) = anEdges%verts(2,childEdges(1))
		aPanels%vertices(2,nPanels+3) = anEdges%verts(2,childEdges(1))
	else
		!
		! divide parent edge 2
		!
		anEdges%hasChildren(edgeIndices(2)) = .TRUE.
		anEdges%children(:,edgeIndices(2)) = [nEdges+1,nEdges+2]
		if ( edgeOrientation(2) ) then
			aPanels%edges(2,nPanels+2) = nEdges+1
			aPanels%edges(2,nPanels+3) = nEdges+2

			anEdges%verts(:,nEdges+1) = [vertIndices(2),nParticles+1]
			anEdges%leftPanel(nEdges+1) = nPanels+2
			anEdges%rightPanel(nEdges+1) = anEdges%rightPanel(edgeIndices(2))

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertIndices(3)]
			anEdges%leftPanel(nEdges+2) = nPanels+3
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeIndices(2))
		else
			aPanels%edges(2,nPanels+2) = nEdges+2
			aPanels%edges(2,nPanels+3) = nEdges+1

			anEdges%verts(:,nEdges+1) = [vertIndices(3),nParticles+1]
			anEdges%leftPanel(nEdges+1) = anEdges%leftPanel(edgeIndices(2))
			anEdges%rightPanel(nEdges+1) = nPanels+3

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertIndices(2)]
			anEdges%leftPanel(nEdges+2) = anEdges%leftPanel(edgeIndices(2))
			anEdges%rightPanel(nEdges+2) = nPanels+2
		endif

		aParticles%x(:,nParticles+1) = Midpoint(aParticles%x(:,vertIndices(2)),aParticles%x(:,vertIndices(3)))
		aParticles%x0(:,nParticles+1) = Midpoint(aParticles%x0(:,vertIndices(2)),aParticles%x0(:,vertIndices(3)))

		aPanels%vertices(3,nPanels+2) = nParticles+1
		aPanels%vertices(2,nPanels+3) = nParticles+1

		nParticles = nParticles+1
		nEdges = nEdges+2
	endif
	if ( alreadyDivided(3) ) then
		!
		! connect to existing child edges
		!
		childEdges = anEdges%children(:,edgeIndices(3))
		if ( edgeOrientation(3) ) then
			aPanels%edges(3,nPanels+3) = childEdges(1)
			aPanels%edges(3,nPanels+4) = childEdges(2)
			anEdges%leftPanel(childEdges(1)) = nPanels+3
			anEdges%leftPanel(childEdges(2)) = nPanels+4
		else
			aPanels%edges(3,nPanels+3) = childEdges(2)
			aPanels%edges(3,nPanels+4) = childEdges(1)
			anEdges%rightPanel(childEdges(1)) = nPanels+4
			anEdges%rightPanel(childEdges(2)) = nPanels+3
		endif
		aPanels%vertices(4,nPanels+3) = anEdges%verts(2,childEdges(1))
		aPanels%vertices(3,nPanels+4) = anEdges%verts(2,childEdges(1))
	else
		!
		! divide parent edge 3
		!
		anEdges%hasChildren(edgeIndices(3)) = .TRUE.
		anEdges%children(:,edgeIndices(3)) = [nEdges+1,nEdges+2]
		if ( edgeOrientation(3)) then
			aPanels%edges(3,nPanels+3) = nEdges+1
			aPanels%edges(3,nPanels+4) = nEdges+2

			anEdges%verts(:,nEdges+1) = [vertIndices(3),nParticles+1]
			anEdges%leftPanel(nEdges+1) = nPanels+3
			anEdges%rightPanel(nEdges+1) = anEdges%rightPanel(edgeIndices(3))

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertIndices(4)]
			anEdges%leftPanel(nEdges+2) = nPanels+4
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeindices(3))
		else
			aPanels%edges(3,nPanels+3) = nEdges+2
			aPanels%edges(3,nPanels+4) = nEdges+1

			anEdges%verts(:,nEdges+1) = [vertIndices(4),nParticles+1]
			anEdges%leftPanel(nEdges+1) = anEdges%leftPanel(edgeIndices(3))
			anEdges%rightPanel(nEdges+1) = nPanels+4

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertIndices(3)]
			anEdges%leftPanel(nEdges+2) = anEdges%leftPanel(edgeIndices(3))
			anEdges%rightPanel(nEdges+2) = nPanels+3
		endif
		aParticles%x(:,nParticles+1) = Midpoint(aParticles%x(:,vertIndices(3)),aparticles%x(:,vertIndices(4)))
		aParticles%x0(:,nParticles+1) = Midpoint(aParticles%x0(:,vertIndices(3)),aparticles%x0(:,vertIndices(4)))

		aPanels%vertices(4,nPanels+3) = nParticles+1
		aPanels%vertices(3,nPanels+4) = nParticles+1

		nParticles = nParticles+1
		nEdges = nEdges+2
	endif
	if ( alreadyDivided(4) ) then
		!
		! connect to existing child edges
		!
		childEdges = anEdges%children(:,edgeIndices(4))
		if ( edgeorientation(4) ) then
			aPanels%edges(4,nPanels+4) = childEdges(1)
			aPanels%edges(4,nPanels+1) = childEdges(2)
			anEdges%leftPanel(childEdges(1)) = nPanels+4
			anEdges%leftPanel(childEdges(2)) = nPanels+1
		else
			aPanels%edges(4,nPanels+4) = childEdges(2)
			aPanels%edges(4,nPanels+1) = childEdges(1)
			anEdges%rightPanel(childEdges(1)) = nPanels+1
			anEdges%rightPanel(childEdges(2)) = nPanels+4
		endif
		aPanels%vertices(1,nPanels+4) = anEdges%verts(2,childEdges(1))
		aPanels%vertices(4,nPanels+1) = anEdges%verts(2,childEdges(1))
	else
		!
		! divide parent edge 4
		!
		anEdges%hasChildren(edgeIndices(4)) = .TRUE.
		anEdges%children(:,edgeIndices(4)) = [nEdges+1, nEdges+2]
		if ( edgeOrientation(4)) then
			aPanels%edges(4,nPanels+4) = nEdges+1
			aPanels%edges(4,nPanels+1) = nEdges+2

			anEdges%verts(:,nEdges+1) = [vertIndices(4),nParticles+1]
			anEdges%leftPanel(nEdges+1) = nPanels+4
			anEdges%rightPanel(nEdges+1) = anEdges%rightPanel(edgeIndices(4))

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertIndices(1)]
			anEdges%leftPanel(nEdges+2) = nPanels+1
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeIndices(4))
		else
			aPanels%edges(4,nPanels+4) = nEdges+2
			aPanels%edges(4,nPanels+1) = nEdges+1

			anEdges%verts(:,nEdges+1) = [vertIndices(1),nParticles+1]
			anEdges%leftPanel(nEdges+1) = anEdges%leftPanel(edgeIndices(4))
			anEdges%rightPanel(nEdges+1) = nPanels+1

			anEdges%verts(:,nEdges+2) = [nParticles+1,vertIndices(4)]
			anEdges%leftPanel(nEdges+2) = anEdges%leftPanel(edgeIndices(4))
			anEdges%rightPanel(nEdges+2) = nPanels+4
		endif
		aParticles%x(:,nParticles+1) = Midpoint(aParticles%x(:,vertIndices(4)), aParticles%x(:,vertIndices(1)))
		aParticles%x0(:,nParticles+1) = Midpoint(aParticles%x0(:,vertIndices(4)), aParticles%x0(:,vertIndices(1)))

		aPanels%vertices(1,nPanels+4) = nParticles+1
		aPanels%vertices(4,nPanels+1) = nParticles+1

		nParticles = nParticles + 1
		nEdges = nEdges+2
	endif
	!
	! change parent panel's active particle to a passive particle
	!
	aParticles%x(:,nParticles+1) = aPanels%x(:,panelIndex)
	aParticles%x0(:,nParticles+1) = aPanels%x0(:,panelIndex)

	aPanels%vertices(3,nPanels+1) = nParticles+1
	aPanels%vertices(4,nPanels+2) = nParticles+1
	aPanels%vertices(1,nPanels+3) = nParticles+1
	aPanels%vertices(2,nPanels+4) = nParticles+1

	anEdges%verts(:,nEdges+1) = [aPanels%vertices(2,nPanels+1),nParticles+1]
	anEdges%leftPanel(nEdges+1) = nPanels+1
	anEdges%rightPanel(nEdges+1) = nPanels+2

	anEdges%verts(:,nEdges+2) = [nParticles+1,aPanels%vertices(3,nPanels+4)]
	anEdges%leftPanel(nEdges+2) = nPanels+4
	anEdges%rightPanel(nEdges+2) = nPanels+3

	anEdges%verts(:,nEdges+3) = [aPanels%vertices(3,nPanels+2),nParticles+1]
	anEdges%leftPanel(nEdges+3) = nPanels+2
	anEdges%rightPanel(nEdges+3) = nPanels+3

	anEdges%verts(:,nEdges+4) = [nParticles+1, aPanels%vertices(1,nPanels+4)]
	anEdges%leftPanel(nEdges+4) = nPanels+1
	anEdges%rightPanel(nEdges+4) = nPanels+4

	aPanels%edges(2,nPanels+1) = nEdges+1
	aPanels%edges(3,nPanels+1) = nEdges+4
	aPanels%edges(4,nPanels+2) = nEdges+1
	aPanels%edges(3,nPanels+2) = nEdges+3
	aPanels%edges(1,nPanels+3) = nEdges+3
	aPanels%edges(4,nPanels+3) = nEdges+2
	aPanels%edges(1,nPanels+4) = nEdges+4
	aPanels%edges(2,nPanels+4) = nEdges+2

	nParticles = nParticles+1
	nEdges = nEdges+4
! DEBUG
!	do j=1,4
!		vertIndices = aPanels%vertices(:,nPanels+j)
!		edgeIndices = apanels%edges(:,nPanels+j)
!		write(logstring,'(A,I2)') 'child panel ', j
!		call LogMessage(log,DEBUG_LOGGING_LEVEL,'  ',logstring)
!		write(logstring,'(A,4I8)') 'vertices = ', vertIndices
!		call LogMessage(log, DEBUG_LOGGING_LEVEL, '  ', logstring)
!		write(logstring,'(A,4I8)') 'edges    = ', edgeIndices
!		call LogMessage(log, DEBUG_LOGGING_LEVEL, '  ', logstring)
!		do k=1,4
!			write(logString,'(A,I1,A,2I8)') 'edge ',k, ' verts = ', anedges%verts(:,edgeIndices(k))
!			call LogMessage(log,DEBUG_LOGGING_LEVEL,'  ',logstring)
!		enddo
!	enddo


	!
	! subpanel centers, nest, and area
	!
	do j=1,4
		aPanels%x(:,nPanels+j) = QuadCentroid( aParticles%x(:,aPanels%vertices(1,nPanels+j)), &
											   aParticles%x(:,aPanels%vertices(2,nPanels+j)), &
											   aParticles%x(:,aPanels%vertices(3,nPanels+j)), &
											   aParticles%x(:,aPanels%vertices(4,nPanels+j)) )
		aPanels%x0(:,nPanels+j) = QuadCentroid( aParticles%x0(:,aPanels%vertices(1,nPanels+j)), &
											    aParticles%x0(:,aPanels%vertices(2,nPanels+j)), &
											    aParticles%x0(:,aPanels%vertices(3,nPanels+j)), &
											    aParticles%x0(:,aPanels%vertices(4,nPanels+j)) )
	enddo
	!DEBUG
	!call LogMessage(log,DEBUG_LOGGING_LEVEL,' ','panel center set')
	do j=1,4
		aPanels%area(nPanels+j) = QuadPanelArea(self,nPanels+j)
		aPanels%nest(nPanels+j) = nestLevel+1
	enddo
	! DEBUG
	!call EndSection(log)

	!
	! update parent panel
	!
	aPanels%area(panelIndex) = 0.0_kreal
	aPanels%hasChildren(panelIndex) = .TRUE.
	aPanels%children(:,panelIndex) = [1,2,3,4] + nPanels


	aPanels%N = nPanels+4
	aPanels%N_Active = aPanels%N_Active + 3

	anEdges%N = nEdges

	aParticles%N = nParticles
end subroutine

function QuadPanelArea(self, panelIndex)
	real(kreal) :: QuadPanelArea
	type(PlaneMesh), intent(in) :: self
	integer(kint), intent(in) :: panelIndex
	!
	type(Particles), pointer :: aParticles
	type(Edges), pointer :: anEdges
	type(Panels), pointer :: aPanels
	integer(kint) :: vertList(8), edgeList(8), j, nVerts
	real(kreal) :: centerX(2)

	aParticles => self%particles
	anEdges => self%edges
	aPanels => self%panels

! 	! Error checking
!	if ( panelIndex > aPanels%N) then
!		call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'QuadPanelArea ERROR : out of bounds.')
!		return
!	endif

!	call StartSection(log,'entering QuadPanelArea')

	if ( aPanels%hasChildren(panelIndex)) then
		write(logString,'(A,I4,A)') 'panel ',panelIndex,' has children; area = 0.'
		call LogMessage(log,WARNING_LOGGING_LEVEL,'AREA WARNING : ',trim(logString))
		QuadPanelArea = 0.0_kreal
	else
		QuadPanelArea = 0.0_kreal
		centerX = aPanels%x(:,panelIndex)

		call CCWEdgesAndParticlesAroundPanel(edgeList,vertList,nVerts,self,panelIndex)

		do j=1,nVerts-1
			QuadPanelArea = QuadPanelArea + TriArea(aParticles%x(:,vertList(j)), centerX, aParticles%x(:,vertList(j+1)) )
		enddo
		QuadPanelArea = QuadPanelArea + TriArea(aParticles%x(:,vertList(nVerts)), centerX, aParticles%x(:,vertList(1)) )
	endif

!	call EndSection(log)
end function

function NearestRootPanel(self, xy)
! used to initialize tree search
	integer(kint) :: NearestRootPanel
	type(PlaneMesh), intent(in) :: self
	real(kreal), intent(in) :: xy(2)
	!
	integer(kint) :: j
	type(Panels), pointer :: aPanels
	real(kreal) :: currentMin, testDist

	NearestRootPanel = 0
	aPanels => self%panels
	currentMin = 1.0d15
	do j=1,self%nRootPanels
		testDist = Distance(aPanels%x(:,j), xy)
		if ( testDist < currentMin ) then
			currentMin = testDist
			NearestRootPanel = j
		endif
	enddo
end function

recursive subroutine LocatePointTreeSearch(inPanel, self, xy, rootPanel)
	integer(kint), intent(out) :: inPanel
	type(PlaneMesh), intent(in) :: self
	real(kreal), intent(in) :: xy(2)
	integer(kint), intent(in) :: rootPanel
	!
	integer(kint) :: j, currentPanel
	real(kreal) :: currentMin, testDist, centroid(2)
	type(Panels), pointer :: apanels
	type(Particles), pointer :: aparticles

	apanels => self%panels
	aparticles => self%particles

	if ( apanels%hasChildren(rootPanel) ) then
		currentMin = 1.0d15
		do j=1,4
			centroid = QuadCentroid( aParticles%x(:,aPanels%vertices(1,aPanels%children(j,rootPanel))), &
									 aParticles%x(:,aPanels%vertices(2,aPanels%children(j,rootPanel))), &
									 aParticles%x(:,aPanels%vertices(3,aPanels%children(j,rootPanel))), &
									 aParticles%x(:,aPanels%vertices(4,aPanels%children(j,rootPanel))) )
			testDist = Distance(xy, centroid)
			if ( testDist < currentMin ) then
				currentMin = testDist
				currentPanel = aPanels%children(j,rootPanel)
			endif
		enddo
		call LocatePointTreeSearch(inPanel, self, xy, currentPanel)
	else
		inPanel = rootPanel
	endif
end subroutine

recursive subroutine LocatePointWalkSearch(inPanel, self, xy, startPanel)
	integer(kint), intent(out) :: inPanel
	type(PlaneMesh), intent(in) :: self
	real(kreal), intent(in) :: xy(2)
	integer(kint), intent(in) :: startPanel
	!
	real(kreal) :: currentMin, testDist, centroid(2)
	integer(kint) :: adjPanels(8), nAdj, currentPanel, j
	type(Particles), pointer :: aparticles
	type(Panels), pointer :: apanels

	! Error checking
	if ( startPanel < 0 .OR. startPanel > self%panels%N) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logkey//'WalkSearch ERROR : ','invalid initial guess.')
		return
	endif

	apanels => self%panels
	aparticles => self%particles

	if (aPanels%hasChildren(startPanel) ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//'WalkSearch ERROR : ','divided panel.')
		return
	endif

	centroid = QuadCentroid( aParticles%x(:, aPanels%vertices(1,startPanel) ), &
						  	 aParticles%x(:, aPanels%vertices(2,startPanel) ), &
						  	 aParticles%x(:, aPanels%vertices(3,startPanel) ), &
						  	 aParticles%x(:, aPanels%vertices(4,startPanel) ) )
	currentMin = Distance(centroid, xy)
	currentPanel = startPanel

	call FindAdjacentPanels(adjPanels, nAdj, self, startPanel)

	do j=1, nAdj
		centroid = QuadCentroid( aParticles%x(:, aPanels%vertices(1,adjPanels(j)) ), &
								 aParticles%x(:, aPanels%vertices(2,adjPanels(j)) ), &
								 aParticles%x(:, aPanels%vertices(3,adjPanels(j)) ), &
								 aParticles%x(:, aPanels%vertices(4,adjPanels(j)) ) )
		testDist = Distance(centroid, xy)
		if ( testDist < currentMin ) then
			currentMin = testDist
			currentPanel = adjPanels(j)
		endif
	enddo

	if ( currentPanel == startPanel) then
		inPanel = currentPanel
		return
	else
		call LocatePointWalkSearch(inPanel, self, xy, currentPanel)
	endif
end subroutine

function PlaneParticleMax(maxnest)
	integer(kint) :: PlaneParticleMax
	integer(kint), intent(in) :: maxnest
	!
	integer(kint) :: j

	PlaneParticleMax = 3
	do j=1,maxnest
		PlaneParticleMax = PlaneParticleMax + 2**j
	enddo
	PlaneParticleMax = PlaneParticleMax*PlaneParticleMax
end function

function PlanePanelMax(maxNest)
	integer(kint) :: PlanePanelMax
	integer(kint), intent(in) :: maxNest
	!
	integer(kint) :: j
	PlanePanelMax = 4
	do j=1,maxNest
		PlanePanelMax = PlanePanelMax + 4**(j+1)
	enddo
end function

function PlaneEdgeMax(maxNest)
	integer(kint) :: PlaneEdgeMax
	integer(kint), intent(in) ::maxNest
	!
	PlaneEdgeMax = 4**(maxNest+2) - 4**maxnest
	! TO DO : this function is inexact.  It should match the following sequence:
	!	nest   : 0		1		2		3		4		5		...
	!	nEdges : 12		52		196		740		2852	11172	...
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

