module EdgesModule

use NumberKindsModule
use LoggerModule
use STDIntVectorModule
use ParticlesModule
use SphereGeomModule, only : crossProduct, SphereMidpoint, SphereTriArea, SphereDistance, ChordDistance
use PlaneGeomModule, only : Midpoint, TriArea
use ParticlesModule

implicit none

private
public Edges
public New, Delete, Copy
public InsertEdge, DivideEdge
public RecordIncidentEdgeAtParticles
public positiveEdge
public onBoundary
public ReplaceIncidentEdgeWithChild
public GetLeafEdgesFromParent, AreaFromLeafEdges
public EdgeLength

type Edges
	integer(kint), pointer :: orig(:) => null()
	integer(kint), pointer :: dest(:) => null()
	integer(kint), pointer :: rightFace(:) => null()
	integer(kint), pointer :: leftFace(:) => null()
	integer(kint), pointer :: child1(:) => null()
	integer(kint), pointer :: child2(:) => null()
	logical(klog), pointer :: hasChildren(:) => null()
	integer(kint), pointer :: parent(:) => null()
	integer(kint) :: N = 0
	integer(kint) :: N_Max = 0
end type

interface New
	module procedure newPrivate
end interface

interface Delete
	module procedure deletePrivate
end interface

interface Copy
	module procedure copyPrivate
end interface

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'Edges'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains

subroutine NewPrivate(self, nMax )
	type(Edges), intent(out) :: self
	integer(kint), intent(in) :: nMax
	
	if ( .NOT. logInit ) call InitLogger(log, procRank )
	
	if ( nMax <= 0 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, "New Edges ERROR : invalid nMax.")
		return
	endif
	
	self%N_Max = nMax
	self%N = 0
	
	allocate(self%orig(nMax))
	self%orig = 0
	allocate(self%dest(nMax))
	self%dest = 0
	allocate(self%leftFace(nMax))
	self%leftFace = 0
	allocate(self%rightFace(nMax))
	self%rightFace = 0
	allocate(self%hasChildren(nMax))
	self%hasChildren = .FALSE.
	allocate(self%child1(nMax))
	self%child1 = 0
	allocate(self%child2(nMax))
	self%child2 = 0
	allocate(self%parent(nMax))
	self%parent = 0
end subroutine

subroutine DeletePrivate(self)
	type(Edges), intent(inout) :: self
	deallocate(self%orig)
	deallocate(self%dest)
	deallocate(self%leftFace)
	deallocate(self%rightFace)
	deallocate(self%hasChildren)
	deallocate(self%child1)
	deallocate(self%child2)
	deallocate(self%parent)
	self%N = 0
	self%N_Max = 0
end subroutine

subroutine copyPrivate( self, other )
	type(Edges), intent(inout) :: self
	type(Edges), intent(in) :: other
	!
	integer(kint) :: j
	
	if ( self%N_Max < other%N ) then
		call LogMessage( log, ERROR_LOGGING_LEVEL, logkey, "CopyEdges ERROR : not enough memory.")
		return
	endif
	
	do j = 1, other%N
		self%orig(j) = other%orig(j)
		self%dest(j) = other%dest(j)
		self%leftFace(j) = other%leftFace(j)
		self%rightFace(j) = other%rightFace(j)
		self%hasChildren(j) = other%hasChildren(j)
		self%child1(j) = other%child1(j)
		self%child2(j) = other%child2(j)
		self%parent(j) = other%parent(j)
	enddo
	self%N = other%N
end subroutine

subroutine InsertEdge( self, aParticles, origIndex, destIndex, leftFace, rightFace )
	type(Edges), intent(inout) :: self
	type(Particles), intent(inout) :: aParticles
	integer(kint), intent(in) :: origIndex, destIndex, leftFace, rightFace
	!
	integer(kint) :: n
	
	if ( self%N >= self%N_Max ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logKey, " InsertEdge : out of memory. ")
		return
	endif
	
	n = self%N
	
	self%orig( n + 1 ) = origIndex
	self%dest( n + 1 ) = destIndex
	self%leftFace( n + 1) = leftFace
	self%rightFace(n + 1) = rightFace
	
	call RecordIncidentEdgeAtParticles( self, n + 1, 	aParticles )
	
	self%N = n + 1
end subroutine

function EdgeLength(self, edgeIndex, aParticles )
	real(kreal) :: EdgeLength
	integer(kint), intent(in) :: edgeIndex
	type(Edges), intent(in) :: self
	type(Particles), intent(in) :: aParticles
	!
	real(kreal) :: v0(3), v1(3)
	v0 = PhysCoord(aParticles, self%orig(edgeIndex))
	v1 = PhysCoord(aParticles, self%dest(edgeIndex))
	EdgeLength = 0.0_kreal
	if ( aParticles%geomKind == SPHERE_GEOM ) then
		EdgeLength = SphereDistance( v0, v1 )
	else
		EdgeLength = ChordDistance(v0, v1)
	endif
end function

subroutine RecordIncidentEdgeAtParticles( self, edgeIndex, aParticles )
	type(Edges), intent(in) :: self
	integer(kint), intent(in) :: edgeIndex
	type(Particles), intent(inout) :: aParticles
	!
	logical :: duplicateEdge
	integer(kint) :: j, origParticle, destParticle
	real(kreal) :: angleVal
	
	if ( self%hasChildren(edgeIndex) ) then
		return
	else
		origParticle = self%orig(edgeIndex)
		destParticle = self%dest(edgeIndex)
		!
		! origin vertex
		!
		duplicateEdge = .FALSE.
		do j = 1, aParticles%nEdges( origParticle )
			if ( aParticles%incidentEdges( j, origParticle ) == edgeIndex ) duplicateEdge = .TRUE.
		enddo
		if ( .NOT. duplicateEdge ) then
			if ( aParticles%nEdges( origParticle ) >= MAX_VERTEX_DEGREE ) then
				call LogMessage(log,ERROR_LOGGING_LEVEL,logkey, " recordEdgeAtOrigin : out of memory.")
				return
			endif
			angleVal = edgeAngleAtOrig( self,  edgeIndex, aParticles)
			aParticles%incidentEdges( aParticles%nEdges( origParticle ) + 1, origParticle ) = edgeIndex
			aParticles%incidentAngles(aParticles%nEdges( origParticle ) + 1, origParticle ) = angleVal
			aParticles%nEdges(origParticle) = aParticles%nEdges(origParticle) + 1
		endif
		
		!
		! destination vertex
		!
		duplicateEdge = .FALSE.
		do j = 1, aParticles%nEdges( self%dest( edgeIndex ))
			if ( aParticles%incidentEdges( j, self%dest(edgeIndex)) == edgeIndex ) duplicateEdge = .TRUE.
		enddo
		if ( .NOT. duplicateEdge ) then
			if ( aParticles%nEdges( destParticle) >= MAX_VERTEX_DEGREE ) then
				call LogMessage(log,ERROR_LOGGING_LEVEL,logkey, " recordEdgeAtDestination : out of memory.")
				return
			endif
			angleVal = edgeAngleAtDest( self, edgeIndex, aParticles )
			aParticles%incidentEdges( aParticles%nEdges(destParticle) + 1, destParticle ) = edgeIndex
			aParticles%incidentAngles(aParticles%nEdges(destParticle) + 1, destParticle ) = angleVal
			aParticles%nEdges(destParticle) = aParticles%nEdges(destParticle) + 1
		endif
	endif
end subroutine

function edgeAngleAtOrig( self, edgeIndex, aParticles )
	real(kreal) :: edgeAngleAtOrig
	type(Edges), intent(in) :: self
	integer(kint), intent(in) :: edgeIndex
	type(Particles), intent(in) :: aParticles
	!
	real(kreal) :: v0(3), v1(3), edge1Vec(3), newEdgeVec(3), cp(3), dp
	
	edgeAngleAtOrig = 0.0_kreal
	if ( aParticles%geomKind == PLANAR_GEOM ) then
		!
		!	define angle relative to positive real axis
		!
		edgeAngleAtOrig = atan2( aParticles%y( self%dest(edgeIndex)) - aParticles%y( self%orig(edgeIndex)), &
								 aParticles%x( self%dest(edgeIndex)) - aParticles%x( self%orig(edgeIndex)) )
	elseif ( aParticles%geomKind == SPHERE_GEOM ) then
		!
		!	define angle relative to first edge at particle
		!
		if ( aParticles%nEdges( self%orig(edgeIndex) ) > 1 ) then
			if ( self%orig(aParticles%incidentEdges(1,self%orig(edgeIndex))) == self%orig(edgeIndex) ) then
				edge1Vec = edgeVector( self, aParticles%incidentEdges(1,self%orig(edgeIndex)), aParticles )
			elseif ( self%dest(aParticles%incidentEdges(1,self%orig(edgeIndex))) == self%orig(edgeIndex) ) then
				edge1Vec = -edgeVector(self, aParticles%incidentEdges(1,self%orig(edgeIndex)), aParticles)
			else
				call LogMessage(log, ERROR_LOGGING_LEVEL, "edgeAngleAtOrig : ", "connectivity error.")
				return
			endif
			newEdgeVec = edgeVector( self, edgeIndex, aParticles)
			edge1Vec = edge1Vec/sqrt(sum( edge1Vec*edge1Vec))
			newEdgeVec = newEdgeVec/sqrt(sum( newEdgeVec*newEdgeVec))
			
			cp = crossProduct(edge1Vec, newEdgeVec)
			edgeAngleAtOrig = atan2( sqrt(sum(cp*cp)), sum(edge1Vec*newEdgeVec))
		endif
	else
		call LogMessage(log, WARNING_LOGGING_LEVEL, logkey, " geomKind not implemented yet.")
	endif
end function

function edgeAngleAtDest( self, edgeIndex, aParticles )
	real(kreal) :: edgeAngleAtDest
	type(Edges), intent(in) :: self
	integer(kint), intent(in) :: edgeIndex
	type(Particles), intent(in) :: aParticles
	!
	real(kreal) :: edge1Vec(3), newEdgeVec(3), cp(3), dp
	
	edgeAngleAtDest = 0.0_kreal
	if ( aParticles%geomKind == PLANAR_GEOM ) then
		!
		!	define angle relative to positive real axis
		!
		edgeAngleAtDest = atan2( aParticles%y( self%orig(edgeIndex)) - aParticles%y( self%dest(edgeIndex)), &
								 aParticles%x( self%orig(edgeIndex)) - aParticles%x( self%dest(edgeIndex)) )
	elseif ( aParticles%geomKind == SPHERE_GEOM ) then
		!
		!	define angle relative to first edge at particle
		!
		if ( aParticles%nEdges( self%dest(edgeIndex) ) > 1 ) then
			if ( self%orig(aParticles%incidentEdges(1,self%dest(edgeIndex))) == self%dest(edgeIndex) ) then
				edge1Vec = edgeVector(self, aParticles%incidentEdges(1,self%dest(edgeIndex)), aParticles)
			elseif ( self%dest(aParticles%incidentEdges(1,self%dest(edgeIndex))) == self%dest(edgeIndex)) then
				edge1Vec = -edgeVector(self, aParticles%incidentEdges(1,self%dest(edgeIndex)), aParticles)
			else
				call LogMessage(log, ERROR_LOGGING_LEVEL, "edgeAngleAtDest : ", "connectivity error.")
				return
			endif
		endif
		newEdgeVec = -edgeVector(self,edgeIndex,aParticles)
		
		edge1Vec = edge1Vec/sqrt(sum( edge1Vec*edge1Vec))
		newEdgeVec = newEdgeVec/sqrt(sum( newEdgeVec*newEdgeVec))
		
		cp = crossProduct(edge1Vec, newEdgeVec)
		edgeAngleAtDest = atan2( sqrt(sum(cp*cp)), sum(edge1Vec*newEdgeVec))
	else
		call LogMessage(log, WARNING_LOGGING_LEVEL, logkey, " geomKind not implemented yet.")
	endif
end function 

function edgeVector( self,  edgeIndex, aParticles )
	real(kreal), dimension(3) :: edgeVector
	type(Edges), intent(in) :: self
	integer(kint), intent(in) :: edgeIndex
	type(Particles), intent(in) :: aParticles
	edgeVector(1) = aParticles%x( self%dest(edgeIndex)) - aParticles%x( self%orig(edgeIndex))
	edgeVector(2) = aParticles%y( self%dest(edgeIndex)) - aParticles%y( self%orig(edgeIndex))
	edgeVector(3) = aParticles%z( self%dest(edgeIndex)) - aParticles%z( self%orig(edgeIndex))
end function

subroutine DivideEdge( self, edgeIndex, aParticles )
	type(Edges), intent(inout) :: self
	integer(kint), intent(in) :: edgeIndex
	type(Particles), intent(inout) :: aParticles
	!
	real(kreal) :: midPt(3), lagMidPt(3), v0(3), v1(3), lV0(3), lV1(3)
	integer(kint) :: pInsertIndex
	
	if ( self%N + 2 > self%N_Max ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, " DivideEdge : out of memory.")
		return
	endif
	
	v0 = PhysCoord(aParticles, self%orig(edgeIndex))
	lV0 = LagCoord(aParticles, self%orig(edgeIndex))
	v1 = PhysCoord(aParticles, self%dest(edgeIndex))
	lV1 = LagCoord(aParticles, self%dest(edgeIndex))
		
	if ( aParticles%geomKind == PLANAR_GEOM ) then
		midPt(1:2) = Midpoint( v0(1:2), v1(1:2) )
		lagMidPt(1:2) = Midpoint( lv0(1:2), lv1(1:2))
	elseif ( aParticles%geomKind == SPHERE_GEOM ) then
		midPt = SphereMidpoint( v0, v1 )
		lagMidPt = SphereMidpoint( lv0, lv1 )
	else
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, "DivideEdge : geomKind not implemented.")
		return
	endif
	
	pInsertIndex = aParticles%N+1
	call InsertParticle( aParticles, midPt, lagMidPt )
	
	self%hasChildren(edgeIndex) = .TRUE.
	self%child1(edgeIndex) = self%N + 1
	self%child2(edgeIndex) = self%N + 2
	self%parent(self%N+1) = edgeIndex
	self%parent(self%N+2) = edgeIndex
	
	self%orig( self%N + 1 ) = self%orig(edgeIndex)
	self%dest( self%N + 1 ) = pInsertIndex
	self%leftFace( self%N + 1 ) = self%leftFace( edgeIndex ) 
	self%rightFace(self%N + 1 ) = self%rightFace( edgeIndex)
	
	self%orig( self%N + 2 ) = pInsertIndex
	self%dest( self%N + 2 ) = self%dest(edgeIndex)
	self%leftFace( self%N + 2 ) = self%leftFace( edgeIndex ) 
	self%rightFace(self%N + 1 ) = self%rightFace( edgeIndex)
	
	call replaceIncidentEdgeWithChild( self, edgeIndex, aParticles)
	call RecordIncidentEdgeAtParticles(self, self%N + 1, aParticles)
	call RecordIncidentEdgeAtParticles(self, self%N + 2, aParticles)
	
	self%N = self%N + 2
end subroutine

subroutine replaceIncidentEdgeWithChild( self, parentIndex, aParticles )
	type(Edges), intent(in) :: self
	integer(kint), intent(in) :: parentIndex
	type(Particles), intent(inout) :: aParticles
	!
	integer(kint) :: j, pEdgeIndex
	
	!
	! replace parent edge at origin vertex
	!
	pEdgeIndex = 0
	do j = 1, aParticles%nEdges( self%orig(parentIndex) )
		if ( aParticles%incidentEdges( j, self%orig(parentIndex) ) == parentIndex ) then
			pEdgeIndex = j
			exit
		endif
	enddo
	
	if ( pEdgeIndex == 0 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, "replaceIncidentEdgeWithChild ERROR : Parent not found.")
		return
	endif
	
	aParticles%incidentEdges( pEdgeIndex, self%orig(parentIndex)) = self%child1(parentIndex)
	
	!
	! replace parent edge at destination vertex
	!
	pEdgeIndex = 0
	do j = 1, aParticles%nEdges( self%dest(parentIndex) )
		if ( aParticles%incidentEdges(j, self%dest(parentIndex)) == parentIndex ) then
			pEdgeIndex = j
			exit
		endif
	enddo
	
	if ( pEdgeIndex == 0 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey, "replaceIncidentEdgeWithChild ERROR : Parent not found.")
		return
	endif
	
	aParticles%incidentEdges( pEdgeIndex, self%dest(parentIndex)) = self%child2(parentIndex)
end subroutine

function positiveEdge( anEdges, faceIndex, edgeIndex )
	logical(klog) :: positiveEdge
	type(Edges), intent(in) :: anEdges
	integer(kint), intent(in) :: faceIndex
	integer(kint), intent(in) :: edgeIndex
	
	positiveEdge = ( faceIndex == anEdges%leftFace(edgeIndex))
end function

function onBoundary( anEdges, edgeIndex )
	logical(klog) :: onBoundary
	type(Edges), intent(in) :: anEdges
	integer(kint), intent(in) :: edgeIndex
	onBoundary = ( anEdges%leftFace(edgeIndex) == 0 .OR. anEdges%rightFace(edgeIndex) == 0 )
end function

subroutine GetLeafEdgesFromParent( anEdges, parentIndex, leafEdges )
	type(Edges), intent(in) :: anEdges
	integer(kint), intent(in) :: parentIndex
	type(STDIntVector), intent(out) :: leafEdges
	!
	integer(kint) :: i, nLeaves
	logical(klog) :: keepGoing
	
	call initialize(leafEdges)
	call leafEdges%pushBack(parentIndex)
	nLeaves = 1
	
	keepGoing = .FALSE.
	if ( anEdges%hasChildren(parentIndex) )	keepGoing = .TRUE.
	
	do while (keepGoing)
		do i = 1, nLeaves
			if ( anEdges%hasChildren( leafEdges%int(i) ) ) then
				call leafEdges%replace(i, anEdges%child1(i) )
				call leafEdges%insert(i+1, anEdges%child2(i))
			endif
		enddo
		nLeaves = leafEdges%N
		keepGoing = .FALSE.
		do i = 1, nLeaves
			if ( anEdges%hasChildren(leafEdges%int(i)) ) keepGoing = .TRUE.
		enddo
	enddo
end subroutine

function AreaFromLeafEdges( self, aParticles, centerParticle, leafEdges, nLeaves )
	real(kreal) :: AreaFromLeafEdges
	type(Edges), intent(in) :: self
	type(Particles), intent(in) :: aParticles
	integer(kint), intent(in) :: centerParticle
	integer(kint), intent(in) :: leafEdges(:)
	integer(kint), intent(in) :: nLeaves
	!
	integer(kint) :: i
	real(kreal) :: centerVec(3), v1Vec(3), v2Vec(3)
	
	AreaFromLeafEdges = 0.0_kreal
	centerVec = PhysCoord( aParticles, centerParticle )
	if ( aParticles%geomKind == PLANAR_GEOM ) then
		do i = 1, nLeaves
			v1Vec = PhysCoord( aParticles, self%orig(leafEdges(i)) )
			v2Vec = PhysCoord( aParticles, self%dest(leafEdges(i)) )
			AreaFromLeafEdges = AreaFromLeafEdges + TriArea( v1Vec(1:2), centerVec(1:2), v2Vec(1:2) )
		enddo
	elseif ( aParticles%geomKind == SPHERE_GEOM ) then
		do i = 1, nLeaves
			v1Vec = PhysCoord( aParticles, self%orig(leafEdges(i)) )
			v2Vec = PhysCoord( aParticles, self%dest(leafEdges(i)) )
			AreaFromLeafEdges = AreaFromLeafEdges + SphereTriArea( v1Vec, centerVec, v2Vec )
		enddo
	else
		call LogMessage(log,ERROR_LOGGING_LEVEL,logkey//" AreaFromLeafEdges ERROR : ", "geomKind not implemented.")
		return
	endif
end function

subroutine InitLogger(aLog,rank)
! Initialize a logger for this module and processor
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I3,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

end module