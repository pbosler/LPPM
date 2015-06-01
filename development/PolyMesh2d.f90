module PolyMesh2dModule

use NumberKindsModule
use STDIntVectorModule
use LoggerModule
use ParticlesModule
use EdgesModule
use FacesModule
use PlaneGeomModule
use SphereGeomModule

implicit none

private
public PolyMesh2d
public New, Delete, Copy
public locateFaceContainingPoint
!public pointIsOutsideMesh
public CCWEdgesAroundFace, CCWVerticesAroundFace, CCWAdjacentFaces
public CCWFacesAroundVertex
!public ResetSurfaceArea
!public VertPhysCoord, VertLagCoord
public LogStats, PrintDebugInfo
public WriteMeshToMatlab

type PolyMesh2d
	type(Particles) :: particles
	type(Edges) :: edges
	type(Faces) :: faces
	integer(kint) :: faceKind = 0
	integer(kint) :: geomKind = 0
	integer(kint) :: meshSeed = 0
	integer(kint) :: initNest = -1
	integer(kint) :: amrLimit = 0
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

interface LogStats
	module procedure LogStatsPrivate
end interface

interface PrintDebugInfo
	module procedure PrintDebugPrivate
end interface

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'PolyMesh2d'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logstring

contains

subroutine newPrivate(self, meshSeed, initNest, maxNest, amrLimit, ampFactor )
	type(PolyMesh2d), intent(out) :: self
	integer(kint), intent(in) :: meshSeed
	integer(kint), intent(in) :: maxNest
	integer(kint), intent(in) :: initNest
	integer(kint), intent(in) :: amrLimit
	real(kreal), intent(in) :: ampFactor
	!
	integer(kint) :: nMaxParticles, nMaxVertices, nMaxFaces, nMaxEdges
	integer(kint) :: i, j, startIndex, nFacesOld

	if ( .NOT. logInit) call InitLogger(log, procRank)
	
	self%meshSeed = meshSeed
	self%initNest = initNest
	self%amrLimit = amrLimit
	if ( meshSeed < TRI_HEX_SEED .OR. meshSeed > CUBED_SPHERE_SEED ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL, logkey," NewMesh ERROR : invalid meshSeed.")
		return
	endif
	self%faceKind = FaceKindFromSeed(meshSeed)
	self%geomKind = GeomKindFromSeed(meshSeed)
	
	print *, "DEBUG : nVertices = ", nVerticesInMesh(self, maxNest)
	print *, "DEBUG : nFaces = ", nFacesInMesh(self, maxNest)
		
	nMaxParticles = nVerticesInMesh(self, maxNest) + nFacesInMesh(self, maxNest)
	nMaxVertices = 0
	nMaxFaces = 0
	nMaxEdges = 0
	do i = 0, initNest
		nMaxFaces = nMaxFaces + nFacesInMesh(self, i)
		nMaxEdges = nMaxEdges + nEdgesInMesh(self, nVerticesInMesh(self,i), nFacesInMesh(self,i) )
	enddo
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL,logkey//" allocating memory for nParticles = ", nMaxParticles )
	call LogMessage(log, DEBUG_LOGGING_LEVEL,logkey//" allocating memory for nEdges = ", nMaxEdges )
	call LogMessage(log, DEBUG_LOGGING_LEVEL,logkey//" allocating memory for nFaces = ", nMaxFaces )
	call LogMessage(log, DEBUG_LOGGING_LEVEL,logkey//" faceKind = ", self%faceKind )
	call LogMessage(log, DEBUG_LOGGING_LEVEL,logkey//" geomKind = ", self%geomKind )
	
	call New(self%particles, nMaxParticles, self%geomKind )
	call New(self%edges, nMaxEdges)
	call New(self%faces, self%faceKind, nMaxFaces )
	
	call initializeMeshFromSeed(self, ampFactor)
	
	startIndex = 1
	do i = 1, initNest
		nFacesOld = self%faces%N
		do j = startIndex, nFacesOld
			if ( self%faceKind == TRI_PANEL ) then
				call DivideTriFace(self%faces, j, self%particles, self%edges)
			elseif ( self%faceKind == QUAD_PANEL ) then
				call DivideQuadFace(self%faces,j, self%particles, self%edges)
			endif
		enddo
		startIndex = nFacesOld + 1
	enddo
	
	do i = 1, self%particles%N
		if ( self%particles%isPassive(i) ) call SortIncidentEdgesAtParticle( self%particles, i )
	enddo
end subroutine

subroutine deletePrivate(self)
	type(PolyMesh2d), intent(inout) :: self
	call Delete(self%faces)
	call Delete(self%edges)
	call Delete(self%particles)
	self%meshSeed = 0
	self%initNest = 0
	self%amrLimit = 0
	self%geomKind = 0
	self%faceKind = 0
end subroutine

subroutine copyPrivate( self, other )
	type(PolyMesh2d), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: other
	self%meshSeed = other%meshSeed
	self%initNest = other%initNest
	self%amrLimit = other%amrLimit
	self%geomKind = other%geomKind
	self%faceKind = other%faceKind
	call Copy(self%particles, other%particles)
	call Copy(self%edges, other%edges)
	call Copy(self%faces, other%faces)
end subroutine

function LocateFaceContainingPoint(self, queryPt)
	integer(kint) :: locateFaceContainingPoint
	type(PolyMesh2d), intent(in) :: self
	real(kreal), intent(in) :: queryPt(:)
	!
	integer(kint) :: treeStart
	integer(kint) :: walkStart
	treeStart = nearestRootFace(self, queryPt)
	walkStart = locatePointTreeSearch(self, queryPt, treeStart)
	locateFaceContainingPoint = locatePointWalkSearch(self, queryPt, walkStart)
end function

subroutine CCWEdgesAroundFace( self, leafEdges, faceIndex )
	type(PolyMesh2d), intent(in) :: self
	type(STDIntVector), intent(out) :: leafEdges
	integer(kint), intent(in) :: faceIndex
	!
	integer(kint) :: i, j, nParentEdges
	type(STDIntVector) :: edgeLeaves(4)
	
	if ( self%faceKind == TRI_PANEL ) then
		nParentEdges = 3
	elseif ( self%faceKind == QUAD_PANEL ) THEN
		nParentEdges = 4
	endif
	
	do i = 1, nParentEdges
		call GetLeafEdgesFromParent( self%edges, self%faces%edges(i,faceIndex), edgeLeaves(i))
	enddo
	
	call initialize(leafEdges)
	do i = 1, nParentEdges
		do j = 1, edgeLeaves(i)%N
			call leafEdges%pushBack(edgeLeaves(i)%integers(j))
		enddo
	enddo
end subroutine

subroutine CCWVerticesAroundFace( self, verts, faceIndex )
	type(PolyMesh2d), intent(in) :: self
	type(STDIntVector), intent(out) :: verts
	integer(kint), intent(in) :: faceIndex
	!
	type(STDIntVector) :: leafEdges
	integer(kint) :: i
	
	call CCWEdgesAroundFace(self, leafEdges, faceIndex)
	call initialize(verts)
	
	if ( positiveEdge(self%edges, faceIndex, leafEdges%int(1)) ) then
		call verts%pushBack(self%edges%orig(leafEdges%int(1)))
	else
		call verts%pushBack(self%edges%dest(leafEdges%int(1)))
	endif
	do i = 1, leafEdges%N
		if ( positiveEdge(self%edges, faceIndex, leafEdges%int(i)) ) then
			call verts%pushBack(self%edges%dest(leafEdges%int(i)))
		else
			call verts%pushBack(self%edges%orig(leafEdges%int(i)))
		endif
	enddo	
end subroutine

subroutine CCWAdjacentFaces( self, adjFaces, faceIndex )
	type(PolyMesh2d), intent(in) :: self
	type(STDIntVector), intent(out) :: adjFaces
	integer(kint), intent(in) :: faceIndex
	!
	type(STDIntVector) :: leafEdges
	integer(kint) :: i
	
	call CCWEdgesAroundFace(self, leafEdges, faceIndex )
	
	call initialize(adjFaces)
	
	do i = 1, leafEdges%N
		if ( positiveEdge(self%edges, faceIndex, leafEdges%int(i)) ) then
			call adjFaces%pushBack(self%edges%rightFace(leafEdges%int(i)))
		else
			call adjFaces%pushBack(self%edges%leftFace(leafEdges%int(i)))
		endif
	enddo
end subroutine

subroutine CCWFacesAroundVertex( self, adjFaces, vertexIndex )
	type(PolyMesh2d), intent(in) :: self
	type(STDIntVector), intent(out) :: adjFaces
	integer(kint), intent(in) :: vertexIndex
	!
	integer(kint) :: i
	
	if ( self%particles%area(vertexIndex) > 0.0_kreal ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" CCWFacesAroundVertex ERROR : ", " vertices should have 0 area.")
		return
	endif
	
	call initialize(adjFaces)
	do i = 1, self%particles%nEdges(vertexIndex)
		call adjFaces%pushBack( self%particles%incidentEdges(i, vertexIndex) )
	enddo
end subroutine

function nearestRootFace(self, queryPt)
	integer(kint) :: nearestRootFace
	type(PolyMesh2d), intent(in) :: self
	real(kreal), intent(in) :: queryPt(:)
	!
	integer(kint) :: i, nRootFaces
	real(kreal) :: dist, testDist, cntd(3)
	
	nearestRootFace = 0
	
	select case (self%meshSeed)
		case ( TRI_HEX_SEED, CUBED_SPHERE_SEED )
			nRootFaces = 6
		case ( QUAD_RECT_SEED, QUAD_RECT_PERIODIC_SEED, POLAR_DISC_SEED)
			nRootFaces = 4
		case ( ICOS_TRI_SPHERE_SEED)
			nRootFaces = 20
	end select
	
	if ( self%geomKind == PLANAR_GEOM ) then
		cntd = faceCentroid( self%faces, 1, self%particles)
		dist = Distance( cntd(1:2), queryPt(1:2) )
		do i = 2, nRootFaces
			cntd = faceCentroid( self%faces, i, self%particles )
			testDist = Distance(cntd(1:2), queryPt(1:2))
			if ( testDist < dist ) then
				dist = testDist
				nearestRootFace = i
			endif
		enddo
	elseif ( self%geomKind == SPHERE_GEOM ) then
		cntd = faceCentroid(self%faces, 1, self%particles)
		dist = SphereDistance( cntd, queryPt)
		do i = 2, nRootFaces
			cntd = faceCentroid(self%faces, i, self%particles)
			testDist = SphereDistance( cntd, queryPt )
			if ( testDist < dist ) then
				dist = testDist
				nearestRootFace = i
			endif
		enddo
	else
		call LogMessage(log, ERROR_LOGGING_LEVEL,logkey//" nearestRootFace ERROR: "," geomKind not implemented.")
	endif
end function

function FaceKindFromSeed( seed )
	integer(kint) :: FaceKindFromSeed
	integer(kint), intent(in) :: seed
	FaceKindFromSeed = 0
	select case (seed)
		case ( TRI_HEX_SEED, ICOS_TRI_SPHERE_SEED)
			FaceKindFromSeed = TRI_PANEL
		case ( QUAD_RECT_SEED, QUAD_RECT_PERIODIC_SEED, POLAR_DISC_SEED, CUBED_SPHERE_SEED)
			FaceKindFromSeed = QUAD_PANEL
	end select
end function

function GeomKindFromSeed( seed )
	integer(kint) :: GeomKindFromSeed
	integer(kint), intent(in) :: seed
	GeomKindFromSeed = 0
	select case (seed)
		case (TRI_HEX_SEED, QUAD_RECT_SEED, QUAD_RECT_PERIODIC_SEED, POLAR_DISC_SEED)
			GeomKindFromSeed = PLANAR_GEOM
		case (ICOS_TRI_SPHERE_SEED, CUBED_SPHERE_SEED)
			GeomKindFromSeed = SPHERE_GEOM
	end select
end function

recursive function locatePointTreeSearch( self, queryPt, index ) result(faceIndex)
	integer(kint) :: faceIndex
	type(PolyMesh2d), intent(in) :: self
	real(kreal), intent(in) :: queryPt(:)
	integer(kint), intent(in) :: index
	!
	real(kreal) :: dist, testDist, cntd(3)
	integer(kint) :: i, nextIndex
	
	faceIndex = 0
	
	dist = 1.0d20
	
	if ( self%faces%hasChildren(index) ) then
		if ( self%geomKind == PLANAR_GEOM ) then
			do i = 1, 4
				cntd = FaceCentroid(self%faces, self%faces%children(i,index), self%particles)
				testDist = Distance( cntd(1:2), queryPt(1:2) )
				if ( testDist < dist ) then
					dist = testDist
					nextIndex = self%faces%children(i,index)
				endif
			enddo
		elseif ( self%geomKind == SPHERE_GEOM ) then
			do i = 1, 4
				cntd = FaceCentroid(self%faces, self%faces%children(i,index), self%particles )
				testDist = SphereDistance( cntd, queryPt )
				if ( testDist < dist ) then
					dist = testDist
					nextIndex = self%faces%children(i,index)
				endif
			enddo
		else
			call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" locatePointTreeSearch ERROR : ", " geomKind not implemented.")
			return
		endif
		faceIndex = locatePointTreeSearch( self, queryPt, nextIndex )
	else
		faceIndex = index
	endif	
end function

recursive function locatePointWalkSearch(self, queryPt, index) result(faceIndex)
	integer(kint) :: faceIndex
	type(PolyMesh2d), intent(in) :: self
	real(kreal), intent(in) :: queryPt(:)
	integer(kint), intent(in) :: index
	!
	real(kreal) :: dist, testDist, cntd(3)
	
	faceIndex = 0
	if ( self%faces%hasChildren(index) ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" locatePointWalkSearch ERROR:"," expected a leaf face.")
		return
	endif
	
	cntd = FaceCentroid(self%faces, index, self%particles)
	if (self%geomKind == PLANAR_GEOM ) then
		
	elseif ( self%geomKind == SPHERE_GEOM ) then
	else
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" locatePointWalkSearch ERROR : ", " geomKind not implemented.")
		return
	endif
end function

function nVerticesInMesh( self, initNest )
	integer(kint) :: nVerticesInMesh
	type(PolyMesh2d), intent(in) :: self
	integer(kint), intent(in) :: initNest
	!
	integer(kint) :: i, steps, result
	
	result = 0
 	if ( self%meshSeed == TRI_HEX_SEED ) then
		do i = 2**initNest + 1, 2**(initNest+1)
			result = result + i
		enddo
		result = 2*result + 2**(initNest+1)+1
	elseif ( self%meshSeed == QUAD_RECT_SEED .OR. self%meshSeed ==  QUAD_RECT_PERIODIC_SEED) then
		result = 3
		do i = 1, initNest
			result = result + 2**i
		enddo
		result = result*result
	elseif ( self%meshSeed == ICOS_TRI_SPHERE_SEED) then
			result = 2 + 10 * 4 ** initNest
	elseif ( self%meshSeed == CUBED_SPHERE_SEED) then
			result = 2 + 6 * 4 ** initNest
	endif
	nVerticesInMesh = result
end function

function nFacesInMesh( self, initNest )
	integer(kint) :: nFacesInMesh
	type(PolyMesh2d), intent(in) :: self
	integer(kint), intent(in) :: initNest
	
	nFacesInMesh = 0
	select case (self%meshSeed)
		case (TRI_HEX_SEED)
			nFacesInMesh = 6 * 4 ** initNest
		case (QUAD_RECT_SEED, QUAD_RECT_PERIODIC_SEED)
			nFacesInMesh = 4 * 4 ** initNest
		case (ICOS_TRI_SPHERE_SEED)
			nFacesInMesh = 20 * 4 ** initNest
		case (CUBED_SPHERE_SEED)
			nFacesInMesh = 6 * 4 ** initNest
	end select
end function

function nEdgesInMesh( self, nVertices, nFaces )
	integer(kint) :: nEdgesInMesh
	type(PolyMesh2d), intent(in) :: self
	integer(kint), intent(in) :: nVertices
	integer(kint), intent(in) :: nFaces
	nEdgesInMesh = 0
	select case (self%geomKind)
		case (PLANAR_GEOM)
			nEdgesInMesh = nFaces + nVertices - 1
		case (SPHERE_GEOM)
			nEdgesInMesh = nFaces + nVertices - 2
	end select
end function

subroutine initializeMeshFromSeed( self, ampFactor )
	type(PolyMesh2d), intent(inout) :: self
	real(kreal), intent(in) :: ampFactor
	!
	character(len=56) :: seedFilename
	integer(kint) :: i
	integer(kint) :: nSeedParticles, nSeedEdges, nSeedFaces, nSeedVerts
	real(kreal), allocatable :: seedXYZ(:,:)
	integer(kint), allocatable :: seedEdgeOrigs(:), seedEdgeDests(:), seedEdgeLefts(:), seedEdgeRights(:)
	integer(kint), allocatable :: seedFaceVerts(:,:), seedFaceEdges(:,:), seedVertexDegree(:)
	
	select case ( self%meshSeed )
		case (TRI_HEX_SEED)
			nSeedParticles = 13
			nSeedEdges = 12
			nSeedFaces = 6
			nSeedVerts = 7
!			self%faceKind = TRI_PANEL
!			self%geomKind = PLANAR_GEOM
			seedFilename = "triHexSeed.dat"
		case (QUAD_RECT_SEED)
			nSeedParticles = 13
			nSeedEdges = 12
			nSeedFaces = 4
			nSeedVerts = 9
!			self%faceKind = QUAD_PANEL
!			self%geomKind = PLANAR_GEOM
			seedFilename = "quadRectSeed.dat"
		case (POLAR_DISC_SEED)
			call LogMessage(log, ERROR_LOGGING_LEVEL,logkey," initMeshFromSeed ERROR : seed not implemented.")
			return
		case (QUAD_RECT_PERIODIC_SEED)
			nSeedParticles = 13
			nSeedEdges = 12
			nSeedFaces = 4
			nSeedVerts = 9
			call LogMessage(log, ERROR_LOGGING_LEVEL,logkey," initMeshFromSeed ERROR : seed not implemented.")
			return
		case (ICOS_TRI_SPHERE_SEED)
			nSeedParticles = 32
			nSeedEdges = 30
			nSeedFaces = 20
			nSeedVerts = 12
!			self%faceKind = TRI_PANEL
!			self%geomKind = SPHERE_GEOM
			seedFilename = "icosTriSeed.dat"
		case (CUBED_SPHERE_SEED)
			nSeedParticles = 14
			nSeedEdges = 12
			nSeedFaces = 6
			nSeedVerts = 8
!			self%faceKind = QUAD_PANEL
!			self%geomKind = SPHERE_GEOM
			seedFilename = "cubedSphereSeed.dat"
	end select
	
	
	allocate(seedXYZ(3,nSeedParticles))
	allocate(seedEdgeOrigs(nSeedEdges))
	allocate(seedEdgeDests(nSeedEdges))
	allocate(seedEdgeLefts(nSeedEdges))
	allocate(seedEdgeRights(nSeedEdges))
	allocate(seedVertexDegree(nSeedVerts))
	if ( self%faceKind == TRI_PANEL ) then
		allocate(seedFaceVerts(3,nSeedFaces))
		allocate(seedFaceEdges(3,nSeedFaces))
	elseif ( self%faceKind == QUAD_PANEL ) then
		allocate(seedFaceVerts(4,nSeedFaces))
		allocate(seedFaceEdges(4,nSeedFaces))
	endif
	
	call readSeedFile(self, seedFilename, nSeedParticles, nSeedEdges, nSeedFaces, seedXYZ, seedEdgeOrigs, seedEdgeDests, &
					  seedEdgeLefts, seedEdgeRights, seedFaceVerts, seedFaceEdges, seedVertexDegree)
					  
	if ( self%meshSeed == ICOS_TRI_SPHERE_SEED ) then
		do i = 1, 20
			seedXYZ(:,12+i) = SphereTriCenter( seedXYZ(:,seedFaceVerts(1,i)), seedXYZ(:,seedFaceVerts(2,i)), &
											   seedXYZ(:,seedFaceVerts(3,i)) )
		enddo
	endif
	
	seedXYZ = ampFactor * seedXYZ
	
	!
	!	initialize mesh with seed data
	!
	do i = 1, nSeedParticles
		call InsertParticle( self%particles, seedXYZ(:,i), seedXYZ(:,i))
	enddo
	if ( self%particles%N /= nSeedParticles ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" initMeshFromSeed ERROR : "," particles%N.")
	endif
	do i = 1, nSeedEdges
		call InsertEdge( self%Edges, self%particles, seedEdgeOrigs(i), seedEdgeDests(i), seedEdgeLefts(i), seedEdgeRights(i))
	enddo
	if ( self%edges%N /= nSeedEdges ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" initMeshFromSeed ERROR : "," edges%N.")
	endif
	do i = 1, nSeedFaces
		call InsertFace( self%Faces, nSeedVerts + i, seedFaceVerts(:,i), seedFaceEdges(:,i))
		call MakeParticleActive( self%particles, nSeedVerts +i)
	enddo
	if ( self%faces%N /= nSeedFaces ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, logkey//" initMeshFromSeed ERROR : "," faces%N.")
	endif
	self%faces%N_Active = nSeedFaces
	if ( self%faceKind == TRI_PANEL ) then
		do i = 1, nSeedFaces
			self%faces%area(i) = TriFaceArea( self%faces, i, self%particles)
			self%particles%area( self%faces%centerParticle(i) ) = self%faces%area(i)
		enddo
	elseif ( self%faceKind == QUAD_PANEL ) then
		do i = 1, nSeedFaces
			self%faces%area(i) = QuadFaceArea( self%faces, i, self%particles)
			self%particles%area( self%faces%centerParticle(i) ) = self%faces%area(i)
		enddo
	endif
	!
	!	initialize dual information
	!
	do i = 1, nSeedVerts
		if ( self%particles%nEdges(i) /= seedVertexDegree(i) ) then
			write(logstring,*) " dual connectivity error at vertex ", i
			call LogMessage(log, WARNING_LOGGING_LEVEL,logkey//" initFromSeed WARNING: ",logstring)
		endif
		call SortIncidentEdgesAtParticle( self%particles, i )
	enddo

	deallocate(seedXYZ)
	deallocate(seedEdgeOrigs)
	deallocate(seedEdgeDests)
	deallocate(seedEdgeLefts)
	deallocate(seedEdgeRights)
	deallocate(seedFaceVerts)
	deallocate(seedFaceEdges)
end subroutine

subroutine readSeedFile( self, seedFilename, nParticles, nEdges, nFaces, seedXYZ, &
						 seedEdgeOrigs, seedEdgeDests, seedEdgeLefts, seedEdgeRights, &
						 seedFaceVerts, seedFaceEdges, seedVertexDegree )
	type(PolyMesh2d), intent(in) :: self
	character(len=*), intent(in) :: seedFilename
	integer(kint), intent(in) :: nParticles, nEdges, nFaces
	real(kreal), intent(out) :: seedXYZ(3,nParticles)
	integer(kint), intent(out) :: seedEdgeOrigs(nEdges)
	integer(kint), intent(out) :: seedEdgeDests(nEdges)
	integer(kint), intent(out) :: seedEdgeLefts(nEdges)
	integer(kint), intent(out) :: seedEdgeRights(nEdges)
	integer(kint), intent(inout) :: seedFaceVerts(:,:)
	integer(kint), intent(inout) :: seedFaceEdges(:,:)
	integer(kint), intent(out) :: seedVertexDegree(nParticles - nFaces)
	!
	character(len=56) :: tempString							 
	integer(kint) :: i, readStat
	
	open(unit=READ_UNIT,file=seedFilename,STATUS='OLD',ACTION='READ',iostat=readStat)
	if ( readStat /= 0 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey," cannot open seed file.")
		return
	endif
	!
	!	read root particle positions
	!
	read(READ_UNIT,*) tempString
	if ( self%geomKind == PLANAR_GEOM ) then
		do i = 1, nParticles
			read(READ_UNIT,*) seedXYZ(1,i), seedXYZ(2,i)
		enddo
	else
		if ( self%meshSeed == ICOS_TRI_SPHERE_SEED ) then
			do i = 1, 12
				read(READ_UNIT,*) seedXYZ(1,i), seedXYZ(2,i), seedXYZ(3,i)
			enddo
		else
			do i = 1, nParticles
				read(READ_UNIT,*) seedXYZ(1,i), seedXYZ(2,i), seedXYZ(3,i)
			enddo
		endif
	endif
	!
	!	read root edges
	!
	read(READ_UNIT,*) tempString
	do i = 1, nEdges
		read(READ_UNIT,*) seedEdgeOrigs(i), seedEdgeDests(i), seedEdgeLefts(i), seedEdgeRights(i)
	enddo
	! 	adjust to fortran (base index = 1, not 0)
	seedEdgeOrigs = seedEdgeOrigs + 1
	seedEdgeDests = seedEdgeDests + 1
	seedEdgeLefts = seedEdgeLefts + 1
	seedEdgeRights = seedEdgeRights + 1
	!
	!	read root faces
	!
	read(READ_UNIT,*) tempString
	if ( self%faceKind == TRI_PANEL ) then
		do i = 1, nFaces
			read(READ_UNIT,*) seedFaceVerts(1,i), seedFaceVerts(2,i), seedFaceVerts(3,i)
		enddo
		read(READ_UNIT,*) tempString
		do i =1, nFaces
			read(READ_UNIT,*) seedFaceEdges(1,i), seedFaceEdges(2,i), seedFaceEdges(3,i)
		enddo
	elseif ( self%faceKind == QUAD_PANEL ) then
		do i = 1, nFaces
			read(READ_UNIT,*) seedFaceVerts(1,i), seedFaceVerts(2,i), seedFaceVerts(3,i), seedFaceVerts(4,i)
		enddo
		read(READ_UNIT,*) tempString
		do i =1, nFaces
			read(READ_UNIT,*) seedFaceEdges(1,i), seedFaceEdges(2,i), seedFaceEdges(3,i), seedFaceEdges(4,i)
		enddo
	endif
	!	adjust to fortran (base index = 1)
	seedFaceEdges = seedFaceEdges + 1
	seedFaceVerts = seedFaceVerts + 1
	!
	!	read vertex degree
	!
	read(READ_UNIT,*) tempString
	do i = 1, nParticles - nFaces
		read(READ_UNIT,*) seedVertexDegree(i)
	enddo
	
	close(READ_UNIT)
end subroutine

subroutine WriteMeshToMatlab( self, fileunit )
	type(PolyMesh2d), intent(in) :: self
	integer(kint), intent(in) :: fileunit
	
	call WriteParticlesToMatlab(self%particles, fileunit)
	call WriteEdgesToMatlab(self%edges, fileunit)
	call WriteFacesToMatlab(self%faces, fileunit)
end subroutine

subroutine LogStatsPrivate(self, aLog)
	type(PolyMesh2d), intent(in) :: self
	type(Logger), intent(inout) :: aLog
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, logKey, "PolyMesh Stats:" )
	call StartSection(aLog)
	call MeshSeedString(logString, self%meshSeed)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "meshSeed = ", trim(logString))
	if ( self%geomKind == PLANAR_GEOM ) then
		call LogMessage(aLog, TRACE_LOGGING_LEVEL, "geomKind = ", "PLANAR_GEOM")
	elseif ( self%geomKind == SPHERE_GEOM ) then
		call LogMessage(aLog, TRACE_LOGGING_LEVEL, "geomKind = ", "SPHERE_GEOM")
	elseif ( self%geomKind == EUCLIDEAN_3D ) then
		call LogMessage(aLog, TRACE_LOGGING_LEVEL, "geomKind = ", "EUCLIDEAN_3D")
	else
		call LogMessage(aLog, TRACE_LOGGING_LEVEL, "geomKind = ", "invalid geomKind.")
	endif
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "initNest = ", self%initNest)
	call LogMessage(aLog, TRACE_LOGGING_LEVEL, "amrLimit = ", self%amrLimit)
	call LogStats(self%particles, aLog)
	call LogStats(self%edges, aLog)
	call LogStats(self%faces, aLog)
	call EndSection(aLog)
end subroutine

subroutine PrintDebugPrivate( self ) 
	type(PolyMesh2d), intent(in) :: self
	print *, "PolyMesh2d DEBUG info : "
	print *, "meshSeed = ", self%meshSeed
	print *, "geomKind = ", self%geomKind
	print *, "initNest = ", self%initNest
	print *, "amrLimit = ", self%amrLimit
	call PrintDebugInfo( self%particles )
	call PrintDebugInfo( self%edges )
	call PrintDebugInfo( self%faces )
end subroutine

subroutine MeshSeedString( seedString, msInt )
	character(len=*), intent(inout) :: seedString
	integer(kint), intent(in) :: msInt
	select case (msInt)
		case (TRI_HEX_SEED )
			write(seedString,*) "TRI_HEX_SEED"
		case (QUAD_RECT_SEED)
			write(seedString,*) "QUAD_RECT_SEED"
		case (QUAD_RECT_PERIODIC_SEED)
			write(seedString,*) "QUAD_RECT_PERIODIC_SEED"
		case (POLAR_DISC_SEED)
			write(seedString,*) "POLAR_DISC_SEED"
		case (ICOS_TRI_SPHERE_SEED)
			write(seedString,*) "ICOS_TRI_SPHERE_SEED"
		case (CUBED_SPHERE_SEED)
			write(seedString,*) "CUBED_SPHERE_SEED"
		case default
			write(seedString,*) "invalid seed"
	end select
end subroutine

subroutine InitLogger(aLog,rank)
! Initialize a logger for this module and processor
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.3,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

end module