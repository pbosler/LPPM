module RefinementModule

use NumberKindsModule
use STDIntVectorModule
use OutputWriterModule
use LoggerModule
use ParticlesModule
use FieldModule
use EdgesModule
use FacesModule
use PlaneGeomModule
use SphereGeomModule
use PolyMesh2dModule

implicit none

private
public RefineSetup, New, Delete
!----------------
! types and module variables
!----------------
!
type RefineSetup
	logical(klog), pointer :: refineFlag(:) => null()
	contains
		final :: deletePrivate
end type
!
!----------------
! interfaces
!----------------
!
interface New
	module procedure newPrivate
end interface

interface Delete
	module procedure deletePrivate
end interface

interface
	function FlagFunction( mesh, dataField, tolerance, faceIndex )
		use NumberKindsModule
		use PolyMesh2dModule
		use FieldModule
		logical(klog) :: FlagFunction
		type(PolyMesh2d), intent(in) :: mesh
		type(Field), intent(in) :: dataField
		real(kreal), intent(in) :: tolerance
		integer(kint), intent(in) :: faceIndex
	end function 
end interface

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'Refine'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains

!
!----------------
! public methods
!----------------
!
subroutine newPrivate(self, nFaces )
	type(RefineSetup), intent(out) :: self
	integer(kint), intent(in) :: nFaces
	
	if (.NOT. logInit ) call InitLogger(log, procRank)
	
	allocate(self%refineFlag(nFaces))
	self%refineFlag = .FALSE.
	
end subroutine

subroutine deletePrivate(self)
	type(RefineSetup), intent(inout) :: self
	if ( associated(self%refineFlag)) deallocate(self%refineFlag)
end subroutine

subroutine DivideFlaggedFaces( self, mesh, indexStart, indexEnd, nParticlesBefore, nParticlesAfter )
	type(RefineSetup), intent(in) :: self
	type(PolyMesh2d), intent(inout) :: mesh
	integer(kint), intent(in) :: indexStart, indexEnd
	integer(kint), intent(out) :: nParticlesBefore, nParticlesAfter
	!
	integer(kint) :: i
	
	nParticlesBefore = mesh%particles%N
	do i = indexStart, indexEnd
		if ( self%refineFlag(i) .AND. mesh%initNest + mesh%amrLimit > CountParents(mesh%faces, i) ) then
			if ( mesh%faceKind == TRI_PANEL) then
				call DivideTriFace( mesh%faces, i, mesh%particles, mesh%edges)
			elseif ( mesh%faceKind == QUAD_PANEL ) then
				call DivideQuadFace( mesh%faces, i, mesh%particles, mesh%edges)
			endif
			self%refineFlag(i) = .FALSE.
		endif
	enddo
	nParticlesAfter = mesh%particles%N
end subroutine

subroutine FlagFaces( self, mesh, dataField, flagFn, tol, startIndex, endIndex, counter )
	type(RefineSetup), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: mesh
	type(Field), intent(in) :: dataField
	procedure(FlagFunction) :: flagFn
	real(kreal), intent(in) :: tol
	integer(kint), intent(in) :: startIndex, endIndex
	integer(kint), intent(inout) :: counter
	!
	integer(kint) :: i 
	
	do i = startIndex, endIndex
		if ( .NOT. mesh%faces%hasChildren(i) ) then
			if ( flagFn( mesh, dataField, tol, i) ) then	
				self%refineFlag(i) = .TRUE.
				counter = counter + 1
			endif
		endif
	enddo
end subroutine

function ScalarIntegralRefinement( mesh, dataField, tol, faceIndex )
	logical(klog) :: ScalarIntegralRefinement
	type(PolyMesh2d), intent(in) :: mesh
	type(Field), intent(in) :: dataField
	real(kreal), intent(in) :: tol
	integer(kint), intent(in) :: faceIndex
	!
	integer(kint) :: pIndex
	pIndex = mesh%faces%centerParticle(faceIndex)
	ScalarIntegralRefinement = ( abs(dataField%scalar(pIndex)) * mesh%particles%area(pIndex) > tol )
end function

function ScalarVariationRefinement( mesh, dataField, tol, faceIndex )
	logical(klog) :: ScalarVariationRefinement
	type(PolyMesh2d), intent(in) :: mesh
	type(Field), intent(in) :: dataField
	real(kreal), intent(in) :: tol
	integer(kint), intent(in) :: faceIndex
	!
	integer(kint) :: pIndex, i
	real(kreal) :: maxScalar, minScalar
	type(STDIntVector) :: faceVerts
	
	pIndex = mesh%faces%centerParticle(faceIndex)
	call CCWVerticesAroundFace( mesh, faceVerts, faceindex)
	maxScalar = dataField%scalar(pIndex)
	minScalar = maxScalar
	do i = 1, faceVerts%N
		if ( dataField%scalar( faceVerts%int(i) ) > maxScalar ) maxScalar = dataField%scalar( faceVerts%int(i))
		if ( dataField%scalar( faceVerts%int(i) ) < minScalar ) minScalar = dataField%scalar( faceVerts%int(i))
	enddo
	ScalarVariationRefinement = ( maxScalar - minScalar > tol )
end function

!
!----------------
! private methods
!----------------
!
subroutine InitLogger(aLog,rank)
! Initialize a logger for this module and processor
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	if ( rank == 0 ) then
		call New(aLog,logLevel)
	else
		call New(aLog,ERROR_LOGGING_LEVEL)
	endif
	logInit = .TRUE.
end subroutine

end module