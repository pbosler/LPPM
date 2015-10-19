module TrivariateMLSQModule

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use STDIntVectorModule
use ParticlesModule
use EdgesModule
use FacesModule
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use SphereGeomModule

implicit none

include 'mpif.h'

private
public QuadMLSQ, New, Delete, LogStats
public InterpolateScalar, InterpolateVector !, InterpolateScalarGradient, InterpolateScalarLaplacian
!public MLSQScalarGradientAtParticles, MLSQLaplacianAtParticles, MLSQDoubleDotAtParticles
public GetParticlesNearFace, GetParticlesNearVertex
!
!----------------
! types and module variables
!----------------
!
type QuadMLSQ
	real(kreal), pointer :: coeffs1(:,:) => null()
	real(kreal), pointer :: coeffs2(:,:) => null()
	real(kreal), pointer :: coeffs3(:,:) => null()
	type(STDIntVector), pointer :: nearbyParticles(:) => null()
	contains
		final :: deletePrivate
end type

integer(kint), parameter :: nCoeffs = 10 ! trivariate quadratic polynomials have 10 coefficients

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

interface LogStats
	module procedure logStatsPrivate
end interface

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'MLSQ'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString

contains
!
!----------------
! Public methods
!----------------
!
subroutine newPrivate(self, aMesh, dataField)
	type(QuadMLSQ), intent(out) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: dataField
	integer(kint) :: i
	
	if ( .NOT. logInit ) call InitLogger(log, procRank)
	
	allocate(self%coeffs1(nCoeffs, aMesh%particles%N))
	self%coeffs1 = 0.0_kreal
	if ( dataField%nDim > 1 ) then
		allocate(self%coeffs2(nCoeffs, aMesh%particles%N))
		self%coeffs2 = 0.0_kreal
		if ( dataField%nDim == 3 ) then
			allocate(self%coeffs3(nCoeffs, aMesh%particles%N))
			self%coeffs3 = 0.0_kreal
		endif
	endif
	allocate(self%nearbyParticles(amesh%particles%N))
	
	do i = 1, aMesh%particles%N
		if ( .NOT. aMesh%particles%isActive(i) ) then
			call GetParticlesNearVertex( self%nearbyParticles(i), aMesh, i)
		endif
	enddo
	
	do i = 1, aMesh%faces%N
		if ( .NOT. aMesh%faces%hasChildren(i) ) then
			call GetParticlesNearFace( self%nearbyParticles( aMesh%faces%centerParticle(i)), aMesh, i)
		endif
	enddo
	
	call SetCoefficients(self, aMesh, dataField)
end subroutine

subroutine deletePrivate(self)
	type(QuadMLSQ), intent(inout) :: self
	integer(kint) :: i
	if ( associated(self%nearbyParticles) ) then
		deallocate(self%nearbyParticles)
	endif
	if ( associated(self%coeffs1) ) deallocate(self%coeffs1)
	if ( associated(self%coeffs2) ) deallocate(self%coeffs2)
	if ( associated(self%coeffs3) ) deallocate(self%coeffs3)
end subroutine

subroutine logStatsPrivate(self, aLog)
	type(QuadMLSQ), intent(in) :: self
	type(Logger), intent(inout) :: aLog
	!
	integer(kint) :: i
	
	call StartSection(aLog, "QuadMLSQ nearbyParticles:")
	do i = 1, size(self%nearbyParticles)
		write(logString,'(I6,A)') i, " :"
		call LogMessage(aLog, TRACE_LOGGING_LEVEL, trim(logString), self%nearbyParticles(i) )
	enddo
	call EndSection(aLog)
end subroutine

subroutine SetCoefficients(self, aMesh, aField)
	type(QuadMLSQ), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: aField
	!
	real(kreal) :: x0(3)
	real(kreal) :: dataLocs(3,31)
	real(kreal) :: A(31,10) ! GetParticlesNear* will return vectors of size <= 29
	real(kreal) :: b(31,3)
	integer(kint) :: mm
	integer(kint) :: i, j
	real(kreal) :: lpkWork(20)
	integer(kint) :: lpkLWork, lpkInfo
	
	lpkLWork = 20
	
	do i = 1, aMesh%particles%N
		dataLocs(:,1) = PhysCoord(aMesh%particles, self%nearbyParticles(i)%int(1))
		mm = self%nearbyParticles(i)%N + 2
		if ( aField%nDim == 1 ) then
			b(1,1) = aField%scalar( self%nearbyParticles(i)%int(1))
		else
			b(1,1) = aField%xComp( self%nearbyParticles(i)%int(1))
			b(1,2) = aField%yComp( self%nearbyParticles(i)%int(1))
			b(1,3) = aField%zComp( self%nearbyParticles(i)%int(1))
		endif
		
		do j = 2, self%nearbyParticles(i)%N
			dataLocs(:,j) = PhysCoord(aMesh%particles, self%nearbyParticles(i)%int(j))
			if ( aField%nDim == 1 ) then
				b(j,1) = aField%scalar(self%nearbyParticles(i)%int(j))
			else
				b(j,1) = aField%xComp(self%nearbyParticles(i)%int(j))
				b(j,2) = aField%yComp(self%nearbyParticles(i)%int(j))
				b(j,3) = aField%zComp(self%nearbyParticles(i)%int(j))
			endif
		enddo
		
		dataLocs(:, mm-1) = 0.0_kreal
		b(mm-1,:) = b(1,:)
		dataLocs(:,mm) = 2.0_kreal * dataLocs(:,1)
		b(mm,:) = b(1,:)
		
		do j = 1, mm
			dataLocs(:,j) = dataLocs(:,j) - dataLocs(:,1)
		enddo
				
		A(1:mm,1) = 1.0_kreal
		A(1:mm,2) = dataLocs(1,1:mm)
		A(1:mm,3) = dataLocs(2,1:mm)
		A(1:mm,4) = dataLocs(3,1:mm)
		A(1:mm,5) = dataLocs(1,1:mm) * dataLocs(2,1:mm)
		A(1:mm,6) = dataLocs(1,1:mm) * dataLocs(3,1:mm)
		A(1:mm,7) = dataLocs(2,1:mm) * dataLocs(3,1:mm)
		A(1:mm,8) = dataLocs(1,1:mm) * dataLocs(1,1:mm)
		A(1:mm,9) = dataLocs(2,1:mm) * dataLocs(2,1:mm)
		A(1:mm,10) = dataLocs(3,1:mm) * dataLocs(3,1:mm)

		if ( aField%nDim == 1 ) then
			!call DGELSY( mm, 10, 1, A, mm, b(1:mm,1), mm, lpkWork, lpkLWork, lpkInfo	)
			call wrapped_scalar_dgelsy( mm, 10, A(1:mm,:), b(1:mm,1) )
			self%coeffs1(:,i) = b(1:10,1)
		else
			call wrapped_vector_dgelsy(mm, 10, A(1:mm,:), b(1:mm,:))
			self%coeffs1(:,i) = b(1:10,1)
			self%coeffs2(:,i) = b(1:10,2)
			self%coeffs3(:,i) = b(1:10,3)
		endif
	enddo
end subroutine

!
!----------------
! Private methods
!----------------
!
subroutine wrapped_scalar_dgelsy( m, n, A, b )
	integer(kint), intent(in) :: m, n
	real(kreal), intent(inout), dimension(m,n) :: A
	real(kreal), intent(inout), dimension(m) :: b
	!
	integer(kint) :: nRhs, lda, ldb, rank, lwork, info
	integer(kint), dimension(n) :: jpvt
	real(kreal) :: rcond
	real(kreal), dimension(100) :: work
	real(kreal), dimension(n) :: svals
	
	nRhs = 1
	lda = m
	ldb = m
	lwork = 100
	jpvt = 0
	rcond = -1.0_kreal
	
!	call DGELSY(m, n, nRhs, A, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
	call DGELSS(m, n, nRhs, A, lda, b, ldb, svals, rcond, rank, work, lwork, info)
	
	if ( info < 0 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL,trim(logkey)//" DGELSY ERROR : bad argument ", info)
	elseif (info > 0 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL,trim(logKey)//" DGELSS ERROR : failed to converge", info)
	endif
end subroutine

subroutine wrapped_vector_dgelsy( m, n, A, b )
	integer(kint), intent(in) :: m, n
	real(kreal), intent(inout), dimension(m,n) :: A
	real(kreal), intent(inout), dimension(m,3) :: b
	!	
	integer(kint) :: nRhs, lda, ldb, rank, lwork, info
	integer(kint), dimension(n) :: jpvt
	real(kreal) :: rcond
	real(kreal), dimension(4*n+1) :: work
	
	nRhs = 3
	lda = m
	ldb = m
	lwork = 4*n+1
	jpvt = 0
	
	call DGELSY(m, n, nRhs, A, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
	
	if ( info < 0 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL,trim(logkey)//" DGELSY ERROR : bad argument ", info)
	endif
end subroutine

function InterpolateScalar( self, aMesh, queryPt )
	real(kreal) :: InterpolateScalar
	type(QuadMLSQ), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	real(kreal), dimension(3), intent(in) :: queryPt
	!
	integer(kint) :: np 
	real(kreal), dimension(10) :: basis
	real(kreal), dimension(3) :: x0
	
	np = nearestParticle( aMesh, queryPt )
	
	x0 = PhysCoord(aMesh%particles, nP)
	
	basis(1) = 1.0_kreal
	basis(2) = queryPt(1) - x0(1)
	basis(3) = queryPt(2) - x0(2)
	basis(4) = queryPt(3) - x0(3)
	basis(5) = (queryPt(1) - x0(1)) * (queryPt(2) - x0(2))
	basis(6) = (queryPt(1) - x0(1)) * (queryPt(3) - x0(3))
	basis(7) = (queryPt(2) - x0(2)) * (queryPt(3) - x0(3))
	basis(8) = (queryPt(1) - x0(1)) * (queryPt(1) - x0(1))
	basis(9) = (queryPt(2) - x0(2)) * (queryPt(2) - x0(2))
	basis(10) = (queryPt(3) - x0(3)) * (queryPt(3) - x0(3))
	
	InterpolateScalar = sum( self%coeffs1(:,nP) * basis)
end function

function InterpolateVector( self, aMesh, queryPt )
	real(kreal), dimension(3) :: InterpolateVector
	type(QuadMLSQ), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	real(kreal), dimension(3), intent(in) :: queryPt
	!
	integer(kint) :: np 
	real(kreal), dimension(10) :: basis
	real(kreal), dimension(3) :: x0
	
	np = nearestParticle( aMesh, queryPt )
	
	x0 = PhysCoord(aMesh%particles, nP)
	
	basis(1) = 1.0_kreal
	basis(2) = queryPt(1) - x0(1)
	basis(3) = queryPt(2) - x0(2)
	basis(4) = queryPt(3) - x0(3)
	basis(5) = (queryPt(1) - x0(1)) * (queryPt(2) - x0(2))
	basis(6) = (queryPt(1) - x0(1)) * (queryPt(3) - x0(3))
	basis(7) = (queryPt(2) - x0(2)) * (queryPt(3) - x0(3))
	basis(8) = (queryPt(1) - x0(1)) * (queryPt(1) - x0(1))
	basis(9) = (queryPt(2) - x0(2)) * (queryPt(2) - x0(2))
	basis(10) = (queryPt(3) - x0(3)) * (queryPt(3) - x0(3))
	
	InterpolateVector(1) = sum( self%coeffs1(:,np) * basis)
	InterpolateVector(2) = sum( self%coeffs2(:,nP) * basis)
	InterpolateVector(3) = sum( self%coeffs3(:,nP) * basis)
end function

subroutine GetParticlesNearFace( nearbyParticles, aMesh, faceIndex)
	type(STDIntVector), intent(out) :: nearbyParticles
	type(PolyMesh2d), intent(in) :: aMesh
	integer(kint), intent(in) :: faceIndex
	!
	type(STDIntVector) :: adjFaces, faceVerts
	integer(kint) :: i, j
	logical(klog) :: duplicateFound
	
	call initialize(nearbyParticles)
	call nearbyParticles%pushBack( aMesh%faces%centerParticle(faceIndex) )
	
	call CCWAdjacentFaces( aMesh, adjFaces, faceIndex )
	! DEBUG 
	duplicateFound = .FALSE.
	do i = 1, adjFaces%N
		do j = i+1, adjFaces%N
			if ( adjFaces%int(i) == adjFaces%int(j) ) duplicateFound = .TRUE.
		enddo
	enddo
	
	if ( duplicateFound ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey)//"GetParticlesNearFace adjFaces duplicate ERROR at face : ", faceIndex )
		call LogMessage(log, ERROR_LOGGING_LEVEL, "adjFaces = ", adjFaces )
	endif
	! END DEBUG
	
	call CCWVerticesAroundFace( aMesh, faceVerts, faceIndex )
	! DEBUG 
	duplicateFound = .FALSE.
	do i = 1, faceVerts%N
		do j = i+1, faceVerts%N
			if ( faceVerts%int(i) == faceVerts%int(j) ) duplicateFound = .TRUE.
		enddo
	enddo
	
	if ( duplicateFound ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey)//" GetParticlesNearFace faceVerts duplicate ERROR at face : ", faceIndex )
		call LogMessage(log, ERROR_LOGGING_LEVEL, "faceVerts = ", faceVerts )
	endif
	! END DEBUG

	do i = 1, adjFaces%N
		call nearbyParticles%pushBack( aMesh%faces%centerParticle( adjFaces%int(i) ))
	enddo
	do i = 1, faceVerts%N
		call nearbyParticles%pushBack( faceVerts%int(i) )
	enddo
	
	do i = 1, adjFaces%N
		call initialize(faceVerts)
		call CCWVerticesAroundFace( aMesh, faceVerts, adjFaces%int(i) )
		do j = 1, faceVerts%N
			call nearbyParticles%pushBackUnique( faceVerts%int(j))
		enddo
	enddo
end subroutine

subroutine GetParticlesNearVertex( nearbyParticles, aMesh, vertIndex)
	type(STDIntVector), intent(out) :: nearbyParticles
	type(PolyMesh2d), intent(in) :: aMesh
	integer(kint), intent(in) :: vertIndex
	!
	type(STDIntVector) :: adjFaces, faceVerts
	integer(kint) :: i, j
	logical(klog) :: duplicateFound
	
	call initialize(nearbyParticles)
	call nearbyParticles%pushBack(vertIndex)	
	
	call CCWFacesAroundVertex( aMesh, adjFaces, vertIndex )
	! DEBUG
	duplicateFound = .FALSE.
	do i = 1, adjFaces%N
		do j = i + 1, adjFaces%N
			if ( adjFaces%int(i) == adjFaces%int(j) ) duplicateFound = .TRUE.
		enddo
	enddo
	if ( duplicateFound ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, &
			trim(logKey)//" GetParticlesNearVertex adjFaces duplicate ERROR at vert", vertIndex)
		call LogMessage(log, ERROR_LOGGING_LEVEL, "adjFaces = ", adjFaces )
	endif
	! END DEBUG
	
	do i = 1, adjFaces%N
		call nearbyParticles%pushBack( aMesh%faces%centerParticle(adjFaces%int(i)) )
	enddo
	
	do i = 1, adjFaces%N
		call initialize(faceVerts)
		call CCWVerticesAroundFace( aMesh, faceVerts, adjFaces%int(i))
		do j = 1, faceVerts%N
			call nearbyParticles%pushBackUnique( faceVerts%int(j))
		enddo
	enddo
end subroutine


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