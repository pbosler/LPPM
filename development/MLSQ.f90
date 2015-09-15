module MLSQModule

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
!use LAPACK95

implicit none

include 'mpif.h'

private
public MLSQ, New, Delete
public InitializeMLSQ
public InterpolateScalar, InterpolateVector, InterpolateScalarGradient, InterpolateScalarLaplacian
public MLSQGradientAtParticles, MLSQLaplacianAtParticles, MLSQDoubleDotProductAtParticles
!
!----------------
! types and module variables
!----------------
!
type MLSQ
	real(kreal), pointer :: basisVectors(:,:,:) => null()
	real(kreal), pointer :: coeffs1(:,:) => null()
	real(kreal), pointer :: coeffs2(:,:) => null()
	real(kreal), pointer :: coeffs3(:,:) => null()

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

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'MLSQ'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains
!
!----------------
! public methods
!----------------
!
subroutine newPrivate( self, aMesh, dataField )
	type(MLSQ), intent(out) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: dataField
	!
	integer(kint) :: nCoeffs
	
	if ( .NOT. logInit ) call InitLogger(log, procRank)
	
	allocate(self%basisVectors(3,3,aMesh%particles%N))
	self%basisVectors = 0.0_kreal
	
	nCoeffs = 10 ! cubic bivariate polynomials
	
	allocate(self%coeffs1( nCoeffs, aMesh%particles%N))
	self%coeffs1 = 0.0_kreal
	if ( dataField%nDim == 2 .OR. dataField%nDim == 3) then
		allocate(self%coeffs2(nCoeffs, aMesh%particles%N))
		self%coeffs2 = 0.0_kreal
	endif
	if ( dataField%nDim == 3 ) then
		allocate(self%coeffs3(nCoeffs, aMesh%particles%N))
		self%coeffs3 = 0.0_kreal
	endif
	
	call InitializeMLSQ( self, aMesh, dataField)
end subroutine

subroutine deletePrivate(self)
	type(MLSQ), intent(inout) :: self
	if ( associated(self%basisVectors)) deallocate(self%basisVectors)
	if ( associated(self%coeffs1)) deallocate(self%coeffs1)
	if ( associated(self%coeffs2)) deallocate(self%coeffs2)
	if ( associated(self%coeffs3)) deallocate(self%coeffs3)
end subroutine

subroutine InitializeMLSQ(self, aMesh, dataField)
	type(MLSQ), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: dataField
	!
	type(STDIntVector) :: nearbyParticles
	integer(kint) :: i
	real(kreal) :: localCoeffs(10,3)
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL, trim(logKey), " entering InitializeMLSQ ... ")
	call StartSection(log, "DEBUGGING MLSQ:")
	do i = 1, aMesh%particles%N
		if ( .NOT. aMesh%particles%isActive(i) ) then
			
			call LogMessage(log, DEBUG_LOGGING_LEVEL, "particle i = ", i )
			
			call GetParticlesNearVertex(nearbyParticles, aMesh, i )
			
			call nearbyParticles%print("nearbyParticles")
			
			call SetBasisVectorsAndCoefficients( self%basisVectors(:,:,i), localCoeffs, nearbyParticles, aMesh, dataField)
			
			if ( dataField%nDim == 1) then
				self%coeffs1(:,i) = localCoeffs(:,1)
			elseif ( dataField%nDim == 2) then
				self%coeffs1(:,i) = localCoeffs(:,1)
				self%coeffs2(:,i) = localCoeffs(:,2)
			elseif ( dataField%nDim == 3) then
				self%coeffs1(:,i) = localCoeffs(:,1)
				self%coeffs2(:,i) = localCoeffs(:,2)
				self%coeffs3(:,i) = localCoeffs(:,3)
			endif		
		endif
	enddo
	do i = 1, aMesh%faces%N
		if ( .NOT. aMesh%faces%hasChildren(i) ) then
			call LogMessage(log, DEBUG_LOGGING_LEVEL, "face i = ", i )
		
			call GetParticlesNearFace(nearbyParticles, aMesh, i)

			call SetBasisVectorsAndCoefficients( self%basisVectors(:,:,aMesh%faces%centerParticle(i)), localCoeffs, &
					 nearbyParticles, aMesh, dataField)

			if ( dataField%nDim == 1) then
				self%coeffs1(:,aMesh%faces%centerParticle(i)) = localCoeffs(:,1)
			elseif ( dataField%nDim == 2) then
				self%coeffs1(:,aMesh%faces%centerParticle(i)) = localCoeffs(:,1)
				self%coeffs2(:,aMesh%faces%centerParticle(i)) = localCoeffs(:,2)
			elseif ( dataField%nDim == 3) then
				self%coeffs1(:,aMesh%faces%centerParticle(i)) = localCoeffs(:,1)
				self%coeffs2(:,aMesh%faces%centerParticle(i)) = localCoeffs(:,2)
				self%coeffs3(:,aMesh%faces%centerParticle(i)) = localCoeffs(:,3)
			endif		
		endif
	enddo
end subroutine

function InterpolateScalar( self, aMesh, interpLoc )
	real(kreal) :: InterpolateScalar
	type(MLSQ), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	real(kreal), intent(in) :: interpLoc(3)
	!
	integer(kint) :: faceIndex, particleIndex
	real(kreal) :: qt(3,3), locCoords(3), locBasis(10), projMat(3,3), pLoc(3)
	real(kreal) :: p, q
	
	faceIndex = LocateFaceContainingPoint(aMesh, interpLoc)
	
	particleIndex = aMesh%faces%centerParticle(faceindex) ! find nearest particle? or just nearest face?
	pLoc = PhysCoord(amesh%particles, particleIndex)
	projMat = SphereProjection( pLoc )
	
	pLoc = MATMUL( projMat, interpLoc)
	
	qt = Transpose(self%basisVectors(:,:,particleIndex))
	locCoords = MATMUL( qt, pLoc - PhysCoord(aMesh%particles, particleIndex))
	
	p = locCoords(1)
	q = locCoords(2)
	locBasis = [ 1.0_kreal, p, p*p, p*p*p, q, p*q, p*p*q, q*q, p*q*q, q*q*q]
	InterpolateScalar = sum( self%coeffs1(:,particleIndex) * locBasis)
end function

function InterpolateVector( self, aMesh, interpLoc )
	real(kreal) :: InterpolateVector(3)
	type(MLSQ), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	real(kreal), intent(in) :: interpLoc(3)
	!
	integer(kint) :: faceIndex, particleIndex
	real(kreal) :: qt(3,3), locCoords(3), locBasis(10), projMat(3,3), pLoc(3)
	real(kreal) :: p, q
	
	faceIndex = LocateFaceContainingPoint(aMesh, interpLoc)
	
	particleIndex = aMesh%faces%centerParticle(faceIndex) ! find nearest particle? or just nearest face?
	pLoc = PhysCoord(amesh%particles, particleIndex)
	projMat = SphereProjection( pLoc )
	
	pLoc = MATMUL( projMat, interpLoc)
	qt = Transpose(self%basisVectors(:,:,particleIndex))
	locCoords = MATMUL( qt, pLoc - PhysCoord(aMesh%particles, particleIndex))
	
	p = locCoords(1)
	q = locCoords(2)
	locBasis = [ 1.0_kreal, p, p*p, p*p*p, q, p*q, p*p*q, q*q, p*q*q, q*q*q]
	
	InterpolateVector(1) = sum( self%coeffs1(:,particleIndex) * locBasis)
	InterpolateVector(2) = sum( self%coeffs2(:,particleIndex) * locBasis)
	if ( associated(self%coeffs3) ) then
		InterpolateVector(3) = sum( self%coeffs3(:,particleIndex) * locBasis)	
	else
		InterpolateVector(3) = 0.0_kreal
	endif
end function

function InterpolateScalarGradient( self, aMesh, interpLoc )
	real(kreal), dimension(3) :: InterpolateScalarGradient
	type(MLSQ), intent(in) :: self
	type(PolyMesh2d), intent(in) :: amesh
	real(kreal), intent(in), dimension(3) :: interpLoc
	!
	integer(kint) :: faceIndex, particleIndex
	real(kreal) :: qt(3,3), locCoords(3), locBasisDp(10), locBasisDq(10), projMat(3,3), pLoc(3)
	real(kreal) :: p, q, dsdp, dsdq
	
	faceIndex = LocateFaceContainingPoint(aMesh, interpLoc)
	
	particleIndex = aMesh%faces%centerParticle(faceIndex)
	
	pLoc = PhysCoord(amesh%particles, particleIndex)
	projMat = SphereProjection( pLoc )
	
	pLoc = MATMUL( projMat, interpLoc)
	qt = Transpose(self%basisVectors(:,:,particleIndex))
	locCoords = MATMUL( qt, pLoc - PhysCoord(aMesh%particles, particleIndex))
	
	p = locCoords(1)
	q = locCoords(2)
	locBasisDp = [0.0_kreal, &
				  1.0_kreal, &
				  2.0_kreal * p, &
				  3.0_kreal * p * p, &
				  0.0_kreal, &
				  q, &
				  2.0_kreal * p * q, &
				  0.0_kreal, &
				  q * q, &
				  0.0_kreal ]
	locBasisDq = [0.0_kreal, &
				  0.0_kreal, &
				  0.0_kreal, &
				  0.0_kreal, &
				  1.0_kreal, &
				  p, &
				  p*p, &
				  2.0_kreal * q, &
				  2.0_kreal * p * q, &
				  3.0_kreal * q * q ] 
	dsdp = sum( self%coeffs1(:, particleIndex) * locBasisDp )
	dsdq = sum( self%coeffs1(:, particleIndex) * locBasisDq )
	InterpolateScalarGradient(1) = dsdp * self%basisVectors(1,1,particleIndex) + dsdq * self%basisVectors(1,2,particleIndex)
	InterpolateScalarGradient(2) = dsdp * self%basisVectors(2,1,particleIndex) + dsdq * self%basisVectors(2,2,particleIndex)
	InterpolateScalarGradient(3) = dsdp * self%basisVectors(3,1,particleIndex) + dsdq * self%basisVectors(3,2,particleIndex)
end function

function InterpolateScalarLaplacian(self, aMesh, interpLoc)
	real(kreal) :: InterpolateScalarLaplacian
	type(MLSQ), intent(in) :: self
	type(PolyMesh2d), intent(in) :: amesh
	real(kreal), intent(in), dimension(3) :: interpLoc
	!
	integer(kint) :: faceIndex, particleIndex
	real(kreal) :: qt(3,3), locCoords(3), locBasis(10), projMat(3,3), pLoc(3)
	real(kreal) :: p, q
	
	faceIndex = LocateFaceContainingPoint(aMesh, interpLoc)
	
	particleIndex = aMesh%faces%centerParticle(faceIndex)
	
	pLoc = PhysCoord(amesh%particles, particleIndex)
	projMat = SphereProjection( pLoc )
	
	pLoc = MATMUL( projMat, interpLoc)
	qt = Transpose(self%basisVectors(:,:,particleIndex))
	locCoords = MATMUL( qt, pLoc - PhysCoord(aMesh%particles, particleIndex))
	
	p = locCoords(1)
	q = locCoords(2)
	
	locBasis = [ 0.0_kreal, &
				 0.0_kreal, &
				 2.0_kreal, &
				 6.0_kreal * p, &
				 0.0_kreal, &
				 0.0_kreal, &
				 2.0_kreal * q, & 
				 2.0_kreal, &
				 2.0_kreal * p, &
				 6.0_kreal * q ]
	InterpolateScalarLaplacian = sum( locBasis * self%coeffs1(:,particleIndex))
end function

subroutine MLSQGradientAtParticles( self, aMesh, scalarGrad )
	type(MLSQ), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(inout) :: scalarGrad
	!
	integer(kint) :: i
	
	call SetFieldToZero(scalarGrad)
	scalarGrad%N = aMesh%particles%N
	
	do i = 1, aMesh%particles%N
		scalarGrad%xComp(i) = self%coeffs1(2,i) * self%basisVectors(1,1,i) + self%coeffs1(5,i) * self%basisVectors(1,2,i)
		scalarGrad%yComp(i) = self%coeffs1(2,i) * self%basisVectors(2,1,i) + self%coeffs1(5,i) * self%basisVectors(2,2,i)
		scalarGrad%zComp(i) = self%coeffs1(2,i) * self%basisVectors(3,1,i) + self%coeffs1(5,i) * self%basisVectors(3,2,i)
	enddo
end subroutine

subroutine MLSQLaplacianAtParticles( self, aMesh, scalarLap )
	type(MLSQ), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(inout) :: scalarLap
	!
	integer(kint) :: i
	
	call SetFieldToZero(scalarLap)
	scalarLap%N = aMesh%particles%N
	do i = 1, aMesh%particles%N
		scalarLap%scalar(i) = 2.0_kreal * self%coeffs1(3,i) + 2.0_kreal * self%coeffs1(8,i)
	enddo
end subroutine

subroutine MLSQDoubleDotProductAtParticles( self, aMesh, doubleDot )
	type(MLSQ), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(inout) :: doubleDot
	!
	integer(kint) :: i
	real(kreal) :: ux, uy, uz, vx, vy, vz, wx, wy, wz 
	
	call SetFieldToZero(doubleDot)
	doubleDot%N = aMesh%particles%N
	do i = 1, aMesh%particles%N
		ux = self%coeffs1(2,i) * self%basisVectors(1,1,i) + self%coeffs1(5,i) * self%basisVectors(1,2,i) 
		uy = self%coeffs1(2,i) * self%basisVectors(2,1,i) + self%coeffs1(5,i) * self%basisVectors(2,2,i)
		uz = self%coeffs1(2,i) * self%basisVectors(3,1,i) + self%coeffs1(5,i) * self%basisVectors(3,2,i)
		vx = self%coeffs2(2,i) * self%basisVectors(1,1,i) + self%coeffs2(5,i) * self%basisVectors(1,2,i)
		vy = self%coeffs2(2,i) * self%basisVectors(2,1,i) + self%coeffs2(5,i) * self%basisVectors(2,2,i)
		vz = self%coeffs2(2,i) * self%basisVectors(3,1,i) + self%coeffs2(5,i) * self%basisVectors(3,2,i)
		wx = self%coeffs3(2,i) * self%basisVectors(1,1,i) + self%coeffs3(5,i) * self%basisVectors(1,2,i)
		wy = self%coeffs3(2,i) * self%basisVectors(2,1,i) + self%coeffs3(5,i) * self%basisVectors(2,2,i)
		wz = self%coeffs3(2,i) * self%basisVectors(3,1,i) + self%coeffs3(5,i) * self%basisVectors(3,2,i)
		doubleDot%scalar(i) = ux*ux + vy*vy + wz*wz + 2.0_kreal * ( uy*vx + uz * wx + vz * wy )
	enddo
end subroutine

!
!----------------
! private methods
!----------------
!
subroutine SetBasisVectorsAndCoefficients( basisVectors, coeffsOut, nearbyParticles, aMesh, dataField)
	real(kreal), intent(out) :: basisVectors(3,3)
	real(kreal), intent(out) :: coeffsOut(10,3)
	type(STDIntVector), intent(in) :: nearbyParticles
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: dataField
	!
	integer(kint) :: i	
	real(kreal) :: projMatrix(3,3), qt(3,3)
	real(kreal) :: localVectors(3,21)
	real(kreal) :: lsqA(21,10)
	real(kreal) :: scalarB(21), vec2B(21,2), vec3B(21,3)
	real(kreal) :: x0(3), xi(3)
	real(kreal) :: p, q
	real(kreal) :: work(20)
	integer(kint) :: lpkInfo
	
	x0 = PhysCoord(aMesh%particles, nearbyParticles%int(1) )

	!
	!	projection matrix to plane tangent to sphere at x0
	!			
	projMatrix(1,1) = 1.0_kreal - x0(1)*x0(1)
	projMatrix(2,1) = - x0(1) * x0(2)
	projMatrix(3,1) = - x0(1) * x0(3)
	projMatrix(1,2) = - x0(1) * x0(2)
	projMatrix(2,2) = 1.0_kreal - x0(2) * x0(2)
	projMatrix(3,2) = - x0(2) * x0(3)
	projMatrix(1,3) = - x0(1) * x0(3)
	projMatrix(2,3) = - x0(2) * x0(3)
	projMatrix(3,3) = 1.0_kreal - x0(3) * x0(3)
	
	localVectors(:,1) = 0.0_kreal
	do i = 2, nearbyParticles%N
		xi = PhysCoord(aMesh%particles, nearbyParticles%int(i) )
		localVectors(:,i) = MATMUL( projMatrix, xi) ! project nearby particles into tangent plane
		localVectors(:,i) = localVectors(:,i) - x0	! set origin to x0
	enddo
	
	!
	!	construct orthonormal basis for tangent plane
	!
	basisVectors(:,1) = localVectors(:,2) / sqrt(sum( localVectors(:,2)*localVectors(:,2)))
	basisVectors(:,2) = localVectors(:,3) - sum( localVectors(:,3) * basisVectors(:,1)) * basisVectors(:,1)
	basisVectors(:,2) = basisVectors(:,2) / sqrt(sum( basisVectors(:,2) * basisVectors(:,2)))
	basisVectors(:,3) = x0 / sqrt(sum(x0*x0))
	
	!
	!	find coordinates in new basis
	!
	qt = Transpose(basisVectors)
	do i = 2, nearbyParticles%N
		localVectors(:,i) = MATMUL(qt, localVectors(:,i))
	enddo
	
	!
	!	build lsq matrix A
	!
	lsqA(1:nearbyParticles%N,1) = 1.0_kreal
	do i = 1, nearbyParticles%N
		p = localVectors(1,i)
		q = localVectors(2,i)
		lsqA(i,1) = 1.0_kreal
		lsqA(i,2) = p
		lsqA(i,3) = p*p
		lsqA(i,4) = p*p*p
		lsqA(i,5) = q
		lsqA(i,6) = p*q
		lsqA(i,7) = p*p*q
		lsqA(i,8) = q*q
		lsqA(i,9) = p*q*q
		lsqA(i,10) = q*q*q
	enddo
	
	!! DEBUG !!
	if ( procRank == 0 ) then
		write(6,'(A)',advance='NO') "A = "
		do i = 1, nearbyParticles%N
			write(6,*) lsqA(i,:)
		enddo
	endif		
	
	!
	!	build rhs B
	!
	if ( dataField%nDim == 1) then
		do i = 1, nearbyParticles%N
			scalarB(i) = dataField%scalar(nearbyParticles%int(i))
		enddo
	elseif ( dataField%nDim == 2) then
		do i = 1, nearbyParticles%N
			vec2B(i,1) = dataField%xComp(nearbyParticles%int(i))
		enddo
		do i = 1, nearbyParticles%N
			vec2B(i,2) = dataField%yComp(nearbyParticles%int(i))
		enddo
	elseif ( dataField%nDim == 3) then
		do i = 1, nearbyParticles%N
			vec3B(i,1) = dataField%xComp(nearbyParticles%int(i))
		enddo
		do i = 1, nearbyParticles%N
			vec3B(i,2) = dataField%yComp(nearbyParticles%int(i))
		enddo
		do i = 1, nearbyParticles%N
			vec3B(i,3) = dataField%zComp(nearbyParticles%int(i))
		enddo
	endif
	
	!
	!	solve lsq system with LAPACK QR factorization
	!
	if ( dataField%nDim == 1 ) then
		call DGELS('N', nearbyParticles%N, 10, 1, lsqA, nearbyParticles%N, scalarB, nearbyParticles%N,&
			 work, 20, lpkInfo)
		coeffsOut(:,1) = scalarB(1:10)
	elseif ( dataField%nDim == 2) then
		call DGELS('N', nearbyParticles%N, 10, 2, lsqA, nearbyParticles%N, vec2B, nearbyParticles%N,&
			 work, 20, lpkInfo)
		coeffsOut(:,1) = vec2B(1:10,1)
		coeffsOut(:,2) = vec2B(1:10,2)
	elseif (dataField%nDIm == 3) then
		call DGELS('N', nearbyParticles%N, 10, 3, lsqA, nearbyParticles%N, vec3B, nearbyParticles%N,&
			 work, 20, lpkInfo)
		coeffsOut(:,1) = vec3B(1:10,1)
		coeffsOut(:,2) = vec3B(1:10,2)
		coeffsOut(:,3) = vec3B(1:10,3)
	endif

end subroutine


subroutine GetParticlesNearVertex( nearbyParticles, aMesh, particleIndex )
	type(STDIntVector), intent(out) :: nearbyParticles
	type(PolyMesh2d), intent(in) :: aMesh
	integer(kint), intent(in) :: particleIndex
	!
	integer(kint) :: i, j, k
	type(STDIntVector) :: adjFaces, faceVerts, adjFaces2
	logical(klog) :: isNew
	
	call initialize(nearbyParticles)
	call nearbyParticles%pushBack( particleIndex)
	
	call CCWFacesAroundVertex( aMesh, adjFaces, particleIndex )
	do i = 1, adjFaces%N
		call nearbyParticles%pushBack( aMesh%faces%centerParticle(adjFaces%int(i)))
	enddo
	do i = 1, adjFaces%N
		call CCWVerticesAroundFace(amesh, faceVerts, adjFaces%int(i))
		do j = 1, faceVerts%N
			isNew = .TRUE.
			do k = 1, nearbyParticles%N
				if ( faceVerts%int(j) == nearbyParticles%int(k) ) isNew = .FALSE.
			enddo
			if ( isNew ) call nearbyParticles%pushBack( faceVerts%int(j))
		enddo
		call CCWAdjacentFaces( aMesh, adjFaces2, adjFaces%int(i) )
		do j = 1, adjFaces2%N
			isNew = .TRUE.
			do k = 1, nearbyParticles%N
				if ( adjFaces2%int(j) == nearbyParticles%int(k) ) isNew = .FALSE.
			enddo
			if ( isNew ) call nearbyParticles%pushBackUnique( aMesh%faces%centerParticle(adjFaces2%int(i)) )
		enddo
	enddo
end subroutine

subroutine GetParticlesNearFace( nearbyParticles, aMesh, faceIndex )
	type(STDIntVector), intent(out) :: nearbyParticles
	type(PolyMesh2d), intent(in) :: aMesh
	integer(kint), intent(in) :: faceIndex
	!
	integer(kint) :: i, j, k 
	type(STDIntVector) :: adjFaces, faceEdges, faceVerts
	logical(klog) :: isNew
	
	call initialize(nearbyParticles)
	call nearbyParticles%pushBack( aMesh%faces%centerParticle(faceIndex) )
	
	call CCWVerticesAroundFace(aMesh, faceVerts, faceIndex)
	do i = 1, faceVerts%N
		call nearbyParticles%pushBack( faceVerts%int(i))
	enddo
	
	call CCWAdjacentFaces( aMesh, adjFaces, faceIndex)
	do i = 1, adjFaces%N
		call nearbyParticles%pushBack( aMesh%faces%centerParticle( adjFaces%int(i) ))
	enddo
	
	do i = 1, adjFaces%N
		call CCWVerticesAroundFace(aMesh, faceVerts, adjFaces%int(i))
		do j = 1, faceVerts%N
			isNew = .TRUE.
			do k = 1, nearbyParticles%N
				if ( faceVerts%int(j) == nearbyParticles%int(k) ) isNew = .FALSE.
			enddo
			if ( isNew ) call nearbyParticles%pushBackUnique(faceVerts%int(j))
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
