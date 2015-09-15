module TriScatteredInterpModule
!> @file ScatteredTriInterp.f90
!> Implements the BIVAR algorithm (H. Akima, 1978, _ACM TOMS_; H. Akima, 1984, _Rocky Mtn. J. Math._)
!> on an existing LPM triangular mesh.
!> 
!> @author Peter Bosler, Sandia National Laboratories Center for Computing Research
!>
!>
!> @defgroup TriScatteredInterp TriScatteredInterp module
!> @brief Provides data structures and methods for computing bivariate interpolation of scattered data on a 2d mesh
!> of triangular faces. On each face a bivariate quintic Hermite polynomial is constructed to provide a C1 interpolating surface.
!> This method is due to H. Akima, 1974 and H. Akima, 1984.
!> @{
use NumberKindsModule
use LoggerModule
use SphereGeomModule
use PlaneGeomModule
use ParticlesModule
use EdgesModule
use FacesModule
use PolyMesh2dModule
use FieldModule
use STDIntVectorModule

implicit none

private
public TriScatteredInterp, New, Delete, Copy
public generatefaceUVs, estimatePartialDerivatives!, BIVAREstimatePartials
public SetPlanarCoefficients
public InterpolateScalar, InterpolateGradient, InterpolateLaplacian

!> @class TriScatteredInterp
!> @brief Container class for holding the coefficients of local quintic polynomials to each face, and the loca u-v map
!> used to create that polynomial.
type TriScatteredInterp
	real(kreal), pointer :: faceU(:,:) => null()
	real(kreal), pointer :: faceV(:,:) => null()
	real(kreal), pointer :: faceCoeffs(:,:) => null()
	logical(klog) :: mapsReady = .FALSE.
	logical(klog) :: coeffsReady = .FALSE.
	
	contains
		final :: deletePrivate
end type

interface New
	module procedure newPrivate
!	module procedure newFromCopy
end interface

interface Copy
	module procedure copyPrivate
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
character(len=28), save :: logKey = 'TriScatteredInterpLog'
integer(kint), parameter :: logLevel = TRACE_LOGGING_LEVEL
character(len=MAX_STRING_LENGTH) :: logString

contains

subroutine newPrivate( self, aMesh )
	type(TriScatteredInterp), intent(out) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	
	if (.NOT. logInit) call InitLogger(log, procRank)
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL,trim(logKey)//"TriScatteredInterp::New"," entering.")
	
	if ( aMesh%geomKind == PLANAR_GEOM ) then
		allocate(self%faceU(2,aMesh%faces%N))
		self%faceU = 0.0_kreal
		allocate(self%faceV(2,aMesh%faces%N))
		self%faceV = 0.0_kreal
	else
		allocate(self%faceU(3,aMesh%faces%N))
		self%faceU = 0.0_kreal
		allocate(self%faceV(3,aMesh%faces%N))
		self%faceV = 0.0_kreal
	endif
	
	allocate(self%faceCoeffs(21,aMesh%faces%N))
	self%faceCoeffs = 0.0_kreal
	
	call generatefaceUVs(self, aMesh)
end subroutine

subroutine copyPrivate(self, other )
	type(TriScatteredInterp), intent(out) :: self
	type(TriScatteredInterp), intent(in) :: other
	!
	integer(kint) :: i
	call LogMessage(log, DEBUG_LOGGING_LEVEL,trim(logKey)//"TriScatteredInterp::Copy"," entering.")
	if ( size(self%faceU,1) /= size(other%faceU,1) .OR. size(self%faceU,2) /= size(other%faceU,2) ) then
		call LogMessage( log, ERROR_LOGGING_LEVEL, trim(logkey)//" newFromCopy ERROR : ", "size mismatch.")
		return
	endif
	
	do i = 1, size(other%faceU,2)
		self%faceU(:,i) = other%faceU(:,i)
		self%faceV(:,i) = other%faceV(:,i)
		self%faceCoeffs(:,i) = other%faceCoeffs(:,i)
	enddo
end subroutine

subroutine deletePrivate(self)
	type(TriScatteredInterp), intent(inout) :: self

	if ( associated(self%faceU)) deallocate(self%faceU)
	if (associated(self%faceV)) deallocate(self%faceV)
	if ( associated(self%faceCoeffs)) deallocate(self%faceCoeffs)
end subroutine

subroutine generatefaceUVs( self, aMesh )
	type(TriScatteredInterp), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	!
	real(kreal) :: v0(3), v1(3), v2(3)
	integer(kint) :: i
	call LogMessage(log, DEBUG_LOGGING_LEVEL,trim(logKey)//"TriScatteredInterp::generateFaceUVs"," entering.")
	if ( aMesh%geomKind == PLANAR_GEOM ) then
		do i = 1, aMesh%Faces%N
			if (.NOT. aMesh%faces%hasChildren(i) ) then
				v0 = PhysCoord(aMesh%particles, aMesh%faces%vertices(1,i))
				v1 = PhysCoord(aMesh%particles, aMesh%faces%vertices(2,i))
				v2 = PhysCoord(aMesh%particles, aMesh%faces%vertices(3,i))

				self%faceU(:,i) = v1(1:2) - v0(1:2)
				self%faceV(:,i) = v2(1:2) - v0(1:2)
			endif
		enddo
	else
		do i = 1, aMesh%Faces%N
			if (.NOT. aMesh%faces%hasChildren(i) ) then
				v0 = PhysCoord(aMesh%particles, aMesh%faces%vertices(1,i))
				v1 = PhysCoord(aMesh%particles, aMesh%faces%vertices(2,i))
				v2 = PhysCoord(aMesh%particles, aMesh%faces%vertices(3,i))

				self%faceU(:,i) = v1 - v0
				self%faceV(:,i) = v2 - v0
			endif
		enddo	
	endif
	self%mapsReady = .TRUE.
end subroutine

!> @brief Sets the coefficients of a quintic Hermite polynomial on each triangular face of a 
!> PolyMesh2d object in planar geometry.
!> 
!> Assumptions : These coefficients will interpolate a surface z(x,y) defined on a planar domain.
!> Scalar values z_i = z(x_i,y_i) are stored in scalarField and correspond to particle i. @n
!> Gradient vectors [ dz/dx(x_i,y_i), dz/dy(x_i,y_i)] are stored in gradField. @n
!> Second derivatives [ d^2 z/dx^2(x_i,y_i), d^2 z/dxdy(x_i,y_i), d^2 z/dy^2(x_i,y_i)] are stored in secondDerivs.
!> 
!> Coefficients are ordered at each face to correspond with the double summation notation 
!> @f$ p(x,y) = \sum_{j=0}^5\sum_{k=0}^{5-j} p_{jk}u^jv^k @f$
!> where _u_ and _v_ are the local coordinates of a point in a triangular face. 
!>
!> @param self
!> @param aMesh planar Mesh object with triangular faces
!> @param scalarField scalar data to be interpolated
!> @param gradient of scalarField
!> @param secondDerivs second derivatives of scalarField
subroutine SetPlanarCoefficients(self, aMesh, scalarField, gradField, secondDerivs )
	type(TriScatteredInterp), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: scalarField
	type(Field), intent(in) :: gradField
	type(Field), intent(in) :: secondDerivs
	!
	integer(kint) :: i, j
	real(kreal) :: a, b, c, d, zu(3), zv(3), zuu(3), zuv(3), zvv(3)
	real(kreal) :: p00, p01, p02, p03, p04, p05, p10, p11, p12, p13, p14, p20, p21, p22, p23, p30, p31, p32, &
				   p40, p41, p50
	real(kreal) :: h1, h2, h3, lu, lv, thuv, thus, aa, bb, cc, dd, g1, g2
	real(kreal) :: testScalar, center(3), v0(3), u, v
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL,trim(logKey)//"TriScatteredInterp::SetPlanarCoefficients"," entering.")
	
	if ( aMesh%geomKind /= PLANAR_GEOM ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" SetPlanarCoefficients ERROR : ", " invalid geomKind.")
		return
	endif
	if ( gradField%nDim /= 2 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" SetPlanarCoefficients ERROR : ", " gradField nDim.")
		return
	endif
	if ( secondDerivs%nDim /= 3 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" SetPlanarCoefficients ERROR : ", " secondDerivs nDim.")
		return
	endif
	
	if ( .NOT. self%mapsReady ) then
		!call LogMessage(log, WARNING_LOGGING_LEVEL, trim(logKey)//" SetPlanarCoefficients WARNING: ", " face maps not ready.")
		call generatefaceUVs(self, aMesh)
	endif
	
	do i = 1, aMesh%faces%N
		if ( .NOT. aMesh%faces%hasChildren(i) ) then
			a = self%faceU(1,i)
			b = self%faceV(1,i)
			c = self%faceU(2,i)
			d = self%faceV(2,i)
			
			!
			! change partial derivatives to local basis
			!
			do j = 1, 3
				zu(j) = a * gradField%xComp( aMesh%faces%vertices(j,i) ) + c * gradField%yComp( aMesh%faces%vertices(j,i) )
				zv(j) = b * gradField%xComp( aMesh%faces%vertices(j,i) ) + d * gradField%yComp( aMesh%faces%vertices(j,i) )
				zuu(j)= a * a * secondDerivs%xComp( aMesh%faces%vertices(j,i) ) + &
					    2.0_kreal * a * c * secondDerivs%yComp( aMesh%faces%vertices(j,i) ) + &
					    c * c * secondDerivs%zComp( aMesh%faces%vertices(j,i) )
				zuv(j) = a * b * secondDerivs%xComp( aMesh%faces%vertices(j,i) ) + &
						 (a*d + b*c) * secondDerivs%yComp( aMesh%faces%vertices(j,i) ) + &
						 c * d * secondDerivs%zComp( aMesh%faces%vertices(j,i) )
				zvv(j) = b * b * secondDerivs%xComp( aMesh%faces%vertices(j,i) ) + &
				         2.0_kreal * b * d * secondDerivs%yComp( aMesh%faces%vertices(j,i) ) + &
				         d * d * secondDerivs%zComp( aMesh%faces%vertices(j,i) )
			enddo
			
			p00 = scalarField%scalar( aMesh%faces%vertices(1,i) )
			p10 = zu(1)
			p01 = zv(1)
			p20 = 0.5_kreal * zuu(1)
			p11 = zuv(1)
			p02 = 0.5_kreal * zvv(1)
			p30 = 10.0_kreal * scalarField%scalar(aMesh%faces%vertices(2,i)) - 4.0_kreal*zu(2) + 0.5_kreal * zuu(2) - &
				  10.0_kreal * p00 - 6.0_kreal * p10 - 3.0_kreal * p20
			p40 = -15.0_kreal * scalarField%scalar(aMesh%faces%vertices(2,i)) + 7.0_kreal * zu(2) - zuu(2) + &
				  15.0_kreal * p00 + 8.0_kreal * p10 + 3.0_kreal * p20
			p50 = 6.0_kreal * scalarField%scalar(aMesh%faces%vertices(2,i)) - 3.0_kreal * zu(2) + 0.5_kreal * zuu(2) - &
				  6.0_kreal * p00 - 3.0_kreal * p10 - p20
			p03 = 10.0_kreal * scalarField%scalar(aMesh%faces%vertices(3,i)) - 4.0_kreal * zv(3) + 0.5_kreal * zvv(3) - &
				  10.0_kreal * p00 - 6.0_kreal * p01 - 3.0_kreal * p02
			p04 = -15.0_kreal * scalarField%scalar(aMesh%faces%vertices(3,i)) + 7.0_kreal * zv(3) - zvv(3) + &
				  15.0_kreal * p00 + 8.0_kreal * p01 + 3.0_kreal * p02
			p05 = 6.0_kreal * scalarField%scalar(aMesh%faces%vertices(3,i)) - 3.0_kreal * zv(3) + 0.5_kreal * zvv(3) - &
				  6.0_kreal * p00 - 3.0_kreal * p01 - p02
			lu = a*a + c*c
			lv = b*b + d*d
			thuv = atan2(d,b) - atan2(c,a)
			p41 = 5.0_kreal * lv * cos(thuv) * p50 / lu
			p14 = 5.0_kreal * lu * cos(thuv) * p05 / lv
			p21 = 3.0_kreal * zv(2) - zuv(2) - 3.0_kreal * p01 - 2.0_kreal * p11 + p41
			p31 = -2.0_kreal * zv(2) + zuv(2) + 2.0_kreal * p01 + p11 - 2.0_kreal * p41
			p12 = 3.0_kreal * zu(3) - zuv(3) - 3.0_kreal * p10 - 2.0_kreal * p11 + p14
			p13 = -2.0_kreal * zu(3) + zuv(3) + 2.0_kreal * p10 + p11 - 2.0_kreal * p14
			thus = atan2( d-c, b-a) - atan2( c, a)
			aa = sin( thuv - thus) / ( lu * sin(thuv))
			bb = - cos( thuv - thus) / (lu* sin(thuv))
			cc = sin(thus) / ( lv * sin(thuv) )
			dd = cos(thus) / ( lv * sin(thuv) )
			g1 = aa * aa * cc * ( 3.0_kreal * bb * cc - 2.0_kreal * aa * dd)
			g2 = aa * cc * cc * ( 2.0_kreal * bb * cc + 3.0_kreal * aa * dd)
			h1 = -5.0_kreal * aa**4 * bb * p50 - aa**3 * (4.0_kreal * bb * cc + aa * dd) * p41 - &
				 cc**3 * (bb * cc + 4.0_kreal * aa * dd) * p14 - 5.0_kreal * cc**4 * dd * p05
			h2 = 0.5_kreal * zvv(2) - p02 - p12
			h3 = 0.5_kreal * zuu(3) - p20 - p21
			p22 = ( g1 * h2 + g2 * h3 - h1)/(g1+g2)
			p32 = h2 - p22
			p23 = h3 - p22
			
			self%faceCoeffs(:,i) = [ p00, p01, p02, p03, p04, p05, &
			 						 p10, p11, p12, p13, p14, &
			 						 p20, p21, p22, p23, &
			 						 p30, p31, p32, &
			 						 p40, p41, &
			 						 p50 ]
			 						 
			if ( logLevel == DEBUG_LOGGING_LEVEL ) then
				center = FaceCenterPhysCoord(aMesh%faces, i, aMesh%particles )
				v0 = PhysCoord(aMesh%particles, aMesh%faces%vertices(1,i))
				u = d * ( center(1) - v0(1)) - b * ( center(2)-v0(2))
				u = u / (a*d - b*c)
				v = -c *( center(1) - v0(1)) + a * ( center(2)-v0(2))
				v = v / (a*d - b*c)
				testScalar = p00 + p01 * v + p02 * v**2 + p03 * v**3 + p04 * v**4 + p05 * v**5 + &
							 p10 * u + p11 * u * v + p12 * u * v**2 + p13 * u * v**3 + p14 * u * v**4 + &
							 p20 * u**2 + p21 * u**2 * v + p22 * u**2 * v**2 + p23 * u**2 * v**3 + &
							 p30 * u**3 + p31 * u**3 * v + p32 * u**3 * v**2 + &
							 p40 * u**4 + p41 * u**4 * v + &
							 p50 * u**5
				
				if ( abs(testScalar - scalarField%scalar(aMesh%faces%centerParticle(i))) > 1.0E-6 ) then
					call StartSection(log)
					call LogMessage(log, WARNING_LOGGING_LEVEL, "DEBUG ", " Coefficients Interpolation ERROR ")
					call LogMessage(log, WARNING_LOGGING_LEVEL, "at face ", i)
					call LogMessage(log, WARNING_LOGGING_LEVEL, "testScalar = ", testScalar )
					call LogMessage(log, WARNING_LOGGING_LEVEL, "carried scalar = ", &
									scalarField%scalar(aMesh%faces%centerParticle(i)))
					call EndSection(log)
				else
					call LogMessage(log, DEBUG_LOGGING_LEVEL,trim(logkey)//" SetPlanarCoefficients :", " interp test passed.")
				endif
			endif
			
			
		endif
	enddo
	
	self%coeffsReady = .TRUE.
end subroutine

subroutine estimatePartialDerivatives(self, aMesh, scalarField, scalarGrad, secondPartials )
	type(TriScatteredInterp), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: scalarField
	type(Field), intent(inout) :: scalarGrad
	type(Field), intent(inout) :: secondPartials
	!
	real(kreal) :: vectorSum(3), vectorSumXX(3), vectorSumYY(3)
	integer(kint) :: i, j, piIndex, pjIndex
	real(kreal) :: pi(3), pj(3), p0(3), p0x(3), p0y(3), pix(3), piy(3), pjx(3), pjy(3)
	
	call SetFieldToZero(scalarGrad)
	call SetFieldToZero(secondPartials)
	
	call LogMessage(log, DEBUG_LOGGING_LEVEL,trim(logKey)//"TriScatteredInterp::estimatePartialDerivatives"," entering.")
	
	scalarGrad%N = aMesh%particles%N
	secondPartials%N = aMesh%particles%N
	
	if ( aMesh%geomKind == PLANAR_GEOM) then
		!
		!	estimate gradient vector at vertex particles
		!		
		do i = 1, amesh%particles%N
			if ( aMesh%particles%isPassive(i) ) then
				
				piIndex = 0
				pjIndex = 0
				vectorSum = 0.0_kreal	
				p0 = PhysCoord(aMesh%particles, i)
				p0(3) = scalarField%scalar(i)
			
				do j = 1, aMesh%particles%nEdges(i)
					if ( aMesh%edges%orig( aMesh%particles%incidentEdges(j,i) ) == i ) then
						piIndex = aMesh%edges%dest( aMesh%particles%incidentEdges(j,i))
					elseif ( aMesh%edges%dest( aMesh%particles%incidentEdges(j,i)) == i ) then
						piIndex = aMesh%edges%orig( aMesh%particles%incidentEdges(j,i))
					else
						call LogMessage(log,ERROR_LOGGING_LEVEL,trim(logKey)//&
										" estimatePartialDerivatives ERROR : bad edge connectivity at particle ", i)
					endif
				
					pi = PhysCoord(aMesh%particles, piIndex)
					pi(3) = scalarField%scalar(piIndex)
					pi = pi - p0
					
					if ( aMesh%edges%orig( aMesh%particles%incidentEdges( mod(j, aMesh%particles%nEdges(i)) + 1, i) ) == i ) then
						pjIndex = aMesh%edges%dest(aMesh%particles%incidentEdges(mod(j,aMesh%particles%nEdges(i)) + 1, i))
					elseif ( aMesh%edges%dest(aMesh%particles%incidentEdges(mod(j,aMesh%particles%nEdges(i)) +1, i) ) == i ) then
						pjIndex = aMesh%edges%orig(aMesh%particles%incidentEdges(mod(j,aMesh%particles%nEdges(i)) + 1, i))
					else
						call LogMessage(log,ERROR_LOGGING_LEVEL,trim(logKey)//&
										" estimatePartialDerivatives ERROR : bad edge connectivity at particle ", i)
					endif
				
					pj = PhysCoord(aMesh%particles, pjIndex)
					pj(3) = scalarField%scalar(pjIndex)
					pj = pj - p0
							
					vectorSum(1) = vectorSum(1) + pi(2)*pj(3) - pi(3)*pj(2)
					vectorSum(2) = vectorSum(2) + pi(3)*pj(1) - pi(1)*pj(3)
					vectorSum(3) = vectorSum(3) + pi(1)*pj(2) - pi(2)*pj(1)
				enddo	
				
				scalarGrad%xComp(i) = - vectorSum(1) / vectorSum(3)
				scalarGrad%yComp(i) = - vectorSum(2) / vectorSum(3)
			else
				if ( .NOT. aMesh%particles%isActive(i) ) then
					call LogMessage(log,ERROR_LOGGING_LEVEL,trim(logKey)//&
						" ERROR : neither isActive(i) or isPassive(i) at i = ", i)
				endif			
			endif
		enddo
		
		!
		!	estimate gradient vector at center particles
		!
		do i = 1, aMesh%faces%N
			if ( .NOT. aMesh%faces%hasChildren(i) ) then
				vectorSum = 0.0_kreal
				piIndex = 0
				pjIndex = 0
				
				p0 = FaceCenterPhysCoord(aMesh%faces, i, aMesh%particles)
				p0(3) = scalarField%scalar( aMesh%faces%centerParticle(i) )
				
				do j = 1, 3
					piIndex = aMesh%faces%vertices(j,i)
					pjIndex = aMesh%faces%vertices(mod(j,3)+1,i)
					
					pi = PhysCoord(aMesh%particles, piIndex)
					pi(3) = scalarField%scalar(piIndex)
					pi = pi - p0
					
					pj = PhysCoord(aMesh%particles, pjIndex)
					pj(3) = scalarField%scalar(pjIndex)
					pj = pj - p0
					
					vectorSum(1) = vectorSum(1) + pi(2)*pj(3) - pj(2)*pi(3)
					vectorSum(2) = vectorSum(2) + pj(1)*pi(3) - pi(1)*pj(3)
					vectorSum(3) = vectorSum(3) + pi(1)*pj(2) - pj(1)*pi(2)
				enddo
							
				scalarGrad%xComp(amesh%faces%centerParticle(i)) = - vectorSum(1)/vectorSum(3)
				scalarGrad%yComp(aMesh%faces%centerParticle(i)) = - vectorSum(2)/vectorSum(3)
			endif
		enddo
		
		!
		!	estimate second derivatives
		!
		do i = 1, aMesh%particles%N
			if ( aMesh%particles%isPassive(i) ) then
				vectorSumXX = 0.0_kreal
				vectorSumYY = 0.0_kreal
				
				p0x = PhysCoord(aMesh%particles,i)
				p0y = p0x
				p0x(3) = scalarGrad%xComp(i)
				p0y(3) = scalarGrad%yComp(i)
								
				do j = 1, aMesh%particles%nEdges(i)
					if ( aMesh%edges%orig( aMesh%particles%incidentEdges(j,i)) == i ) then
						piIndex = aMesh%edges%dest( aMesh%particles%incidentEdges(j,i))
					else
						piIndex = aMesh%edges%orig(aMesh%particles%incidentEdges(j,i))
					endif
					if ( aMesh%edges%orig( aMesh%particles%incidentEdges(mod(j,aMesh%particles%nEdges(i))+1,i)) == i ) then
						pjIndex = aMesh%edges%dest( aMesh%particles%incidentEdges(mod(j,aMesh%particles%nEdges(i))+1,i))
					else
						pjIndex = aMesh%edges%orig( aMesh%particles%incidentEdges(mod(j,aMesh%particles%nEdges(i))+1,i))
					endif
					
					pix = PhysCoord(aMesh%particles, piIndex)
					piy = pix
					pix(3) = scalarGrad%xComp(piIndex)
					piy(3) = scalarGrad%yComp(piIndex)
					
					pjx = PhysCoord(aMesh%particles, pjIndex)
					pjy = pjx
					pjx(3) = scalarGrad%xComp(pjIndex)
					pjy(3) = scalarGrad%yComp(pjIndex)	
					
					pix = pix - p0x
					piy = piy - p0y
					pjx = pjx - p0x
					pjy = pjy - p0y
					
					vectorSumXX(1) = vectorSumXX(1) + pix(2)*pjx(3) - pjx(2)*pix(3)
					vectorSumXX(2) = vectorSumXX(2) + pjx(1)*pix(3) - pix(1)*pjx(3)
					vectorSumXX(3) = vectorSumXX(3) + pix(1)*pjx(2) - pjx(1)*pix(2)
					
					vectorSumYY(1) = vectorSumYY(1) + piy(2)*pjy(3) - pjy(2)*piy(3)
					vectorSumYY(2) = vectorSumYY(2) + pjy(1)*piy(3) - piy(1)*pjy(3)
					vectorSumYY(3) = vectorSumYY(3) + piy(1)*pjy(2) - pjy(1)*piy(2)
				enddo
				
				! store d^2z/dx^2 in xComp
				secondPartials%xComp(i) = - vectorSumXX(1) / vectorSumXX(3)
				! store mixed partial in yComp
				secondPartials%yComp(i) = -0.5_kreal * (vectorSumXX(2)/vectorSumXX(3)  + vectorSumYY(1)/vectorSumYY(3) )
				! store d^z/dy^2 in zComp
				secondPartials%zComp(i) = - vectorSumYY(2) / vectorSumYY(3)

			endif
		enddo
		
		do i = 1, aMesh%faces%N
			if ( .NOT. aMesh%faces%hasChildren(i) ) then
				vectorSumXX = 0.0_kreal
				vectorSumYY = 0.0_kreal
				
				p0x = FaceCenterPhysCoord(aMesh%faces, i, aMesh%particles)
				p0y = p0x
				p0x(3) = scalarGrad%xComp( aMesh%faces%centerParticle(i))
				p0y(3) = scalarGrad%yComp( aMesh%faces%centerParticle(i))
				
				do j = 1, 3
					piIndex = aMesh%faces%vertices(j,i)
					pjIndex = aMesh%faces%vertices(mod(j,3)+1,i)
					
					pix = PhysCoord(aMesh%particles, piIndex)
					piy = pix
					pix(3) = scalarGrad%xComp(piIndex)
					piy(3) = scalarGrad%yComp(piIndex)
					
					pjx = PhysCoord(aMesh%particles, pjIndex)
					pjy = pjx
					pjx(3) = scalarGrad%xComp(pjIndex)
					pjy(3) = scalarGrad%yComp(pjIndex)
					
					vectorSumXX(1) = vectorSumXX(1) + pix(2)*pjx(3) - pjx(2)*pix(3)
					vectorSumXX(2) = vectorSumXX(2) + pjx(1)*pix(3) - pix(1)*pjx(3)
					vectorSumXX(3) = vectorSumXX(3) + pix(1)*pjx(2) - pjx(1)*pix(2)
					
					vectorSumYY(1) = vectorSumYY(1) + piy(2)*pjy(3) - pjy(2)*piy(3)
					vectorSumYY(2) = vectorSumYY(2) + pjy(1)*piy(3) - piy(1)*pjy(3)
					vectorSumYY(3) = vectorSumYY(3) + piy(1)*pjy(2) - pjy(1)*piy(2)
				enddo
				
				secondPartials%xComp(aMesh%faces%centerParticle(i)) = - vectorSumXX(1) / vectorSumXX(3)
				secondPartials%yComp(aMesh%faces%centerParticle(i)) = -0.5_kreal * &
					( vectorSumXX(2)/vectorSumXX(3) + vectorSumYY(1)/vectorSumYY(3))
				secondPartials%zComp(aMesh%faces%centerParticle(i)) = -vectorSumYY(2) / vectorSumYY(3)
			endif
		enddo
		
	else
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" estimatePartialDerivatives ERROR : ", &
			 			" geomKind not implemented yet.")
		return
	endif ! PLANAR_GEOM
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
			if ( isNew ) call nearbyParticles%pushBack(faceVerts%int(j))
		enddo
	enddo
end subroutine

subroutine GetParticlesNearVertex( nearbyParticles, aMesh, particleIndex )
	type(STDIntVector), intent(out) :: nearbyParticles
	type(PolyMesh2d), intent(in) :: aMesh
	integer(kint), intent(in) :: particleIndex
	!
	integer(kint) :: i, j, k
	type(STDIntVector) :: adjFaces, faceVerts
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
	enddo
end subroutine

subroutine estimatePartialDerivativesPSE( self, aMesh, scalarField, scalarGrad, secondPartials )
	type(TriScatteredInterp), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: scalarField
	type(Field), intent(inout) :: scalarGrad
	type(Field), intent(inout) :: secondPartials
	!
	integer(kint) :: i, j
	real(kreal) :: meshSize, eps
	real(kreal) :: xi, yi, xj, yj
	
	meshSize = MaxEdgeLength(aMesh%edges, aMesh%particles)
	eps = 3.0_kreal * meshSize

	call SetFieldToZero(scalarGrad)
	call SetFieldToZero(secondPartials)

	do i = 1, aMesh%particles%N
		xi = aMesh%particles%x(i)
		yi = aMesh%particles%y(i)
		
		do j = 1, aMesh%particles%N
			if ( aMesh%particles%isActive(j) ) then
				xj = aMesh%particles%x(j)
				yj = aMesh%particles%y(j)
				
				scalarGrad%xComp(i) = scalarGrad%xComp(i) + aMesh%particles%area(j) * & 
									  (scalarField%scalar(j) + scalarField%scalar(i)) * &
									  dfdxKernel4( xi, yi, xj, yj, eps)
				scalarGrad%yComp(i) = scalarGrad%yComp(i) + aMesh%particles%area(j) * &
									  (scalarField%scalar(j) + scalarField%scalar(i)) * &
									  dfdyKernel4( xi, yi, xj, yj, eps)
			endif			
		enddo
		scalarGrad%xComp(i) = scalarGrad%xComp(i) / eps
		scalarGrad%yComp(i) = scalarGrad%yComp(i) / eps 
	enddo
	
	do i = 1, aMesh%particles%N
		xi = aMesh%particles%x(i)
		yi = aMesh%particles%y(i)
		
		do j = 1, aMesh%particles%N
			if ( aMesh%particles%isActive(j) ) then
				xj = aMesh%particles%x(j)
				yj = aMesh%particles%y(j)
				
				secondPartials%xComp(i) = secondPartials%xComp(i) + aMesh%particles%area(j) * &
						(scalarGrad%xComp(j) + scalarGrad%xComp(i)) * dfdxKernel4( xi, yi, xj, yj, eps )
				secondPartials%yComp(i) = secondPartials%yComp(i) + 0.5_kreal * aMesh%particles%area(j) * &
						(scalarGrad%xComp(j) + scalarGrad%xComp(i)) * dfdyKernel4( xi, yi, xj, yj, eps ) + &
						0.5_kreal * aMesh%particles%area(j) * (scalarGrad%yComp(j)-scalarGrad%yComp(i)) * &
						dfdxKernel4(xi,yi, xj, yj, eps)
				secondPartials%zComp(i) = secondPartials%zComp(i) + aMesh%particles%area(j) * &
						(scalarGrad%yComp(j) + scalarGrad%yComp(i)) * dfdyKernel4(xi, yi, xj, yj, eps )
			endif
		enddo
		secondPartials%xComp(i) = secondPartials%xComp(i) / eps
		secondPartials%yComp(i) = secondPartials%yComp(i) / eps
		secondPartials%zComp(i) = secondPartials%zComp(i) / eps
	enddo
end subroutine

pure function dfdxKernel4( xp, yp, xq, yq, eps)
	real(kreal) :: dfdxKernel4
	real(kreal), intent(in) :: xp, yp, xq, yq, eps
	!
	real(kreal) :: kernelIn(2)
	
	kernelIn = [xp-xq, yp-yq]/eps

	dfdxKernel4 = kernelIn(1) / PI * exp( - sum( kernelIn*kernelIn)) * (-6.0_kreal + 2.0_kreal * sum(kernelIn*kernelIn))
	
	dfdxKernel4 = dfdxKernel4 / eps / eps
end function

pure function dfdyKernel4( xp, yp, xq, yq, eps )
	real(kreal) :: dfdyKernel4
	real(kreal), intent(in) :: xp, yp, xq, yq, eps
	!
	real(kreal) :: kernelIn(2)
	
	kernelIn = [xp-xq, yp-yq]/eps
	
	dfdyKernel4 = kernelIn(2) / PI * exp( -sum(kernelIn*kernelIn)) * (-6.0_kreal + 2.0_kreal * sum(kernelIn*kernelIn))
end function

subroutine estimatePartialDerivativesMLSQ( self, aMesh, scalarField, scalarGrad, secondPartials )
	type(TriScatteredInterp), intent(inout) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	type(Field), intent(in) :: scalarField
	type(Field), intent(inout) :: scalarGrad
	type(Field), intent(inout) :: secondPartials
	!
	type(STDIntVector) :: nearbyParticles
	integer(kint) :: i, j
	integer(kint) :: nn
	integer(kint), parameter :: maxLSQM = 80
	real(kreal) :: A(maxLSQM,10), b(maxLSQM), x, y, coeff(10)
	integer(kint) :: centerParticleIndex
	
	call SetFieldToZero(scalarGrad)
	call SetFieldToZero(secondPartials)
	scalarGrad%N = aMesh%particles%N
	secondPartials%N = aMesh%particles%N

	nn = 10
	
	if ( aMesh%geomKind == PLANAR_GEOM ) then
		!
		! estimate gradient at vertices of mesh triangles
		!

		do i = 1, aMesh%particles%N
			if ( aMesh%particles%isPassive(i) ) then
				A = 0.0_kreal
				b = 0.0_kreal
				call GetParticlesNearVertex( nearbyParticles, aMesh, i )
				if ( nearbyParticles%N > maxLSQM ) then
					call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" MLSQ ERROR: ", "too many equations at vertex.")
					return
				endif
				do j = 1, nearbyParticles%N
					x = aMesh%particles%x(nearbyParticles%int(j))
					y = aMesh%particles%y(nearbyParticles%int(j))
									
					A(j,1) = 1.0_kreal
					A(j,2) = y
					A(j,3) = y * y
					A(j,4) = y * y * y
					A(j,5) = x
					A(j,6) = x * y
					A(j,7) = x * y * y
					A(j,8) = x * x
					A(j,9) = x * x * y
					A(j,10)= x * x * x
					
					b(j) = scalarField%scalar(nearbyParticles%int(j))
				enddo
				
				call lapack_qrSolve( nearbyParticles%N, nn, A, b, coeff )
				
				x = aMesh%particles%x(i)
				y = aMesh%particles%y(i)
				
				scalarGrad%xComp(i) = coeff(5) + coeff(6) * y + coeff(7) * y * y + &
									  2.0_kreal * coeff(8) * x + 2.0_kreal * coeff(9) * x * y + &
									  3.0_kreal * coeff(10) * x * x
				scalarGrad%yComp(i) = coeff(2) + 2.0_kreal * coeff(3) * y + 3.0_kreal * coeff(4) * y * y + &
									  coeff(6) * x + 2.0_kreal * coeff(7) * x * y + &
									  coeff(9) * x * x

				secondPartials%xComp(i) = 2.0_kreal * coeff(8)	+ 2.0_kreal * coeff(9) * y + 6.0_kreal * coeff(10) *x 
				secondPartials%yComp(i) = coeff(6) * 2.0_kreal * coeff(7) * y + 2.0_kreal * coeff(9) * x
				secondPartials%zComp(i) = 2.0_kreal * coeff(3) + 6.0_kreal * coeff(4) * y + 2.0_kreal * coeff(7) * x
			endif ! passive
		enddo ! vertices loop
		
		do i = 1, aMesh%faces%N
			if (.NOT. aMesh%faces%hasChildren(i) ) then
				A = 0.0_kreal
				b = 0.0_kreal
				call GetParticlesNearFace( nearbyParticles, aMesh, i)
				if ( nearbyParticles%N > maxLSQM ) then
					call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//"MLSQ ERROR: ", " too many equations at face.")
					return
				endif
				do j = 1, nearbyParticles%N
					x = aMesh%particles%x(nearbyParticles%int(j))
					y = aMesh%particles%y(nearbyParticles%int(j))
					
					b(j) = scalarField%scalar(nearbyParticles%int(j))
					
					A(j,1) = 1.0_kreal
					A(j,2) = y
					A(j,3) = y * y
					A(j,4) = y * y * y
					A(j,5) = x
					A(j,6) = x * y
					A(j,7) = x * y * y
					A(j,8) = x * x
					A(j,9) = x * x * y
					A(j,10)= x * x * x
				enddo
				
				call lapack_qrSolve( nearbyParticles%N, nn, A, b, coeff )
				
				centerParticleIndex = aMesh%faces%centerParticle(i)
				x = aMesh%particles%x( centerParticleIndex)
				y = aMesh%particles%y( centerParticleIndex)
				
				scalarGrad%xComp(centerParticleIndex) = coeff(5) + coeff(6) * y + coeff(7) * y * y + &
									  2.0_kreal * coeff(8) * x + 2.0_kreal * coeff(9) * x * y + &
									  3.0_kreal * coeff(10) * x * x	
				scalarGrad%yComp(centerParticleIndex) = coeff(2) + 2.0_kreal * coeff(3) * y + &
									  3.0_kreal * coeff(4) * y * y + coeff(6) * x + 2.0_kreal * coeff(7) * x * y + &
									  coeff(9) * x * x	
				secondPartials%xComp(centerParticleIndex) = 2.0_kreal * coeff(8) + 2.0_kreal * coeff(9) * y +&
													        6.0_kreal * coeff(10) *x 
				secondPartials%yComp(centerParticleIndex) = coeff(6) * 2.0_kreal * coeff(7) * y +&
													        2.0_kreal * coeff(9) * x
				secondPartials%zComp(centerParticleIndex) = 2.0_kreal * coeff(3) + 6.0_kreal * coeff(4) * y + &
															2.0_kreal * coeff(7) * x
			endif ! faces%hasChildren(i)
		enddo ! centers loop
	endif ! PLANAR_GEOM	
end subroutine

subroutine lapack_qrSolve( m, n, A, b, x )
	integer(kint), intent(in) :: m
	integer(kint), intent(in) :: n
	real(kreal), intent(inout) :: A(m, n)
	real(kreal), intent(inout) :: b(m)
	real(kreal), intent(out) :: x(n)
	!
	integer(kint) :: nRHS, lda, ldb, jPvt(n), rank, lwork, info
	real(kreal) :: rcond
	real(kreal) :: work(4*n+1)
	
	nRHS = 1
	lda = m
	ldb = m
	lwork = 4*n+1
	jPvt = 0
	
	call DGELSY( m, n, nRHS, A, lda, b, ldb, jPvt, rcond, rank, work, lwork, info)
	
!	if ( info < 0 ) then
!		call LogMessage(log, ERROR_LOGGING_LEVEL, "lapack_qrSolve ERROR in arg: ", -info)
!		return
!	endif
	x = b(1:n)
end subroutine

function InterpolateScalar( self, aMesh, xyIn )
	real(kreal) :: InterpolateScalar
	type(TriScatteredInterp), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	real(kreal), intent(in) :: xyIn(:)
	!
	real(kreal) :: uv(2), u, v
	integer(kint) :: faceIndex
	!real(kreal) :: vmBasis(21)
	real(kreal) :: p0, p1, p2, p3, p4, p5
	
	
	InterpolateScalar = 0.0_kreal
	faceIndex = 0
	if ( .NOT. ( self%mapsReady .AND. self%coeffsReady ) ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" InterpolateScalar ERROR : ",&
										" coefficients and/or maps not ready.")
		return
	endif
	
	faceIndex = LocateFaceContainingPoint(aMesh, xyIn)
	if ( .NOT. pointIsOutsideMesh(aMesh, xyIn) ) then
		
		call convertXYToLocalUV( self, aMesh, faceIndex, xyIn, uv)
		u = uv(1)
		v = uv(2)
!		vmBasis = [1.0_kreal, v, v**2, v**3, v**4, v**5, &
!					u, u*v, u*v**2, u*v**3, u*v**4, &
!					u**2, u**2*v, u**2*v**2, u**2*v**3, &
!					u**3, u**3*v, u**3*v**2, &
!					u**4, u**4*v, &
!					u**5]
		p0 = self%faceCoeffs(1,faceIndex) + v * ( self%faceCoeffs(2,faceIndex) + v * (self%faceCoeffs(3,faceIndex) + &
			 v * (self%faceCoeffs(4,faceIndex) + v * (self%faceCoeffs(5,faceIndex) + v * self%faceCoeffs(6,faceIndex)))))
		p1 = self%faceCoeffs(7,faceIndex) + v * (self%faceCoeffs(8,faceIndex) + v * (self%faceCoeffs(9,faceIndex) + &
			 v * (self%faceCoeffs(10,faceIndex) + v * self%faceCoeffs(11,faceIndex))))
		p2 = self%faceCoeffs(12,faceIndex) + v * (self%faceCoeffs(13,faceIndex) + v * &
			 (self%faceCoeffs(14,faceIndex) + v * self%faceCoeffs(15,faceIndex)))
		p3 = self%faceCoeffs(16,faceIndex) + v * (self%faceCoeffs(17,faceIndex) + v * self%faceCoeffs(18,faceIndex))
		p4 = self%faceCoeffs(19,faceIndex) + v * self%faceCoeffs(20,faceIndex)
		p5 = self%faceCoeffs(21,faceIndex)
		
		InterpolateScalar = p0 + u * (p1 + u*(p2 + u * (p3 + u * (p4 + u*p5))))
		
!		InterpolateScalar = sum( vmBasis * self%faceCoeffs(:,faceIndex))
	endif
end function 

function InterpolateGradient( self, aMesh, xyIn)
	real(kreal) :: InterpolateGradient(2)
	type(TriScatteredInterp), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	real(kreal), intent(in) :: xyIn(:)
	!
	real(kreal) :: uv(2), u, v, vmBasisU(21), vmBasisV(21)
	real(kreal) :: a, b, c, d
	integer(kint) :: faceIndex
	
	InterpolateGradient = 0.0_kreal
	
	if ( .NOT. ( self%mapsReady .AND. self%coeffsReady ) ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" InterpolateGradient ERROR : ",&
										" coefficients and/or maps not ready.")
		return
	endif
	
	faceIndex = LocateFaceContainingPoint(aMesh, xyIn)
	if ( .NOT. pointIsOutsideMesh(aMesh, xyIn) ) then
		
		call convertXYToLocalUV( self, aMesh, faceIndex, xyIn, uv)
		u = uv(1)
		v = uv(2)
		
		vmBasisU(1:6) = 0.0_kreal
		vmBasisU(7:21) = [1.0_kreal, v, v**2, v**3, v**4, &
						  2.0_kreal * u, 2.0_kreal * u * v, 2.0_kreal * u * v**2, 2.0_kreal * u * v**3, &
						  3.0_kreal * u**2, 3.0_kreal * u**2 * v, 3.0_kreal * u**2 * v**2, &
						  4.0_kreal * u**3, 4.0_kreal * u**3*v, &
						  5.0_kreal * u**4]
		vmBasisV = [0.0_kreal, 1.0_kreal, 2.0_kreal * v, 3.0_kreal * v**2, 4.0_kreal * v**3, 5.0_kreal * v**4, &
		            0.0_kreal, u, 2.0_kreal * u * v, 3.0_kreal * u * v**2, 4.0_kreal * u * v**3, &
		            0.0_kreal, u**2, 2.0_kreal * u**2 * v, 3.0_kreal * u**2 * v**2, &
		            0.0_kreal, u**3, 2.0_kreal * u**3 * v, &
		            0.0_kreal, u**4, &
		            0.0_kreal]
		
		
		a = self%faceU(1, faceIndex)
		b = self%faceV(1, faceIndex)
		c = self%faceU(2, faceIndex)
		d = self%faceV(2, faceIndex)
		
		InterpolateGradient(1) = d * sum( vmBasisU * self%faceCoeffs(:,faceIndex)) / (a*d - b*c ) - &
								 c * sum( vmBasisV * self%faceCoeffs(:,faceIndex)) / (a*d - b*c )
		InterpolateGradient(2) =-b * sum( vmBasisU * self%faceCoeffs(:,faceIndex)) / (a*d - b*c ) + &
		 						 a * sum( vmBasisV * self%faceCoeffs(:,faceIndex)) / (a*d - b*c )
	endif
end function

function InterpolateLaplacian( self, aMesh, xyIn )
	real(kreal) :: InterpolateLaplacian
	type(TriScatteredInterp), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	real(kreal), intent(in) :: xyIn(:)
	!
	real(kreal) :: uv(2), u, v, vmBasisUU(21), vmBasisUV(21), vmBasisVV(21)
	real(kreal) :: a, b, c, d, det2
	integer(kint) :: faceIndex
	
	InterpolateLaplacian = 0.0_kreal
	
	if ( .NOT. ( self%mapsReady .AND. self%coeffsReady ) ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logkey)//" InterpolateLaplacian ERROR : ",&
										" coefficients and/or maps not ready.")
		return
	endif
	
	faceIndex = LocateFaceContainingPoint(aMesh, xyIn)
	if ( .NOT. pointIsOutsideMesh(aMesh, xyIn) ) then
	
		call convertXYToLocalUV( self, aMesh, faceIndex, xyIn, uv)
		
		vmBasisUU(1:11) = 0.0_kreal
		vmBasisUU(12:21) = [2.0_kreal, 2.0_kreal * v, 2.0_kreal * v**2, 2.0_kreal * v**3, &
						    6.0_kreal * u, 6.0_Kreal * u * v, 6.0_kreal * u * v**2, &
						    12.0_kreal * u**2, 12.0_kreal * u**2 * v, &
						    20.0_kreal * u**3]
		vmBasisUV(1:7) = 0.0_kreal
		vmBasisUV(8:21) = [1.0_kreal, 2.0_kreal * v, 3.0_kreal * v**2, 4.0_kreal * v**3, &
						   0.0_kreal, 2.0_kreal * u, 4.0_kreal * u * v, 6.0_kreal * u * v**2, &
						   0.0_kreal, 3.0_kreal * u**2, 6.0_kreal * u**2 * v, &
						   0.0_kreal, 4.0_kreal * u**3, 0.0_kreal]
		vmBasisVV = [ 0.0_kreal, 0.0_kreal, 2.0_kreal, 6.0_kreal * v, 12.0_kreal * v**2, 20.0_kreal * v**3, &
					  0.0_kreal, 0.0_kreal, 2.0_kreal * u, 6.0_kreal * u * v, 12.0_kreal * u * v**2, &
					  0.0_kreal, 0.0_kreal, 2.0_kreal * u**2, 6.0_kreal*u**2 * v, &
					  0.0_kreal, 0.0_kreal, 2.0_kreal * u**3, &
					  0.0_kreal, 0.0_kreal, 0.0_kreal]
					  
		a = self%faceU(1, faceIndex)
		b = self%faceV(1, faceIndex)
		c = self%faceU(2, faceIndex)
		d = self%faceV(2, faceIndex)
		
		det2 = (a*d-b*c)*(a*d-b*c)
		
		InterpolateLaplacian = 	d*d * sum( vmBasisUU * self%faceCoeffs(:,faceIndex))/det2 - &
					2.0_kreal * c*d * sum( vmBasisUV * self%faceCoeffs(:,faceIndex))/det2 + &
								c*c * sum( vmBasisVV * self%faceCoeffs(:,faceIndex))/det2 + &
								b*b * sum( vmBasisUU * self%faceCoeffs(:,faceIndex))/det2 - &
					2.0_kreal * a*b * sum( vmBasisUV * self%faceCoeffs(:,faceIndex))/det2 + &
								a*a * sum( vmBasisVV * self%faceCoeffs(:,faceIndex))/det2
	endif
end function

subroutine convertXYToLocalUV(self, aMesh, faceIndex, xyIn, uvOut)
	type(TriScatteredInterp), intent(in) :: self
	type(PolyMesh2d), intent(in) :: aMesh
	integer(kint), intent(in) :: faceIndex
	real(kreal), intent(in) :: xyIn(:)
	real(kreal), intent(out) :: uvOut(2)
	!
	real(kreal) :: v0(3), a, b, c, d
	
	if ( .NOT. self%mapsReady ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, trim(logKey)//" convertXYToLocalUV ERROR:",&
			" cannot use local maps before they are defined.")
		return
	endif
	
	v0 = PhysCoord(amesh%particles, amesh%faces%vertices(1,faceIndex))
	
	a = self%faceU(1, faceIndex)
	b = self%faceV(1, faceIndex)
	c = self%faceU(2, faceIndex)
	d = self%faceV(2, faceIndex)
	
	uvOut(1) = ( d * ( xyIn(1) - v0(1) ) - b * (xyIn(2) - v0(2)))/(a*d - b*c)
	uvOut(2) = (-c * ( xyIn(1) - v0(1) ) + a * (xyIn(2) - v0(2)))/(a*d - b*c)
end subroutine

subroutine InitLogger(aLog,rank)
! Initialize a logger for this module and processor
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.3,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

!> @}
end module
