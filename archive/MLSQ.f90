module MovingLeastSquaresModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	Defines the moving least squares data structure and methods used by SphereMesh
!	for approximating functions on the sphere.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
!----------------
use NumberKindsModule
use LoggerModule
use SphereGeomModule
use ParticlesModule
use EdgesModule
use PanelsModule
use SphereMeshModule
use LAPACK95

implicit none

include 'mpif.h'

private
public MLSQData
public New, Delete


!
!----------------
! Types and module constants
!----------------
!
type MLSQData
	integer(kint), pointer :: adjacentPanels(:,:) => null(), &
							  nAdjacentPanels(:) => null(), &
							  nearbyParticles(:,:) => null(), &
							  nNearbyParticles(:,:) => null()
	real(kreal), pointer :: xCoords(:,:)=>null(), &
							yCoords(:,:)=>null(), &
							zCoords(:,:)=>null()
	real(kreal), pointer :: uData(:,:)=>null(),&
						 	vData(:,:)=>null(),&
						 	wData(:,:)=>null(),&
						 	hData(:,:)=>null()
	real(kreal), pointer :: uCoeff(:,:)=>null(),&
						 	vCoeff(:,:)=>null(),&
						 	wCoeff(:,:)=>null(),&
						 	hCoeff(:,:)=>null()
end type
!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'MLSQ'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=128) :: logstring
character(len=24) :: formatstring
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
subroutine NewPrivate(self, aMesh)
	type(MLSQData), intent(out) :: self
	type(SphereMesh), intent(in) :: aMesh
	!
	integer(kint) :: nActive, i, j, k, insertPoint
	type(Panels), pointer :: aPanels
	
	aPanels => aMesh%panels
	nPanels = aPanels%N
		
	allocate(self%adjacentPanels(8,nPanels))
	allocate(self%nAdjacentPanels(nPanels))
	allocate(self%nearbyParticles(20,nPanels))
	allocate(self%nNearbyParticles(nPanels))
	self%adjacentPanels = 0
	self%nAdjacentPanels = 0
	self%nearbyParticles = 0
	self%nNearbyParticles = 0
	
!	allocate(self%xCoords(29,nPanels))
!	allocate(self%yCoords(29,nPanels))
!	allocate(self%zCoords(29,nPanels))
!	self%xCoords = 0.0_kreal
!	self%yCoords = 0.0_kreal
!	self%zCoords = 0.0_kreal
!	
!	allocate(self%uData(29,nPanels))
!	allocate(self%vData(29,nPanels))
!	allocate(self%wData(29,nPanels))
!	allocate(self%hData(29,nPanels))
!	self%uData = 0.0_kreal
!	self%vData = 0.0_kreal
!	self%wData = 0.0_kreal
!	self%hData = 0.0_kreal
!		
!	allocate(self%uCoeff(10,nPanels))
!	allocate(self%vCoeff(10,nPanels))
!	allocate(self%wCoeff(10,nPanels))
!	allocate(self%hCoeff(10,nPanels))
!	self%uCoeff = 0.0_kreal
!	self%vCoeff = 0.0_kreal
!	self%wCoeff = 0.0_kreal
!	self%hCoeff = 0.0_kreal
	
	insertPoint = 1
	do i=1,aPanel%N
		if ( .NOT. aPanels%hasChildren(i) ) then
			call FindNearbyParticlesToPanel(self%nearbyParticles(:,insertPoint), self%nNearbyParticles(inserPoint), &
				self%adjacentPanels(:,insertPoint), self%nAdjacentPanels(insertPoint), aMesh, i)
			insertPoint = insertPoint + 1
		endif
	enddo	
end subroutine
!
!----------------
! Public member functions
!----------------
!
subroutine SetMLSQCoefficients(self,particlesXYZIn,panelsXYZIn, particlesVelocityIn,panelsVelocityIn,&
			particlesHIn, panelsHIn, panelsAreaIn)
	type(MLSQData), intent(inout) :: self
	real(kreal), intent(in) :: particlesXYZIn(:,:), activePanelsXYZIn(:,:), particlesVelocityIn(:,:)
	real(kreal), intent(in) :: activePanelsVelocityIn(:,:), particlesHIn(:), activePanelsHIn(:), activePanelsAreaIn(:)
	!
	integer(kint) :: nLocations
	real(kreal) :: uavg, vavg, wavg, havg
	real(kreal) :: A(31,10), dataLocs(3,31), uData(31), vData(31), wData(31), hData(31)
	
	
	
	uavg = sum(activePanelsVelocityIn(1,:)*activePanelsAreaIn/EARTH_RADIUS)/(4.0_kreal*EARTH_RADIUS)
	vavg = sum(activePanelsVelocityIn(2,:)*activePanelsAreaIn/EARTH_RADIUS)/(4.0_kreal*EARTH_RADIUS)
	wavg = sum(activePanelsVelocityIn(3,:)*activePanelsAreaIn/EARTH_RADIUS)/(4.0_kreal*EARTH_RADIUS)
	havg = sum(activePanelsHIn(1,:)*activePanelsAreaIn/EARTH_RADIUS)/(4.0_kreal*EARTH_RADIUS)
	
	do j=1,nPanels
		if ( .NOT. aPanels%hasChildren(j)) then
			dataLocs = 0.0_kreal
			do i=1,self%nNearbyParticles(j)
				dataLocs(:,i) = particlesXYZIn(:,self%nearbyParticles(i,j))
				uData(i) = particlesVelocityIn(1,self%nearbyParticles(i,j))
				vData(i) = particlesVelocityIn(2,self%nearbyParticles(i,j))
				wData(i) = particlesVelocityIn(3,self%nearbyParticles(i,j))
				hData(i) = particlesHIn(self%nearbyParticles(i,j))
			enddo
			do i=1,self%nAdjacentPanels(j)
				dataLocs(:,self%nNearbyParticles(j) + i) = panelsXYZIn(:,self%adjacentPanels(i,j))
				uData(self%nNearbyParticles(j) + i) = panelsVelocity(1,self%adjacentPanels(i,j))
				vData(self%nNearbyParticles(j) + i) = panelsVelocity(2,self%adjacentPanels(i,j))
				wData(self%nNearbyParticles(j) + i) = panelsVelocity(3,self%adjacentPanels(i,j))
				hData(self%nNearbyParticles(j) + i) = panelsH(self%adjacentPanels(i,j))
			enddo
			dataLocs(:,self%nNearbyParticles(j) + self%nAdjacentPanels(j) + 1) = panelsXYZIn(:,j)
			uData(self%nNearbyParticles(j) + self%nAdjacentPanels(j) + 1) = panelsVelocity(1,j)
			vData(self%nNearbyParticles(j) + self%nAdjacentPanels(j) + 1) = panelsVelocity(2,j)
			wData(self%nNearbyParticles(j) + self%nAdjacentPanels(j) + 1) = panelsVelocity(3,j)
			hData(self%nNearbyParticles(j) + self%nAdjacentPanels(j) + 1) = panelsH(1,j)
			
			dataLocs(:,self%nNearbyParticles(j) + self%nAdjacentPanels(j) + 2) = 0.0_kreal
			uData(self%nNearbyParticles(j) + self%nAdjacentPanels + 2 ) = uavg
			vData(self%nNearbyParticles(j) + self%nAdjacentPanels + 2 ) = vavg
			wData(self%nNearbyParticles(j) + self%nAdjacentPanels + 2 ) = wavg
			hData(self%nNearbyParticles(j) + self%nAdjacentPanels + 2 ) = havg
			
			dataLocs(:,self%nNearbyParticles(j) + self%nAdjacentPanels(j) + 3) = 2.0_kreal*panelsXYZ(:,j)
			uData(self%nNearbyParticles(j) + self%nAdjacentPanels + 3 ) = uavg
			vData(self%nNearbyParticles(j) + self%nAdjacentPanels + 3 ) = vavg
			wData(self%nNearbyParticles(j) + self%nAdjacentPanels + 3 ) = wavg
			hData(self%nNearbyParticles(j) + self%nAdjacentPanels + 3 ) = havg
			
			nLocations = self%nNearbyParticles(j) + self%nAdjacentPanels(j) + 3
			
			A = 0.0_kreal
			A(1:nLocations, 1 ) = 1.0_kreal				
			A(1:nLocations, 2 ) = dataLocs(1,1:nPoints) ! x
			A(1:nLocations, 3 ) = dataLocs(2,1:nPoints) ! y 
			A(1:nLocations, 4 ) = dataLocs(3,1:nPoints) ! z
			A(1:nLocations, 5 ) = dataLocs(1,1:nPoints) * dataLocs(2,1:nPoints) ! xy
			A(1:nLocations, 6 ) = dataLocs(1,1:nPoints) * dataLocs(3,1:nPoints) ! xz
			A(1:nLocations, 7 ) = dataLocs(2,1:nPoints) * dataLocs(3,1:nPoints) ! yz
			A(1:nLocations, 8 ) = dataLocs(1,1:nPoints) * dataLocs(1,1:nPoints)	! x^2
			A(1:nLocations, 9 ) = dataLocs(2,1:nPoints) * dataLocs(2,1:nPoints) ! y^2
			A(1:nLocations, 10) = dataLocs(3,1:nPoints) * dataLocs(3,1:nPoints) ! z^2
			
			AT = transpose(A)
			ATA = matmul(AT,A)
			
			! find cholesky factorization with lapack
			call DPOTRF('U',10,ATA,10,errCode)
			if ( errCode < 0 ) then
				call LogMessage(log,ERROR_LOGGING_LEVEL,trim(logKey)//' potrf ERROR : found illegal value at position ',errCode)
				call LogMessage(log,ERROR_LOGGING_LEVEL,trim(logKey)//' error found at panel ',j)
			elseif ( errCode > 0 ) then
				call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'potrf ERROR : non-symmetric positive definite matrix.')
			endif
			
			! invert ATA with lapack
			call DPOTRI('U',10,ATA,10,errCode)
			if ( errCode < 0 ) then
				call LogMessage(log,ERROR_LOGGING_LEVEL,trim(logkey)//' potri ERROR : found illegal value at position ',errCode)
				call LogMessage(log,ERROR_LOGGING_LEVEL,trim(logKey)//' error found at panel ',j)
			elseif (errCode > 0 ) then
				call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'potri ERROR : found zero on diagonal of Cholesky matrix')
				call LogMessage(log,ERROR_LOGGING_LEVEL,trim(logKey)//' error found at panel ',j)
			endif
			
			! Unpack LAPACK'S symmetric storage
			do m=1,10
				do k=m+1,10
					ata(k,m) = ata(m,k)
				enddo
			enddo
			
			ATAInvAT = matmul(ATA,AT)
			
			uCoeff = matmul(ATAInvAT(:,1:nLocations),uData(1:nLocations))
			vCoeff = matmul(ATAInvAT(:,1:nLocations),vData(1:nLocations))
			wCoeff = matmul(ATAInvAT(:,1:nLocations),wData(1:nLocations))
			hCoeff = matmul(ATAInvAT(:,1:nLocations),hData(1:nLocations))
			
		endif
	enddo
	
end subroutine

!
!----------------
! Module methods : module- or type-specific private functions
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