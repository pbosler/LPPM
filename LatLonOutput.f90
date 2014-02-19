module LatLonOutputModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	Defines output data structures and methods for plotting Lagrangian meshes of the sphere.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
!----------------
use NumberKindsModule
use SphereGeomModule
use IntegerListModule
use LoggerModule
use ParticlesModule
use EdgesModule
use PanelsModule
use SphereMeshModule
use STRIPACKInterfaceModule
use SSRFPACKInterfaceModule


implicit none
private
public LLSource
public New, Delete
public LLOutput, UpdateFilename

type LLSource
	character(len = 256) :: filename
	character(len = 128) :: title
	type(STRIPACKData) :: delTri
	type(SSRFPACKData) :: vectorSource
	type(SSRFPACKData) :: scalarSource
	integer(kint) :: nLon
	real(kreal), pointer :: lats(:) => null(),&
				   			lons(:) => null()
end type
!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'LatLonOutput'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
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

interface UpdateFilename
	module procedure UpdateFilenameLL
end interface

contains
!
!----------------
! Standard methods : Constructor / Destructor, Copy
!----------------
!

subroutine NewPrivate(self,aMesh,filename,nLon)
	type(LLSource), intent(out) :: self
	type(SphereMesh), intent(in) :: aMesh
	character(len=*), intent(in) :: filename
	integer(kint), intent(in), optional :: nLon
	! local variables
	real(kreal) :: dLambda
	integer(kint) :: j

	if (.not. logInit) call InitLogger(log,procRank)
	if ( present(nLon) ) then
		dLambda = 2.0_kreal*PI/real(nLon,kreal)
	else
		dLambda = PI/180.0_kreal
	endif
	allocate(self%lons(nLon))
	allocate(self%lats(nLon/2 + 1))
	do j=1,nLon
		self%lons(j) = (j-1)*dLambda
	enddo
	do j=1,nLon/2+1
		self%lats(j) = -PI/2.0_kreal + (j-1)*dLambda
	enddo
	call New(self%delTri,aMesh)
	call DelaunayTriangulation(delTri)
end subroutine

subroutine DeletePrivate(self)
	type(LLSource), intent(inout) :: self
	deallocate(self%lats)
	deallocate(self%lons)
	call Delete(self%deltri)
	call Delete(self%scalarSource)
	call Delete(self%vectorSource)
end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine UpdateFileNameLL(self,filename)
	type(LLSource), intent(inout) ::: self
	character(len=*), intent(in) :: filename
	self%filename = trim(filename)
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
