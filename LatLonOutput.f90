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
public LLOutputMatlab, UpdateFilename

type LLSource
	character(len = 256) :: filename
	character(len = 128) :: title
	type(STRIPACKData) :: delTri
	type(SSRFPACKData) :: vectorSource
	type(SSRFPACKData), pointer :: scalarSource => null()
	integer(kint) :: nTracer
	integer(kint) :: problemKind
	type(SSRFPACKData), pointer :: tracerSource => null()
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

	self%filename = filename
	self%nLon = nLon

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
	call DelaunayTriangulation(self%delTri)

	call New(self%vectorSource,self%delTri,.TRUE.)

	self%problemKind = aMesh%problemKind
	if ( self%problemKind /= ADVECTION_SOLVER ) then
		allocate(self%scalarSource)
		call New(self%scalarSource,self%delTri,.FALSE.)
	endif

	self%nTracer = aMesh%nTracer
	if ( self%nTracer > 0 ) then
		allocate(self%tracerSource)
		call New(self%tracerSource, self%delTri, self%nTracer)
	endif

end subroutine

subroutine DeletePrivate(self)
	type(LLSource), intent(inout) :: self
	deallocate(self%lats)
	deallocate(self%lons)
	call Delete(self%deltri)
	call Delete(self%scalarSource)
	call Delete(self%vectorSource)
	if ( self%nTracer > 0 ) then
		call Delete(self%tracerSource)
		deallocate(self%tracerSource)
	endif
	if ( associated(self%scalarSource) ) then
		call Delete(self%scalarSource)
		deallocate(self%scalarSource)
	endif
end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine UpdateFileNameLL(self,filename)
	type(LLSource), intent(inout) :: self
	character(len=*), intent(in) :: filename
	self%filename = trim(filename)
end subroutine

subroutine LLOutputMatlab(self,aMesh)
	type(LLSource), intent(inout) :: self
	type(SphereMesh), intent(in) :: aMesh
	!
	integer(kint) :: i, j, k, writeStat
	real(kreal) :: interpVector(3), interpScalar, xyzLoc(3), lonUnitVector(3), latUnitVector(3)
	real(kreal), allocatable :: llgrid1(:,:), llgrid2(:,:)

	! open output file
	open(unit = WRITE_UNIT_1, file = self%filename, status = 'REPLACE', action = 'WRITE', iostat = writeStat)
	if ( writeStat /= 0 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'LLOutputMatlab : ', 'ERROR opening output file.')
		return
	endif

	! write lon/lat vectors
	write(WRITE_UNIT_1,'(A)',ADVANCE = 'NO') 'lon = [ '
	do i=1,self%nlon-1
		write(WRITE_UNIT_1,'(F24.15,A)') self%lons(i), '; ...'
	enddo
	write(WRITE_UNIT_1,'(F24.15,A)') self%lons(self%nLon), ' ] ; '

	write(WRITE_UNIT_1,'(A)',advance = 'NO') 'lat = [ '
	do i=1,self%nlon/2
		write(WRITE_UNIT_1,'(F24.15,A)') self%lats(i), '; ...'
	enddo
	write(WRITE_UNIT_1,'(F24.15,A)') self%lats(self%nlon/2+1), ' ] ;'

	allocate(llgrid1(self%nLon/2+1,self%nLon))
	llgrid1 = 0.0_kreal
	allocate(llgrid2(self%nLon/2+1,self%nLon))
	llgrid2 = 0.0_kreal
	!
	! Lagrangian parameter
	!
	call SetSourceLagrangianParameter(self%vectorSource,self%delTri)
	do j=1,self%nLon
		do i=1,self%nLon/2+1
			xyzLoc = EARTH_RADIUS*[ cos(self%lats(i))*cos(self%lons(j)), &
						cos(self%lats(i))*sin(self%lons(j)), sin(self%lats(i)) ]
			interpVector = InterpolateVector(xyzLoc, self%vectorSource, self%delTri)
			llgrid1(i,j) = Longitude(interpVector)
			llgrid2(i,j) = Latitude(interpVector)
		enddo
	enddo

	write(WRITE_UNIT_1,'(A)',advance = 'NO') ' alphaLongitude = [ '
	do i=1,self%nLon/2
		do j=1,self%nLon -1
			write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid1(i,j)
		enddo
		write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid1(i,self%nLon), ' ; ...'
	enddo
	do j=1,self%nLon -1
		write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid1(self%nLon/2+1,j)
	enddo
	write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid1(self%nLon/2+1,self%nLon), ' ] ; '

	write(WRITE_UNIT_1,'(A)',advance = 'NO') ' alphaLatitude = [ '
	do i=1,self%nLon/2
		do j=1,self%nLon -1
			write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid2(i,j)
		enddo
		write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid2(i,self%nLon), ' ; ...'
	enddo
	do j=1,self%nLon -1
		write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid2(self%nLon/2+1,j)
	enddo
	write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid2(self%nLon/2+1,self%nLon), ' ] ; '

	!
	!	velocity
	!
	call SetSourceVelocity(self%vectorSource,self%delTri)
	llgrid1 = 0.0_kreal
	llgrid2 = 0.0_kreal
	do j=1,self%nLon
		do i=2,self%nLon/2
			xyzLoc = EARTH_RADIUS*[ cos(self%lats(i))*cos(self%lons(j)), &
						cos(self%lats(i))*sin(self%lons(j)), sin(self%lats(i)) ]
			interpVector = InterpolateVector(xyzLoc, self%vectorSource, self%delTri)

			lonUnitVector = [ - xyzLoc(2), &
								xyzLoc(1), 0.0_kreal ]
			lonUnitVector = lonUnitVector/(EARTH_RADIUS*sqrt(EARTH_RADIUS*EARTH_RADIUS - xyzLoc(3)*xyzLoc(3)))

			latUnitVector = [ -xyzLoc(1)*xyzLoc(3), &
							  -xyzLoc(2)*xyzLoc(3), &
							  EARTH_RADIUS*EARTH_RADIUS - xyzLoc(3)*xyzLoc(3) ]
			latUnitVector = latUnitVector/(EARTH_RADIUS*EARTH_RADIUS*sqrt(EARTH_RADIUS*EARTH_RADIUS - xyzLoc(3)*xyzLoc(3)))

			llgrid1(i,j) = sum(interpVector*lonUnitVector)

			llgrid2(i,j) = sum(interpVector*latUnitVector)

		enddo
	enddo

	write(WRITE_UNIT_1,'(A)',advance = 'NO') ' zonalVelocityU = [ '
	do i=1,self%nLon/2
		do j=1,self%nLon -1
			write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid1(i,j)
		enddo
		write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid1(i,self%nLon), ' ; ...'
	enddo
	do j=1,self%nLon -1
		write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid1(self%nLon/2+1,j)
	enddo
	write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid1(self%nLon/2+1,self%nLon), ' ] ; '

	write(WRITE_UNIT_1,'(A)',advance = 'NO') ' meridionalVelocityV = [ '
	do i=1,self%nLon/2
		do j=1,self%nLon -1
			write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid2(i,j)
		enddo
		write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid2(i,self%nLon), ' ; ...'
	enddo
	do j=1,self%nLon -1
		write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid2(self%nLon/2+1,j)
	enddo
	write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid2(self%nLon/2+1,self%nLon), ' ] ; '

	!
	!	tracers
	!
	if ( self%nTracer > 0 ) then
		call SetSourceTracer(self%tracerSource, self%delTri)
		llgrid1 = 0.0_kreal
		do k = 1, self%nTracer
			do j=1,self%nLon
				do i=1,self%nLon/2 + 1
					xyzLoc = EARTH_RADIUS*[ cos(self%lats(i))*cos(self%lons(j)), &
						cos(self%lats(i))*sin(self%lons(j)), sin(self%lats(i)) ]

					llgrid1(i,j) = InterpolateTracer(xyzLoc, self%tracerSource, self%delTri, k)
				enddo
			enddo

			write(WRITE_UNIT_1,'(A,I1,A)',advance = 'NO') ' tracer', k,' = [ '
			do i=1,self%nLon/2
				do j=1,self%nLon -1
					write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid1(i,j)
				enddo
				write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid1(i,self%nLon), ' ; ...'
			enddo
			do j=1,self%nLon -1
				write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid1(self%nLon/2+1,j)
			enddo
			write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid1(self%nLon/2+1,self%nLon), ' ] ; '
		enddo
	endif

	!
	!	bve variables
	!
	if ( self%problemKind == BVE_SOLVER) then
		!
		!	absolute vorticity
		!
		llgrid1 = 0.0_kreal
		call SetSourceAbsVort(self%scalarSource,self%delTri)
		do j=1,self%nLon
			do i=1,self%nLon/2 + 1
				xyzLoc = EARTH_RADIUS*[ cos(self%lats(i))*cos(self%lons(j)), &
					cos(self%lats(i))*sin(self%lons(j)), sin(self%lats(i)) ]

				llgrid1(i,j) = InterpolateScalar(xyzLoc, self%scalarSource, self%delTri)
			enddo
		enddo

		write(WRITE_UNIT_1,'(A)',advance = 'NO') ' absVort = [ '
		do i=1,self%nLon/2
			do j=1,self%nLon -1
				write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid1(i,j)
			enddo
			write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid1(i,self%nLon), ' ; ...'
		enddo
		do j=1,self%nLon -1
			write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid1(self%nLon/2+1,j)
		enddo
		write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid1(self%nLon/2+1,self%nLon), ' ] ; '

		!
		!	relative vorticity
		!
		llgrid1 = 0.0_kreal
		call SetSourceRelVort(self%scalarSource,self%delTri)
		do j=1,self%nLon
			do i=1,self%nLon/2 + 1
				xyzLoc = EARTH_RADIUS*[ cos(self%lats(i))*cos(self%lons(j)), &
					cos(self%lats(i))*sin(self%lons(j)), sin(self%lats(i)) ]

				llgrid1(i,j) = InterpolateScalar(xyzLoc, self%scalarSource, self%delTri)
			enddo
		enddo
		write(WRITE_UNIT_1,'(A)',advance = 'NO') ' relVort = [ '
		do i=1,self%nLon/2
			do j=1,self%nLon -1
				write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid1(i,j)
			enddo
			write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid1(i,self%nLon), ' ; ...'
		enddo
		do j=1,self%nLon -1
			write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid1(self%nLon/2+1,j)
		enddo
		write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid1(self%nLon/2+1,self%nLon), ' ] ; '

		!
		!	kinetic energy
		!
		llgrid1 = 0.0_kreal
		call SetSourceKineticEnergy(self%scalarSource,self%delTri)
		do j=1,self%nLon
			do i=1,self%nLon/2 + 1
				xyzLoc = EARTH_RADIUS*[ cos(self%lats(i))*cos(self%lons(j)), &
					cos(self%lats(i))*sin(self%lons(j)), sin(self%lats(i)) ]

				llgrid1(i,j) = InterpolateScalar(xyzLoc, self%scalarSource, self%delTri)
			enddo
		enddo
		write(WRITE_UNIT_1,'(A)',advance = 'NO') ' ke = [ '
		do i=1,self%nLon/2
			do j=1,self%nLon -1
				write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid1(i,j)
			enddo
			write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid1(i,self%nLon), ' ; ...'
		enddo
		do j=1,self%nLon -1
			write(WRITE_UNIT_1,'(F24.10)',advance='NO') llgrid1(self%nLon/2+1,j)
		enddo
		write(WRITE_UNIT_1,'(F24.10, A)',advance='YES') llgrid1(self%nLon/2+1,self%nLon), ' ] ; '
	endif

	!
	!	TO DO : SWE VARIABLES
	!


	deallocate(llgrid1)
	deallocate(llgrid2)
	close(WRITE_UNIT_1)
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
