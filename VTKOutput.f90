module VTKOutputModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method - Spherical Model
!------------------------------------------------------------------------------
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
!
!> @defgroup VTKOutput VTKOutput - Sphere
!> Provides data structures and methods for outputting data from SphereMesh objects into VTK 2.0 file format.
!
!
! DESCRIPTION:
!> @file
!> Provides data structures and methods for outputting data from SphereMesh objects into VTK 2.0 file format.
!
!------------------------------------------------------------------------------
use NumberKindsModule
use SphereGeomModule
use IntegerListModule
use LoggerModule
use ParticlesModule
use EdgesModule
use PanelsModule
use SphereMeshModule


implicit none
private
public VTKSource
public New, Delete
public VTKOutput, UpdateFilename, UpdateTitle
public VTKOutputMidpointRule

type VTKSource
	character(len = 256) :: filename
	character(len = 128) :: title
	integer(kint) :: nScalars, &
				     nPoints, &
				     nCells, &
				     cellListSize
end type
!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'VTKOutput'
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
	module procedure UpdateFilenameVTK
end interface

contains
!
!----------------
! Standard methods : Constructor / Destructor, Copy
!----------------
!

subroutine NewPrivate(self,aMesh,filename,title)
	! calling parameters
	type(VTKSource), intent(out) :: self
	type(SphereMesh), intent(in) :: aMesh
	character(len=*), intent(in) :: filename
	character(len=*), intent(in), optional :: title
	! local variables
	type(Particles), pointer :: aParticles
	type(Edges), pointer :: anEdges
	type(Panels), pointer :: aPanels
	integer(kint) :: j, panelKind
	integer(kint) :: edgeList(8), vertList(8), nVerts

	aParticles=>aMesh%particles
	anEdges=>aMesh%edges
	aPanels=>aMesh%panels

	self%nScalars = 2	! lagrangian parameter and velocity
	self%nPoints = aParticles%N + aPanels%N_Active
	self%nCells = 0
	self%cellListSize = 0

	panelKind = GetPanelKind(aPanels)

	! Count number of subtriangles
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			call CCWEdgesAndParticlesAroundPanel(edgeList,vertList,nVerts,aMesh,j)
			self%nCells = self%nCells + nVerts
			self%cellListSize = self%cellListSize + 4*nVerts
		endif
	enddo

	self%filename = filename
	if ( present(title) ) then
		self%title = title
	else
		self%title = '_'
	endif
	if ( associated(aPanels%absVort) ) self%nScalars = self%nScalars+1
	if ( associated(aPanels%relVort) ) self%nScalars = self%nScalars+1
	if ( associated(aPanels%potVort) ) self%nScalars = self%nScalars+1
	if ( associated(aPanels%h) ) self%nScalars = self%nScalars+1
	if ( associated(aPanels%ke)) self%nScalars = self%nScalars+1
	if ( associated(aPanels%div) ) self%nScalars = self%nScalars+1
	self%nScalars = self%nScalars + GetNTracer(aPanels)
end subroutine

subroutine DeletePrivate(self)
	type(VTKSource), intent(inout) :: self
	self%nScalars = 0
	self%nPoints = 0
	self%nCells = 0
	self%title = ' '
	self%filename = ' '
end subroutine

!
!----------------
! Public functions
!----------------
!

subroutine UpdateFilenameVTK(self,newfilename)
	type(VTKSource), intent(inout) :: self
	character(len=*), intent(in) :: newFilename
	self%filename = trim(newFilename)
end subroutine

subroutine UpdateTitle(self, newtitle)
	type(VTKSource), intent(inout) :: self
	character(len=*), intent(in) :: newTitle
	self%title = newtitle
end subroutine

subroutine vtkOutput(self,aMesh)
	type(VTKSource), intent(in) :: self
	type(SphereMesh), intent(in) :: aMesh
	integer(kint) :: writeStat, j, k
	type(Particles), pointer :: aParticles
	type(Edges), pointer :: anEdges
	type(Panels), pointer :: aPanels
	integer(kint) :: edgeList(8), vertList(8), nVerts, panelKind, panelK, nTracer
	character(len=28) :: dataString
	real(kreal), allocatable :: tempArea(:)

	open(unit = WRITE_UNIT_1,file=self%filename,status='REPLACE',action='WRITE',iostat=writeStat)
	if ( writeStat /= 0 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'VTKOutput ERROR opening output file.')
		return
	endif

	aParticles=>aMesh%particles
	anEdges=>aMesh%edges
	aPanels=>aMesh%panels
	vertList = -1
	panelKind = GetPanelKind(aPanels)
	!
	! VTK File header
	!
	write(WRITE_UNIT_1,'(A)') '# vtk DataFile Version 2.0'
	write(WRITE_UNIT_1,'(A)') self%title
	write(WRITE_UNIT_1,'(A)') 'ASCII'
	write(WRITE_UNIT_1,'(A)') 'DATASET POLYDATA'
	!
	! VTK points
	!
	write(WRITE_UNIT_1,'(A,I8,A)') 'POINTS ',self%nPoints,' double'
	do j=1,aparticles%N
		write(WRITE_UNIT_1,'(3F24.8)') aParticles%x(1,j)/EARTH_RADIUS,aParticles%x(2,j)/EARTH_RADIUS,aParticles%x(3,j)/EARTH_RADIUS
	enddo
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			write(WRITE_UNIT_1,'(3F24.8)') aPanels%x(1,j)/EARTH_RADIUS, aPanels%x(2,j)/EARTH_RADIUS, aPanels%x(3,j)/EARTH_RADIUS
		endif
	enddo

	write(WRITE_UNIT_1,'(A,I8,I8)') 'POLYGONS ',self%nCells, self%cellListSize
	! Draw subtriangles
	panelK = 0
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			call CCWEdgesAndParticlesAroundPanel(edgeList,vertList,nVerts,aMesh,j)
			vertList = vertList - 1
			do k=1,nVerts-1
				write(WRITE_UNIT_1,'(4I8)') 3, vertList(k),aParticles%N + panelK, vertList(k+1)
			enddo
			write(WRITE_UNIT_1,'(4I8)') 3, vertList(nverts), aparticles%N + panelK, vertList(1)
			panelK = panelK + 1
		endif
	enddo

	if ( self%nScalars > 0 ) then
		write(WRITE_UNIT_1,'(A,I8)') 'POINT_DATA  ',self%nPoints
		datastring = 'SCALARS lagParam double 3'
		write(WRITE_UNIT_1,'(A)') trim(dataString)
		write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
		do j=1,aParticles%N
			write(WRITE_UNIT_1,'(3F24.8)') aParticles%x0(:,j)/EARTH_RADIUS
		enddo
		do j=1,aPanels%N
			if (.NOT. aPanels%hasChildren(j) ) then
				write(WRITE_UNIT_1,'(3F24.8)') aPanels%x0(:,j)/EARTH_RADIUS
			endif
		enddo

		if ( associated(aPanels%absVort)) then
			dataString = 'SCALARS  AbsVort  double 1'
			write(WRITE_UNIT_1,'(A)') dataString
			write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
			do j=1,aParticles%N
				write(WRITE_UNIT_1,'(F24.15)') aParticles%absVort(j)
			enddo
			do j=1,aPanels%N
				if ( .NOT. aPanels%hasChildren(j)) then
					write(WRITE_UNIT_1,'(F24.15)') aPanels%absVort(j)
				endif
			enddo
		endif
		if ( associated(aPanels%relVort)) then
			dataString = 'SCALARS  RelVort  double 1'
			write(WRITE_UNIT_1,'(A)') dataString
			write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
			do j=1,aParticles%N
				write(WRITE_UNIT_1,'(F24.15)') aParticles%relVort(j)
			enddo
			do j=1,aPanels%N
				if ( .NOT. aPanels%hasChildren(j)) then
					write(WRITE_UNIT_1,'(F24.15)') aPanels%relVort(j)
				endif
			enddo
		endif

		if ( associated(aPanels%potVort)) then
			dataString = 'SCALARS  PotVort  double 1'
			write(WRITE_UNIT_1,'(A)') dataString
			write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
			do j=1,aParticles%N
				write(WRITE_UNIT_1,'(F24.15)') aParticles%potVort(j)
			enddo
			do j=1,aPanels%N
				if ( .NOT. aPanels%hasChildren(j)) then
					write(WRITE_UNIT_1,'(F24.15)') aPanels%potVort(j)
				endif
			enddo
			allocate(tempArea(aParticles%N))
			tempArea = 0.0_kreal
			do j=1,aPanels%N
				if (.NOT. aPanels%hasChildren(j) ) then
					call CCWEdgesAndParticlesAroundPanel(edgeList,vertList,nVerts,aMesh,j)
					do k=1,nVerts
						tempArea(vertList(k)) = aPanels%area(j)
					enddo
				endif
			enddo
			dataString = 'SCALARS Area double 1'
			write(WRITE_UNIT_1,'(A)') dataString
			write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
			do j=1,aParticles%N
				write(WRITE_UNIT_1,'(F24.5)') tempArea(j)
			enddo
			do j=1,aPanels%N
				if ( .NOT. aPanels%hasChildren(j)) then
					write(WRITE_UNIT_1,'(F24.5)') aPanels%area(j)
				endif
			enddo
			deallocate(tempArea)
		endif

		if ( associated(aPanels%h)) then
			dataString = 'SCALARS  h  double 1'
			write(WRITE_UNIT_1,'(A)') dataString
			write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
			do j=1,aParticles%N
				write(WRITE_UNIT_1,'(F24.15)') aParticles%h(j)
			enddo
			do j=1,aPanels%N
				if ( .NOT. aPanels%hasChildren(j)) then
					write(WRITE_UNIT_1,'(F24.15)') aPanels%h(j)
				endif
			enddo
		endif

		if ( associated(aPanels%div)) then
			dataString = 'SCALARS  div  double 1'
			write(WRITE_UNIT_1,'(A)') dataString
			write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
			do j=1,aParticles%N
				write(WRITE_UNIT_1,'(F24.15)') aParticles%div(j)
			enddo
			do j=1,aPanels%N
				if ( .NOT. aPanels%hasChildren(j)) then
					write(WRITE_UNIT_1,'(F24.15)') aPanels%div(j)
				endif
			enddo
		endif

		if ( associated(aPanels%ke)) then
			dataString = 'SCALARS KE double 1'
			write(WRITE_UNIT_1,'(A)') datastring
			write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
			do j=1,aParticles%N
				write(WRITE_UNIT_1,'(F24.15)') aParticles%ke(j)
			enddo
			do j=1,aPanels%N
				if (.NOT. aPanels%hasChildren(j)) then
					write(WRITE_UNIT_1,'(F24.15)') aPanels%ke(j)
				endif
			enddo
		endif

		nTracer = GetNTracer(aPanels)
		if ( nTracer > 0 ) then
			do k=1,nTracer
				write(dataString,'(A6,I1)') trim('tracer'),k
				write(WRITE_UNIT_1,'(A,A,A)') 'SCALARS  ',trim(dataString),'  double 1'
				write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
				do j=1,aParticles%N
					write(WRITE_UNIT_1,'(F24.15)') aParticles%tracer(j,k)
				enddo
				do j=1,aPanels%N
					if ( .NOT. aPanels%hasChildren(j)) then
						write(WRITE_UNIT_1,'(F24.15)') aPanels%tracer(j,k)
					endif
				enddo
			enddo
		endif
	endif
	if ( associated(aPanels%u) ) then
		dataString = 'SCALARS velocity double 3'
		write(WRITE_UNIT_1,'(A)') dataString
		write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
		do j=1,aParticles%N
			write(WRITE_UNIT_1,'(3F24.8)') aParticles%u(:,j)
		enddo
		do j=1,aPanels%N
			if ( .NOT. aPanels%hasChildren(j) ) then
				write(WRITE_UNIT_1,'(3F24.8)') aPanels%u(:,j)
			endif
		enddo
	endif
	close(WRITE_UNIT_1)
end subroutine

subroutine VTKOutputMidpointRule(self,aMesh)
	type(VTKSource), intent(in) :: self
	type(SphereMesh), intent(in) :: aMesh
	!
	integer(kint) :: writeStat, j, k
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: nVerts

	open(unit=WRITE_UNIT_1, file=self%filename, status='REPLACE', action='WRITE', iostat=writeStat)
	if ( writeStat /= 0 ) then
		call LogMessage(log, ERROR_LOGGING_LEVEL, 'vtkOutput ERROR : ', 'cannot open vtk output file.')
		return
	endif

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	!
	! VTK File header
	!
	write(WRITE_UNIT_1,'(A)') '# vtk DataFile Version 2.0'
	write(WRITE_UNIT_1,'(A)') self%title
	write(WRITE_UNIT_1,'(A)') 'ASCII'
	write(WRITE_UNIT_1,'(A)') 'DATASET POLYDATA'

	!
	! VTK points
	!
	write(WRITE_UNIT_1,'(A,I8,A)') 'POINTS ', aParticles%N, ' double'
	do j = 1, aParticles%N
		write(WRITE_UNIT_1,'(3F24.8)') aParticles%x(1,j)/EARTH_RADIUS, aParticles%x(2,j)/EARTH_RADIUS, aParticles%x(3,j)/EARTH_RADIUS
	enddo

	!
	! VTK Cells
	!
	if ( aMesh%panelKind == TRI_PANEL ) then
		nVerts = 3
	elseif ( aMesh%panelKind == QUAD_PANEL ) then
		nVerts = 4
	endif

	write(WRITE_UNIT_1,'(A, I8, I8)') 'POLYGONS ', aPanels%N_Active, aPanels%N_Active*(nVerts+1)
	if ( aMesh%panelKind == TRI_PANEL ) then
		do j = 1, aPanels%N
			if ( .NOT. aPanels%hasChildren(j) ) then
				write(WRITE_UNIT_1,'(4I8)') 3, aPanels%vertices(:,j) - 1
			endif
		enddo
	elseif ( aMesh%panelKind == QUAD_PANEL ) then
		do j = 1, aPanels%N
			if ( .NOT. aPanels%hasChildren(j) ) then
				write(WRITE_UNIT_1,'(5I8)') 4, aPanels%vertices(:,j) - 1
			endif
		enddo
	endif


	!
	! VTK Scalar Data
	!
	write(WRITE_UNIT_1,'(A, I8)') 'POINT_DATA ', aParticles%N
	write(WRITE_UNIT_1,'(A)') 'SCALARS lagParam double 3'
	write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
	do j = 1, aparticles%N
		write(WRITE_UNIT_1,'(3F24.8)') aParticles%x0(:,j)/EARTH_RADIUS
	enddo
	write(WRITE_UNIT_1,'(A)') 'SCALARS velocity double 3'
	write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
		do j = 1, aparticles%N
		write(WRITE_UNIT_1,'(3F24.8)') aParticles%u(:,j)
	enddo
	if ( associated(aParticles%absVort) ) then
		write(WRITE_UNIT_1,'(A)') 'SCALARS absVort double 1'
		write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
		do j = 1, aParticles%N
			write(WRITE_UNIT_1, '(F24.8)') aParticles%absVort(j)
		enddo
	endif
	if ( associated(aParticles%relVort) ) then
		write(WRITE_UNIT_1,'(A)') 'SCALARS relVort double 1'
		write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
		do j = 1, aParticles%N
			write(WRITE_UNIT_1, '(F24.8)') aParticles%relVort(j)
		enddo
	endif
	if ( associated(aParticles%potVort) ) then
		write(WRITE_UNIT_1,'(A)') 'SCALARS potVort double 1'
		write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
		do j = 1, aParticles%N
			write(WRITE_UNIT_1, '(F24.8)') aParticles%potVort(j)
		enddo
	endif
	if ( associated(aPanels%h)) then
		write(WRITE_UNIT_1,'(A)') 'SCALARS  h  double 1'
		write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
		do j=1,aParticles%N
			write(WRITE_UNIT_1,'(F24.15)') aParticles%h(j)
		enddo
	endif
	if ( associated(aPanels%div)) then
		write(WRITE_UNIT_1,'(A)') 'SCALARS  div  double 1'
		write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
		do j=1,aParticles%N
			write(WRITE_UNIT_1,'(F24.8)') aParticles%div(j)
		enddo
	endif
	if ( aMesh%nTracer > 0 ) then
		do k = 1, aMesh%nTracer
			write(WRITE_UNIT_1,'(A,I1,A)') 'SCALARS tracer', k, '  double 1'
			write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
			do j = 1, aParticles%N
				write(WRITE_UNIT_1,'(F24.8)') aParticles%tracer(j,k)
			enddo
		enddo
	endif


	!
	! VTK Cell data
	!
	write(WRITE_UNIT_1,'(A,I8)') 'CELL_DATA ', aPanels%N_Active
	write(WRITE_UNIT_1,'(A)') 'SCALARS lagParamPanel double 3'
	write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
	do j = 1, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			write(WRITE_UNIT_1,'(3F24.8)') aPanels%x0(:,j)/EARTH_RADIUS
		endif
	enddo
	write(WRITE_UNIT_1,'(A)') 'SCALARS velocityPanel double 3'
	write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
	do j = 1, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			write(WRITE_UNIT_1,'(3F24.8)') aPanels%u(:,j)
		endif
	enddo
	if ( associated(aPanels%absVort) ) then
		write(WRITE_UNIT_1,'(A)') 'SCALARS absVortPanel double 1'
		write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
		do j = 1, aPanels%N
			if ( .NOT. aPanels%hasChildren(j) ) then
				write(WRITE_UNIT_1,'(3F24.8)') aPanels%absVort(j)
			endif
		enddo
	endif
	if ( associated(aPanels%relVort) ) then
		write(WRITE_UNIT_1,'(A)') 'SCALARS relVortPanel double 1'
		write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
		do j = 1, aPanels%N
			if ( .NOT. aPanels%hasChildren(j) ) then
				write(WRITE_UNIT_1,'(3F24.8)') aPanels%relVort(j)
			endif
		enddo
	endif
	if ( associated(aPanels%potVort) ) then
		write(WRITE_UNIT_1,'(A)') 'SCALARS potVortPanel double 1'
		write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
		do j = 1, aPanels%N
			if ( .NOT. aPanels%hasChildren(j) ) then
				write(WRITE_UNIT_1,'(3F24.8)') aPanels%potVort(j)
			endif
		enddo
	endif
	if ( associated(aPanels%h) ) then
		write(WRITE_UNIT_1,'(A)') 'SCALARS hPanel double 1'
		write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
		do j = 1, aPanels%N
			if ( .NOT. aPanels%hasChildren(j) ) then
				write(WRITE_UNIT_1,'(3F24.8)') aPanels%h(j)
			endif
		enddo
	endif
	if ( associated(aPanels%div) ) then
		write(WRITE_UNIT_1,'(A)') 'SCALARS divPanel double 1'
		write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
		do j = 1, aPanels%N
			if ( .NOT. aPanels%hasChildren(j) ) then
				write(WRITE_UNIT_1,'(3F24.8)') aPanels%div(j)
			endif
		enddo
	endif
	if ( aMesh%nTracer > 0 ) then
		do k = 1, aMesh%nTracer
			write(WRITE_UNIT_1,'(A,I1,A)') 'SCALARS tracerPanel', k, ' double 1'
			write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
			do j = 1, aPanels%N
				if ( .NOT. aPanels%hasChildren(j) ) then
					write(WRITE_UNIT_1,'(F24.8)') aPanels%tracer(j,k)
				endif
			enddo
		enddo
	endif

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
