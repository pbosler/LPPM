module VTKOutputModule
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


implicit none
private
public VTKSource
public New, Delete
public VTKOutput, UpdateFilename

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

subroutine vtkOutput(self,aMesh)
	type(VTKSource), intent(in) :: self
	type(SphereMesh), intent(in) :: aMesh
	integer(kint) :: writeStat, j, k
	type(Particles), pointer :: aParticles
	type(Edges), pointer :: anEdges
	type(Panels), pointer :: aPanels
	integer(kint) :: edgeList(8), vertList(8), nVerts, panelKind, panelK, nTracer
	character(len=28) :: dataString

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
		write(WRITE_UNIT_1,'(3F24.8)') aParticles%x(1,j),aParticles%x(2,j),aParticles%x(3,j)
	enddo
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			write(WRITE_UNIT_1,'(3F24.8)') aPanels%x(1,j), aPanels%x(2,j), aPanels%x(3,j)
		endif
	enddo

! TO DO : AMR / Variable size polygons

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
			write(WRITE_UNIT_1,'(3F24.8)') aParticles%x0(:,j)
		enddo
		do j=1,aPanels%N
			if (.NOT. aPanels%hasChildren(j) ) then
				write(WRITE_UNIT_1,'(3F24.8)') aPanels%x0(:,j)
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
		if ( associated(aPanels%absVort)) then
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

		nTracer = GetNTracer(aPanels)
		if ( nTracer > 0 ) then
			do k=1,nTracer
				write(dataString,'(A6,I1)') trim('Tracer'),k
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
