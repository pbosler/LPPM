module PlaneOutputModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	Defines the RK4 data structure used by SphereMesh.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
!----------------
! USAGE :  This module provides methods for integrating the barotropic vorticity equation on the sphere.
!----------------
use NumberKindsModule
use LoggerModule
use PlaneGeomModule
use ParticlesModule
use EdgesModule
use PanelsModule
use PlaneMeshModule

implicit none

include 'mpif.h'

private
public PlaneOutput
public New, Delete
public OutputForVTK, OutputForNCL, OutputForMatlab
public UpdateFilename

!
!----------------
! Types and module constants
!----------------
!
type PlaneOutput
	character(len=MAX_STRING_LENGTH) :: filename
	character(len=128) :: title
	integer(kint) :: vtkNPoints, &
					 vtknCells, &
					 vtkCellListSize
end type
!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'PlaneOutput'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=28) :: formatString
character(len=128) :: logString
!
!----------------
! Interfaces
!----------------
!
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

subroutine NewPrivate(self, aMesh, filename, title)
	type(PlaneOutput), intent(out) :: self
	type(PlaneMesh), intent(in) :: aMesh
	character(len=*), intent(in) :: filename
	character(len=*), intent(in), optional :: title
	!
	type(Particles), pointer :: aparticles
	type(Edges), pointer :: anedges
	type(Panels), pointer :: apanels
	integer(kint) :: j
	integer(kint) :: edgelist(8), vertlist(8), nverts

	if ( .NOT. loginit) call InitLogger(log,procRank)

	aparticles => aMesh%particles
	anedges => aMesh%edges
	apanels => aMesh%panels

	self%filename = trim(filename)
	if ( present(title) ) then
		self%title = trim(title)
	else
		self%title = ' '
	endif

	self%vtkNpoints = aParticles%N + aPanels%N_Active
	self%vtkNCells = 0
	self%vtkcelllistsize = 0

	do j=1,aPanels%N
		if (.not. aPanels%hasChildren(j) ) then
			call CCWEdgesAndParticlesAroundPanel(edgelist, vertlist, nverts, aMesh, j)
			self%vtkncells = self%vtkncells + nverts
			self%vtkcelllistsize = self%vtkcelllistsize + 4*nVerts
		endif
	enddo
end subroutine

subroutine DeletePrivate(self)
	type(PlaneOutput), intent(inout) :: self
	self%vtknpoints = 0
	self%vtkncells = 0
	self%vtkcelllistsize = 0
	self%filename = ' '
	self%title = ' '
end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine UpdateFilename(self, newfilename)
	type(PlaneOutput), intent(inout) :: self
	character(len=*), intent(in) :: newfilename
	self%filename = trim(newfilename)
end subroutine

subroutine OutputForVTK(self, aMesh)
	type(PlaneOutput), intent(in) :: self
	type(PlaneMesh), intent(in) :: aMesh
	!
	integer(kint) :: writeStat, j, k
	type(Particles), pointer :: aparticles
	type(Edges), pointer :: anedges
	type(Panels), pointer :: apanels
	integer(kint) :: edgelist(8), vertlist(8), nverts, panelk
	character(len=28) :: datastring

	open(unit = WRITE_UNIT_1,file=self%filename,status='REPLACE',action='WRITE',iostat=writeStat)
	if ( writeStat /= 0 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'VTKOutput ERROR opening output file.')
		return
	endif

	aParticles=>aMesh%particles
	anEdges=>aMesh%edges
	aPanels=>aMesh%panels
	vertList = -1

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
	write(WRITE_UNIT_1,'(A,I8,A)') 'POINTS ',self%vtknPoints,' double'
	do j=1,aparticles%N
		write(WRITE_UNIT_1,'(3F24.8)') aParticles%x(1,j),aParticles%x(2,j), 0.0_kreal
	enddo
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			write(WRITE_UNIT_1,'(3F24.8)') aPanels%x(1,j), aPanels%x(2,j), 0.0_kreal
		endif
	enddo
	!
	! VTK cells
	!
	write(WRITE_UNIT_1,'(A,I8,I8)') 'POLYGONS ',self%vtknCells, self%vtkcellListSize
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
	!
	! VTK scalars
	!
	write(WRITE_UNIT_1,'(A,I8)') 'POINT_DATA  ',self%vtknPoints
	datastring = 'SCALARS lagParam double 3'
	write(WRITE_UNIT_1,'(A)') trim(dataString)
	write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
	do j=1,aParticles%N
		write(WRITE_UNIT_1,'(3F24.8)') aParticles%x0(1,j), aParticles%x0(2,j), 0.0_kreal
	enddo
	do j=1,aPanels%N
		if (.NOT. aPanels%hasChildren(j) ) then
			write(WRITE_UNIT_1,'(3F24.8)') aPanels%x0(1,j), aPanels%x0(2,j), 0.0_kreal
		endif
	enddo

	dataString = 'SCALARS velocity double 3'
	write(WRITE_UNIT_1,'(A)') dataString
	write(WRITE_UNIT_1,'(A)') 'LOOKUP_TABLE default'
	do j=1,aParticles%N
		write(WRITE_UNIT_1,'(3F24.8)') aParticles%u(1,j), aParticles%u(2,j), 0.0_kreal
	enddo
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			write(WRITE_UNIT_1,'(3F24.8)') aPanels%u(1,j), aPanels%u(2,j), 0.0_kreal
		endif
	enddo

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

	do k=1, amesh%nTracer
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
	close(WRITE_UNIT_1)
end subroutine

subroutine OutputForNCL(self, aMesh)
	type(PlaneOutput), intent(in) :: self
	type(PlaneMesh), intent(in) :: aMesh
	! TO DO
end subroutine

subroutine OutputForMatlab(self, aMesh)
	type(PlaneOutput), intent(in) :: self
	type(PlaneMesh), intent(in) :: aMesh
	!
	integer(kint) :: writeStat, j
	type(Particles), pointer :: aparticles
	type(Edges), pointer :: anedges
	type(Panels), pointer :: apanels

	open(unit = WRITE_UNIT_1,file=self%filename,status='REPLACE',action='WRITE',iostat=writeStat)
	if ( writeStat /= 0 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'Matlab Output ERROR opening output file.')
		return
	endif

	aParticles=>aMesh%particles
	anEdges=>aMesh%edges
	aPanels=>aMesh%panels

	write(WRITE_UNIT_1,'(A,F24.12, A, F24.12, A)') 'x = [ ', aParticles%x(1,1),' ; ', aParticles%x(2,1), ' ,... '
	do j=2,aParticles%N
		write(WRITE_UNIT_1,'(F24.12, A, F24.12, A)') aParticles%x(1,j), ' ; ', aParticles%x(2,j), ' ,... '
	enddo
	do j=1,aPanels%N-1
		if ( .NOT. aPanels%hasChildren(j) ) then
			write(WRITE_UNIT_1,'(F24.12, A, F24.12, A)') aPanels%x(1,j), ' ; ', aPanels%x(2,j), ' ,...'
		endif
	enddo
	write(WRITE_UNIT_1,'(F24.12, A, F24.12, A)') aPanels%x(1,aPanels%N), ' ; ', aPanels%x(2,aPanels%N), ' ];'

	write(WRITE_UNIT_1,'(A)') ' '

	write(WRITE_UNIT_1,'(A,F24.12, A, F24.12, A)') 'x0 = [ ', aParticles%x0(1,1),' ; ', aParticles%x0(2,1), ' ,... '
	do j=2,aParticles%N
		write(WRITE_UNIT_1,'(F24.12, A, F24.12, A)') aParticles%x0(1,j), ' ; ', aParticles%x0(2,j), ' ,... '
	enddo
	do j=1,aPanels%N-1
		if ( .NOT. aPanels%hasChildren(j) ) then
			write(WRITE_UNIT_1,'(F24.12, A, F24.12, A)') aPanels%x0(1,j), ' ; ', aPanels%x0(2,j), ' ,...'
		endif
	enddo
	write(WRITE_UNIT_1,'(F24.12, A, F24.12, A)') aPanels%x0(1,aPanels%N), ' ; ', aPanels%x0(2,aPanels%N), ' ];'

	write(WRITE_UNIT_1,'(A,F24.12, A, F24.12, A)') 'u = [ ', aParticles%u(1,1),' ; ', aParticles%u(2,1), ' ,... '
	do j=2,aParticles%N
		write(WRITE_UNIT_1,'(F24.12, A, F24.12, A)') aParticles%u(1,j), ' ; ', aParticles%u(2,j), ' ,... '
	enddo
	do j=1,aPanels%N-1
		if ( .NOT. aPanels%hasChildren(j) ) then
			write(WRITE_UNIT_1,'(F24.12, A, F24.12, A)') aPanels%u(1,j), ' ; ', aPanels%u(2,j), ' ,...'
		endif
	enddo
	write(WRITE_UNIT_1,'(F24.12, A, F24.12, A)') aPanels%u(1,aPanels%N), ' ; ', aPanels%u(2,aPanels%N), ' ];'

	write(WRITE_UNIT_1,'(A,F24.12, A)') 'vort = [ ', aParticles%relVort(1),' ,... '
	do j=2,aParticles%N
		write(WRITE_UNIT_1,'(F24.12, A )') aParticles%relVort(j), ' ,... '
	enddo
	do j=1,aPanels%N-1
		if ( .NOT. aPanels%hasChildren(j) ) then
			write(WRITE_UNIT_1,'(F24.12, A)') aPanels%relVort(j), ' ,...'
		endif
	enddo
	write(WRITE_UNIT_1,'(F24.12, A)') aPanels%relVort(aPanels%N),  ' ];'


	close(WRITE_UNIT_1)
end subroutine

!
!----------------
! Module methods : module- or type-specific private functions
!----------------
!

subroutine InitLogger(aLog, rank)
	type(Logger), intent(inout) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	if ( rank == 0 ) then
		call New(aLog,logLevel)
	else
		call New(aLog,WARNING_LOGGING_LEVEL)
	endif
	logInit = .TRUE.
end subroutine

end module
