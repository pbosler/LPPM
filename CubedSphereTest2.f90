program CubedSphereTest

use NumberKindsModule
use LoggerModule
use SphereMeshModule
use ParticlesModule
use EdgesModule
use PanelsModule
use SphereMeshOutputModule

implicit none

type(SphereMesh), pointer :: grid
integer(kint) :: panelKind, initNest, AMR, nTracer, problemKind
type(Particles), pointer :: gridParticles
type(Edges), pointer :: gridEdges
type(Panels), pointer :: gridPanels

type(VTKMetaData) :: vtkOut
character(len=56) :: testVTKfile
type(LLContourData) :: nclOut
character(len=56) :: testNCLFile

type(Logger) :: exeLog
character(len=56) :: logString

integer(kint) :: i, j, inpanel
character(len=56) :: tracerFile

real(kreal) :: tracerLL(181,360), lats(181), lons(360), xyz(3)

initNest = 4
write(testVTKFile,'(A,I1,A)') 'EvenMoreBetterPlotVTKFile',initNest,'.vtk'
write(testNCLFile,'(A,I1,A)') 'EvenMoreBetterPlotNCLFile',initNest,'.dat'
panelKind = 4
AMR = 0
nTracer = 1
!problemKind = BVE_SOLVER
problemKind = ADVECTION_SOLVER

call New(exeLog,DEBUG_LOGGING_LEVEL)


	allocate(grid)

	call New(grid,panelKind,initNest,AMR,nTracer,problemKind)
	gridParticles => grid%particles
	gridEdges => grid%edges
	gridPanels => grid%panels

!
	write(logString,'(A,I4,A)') 'Mesh at initNest = ',initNest,' returned : '
	call LogStats(grid,exeLog,trim(logString))
!
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,'Surf. area should be = ',4.0_kreal*PI)

	write(6,'(A,F24.15)') "surf area error = ",abs(4.0_kreal*PI - sum(gridPanels%area))

	do j=1,gridPanels%N
		gridPanels%tracer(j,1) = real(j,kreal)
	enddo
	do j=1,gridParticles%N
		inPanel = LocatePoint(grid,gridParticles%x(:,j))
		gridParticles%tracer(j,1) = real(inPanel,kreal)
	enddo

	do i=1,181
		do j=1,360
			xyz = [ cos(-PI/2.0_kreal + (i-1)*PI/180.0_kreal)*cos( (j-1)*PI/180.0_kreal),&
					cos(-PI/2.0_kreal + (i-1)*PI/180.0_kreal)*sin( (j-1)*PI/180.0_kreal),&
					sin(-PI/2.0_kreal + (i-1)*PI/180.0_kreal ) ]
			tracerLL(i,j) = real(LocatePoint(grid,xyz),kreal)
		enddo
	enddo

	tracerFile = 'panelIDTracer.dat'
	open(unit=WRITE_UNIT_2,file=tracerFile,status='REPLACE',action='WRITE')
		do i=1,181
			do j=1,360
				write(WRITE_UNIT_2,'(F16.3)',advance = 'NO') tracerLL(i,j)
			enddo
			write(WRITE_UNIT_2,'(A)',advance='YES') ' '
		enddo
	close(WRITE_UNIT_2)

!	write(6,'(A12,A24)') 'panelIndex','area'
!	do j=1,gridPanels%N
!		write(6,'(I12,F24.15)') j, gridPanels%area(j)
!	enddo
!
!
!
!	write(6,'(5A8)') 'edge #','vert1','vert2','left','right'
!	do j=1,gridEdges%N
!		write(6,'(5I8)') j,gridEdges%verts(:,j),gridEdges%leftpanel(j), gridEdges%rightPanel(j)
!	enddo
!	write(6,'(A)') ' '
!
!	write(6,'(9A8)') 'panel #','edge1','edge2','edge3','edge4','vert1','vert2','vert3','vert4'
!	do j=1,gridPanels%N
!		write(6,'(9I8)') j,gridPanels%edges(:,j),gridPanels%vertices(:,j)
!	enddo


	call New(vtkOut,grid,testVTKFile,'test1')
	call vtkOutput(vtkOut,grid)
	call Delete(vtkOut)

	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,'VTKOutput : ',' done.')

	call New(nclOut,testnclFile,360)
	call MatlabLLOutputCubicInterp(nclOut,grid)
	call Delete(nclOut)

	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,'MatlabOutput : ',' done.')

	call Delete(grid)

	deallocate(grid)


call Delete(exeLog)

end program
