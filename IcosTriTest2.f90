program IcosTrianglesTest
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	This program is a unit test for a basic cubed sphere sphere, no time integration.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!

!----------------
use NumberKindsModule
use LoggerModule
use SphereMeshModule
use ParticlesModule
use EdgesModule
use PanelsModule
use VTKOutputModule
use AdvectionModule

implicit none

type(SphereMesh) :: sphere
integer(kint) :: panelKind, initNest, AMR, nTracer, problemKind
type(Particles), pointer :: sphereParticles
type(Edges), pointer :: sphereEdges
type(Panels), pointer :: spherePanels

type(VTKSource) :: vtkOut
character(len=56) :: testVTKfile

type(Logger) :: exeLog
character(len=56) :: logString

integer(kint) :: j

!
!	set mesh parameters
!
	initNest = 5
	panelKind = 3
	AMR = 0
	nTracer = 1
	problemKind = ADVECTION_SOLVER
	call New(exeLog,DEBUG_LOGGING_LEVEL)

!
!	Build a mesh
!
	call New(sphere,panelKind,initNest,AMR,nTracer,problemKind)
	sphereParticles => sphere%particles
	sphereEdges => sphere%edges
	spherePanels => sphere%panels
	write(logString,'(A,I4,A)') 'Mesh at initNest = ',initNest,' returned : '
	call LogStats(sphere,exeLog,trim(logString))
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,'Surf. area should be = ',4.0_kreal*PI*EARTH_RADIUS*EARTH_RADIUS)
	write(6,'(A,F24.15)') "surf area error = ",abs(4.0_kreal*PI*EARTH_RADIUS*EARTH_RADIUS - sum(spherePanels%area))

!
!	Define VTKOutput
!
	write(testVTKFile,'(A,I1,A)') 'icosTriMesh',initNest,'.vtk'
	call New(vtkOut,sphere,testVTKFile,'test1')
	call vtkOutput(vtkOut,sphere)
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,'VTKOutput : ',' done.')

!
!	clean up
!
	call Delete(vtkOut)
	call Delete(sphere)
	call Delete(exeLog)

end program
