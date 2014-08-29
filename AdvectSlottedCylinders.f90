program SlottedCylindersAdvection

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use SphereMeshModule
use AdvectionModule
use ParticlesModule
use PanelsModule
use SphereMeshModule
use TracerSetupModule
use VTKOutputModule
use BVESetupModule
use SphereRemeshModule

implicit none

include 'mpif.h'

!
! mesh variables
!
type(SphereMesh) :: sphere
integer(kint) :: panelKind, initNest, AMR, nTracer
type(Particles), pointer :: sphereParticles
type(Panels), pointer :: spherePanels

!
! tracer variables
!
type(TracerSetup) :: slotC
integer(kint) :: tracerID

!
! vorticity placeholder
!
type(BVESetup) :: nullVort

!
! remeshing / refinement variables
!
type(RemeshSetup) :: remesh
integer(kint) :: remeshInterval, resetAlphaInterval, amrLimit, remeshCounter
type(ReferenceSphere), pointer :: reference

!
! timestepping variables
!
type(AdvRK4Data) :: timekeeper
real(kreal) :: t, tfinal, dt
integer(kint) :: timesteps, timeJ

!
! output variables
!


end program
