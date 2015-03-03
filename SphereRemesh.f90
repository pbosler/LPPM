module SphereRemeshModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method - Spherical Model
!------------------------------------------------------------------------------
!
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
!> @defgroup SphereRemesh Sphere Remesh
!> Data types, subroutines, and functions for adaptive remeshing spherical particle/panel meshes.
!
!
! DESCRIPTION:
!> @file
!> Data types, subroutines, and functions for adaptive remeshing spherical particle/panel meshes.
!
!------------------------------------------------------------------------------

use NumberKindsModule
use LoggerModule
use SphereGeomModule
use ParticlesModule
use PanelsModule
use SphereMeshModule
use BVESetupModule
use TracerSetupModule
use STRIPACKInterfaceModule
use SSRFPACKInterfaceModule

implicit none

include 'mpif.h'

private
public RemeshSetup, ReferenceSphere
public New, Delete, SetReferenceValues
public LagrangianRemeshToInitialTime, LagrangianRemeshToReference
public InitialRefinement
public ResetLagrangianParameter
public DirectRemesh

!
!----------------
! Types and module constants
!----------------
!
type RemeshSetup
	logical(klog) :: useAMR ! if .FALSE., uniform meshes will be used
	logical(klog) :: vorticityRefine
	real(kreal) :: maxCircTol
	real(kreal) :: vortVarTol
	logical(klog) :: flowMapRefine
	real(kreal) :: lagVarTol
	logical(klog) :: tracerRefine
	integer(kint) :: tracerID
	real(kreal) :: tracerMassTol
	real(kreal) :: tracerVarTol
	logical(klog) :: useReferenceVal = .FALSE.
	real(kreal) :: refVal
	real(kreal) :: refTol
	integer(kint) :: refinementLimit
end type

type ReferenceSphere
	type(STRIPACKData) :: delTri
	type(SSRFPACKData), pointer :: absVortSource => null()
	type(SSRFPACKData), pointer :: tracerSource => null()
end type


!
!----------------
! Interfaces
!----------------
!
interface New
	module procedure NewPrivateAll
	module procedure NewPrivateVorticity
	module procedure NewPrivateVorticityAndFlowMap
	module procedure NewPrivateTracer
	module procedure NewPrivateTracerAndFlowMap
	module procedure NewPrivateNoAMR
	module procedure NewReference
end interface

interface Delete
	module procedure DeletePrivate
	module procedure DeleteReference
end interface

interface InitialRefinement
	module procedure InitialRefinementPrivate
end interface

interface LagrangianRemeshToInitialTime
	module procedure LagrangianRemeshToInitialTimePrivate
end interface

interface
	subroutine SetTracerOnMesh(genMesh, genTracer)
		use SphereMeshModule
		use TracerSetupModule
		implicit none
		type(SphereMesh), intent(inout) :: genMesh
		type(TracerSetup), intent(in) :: genTracer
	end subroutine
end interface

interface
	subroutine SetVorticityOnMesh(genMesh,genVort, t)
		use SphereMeshModule
		use BVESetupModule
		implicit none
		type(SphereMesh), intent(inout) :: genMesh
		type(BVESetup), intent(in) :: genVort
		real(8), intent(in), optional :: t
	end subroutine
end interface

!
!----------------
! Logging
!----------------
!

logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'SphereRemesh'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=128) :: logString
character(len=24) :: formatString

contains
!
!----------------
! Standard methods : Constructor / Destructor
!----------------
!

subroutine NewPrivateAll(self, maxCircTol, vortVarTol, lagVarTol, tracerID, tracerMassTol, tracerVarTol, limit)
	type(RemeshSetup), intent(out) :: self
	real(kreal), intent(in) :: maxCircTol, vortVarTol, lagVarTol
	integer(kint), intent(in) :: tracerID
	real(kreal), intent(in) :: tracerMassTol, tracerVarTol
	integer(kint), intent(in) :: limit

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	if 	(maxCircTol > 0.0_kreal .OR. vortVarTol > 0.0_kreal) then
		self%vorticityRefine = .TRUE.
		self%maxCircTol = maxCircTol
		self%vortVarTol = vortVarTol
	else
		self%vorticityRefine = .FALSE.
		self%maxCircTol = 0.0_kreal
		self%vortVarTol = 0.0_kreal
	endif
	if ( lagVarTol > 0.0_kreal ) then
		self%flowmapRefine = .TRUE.
		self%lagVarTol = lagVarTol
	else
		self%flowMapREfine = .FALSE.
		self%lagVarTol = 0.0_kreal
	endif
	if ( tracerID > 0 .AND. ( tracerMassTol > 0.0_kreal .OR. tracerVarTol > 0.0_kreal ) ) then
		self%tracerRefine = .TRUE.
		self%tracerID = tracerID
		self%tracerMassTol = tracerMassTol
		self%tracerVarTol = tracerVarTol
	else
		self%tracerRefine = .FALSE.
		self%tracerID = 0
		self%tracerMassTol = 0.0_kreal
		self%tracerVarTol = 0.0_kreal
	endif
	if ( limit > 0 ) then
		if ( (self%vorticityRefine .OR. self%flowMapRefine) .OR. self%tracerRefine ) then
			self%useAMR = .TRUE.
			self%refinementLimit = limit
		else
			self%useAMR = .FALSE.
			self%refinementLimit = 0
		endif
	else
		self%useAMR = .FALSE.
		self%vorticityRefine = .FALSE.
		self%flowmapRefine = .FALSE.
		self%tracerRefine = .FALSE.
		self%refinementLimit = 0
	endif
	self%useReferenceVal = .FALSE.
	self%refTol = 1.0e20
	self%refVal = 0.0
end subroutine

subroutine NewPrivateVorticity(self, maxCircTol, vortVarTol, limit)
	type(RemeshSetup), intent(out) :: self
	real(kreal), intent(in) :: maxCircTol, vortVarTol
	integer(kint), intent(in) :: limit

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	if ( maxCircTol > 0.0_kreal .OR. vortVarTol > 0.0_kreal) then
		if ( limit > 0 ) then
			self%useAMR = .TRUE.
			self%maxCircTol = maxCircTol
			self%vortVarTol = vortVarTol
			self%refinementLimit = limit
		else
			self%useAMR = .FALSE.
			self%refinementLimit = 0
			self%maxCircTol = 0.0_kreal
			self%vortVarTol = 0.0_kreal
		endif
	else
		self%useAMR = .FALSE.
		self%refinementLimit = 0
		self%maxCircTol = 0.0_kreal
		self%vortVarTol = 0.0_kreal
	endif
	self%flowMapRefine = .FALSE.
	self%lagVarTol = 0.0_kreal
	self%tracerRefine = .FALSE.
	self%tracerID = 0
	self%tracerMassTol = 0.0_kreal
	self%tracerVarTol = 0.0_kreal
	self%useReferenceVal = .FALSE.
	self%refTol = 1.0e20
	self%refVal = 0.0
end subroutine

subroutine NewPrivateVorticityAndFlowMap(self, maxCircTol, vortVarTol, lagVarTol, limit)
	type(RemeshSetup), intent(out) :: self
	real(kreal), intent(in) :: maxCircTol, vortVarTol, lagVarTol
	integer(kint), intent(in) :: limit

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	if ( maxCircTol > 0.0_kreal .OR. vortVarTol > 0.0_kreal ) then
		self%vorticityRefine = .TRUE.
		self%maxCircTol = maxCircTol
		self%vortVarTol = vortVarTol
	else
		self%vorticityREfine = .FALSE.
		self%maxCircTol = 0.0_kreal
		self%vortVarTol = 0.0_kreal
	endif
	if ( lagVarTol > 0.0_kreal ) then
		self%flowMapRefine = .TRUE.
		self%lagVarTol = lagVarTol
	else
		self%flowMapRefine = .FALSE.
		self%lagVarTol = 0.0_kreal
	endif
	if ( limit > 0 ) then
		if ( self%vorticityRefine .OR. self%flowMapRefine ) then
			self%useAMR = .TRUE.
			self%refinementLimit = limit
		else
			self%useAMR = .FALSE.
			self%refinementLimit = 0
		endif
	else
		self%useAMR = .FALSE.
		self%refinementLimit = 0
	endif
	self%tracerRefine = .FALSE.
	self%tracerID = 0
	self%tracerMassTol = 0.0_kreal
	self%tracerVarTol = 0.0_kreal
	self%useReferenceVal = .FALSE.
	self%refTol = 1.0e20
	self%refVal = 0.0
end subroutine

subroutine NewPrivateTracer(self, tracerID, tracerMassTol, tracerVarTol, limit)
	type(RemeshSetup), intent(out) :: self
	integer(kint), intent(in) :: tracerID
	real(kreal), intent(in) :: tracerMassTol, tracerVarTol
	integer(kint), intent(in) :: limit

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	if ( tracerID > 0 ) then
		if ( tracerMassTol > 0.0_kreal .OR. tracerVarTol > 0.0_kreal) then
			self%tracerRefine = .TRUE.
			self%tracerID = tracerID
			self%tracerMassTol = tracerMassTol
			self%tracerVarTol = tracerVarTol
		else
			self%tracerRefine = .FALSE.
			self%tracerID = 0
			self%tracerMassTol = 0.0_kreal
			self%tracerVarTol = 0.0_kreal
		endif
	else
		self%tracerRefine = .FALSE.
		self%tracerID = 0
		self%tracerMassTol = 0.0_kreal
		self%tracerVarTol = 0.0_kreal
	endif
	if ( self%tracerRefine .AND. limit > 0 ) then
		self%useAMR = .TRUE.
		self%refinementLimit = limit
	else
		self%useAMR = .FALSE.
		self%tracerRefine = .FALSE.
		self%refinementLimit = 0
	endif
	self%vorticityRefine = .FALSE.
	self%maxCircTol = 0.0_kreal
	self%vortVarTol = 0.0_kreal
	self%flowMapRefine = .FALSE.
	self%lagVarTol = 0.0_kreal
	self%useReferenceVal = .FALSE.
	self%refTol = 1.0e20
	self%refVal = 0.0
end subroutine

subroutine NewPrivateTracerAndFlowMap(self, tracerID, tracerMassTol, tracerVarTol, lagVarTol, limit)
	type(RemeshSetup), intent(out) :: self
	integer(kint), intent(in) :: tracerID
	real(kreal), intent(in) :: tracerMassTol, tracerVarTol, lagVarTol
	integer(kint), intent(in) :: limit

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	if ( tracerID > 0 ) then
		if ( tracerMassTol > 0.0_kreal .OR. tracerVarTol > 0.0_kreal ) then
			self%tracerRefine = .TRUE.
			self%tracerID = tracerId
			self%tracerMassTol = tracerMassTol
			self%tracerVarTol = tracerVarTol
		else
			self%tracerRefine = .FALSE.
			self%tracerID = 0
			self%tracerMassTol = 0.0_kreal
			self%tracerVarTol = 0.0_kreal
		endif
	else
		self%tracerRefine = .FALSE.
		self%tracerID = 0
		self%tracerMassTol = 0.0_kreal
		self%tracerVarTol = 0.0_kreal
	endif
	if ( lagVarTol > 0.0_kreal ) then
		self%flowMapRefine = .TRUE.
		self%lagVarTol = lagVarTol
	else
		self%flowMapRefine = .FALSE.
		self%lagVarTol = 0.0_kreal
	endif
	if ( (self%tracerRefine .OR. self%flowMapRefine ) .AND. limit > 0 ) then
		self%useAMR = .TRUE.
		self%refinementLimit = limit
	else
		self%useAMR = .FALSE.
		self%tracerRefine = .FALSE.
		self%flowMapRefine = .FALSE.
		self%refinementLimit = 0
	endif
	self%vorticityRefine = .FALSE.
	self%maxCircTol = 0.0_kreal
	self%vortVarTol = 0.0_kreal
	self%useReferenceVal = .FALSE.
	self%refTol = 1.0e20
	self%refVal = 0.0
end subroutine

subroutine NewPrivateNoAMR(self)
	type(RemeshSetup), intent(out) :: self

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	self%useAMR = .FALSE.
	self%vorticityRefine = .FALSE.
	self%maxCircTol = 0.0_kreal
	self%vortVarTol = 0.0_kreal
	self%flowMapRefine = .FALSE.
	self%lagVarTol = 0.0_kreal
	self%tracerRefine = .FALSE.
	self%tracerID = 0
	self%tracerMassTol = 0.0_kreal
	self%tracerVarTol = 0.0_kreal
	self%refinementLimit = 0
	self%useReferenceVal = .FALSE.
	self%refTol = 1.0e20
	self%refVal = 0.0
end subroutine

subroutine DeletePrivate(self)
	type(RemeshSetup), intent(inout) :: self
	self%useAMR = .FALSE.
	self%vorticityRefine = .FALSE.
	self%maxCircTol = 0.0_kreal
	self%vortVarTol = 0.0_kreal
	self%flowMapRefine = .FALSE.
	self%lagVarTol = 0.0_kreal
	self%tracerRefine = .FALSE.
	self%tracerID = 0
	self%tracerMassTol = 0.0_kreal
	self%tracerVarTol = 0.0_kreal
	self%refinementLimit = 0
	self%useReferenceVal = .FALSE.
	self%refTol = 1.0e20
	self%refVal = 0.0
end subroutine

subroutine NewReference(self, aMesh)
	type(ReferenceSphere), intent(out) :: self
	type(SphereMesh), intent(in) :: aMesh
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	call New(self%delTri, aMesh)

	aParticles => amesh%particles
	aPanels => aMesh%panels

	if ( aMesh%problemKind == BVE_SOLVER ) then
		allocate(self%absVortSource)
		call New(self%absVortSource, self%delTri, .FALSE.)
		call SetSourceAbsVort(self%absVortSource, self%delTri)
	endif

	if ( aMesh%nTracer > 0 ) then
		allocate(self%tracerSource)
		call New(self%tracerSource, self%delTri, amesh%nTracer)
		call SetSourceTracer(self%tracerSource, self%delTri)
	endif
end subroutine

subroutine DeleteReference(self)
	type(ReferenceSphere), intent(inout) :: self
	if ( associated(self%tracerSource) ) then
		call Delete(self%tracerSource)
		deallocate(self%tracerSource)
	endif
	if ( associated(self%absVortSource) ) then
		call Delete(self%absVortSource)
		deallocate(self%absVortSource)
	endif
	call Delete(self%delTri)
end subroutine


!
!----------------
! Public functions
!----------------
!
subroutine SetReferenceValues( refineObj, refVal, refTol )
	type(RemeshSetup), intent(inout) :: refineObj
	real(kreal), intent(in) :: refVal, refTol
	refineObj%refVal = refVal
	refineObj%refTol = refTol
	refineObj%useReferenceVal = .TRUE.
end subroutine


subroutine InitialRefinementPrivate(aMesh, remesh, updateTracerOnMesh, tracerDef, updateVorticityOnMesh, vorticityDef, t)
	type(SphereMesh), intent(inout) :: aMesh
	type(RemeshSetup), intent(in) :: remesh
	procedure(SetTracerOnMesh) :: updateTracerOnMesh
	type(TracerSetup), intent(in) :: tracerDef
	procedure(SetVorticityOnMesh) :: updateVorticityOnMesh
	type(BVESetup), intent(in) :: vorticityDef
	real(kreal), intent(in), optional :: t
	!
	integer(kint) :: refinecount, spaceLeft, counters(6), j
	integer(kint) :: startIndex, nOldPanels, amrLoopCounter
	logical(klog) :: keepGoing
	logical(klog), allocatable :: refineFlag(:)
	type(Panels), pointer :: aPanels

	if ( .NOT. remesh%useAMR ) then
		call LogMessage(log, WARNING_LOGGING_LEVEL, logKey, 'InitialRefinement WARNING : uniform mesh, no AMR.')
		return
	endif

	call StartSection(log, 'InitialRefinement : ')

	aPanels => aMesh%panels
	allocate(refineFlag(aPanels%N_Max))
	refineFlag = .FALSE.
	keepGoing = .FALSE.
	counters = 0
	startIndex = 1
	!
	! apply refinement criteria
	!
	if ( remesh%vorticityRefine) then
		call FlagPanelsForMaxCirculationRefinement(refineFlag, aMesh, remesh, startIndex, counters(1))
		call FlagPanelsForVorticityVariationRefinement(refineFlag, aMesh, remesh, startIndex, counters(2))
	endif
	if ( remesh%flowmapRefine) then
		call FlagPanelsForFlowMapRefinement(refineFlag, aMesh, remesh, startIndex, counters(3))
	endif
	if ( remesh%tracerRefine ) then
		if ( remesh%useReferenceVal ) then
			call FlagPanelsForTracerInterfaceRefinement(refineFlag, aMesh, remesh, startIndex, counters(4))
		else
			call FlagPanelsForTracerMassRefinement(refineFlag, aMesh, remesh, startIndex, counters(4))
		endif
		call FlagPanelsForTracerVariationRefinement(refineFlag, aMesh, remesh, startIndex, counters(5))
	endif

	call FlagPanelsForTransitionRefinement( refineFlag, aMesh, remesh, startIndex, counters(6) )
	
	refineCount = count(refineFlag)
	spaceLeft = aPanels%N_Max - aPanels%N

	if ( refineCount > 0 ) then
		if ( spaceLeft / 4 > refineCount ) then
			keepGoing = .TRUE.
			amrLoopCounter = 0
			do while (keepGoing)
				amrLoopCounter = amrLoopCounter + 1
				write(logstring,'(A,I2,A,I8,A)') 'AMR loop ', amrLoopCounter, ' : refining ', refineCount, ' panels.'
				call LogMessage(log, TRACE_LOGGING_LEVEL, 'InitRefine : ', trim(logstring))
				!
				! divide flagged panels
				!
				nOldPanels = aPanels%N
				do j = 1, nOldPanels
					if ( refineFlag(j) ) then
						call DividePanel(aMesh, j)
						refineFlag(j) = .FALSE.
					endif
				enddo
				!
				!	TO DO : ensure adjacent panels differ by at most one level of refinement
				!

				!
				! set data on new panels and particles
				!
				if ( aMesh%nTracer > 0 ) call UpdateTracerOnMesh(aMesh, tracerDef)
				if ( present(t) ) then
					call UpdateVorticityOnMesh(aMesh, vorticityDef, t)
				else
					call UpdateVorticityOnMesh(aMesh, vorticityDef)
				endif

				!
				! prevent too much refinement
				!
				if ( amrLoopCounter >= remesh%refinementLimit ) then
					keepGoing = .FALSE.
					call LogMessage(log, WARNING_LOGGING_LEVEL, 'InitRefine WARNING : ',' refinement limit reached.')
					call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : nActive = ', aPanels%N_Active)
					write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
					call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
					write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
					call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
					write(logstring,'(A, I8, A)') 'flowmap variation criterion triggered ', counters(3), ' times.'
					call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
					write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
					call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
					write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
					call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
				else ! proceed with next round of amr
					!
					! apply refinement criteria
					!
					startIndex = nOldPanels + 1
					nOldPanels = aPanels%N
					if ( remesh%vorticityRefine) then
						call FlagPanelsForMaxCirculationRefinement(refineFlag, aMesh, remesh, startIndex, counters(1))
						call FlagPanelsForVorticityVariationRefinement(refineFlag, aMesh, remesh, startIndex, counters(2))
					endif
					if ( remesh%flowmapRefine) then
						call FlagPanelsForFlowMapRefinement(refineFlag, aMesh, remesh, startIndex, counters(3))
					endif
					if ( remesh%tracerRefine ) then
						if ( remesh%useReferenceVal ) then
							call FlagPanelsForTracerInterfaceRefinement(refineFlag, aMesh, remesh, startIndex, counters(4))
						else
							call FlagPanelsForTracerMassRefinement(refineFlag, aMesh, remesh, startIndex, counters(4))
						endif
						call FlagPanelsForTracerVariationRefinement(refineFlag, aMesh, remesh, startIndex, counters(5))
					endif
					call FlagPanelsForTransitionRefinement( refineFlag, aMesh, remesh, startIndex, counters(6) )
					
					refineCount = count(refineFlag)
					spaceLeft = aPanels%N_Max - aPanels%N

					!
					! check stopping criteria
					!
					if ( refineCount == 0 ) then
						!
						! refinement has converged
						!
						keepGoing = .FALSE.
						call LogMessage(log, TRACE_LOGGING_LEVEL, 'InitRefine : ', 'refinement converged.')
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : nActive = ', aPanels%N_Active)
						write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') 'flowmap variation criterion triggered ', counters(3), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
					elseif ( spaceLeft / 4 < refineCount ) then
						!
						! not enough memory to proceed
						!
						keepGoing = .FALSE.
						call LogMessage(log, WARNING_LOGGING_LEVEL, 'InitRefine WARNING : ', 'not enough memory to continue AMR.')
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : nActive = ', aPanels%N_Active)
						write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') 'flowmap variation criterion triggered ', counters(3), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
					endif! stopping criteria triggered
				endif! below refinement limit
			enddo ! while keepGoing
		else ! not enough memory
			call LogMessage(log, WARNING_LOGGING_LEVEL, 'InitRefine : ',' not enough memory for AMR.')
		endif
	else ! refine count == 0
		 call LogMessage(log, TRACE_LOGGING_LEVEL, 'InitRefine : ','no refinement necessary.')
	endif
	deallocate(refineFlag)
	call EndSection(log)
end subroutine

subroutine LagrangianRemeshToInitialTimePrivate(aMesh, remesh, setVorticity, vorticityDef, setTracer, tracerDef, t)
	type(SphereMesh), intent(inout) :: aMesh
	type(RemeshSetup), intent(in) :: remesh
	procedure(SetVorticityOnMesh) :: setVorticity
	type(BVESetup), intent(in) :: vorticityDef
	procedure(SetTracerOnMesh) :: setTracer
	type(TracerSetup) :: tracerDef
	real(kreal), intent(in), optional :: t
	!
	type(STRIPACKData) :: delTri
	type(SSRFPACKData) :: lagSource
	type(SphereMesh) :: newMesh
	type(Particles), pointer :: newParticles
	type(Panels), pointer :: newPanels
	integer(kint) :: j, amrLoopCounter, counters(5)
	logical(klog), allocatable :: refineFlag(:)
	logical(klog) :: keepGoing
	integer(kint) :: startIndex, nOldPanels, nOldParticles, refineCount, spaceLeft

	nullify(newParticles)
	nullify(newPanels)
	keepGoing = .FALSE.
	amrLoopCounter = 0
	startIndex = 1
	counters = 0

	call LogMessage(log,DEBUG_LOGGING_LEVEL,'LagRemeshInit : ',' entering Lagrangian remesh.')

	!
	! build a new uniform mesh
	!
	call New(newMesh, aMesh%panelKind, aMesh%initNest, aMesh%AMR, aMesh%nTracer, aMesh%problemKind)
	newParticles => newMesh%particles
	newPanels => newMesh%panels

	!
	! setup source data for interpolation from old mesh
	!
	call New(delTri, aMesh)
	call New(lagSource, delTri, .TRUE.)
	call SetSourceLagrangianParameter(lagSource, delTri)

	!
	! interpolate Lagrangian parameter from old mesh to new mesh
	!
	do j = 1, newParticles%N
		newParticles%x0(:,j) = InterpolateVector(newParticles%x(:,j), lagSource, delTri)
		! renormalize to spherical surface
		newParticles%x0(:,j) = newParticles%x0(:,j) / &
			sqrt(sum(newParticles%x0(:,j)*newParticles%x0(:,j))) * EARTH_RADIUS
	enddo
	do j = 1, newPanels%N
		newPanels%x0(:,j) = InterpolateVector(newPanels%x0(:,j), lagSource, delTri)
		newPanels%x0(:,j) = newPanels%x0(:,j) / &
			sqrt( sum(newPanels%x0(:,j) * newPanels%x0(:,j))) * EARTH_RADIUS
	enddo

	!
	! set flow data on new mesh
	!
	if ( aMesh%nTracer > 0 ) call SetTracer(newMesh, tracerDef)
	if ( present(t) ) then
		call SetVorticity(newMesh, vorticityDef, t)
	else
		call SetVorticity(newMesh, vorticityDef)
	endif

	!
	! AMR
	!
	if ( aMesh%AMR > 0 .AND. remesh%useAMR ) then
		call StartSection(log, 'LagRemeshInit, AMR : ')
		allocate(refineFlag(newPanels%N_Max))
		refineFlag = .FALSE.
		startIndex = 1
		amrLoopCounter = 0
		if ( remesh%vorticityRefine ) then
			call FlagPanelsForMaxCirculationRefinement(refineFlag, newMesh, remesh, startIndex, counters(1))
			call FlagPanelsForVorticityVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(2))
		endif
		if ( remesh%flowMapRefine) then
			call FlagPanelsForFlowMapRefinement(refineFlag, newMesh, remesh, startIndex, counters(3))
		endif
		if ( remesh%tracerRefine ) then
			if ( remesh%useReferenceVal ) then
				call FlagPanelsForTracerInterfaceRefinement(refineFlag, newMesh, remesh, startIndex, counters(4))
			else
				call FlagPanelsForTracerMassRefinement(refineFlag, newMesh, remesh, startIndex, counters(4))
			endif
			call FlagPanelsForTracerVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(5))
		endif
		call FlagPanelsForTransitionRefinement( refineFlag, newMesh, remesh, startIndex, counters(6) )
		
		refineCount = count(refineFlag)
		spaceLeft = newPanels%N_Max - newPanels%N

		if ( refineCount > 0 ) then
			if ( spaceleft / 4 > refineCount ) then
				keepGoing = .TRUE.
				do while (keepGoing)
					amrLoopCounter = amrLoopCounter + 1
					write(logString,*) 'AMR Loop ', amrLoopCounter, ' : refining ', refineCount, ' of ', newMesh%panels%N_Active, ' panels.'
					call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshInitTime : ',trim(logstring))
					!
					! refine flagged panels
					!
					nOldPanels = newPanels%N
					nOldParticles = newParticles%N
					do j = startIndex, newPanels%N
						if ( refineFlag(j) ) then
							call DividePanel(newMesh, j)
							refineFlag(j) = .FALSE.
						endif
					enddo
					!
					! 	TO DO : ensure adjacent panels differ by at most one level of refinement
					!

					!
					! interpolate Lagrangian parameter to new particles and panels
					!
					do j = nOldParticles + 1, newParticles%N
						newParticles%x0(:,j) = InterpolateVector(newParticles%x(:,j), lagSource, delTri)
						newParticles%x0(:,j) = newParticles%x0(:,j) / &
							sqrt( sum( newParticles%x0(:,j)*newParticles%x0(:,j))) * EARTH_RADIUS
					enddo
					do j = nOldPanels+1, newPanels%N
						newPanels%x0(:,j) = InterpolateVector(newPanels%x(:,j), lagSource, delTri)
						newPanels%x0(:,j) = newPanels%x0(:,j) / &
							sqrt( sum( newPanels%x0(:,j) * newPanels%x0(:,j)) ) * EARTH_RADIUS
					enddo

					!
					! set flow data on new particles, panels
					!
					if ( aMesh%nTracer > 0 ) call SetTracer(newMesh, tracerDef)
					if ( present(t) ) then
						call SetVorticity(newMesh, vorticityDef, t)
					else
						call SetVorticity(newMesh, vorticityDef)
					endif

					if ( amrLoopCounter >= remesh%refinementLimit ) then
						!
						! prevent too much refinement
						!
						keepGoing = .FALSE.
						call LogMessage(log, WARNING_LOGGING_LEVEL, 'LagRemeshInitTime WARNING : ', 'refinement limit reached.')
						call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshInitTime : nActive = ', newPanels%N_Active)
						write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitial : ',trim(logstring))
						write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitial : ',trim(logstring))
						write(logstring,'(A, I8, A)') 'flowmap variation criterion triggered ', counters(3), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitial : ',trim(logstring))
						write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitial : ',trim(logstring))
						write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitial : ',trim(logstring))
					else
						!
						! apply refinement critera
						!
						startIndex = nOldPanels + 1
						nOldPanels = newPanels%N
						if ( remesh%vorticityRefine ) then
							call FlagPanelsForMaxCirculationRefinement(refineFlag, newMesh, remesh, startIndex, counters(1) )
							call FlagPanelsForVorticityVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(2) )
						endif
						if ( remesh%flowMapRefine ) then
							call FlagPanelsForFlowMapRefinement(refineFlag, newMesh, remesh, startIndex, counters(3))
						endif
						if ( remesh%tracerRefine ) then
							if ( remesh%useReferenceVal ) then
								call FlagPanelsForTracerInterfaceRefinement(refineFlag, newMesh, remesh, startIndex, counters(4))
							else
								call FlagPanelsForTracerMassRefinement(refineFlag, newMesh, remesh, startIndex, counters(4))
							endif
							call FlagPanelsForTracerVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(5))
						endif
						call FlagPanelsForTransitionRefinement( refineFlag, newMesh, remesh, startIndex, counters(6) )

						refineCount = count(refineFlag)
						spaceLeft = newPanels%N_Max - newPanels%N

						if ( refineCount == 0 ) then
							keepGoing = .FALSE.
							call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshInitTime : ', ' refinement converged.')
							call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshInitTime : nActive = ', newPanels%N_Active)
							write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitial : ',trim(logstring))
							write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitial : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'flowmap variation criterion triggered ', counters(3), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitial : ',trim(logstring))
							write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitial : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitial : ',trim(logstring))
						elseif ( spaceLeft / 4 < refineCount ) then
							keepGoing = .FALSE.
							call LogMessage(log, WARNING_LOGGING_LEVEL, 'LagRemeshInitTime WARNING : ', 'not enough memory to continue AMR.')
							call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshInitTime : nActive = ', newPanels%N_Active)
							write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitial : ',trim(logstring))
							write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitial : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'flowmap variation criterion triggered ', counters(3), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitial : ',trim(logstring))
							write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitial : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitial : ',trim(logstring))
						endif! stopping criteria triggered
					endif ! limit reached
				enddo! while keepgoing
			else ! not enough memory for amr
				call LogMessage(log, WARNING_LOGGING_LEVEL, 'LagRemeshInitTime WARNING : ', 'not enough memory for AMR.')
			endif
		else ! no refinement necessary
			call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshInitTime : ', ' no refinement necessary.')
		endif ! refineCount > 0

		deallocate(refineFlag)
		call EndSection(log)
	endif ! AMR

	!
	! replace old mesh with new mesh
	!
	call Copy(aMesh, newMesh)

	!
	! clean up
	!
	call Delete(newMesh)
	call Delete(delTri)
	call Delete(lagSource)
end subroutine

subroutine LagrangianRemeshToReference(aMesh, reference, remesh, setVorticity, vortDef, t)
	type(SphereMesh), intent(inout) :: aMesh
	type(ReferenceSphere), intent(in) :: reference
	type(RemeshSetup), intent(in) :: remesh
	procedure(SetVorticityOnMesh), optional :: setVorticity
	type(BVESetup), intent(in), optional :: vortDef
	real(kreal), intent(in), optional :: t
	!
	type(SphereMesh) :: newMesh
	type(Particles), pointer :: newParticles
	type(Panels), pointer :: newPanels
	type(STRIPACKData) :: delTri
	type(SSRFPACKData) :: lagSource
	integer(kint) :: j, k, amrLoopCounter, counters(5)
	logical(klog), allocatable :: refineFlag(:)
	logical(klog) :: keepGoing
	integer(kint) :: startIndex, nOldPanels, nOldParticles, refineCount, spaceLeft

	nullify(newParticles)
	nullify(newPanels)
	amrLoopCounter = 0
	counters = 0
	startIndex = 1
	keepGoing = .FALSE.

	call LogMessage(log, DEBUG_LOGGING_LEVEL, logkey, ' entering LagrangianRemeshToReferenceTime')

	!
	! build a new mesh
	!
	call New(newMesh, aMesh%panelKind, aMesh%initNest, aMesh%AMR, aMesh%nTracer, aMesh%problemKind)
	newParticles => newMesh%particles
	newPanels => newMesh%panels
 
	!
	! setup old mesh as source for interpolation of Lagrangian parameters
	!
	call New(deltri, aMesh)
	call New(lagSource, delTri, .TRUE.)
	call SetSourceLagrangianParameter(lagSource, delTri)

	!
	! interpolate Lagrangian parameter from old mesh to new mesh
	!
	do j = 1, newParticles%N
		newParticles%x0(:,j) = InterpolateVector(newParticles%x(:,j), lagSource, delTri)
		newParticles%x0(:,j) = newParticles%x0(:,j) / &
			sqrt( sum( newParticles%x0(:,j) * newParticles%x0(:,j) ) ) * EARTH_RADIUS
	enddo
	do j = 1, newPanels%N
		newPanels%x0(:,j) = InterpolateVector(newPanels%x0(:,j), lagSource, delTri)
		newPanels%x0(:,j) = newPanels%x0(:,j) / &
			sqrt( sum( newPanels%x0(:,j) * newPanels%x0(:,j))) * EARTH_RADIUS
	enddo

	!
	! set flow data on new mesh
	!
	do k = 1, aMesh%nTracer
		do j = 1, newParticles%N
			newParticles%tracer(j,k) = InterpolateTracer(newParticles%x0(:,j),&
					reference%tracerSource, reference%delTri, k)
		enddo

		do j = 1, newPanels%N
			newPanels%tracer(j,k) = InterpolateTracer(newPanels%x0(:,j), &
					reference%tracerSource, reference%delTri, k)
		enddo
	enddo

	if ( aMesh%problemKind == BVE_SOLVER) then
		if ( present(t) ) then
			call SetVorticity(newMesh, vortDef, t)	
		else
			do j = 1, newParticles%N
				newParticles%absVort(j) = InterpolateScalar(newParticles%x0(:,j), &
						reference%absVortSource, reference%delTri)
				newParticles%relVort(j) = newParticles%absVort(j) - &
					2.0_kreal * OMEGA * newParticles%x(3,j) / EARTH_RADIUS
			enddo
			do j = 1, newPanels%N
				newPanels%absVort(j) = InterpolateScalar(newPanels%x0(:,j), &
						reference%absVortSource, reference%delTri)
				newPanels%relVort(j) = newPanels%absVort(j) - &
					2.0_kreal * OMEGA * newPanels%x(3,j) / EARTH_RADIUS
			enddo
		endif
	endif

	call LogMessage(log, DEBUG_LOGGING_LEVEL,'LagRemeshToRef : ',' uniform mesh ready.')

	!
	! AMR
	!
	if ( aMesh%AMR > 0 .AND. remesh%useAMR ) then
		call StartSection(log, 'LagRemeshToRef AMR : ')
		allocate(refineFlag(newPanels%N_Max))
		refineFlag = .FALSE.
		amrLoopCounter = 0

		!
		! apply amr criteria
		!
		if ( remesh%vorticityRefine ) then
			call FlagPanelsForMaxCirculationRefinement(refineFlag, newMesh, remesh, startIndex, counters(1))
			call FlagPanelsForVorticityVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(2))
		endif
		if ( remesh%flowMapRefine ) then
			call FlagPanelsForFlowMapRefinement(refineFlag, newMesh, remesh, startIndex, counters(3))
		endif
		if ( remesh%tracerRefine ) then
			if ( remesh%useReferenceVal ) then
				call FlagPanelsForTracerInterfaceRefinement(refineFlag, newMesh, remesh, startIndex, counters(4))
			else
				call FlagPanelsForTracerMassRefinement(refineFlag, newMesh, remesh, startIndex, counters(4))
			endif
			call FlagPanelsForTracerVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(5))
		endif
		call FlagPanelsForTransitionRefinement( refineFlag, newMesh, remesh, startIndex, counters(6) )

		refineCount = count(refineFlag)
		spaceLeft = newPanels%N_Max - newPanels%N

		if ( refineCount > 0 ) then
			if ( spaceLeft / 4 > refineCount ) then
				keepGoing = .TRUE.
				do while (keepGoing)
					amrLoopCounter = amrLoopCounter + 1
					write(logString,'(A,I2,A,I8,A)') 'AMR loop ', amrLoopCounter, ' : refining ', refineCount,' panels.'
					call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshToRef : ', trim(logstring))

					!
					! divide flagged panels
					!
					nOldPanels = newPanels%N
					nOldParticles = newParticles%N
					do j = startIndex, nOldPanels
						if ( refineFlag(j) ) then
							call DividePanel(newMesh, j)
							refineFlag(j) = .FALSE.
						endif
					enddo

					!
					! TO DO : ensure adjacent panels differ by no more than one level of refinement
					!

					!
					! interpolate Lagrangian parameter to new particles, panels
					!
					do j = nOldParticles + 1, newParticles%N
						newParticles%x0(:,j) = InterpolateVector(newParticles%x(:,j), lagSource, delTri)
						newParticles%x0(:,j) = newParticles%x0(:,j) / &
							sqrt( sum( newParticles%x0(:,j) * newParticles%x0(:,j)) ) * EARTH_RADIUS
					enddo
					do j = nOldPanels + 1, newPanels%N
						newPanels%x0(:,j) = InterpolateVector(newPanels%x(:,j), lagSource, delTri)
						newpanels%x0(:,j) = newPanels%x0(:,j) / &
							sqrt( sum( newPanels%x0(:,j) * newPanels%x0(:,j) ) ) * EARTH_RADIUS
					enddo

					!
					! set flow data on new particles, panels
					!
					do k = 1, aMesh%nTracer
						do j = nOldParticles+1, newParticles%N
							newParticles%tracer(j,k) = INterpolateTracer(newParticles%x0(:,j), &
								reference%tracerSource, reference%delTri, k)
						enddo
						do j = nOldPanels+1, newPanels%N
							newPanels%tracer(j,k) = InterpolateTracer(newPanels%x0(:,j), &
								reference%tracerSource, reference%delTri, k)
						enddo
					enddo

					if ( aMesh%problemKind == BVE_SOLVER) then
						if ( present(t) ) then
							call SetVorticity(newMesh, vortDef, t)	
						else
							do j = nOldParticles+1, newParticles%N
								newParticles%absVort(j) = InterpolateScalar( newParticles%x0(:,j), &
									reference%absVortSource, reference%delTri)
								newParticles%relVort(j) = newParticles%absVort(j) - &
									2.0_kreal * OMEGA * newParticles%x(3,j) / EARTH_RADIUS
							enddo
							do j = nOldPanels + 1, newPanels%N
								newPanels%absVort(j) = InterpolateScalar(newPanels%x0(:,j), &
									reference%absVortSource, reference%delTri)
								newPanels%relVort(j) = newPanels%absVort(j) - &
									2.0_kreal * OMEGA * newPanels%x(3,j) / EARTH_RADIUS
							enddo
						endif
					endif

					!
					! prevent too much refinement
					!
					if ( amrLoopCounter >= remesh%refinementLimit ) then
						keepGoing = .FALSE.
						call LogMessage(log, WARNING_LOGGING_LEVEL, 'LagRemeshToRef WARNING : ', 'refinement limit reached.')
						call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshToRef N_Active = ', newPanels%N_Active)
						write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshToRef : ',trim(logstring))
						write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshToRef : ',trim(logstring))
						write(logstring,'(A, I8, A)') 'flwmap variation criterion triggered ', counters(3), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshToRef : ',trim(logstring))
						write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshToRef : ',trim(logstring))
						write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshToRef : ',trim(logstring))
					else
						!
						! apply refinement criteria to new panels
						!
						startIndex = nOldPanels + 1
						nOldPanels = newPanels%N
						if ( remesh%vorticityRefine ) then
							call FlagPanelsForMaxCirculationRefinement(refineFlag, newMesh, remesh, startIndex, counters(1) )
							call FlagPanelsForVorticityVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(2) )
						endif
						if ( remesh%flowMapRefine ) then
							call FlagPanelsForFlowMapRefinement(refineFlag, newMesh, remesh, startIndex, counters(3))
						endif
						if ( remesh%tracerRefine ) then
							if ( remesh%useReferenceVal ) then
								call FlagPanelsForTracerInterfaceRefinement(refineFlag, newMesh, remesh, startIndex, counters(4))
							else
								call FlagPanelsForTracerMassRefinement(refineFlag, newMesh, remesh, startIndex, counters(4))
							endif
							call FlagPanelsForTracerVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(5))
						endif
						call FlagPanelsForTransitionRefinement( refineFlag, newMesh, remesh, startIndex, counters(6) )

						refineCount = count(refineFlag)
						spaceLeft = newPanels%N_Max - newPanels%N

						!
						! check stopping criteria
						!
						if ( refineCount == 0 ) then
							keepGoing = .FALSE.
							call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshToRef : ', ' refinement converged.')
							call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshToRef N_Active = ', newPanels%N_Active)
							write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshToRef : ',trim(logstring))
							write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshToRef : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'flwmap variation criterion triggered ', counters(3), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshToRef : ',trim(logstring))
							write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshToRef : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshToRef : ',trim(logstring))
						elseif ( spaceLeft / 4 < refineCount ) then
							keepGoing = .FALSE.
							call LogMessage(log, WARNING_LOGGING_LEVEL, 'LagRemeshToRef WARNING: ', ' not enough memory to continue AMR.')
							call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshToRef N_Active = ', newPanels%N_Active)
							write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshToRef : ',trim(logstring))
							write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshToRef : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'flwmap variation criterion triggered ', counters(3), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshToRef : ',trim(logstring))
							write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshToRef : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshToRef : ',trim(logstring))
						endif ! stopping criteria triggered
					endif ! below refinement limit
				enddo ! while keepgoing
			else ! not enough memory for AMR
				call LogMessage(log, WARNING_LOGGING_LEVEL,'LagRemeshToRef WARNING : ', ' not enough memory for AMR.')
			endif
		else ! refince count == 0 , no refinement necessary
			call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshToRef : ', 'no refinement necessary.')
		endif ! refinecoutn > 0
		deallocate(refineFlag)
		call EndSection(log)
	endif! AMR

	!
	! replace old mesh with new mesh
	!
	call Copy(aMesh, newMesh)

	!
	! clean up
	!
	call Delete(delTri)
	call Delete(lagSource)
	call Delete(newMesh)
end subroutine

subroutine DirectRemesh( aMesh, remesh )
	type(SphereMesh), intent(inout) :: aMesh
	type(RemeshSetup), intent(in) :: remesh
	!
	type(STRIPACKData) :: delTri
	type(SSRFPACKData) :: vectorSource, absVortSource, tracerSource
	type(SphereMesh) :: newMesh
	type(Particles), pointer :: newParticles
	type(Panels), pointer :: newPanels
	integer(kint) :: j, k,  amrLoopCounter, counters(4)
	logical(klog), allocatable :: refineFlag(:)
	logical(klog) :: keepGoing
	integer(kint) :: startIndex, nOldPanels, nOldParticles, refineCount, spaceLeft

	nullify(newParticles)
	nullify(newPanels)
	keepGoing = .FALSE.
	startIndex = 1
	amrLoopCounter = 0
	counters = 0

	call LogMessage(log, DEBUG_LOGGING_LEVEL, 'DirectRemesh : ', 'entering.')

	!
	! build a new uniform mesh
	!
	call New(newMesh, aMesh%panelKind, aMesh%initNest, aMesh%AMR, aMesh%nTracer, aMesh%problemKind)
	newParticles => newMesh%particles
	newPanels => newMesh%panels

	!
	! setup source data for interpolation from old mesh
	!
	call New(delTri, aMesh)
!	call New(vectorSource, delTri, .TRUE.)
!	call New(scalarSource, delTri, .FALSE.)

	!
	! set flow data on new mesh
	!
	if ( aMesh%nTracer > 0 ) then
		call New(tracerSource, delTri, aMesh%nTracer)
		call SetSourceTracer(tracerSource, delTri)

		do k = 1, aMesh%nTracer
			do j = 1, newParticles%N
				newParticles%tracer(j, k) = InterpolateTracer(newParticles%x(:,j), tracerSource, delTri, k)
			enddo
			do j = 1, newPanels%N
				if ( newPanels%hasChildren(j) ) then
					newPanels%tracer(j, k) = 0.0_kreal
				else
					newPanels%tracer(j, k) = InterpolateTracer(newPanels%x(:,j), tracerSource, delTri, k)
				endif
			enddo
		enddo
	endif

	if ( aMesh%problemKind == BVE_SOLVER) then
		call New(absVortSource, delTri, .FALSE.)
		call SetSourceAbsVort(absVortSource, delTri)

		do j = 1, newParticles%N
			newParticles%absVort(j) = InterpolateScalar(newParticles%x(:,j), absVortSource, delTri)
			newParticles%relVort(j) = newParticles%absVort(j) - 2.0_kreal * OMEGA * newParticles%x(3,j)/EARTH_RADIUS
		enddo
		do j = 1, newPanels%N
			if ( newPanels%hasChildren(j) ) then
				newPanels%absVort(j) = 0.0_kreal
			else
				newPanels%absVort(j) = InterpolateScalar(newPanels%x(:,j), absVortSource, deltri)
				newPanels%relVort(j) = newPanels%absVort(j) - 2.0_kreal * OMEGA * newPanels%x(3,j)/EARTH_RADIUS
			endif
		enddo

		!
		! TO DO : Stream function !
		!

	endif

	!
	! TO DO : SWE Variables
	!


	!
	! AMR
	!
	if ( aMesh%AMR > 0 .AND. remesh%useAMR ) then
		call StartSection(log, 'DirectRemesh AMR : ')
		allocate(refineFlag(newPanels%N_Max))
		refineFlag = .FALSE.
		startIndex = 1
		amrLoopCounter = 0
		if ( remesh%vorticityRefine ) then
			call FlagPanelsForMaxCirculationRefinement(refineFlag, newMesh, remesh, startIndex, counters(1))
			call FlagPanelsForVorticityVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(2))
		endif
		if ( remesh%tracerRefine ) then
			call FlagPanelsForTracerMassRefinement(refineFlag, newMesh, remesh, startIndex, counters(3))
			call FlagPanelsForTracerVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(4))
		endif
		call FlagPanelsForTransitionRefinement( refineFlag, newMesh, remesh, startIndex, counters(6) )

		refineCount = count(refineFlag)
		spaceLeft = newPanels%N_Max - newPanels%N

		if ( refineCount > 0 ) then
			if ( spaceleft / 4 > refineCount ) then
				keepGoing = .TRUE.
				do while (keepGoing)
					amrLoopCounter = amrLoopCounter + 1
					write(logString,'(A,I2,A,I8,A)') 'AMR Loop ', amrLoopCounter, ' : refining ', refineCount, ' panels.'
					call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshInitTime : ',trim(logstring))
					!
					! refine flagged panels
					!
					nOldPanels = newPanels%N
					nOldParticles = newParticles%N
					do j = startIndex, newPanels%N
						if ( refineFlag(j) ) then
							call DividePanel(newMesh, j)
							refineFlag(j) = .FALSE.
						endif
					enddo

					do k = 1, aMesh%nTracer
						do j = nOldParticles + 1, newParticles%N
							newParticles%tracer(j, k) = InterpolateTracer(newParticles%x(:,j), tracerSource, delTri, k)
						enddo
						do j = nOldPanels + 1, newPanels%N
							newPanels%tracer(j, k) = InterpolateTracer(newPanels%x(:,j), tracerSource, delTri, k)
						enddo
					enddo

					if ( aMesh%problemKind == BVE_SOLVER) then
						do j = nOldParticles + 1, newParticles%N
							newParticles%absVort(j) = InterpolateScalar(newParticles%x(:,j), absVortSource, delTri)
							newParticles%relVort(j) = newParticles%absVort(j) - 2.0_kreal * OMEGA * newParticles%x(3,j)/EARTH_RADIUS
						enddo
						do j = nOldPanels + 1, newPanels%N
							newPanels%absVort(j) = InterpolateScalar(newPanels%x(:,j), absVortSource, delTri)
							newPanels%relVort(j) = newPanels%absVort(j) - 2.0_kreal * OMEGA * newPanels%x(3,j) / EARTH_RADIUS
						enddo
						!
						! TO DO : stream function
						!
					endif

					!
					! TO DO : SWE Variables
					!
					if ( amrLoopCounter >= remesh%refinementLimit ) then
						!
						! prevent too much refinement
						!
						keepGoing = .FALSE.
						call LogMessage(log, WARNING_LOGGING_LEVEL, 'DirectRemesh WARNING : ', 'refinement limit reached.')
						call LogMessage(log, TRACE_LOGGING_LEVEL, 'DirectRemesh : nActive = ', newPanels%N_Active)
						write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'DirectRemesh : ',trim(logstring))
						write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'DirectRemesh : ',trim(logstring))
						write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(3), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'DirectRemesh : ',trim(logstring))
						write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(4), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'DirectRemesh : ',trim(logstring))
					else
						!
						! apply refinement critera
						!
						startIndex = nOldPanels + 1
						nOldPanels = newPanels%N
						if ( remesh%vorticityRefine ) then
							call FlagPanelsForMaxCirculationRefinement(refineFlag, newMesh, remesh, startIndex, counters(1))
							call FlagPanelsForVorticityVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(2))
						endif
						if ( remesh%tracerRefine ) then
							call FlagPanelsForTracerMassRefinement(refineFlag, newMesh, remesh, startIndex, counters(3))
							call FlagPanelsForTracerVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(4))
						endif
						call FlagPanelsForTransitionRefinement( refineFlag, newMesh, remesh, startIndex, counters(6) )

						refineCount = count(refineFlag)
						spaceLeft = newPanels%N_Max - newPanels%N

						if ( refineCount == 0 ) then
							keepGoing = .FALSE.
							call LogMessage(log, TRACE_LOGGING_LEVEL, 'DirectRemesh : ', 'refinement converged.')
							call LogMessage(log, TRACE_LOGGING_LEVEL, 'DirectRemesh : nActive = ', newPanels%N_Active)
							write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'DirectRemesh : ',trim(logstring))
							write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'DirectRemesh : ',trim(logstring))
							write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(3), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'DirectRemesh : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(4), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'DirectRemesh : ',trim(logstring))
						elseif ( spaceLeft / 4 < refineCount ) then
							keepGoing = .FALSE.
							call LogMessage(log, WARNING_LOGGING_LEVEL, 'DirectRemesh WARNING : ', 'not enough memory to continue AMR.')
							call LogMessage(log, TRACE_LOGGING_LEVEL, 'DirectRemesh : nActive = ', newPanels%N_Active)
							write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'DirectRemesh : ',trim(logstring))
							write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'DirectRemesh : ',trim(logstring))
							write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(3), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'DirectRemesh : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(4), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'DirectRemesh : ',trim(logstring))
						endif! stopping criteria triggered
					endif ! limit reached
				enddo! while keepgoing
			else ! not enough memory for amr
				call LogMessage(log, WARNING_LOGGING_LEVEL, 'DirectRemesh WARNING : ', 'not enough memory for AMR.')
			endif
		else ! no refinement necessary
			call LogMessage(log, TRACE_LOGGING_LEVEL, 'DirectRemesh : ', ' no refinement necessary.')
		endif ! refineCount > 0
		deallocate(refineFlag)
		call EndSection(log)
	endif ! amr

	!
	! replace old mesh with new mesh
	!

	call Copy(aMesh, newMesh)

	!
	! clean up
	!
	call Delete(newMesh)
	call Delete(delTri)
	if ( amesh%nTracer > 0 ) call Delete(tracerSource)
	if ( aMesh%problemKind == BVE_SOLVER ) call Delete(absVortSource)

end subroutine

!subroutine ResetLagrangianParameter(ref)
!	type(ReferenceSphere), intent(inout) :: ref
!end subroutine

!
!----------------
! Module methods : module- or type-specific private functions
!----------------
!

subroutine FlagPanelsForMaxCirculationRefinement(refineFlag, aMesh, remesh, startIndex, counter)
	logical(klog), intent(inout) :: refineFlag(:)
	type(SphereMesh), intent(in) :: aMesh
	type(RemeshSetup), intent(in) :: remesh
	integer(kint), intent(in) :: startIndex
	integer(kint), intent(inout) :: counter
	!
	type(Panels), pointer :: aPanels
	integer(kint) :: j

	aPanels => aMesh%panels

	do j = startIndex, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			if ( abs(aPanels%relVort(j)) * aPanels%area(j) > remesh%maxCircTol ) then
				refineFlag(j) = .TRUE.
				counter = counter + 1
			endif
		endif
	enddo
end subroutine

subroutine FlagPanelsForTracerInterfaceRefinement( refineFlag, aMEsh, remesh, startIndex, counter )
	logical(klog), intent(inout) :: refineFlag(:)
	type(SphereMesh), intent(in) :: aMesh
	type(RemeshSetup), intent(in) :: remesh
	integer(kint), intent(in) :: startIndex
	integer(kint), intent(inout) :: counter
	!
	type(Panels), pointer :: aPanels
	integer(kint) :: j
	
	if ( .NOT. remesh%useReferenceVal ) then
		call LogMessage(log, WARNING_LOGGING_LEVEL, 'FlagPanelsTracerInterface WARNING : ', 'reference values not set.')
		return
	endif
	
	aPanels=>aMesh%panels
	do j = 1, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			if ( abs(aPanels%tracer(j, remesh%tracerID) - remesh%refVal)*aPanels%area(j)/EARTH_SURFACE_AREA < remesh%refTol ) then
				refineFlag(j) = .TRUE.
				counter = counter + 1
			endif		
		endif
	enddo
end subroutine

subroutine FlagPanelsForVorticityVariationRefinement(refineFlag, amesh, remesh, startINdex, counter)
	logical(klog), intent(inout) :: refineFlag(:)
	type(SphereMesh), intent(in) :: aMesh
	type(RemeshSetup), intent(in) :: remesh
	integer(kint), intent(in) :: startIndex
	integer(kint), intent(inout) :: counter
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j, k, edgeList(8), vertList(8), nVerts
	real(kreal) :: minVort, maxVort

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	do j = startIndex, aPanels%N
		if (.NOT. aPanels%hasChildren(j) ) then
			maxVort = aPanels%relVort(j)
			minVort = aPanels%relVort(j)
			call CCWEdgesAndParticlesAroundPanel(edgeList, vertList, nVerts, aMesh, j)
			do k = 1, nVerts
				if ( aParticles%relVort(vertList(k)) > maxVort) maxVort = aParticles%relVort(vertList(k))
				if ( aParticles%relVort(vertList(k)) < minVort) minVort = aParticles%relVort(vertList(k))
			enddo

			if ( maxVort - minVort > remesh%vortVarTol ) then
				refineFlag(j) = .TRUE.
				counter = counter + 1
			endif
		endif
	enddo
end subroutine

subroutine FlagPanelsForTransitionRefinement( refineFlag, aMesh, remesh, startIndex, counter )
	logical(klog), intent(inout) :: refineFlag(:)
	type(SphereMesh), intent(inout) :: aMesh
	type(RemeshSetup), intent(in) :: remesh
	integer(kint), intent(in) :: startIndex
	integer(kint), intent(inout) :: counter
	!
	integer(kint) :: adjPanels(8), adjNestLevels(8), nAdj, j, k
	type(Panels), pointer :: aPanels
	
	aPanels=>aMesh%panels
	
	do j = 1, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			call FindAdjacentPanels(adjPanels, nAdj, aMesh, j)
			do k = 1, nAdj
				adjNestLevels(k) = aPanels%nest( adjPanels(k) )
			enddo	
			if ( maxval(adjNestLevels) - aPanels%nest(j) >= 2 ) then
				refineFlag(j) = .TRUE.
				counter = counter + 1
			endif
		endif
	enddo
end subroutine

subroutine FlagPanelsForFLowMapRefinement(refineFlag, amesh, remesh, startINdex, counter)
	logical(klog), intent(inout) :: refineFlag(:)
	type(SphereMesh), intent(in) :: aMesh
	type(RemeshSetup), intent(in) :: remesh
	integer(kint), intent(in) :: startIndex
	integer(kint), intent(inout) :: counter
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j, k, edgeList(8), vertList(8), nVerts
	real(kreal) :: maxx0(3), minx0(3)

		aParticles => aMesh%particles
	aPanels => aMesh%panels

	do j = startIndex, aPanels%N
		if (.NOT. aPanels%hasChildren(j) ) then
			maxx0 = aPanels%x0(:,j)
			minx0 = aPanels%x0(:,j)
			call CCWEdgesAndParticlesAroundPanel(edgeList, vertList, nVerts, aMesh, j)
			do k = 1, nVerts
				if ( aParticles%x0(1,vertLIst(k)) > maxx0(1) ) maxx0(1) = aParticles%x0(1,vertlist(k))
				if ( aParticles%x0(1,vertList(k)) < minx0(1) ) minx0(1) = aParticles%x0(1,vertList(k))
				if ( aParticles%x0(2,vertLIst(k)) > maxx0(2) ) maxx0(2) = aParticles%x0(2,vertlist(k))
				if ( aParticles%x0(2,vertList(k)) < minx0(2) ) minx0(2) = aParticles%x0(2,vertList(k))
				if ( aParticles%x0(3,vertLIst(k)) > maxx0(3) ) maxx0(3) = aParticles%x0(3,vertlist(k))
				if ( aParticles%x0(3,vertList(k)) < minx0(3) ) minx0(3) = aParticles%x0(3,vertList(k))
			enddo

			if ( sum(maxx0 - minx0) > remesh%lagVarTol ) then
				refineFlag(j) = .TRUE.
				counter = counter + 1
			endif
		endif
	enddo
end subroutine

subroutine FlagPanelsForTracerMassRefinement(refineFlag, aMesh, remesh, startIndex, counter)
	logical(klog), intent(inout) :: refineFlag(:)
	type(SphereMesh), intent(in) :: aMesh
	type(RemeshSetup), intent(in) :: remesh
	integer(kint), intent(in) :: startIndex
	integer(kint), intent(inout) :: counter
	!
	type(Panels), pointer :: aPanels
	integer(kint) :: j

	aPanels => aMesh%panels

	do j = startIndex, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			if ( aPanels%tracer(j,remesh%tracerID) * aPanels%area(j) > remesh%tracerMassTol ) then
				refineFlag(j) = .TRUE.
				counter = counter + 1
			endif
		endif
	enddo
end subroutine

subroutine FlagPanelsForTracerVariationRefinement(refineFlag, amesh, remesh, startINdex, counter)
	logical(klog), intent(inout) :: refineFlag(:)
	type(SphereMesh), intent(in) :: aMesh
	type(RemeshSetup), intent(in) :: remesh
	integer(kint), intent(in) :: startIndex
	integer(kint), intent(inout) :: counter
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j, k, edgeList(8), vertList(8), nVerts
	real(kreal) :: minTracer, maxTracer

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	do j = startIndex, aPanels%N
		if (.NOT. aPanels%hasChildren(j) ) then
			maxTracer = aPanels%tracer(j, remesh%tracerID)
			minTracer = aPanels%tracer(j, remesh%tracerID)
			call CCWEdgesAndParticlesAroundPanel(edgeList, vertList, nVerts, aMesh, j)
			do k = 1, nVerts
				if ( aParticles%tracer(vertList(k),remesh%tracerID) > maxTracer) &
						maxTracer = aParticles%tracer(vertList(k), remesh%tracerID)
				if ( aParticles%tracer(vertList(k),remesh%tracerID) < minTracer) &
						minTracer = aParticles%tracer(vertList(k),remesh%tracerID)
			enddo

			if ( maxTracer - minTracer > remesh%tracerVarTol ) then
				refineFlag(j) = .TRUE.
				counter = counter + 1
			endif
		endif
	enddo
end subroutine

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
