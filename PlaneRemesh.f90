module PlaneRemeshModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method
!------------------------------------------------------------------------------
!
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
!> @defgroup PlaneRemesh Plane Remesh
!> Data types, subroutines, and functions for adaptive remeshing planar particle/panel meshes.
!
!
! DESCRIPTION:
!> @file
!> Data types, subroutines, and functions for adaptive remeshing planar particle/panel meshes.
!
!------------------------------------------------------------------------------
use NumberKindsModule
use LoggerModule
use PlaneGeomModule
use ParticlesModule
use EdgesModule
use PanelsModule
use PlaneMeshModule
use PlaneVorticityModule
use PlaneTracerModule
use BIVARInterfaceModule
use BIVARModule

implicit none

private
public RemeshSetup
public New, Delete
public InitialRefinement
public LagrangianRemeshToInitialTime, LagrangianRemeshToReferenceTime


!
!----------------
! Types and module constants
!----------------
!
type RemeshSetup
	logical(klog) :: uniformMesh			! if TRUE, no AMR
	logical(klog) :: vorticityRefine		! if TRUE, remeshing will apply vorticity AMR criteria
	real(kreal) :: maxCircTol
	real(kreal) :: vortVarTol
	logical(klog) :: flowMapRefine			! if TRUE, remeshing will apply flow map AMR criteria
	real(kreal) :: lagVarTol
	logical(klog) :: tracerRefine			! if TRUE, remeshing will apply tracer AMR criteria
	integer(kint) :: tracerID
	real(kreal) :: maxMassTol
	real(kreal) :: massVarTol
	integer(kint) :: refinementLimit
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
	module procedure NewPrivateNoAMR
end interface

interface Delete
	module procedure DeletePrivate
end interface

interface InitialRefinement
	module procedure InitialRefinementPrivate
end interface

interface LagrangianRemeshToInitialTime
	module procedure LagrangianRemeshToInitialTimePrivate
end interface

interface LagrangianRemeshToReferenceTime
	module procedure LagrangianRemeshToReferenceTimePrivate
end interface

interface
	subroutine SetVorticityOnMesh( genMesh, genVort)
		use PlaneMeshModule
		use PlaneVorticityModule
		implicit none
		type(PlaneMesh), intent(inout) :: genMesh
		type(VorticitySetup), intent(in) :: genVort
	end subroutine
end interface

interface
	subroutine SetTracerOnMesh( genMesh, genTracer )
		use PlaneMeshModule
		use PlaneTracerModule
		implicit none
		type(PlaneMesh), intent(inout) :: genMesh
		type(TracerSetup), intent(in) :: genTracer
	end subroutine
end interface

!
!----------------
! Logging
!----------------
!

logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'PlaneRemesh'
integer(kint), parameter :: logLevel = TRACE_LOGGING_LEVEL
character(len=128) :: logString
character(len=24) :: formatString

contains
!
!----------------
! Standard methods : Constructor / Destructor
!----------------
!
subroutine NewPrivateAll(self, maxCircTol, vortVarTol, lagVarTol, tracerID, maxMassTol, massVarTol, limit)
	type(RemeshSetup), intent(out) :: self
	real(kreal), intent(in) :: maxCircTol, vortVarTol, lagVarTol
	integer(kint), intent(in) :: tracerID
	real(kreal), intent(in) :: maxMassTol, massVarTol
	integer(kint), intent(in) :: limit

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	self%uniformMesh = .FALSE.
	self%vorticityRefine = .TRUE.
	self%maxCircTol = maxCircTol
	self%vortVarTol = vortVarTol
	self%flowMapRefine = .TRUE.
	self%lagVarTol = lagVarTol
	self%tracerRefine = .TRUE.
	self%tracerID = tracerID
	self%maxMassTol = maxMassTol
	self%massVarTol = massVarTol
	self%refinementLimit = limit
end subroutine

subroutine NewPrivateVorticity(self, maxCircTol, vortVarTol, limit)
	type(RemeshSetup), intent(out) :: self
	real(kreal), intent(in) :: maxCircTol, vortVarTol
	integer(kint), intent(in) :: limit

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	self%uniformMesh = .FALSE.
	self%vorticityRefine = .TRUE.
	self%maxCircTol = maxCircTol
	self%vortVarTol = vortVarTol
	self%flowMapRefine = .FALSE.
	self%lagVarTol = 0.0_kreal
	self%tracerRefine = .FALSE.
	self%tracerID = 0
	self%maxMassTol = 0.0_kreal
	self%massVarTol = 0.0_kreal
	self%refinementLimit = limit
end subroutine

subroutine NewPrivateVorticityAndFlowMap(self, maxCircTol, vortVarTol, lagVarTol, limit)
	type(RemeshSetup), intent(out) :: self
	real(kreal), intent(in) :: maxCircTol, vortVarTol, lagVarTol
	integer(kint), intent(in) :: limit

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	self%uniformMesh = .FALSE.
	self%vorticityRefine = .TRUE.
	self%maxCircTol = maxCircTol
	self%vortVarTol = vortVarTol
	self%flowMapRefine = .TRUE.
	self%lagVarTol = lagVarTol
	self%tracerRefine = .FALSE.
	self%tracerID = 0
	self%maxMassTol = 0.0_kreal
	self%massVarTol = 0.0_kreal
	self%refinementLimit = limit
end subroutine

subroutine NewPrivateTracer(self, tracerID, maxMassTol, massVarTol, limit)
	type(RemeshSetup), intent(inout) :: self
	integer(kint), intent(in) :: tracerID
	real(kreal), intent(in) :: maxMassTol, massVarTol
	integer(kint), intent(in) :: limit

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	self%uniformMesh = .FALSE.
	self%vorticityRefine = .FALSE.
	self%maxCircTol = 0.0_kreal
	self%vortVarTol = 0.0_kreal
	self%flowMapRefine = .FALSE.
	self%lagVarTol = 0.0_kreal
	self%tracerRefine = .TRUE.
	self%tracerID = tracerID
	self%maxMassTol = maxMassTol
	self%massVarTol = massVarTol
	self%refinementLimit = limit
end subroutine

subroutine NewPrivateNoAMR(self)
	type(RemeshSetup), intent(inout) :: self

	if ( .NOT. logInit ) call InitLogger(log, procRank)

	self%uniformMesh = .TRUE.
	self%vorticityRefine = .FALSE.
	self%maxCircTol = 0.0_kreal
	self%vortVarTol = 0.0_kreal
	self%flowmapRefine = .FALSE.
	self%lagVarTol = 0.0_kreal
	self%tracerRefine = .FALSE.
	self%tracerID = 0
	self%maxMassTol = 0.0_kreal
	self%massVarTol = 0.0_kreal
	self%refinementLimit = 0
end subroutine

subroutine DeletePrivate(self)
	type(RemeshSetup), intent(inout) :: self
	self%uniformMesh = .TRUE.
	self%vorticityRefine = .FALSE.
	self%flowMapRefine = .FALSE.
	self%tracerRefine = .FALSE.
end subroutine


!
!----------------
! Public functions
!----------------
!
subroutine InitialRefinementPrivate( aMesh, remeshData, updateTracerOnMesh, tracerDefinition, updateVorticityOnMesh, vorticityDefinition)
	type(PlaneMesh), intent(inout) :: aMesh
	type(RemeshSetup), intent(in) :: remeshData
	procedure(SetTracerOnMesh) :: updateTracerOnMesh
	type(TracerSetup), intent(in) :: tracerDefinition
	procedure(SetVorticityOnMesh) :: updateVorticityOnMesh
	type(VorticitySetup), intent(in) :: vorticityDefinition
	!
	integer(kint) :: refineCount, spaceLeft, counters(5), j
	integer(kint) :: startIndex, nOldPanels, amrLoopCounter, limit
	logical(klog) :: keepGoing
	logical(klog), allocatable :: refineFlag(:)
	type(Panels), pointer :: aPanels

	if ( remeshData%uniformMesh ) then
		call LogMessage(log, WARNING_LOGGING_LEVEL, logKey, ' InitialRefinement WARNING : uniform mesh = no AMR.')
		return
	endif
	if ( remeshData%refinementLimit == 0 ) then
		call LogMessage(log, WARNING_LOGGING_LEVEL, logKey, ' InitialRefinement WARNING : amr limit = 0.')
		return
	endif

	call StartSection(log, 'InitialRefinement')

	aPanels => aMesh%panels
	allocate(refineFlag(aPanels%N_Max))
	refineFlag = .FALSE.
	keepGoing = .FALSE.
	!
	!	apply refinement criteria
	!
	limit = remeshData%refinementLimit
	startIndex = 1
	counters = 0
	if ( remeshData%vorticityRefine ) then
		!limit = max( remeshData%refinementLimit, limit)
		call FlagPanelsForMaxCirculationRefinement(refineFlag, aMesh, remeshData, startIndex, counters(1))
		call FlagPanelsForVorticityVariationRefinement(refineFlag, aMesh, remeshData, startIndex, counters(2))
	endif
	if ( remeshData%flowMapRefine ) then
		!limit = max( remeshData%refinementLimit, limit)
		call FlagPanelsForFlowMapRefinement(refineFlag, aMesh, remeshData, startIndex, counters(3))
	endif
	if ( remeshData%tracerRefine) then
		!limit = max( remeshData%refinementLimit, limit)
		call FlagPanelsForTracerMassRefinement(refineFlag, aMesh, remeshData, startIndex, counters(4))
		call FlagPanelsForTracerVariationRefinement(refineFlag, aMesh, remeshData, startIndex, counters(5))
	endif

	refineCount = count(refineFlag)
	spaceLeft = aPanels%N_Max - aPanels%N
	if ( refineCount > 0 ) then
		if ( spaceleft/4 > refineCount ) then
			keepGoing = .TRUE.
			amrLoopCounter = 0
			do while (keepGoing)
				amrLoopCounter = amrLoopCounter + 1
				write(logstring,'(A,I2,A,I8,A)') 'AMR loop ', amrLoopCounter, ' : refining ', refineCount,' panels.'
				call LogMessage(log,TRACE_LOGGING_LEVEL,'InitRefine : ',logString)
				!
				! divide flagged panels
				!
				nOldPanels = aPanels%N
				do j = startIndex, nOldPanels
					if ( refineFlag(j) ) then
						call DividePanel(aMesh, j)
						refineFlag(j) = .FALSE.
					endif
				enddo
				!
				!	TO DO : ensure adjacent panels do not differ by more than one level of refinement
				!
				! set data on refined panels
				if ( aMesh%nTracer > 0 ) call UpdateTracerOnMesh(aMesh, tracerDefinition)
				call UpdateVorticityOnMesh(amesh, vorticityDefinition)
				!
				! prevent too much refinement
				!
				if ( amrLoopCounter >= remeshData%refinementLimit ) then
					keepGoing = .FALSE.
					call LogMessage(log, WARNING_LOGGING_LEVEL, 'InitRefine WARNING : ',' refinement limit reached.')
					call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : nActive = ', aPanels%N_Active)
					write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
					call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
					write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
					call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
					write(logstring,'(A, I8, A)') 'flwmap variation criterion triggered ', counters(3), ' times.'
					call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
					write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
					call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
					write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
					call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
				else
					!
					! apply refinement criteria to new panels
					!
					startIndex = nOldPanels + 1
					nOldPanels = aPanels%N
					if ( remeshData%vorticityRefine ) then
						call FlagPanelsForMaxCirculationRefinement(refineFlag, aMesh, remeshData, startIndex, counters(1))
						call FlagPanelsForVorticityVariationRefinement(refineFlag, aMesh, remeshData, startIndex, counters(2))
					endif
					if ( remeshData%flowMapRefine ) then
						call FlagPanelsForFlowMapRefinement(refineFlag, aMesh, remeshData, startIndex, counters(3))
					endif
					if ( remeshData%tracerRefine ) then
						call FlagPanelsForTracerMassRefinement(refineFlag, aMesh, remeshData, startIndex, counters(4))
						call FlagPanelsForTracerVariationRefinement(refineFlag, aMesh, remeshData, startIndex, counters(5))
					endif

					refineCount = count(refineFlag)
					spaceleft = aPanels%N_Max - aPanels%N

					if ( refineCount == 0 ) then
						!
						! stop if refinement is not needed
						!
						call LogMessage(log, TRACE_LOGGING_LEVEL, 'InitRefine : ', 'refinement converged.')
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : nActive = ', aPanels%N_Active)
						write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') 'flwmap variation criterion triggered ', counters(3), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						keepGoing = .FALSE.
					elseif ( spaceLeft / 4 < refineCount ) then
						!
						! stop if not enough memory
						!
						call LogMessage(log, WARNING_LOGGING_LEVEL, 'InitRefine WARNING : ', 'not enough memory to continue AMR.')
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : nActive = ', aPanels%N_Active)
						write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') 'flwmap variation criterion triggered ', counters(3), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'InitRefine : ',trim(logstring))
						keepGoing = .FALSE.
					endif ! stopping criteria met
				endif! below refinement limit
			enddo! while keepgoing
		else ! not enough memory
			call LogMessage(log, WARNING_LOGGING_LEVEL, 'InitRefine WARNING : ', 'not enough memory.')
		endif
	else ! no refinement needed
		call LogMessage(log, TRACE_LOGGING_LEVEL, 'InitRefine : ',' no refinement necessary.')
	endif
	call EndSection(log)
	deallocate(refineFlag)
end subroutine

subroutine LagrangianRemeshToInitialTimePrivate(aMesh, remesh, setVorticity, vorticityDef, setTracer, tracerDef)
	type(PlaneMesh), intent(inout) :: aMesh
	type(RemeshSetup), intent(in) :: remesh
	procedure(SetVorticityOnMesh) :: setVorticity
	type(VorticitySetup), intent(in) :: vorticityDef
	procedure(SetTracerOnMesh) :: setTracer
	type(TracerSetup), intent(in) :: tracerDef
	!
	type(BIVARSetup) :: bivarData
	type(PlaneMesh) :: newMesh
	type(Particles), pointer :: newParticles
	type(Edges), pointer :: newEdges
	type(Panels), pointer :: newPanels
	integer(kint) :: j, amrLoopCounter, counters(5), md
	logical(klog), allocatable :: refineFlag(:)
	logical(klog) :: keepGoing
	integer(kint) :: startIndex, nOldPanels, nOldParticles, refineCount, spaceLeft
	integer(kint), allocatable :: integerWorkParticles(:), integerWorkPanels(:)

	nullify(newParticles)
	nullify(newEdges)
	nullify(newPanels)
	keepGoing = .FALSE.
	amrLoopCounter = 0
	counters = 0

	call LogMessage(log, DEBUG_LOGGING_LEVEL, logkey, ' entering LagrangianRemeshToInitialTime.')

	!
	! build a new uniform mesh
	!
	call New(newMesh, aMesh%initNest, aMesh%amr, aMesh%nTracer)
	call InitializeRectangle(newMesh, MinX(aMesh), MaxX(aMesh), MinY(aMesh), MaxY(aMesh), aMesh%boundaryType)
	newParticles => newMesh%particles
	newEdges => newMesh%edges
	newPanels => newMesh%panels

	!
	! set old mesh as interpolation data source
	!
	call New(bivarData, aMesh)
	allocate(integerWorkParticles( 31*bivarData%n + newParticles%N_Max))
	allocate(integerWorkPanels( 31*bivarData%n + newPanels%N_Max))
	integerWorkParticles = 0
	integerWorkPanels = 0
	!
	! interpolate Lagrangian parameter from old mesh to new mesh
	!
	md = 1
	call IDBVIP( md, bivarData%n, bivarData%x, bivarData%y, bivarData%x0, &
			     newParticles%N, newParticles%x(1,1:newParticles%N), newParticles%x(2,1:newParticles%N), &
			     newParticles%x0(1,1:newParticles%N), integerWorkParticles, bivarData%realWork)
	md = 3
	call IDBVIP( md, bivarData%n, bivarData%x, bivarData%y, bivarData%y0, &
				 newParticles%N, newParticles%x(1,1:newParticles%N), newParticles%x(2,1:newParticles%N), &
				 newParticles%x0(2,1:newParticles%N), integerWorkParticles, bivarData%realWork)
	md = 1
	call IDBVIP( md, bivarData%n, bivarData%x, bivarData%y, bivarData%x0, &
				 newPanels%N, newPanels%x(1,1:newPanels%N), newPanels%x(2,1:newPanels%N), &
				 newPanels%x0(1,1:newPanels%N), integerWorkPanels, bivarData%realWork)
	md = 3
	call IDBVIP( md, bivarData%n, bivarData%x, bivarData%y, bivarData%y0, &
				 newPanels%N, newPanels%x(1,1:newPanels%N), newPanels%x(2,1:newPanels%N), &
				 newPanels%x0(2,1:newPanels%N), integerWorkPanels, bivarData%realWork)

	!
	! set tracer values on new mesh
	!
	if ( aMesh%nTracer > 0 ) call SetTracer(newMesh, tracerDef)
	!
	! set vorticity values on new mesh
	!
	call SetVorticity(newMesh, vorticityDef)

	call LogMessage(log, DEBUG_LOGGING_LEVEL, logkey, ' LagRemeshInitTime : new uniform mesh ready.')

	!
	! AMR
	!
	if ( aMesh%AMR > 0 ) then
		call StartSection(log,'LagrangianRemeshToInitialTime AMR :')
		allocate(refineFlag(newPanels%N_Max))
		refineFlag = .FALSE.
		startIndex = 1
		if ( remesh%vorticityRefine ) then
			call FlagPanelsForMaxCirculationRefinement(refineFlag, newMesh, remesh, startIndex, counters(1) )
			call FlagPanelsForVorticityVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(2) )
		endif
		if ( remesh%flowMapRefine ) then
			call FlagPanelsForFlowMapRefinement(refineFlag, newMesh, remesh, startIndex, counters(3))
		endif
		if ( remesh%tracerRefine ) then
			call FlagPanelsForTracerMassRefinement(refineFlag, newMesh, remesh, startIndex, counters(4))
			call FlagPanelsForTracerVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(5))
		endif

		refineCount = count(refineFlag)
		spaceLeft = newPanels%N_Max - newPanels%N
		if ( refineCount > 0 ) then
			if ( spaceLeft/4 > refineCount ) then
				keepGoing = .TRUE.
				do while (keepGoing)
					amrLoopCounter = amrLoopCounter + 1
					write(logString,'(A,I2,A,I8,A)') 'AMR Loop ', amrLoopCounter, ' : refining ', refineCount, ' panels.'
					call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshInitTime : ',trim(logstring))
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
					! TO DO : ensure adjacent panels differ by no more than 1 level of refinement
					!

					!
					! interpolate Lagrangian parameter to new panels
					!
					md = 1
					call IDBVIP( md, bivarData%n, bivarData%x, bivarData%y, bivarData%x0, &
						 newParticles%N - nOldParticles, newParticles%x(1,nOldParticles+1:newParticles%N), newParticles%x(2,nOldParticles+1:newParticles%N), &
						 newParticles%x0(1,nOldParticles+1:newParticles%N), integerWorkParticles, bivarData%realWork)
					md = 3
					call IDBVIP( md, bivarData%n, bivarData%x, bivarData%y, bivarData%y0, &
						 newParticles%N - nOldParticles, newParticles%x(1,nOldParticles+1:newParticles%N), newParticles%x(2,nOldParticles+1:newParticles%N), &
						 newParticles%x0(2,nOldParticles+1:newParticles%N), integerWorkParticles, bivarData%realWork)
					md = 1
					call IDBVIP( md, bivarData%n, bivarData%x, bivarData%y, bivarData%x0, &
						 newPanels%N - nOldPanels, newPanels%x(1,nOldPanels+1:newPanels%N), newPanels%x(2,nOldPanels+1:newPanels%N), &
						 newPanels%x0(1,nOldPanels+1:newPanels%N), integerWorkPanels, bivarData%realWork)
					md = 3
					call IDBVIP( md, bivarData%n, bivarData%x, bivarData%y, bivarData%y0, &
						 newPanels%N - nOldPanels, newPanels%x(1,nOldPanels+1:newPanels%N), newPanels%x(2,nOldPanels+1:newPanels%N), &
						 newPanels%x0(2,nOldPanels+1:newPanels%N), integerWorkPanels, bivarData%realWork)
					!
					! set data on new mesh
					!
					if ( amesh%nTracer > 0 ) call SetTracer(newMesh, tracerDef)
					call SetVorticity(newMesh, vorticityDef)

					if ( amrLoopCounter >= remesh%refinementLimit ) then
						!
						! prevent too much refinement
						!
						keepGoing = .FALSE.
						call LogMessage(log, WARNING_LOGGING_LEVEL, 'LagRemeshInitTime WARNING : ', 'refinement limit reached.')
						call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshInitTime : nActive = ', newPanels%N_Active)
						write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitTime : ',trim(logstring))
						write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitTime : ',trim(logstring))
						write(logstring,'(A, I8, A)') 'flwmap variation criterion triggered ', counters(3), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitTime : ',trim(logstring))
						write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitTime : ',trim(logstring))
						write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
						call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitTime : ',trim(logstring))
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
							call FlagPanelsForTracerMassRefinement(refineFlag, newMesh, remesh, startIndex, counters(4))
							call FlagPanelsForTracerVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(5))
						endif

						refineCount = count(refineFlag)
						spaceLeft = newPanels%N_Max - newPanels%N

						if ( refineCount == 0 ) then
							keepGoing = .FALSE.
							call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshInitTime : ', ' refinement converged.')
							call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshInitTime : nActive = ', newPanels%N_Active)
							write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitTime : ',trim(logstring))
							write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitTime : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'flwmap variation criterion triggered ', counters(3), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitTime : ',trim(logstring))
							write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitTime : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitTime : ',trim(logstring))
						elseif ( spaceLeft / 4 < refineCount ) then
							keepGoing = .FALSE.
							call LogMessage(log, WARNING_LOGGING_LEVEL, 'LagRemeshInitTime WARNING : ', 'not enough memory to continue AMR.')
							call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshInitTime : nActive = ', newPanels%N_Active)
							write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitTime : ',trim(logstring))
							write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitTime : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'flwmap variation criterion triggered ', counters(3), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitTime : ',trim(logstring))
							write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitTime : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshInitTime : ',trim(logstring))
						endif ! stopping criteria met
					endif ! limit reached
				enddo! while keepGoing
			else ! not enough space left in memory
				call LogMessage(log, WARNING_LOGGING_LEVEL, 'LagRemeshInitTime WARNING : ', 'not enough memory for AMR.')
			endif
		else ! refineCount == 0
			call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshInitTime : ', ' no refinement necessary.')
		endif

		call EndSection(log)
		deallocate(refineFlag)
	endif! AMR
	!
	! replace old mesh with new mesh
	!
	call Copy(aMesh, newMesh)

	!
	! clean up
	!
	deallocate(integerWorkParticles)
	deallocate(integerWorkPanels)
	call Delete(bivarData)
	call Delete(newMesh)
end subroutine

subroutine LagrangianRemeshToReferenceTimePrivate(aMesh, reference, remesh)
	type(PlaneMesh), intent(inout) :: aMesh
	type(BIVARSetup), intent(in) :: reference
	type(RemeshSetup), intent(in) :: remesh
	!
	type(PlaneMesh) :: newMesh
	type(Particles), pointer :: newParticles
	type(Edges), pointer :: newEdges
	type(Panels), pointer :: newPanels
	type(BIVARSetup) :: bivarData
	integer(kint) :: j, k, amrLoopCounter, counters(5), md
	logical(klog), allocatable :: refineFlag(:)
	logical(klog) :: keepGoing
	integer(kint) :: startIndex, nOldPanels, nOldParticles, refineCount, spaceLeft
	integer(kint), allocatable :: integerWorkParticles(:), integerWorkPanels(:)

	nullify(newparticles)
	nullify(newEdges)
	nullify(newPanels)
	amrLoopCounter = 0
	counters = 0
	md = 1
	startIndex = 1
	keepGoing = .FALSE.

	call LogMessage(log, DEBUG_LOGGING_LEVEL, logkey, ' entering LagrangianRemeshToReferenceTime')

	!
	! build a new uniform mesh
	!
	call New(newMesh, aMesh%initNest, aMesh%amr, aMesh%ntracer)
	call InitializeRectangle(newMesh, MinX(aMesh), MaxX(aMesh), MinY(aMesh), MaxY(aMesh), aMesh%boundaryType)
	call LogMessage(log, DEBUG_LOGGING_LEVEL, logkey, ' base mesh returned.')
	newParticles => newMesh%particles
	newEdges => newMesh%edges
	newPanels => newMesh%panels



	!
	! setup old mesh as interpolation source
	!
	call New(bivarData, amesh)
	allocate(integerWorkParticles( 31*bivarData%n + newParticles%N_Max))
	allocate(integerWorkPanels( 31*bivarData%n + newPanels%N_Max))
	integerWorkParticles = 0
	integerWorkPanels = 0
	call LogMessage(log, DEBUG_LOGGING_LEVEL, logkey, ' alpha source data ready.')

	!
	! interpolate Lagrangian parameter from old mesh to new mesh
	!
	md = 1
	call IDBVIP( md, bivarData%n, bivarData%x, bivarData%y, bivarData%x0, &
				 newParticles%N, newParticles%x(1,1:newParticles%N), newParticles%x(2,1:newParticles%N), &
				 newParticles%x0(1,1:newParticles%N), integerWorkParticles, bivarData%realWork)
	md = 3
	call IDBVIP( md, bivarData%n, bivarData%x, bivarData%y, bivarData%y0, &
				 newParticles%N, newParticles%x(1,1:newParticles%N), newParticles%x(2,1:newParticles%N), &
				 newParticles%x0(2,1:newParticles%N), integerWorkParticles, bivarData%realWork)
	md = 1
	call IDBVIP( md, bivarData%N, bivarData%x, bivarData%y, bivarData%x0, &
				 newPanels%N, newPanels%x(1,1:newPanels%n), newPanels%x(2,1:newPanels%n), &
				 newPanels%x0(1,1:newPanels%N), integerWorkPanels, bivarData%realWork)
	md = 3
	call IDBVIP( md, bivarData%N, bivarData%x, bivarData%y, bivarData%y0, &
				 newPanels%N, newPanels%x(1,1:newPanels%n), newPanels%x(2,1:newPanels%n), &
				 newPanels%x0(2,1:newPanels%N), integerWorkPanels, bivarData%realWork)

	call LogMessage(log, DEBUG_LOGGING_LEVEL, logkey, ' base mesh interpolation complete.')


	!
	! set tracer values on new mesh
	!
	do k = 1, aMesh%nTracer
		call AssignTracer(newParticles%tracer(1:newParticles%N,k), newParticles%x0(:,1:newParticles%N), reference, k)
		call AssignTracer(newPanels%tracer(1:newPanels%N, k), newPanels%x0(:,1:newPanels%N), reference, k)
	enddo
	!
	! set vorticity values on new mesh
	!
	call AssignVorticity(newParticles%relVort(1:newParticles%N), newParticles%x0(:,1:newParticles%N), reference)
	call AssignVorticity(newPanels%relVort(1:newPanels%N), newPanels%x0(:,1:newPanels%N), reference)

	call LogMessage(log, DEBUG_LOGGING_LEVEL, logkey, ' LagRemeshReference : new uniform mesh ready.')

	if ( aMesh%AMR > 0 ) then
		call StartSection(log,'LagrangianRemeshToReference AMR :')
			allocate(refineFlag(newPanels%N_Max))
			refineFlag = .FALSE.
			if ( remesh%vorticityRefine ) then
				call FlagPanelsForMaxCirculationRefinement(refineFlag, newMesh, remesh, startIndex, counters(1))
				call FlagPanelsForVorticityVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(2))
			endif
			if ( remesh%flowMapRefine ) then
				call FlagPanelsForFlowMapRefinement(refineFlag, newMesh, remesh, startIndex, counters(3))
			endif
			if ( remesh%tracerRefine ) then
				call FlagPanelsForTracerMassRefinement(refineFlag, newMesh, remesh, startIndex, counters(4))
				call FlagPanelsForTracerVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(5))
			endif

			refineCount = count(refineFlag)
			spaceLeft = newPanels%N_Max - newPanels%N

			if ( refineCount > 0 ) then
				if ( spaceLeft / 4 > refineCount ) then
					keepGoing = .TRUE.
					do while (keepGoing)
						amrLoopCounter = amrLoopCounter + 1
						write(logString,'(A,I2,A,I8,A)') 'AMR loop ', amrLoopCounter, ' : refining ', refineCount,' panels.'
						call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshReference : ', trim(logstring))
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
						! TO DO : ensure adjacent panels differ by no more than 1 level of refinement
						!

						!
						! interpolate Lagrangian parameter to new panels
						!
						md = 1
						call IDBVIP(md, bivarData%n, bivarData%x, bivarData%y, bivarData%x0, &
								    newParticles%N - nOldParticles, &
								    newParticles%x(1,nOldParticles+1:newParticles%N), newParticles%x(2,nOldParticles+1:newParticles%N), &
								    newParticles%x0(1, nOldParticles+1:newParticles%N), integerWorkParticles, bivarData%realwork)
						md = 3
						call IDBVIP(md, bivarData%n, bivarData%x, bivarData%y, bivarData%y0, &
								    newParticles%N - nOldParticles, &
								    newParticles%x(1,nOldParticles+1:newParticles%N), newParticles%x(2,nOldParticles+1:newParticles%N), &
								    newParticles%x0(2, nOldParticles+1:newParticles%N), integerWorkParticles, bivarData%realwork)
						md = 1
						call IDBVIP(md, bivarData%n, bivarData%x, bivarData%y, bivarData%x0, &
						 			newPanels%N - nOldPanels, &
						 			newPanels%x(1,nOldPanels+1:newPanels%N), newPanels%x(2,nOldPanels+1:newPanels%N), &
						 			newPanels%x0(1,nOldPanels+1:newPanels%N), integerWorkPanels, bivarData%realWork)
						md = 3
						call IDBVIP(md, bivarData%n, bivarData%x, bivarData%y, bivarData%y0, &
						 			newPanels%N - nOldPanels, &
						 			newPanels%x(1,nOldPanels+1:newPanels%N), newPanels%x(2,nOldPanels+1:newPanels%N), &
						 			newPanels%x0(2,nOldPanels+1:newPanels%N), integerWorkPanels, bivarData%realWork)
						!
						! set data on new particles, panels
						!
						do k = 1, aMesh%nTracer
							call AssignTracer(newParticles%tracer(nOldParticles+1:newParticles%N,k), &
											  newParticles%x0(:,nOldParticles+1:newParticles%N), reference, k)
							call AssignTracer(newPanels%tracer(nOldPanels+1:newPanels%N,k), &
											  newPanels%x0(:,nOldPanels+1:newPanels%N), reference, k)
						enddo
						call AssignVorticity(newParticles%relVort(nOldParticles+1:newParticles%N), &
											 newParticles%x0(:, nOldParticles+1:newParticles%N), reference)
						call AssignVorticity(newPanels%relVort(nOldPanels+1:newPanels%N), &
											 newPanels%x0(:,nOldPanels+1:newPanels%N), reference)

						!
						! prevent too much refinement
						!
						if ( amrLoopCounter >= remesh%refinementLimit ) then
							keepGoing = .FALSE.
							call LogMessage(log, WARNING_LOGGING_LEVEL, 'LagRemeshReference WARNING : ', 'refinement limit reached.')
							call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshReference N_Active = ', newPanels%N_Active)
							write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshReference : ',trim(logstring))
							write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshReference : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'flwmap variation criterion triggered ', counters(3), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshReference : ',trim(logstring))
							write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshReference : ',trim(logstring))
							write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
							call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshReference : ',trim(logstring))
						else
							!
							! apply refinement criteria to new panels
							!
							startIndex = nOldPanels+1
							nOldPanels = newPanels%N
							if ( remesh%vorticityRefine ) then
								call FlagPanelsForMaxCirculationRefinement(refineFlag, newMesh, remesh, startIndex, counters(1) )
								call FlagPanelsForVorticityVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(2) )
							endif
							if ( remesh%flowMapRefine ) then
								call FlagPanelsForFlowMapRefinement(refineFlag, newMesh, remesh, startIndex, counters(3))
							endif
							if ( remesh%tracerRefine ) then
								call FlagPanelsForTracerMassRefinement(refineFlag, newMesh, remesh, startIndex, counters(4))
								call FlagPanelsForTracerVariationRefinement(refineFlag, newMesh, remesh, startIndex, counters(5))
							endif

							refineCount = count(refineFlag)
							spaceLeft = newPanels%N_Max - newPanels%N

							!
							! check stopping criteria
							!
							if ( refineCount == 0 ) then
								keepGoing = .FALSE.
								call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshReference : ', ' refinement converged.')
								call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshReference N_Active = ', newPanels%N_Active)
								write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
								call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshReference : ',trim(logstring))
								write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
								call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshReference : ',trim(logstring))
								write(logstring,'(A, I8, A)') 'flwmap variation criterion triggered ', counters(3), ' times.'
								call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshReference : ',trim(logstring))
								write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
								call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshReference : ',trim(logstring))
								write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
								call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshReference : ',trim(logstring))
							elseif ( spaceLeft / 4 < refineCount ) then
								keepGoing = .FALSE.
								call LogMessage(log, WARNING_LOGGING_LEVEL, 'LagRemeshReference WARNING: ', ' not enough memory to continue AMR.')
								call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshReference N_Active = ', newPanels%N_Active)
								write(logstring,'(A, I8, A)') ' max circulation criterion triggered ', counters(1), ' times.'
								call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshReference : ',trim(logstring))
								write(logstring,'(A, I8, A)') ' vort. variation criterion triggered ', counters(2), ' times.'
								call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshReference : ',trim(logstring))
								write(logstring,'(A, I8, A)') 'flwmap variation criterion triggered ', counters(3), ' times.'
								call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshReference : ',trim(logstring))
								write(logstring,'(A, I8, A)') '       tracermax criterion triggered ', counters(4), ' times.'
								call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshReference : ',trim(logstring))
								write(logstring,'(A, I8, A)') 'tracer variation criterion triggered ', counters(5), ' times.'
								call LogMessage(log, TRACE_LOGGING_LEVEL,'LagRemeshReference : ',trim(logstring))
							endif ! stopping criteria met
						endif! refinement limit reached
					enddo!while keepgoing
				else
					call LogMessage(log, WARNING_LOGGING_LEVEL,'LagRemeshReference : ',' not enough memory for AMR.')
				endif ! spaceleft
			else ! refinecount == 0
				call LogMessage(log, TRACE_LOGGING_LEVEL, 'LagRemeshReference : ',' no refinement necessary.')
			endif
		call EndSection(log)
		deallocate(refineFlag)
	endif!AMR
	!
	!	replace old mesh with new mesh
	!
	call Copy(aMesh, newMesh)

	!
	! clean up
	!
	deallocate(integerWorkParticles)
	deallocate(integerWorkPanels)
	call Delete(bivarData)
	call Delete(newMesh)
end subroutine


!
!----------------
! Module methods : module- or type-specific private functions
!----------------
!
subroutine FlagPanelsForMaxCirculationRefinement(refineFlag, aMesh, remeshData, startIndex, counter)
	logical(klog), intent(inout) :: refineFlag(:)
	type(PlaneMesh), intent(in) :: aMesh
	type(RemeshSetup), intent(in) :: remeshData
	integer(kint), intent(in) :: startIndex
	integer(kint), intent(inout) :: counter
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j

	aPanels => aMesh%panels
	aParticles => aMesh%particles

	do j = startIndex, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			if ( abs(aPanels%relVort(j)) * aPanels%area(j) > remeshData%maxCircTol ) then
				refineFlag(j) = .TRUE.
				counter = counter + 1
			endif
		endif
	enddo
end subroutine

subroutine FlagPanelsForVorticityVariationRefinement(refineFlag, aMesh, remeshData, startIndex, counter)
	logical(klog), intent(inout) :: refineFlag(:)
	type(PlaneMesh), intent(in) :: aMesh
	type(RemeshSetup), intent(in) :: remeshData
	integer(kint), intent(in) :: startIndex
	integer(kint), intent(inout) :: counter
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: edgeList(8), vertList(8), nverts, j, k
	real(kreal) :: maxVort, minVort, vortVar

	aPanels => aMesh%panels
	aParticles => aMesh%particles

	do j = startIndex, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			maxVort = aPanels%relVort(j)
			minVort = aPanels%relVort(j)
			call CCWEdgesAndParticlesAroundPanel(edgeList, vertLIst, nVerts, aMesh, j)
			do k = 1, nVerts
				if ( aParticles%relVort(vertList(k)) > maxVort ) maxVort = aParticles%relVort(vertList(k))
				if ( aParticles%relVort(vertList(k)) < minVort ) minVort = aParticles%relVort(vertList(k))
			enddo
			vortVar = maxVort - minVort
			if ( vortVar > remeshData%vortVarTol ) then
				refineFlag(j) = .TRUE.
				counter = counter + 1
			endif
		endif
	enddo
end subroutine

subroutine FlagPanelsForFlowMapRefinement(refineFlag, aMesh, remeshData, startIndex, counter)
	logical(klog), intent(inout) :: refineFlag(:)
	type(PlaneMesh), intent(in) :: aMesh
	type(RemeshSetup), intent(in) :: remeshData
	integer(kint), intent(in) :: startIndex
	integer(kint), intent(inout) :: counter
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: edgeList(8), vertList(8), nverts, j, k
	real(kreal) :: maxX0(2), minX0(2), lagVar

	aPanels => aMesh%panels
	aParticles => aMesh%particles

	do j = startIndex, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			maxX0 = aPanels%x0(:,j)
			minX0 = aPanels%x0(:,j)
			call CCWEdgesAndParticlesAroundPanel(edgeList, vertLIst, nVerts, aMesh, j)
			do k = 1, nVerts
				if ( aParticles%x0(1,vertList(k)) > maxX0(1) ) maxX0(1) = aParticles%x0(1,vertList(k))
				if ( aParticles%x0(1,vertList(k)) < minX0(1) ) minX0(1) = aParticles%x0(1,vertList(k))
				if ( aParticles%x0(2,vertList(k)) > maxX0(2) ) maxX0(2) = aParticles%x0(2,vertList(k))
				if ( aParticles%x0(2,vertList(k)) < minX0(2) ) minx0(2) = aParticles%x0(2,vertlist(k))
			enddo
			lagVar = sum( maxX0 - minx0)
			if ( lagVar > remeshData%lagVarTol ) then
				refineFlag(j) = .TRUE.
				counter = counter + 1
			endif
		endif
	enddo
end subroutine

subroutine FlagPanelsForTracerMassRefinement(refineFlag, aMesh, remeshData, startIndex, counter)
	logical(klog), intent(inout) :: refineFlag(:)
	type(PlaneMesh), intent(in) :: aMesh
	type(RemeshSetup), intent(in) :: remeshData
	integer(kint), intent(in) :: startIndex
	integer(kint), intent(inout) :: counter
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j

	aPanels => aMesh%panels
	aParticles => aMesh%particles

	do j = startIndex, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			if ( aPanels%tracer(j,remeshData%tracerID) * aPanels%area(j) > remeshData%maxMassTol ) then
				refineFlag(j) = .TRUE.
				counter = counter + 1
			endif
		endif
	enddo
end subroutine

subroutine FlagPanelsForTracerVariationRefinement(refineFlag, aMesh, remeshData, startIndex, counter)
	logical(klog), intent(inout) :: refineFlag(:)
	type(PlaneMesh), intent(in) :: aMesh
	type(RemeshSetup), intent(in) :: remeshData
	integer(kint), intent(in) :: startIndex
	integer(kint), intent(inout) :: counter
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: edgeList(8), vertList(8), nverts, j, k
	real(kreal) :: maxTracer, minTracer, tracerVar

	aPanels => aMesh%panels
	aParticles => aMesh%particles

	do j = startIndex, aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			maxTracer = aPanels%tracer(j, remeshData%tracerID)
			minTracer = aPanels%tracer(j, remeshData%tracerID)
			call CCWEdgesAndParticlesAroundPanel(edgeList, vertLIst, nVerts, aMesh, j)
			do k = 1, nVerts
				if ( aParticles%tracer(vertList(k), remeshData%tracerID) > maxTracer ) &
							maxTracer = aParticles%tracer(vertList(k), remeshData%tracerID)
				if ( aParticles%tracer(vertList(k), remeshData%tracerID) < minTracer ) &
							minTracer = aParticles%tracer(vertList(k), remeshData%tracerID)
			enddo
			tracerVar = maxTracer - minTracer
			if ( tracerVar > remeshData%massVarTol ) then
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

