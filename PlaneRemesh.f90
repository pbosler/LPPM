module PlaneRemeshModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
!----------------
! USAGE :  This module provides methods for interpolating and remeshing a planar mesh.
!----------------
use NumberKindsModule
use LoggerModule
use PlaneGeomModule
use ParticlesModule
use EdgesModule
use PanelsModule
use PlaneMeshModule
use PlaneVorticityModule
use PlaneTracerModule

implicit none

private 
public RemeshSetup
public New, Delete
public InitialRefinement!, LagrangianRemesh


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
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
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
end subroutine


!
!----------------
! Public functions
!----------------
!
subroutine InitialRefinement( aMesh, remeshData, updateTracerOnMesh, tracerDefinition, updateVorticityOnMesh, vorticityDefinition)
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

end module

