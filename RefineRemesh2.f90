module RefineRemeshModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!	Defines the mesh data structure used by icosahedral triangle and cubed
!	sphere Lagrangian meshes of the sphere.
!
!	Bosler, P.A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis; the University of Michigan, 2013.
!
!----------------
use NumberKindsModule
use SphereGeomModule
use LoggerModule
use ParticlesModule
use EdgesModule
use PanelsModule
use STRIPACKInterfaceModule
use SSRFPACKInterfaceModule
use SphereMeshModule
use TracerSetupModule
use BVESetupModule

implicit none

include 'mpif.h'

private
public RefinementSetup
public New, Delete
!public LagrangianRemesh, DirectRemesh
public InitialRefinement
public NULL_REFINE, TRACER_REFINE, RELVORT_REFINE, FLOWMAP_REFINE

!
!----------------
! Types and module constants
!----------------
!
type RefinementSetup
	real(kreal) :: maxTol			! tolerance for extrema
	real(kreal) :: varTol			! tolerance for variation
	integer(kint) :: type	! identifier for physical data field
	integer(kint) :: tracerID
	integer(kint) :: limit
end type

integer(kint), parameter :: NULL_REFINE = 70, &
							TRACER_REFINE = 71, &
							RELVORT_REFINE = 72, &
							FLOWMAP_REFINE = 73
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
	subroutine SetVorticityOnMesh(genMesh,genVort)
		use SphereMeshModule
		use BVESetupModule
		implicit none
		type(SphereMesh), intent(inout) :: genMesh
		type(BVESetup), intent(in) :: genVort
	end subroutine
end interface
!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'RefineRemesh'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=128) :: logString
character(len=24) :: formatString
contains
!
!----------------
! Standard methods : Constructor / Destructor
!----------------
!
subroutine NewPrivate(self, limit, maxTol, varTol, type, tracerID)
	type(RefinementSetup), intent(out) :: self
	real(kreal), intent(in) :: maxTol, varTol, limit
	integer(kint), intent(in) :: type
	integer(kint), intent(in), optional :: tracerID
	
	if ( .NOT. logInit ) call InitLogger(log, procRank)
	
	if ( type < NULL_REFINE .OR. type > FLOWMAP_REFINE ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,' invalid refinement type.')
		return
	elseif ( type == TRACER_REFINE) then
		if ( .NOT. present(tracerID) ) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,' must specify which tracer to refine.')
			return
		else 
			self%tracerID = tracerID	
		endif	
	endif
	self%limit = limit
	self%maxTol = maxTOl
	self%varTol = varTol
	self%type = type
end subroutine

subroutine DeletePrivate(self)
	type(RefinementSetup), intent(inout) :: self
	self%type = NULL_REFINE
end subroutine
!
!----------------
! Public functions
!----------------
!
subroutine InitialRefinement(aMesh, refineTracer, updateTracerOnMesh, tracerDef, &
									refineRelVort, updateVorticityOnMesh, vortDef)
	type(SphereMesh), intent(inout) :: aMesh
	type(RefinementSetup), intent(in) :: refineTracer
	procedure(SetTracerOnMesh) :: updateTracerOnMesh
	type(TracerSetup), intent(in) :: tracerDef
	type(RefinementSetup), intent(in) :: refineRelVort
	procedure(SetVorticityOnMesh) :: updateVorticityOnMesh
	type(BVESetup), intent(in) :: vortDef
	! local variables
	integer(kint) :: refineCount, spaceLeft, j, counter1, counter2
	integer(kint) :: startIndex, nOldPanels, amrLoopCounter, limit
	logical(klog) :: keepGoing
	type(Panels), pointer :: aPanels
	logical(klog), allocatable :: refineFlag(:)
	
	! check for invalid states
	if ( refineTracer%type == TRACER_REFINE .AND. ( GetNTracer(aMesh%panels) < refineTracer%tracerID ) ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'InitialRefinement ERROR : ','invalid tracer number.')
		return
	endif
	if ( refineTracer%type /= NULL_REFINE .AND. refineTracer%type /= TRACER_REFINE ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'InitialRefinement ERROR : ','invalid tracer refinement type.')
		return
	endif
	if ( refineRelVort%type /= NULL_REFINE .AND. refineRelVort%type /= RELVORT_REFINE ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'InitialRefinement ERROR : ','invalid relVort refinement type.')
		return
	endif
	if ( refineRelVort%type == NULL_REFINE .AND. refineTracer%type == NULL_REFINE ) then
		call LogMessage(log,WARNING_LOGGING_LEVEL,'InitialRefinement WARNING : ','NULL refinement data.')
		return
	endif

	aPanels => aMesh%panels
	allocate(refineFlag(aPanels%N_Max))
	refineFlag = .FALSE.
	keepGoing = .FALSE.
	!
	! 	Apply refinement criteria
	!
	limit = 0
	startIndex = 1
	if ( refineTracer%type /= NULL_REFINE ) then
		limit = max(limit,refineTracer%limit)
		call FlagPanelsForTracerMaxRefinement(refineFlag,aMesh,refineTracer,startIndex)
		counter1 = count(refineFlag)
		call FlagPanelsForTracerVariationRefinement(refineFlag,aMesh,refineTracer,startIndex)
		counter2 = count(refineFlag) - counter1
		write(formatString,'(A)') '(A,I8,A)'
		write(logString,formatString) 'tracerMax criterion triggered ', counter1, ' times.'
		call LogMessage(log,TRACE_LOGGING_LEVEL,'InitRefine : ',logString)
		write(logString,formatString) 'tracerVar criterion triggered ', counter2, ' times.'
		call LogMessage(log,TRACE_LOGGING_LEVEL,'InitRefine : ',logString)
	endif	
	if ( refineRelVort%type /= NULL_REFINE) then
		limit = max(limit,refineRelVort%limit)
		counter1 = count(refineFlag)
		call FlagPanelsForCirculationRefinement(refineFlag,aMesh,refineRelVort,startIndex)
		counter1 = count(refineFlag) - counter1
		call FlagPanelsForRelVortVariationRefinement(refineFlag,aMesh,refineRelVort,startIndex)
		counter2 = count(refineFlag) - counter2
		write(formatString,'(A)') '(A,I8,A)'
		write(logString,formatString) 'circMax criterion triggered ', counter1, ' times.'
		call LogMessage(log,TRACE_LOGGING_LEVEL,'InitRefine : ',logString)
		write(logString,formatString) 'relVortVar criterion triggered ', counter2, ' times.'
		call LogMessage(log,TRACE_LOGGING_LEVEL,'InitRefine : ',logString)
	endif
	
	refineCount = count(refineFlag)
	spaceLeft = aPanels%N_Max - aPanels%N
	
	!
	!	exit if refinement is not needed
	!
	if ( refineCount == 0 ) then
		call LogMessage(log,TRACE_LOGGING_LEVEL,'InitRefine : ','no refinement necessary.')
		deallocate(refineFlag)
		return
	endif
		
	!
	!	check for memory, exit if insufficient
	!
	if ( spaceLeft/4 < refineCount ) then
		call LogMessage(log,WARNING_LOGGING_LEVEL,'InitRefine : ',' insufficient memory for AMR.')
		deallocate(refineFlag)
		return		
	endif
	
	keepGoing = .TRUE.
	amrLoopCounter = 0
	
	do while (keepGoing)
		amrLoopCounter = amrLoopCounter + 1
		write(logString,formatString) 'AMR loop ',amrLoopCounter,' : refining ',refineCount,' panels.'
		call LogMessage(log,TRACE_LOGGING_LEVEL,'InitRefine : ',logString)
		!
		!	Divide flagged panels
		!
		nOldPanels = aPanels%N
		do j=startIndex, aPanels%N
			if ( refineFlag(j) ) then
				call DividePanel(aMesh,j)
				refineFlag(j) = .FALSE.
			endif
		enddo
		!
		!	Ensure adjacent panels differ by no more than one mesh level
		!
		call FlagPanelsAtRefinementBoundaries(refineFlag,aMesh)
		do j=startIndex,aPanels%N
			if ( refineFlag(j) ) then
				call DividePanel(aMesh,j)
				refineFlag(j) = .FALSE.
			endif
		enddo
		!
		!	Set data on refined mesh
		!
		call UpdateTracerOnMesh(aMesh,tracerDef)
		call UpdateVorticityOnMesh(aMesh,vortDef)
		!
		!	Prevent too much refinement
		!
		if ( amrLoopCounter + 1 >= limit ) then
			keepGoing = .FALSE.
			call LogMessage(log,WARNING_LOGGING_LEVEL,'InitRefine WARNING : ',' refinement limit reached.')
		endif
		!
		!	Apply refinement criteria
		!
		startIndex = nOldPanels+1
		nOldPanels = aPanels%N
		if ( refineTracer%type /= NULL_REFINE ) then
			call FlagPanelsForTracerMaxRefinement(refineFlag,aMesh,refineTracer,startIndex)
			counter1 = count(refineFlag)
			call FlagPanelsForTracerVariationRefinement(refineFlag,aMesh,refineTracer,startIndex)
			counter2 = count(refineFlag) - counter1
			write(formatString,'(A)') '(A,I8,A)'
			write(logString,formatString) 'tracerMax criterion triggered ', counter1, ' times.'
			call LogMessage(log,TRACE_LOGGING_LEVEL,'InitRefine : ',logString)
			write(logString,formatString) 'tracerVar criterion triggered ', counter2, ' times.'
			call LogMessage(log,TRACE_LOGGING_LEVEL,'InitRefine : ',logString)
		endif	
		if ( refineRelVort%type /= NULL_REFINE) then
			counter1 = count(refineFlag)
			call FlagPanelsForCirculationRefinement(refineFlag,aMesh,refineRelVort,startIndex)
			counter1 = count(refineFlag) - counter1
			call FlagPanelsForRelVortVariationRefinement(refineFlag,aMesh,refineRelVort,startIndex)
			counter2 = count(refineFlag) - counter2
			write(formatString,'(A)') '(A,I8,A)'
			write(logString,formatString) 'circMax criterion triggered ', counter1, ' times.'
			call LogMessage(log,TRACE_LOGGING_LEVEL,'InitRefine : ',logString)
			write(logString,formatString) 'relVortVar criterion triggered ', counter2, ' times.'
			call LogMessage(log,TRACE_LOGGING_LEVEL,'InitRefine : ',logString)
		endif
		!
		!	exit if refinement is not needed
		!
		if ( refineCount == 0 ) then
			call LogMessage(log,TRACE_LOGGING_LEVEL,'InitRefine : ','refinement converged.')
			keepGoing = .FALSE.
		endif
		
		!
		!	check for memory, exit if insufficient
		!
		if ( spaceLeft/4 < refineCount ) then
			call LogMessage(log,WARNING_LOGGING_LEVEL,'InitRefine : ',' insufficient memory for AMR.')
			keepGoing = .FALSE.
		endif
	enddo	
	 deallocate(refineFlag)
end subroutine

!
!----------------
! Module methods : module- or type-specific private functions
!----------------
!
subroutine FlagPanelsForTracerMaxRefinement(refineFlag,aMesh,refineTracer,startIndex)
	logical(klog), intent(inout) :: refineFlag(:)
	type(SphereMesh), intent(in) :: aMesh
	type(RefinementSetup), intent(in) :: refineTracer
	integer(kint), intent(in) :: startIndex
	! local variables
	type(Panels), pointer :: aPanels
	integer(kint) :: j
	
	if ( refineTracer%type /= TRACER_REFINE ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'FlagPanelsTracer ERROR :',' invalid refinement type.')
		return
	endif
	aPanels => aMesh%panels	
	do j=startIndex,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			if ( abs(aPanels%tracer(j,refineTracer%tracerID))*aPanels%area(j) > refineTracer%maxTol ) refineFlag(j) = .TRUE.
		endif
	enddo
end subroutine

subroutine FlagPanelsForCirculationRefinement(refineFlag,aMesh,refineRelVort,startIndex)
	logical(klog), intent(inout) :: refineFlag(:)
	type(SphereMesh), intent(in) :: aMesh
	type(RefinementSetup), intent(in) :: refineRelVort
	integer(kint), intent(in) :: startIndex
	! local variables
	type(Panels), pointer :: aPanels
	integer(kint) :: j

	if ( refineRelVort%type /= RELVORT_REFINE ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'FlagPanelsCirc ERROR :',' invalid refinement type.')
		return
	endif	
	aPanels => aMesh%panels
	do j=startIndex,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			if ( abs(aPanels%relVort(j))*aPanels%area(j) > refineRelVort%maxTol ) refineFlag(j) = .TRUE.
		endif
	enddo
end subroutine


subroutine FlagPanelsForTracerVariationRefinement(refineFlag,aMesh,refineTracer,startIndex)
	logical(klog), intent(inout) :: refineFlag(:)
	type(SphereMesh), intent(in) :: aMesh
	type(RefinementSetup), intent(in) :: refineTracer
	integer(kint), intent(in) :: startIndex
	! local variables
	type(Panels), pointer :: aPanels
	type(Particles), pointer :: aParticles
	integer(kint) :: edgeList(8), vertList(8), nVerts
	integer(kint) :: j, k
	real(kreal) :: maxTracer, minTracer, tracerVar

	if ( refineTracer%type /= TRACER_REFINE ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'FlagPanelsTracer ERROR :',' invalid refinement type.')
		return
	endif	
	aParticles => aMesh%particles
	aPanels => aMesh%panels
	
	do j=startIndex,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			maxTracer = aPanels%tracer(j,refineTracer%tracerID)
			minTracer = maxTracer
			call CCWEdgesAndParticlesAroundPanel(edgeList,vertList,nVerts,aMesh,j)
			do k=1,nVerts
				if ( aParticles%tracer(vertList(k),refineTracer%tracerID) > maxTracer) &
					maxTracer = aParticles%tracer(vertList(k),refineTracer%tracerID)
				if ( aParticles%tracer(vertList(k),refineTracer%tracerID) < minTracer) &
					minTracer = aParticles%tracer(vertList(k),refineTracer%tracerID)	
			enddo
			tracerVar = maxTracer - minTracer
			if ( tracerVar > refineTracer%varTol ) refineFlag(j) = .TRUE.
		endif
	enddo
end subroutine

subroutine FlagPanelsForRelVortVariationRefinement(refineFlag,aMesh,refineRelVort,startIndex)
	logical(klog), intent(inout) :: refineFlag(:)
	type(SphereMesh), intent(in) :: aMesh
	type(RefinementSetup), intent(in) :: refineRelVort
	integer(kint), intent(in) :: startIndex
	! local variables
	type(Panels), pointer :: aPanels
	type(Particles), pointer :: aParticles
	integer(kint) :: edgeList(8), vertList(8), nVerts
	integer(kint) :: j, k
	real(kreal) :: maxRelvort, minRelvort, relVortVar

	if ( refineRelVOrt%type /= RELVORT_REFINE ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'FlagPanelsRelVortVar ERROR :',' invalid refinement type.')
		return
	endif	
	aParticles => aMesh%particles
	aPanels => aMesh%panels
	
	do j=startIndex,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			maxrelVort = aPanels%relVort(j)
			minRelVort = maxRelVort
			call CCWEdgesAndParticlesAroundPanel(edgeList,vertList,nVerts,aMesh,j)
			do k=1,nVerts
				if ( aParticles%relVort(vertList(k)) > maxrelVort) &
					maxrelVort = aParticles%relVort(vertList(k))
				if ( aParticles%relVort(vertList(k)) < minrelVort) &
					minRelVort = aParticles%relVort(vertList(k))			
			enddo
			relVortVar = maxRelVort - minRelVort
			if ( relVortVar > refineRelVort%varTol ) refineFlag(j) = .TRUE.
		endif
	enddo
end subroutine

subroutine InitLogger(aLog,rank)
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
