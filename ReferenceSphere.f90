module ReferenceSphereModule
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
use STRIPACKInterfaceModule
use SSRFPACKInterfaceModule
use RefineRemeshModule


implicit none
private
public ReferenceSphere
public New, Delete
public LagrangianRemesh

type ReferenceSphere
	type(STRIPACKData) :: delTri
	type(SSRFPACKData), pointer :: absVortSource => null(), &
								   tracerSource => null()
	integer(kint) :: startSearch = 1
end type

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

interface LagrangianRemesh
	module procedure LagrangianRemeshPrivate
end interface


!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'ReferenceSphere'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=28) :: formatString
character(len=128) :: logString

contains
!
!----------------
! Standard methods : Constructor / Destructor, Copy
!----------------
!
subroutine NewPrivate(self, oldSphere)
	type(ReferenceSphere), intent(out) :: self
	type(SphereMesh), intent(in) :: oldSphere
	!
	integer(kint) :: i, j, k, errCode
	type(Panels), pointer :: oldPanels
	type(Particles), pointer :: oldParticles
	real(kreal) :: dSig
	logical(klog) :: false

	if ( .NOT. logInit) call InitLogger(log,procRank)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,' creating new reference sphere.')

	call New(self%delTri,oldSphere)

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,' reference sphere delaunay graph ready.')
	oldParticles => oldSphere%particles
	oldPanels => oldSphere%panels

	if ( oldSphere%problemKind == BVE_SOLVER) then
		false = .FALSE.
		allocate(self%absVortSource)
		call New(self%absVortSource, self%delTri, false)
		call SetSourceAbsVort(self%absVortSource, self%delTri)
	endif

	if ( oldSphere%nTracer > 0 ) then
		allocate(self%tracerSource)
		call New(self%tracerSource, self%delTri, oldSphere%nTracer)
		call SetSourceTracer(self%tracerSource, self%delTri)
	endif

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,' reference sphere ready.')

end subroutine

subroutine DeletePrivate(self)
	type(ReferenceSphere), intent(inout) :: self

	call Delete(self%delTri)
	if ( associated(self%absVortSource) ) then
		call Delete(self%absVortSource)
		deallocate(self%absVortSource)
	endif
	if ( associated(self%tracerSource) ) then
		call Delete(self%tracerSource)
		deallocate(self%tracerSource)
	endif

end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine LagrangianRemeshPrivate(aMesh, reference, vortRefine, tracerRefine, flowMapRefine)
	type(SphereMesh), intent(inout) :: aMesh
	type(ReferenceSphere), intent(inout) :: reference
	type(RefinementSetup), intent(in) :: vortRefine
	type(RefinementSetup), intent(in) :: tracerRefine
	type(RefinementSetup), intent(in) :: flowMapRefine
	!
	type(STRIPACKData) :: delTri
	type(SSRFPACKData) :: lagSource
	logical(klog) :: vectorInterp
	type(SphereMesh) :: newMesh
	type(Particles), pointer :: newParticles
	type(Edges), pointer :: newEdges
	type(Panels), pointer :: newPanels
	integer(kint) :: j, k, amrLoopCounter, counter1, counter2
	logical(klog), allocatable :: refineFlag(:)
	logical(klog) :: keepGoing, refineTracer, refineVort, refineFlowMap
	integer(kint) :: startIndex, nOldPanels, nOldParticles, refineCount, spaceLeft, limit

	nullify(newParticles)
	nullify(newEdges)
	nullify(newPanels)
	vectorInterp = .TRUE.
	refineFlowMap = .FALSE.
	refineTracer = .FALSE.
	refineVort = .FALSE.

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,' entering Lagrangian remesh.')
	!
	!	determine what types of AMR to use
	!
	if ( tracerRefine%type == TRACER_REFINE .AND. &
		GetNTracer(aMesh%panels) <= tracerRefine%tracerID ) refineTracer = .TRUE.
	if ( vortRefine%type == RELVORT_REFINE ) refineVort = .TRUE.
	if ( flowMapRefine%type == FLOWMAP_REFINE) refineFlowMap = .TRUE.
	!
	!	set existing mesh as data source for interpolation
	!
	call New(delTri,aMesh)
	call DelaunayTriangulation(delTri)
	call New(lagSource,delTri,vectorInterp)
	!if ( present(interpSmoothTol) ) call SetSigmaTol(lagSource,interpSmoothTol)
	call SetSourceLagrangianParameter(lagSource,delTri)

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,' remesh source data ready.')
	!
	!	Build a new mesh
	!
	call New(newMesh,aMesh%panelKind,aMesh%initNest,aMesh%AMR,aMesh%nTracer,aMesh%problemKind)
	newParticles => newMesh%particles
	newEdges => newMesh%edges
	newPanels => newMesh%panels

	!
	!	interpolate lagrangian parameter from old mesh to new mesh
	!
	do j=1,newParticles%N
		newParticles%x0(:,j) = InterpolateVector(newParticles%x(:,j),lagSource,delTri)
		! renormalize to spherical surface
		newParticles%x0(:,j) = newParticles%x0(:,j) / &
			sqrt(sum(newParticles%x0(:,j)*newParticles%x0(:,j)))*EARTH_RADIUS
	enddo
	do j=1,newPanels%N
		newPanels%x0(:,j) = InterpolateVector(newPanels%x(:,j),lagSource,delTri)
		! renormalize
		newPanels%x0(:,j) = newPanels%x0(:,j) / &
			sqrt(sum(newPanels%x0(:,j)*newPanels%x0(:,j)))*EARTH_RADIUS
	enddo
	!
	!	set tracer values on new mesh
	!
	if ( aMesh%nTracer > 0 ) then
		do j=1,newParticles%N
			do k=1,aMesh%nTracer
				newParticles%tracer(j,k) = GetTracer(reference,newParticles%x0(:,j),k)
			enddo
		enddo
		do j=1,newPanels%N
			if (.NOT. newPanels%hasChildren(j) ) then
				do k=1,aMesh%nTracer
					newPanels%tracer(j,k) = GetTracer(reference,newPanels%x0(:,j),k)
				enddo
			else
				newPanels%tracer(j,:) = 0.0_kreal
			endif
		enddo
	endif
	!
	!	set vorticity values on new mesh
	!
	if ( aMesh%problemKind == BVE_SOLVER ) then
		do j=1,newParticles%N
			newParticles%absVort(j) = GetAbsVort(reference,newParticles%x0(:,j))
			newParticles%relVort(j) = newParticles%absVort(j) - 2.0_kreal*OMEGA*newParticles%x(3,j)/EARTH_RADIUS
		enddo
		do j=1,newPanels%N
			newPanels%absVort(j) = GetAbsVort(reference,newPanels%x0(:,j))
			newPanels%relVort(j) = newPanels%absVort(j) - 2.0_kreal*OMEGA*newPanels%x(3,j)/EARTH_RADIUS
		enddo
	endif
	!
	!	AMR
	!
	if ( aMesh%AMR > 0 )  then
		allocate(refineFlag(newPanels%N_Max))
		refineFlag = .FALSE.
		startIndex = 1
		keepGoing = .FALSE.
		limit = 0
		!
		!	Apply refinement criteria
		!
		if ( refineTracer ) then
			limit = max(limit,tracerRefine%limit)
			call FlagPanelsForTracerMaxRefinement(refineFlag,newMesh,tracerRefine,startIndex)
			counter1 = count(refineFlag)
			call FlagPanelsForTracerVariationRefinement(refineFlag,newMesh,tracerRefine,startIndex)
			counter2 = count(refineFlag) - counter1
			write(formatString,'(A)') '(A,I8,A)'
			write(logString,formatString) 'LagRemesh AMR : tracerMax criterion triggered ', counter1, ' times.'
			call LogMessage(log,TRACE_LOGGING_LEVEL,logkey,logString)
			write(logString,formatString) 'LagRemesh AMR : tracerVar criterion triggered ', counter2, ' times.'
			call LogMessage(log,TRACE_LOGGING_LEVEL,logkey,logString)
		endif
		if ( refineVort) then
			limit = max(limit,vortRefine%limit)
			counter1 = count(refineFlag)
			call FlagPanelsForCirculationRefinement(refineFlag,newMesh,vortRefine,startIndex)
			counter1 = count(refineFlag) - counter1
			call FlagPanelsForRelVortVariationRefinement(refineFlag,newMesh,vortRefine,startIndex)
			counter2 = count(refineFlag) - counter1
			write(formatString,'(A)') '(A,I8,A)'
			write(logString,formatString) 'LagRemesh AMR : circMax criterion triggered ', counter1, ' times.'
			call LogMessage(log,TRACE_LOGGING_LEVEL,logkey,logString)
			write(logString,formatString) 'LagRemesh AMR : relVortVar criterion triggered ', counter2, ' times.'
			call LogMessage(log,TRACE_LOGGING_LEVEL,logkey,logString)
		endif
		if ( refineFlowMap ) then
			limit = max(limit,flowMapRefine%limit)
			counter1 = count(refineFlag)
			call FlagPanelsForFlowMapRefinement(refineFlag,newMesh,flowMapRefine,startIndex)
			counter1 = count(refineFlag) - counter1
			write(formatString,'(A)') '(A,I8,A)'
			write(logString,formatString) 'LagRemesh AMR : flowMap variation criterion triggered ', counter1, ' times.'
			call LogMessage(log,TRACE_LOGGING_LEVEL,logkey,logString)
		endif

		refineCount = count(refineFlag)
		spaceLeft = newPanels%N_Max - newPanels%N

		!
		!	exit if refinement is not needed, or insufficient memory
		!
		if ( refineCount == 0 ) then
			call LogMessage(log,TRACE_LOGGING_LEVEL,logkey,'LagRemesh : no refinement necessary.')
			keepGoing = .FALSE.
		elseif ( spaceLeft/4 < refineCount ) then
			call LogMessage(log,WARNING_LOGGING_LEVEL,logkey,'LagRemesh : insufficient memory for AMR.')
			keepGoing = .FALSE.
		else
			keepGoing = .TRUE.
		endif

		amrLoopCounter = 0

		do while (keepGoing)
			amrLoopCounter = amrLoopCounter + 1

			write(logString,'(A,I3,A,I8,A)') 'AMR loop ',amrLoopCounter,' : refining ',refineCount,' panels.'
			call LogMessage(log,TRACE_LOGGING_LEVEL,logkey,logString)
			!
			!	divide flagged panels
			!
			nOldPanels = newPanels%N
			nOldParticles = newParticles%N
			do j=startIndex,newPanels%N
				if ( refineFlag(j) ) then
					call DividePanel(newMesh,j)
					refineFlag(j) = .FALSE.
				endif
			enddo
			!
			!	ensure adjacent panels differ by no more than one level
			!
			call FlagPanelsAtRefinementBoundaries(refineFlag,newMesh)
			do j=1,newPanels%N
				if ( refineFlag(j) ) then
					call DividePanel(newMesh,j)
					refineFlag(j) = .FALSE.
				endif
			enddo

			!
			!	set problem data on mesh
			!
			do j=nOldParticles+1,newParticles%N
				newParticles%x0(:,j) = InterpolateVector(newParticles%x(:,j),lagSource,delTri)
				newParticles%x0(:,j) = newParticles%x0(:,j) / &
					sqrt(sum(newParticles%x0(:,j)*newParticles%x0(:,j)))*EARTH_RADIUS
			enddo
			do j=nOldPanels+1,newPanels%N
				newPanels%x0(:,j) = InterpolateVector(newPanels%x(:,j),lagSource,delTri)
				newPanels%x0(:,j) = newPanels%x0(:,j) / &
					sqrt(sum(newPanels%x0(:,j)*newPanels%x0(:,j)))*EARTH_RADIUS
			enddo
			if ( aMesh%nTracer > 0 ) then
				do j=nOldParticles+1,newParticles%N
					do k=1,aMesh%nTracer
						newParticles%tracer(j,k) = GetTracer(reference,newParticles%x0(:,j),k)
					enddo
				enddo
				do j=nOldPanels+1,newPanels%N
					do k=1,aMesh%nTracer
						newPanels%tracer(j,k) = GetTracer(reference,newPanels%x0(:,j),k)
					enddo
				enddo
			endif
			if ( aMesh%problemKind == BVE_SOLVER) then
				do j=nOldParticles+1,newParticles%N
					newParticles%absVort(j) = GetAbsVort(reference,newParticles%x0(:,j))
					newParticles%relVort(j) = newParticles%absVort(j) - 2.0_kreal*OMEGA*newParticles%x(3,j)/EARTH_RADIUS
				enddo
				do j=nOldPanels+1,newPanels%N
					newPanels%absVort(j) = GetAbsVort(reference,newPanels%x0(:,j))
					newPanels%relVort(j) = newPanels%absVort(j) - 2.0_kreal*OMEGA*newPanels%x(3,j)/EARTH_RADIUS
				enddo
			endif

			!
			!	prevent too much refinement
			!
			if ( amrLoopCounter >= limit ) then
				keepGoing = .FALSE.
				call LogMessage(log,WARNING_LOGGING_LEVEL,logkey,'LagRemesh WARNING : refinement limit reached.')
			endif

			!
			!	apply refinement criteria
			!
			startIndex = nOldPanels + 1
			nOldPanels = newPanels%N
			if ( refineTracer ) then
				limit = max(limit,tracerRefine%limit)
				call FlagPanelsForTracerMaxRefinement(refineFlag,newMesh,tracerRefine,startIndex)
				counter1 = count(refineFlag)
				call FlagPanelsForTracerVariationRefinement(refineFlag,newMesh,tracerRefine,startIndex)
				counter2 = count(refineFlag) - counter1
				write(formatString,'(A)') '(A,I8,A)'
				write(logString,formatString) 'LagRemesh AMR : tracerMax criterion triggered ', counter1, ' times.'
				call LogMessage(log,TRACE_LOGGING_LEVEL,logkey,logString)
				write(logString,formatString) 'LagRemesh AMR : tracerVar criterion triggered ', counter2, ' times.'
				call LogMessage(log,TRACE_LOGGING_LEVEL,logkey,logString)
			endif
			if ( refineVort) then
				limit = max(limit,vortRefine%limit)
				counter1 = count(refineFlag)
				call FlagPanelsForCirculationRefinement(refineFlag,newMesh,vortRefine,startIndex)
				counter1 = count(refineFlag) - counter1
				call FlagPanelsForRelVortVariationRefinement(refineFlag,newMesh,vortRefine,startIndex)
				counter2 = count(refineFlag) - counter1
				write(formatString,'(A)') '(A,I8,A)'
				write(logString,formatString) 'LagRemesh AMR : circMax criterion triggered ', counter1, ' times.'
				call LogMessage(log,TRACE_LOGGING_LEVEL,logkey,logString)
				write(logString,formatString) 'LagRemesh AMR : relVortVar criterion triggered ', counter2, ' times.'
				call LogMessage(log,TRACE_LOGGING_LEVEL,logkey,logString)
			endif
			if ( refineFlowMap ) then
				limit = max(limit,flowMapRefine%limit)
				counter1 = count(refineFlag)
				call FlagPanelsForFlowMapRefinement(refineFlag,newMesh,flowMapRefine,startIndex)
				counter1 = count(refineFlag) - counter1
				write(formatString,'(A)') '(A,I8,A)'
				write(logString,formatString) 'LagRemesh AMR : flowMap variation criterion triggered ', counter1, ' times.'
				call LogMessage(log,TRACE_LOGGING_LEVEL,logkey,logString)
			endif

			refineCount = count(refineFlag)
			spaceLeft = newPanels%N_Max - newPanels%N

			!
			!	exit if refinement is not needed, or insufficient memory
			!
			if ( refineCount == 0 ) then
				call LogMessage(log,TRACE_LOGGING_LEVEL,logkey,'LagRemesh : refinement comverged.')
				keepGoing = .FALSE.
			elseif ( spaceLeft/4 < refineCount ) then
				call LogMessage(log,WARNING_LOGGING_LEVEL,logkey,'LagRemesh : WARNING insufficient memory to continue AMR.')
				keepGoing = .FALSE.
			else
				keepGoing = .TRUE.
			endif
		enddo ! while keepgoing
		deallocate(refineFlag)
	endif ! AMR
	!
	!	replace old mesh with new mesh
	!
	call Copy(aMesh,newMesh)

	!
	!	clean up
	!
	call Delete(newMesh)
	call Delete(delTri)
	call Delete(lagSource)
end subroutine
!
!----------------
! Module methods : type-specific functions
!----------------
!
function GetAbsVort(self, alphaIn)
	real(kreal) :: GetAbsVort
	type(ReferenceSphere), intent(inout) :: self
	real(kreal), intent(in) :: alphaIn(3)
	!
	GetAbsVort = InterpolateScalar(alphaIn, self%absVortSource, self%delTri)
end function

function GetTracer(self, alphaIn, tracerID)
	real(kreal) :: GetTracer
	type(ReferenceSphere), intent(inout) :: self
	real(kreal), intent(in) :: alphaIn(3)
	integer(kint), intent(in) :: tracerID
	!
	GetTracer = InterpolateTracer(alphaIn, self%tracerSource, self%delTri, tracerID)
end function



subroutine InitLogger(aLog,rank)
	type(Logger), intent(inout) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	call New(aLog,logLevel)
	logInit = .TRUE.
end subroutine

end module
