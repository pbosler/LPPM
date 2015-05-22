program MovingVorticesAdvectionAMRVorticity

use NumberKindsModule
use OutputWriterModule
use LoggerModule
use SphereGeomModule
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
type(TracerSetup) :: testCaseTracer
integer(kint) :: tracerID
real(kreal) :: vortStartLon, vortStartLat

!
! vorticity placeholder
!
type(BVESetup) :: nullVort
real(kreal) :: maxCircTol, vortVarTol
!
! remeshing / refinement variables
!
type(RemeshSetup) :: remesh
integer(kint) :: remeshInterval, resetAlphaInterval, amrLimit, remeshCounter
real(kreal) :: tracerMassTol, tracerVarTol, lagVarTol
type(ReferenceSphere), pointer :: reference

!
! time stepping variables
!
type(AdvRK4Data) :: timekeeper
real(kreal) :: t, tfinal, dt
integer(kint) :: timesteps, timeJ

!
! output variables
!
type(VTKSource) :: vtkOut, vtkMeshOut
character(len = MAX_STRING_LENGTH) :: vtkRoot, vtkFile, vtkMeshFile, outputDir, jobPrefix, dataFile, summaryFile
character(len = 56) :: amrString
integer(kint) :: frameCounter, frameOut, readWriteStat
type(OutputWriter) :: writer

!
! test case variables
!
real(kreal), allocatable :: totalMasstestCaseTracer(:), sphereL1(:), sphereL2(:), sphereLinf(:), panelsLinf(:),&
						 	particlesLinf(:), phiMax(:), phiMin(:), tracerVar(:), panelTracer0(:), panelArea0(:)
real(kreal) :: deltaPhi, phimax0, phimin0
real(kreal) :: mass0, var0

!
! logging
!
type(Logger) :: exeLog
character(len=28) :: logkey
character(len=MAX_STRING_LENGTH) :: logstring

!
! mpi / computing environment / general variables
!
integer(kint) :: errCode
real(kreal) :: wallclock
integer(kint) :: j

!
! namelists and user input
!
character(len=MAX_STRING_LENGTH) :: namelistFile = 'MovingVortices2Direct.namelist'
namelist /meshDefine/ initNest, AMR, panelKind, amrLimit, maxCircTol, vortVarTol,&
	 tracerMassTol, tracerVarTol, lagVarTol
namelist /timestepping/ tfinal, dt, remeshInterval, resetAlphaInterval
namelist /fileIO/ outputDir, jobPrefix, frameOut


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	INITIALIZE COMPUTER, MESH, TEST CASE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

call MPI_INIT(errCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, errCode)
call MPI_COMM_RANK(MPI_COMM_WORLD, procRank, errCode)

call InitLogger(exeLog, procRank)

wallclock = MPI_WTIME()

nTracer = 3

!
! get user input
!
call ReadNamelistFile(procRank)

!
! define tracer
!
tracerID = 1
vortStartLon = 0.0_kreal
vortStartLat = 0.0_kreal
call New(testCaseTracer, 1, 2)
call InitMovingVortsTracer(testCaseTracer, vortStartLon, vortStartLat, tracerID)

!
! build initial mesh
!
call New(sphere, panelKind, initNest, AMR, nTracer, BVE_SOLVER)

sphereParticles => sphere%particles
spherePanels => sphere%panels



call SetTestCaseVorticityOnMesh(sphere, nullVort, 0.0_kreal)
call SetMovingVortsTracerOnMesh(sphere, testCaseTracer)

!
! initialize remeshing and refinement
!
call ConvertFromRelativeTolerances(sphere, maxCircTol, vortVarTol, tracerMassTol, tracerVarTol, tracerID, lagVarTol)
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'maxCircTol = ', maxCircTol )
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'vortVarTol  = ', vortVarTol )
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'tracerMassTol = ', tracerMassTol )
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, 'tracerVarTol  = ', tracerVarTol )
	call LogMessage(exeLog, TRACE_LOGGING_LEVEL, '   lagVarTol  = ', lagVarTol )
call New(remesh, maxCircTol, vortVarTol, lagVarTol, tracerID, tracerMassTol, tracerVarTol, amrLimit)
nullify(reference)
if ( AMR > 0 ) then
	!call InitialRefinement(sphere, remesh, SetMovingVortsTracerOnMesh,&
	!	testCaseTracer,SetTestCaseVorticityOnMesh,nullvort,0.0_kreal)
	if ( panelKind == QUAD_PANEL ) &
		write(amrstring,'(A,I1,A,I0.2,A)') 'quadAMR_', initNest, 'to', initNest+amrLimit, '_'
	if ( panelKind == TRI_PANEL ) &
		write(amrstring,'(A,I1,A,I0.2,A)') 'triAMR_', initNest, 'to', initNest+amrLimit, '_'
else
	if ( panelKind == QUAD_PANEL ) &
		write(amrstring,'(A,I1,A)') 'quadUnif_', initNest, '_'
	if ( panelKind == TRI_PANEL ) &
		write(amrstring,'(A,I1,A)') 'triUnif_', initNest, '_'
endif

if ( AMR <= 0 ) then 
	allocate(panelTracer0(spherePanels%N))
	allocate(panelArea0(spherePanels%N))
	do j = 1, spherePanels%N
		if ( spherePanels%hasChildren(j) ) then
			panelTracer0(j) = 0.0_kreal
			panelArea0(j) = 0.0_kreal
		else
			panelTracer0(j) = spherePanels%tracer(j,1)
			panelArea0(j) = spherePanels%area(j)
		endif
	enddo
endif

do j = 1, sphereParticles%N
	sphereParticles%tracer(j,3) = testCaseTracerExact(sphereParticles%x(:,j), 0.0_kreal,testCaseTracer)
enddo
do j = 1, spherePanels%N
	if ( spherePanels%hasChildren(j) ) then
		spherePanels%tracer(j,3) = 0.0_kreal
	else
		spherePanels%tracer(j,3) = testCaseTracerExact( spherePanels%x(:,j), 0.0_kreal,testCaseTracer)
	endif
enddo

!
! initialize output
!
if ( procrank == 0 ) then
	call LogStats( sphere, exeLog)

	write(vtkRoot,'(A,A,A,A,A)') trim(outputDir), '/vtkOut/',trim(jobPrefix),trim(amrString),'_'
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),0,'.vtk'
	
	write(vtkMeshFile,'(A,A,I0.4,A)') trim(vtkRoot), '_mesh_',0,'.vtk'
	
	write(summaryFile,'(A,A,A,A)') trim(outputDir), trim(jobPrefix), trim(amrString), '_summary.txt'
	write(datafile,'(A,A,A,A)') trim(outputDir), trim(jobPrefix), trim(amrstring), '_calculatedData.m'
	call New(vtkOut, sphere, vtkFile, 'moving vortices')
	call New(vtkMeshOut, sphere, vtkMeshFile, 'moving vortices')
	call VTKOutput(vtkOut, sphere)
	call VTKOutputMidpointRule(vtkMeshOut,sphere)
endif

!
! initialize time stepping
!
call New(timekeeper, sphere, numProcs)
timesteps = floor(tfinal / dt)
t = 0.0_kreal
remeshCounter = 0
frameCounter = 1
allocate(totalMasstestCaseTracer(0:timesteps))
totalMasstestCaseTracer = 0.0_kreal
mass0 = sum(panelTracer0*panelArea0)
allocate(sphereL2(0:timesteps))
sphereL2 = 0.0_kreal
allocate(sphereLinf(0:timesteps))
sphereLinf = 0.0_kreal
allocate(particlesLinf(0:timesteps))
particlesLinf = 0.0_kreal
allocate(panelsLinf(0:timesteps))
panelsLinf = 0.0_kreal
allocate(phiMax(0:timesteps))
phiMax = 0.0_kreal
allocate(phiMin(0:timesteps))
phiMin = 0.0_kreal
allocate(tracerVar(0:timesteps))
tracerVar = 0.0_kreal
var0 = TracerVariance(sphere, tracerID)
allocate(sphereL1(0:timesteps))
sphereL1 = 0.0_kreal


!phimax0 = max( maxval(sphereParticles%tracer(1:sphereParticles%N,1)), maxval(spherePanels%tracer(1:spherePanels%N,1)) )
!phimin0 = min( minval(sphereParticles%tracer(1:sphereParticles%N,1)), minval(spherePanels%tracer(1:spherePanels%N,1)) )
phimax0 = maxval(sphereParticles%tracer(1:sphereParticles%N,1))
phimin0 = maxval(sphereParticles%tracer(1:sphereParticles%N,1))
do j = 1, spherePanels%N
	if ( .NOT. spherePanels%hasChildren(j)) then
		if ( spherePanels%tracer(j,1) > phimax0 ) phimax0 = spherePanels%tracer(j,1)
		if ( spherePanels%tracer(j,1) < phimin0 ) phimin0 = spherePanels%tracer(j,1)
	endif
enddo
deltaPhi = phimax0 - phimin0


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	RUN THE PROBLEM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do timeJ = 0, timesteps - 1
	if ( mod( timeJ+1, remeshInterval) == 0 ) then
		!
		! remesh before timestep
		!
		remeshCounter = remeshCounter + 1
		call DirectRemesh(sphere, remesh)
		!
		! delete objects associated with old mesh
		!
		call Delete(timekeeper)
		if ( procrank == 0 ) then
			call Delete(vtkOUt)
			call Delete(vtkMeshOut)
		endif
		!
		! create new associated objects for new mesh
		!
		call New(timekeeper, sphere, numProcs)
		if ( procRank == 0 ) then	
			call New(vtkOut, sphere, vtkFile, 'moving vortices')
			call New(vtkMeshOut, sphere, vtkMeshFile, 'moving vortices')
		endif
		sphereParticles => sphere%particles
		spherePanels => sphere%panels
	endif ! remesh

	!
	! advance time
	!
	call AdvectionRK4Timestep(timekeeper, sphere, dt, t, procRank, numProcs, MovingVorticesVelocity)
	t = real( timeJ+1, kreal) * dt
	call SetTestCaseVorticityOnMesh(sphere, nullVort, t)
	

	do j = 1, sphereParticles%N
		sphereParticles%tracer(j,3) = testCaseTracerExact(sphereParticles%x(:,j), t, testCaseTracer)
	enddo
	do j = 1, spherePanels%N
		if ( spherePanels%hasChildren(j) ) then
			spherePanels%tracer(j,3) = 0.0_kreal
		else
			spherePanels%tracer(j,3) = testCaseTracerExact( spherePanels%x(:,j), t, testCaseTracer)
		endif
	enddo

	!
	! calculate error
	!
	do j = 1, sphereParticles%N
		sphereParticles%tracer(j,2) = (sphereParticles%tracer(j,1) - sphereParticles%tracer(j,3)) /&
			 maxval(abs(sphereParticles%tracer(1:sphereParticles%N,1)))
	enddo
	do j = 1, spherePanels%N
		if ( spherePanels%hasChildren(j) ) then
			spherePanels%tracer(j,2) = 0.0_kreal
		else
			spherePanels%tracer(j,2) = (spherePanels%tracer(j,1) - spherePanels%tracer(j,3))/&
				 maxval(abs(sphereParticles%tracer(1:sphereParticles%N,1)))
		endif
	enddo
	
	if ( AMR <= 0 ) then
		totalMasstestCaseTracer(timeJ + 1) = 0.0_kreal
		do j = 1, spherePanels%N
			if ( .NOT. spherePanels%hasChildren(j) ) then
				totalMasstestCaseTracer(timeJ + 1) = totalMasstestCaseTracer(timeJ + 1) + &
					(spherePanels%tracer(j,1) - panelTracer0(j)) * panelArea0(j)
			endif
		enddo	
		totalMasstestCaseTracer(timeJ+1) = totalMasstestCaseTracer(timeJ+1)/mass0
	else
		totalMasstestCaseTracer(timeJ+1) = TotalMass(sphere, tracerID) / mass0 - 1.0_kreal
	endif
	
	tracerVar(timeJ+1) = ( TracerVariance(sphere, tracerID) - var0 ) / var0

	particlesLinf(timeJ+1) = maxval(sphereParticles%tracer(1:sphereParticles%N,2)) /&
		 maxval(abs(sphereParticles%tracer(1:sphereParticles%N,1)))
	panelsLinf(timeJ+1) = maxval( spherePanels%tracer(1:spherePanels%N,2) ) / maxval( abs(spherePanels%tracer(1:spherePanels%N,1) ))

	sphereLinf(timeJ+1) = max( particlesLinf(timeJ+1), panelsLinf(timeJ+1) )
	sphereL2(timeJ+1) = sum( spherePanels%tracer(1:spherePanels%N,2) * &
		spherePanels%tracer(1:spherePanels%N,2) * spherePanels%area(1:spherePanels%N) )
	sphereL2(timeJ+1) = sphereL2(timeJ+1) / sum( spherePanels%tracer(1:spherePanels%N,1) *&
		 spherePanels%tracer(1:spherePanels%N,1) * spherePanels%area(1:spherePanels%N) )
	sphereL2(timeJ+1) = sqrt(sphereL2(timeJ+1))

	sphereL1(timeJ+1) = sum( abs(spherePanels%tracer(1:spherePanels%N,2)) *&
		 spherePanels%area(1:spherePanels%N) )
	sphereL1(timeJ+1) = sphereL1(timeJ+1) / &
		sum( abs(spherePanels%tracer(1:spherePanels%N,1)) * spherePanels%area(1:spherePanels%N) )

!	phimax(timeJ+1) = ( max( maxval(sphereParticles%tracer(1:sphereParticles%N,1)),&
!		 maxval( spherePanels%tracer(1:spherePanels%N,1)) ) - phimax0) / deltaPhi
!	phimin(timeJ+1) = ( min( minval(sphereParticles%tracer(1:sphereParticles%N,1)), &
!		 minval( spherePanels%tracer(1:spherePanels%N,1)) ) - phimin0)/ deltaPhi

	phimax(timeJ+1) = maxval(sphereParticles%tracer(1:sphereParticles%N,1))
	phimin(timeJ+1) = minval(sphereParticles%tracer(1:sphereParticles%N,1))
	do j = 1, spherePanels%N
		if ( .NOT. spherePanels%hasChildren(j) ) then
			if ( spherePanels%tracer(j,1) > phimax(timeJ+1)) phimax(timeJ+1) = spherePanels%tracer(j,1)
			if ( spherePanels%tracer(j,1) < phimin(timeJ+1)) phimin(timeJ+1) = spherePanels%tracer(j,1)
		endif
	enddo
	phimax(timeJ+1) = (phiMax(timeJ+1)-phimax0)/deltaPhi
	phimin(timeJ+1) = (phimin(timeJ+1)-phimin0)/deltaPhi

	!
	! output data
	!
	if ( procRank == 0 .AND. mod( timeJ+1, frameOut) == 0 ) then
		call LogMessage(exelog, TRACE_LOGGING_LEVEL, 'day = ', t/ONE_DAY)

		write(vtkFile, '(A,I0.4,A)') trim(vtkRoot), frameCounter, '.vtk'
		write(vtkMeshFile, '(A,A,I0.4,A)') trim(vtkRoot),'_mesh_',frameCounter,'.vtk'
		
		call UpdateFilename(vtkOut, vtkFile)
		call UpdateFilename(vtkMeshOut,vtkMeshFile)
		call VTKOutput(vtkOut, sphere)
		call VTKOutputMidpointRule(vtkMeshOut,sphere)

		frameCounter = frameCounter + 1
	endif
enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	OUTPUT FINAL DATA
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	if ( procRank == 0 ) then
		open( unit = WRITE_UNIT_1, file = datafile, status = 'REPLACE', action = 'WRITE', iostat = readwritestat)
		if ( readwritestat /= 0 ) then
			call LogMessage(exeLog, ERROR_LOGGING_LEVEL, 'data file ERROR : ', ' failed to open data file.')
		else
			write(WRITE_UNIT_1,'(A,F24.15,A)') 'passiveLinf = [', particlesLinf(0), ' ;'
			do j = 1, timesteps -1
				write(WRITE_UNIT_1,'(F24.15,A)') particlesLinf(j), ' ;'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') particlesLinf(timesteps), ' ];'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'activeLinf = [', panelsLinf(0), ' ;'
			do j = 1, timesteps -1
				write(WRITE_UNIT_1,'(F24.15,A)') panelsLinf(j), ' ;'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') panelsLinf(timesteps), ' ];'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'sphereLinf = [', sphereLinf(0), ' ;'
			do j = 1, timesteps -1
				write(WRITE_UNIT_1,'(F24.15,A)') sphereLinf(j), ' ;'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') sphereLinf(timesteps), ' ];'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'sphereL2 = [', sphereL2(0), ' ;'
			do j = 1, timesteps -1
				write(WRITE_UNIT_1,'(F24.15,A)') sphereL2(j), ' ;'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') sphereL2(timesteps), ' ];'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'sphereL1 = [', sphereL1(0), ' ;'
			do j = 1, timesteps -1
				write(WRITE_UNIT_1,'(F24.15,A)') sphereL1(j), ' ;'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') sphereL1(timesteps), ' ];'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'phi_max = [', phimax(0), ' ;'
			do j = 1, timesteps -1
				write(WRITE_UNIT_1,'(F24.15,A)') phimax(j), ' ;'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') phimax(timesteps), ' ];'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'phi_min = [', phimin(0), ' ;'
			do j = 1, timesteps -1
				write(WRITE_UNIT_1,'(F24.15,A)') phimin(j), ' ;'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') phimin(timesteps), ' ];'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'dt_day = ', dt / ONE_DAY, ' ;'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'tfinal_day = ', tfinal / ONE_DAY, ' ;'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'mass = [ ', totalMasstestCaseTracer(0), ' ; ...'
			do j = 1, timesteps-1
				write(WRITE_UNIT_1,'(F24.15,A)') totalMasstestCaseTracer(j), ' ; ...'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') totalMasstestCaseTracer(timesteps), ' ] ;'

			write(WRITE_UNIT_1,'(A,F24.15,A)') 'tracerVar = [ ', tracerVar(0), ' ; ...'
			do j = 1, timesteps-1
				write(WRITE_UNIT_1,'(F24.15,A)') tracerVar(j), ' ; ...'
			enddo
			write(WRITE_UNIT_1,'(F24.15,A)') tracerVar(timesteps), ' ] ;'
		endif
		close(WRITE_UNIT_1)

		write(logstring,'(A, F8.2,A)') 'elapsed time = ', (MPI_WTIME() - wallClock)/60.0, ' minutes.'
		call LogMessage(exelog,TRACE_LOGGING_LEVEL,'PROGRAM COMPLETE : ',trim(logstring))

	endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	FREE MEMORY, CLEAN UP, FINALIZE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (associated(reference)) then
	call Delete(reference)
	deallocate(reference)
endif
deallocate(totalMasstestCaseTracer)
deallocate(tracerVar)
deallocate(sphereL1)
deallocate(sphereL2)
deallocate(sphereLinf)
deallocate(particlesLinf)
deallocate(panelsLinf)
deallocate(phiMax)
deallocate(phiMin)
call Delete(timekeeper)
call Delete(remesh)
if ( procrank == 0 ) then
	call Delete(vtkOut)
	call Delete(vtkMeshOUt)
endif
call Delete(sphere)
call Delete(testCaseTracer)
call Delete(exeLog)
if ( allocated(panelTracer0)) deallocate(panelTracer0)
if ( allocated(panelArea0)) deallocate(panelArea0)
call MPI_FINALIZE(errCode)

contains

function testCaseTracerExact(xyz, t, mvTracer)
	real(kreal) :: testCaseTracerExact
	real(kreal), intent(in) :: xyz(3), t
	type(TracerSetup), intent(in) :: mvTracer
	!
	real(kreal) :: lat, lon, wr, rho, vortCenterLon, vortCenterLat, vortStartingLon, vortStartingLat
	real(kreal) :: lonPrime, latPrime
	real(kreal), parameter :: u0 = 2.0_kreal * PI * EARTH_RADIUS /  (12.0_kreal * ONE_DAY)

	lat = Latitude(xyz)
	lon = Longitude(xyz)

	vortStartingLon = mvTracer%reals(1)
	vortStartingLat = mvTracer%reals(2)
	!
	! find position of vortex center at time t
	!
	vortCenterLon = 1.5_kreal*PI + OMEGA * t / 12.0_kreal
	vortCenterLat = 0.0_kreal
	!
	! Find coordinates of xyz in a coordinate system whose north pole is at the vortex location
	!
	lonPrime = atan4( cos(lat)*sin( lon - vortCenterLon),  &
		cos(lat)*sin(vortCenterLat)*cos( lon - vortCenterLon) - cos(vortCenterLat)*sin(lat) )
	latPrime = asin( sin(lat)*sin(vortCenterLat) + cos(lat)*cos(vortCenterLat)*cos( lon - vortCenterLon ) )
	!
	! Determine angular tangential velocity induced by vortex about its center
	!
	rho = 3.0_kreal * cos( latPrime )
	wr = u0 * 1.5_kreal * sqrt(3.0_kreal) * tanh(rho) * rho /&
		 ( EARTH_RADIUS * cosh(rho) * cosh(rho) * (rho * rho + ZERO_TOL*ZERO_TOL))

	testCaseTracerExact = 1.0_kreal - tanh( 0.2_kreal * rho * sin(lonPrime - wr*t) )
end function

function MovingVortsVorticity(xyz, t)
	real(kreal) :: MovingVortsVorticity
	real(kreal), intent(in) :: xyz(3), t
	!
	real(kreal) :: lat, lon, rho, omg, rhoDenom, rho_lam, rho_theta, omg_rho, v_lam, ucostheta_theta
	real(kreal) :: cosDenom
	real(kreal), parameter :: u0 = 2.0_kreal * PI * EARTH_RADIUS /  (12.0_kreal * ONE_DAY)
		
	lat = Latitude(xyz)
	lon = Longitude(xyz)
	
	rho = 3.0_kreal*sqrt( 1.0_kreal - cos(lat)*cos(lat)*sin(lon - u0 * t / EARTH_RADIUS) * sin(lon - OMEGA*t/12.0_kreal) )
	rhoDenom = rho / (rho*rho + ZERO_TOL*ZERO_TOL)
	
	omg = u0 * 1.5_kreal * sqrt(3.0_kreal) * tanh( rho ) * rhoDenom / cosh(rho) /cosh(rho)
	omg_rho = u0 * 1.5_kreal * sqrt(3.0_kreal) * &
		( rho - tanh(rho)*(cosh(rho)*cosh(rho) + 2.0_kreal*rho*cosh(rho)*sinh(rho))) * &
		rhoDenom*rhoDenom / (cosh(rho)**4)
	
	rho_lam = -3.0_kreal*cos(lat)*cos(lat)*sin(lon-u0 * t / EARTH_RADIUS)*cos(lon-u0 * t / EARTH_RADIUS) * rhoDenom
	rho_theta = 3.0_kreal*cos(lat)*sin(lat)*sin(lon-u0 * t / EARTH_RADIUS)*sin(lon-u0 * t / EARTH_RADIUS) * rhoDenom
	
	v_lam = omg_rho * rho_lam * cos(lon-u0 * t / EARTH_RADIUS) - omg * sin(lon-OMEGA*t/12.0_kreal)
	ucostheta_theta = -omg * sin(lat)*sin(lat)*sin(lon-u0 * t / EARTH_RADIUS) + &
	  omg_rho*rho_theta*sin(lat)*sin(lon-u0 * t / EARTH_RADIUS) + omg*cos(lat)*sin(lon-u0 * t / EARTH_RADIUS)
	
	cosDenom = cos(lat)/( cos(lat)*cos(lat) + ZERO_TOL * ZERO_TOL)
	
	MovingVortsVorticity = v_lam * cosDenom / EARTH_RADIUS - ucostheta_theta * cosDenom / EARTH_RADIUS
end function

subroutine StoreVorticityInTracer(aMesh, t, tracerID)
	type(SphereMesh), intent(inout) :: aMesh
	real(kreal), intent(in) :: t
	integer(kint), intent(in) :: tracerID
	!
	integer(kint) :: j
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	
	aParticles => aMesh%particles
	aPanels => aMesh%panels
	
	do j = 1, aParticles%N
		aParticles%tracer(j,tracerID) = MovingVortsVorticity(aParticles%x(:,j), t)
	enddo
	do j = 1, aPanels%N
		if ( aPanels%hasCHildren(j) ) then
			aPanels%tracer(j,tracerID) = 0.0_kreal
		else
			aPanels%tracer(j,tracerID) = MovingVortsVorticity(aPanels%x(:,j), t)
		endif
	enddo
end subroutine

subroutine SetTestCaseVorticityOnMesh(aMesh, aVorticity, time)
	type(SphereMesh), intent(inout) :: aMesh
	type(BVESetup), intent(in) :: aVorticity
	real(kreal), intent(in) :: time
	!
	integer(kint) :: j
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	aParticles => aMesh%particles
	aPanels => aMesh%panels
	
	do j = 1, aParticles%N
		aParticles%absVort(j) = 0.0_kreal
		aParticles%relvort(j) = MovingVortsVorticity(aParticles%x(:,j), time)
	enddo
	do j = 1, aPanels%N
		if ( aPanels%hasCHildren(j) ) then
			aPanels%relvort(j) = 0.0_kreal
			aPanels%absVort(j) = 0.0_kreal
		else
			aPanels%relvort(j) = MovingVortsVorticity(aPanels%x(:,j), time)
			aPanels%absVort(j) = 0.0_kreal
		endif
	enddo
	
end subroutine

subroutine ConvertFromRelativeTolerances(aMesh, maxCircTol, vortVarTol, tracerMassTol, tracerVarTol, tracerID, lagVarTol)
	type(SphereMesh), intent(in) :: amesh
	real(kreal), intent(inout) :: maxCircTol, vortVarTol, tracerMassTol, tracerVarTol, lagVarTol
	integer(kint), intent(in) :: tracerID
	maxCircTol = maxCircTol * MaximumCirculation(aMesh)
	vortVarTol = vortVarTol * MaximumVorticityVariation(aMesh)
	tracerMassTol = tracerMassTol * MaximumTracerMass(aMesh, tracerID)
	tracerVarTol = tracerVarTol * MaximumTracerVariation(aMesh, tracerID)
	lagVarTol = lagVarTol * MaximumLagrangianParameterVariation(aMesh)
end subroutine

subroutine ReadNamelistFile(rank)
	integer(kint), intent(in) :: rank
	integer(kint), parameter :: BCAST_INT_SIZE = 6, BCAST_REAL_SIZE= 7
	integer(kint) :: broadcastIntegers(BCAST_INT_SIZE)
	real(kreal) :: broadcastReals(BCAST_REAL_SIZE)

	if ( rank == 0 ) then
		open(unit=READ_UNIT, file=namelistfile, status='OLD', action='READ', iostat=readWriteStat)
			if ( readWriteStat /= 0 ) stop 'cannot read namelist file.'
			read(READ_UNIT, nml=meshDefine)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=timestepping)
			rewind(READ_UNIT)
			read(READ_UNIT, nml=fileIO)
			rewind(READ_UNIT)
		close(READ_UNIT)

		broadcastIntegers(1) = panelKind
		broadcastIntegers(2) = initNest
		broadcastIntegers(3) = AMR
		broadcastIntegers(4) = amrLimit
		broadcastIntegers(5) = remeshInterval
		broadcastIntegers(6) = resetAlphaInterval

		broadcastReals(1) = tracerMassTol
		broadcastReals(2) = tracerVarTol
		broadcastReals(3) = dt
		broadcastReals(4) = tfinal
		broadcastReals(5) = lagVarTol
		broadcastReals(6) = maxCircTol
		broadcastReals(7) = vortVarTol
	endif

	call MPI_BCAST(broadcastIntegers, BCAST_INT_SIZE, MPI_INTEGER, 0, MPI_COMM_WORLD, errCode)
	panelKind = broadcastIntegers(1)
	initNest = broadcastIntegers(2)
	AMR = broadcastIntegers(3)
	amrLimit = broadcastIntegers(4)
	remeshInterval = broadcastIntegers(5)
	resetAlphaInterval = broadcastIntegers(6)

	call MPI_BCAST(broadcastReals, BCAST_REAL_SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errCode)
	tracerMassTol = broadcastReals(1)
	tracerVarTol = broadcastReals(2)
	dt = broadcastReals(3) * ONE_DAY		! convert time to seconds
	tfinal = broadcastReals(4) * ONE_DAY	! convert time to seconds
	lagVarTol = broadcastReals(5)
	maxCircTol = broadcastReals(6)
	vortVarTol = broadcastReals(7)
end subroutine

subroutine InitLogger(alog,rank)
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	if ( rank == 0 ) then
		call New(aLog,DEBUG_LOGGING_LEVEL)
	else
		call New(aLog,WARNING_LOGGING_LEVEL)
	endif
    write(logKey,'(A,I0.2,A)') 'EXE_LOG',rank,' : '
end subroutine

end program
