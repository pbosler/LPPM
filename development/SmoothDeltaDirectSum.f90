module  SmoothDeltaDirectSumModule

use NumberKindsModule
use LoggerModule
use ParticlesModule
use EdgesModule, only : MaxEdgeLength
use PolyMesh2dModule
use FieldModule
use MPISetupModule
use SphereGeomModule, only : ChordDistance, SphereDistance, SphereProjection

implicit none

include 'mpif.h'

private

public SphereDelta, New, Delete
public InterpolateScalar

type SphereDelta
	real(kreal) :: eps
end type

interface New
	module procedure newPrivate
end interface

interface Delete
	module procedure deletePrivate
end interface

interface InterpolateScalar
	module procedure interpolateScalarPrivate
end interface

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'SphereDelta'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains
!----------------
!
! Public functions
!
!----------------
subroutine newPrivate(self, sphereMesh, smoothRadiusMultplier)
	type(SphereDelta), intent(out) :: self
	type(PolyMesh2D), intent(in) :: sphereMesh
	real(kreal), intent(in) :: smoothRadiusMultplier
	real(kreal) :: mult
	
	if ( .NOT. logInit ) call InitLogger(log, procRank)
	
	if ( present(smoothRadiusMultplier) ) then
		mult = smoothRadiusMultplier
	else
		mult = 2.0_kreal
	endif
	
	self%eps = mult * MaxEdgeLength( sphereMesh%edges, sphereMesh%particles)
end subroutine

subroutine deletePrivate(self)
	type(SphereDelta), intent(inout) :: self
end subroutine

function interpolateScalarPrivate(self, sphereMesh, scalarField, interpLoc )
	real(kreal) :: interpolateScalarPrivate
	type(SphereDelta), intent(in) :: self
	type(PolyMesh2D), intent(in) :: sphereMesh
	type(Field), intent(in) :: scalarField
	real(kreal), dimension(3), intent(in) :: interpLoc
	!
	integer(kint) :: j
	real(kreal) :: dotProd
	
	interpolateScalarPrivate = 0.0_kreal
	do j = 1, sphereMesh%particles%N
		if ( sphereMEsh%particles%isActive(j) then
			dotProd = sphereMesh%particles%x(j) * interpLoc(1) + sphereMesh%particles%y(j) * interpLoc(2) + &
					  sphereMesh%particles%z(j) * interpLoc(3)
			interpolateScalarPrivate = interpolateScalarPrivate + scalarField%scalar(j) * sphereMesh%particles%area(j) * &
				(-(1.0_kreal - dotProd)*(1.0_kreal - dotProd) + 2.0_kreal * self%eps*self%eps* dotProd) / &
				( 4.0_kreal * PI * ( 1.0_kreal - dotProd + self%eps*self%eps) * (1.0_kreal - dotProd + self%eps*self%eps))
		endif
	enddo
end function

subroutine interpolateScalarParallel(self, sphereMesh, scalarField, scalarInterp, xPts, yPts, zPts, interpMPI)
	type(SphereDelta), intent(in) :: self
	type(PolyMesh2d), intent(in) :: sphereMesh
	type(Field), intent(in) :: scalarField
	real(kreal), dimension(:), intent(out) :: scalarInterp
	real(kreal), dimension(:), intent(in) :: xPts, yPts, zPts
	type(MPISetup), intent(in) :: interpMPI
	!
	integer(kint) :: i, j, mpiErrCode
	real(kreal) :: dotProd, avg
	
	avg = ScalarAverage(scalarField, sphereMesh%particles)
	
	do i = interpMPI%indexStart(procRank), interpMPI%indexEnd(procRank)
		scalarInterp(i) = 0.0_kreal
		do j = 1, sphereMesh%particles%N
			if ( sphereMesh%particles%isActive(j) ) then
				dotProd = xPts(i) * sphereMesh%particles%x(j) + yPts(i) * sphereMesh%particles%y(j) + &
						  zPts(i) * sphereMesh%particles%z(j)
				scalarInterp(i) = scalarInterp(i) + scalarField%scalar(j) * sphereMesh%particles%area(j) * &
					(-(1.0_kreal - dotProd)*(1.0_kreal - dotProd) + 2.0_kreal * self%eps * self%eps * dotProd ) / &
					( 4.0_kreal * PI * (1.0_kreal - dotProd + self%eps*self%eps)*(1.0_kreal - dotProd + self%eps*self%eps))
			endif
		enddo	
		scalarInterp(i) = scalarInterp(i) - avg
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST(scalarInterp(interpMPI%indexStart(i):interpMPI%indexEnd(i)), interpMPI%messageLength(i), &
				MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	enddo	
end subroutine

subroutine interpolateScalarToLatLonParallel(self, sphereMesh, scalarField, scalarInterp, lons, lats, interpMPI )
	type(SphereDelta), intent(in) :: self
	type(PolyMesh2d), intent(in) :: sphereMesh
	type(Field), intent(in) :: scalarField
	real(kreal), dimension(:,:), intent(out) :: scalarInterp
	real(kreal), dimension(:), intent(in) :: lons
	real(kreal), dimension(:), intent(in) :: lats
	type(MPISetup), intent(in) :: interpMPI
	!
	integer(kint) :: i, j, k, mpiErrCode
	real(kreal) :: dotProd, avg
	
	avg = ScalarAverage(scalarField, sphereMesh%particles)
	
	do j = interpMPI%indexStart(procRank), interpMPI%indexEnd(procRank)
		do i = 1, size(lats)
			scalarInterp(i,j) = 0.0_kreal
			do k = 1, sphereMesh%particles%N
				if ( sphereMesh%particles%isActive(k)) then
					dotProd = cos(lons(j))*cos(lats(i)) * sphereMesh%particles%x(k) + &
							  sin(lons(j))*cos(lats(i)) * sphereMesh%particles%y(k) + &
							  sin(lats(i)) * sphereMesh%particles%z(k)
					scalarInterp(i,j) = scalarInterp(i,j) + scalarField%scalar(k) * sphereMesh%particles%area(k) * &
					(-(1.0_kreal - dotProd)*(1.0_kreal - dotProd) + 2.0_kreal * self%eps * self%eps * dotProd ) / &
					( 4.0_kreal * PI * (1.0_kreal - dotProd + self%eps*self%eps)*(1.0_kreal - dotProd + self%eps*self%eps))
				endif
			enddo
			scalarInterp(i,j) = scalarInterp(i,j) - avg
		enddo
	enddo
	
	do i = 0, numProcs - 1
		call MPI_BCAST( scalarInterp(:,interpMPI%indexStart(i):interpMPI%indexEnd(i)), size(lats)*interpMPI%messageLength(i), &
			MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiErrCode)
	enddo
end subroutine	

!----------------
!
! Private functions
!
!----------------

pure function deltaKernel( x, y, z, xt, yt, zt, eps )
	real(kreal) :: deltaKernel
	real(kreal), intent(in) :: x, y, z
	real(kreal), intent(in) :: xt, yt, zt
	real(kreal), intent(in) :: eps
	!
	real(kreal) :: dotProd
	
	dotProd = x * xt + y * yt + z * zt
	
	deltaKernel = ( -(1.0_kreal - dotProd) * (1.0_kreal - dotProd) + 2.0_kreal * eps * dotProd ) / &
				  ( 4.0_kreal * PI * ( 1.0_kreal - dotProd + eps*eps) * (1.0_kreal - dotProd + eps*eps))
end function


subroutine InitLogger(aLog,rank)
! Initialize a logger for this module and processor
	type(Logger), intent(out) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	if ( rank == 0 ) then
		call New(aLog,logLevel)
	else
		call New(aLog,ERROR_LOGGING_LEVEL)
	endif
	logInit = .TRUE.
end subroutine

end module
