module NetCDFOutputModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method - Spherical Model
!------------------------------------------------------------------------------
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
!
!> @defgroup NetCDFOut NetCDF Output
!> Provides data structures and methods for outputting data from SphereMesh objects into NetCDF file format.
!
!
! DESCRIPTION:
!> @file
!> Provides data structures and methods for outputting data from SphereMesh objects into NetCDF file format.
!
!------------------------------------------------------------------------------
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
use netcdf

implicit none
private
public NetCDFSource
public New, Delete
public OutputNetCDF, UpdateFilename

type NetCDFSource
	character(len=MAX_STRING_LENGTH) :: filename
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

interface UpdateFilename
	module procedure UpdateFilenamePrivate
end interface

contains
!
!----------------
! Standard methods : Constructor / Destructor, Copy
!----------------
!
subroutine NewPrivate(self, aMesh, filename )
	type(NetCDFSource), intent(out) :: self
	type(SphereMesh), intent(in) :: aMesh	
	character(len=*), intent(in) :: filename
	self%filename = trim(filename)
end subroutine

subroutine DeletePrivate(self)
	type(NetCDFSource), intent(inout) :: self
end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine UpdateFilenamePrivate(self, filename)
	type(NetCDFSource), intent(inout) :: self
	character(len=*), intent(in) :: filename
	self%filename = trim(filename)
end subroutine

subroutine OutputNativeDataToNetCDF(self, aMesh)
	type(NetCDFSource), intent(in) :: self
	type(SphereMesh), intent(in) :: aMesh
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: ncID
	integer(kint) :: nParticles
	integer(kint) :: nActive
	integer(kint) :: passiveXID, passiveYID, passiveZID, activeXId, activeYId, activeZID
	character(len=*), parameter :: latName = "latitude"	
	character(len=*), parameter :: lonName = "longitude"	
	character(len=*), parameter :: xName = "x"	
	character(len=*), parameter :: yName = "y"	
	character(len=*), parameter :: zName = "z"
	character(len=*), parameter :: alpha_xName = "alpha_x"	
	character(len=*), parameter :: alpha_yName = "alpha_y"	
	character(len=*), parameter :: alpha_zName = "alpha_z"	
	character(len=*), parameter :: units = "units"
	character(len=*), parameter :: relVortUnits = "seconds^{-1}"
	character(len=*), parameter :: absVortUnits = "seconds^{-1}"
	character(len=*), parameter :: streamUnits = "metersSquared_per_seconds"
	character(len=*), parameter :: velocityUnits = "metersPerSecond"
	character(len=*), parameter :: tracerUnits = "tracerUnits"
	character(len=*), parameter :: xUnits = "earthRadii"
	character(len=*), parameter :: yUnits = "earthRadii"
	character(len=*), parameter :: zUnits = "earthRadii"
	character(len=*), parameter :: latUnits = "degrees_north"
	character(len=*), parameter :: lonUnits = "degrees_east"
	real(kreal), allocatable :: passiveLon(:), passiveLat(:), activeLon(:), activeLat(:)
	
	aParticles => aMesh%particles
	aPanels=>aMesh%panels
	nParticles = aParticles%n
	nActive = aPanels%N_Active
	
	call check( NF90_Create(self%filename, NF90_NETCDF4, ncID) )

	!
	! Define the coordinate variables
	!
!	call check( NF90_Def_Dim(ncID, "passive_x ", nParticles, passiveXID)
!	call check( NF90_Def_Dim(ncID, "passive_y ", nParticles, passiveYID)
!	call check( NF90_Def_Dim(ncID, "passive_z ", nParticles, passiveZID)
!	call check( NF90_Def_Dim(ncID, "active_x ", nActive, activeXID)
!	call check( NF90_Def_Dim(ncID, "active_y ", nActive, activeYID)
!	call check( NF90_Def_Dim(ncID, "active_z ", nActive, activeZID)
	call check( NF90_Def_Dim(ncID, "passive_Lon ", nParticles, passiveLonID)
	call check( NF90_Def_Dim(ncID, "passive_Lat ", nParticles, passiveLatID)
	call check( NF90_Def_Dim(ncID, "active_Lon ", nActive, activeLonID)
	call check( NF90_Def_Dim(ncID, "active_Lat ", nActive, activeLatID)
	
	!
	! Define the field variables
	!

	
		
end subroutine

subroutine check(status)
	 integer, intent ( in) :: status
	 if(status /= nf90_noerr) then
	   !print *, trim(nf90_strerror(status))
	   call LogMessage(log, ERROR_LOGGING_LEVEL,"NetCDF ERROR : ", status)
	 end if
end subroutine check


end module