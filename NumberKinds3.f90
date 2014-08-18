module NumberKindsModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method - Spherical Model
!------------------------------------------------------------------------------
!
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
!> @defgroup NumberKinds Number kinds module
!> Defines numerical constants for use by all other LPPM modules and executables.
!> 
!
!
! DESCRIPTION:
!> @file
!> Defines numerical constants for use by all other LPPM modules and executables.
!
!------------------------------------------------------------------------------
	implicit none
	public

	! Type Constants
	integer, parameter :: KREAL = kind(0.d0)   !> @var compiler generated kind for double precision reals 
	integer, parameter :: KINT = kind(1)	   !> @var compiler generated kind for default integer
	integer, parameter :: KLOG = kind(.TRUE.)  !> @var compiler generated kind logical

	real(KREAL), parameter :: ZERO_TOL = 1.0d-14 !> @var zero tolerance for real numbers

	! Physical constants
	real(KREAL), parameter :: PI = 3.1415926535897932384626433832795027975_KREAL, & !> @var Pi
							  RAD_2_DEG = 180.0_kreal/PI, & !> @var convert radians to degrees
							  GRAV = 9.80616, & !> @var {m s^(-2)} acceleration due to gravity
							  ONE_DAY = 86140.0_kreal, & !> @var {s} Earth's sidereal day
							  EARTH_RADIUS = 6371220_kreal !> @var { m }

    real(KREAL), save ::	  OMEGA = 2.0_kreal*PI / ONE_DAY  !>@var {s^(-1)} rotation rate of sphere


	! I/O Constants
	integer(KINT), parameter :: STD_ERR = 0, &
								STD_IN  = 5, &
							    STD_OUT = 6, &
							    READ_UNIT = 11, &
							    WRITE_UNIT_1 = 12, &
							    WRITE_UNIT_2 = 13, &
							    WRITE_UNIT_3 = 14, &
							    MAX_STRING_LENGTH = 256
	! Mesh & application constants
	integer(KINT), parameter ::	QUAD_PANEL = 4, &
								TRI_PANEL = 3, &
								ADVECTION_SOLVER = 90, &
								BVE_SOLVER = 91, &
								SWE_SOLVER = 92, &
								PLANE_SOLVER = 99, &
								MAX_POLYGON_SIDES = 20, &
								FREE_BOUNDARIES = 101, &
								PERIODIC_BOUNDARIES = 102


	! MPI variables
	integer(KINT), save :: numProcs = 1, &
						   procRank = 0


end module NumberKindsModule

