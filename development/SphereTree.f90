module SphereTreeModule

use NumberKindsModule
use LoggerModule
use ParticlesModule
use FieldModule
use MPISetupModule
use SphereGeomModule, only : ChordDistance, SphereDistance, SphereProjection

implicit none

include 'mpif.h'

private

type kdPriorityQueue
	integer(kint) :: n = 0
	integer(kint), pointer :: elems(:) => null()
	real(kreal), pointer :: pri(:) => null()
	
	contains
		procedure, public :: newPriorityQueue
		final :: deletePriorityQueue
end type

interface New
	module procedure newPriorityQueue
end interface

interface Delete
	module procedure deletePriorityQueue
end interface

!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'SphereTree'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL

contains
!----------------
!
! Public functions
!
!----------------

subroutine newPriorityQueue( q, nmax )
	type(kdPriorityQueue), intent(out) :: q
	integer(kint), intent(in) :: nmax 
	allocate(q%elems(nmax))
	q%elems = 0
	allocate(q%pri(nmax))
	q%pri = 1.0d20
	q%n = 0
end subroutine

subroutine deletePriorityQueue( q )
	type(kdPriorityQueue), intent(inout) :: q
	if ( associated(q%elems) ) deallocate(q%elems)
	if ( associated(q%pri) ) deallocate(q%pri)
end subroutine


!----------------
!
! Private functions
!
!----------------

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