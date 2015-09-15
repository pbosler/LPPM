module MPISetupModule

use NumberKindsModule

implicit none

private
public MPISetup, New, Delete, Copy

type MPISetup
	integer(kint), pointer :: indexStart(:) => null()
	integer(kint), pointer :: indexEnd(:) => null()
	integer(kint), pointer :: messageLength(:) => null()
	integer(kint) :: n = 0
	contains
		final :: deletePrivate
end type

interface New
	module procedure newPrivate
end interface

interface Copy
	module procedure copyPrivate
end interface

interface Delete
	module procedure deletePrivate
end interface

contains

subroutine newPrivate( self, nItems, nProcs )
	type(MPISetup), intent(out) :: self
	integer(kint), intent(in) :: nItems
	integer(kint), intent(in) :: nProcs
	!
	integer(kint) :: i, chunkSize 
	
	allocate(self%indexStart(0:nProcs-1))
	allocate(self%indexEnd(0:nProcs-1))
	allocate(self%messageLength(0:nProcs-1))
	
	chunkSize = nItems / nProcs
	do i = 0, nProcs - 1
		self%indexStart(i) = i * chunkSize + 1
		self%indexEnd(i) = (i+1) * chunkSize
	enddo
	self%indexEnd(nProcs-1) = nItems
	self%messageLength = self%indexEnd - self%indexStart + 1
	self%n = nItems
end subroutine

subroutine deletePrivate(self)
	type(MPISetup), intent(inout) :: self
	if ( associated(self%indexStart) ) then
		deallocate(self%indexStart)
	 	deallocate(self%indexEnd)
	 	deallocate(self%messageLength)
	endif
end subroutine

subroutine copyPrivate(self, other)
	type(MPISetup), intent(inout) :: self
	type(MPISetup), intent(in) :: other
	call deletePrivate(self)
	call newPrivate(self, other%n, numProcs)
end subroutine

end module