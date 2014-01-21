module IntegerListModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!
! This module defines a simple linked-list data structure for dynamically sized integer data.
!
!
use NumberKindsModule

implicit none
private

public IntegerListNode
public New, Delete
public AddIntegerToList, ReturnIntegerArrayFromList
public ConcatenateLists
public PrintList
public AddUniqueIntegerToList

!
!----------------
! Types and module constants
!----------------
!

type IntegerListNode
	integer(kint) :: value
	integer(kint) :: listSize
	type(IntegerListNode), pointer :: Next
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

contains
!
!----------------
! Standard methods : Constructor / Destructor
!----------------
!

subroutine NewPrivate(listRoot,number)
! allocates a new listRoot.
	type(IntegerListNode), pointer, intent(out) :: listRoot
	integer(kint), intent(in), optional :: number
	allocate(listRoot)
	if ( present(number) ) then
		listRoot%value = number
		listRoot%listSize = 1_kint
	else
		listRoot%listSize = 0_kint
	endif
	nullify(listRoot%next)
end subroutine


recursive subroutine DeletePrivate(listRoot)
! Deletes an entire list
	type(IntegerListNode), pointer, intent(inout) :: listRoot
	type(IntegerListNode), pointer :: next
	next => listRoot%next
	deallocate(listRoot)
	if ( associated(next) ) then
		call DeletePrivate(next)
	endif
end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine AddIntegerToList(listRoot,number)
! increases size of list by 1, adds new integer to the end of current list.
	type(IntegerListNode), pointer, intent(inout) :: listRoot
	integer(kint), intent(in) :: number
	type(IntegerListNode), pointer :: current, next
	! check special case of empty list
	if ( listRoot%listSize == 0 ) then
		listRoot%value = number
		listRoot%listSize = 1_kint
		nullify(listRoot%next)
	else
		listRoot%listSize = listRoot%listSize + 1
		current => listRoot
		next => listRoot%next
		! Find end of list
		do while ( associated(next))
			current => current%next
			next => current%next
		enddo
		! allocate a new node
		allocate(next)
		! connect new node to end of list
		current%next => next
		next%value = number
		nullify(next%next)
	endif
end subroutine

subroutine AddUniqueIntegerToList(listRoot,number)
	type(IntegerListNode), pointer, intent(inout) :: listRoot
	integer(kint), intent(in) :: number
	type(IntegerListNode), pointer :: current, next
	! check special case of empty list
	if ( listRoot%listSize == 0 ) then
		listRoot%value = number
		listRoot%listSize = 1_kint
		nullify(listRoot%next)
	else
		! check root node value
		current => listRoot
		next => listRoot%next
		if ( current%value == number ) return
		! check current node
		do while ( associated(next) )
			current => current%next
			next => current%next
			if ( current%value == number ) return
		enddo
		! no equal values found; add number to end of list
		allocate(next)
		current%next => next
		next%value = number
		nullify(next%next)		
	endif
end subroutine

subroutine ConcatenateLists(list1,list2)
! Joins list 2 to the end of list1.
	type(IntegerListNode), pointer, intent(inout) :: list1
	type(IntegerListNode), pointer, intent(in) :: list2
	type(IntegerListNode), pointer :: current2
	! special case where list 2 is empty
	if ( list2%listSize == 0 ) then
		write(6,'(A)') 'ERROR : cannot concatenate empty lists.'
		return
	endif
	current2=>list2
	do while (associated(current2))
		call AddIntegerToList(list1,current2%value)
		current2=>current2%next
	enddo

end subroutine

subroutine ReturnIntegerArrayFromList(intArray,listRoot)
! Converts linked list data structure into an integer array.
! Integer array must be preallocated and correctly sized external to this program.
	integer(kint), intent(out) :: intArray(:)
	type(IntegerListNode), pointer, intent(in) :: listRoot
	type(IntegerListNode), pointer :: current
	integer(kint) :: j
	if ( size(intArray) /= listRoot%listSize ) then
		print *, "ReturnIntegerArrayFromList ERROR : size mismatch."
		return
	endif
	current => listRoot
	j = 1
	do while (associated(current))
		intArray(j) = current%value
		current=>current%next
		j = j+1
	enddo
end subroutine


subroutine PrintList(listRoot)
! Prints list data to standard out.
	type(IntegerListNode), pointer, intent(in) :: listRoot
	type(IntegerListNode), pointer :: current
	integer(kint) :: j
	current => listRoot
	j = 1
	write(6,'(A,I8)') 'ListSize = ',listRoot%listSize
	write(6,'(2A12)') 'List Entry','Value'
	do while ( associated(current) )
		write(6,'(2I12)') j, current%value
		current=>current%next
		j = j+1
	enddo
end subroutine



!
!----------------
! Module methods : type-specific functions
!----------------
!


end module
