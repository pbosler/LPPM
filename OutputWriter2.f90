module OutputWriterModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
!
! This module defines methods for formatting basic output to the screen or to a file unit, typically for logging purposes.
! Adapted from
!		McCormack, D. "Scientific software development with Fortran." www.mentalfaculty.com, 2009.
!
use NumberKindsModule

implicit none
private

public OutputWriter
public New, Delete
public StartSection, EndSection
public Write

!
!----------------
! Types and module constants
!----------------
!

type OutputWriter
	private
	integer(kint) :: fileUnit		! defines unit where output will be written
	integer(kint) :: indentLevel	! defines number of indentations to apply before writing first character
	character(len=MAX_STRING_LENGTH) :: fileName
end type

integer(KINT), parameter :: TAB_SPACE = 4

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

interface StartSection
	module procedure StartSectionWriter
end interface

interface EndSection
	module procedure EndSectionWriter
end interface

interface Write
	module procedure WriteString
	module procedure WriteInteger
	module procedure WriteReal
end interface


contains
!
!----------------
! Standard methods : Constructor / Destructor
!----------------
!

subroutine NewPrivate(self,fileUnit, fileName)
	type(OutputWriter), intent(out) :: self
	integer(kint), intent(in), optional :: fileUnit
	character(len=*), intent(in), optional :: fileName
	integer(kint) :: openstat

	self%indentLevel = 0_KINT

	if (present(fileUnit)) then
		if ( fileUnit /= STD_OUT) then
			if ( .NOT. present(fileName)) then
				self%fileUnit = STD_OUT
			else
				self%fileName = fileName
				self%fileUnit = fileUnit
			endif
		endif
	else
		self%fileUnit = STD_OUT
	endif

	if ( self%fileUnit /= STD_OUT .AND. self%fileUnit /= STD_ERR ) then
		open(unit=self%fileUnit,file=self%fileName,status='REPLACE',action='WRITE',iostat=openStat)
			if (openstat /= 0 ) then
				print *, "OutputWriter ERROR opening output file."
				print *, "Redirecting output to STD_OUT."
				self%fileUnit = STD_OUT
			endif
	endif
end subroutine


subroutine DeletePrivate(self)
	type(OutputWriter), intent(inout) :: self
	if ( self%fileUnit /= STD_OUT .AND. self%fileUnit /= STD_ERR ) then
		close(self%fileUnit)
	endif
end subroutine

!
!----------------
! Public functions
!----------------
!

subroutine WriteString(self,key,str)
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: key
	character(len=*), intent(in), optional :: str
	character(len=32) :: form
	if (present(str)) then
		form = FormatWithIndent(self,'(A,2X,A)')
		write(self%fileUnit,form) trim(key), trim(str)
	else
		form = FormatWithIndent(self,'(A)')
		write(self%fileUnit,form) trim(key)
	endif
end subroutine


subroutine WriteInteger(self,key,val)
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: key
	integer(kint), intent(in) :: val
	character(len=32) :: form
	form = FormatWithIndent(self,'(A,2X,I8)')
	write(self%fileUnit,form), trim(key),val
end subroutine

subroutine WriteReal(self,key,val)
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: key
	real(kreal), intent(in) :: val
	character(len=32) :: form
	if ( val >= 1.0d9 ) then
		form = FormatWithIndent(self,'(A,2X,E24.8E2)')
	else
		form = FormatWithIndent(self,'(A,2X,F24.15)')
	endif
	write(self%fileUnit,form) trim(key), val
end subroutine


subroutine StartSectionWriter(self,sectionName,description)
	type(OutputWriter), intent(inout) :: self
	character(len=*), intent(in) :: sectionName
	character(len=*), intent(in), optional :: description
	character(len=32) :: form
	form = FormatWithIndent(self,'(A)')
	write(self%fileUnit,form) sectionName
	self%indentLevel = self%indentLevel + 1
	if ( present(description) ) then
		form = FormatWithIndent(self,'(A)')
		write(self%fileUnit,form) description
	endif
end subroutine


subroutine EndSectionWriter(self)
	type(OutputWriter), intent(inout) :: self
	if ( self%indentLevel == 0 ) then
		print *, "EndSection WARNING : Indentation level already 0."
	else
		self%indentLevel = self%indentLevel - 1
	endif
end subroutine

!
!----------------
! Module methods : type-specific functions
!----------------
!

function FormatWithIndent(self,formatString)
	! Returns a string suitable for formatting write statements, with an indent level appropriate to
	! current state of OutputWriter::self
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: formatString
	character(len=len(formatString)+10) :: FormatWithIndent
	if ( self%indentLevel>0) then
		write(FormatWithIndent,'(A,I2,2A)') '(',TAB_SPACE*self%indentLevel,'X,',formatString(2:)
	else
		FormatWithIndent = formatString
	endif
end function


end module
