module LoggerModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method - Spherical Model
!------------------------------------------------------------------------------
!
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
!> @defgroup Logger Logger module
!> A Logger object for writing output to console or to files
!
!
! DESCRIPTION:
!> @file
!> A Logger object for writing output to console or to files
!
!------------------------------------------------------------------------------

use NumberKindsModule
use OutputWriterModule

implicit none
private

public Logger
public New, Delete
public DEBUG_LOGGING_LEVEL, TRACE_LOGGING_LEVEL, WARNING_LOGGING_LEVEL, ERROR_LOGGING_LEVEL
public LogMessage
public StartSection, EndSection

!
!----------------
! Types and module constants
!----------------
!

integer(kint), parameter :: DEBUG_LOGGING_LEVEL = 1, &
							TRACE_LOGGING_LEVEL = 2, &
							WARNING_LOGGING_LEVEL = 3, &
							ERROR_LOGGING_LEVEL = 4


type Logger
	integer(kint) :: level
	type(OutputWriter) :: writer
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

interface LogMessage
	module procedure LogString
	module procedure LogInteger
	module procedure LogReal
end interface

interface StartSection
	module procedure StartSectionLogger
end interface

interface EndSection
	module procedure EndSectionLogger
end interface

contains

!
!----------------
! Standard methods : Constructor / Destructor
!----------------
!

subroutine NewPrivate(self,level,fileUnit,fileName)
! Initializes a new Logger object, defines unit for chosen output.
	type(Logger), intent(out) :: self
	integer(kint), intent(in) :: level
	integer(kint), intent(in), optional :: fileUnit
	character(len=*), intent(in), optional :: fileName
	self%level = level

	if ( present(fileUnit) ) then
		if ( present(fileName) ) then
			call New(self%writer,fileUnit,fileName)
		else
			print *,"Logger ERROR : must supply filename for non-STD_OUT logging."
			print *,"Redirecting output to STD_OUT"
			call New(self%writer)
		endif
	else
		call New(self%writer)
	endif


end subroutine


subroutine DeletePrivate(self)
! Placeholder for future development.
	type(Logger), intent(inout) :: self
	call Delete(self%writer)
end subroutine


!
!----------------
! Public functions
!----------------
!

subroutine LogString(self,msgLevel,key,string)
	type(Logger), intent(in) :: self
	integer(kint), intent(in) :: msgLevel
	character(len=*), intent(in) :: key
	character(len=*), intent(in), optional :: string
	if ( msgLevel >= self%level ) then
		if ( present(string) ) then
			call Write(self%writer,key,string)
		else
			call Write(self%writer,key)
		endif
	endif
end subroutine

subroutine LogInteger(self,msgLevel,key,val)
	type(Logger), intent(in) :: self
	integer(kint), intent(in) :: msgLevel
	character(len=*), intent(in) :: key
	integer(kint), intent(in) :: val
	if ( msgLevel >= self%level) then
		call Write(self%writer,key,val)
	endif
end subroutine


subroutine LogReal(self,msgLevel,key,val)
	type(Logger), intent(in) :: self
	integer(kint), intent(in) :: msgLevel
	character(len=*), intent(in) :: key
	real(kreal), intent(in) :: val
	if (msgLevel>=self%level) then
		call Write(self%writer,key,val)
	endif
end subroutine

subroutine StartSectionLogger(self,sectionName,description)
	type(Logger), intent(inout) :: self
	character(len=*), intent(in), optional :: sectionName
	character(len=*), intent(in), optional :: description
	
	if ( present(sectionName) ) then
		if ( present(description) ) then
			call StartSection(self%writer,sectionName,description)
		else
			call StartSection(self%writer,sectionName)
		endif
	else
		call StartSection(self%writer)
	endif
end subroutine

subroutine EndSectionLogger(self)
	type(Logger), intent(inout) :: self
	call EndSection(self%writer)
end subroutine


!
!----------------
! Module methods : type-specific functions
!----------------
!



end module
