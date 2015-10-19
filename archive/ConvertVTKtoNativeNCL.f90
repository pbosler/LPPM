program ConvertNativeDataFromVTKtoNCL

use NumberKindsModule
use SphereGeomModule

implicit none

! file IO
character(len=MAX_STRING_LENGTH) :: namelistfile = 'convertNativeVTK.namelist'
character(len=MAX_STRING_LENGTH) :: inputfile
character(len=MAX_STRING_LENGTH) :: outputPrefix, outputFile
real(kreal) :: programStart, programEnd

namelist /filenames/ inputfile, outputPrefix

integer(kint) :: readStat, writeStat
! string processing
character(len=MAX_STRING_LENGTH) :: lineIn, scratch1, scratch2, formatString, foundName

integer(kint) :: lineNumber, scalarCount, polygonHeaderLine, pointDataLine
integer(kint) :: i, j, k, i1, i2, i3, ndim
logical(klog) :: keepGoing
! data conversion
real(kreal) :: xyz(3)
real(kreal), allocatable :: lon(:), lat(:), scalars(:), vectors(:,:)
integer(kint) :: nTotal, nTri

open(unit=READ_UNIT, file=namelistfile, status='OLD',action='READ',iostat = readstat)
	if ( readstat /= 0 ) stop "ERROR : cannot open namelist file"
	read(READ_UNIT,nml=filenames)
close(READ_UNIT)

! read vtk metadata
open(unit=READ_UNIT, file=inputfile, status='old', action='READ',iostat=readStat)
	if ( readStat /= 0 ) stop "ERROR opening VTK input file."
keepGoing = .TRUE.
lineNumber = 0
scalarCount = 0
polygonHeaderLine = -1
pointDataLine = -1 
print '("Getting VTK File MetaData ...")'
print '("VTK FILE INFO :")'
do while ( keepGoing )
	! read next line into lineIn
	read(READ_UNIT,'(A)',iostat = readStat) lineIn
	lineNumber = lineNumber + 1
	if ( readStat > 0 ) then
		print '("ERROR : problem encountered at line ",I8," .")' , lineNumber
		keepGoing = .FALSE.
		stop
	elseif ( readStat < 0 ) then
		print '("End of file found at line ",I8,".")', lineNumber
		keepGoing = .FALSE.
	else
		if ( lineNumber == 5) then 
			! get number of particles
			lineIn = adjustl(lineIn)
			lineIn = trim(lineIn)
			i1 = scan(lineIn,' ')
			i2 = scan(lineIn,' ',back=.TRUE.)
			scratch1 = lineIn(i1+1:i2)
			read(scratch1,'(I8)') nTotal
			print '("... found ",I8," total particles ...")', nTotal
			polygonHeaderLine = 5 + nTotal + 1
		elseif ( lineNumber == polygonHeaderLine ) then
			! count the number of triangles
			lineIn = adjustl(lineIn)
			lineIn = trim(lineIn)
			i1 = scan(lineIn,' ')
			i2 = scan(trim(lineIn),' ',back=.TRUE.)
			scratch1 = lineIn(i1:i2)
			scratch1 = trim(adjustl(scratch1))
			read(scratch1,'(I8)') nTri
			print '("... found ",I8," triangles ...")', nTri
			print '("... so there are = ",I8," active panels ...")', nTri/4
			pointDataLine = polygonHeaderLine + nTri + 1			
		elseif ( lineNumber > pointDataLine .AND. lineNumber > 5 ) then
			i1 = 0 
			lineIn = adjustl(lineIn)
			i1 = scan(lineIn,'S')
			if ( i1 == 1 ) scalarCount = scalarCount + 1			
		endif
	endif
enddo
print '("... found ", I8, " data fields.")', scalarCount

print '("Allocating memory ...")'


allocate(scalars(nTotal))
scalars = 0.0_kreal
allocate(vectors(3,nTotal))
vectors = 0.0_kreal
allocate(lon(nTotal))
lon = 0.0_kreal
allocate(lat(nTotal))
lat = 0.0_kreal


print '("Reading particle and field data ...")'

rewind(READ_UNIT)
keepGoing = .TRUE.
lineNumber = 0
j = 1
k = 1
nDim = 1
do while (keepGoing)
	read(READ_UNIT,'(A)',iostat = readStat) lineIN
	lineNumber = lineNumber + 1
	if ( readStat > 0 ) then
		print '("ERROR : problem encountered at line ",I8," .")' , lineNumber
		keepGoing = .FALSE.
		stop
	elseif ( readStat < 0 ) then
		print '("End of file found at line ",I8,".")', lineNumber
		keepGoing = .FALSE.
	else
		if ( lineNumber > 5 .AND. lineNumber < polygonHeaderLine ) then
			! read vertex coordinate data
			read ( lineIn, '(3D)') xyz
			lon(j) = Longitude(xyz)
			lat(j) = Latitude(xyz)
			j = j + 1
		endif
		if ( lineNumber == polygonHeaderLine ) then
			j = 1 ! reset counter
		! DEBUG	print '("Max lat = ",F5.3)', maxval(lat)
		endif
		if ( lineNumber > polygonHeaderLine .AND. lineNumber > pointDataLine ) then 
			lineIn = adjustl(lineIn)
			i1 = scan(lineIn,'S')
			if ( i1 == 1 ) then
				print '("... processing data field ",I2," of ", I2,".")', k, scalarCount
				k = k + 1
				j = 1
				i2 = scan(lineIn,' ')
				scratch1 = lineIn(i2+1:len_trim(lineIn))
				scratch1 = adjustl(scratch1)
				i3 = scan(scratch1,' ')
				foundName = scratch1(1:i3)				
				scratch2 = lineIn(len_trim(lineIN)-1:len_trim(lineIn))
				read(scratch2,'(I)') nDim
				write(6,'(A,A,A,I1)') "... data field ", trim(foundName), " has dimension ", ndim				
			endif
			i1 = 0
			i1 = scan(lineIn,'L')
			if ( i1 == 0 ) then
				! read scalar data	
				lineIn = trim(lineIn)
				if ( nDim == 1 ) then
					read(lineIn,'(D)') scalars(j)
					j = j + 1
				elseif ( nDim == 3) then
					read(lineIN,'(3F24.8)') vectors(:,j)
					j = j + 1
				endif
			endif
			if ( j == nTotal + 1 ) then
				write(outputFile,'(A,A,A,A)') trim(outputPrefix),'_',trim(foundName),'.txt'
				open(unit=WRITE_UNIT_1,file=outputfile,status='REPLACE',action='WRITE',iostat=writestat)
					if ( writeStat /= 0 ) stop "ERROR writing data file."
				write(WRITE_UNIT_1,'(3A24)') "lon", "lat", trim(foundname)
				if ( ndim == 1 ) then
					do i = 1, nTotal
						write(WRITE_UNIT_1,'(3F24.8)') lon(i), lat(i), scalars(i)
					enddo		
				elseif ( ndim == 3 ) then
					do i = 1, nTotal
						write(WRITE_UNIT_1,'(5F24.8)') lon(i), lat(i), vectors(:,i)
					enddo
				endif
				close(WRITE_UNIT_1)
			endif			
		endif
	endif
enddo
close(READ_UNIT)
call cpu_time(programEnd)
print *, " PROGRAM COMPLETE "
print '("Elapsed time : ", f7.3, " minutes.")', (programEnd - programStart)/60.0

deallocate(scalars)
deallocate(vectors)
deallocate(lat)
deallocate(lon)

end program