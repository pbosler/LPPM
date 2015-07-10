program lsqTest

implicit none

integer, parameter :: m = 20, n = 10
real(8) :: a(m,n), b(m), x(n)
integer :: i, j

character(len=56) :: formatstring

do j = 1, n
	do i = 1, m
		!a(i,j) = real( (i+j)/real(i*j,8) + i, 8)
		call random_number(a(i,j))
	enddo
enddo

do i = 1, m
	b(i) = 1.0d0
enddo

write(6,'(A)',advance='NO') "A = ["
do i = 1, m-1
	do j = 1, n-1
		write(6,'(F18.6,2X,A)',advance='NO') a(i,j), ", "
	enddo
	write(6,'(F18.6,2X,A)') a(i,n), "; ..."
enddo
do j = 1, n-1
	write(6,'(F18.6,2X,A)',advance='NO') a(m,j), ", "
enddo
write(6,'(F18.6,2X,A)') a(m,n), "];"

call wrapped_dgelsy( m, n, A, b, x )

write(6,'(A)', advance = 'NO') "coefficients = "
write(6,'(4(E18.6,2X))') x(1), x(2), x(3), x(4)
write(6,'(3(E18.6,2X))') x(5), x(6), x(7)
write(6,'(2(E18.6,2X))') x(8), x(9)
write(6,'(E18.6)') x(10)

contains

subroutine wrapped_dgelsy( m, n, A, b, x )
	integer, intent(in) :: m
	integer, intent(in) :: n
	real(8), intent(inout) :: A(m,n)
	real(8), intent(inout) :: b(m)
	real(8), intent(out) :: x(n)
	!
	integer :: nRhs, lda, ldb, jpvt(n), rank, lwork, info
	real(8) :: rcond
	real(8) :: work(4*n+1)
	
	nrhs = 1
	lda = m
	ldb = m
	lwork = 4*n+1
	jpvt = 0
	
	call DGELSY( m, n, nrhs, A, lda, b, ldb, jpvt, rcond, rank, work, lwork, info )
	
	if ( info == 0 ) then
		print *, "SUCCESS."
	elseif ( info < 0 ) then
		print *, "illegal value in argument ", -i
	else
		print *, "undefined error"
	endif
	x = b(1:n)
end subroutine

end program