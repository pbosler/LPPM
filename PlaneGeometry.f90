module PlaneGeomModule
!******************************************************************************
!	Peter A. Bosler
!	Department of Mathematics
!	University of Michigan
!	pbosler@umich.edu
!
!******************************************************************************
!
! This module defines functions for performing geometric calculations
! in the plane.
!
use NumberKindsModule

implicit none

public

contains
!----------------
! Basic geometry : length, area, centers of mass, etc.
!----------------

function Distance( xyA, xyB)
	real(kreal) :: Distance
	real(kreal), intent(in) :: xyA(2), xyB(2)
	Distance = sqrt( sum( (xyB - xyA)*(xyB - xyA) ) )
end function

function Midpoint( xyA, xyB )
	real(kreal) :: Midpoint(2)
	real(kreal), intent(in) :: xyA(2), xyB(2)
	Midpoint = 0.5_kreal*(xyA + xyB)
end function

function TriCentroid( xyA, xyB, xyC )
	real(kreal) :: TriCentroid(2)
	real(kreal), intent(in) :: xyA(2), xyB(2), xyC(2)
	TriCentroid = (xyA + xyB + xyC)/3.0_kreal
end function

function QuadCentroid( xyA, xyB, xyC, xyD )
	real(kreal) :: QuadCentroid(2)
	real(kreal), intent(in) :: xyA(2), xyB(2), xyC(2), xyD(2)
	QuadCentroid = 0.25_kreal*(xyA + xyB + xyC + xyD)
end function

function TriArea( xA, xB, xC )
	real(kreal) :: TriArea
	real(kreal), intent(in) :: xA(2), xB(2), xC(2)
	TriArea = 0.5_kreal*abs( -xB(1)*xA(2) + xC(1)*xA(2) + xA(1)*xB(2) - xC(1)*xB(2) - xA(1)*xC(2) + xB(1)*xC(2))
end function



end module
