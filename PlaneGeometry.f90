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

!----------------
! Misc. functions
!----------------
!
!    FUNCTION BESSJ (N,X)
!
!!     This subroutine calculates the first kind modified Bessel function
!!     of integer order N, for any REAL X. We use here the classical
!!     recursion formula, when X > N. For X < N, the Miller's algorithm
!!     is used to avoid overflows.
!!     REFERENCE:
!!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!!     MATHEMATICAL TABLES, VOL.5, 1962.
!
!      PARAMETER (IACC = 40,BIGNO = 1.D10, BIGNI = 1.D-10)
!      REAL *8 X,BESSJ,BESSJ0,BESSJ1,TOX,BJM,BJ,BJP,SUM
!      IF (N.EQ.0) THEN
!      BESSJ = BESSJ0(X)
!      RETURN
!      ENDIF
!      IF (N.EQ.1) THEN
!      BESSJ = BESSJ1(X)
!      RETURN
!      ENDIF
!      IF (X.EQ.0.) THEN
!      BESSJ = 0.
!      RETURN
!      ENDIF
!      print*,'here'
!      TOX = 2./X
!      IF (X.GT.FLOAT(N)) THEN
!      BJM = BESSJ0(X)
!      BJ  = BESSJ1(X)
!      DO 11 J = 1,N-1
!      BJP = J*TOX*BJ-BJM
!      BJM = BJ
!      BJ  = BJP
!   11 CONTINUE
!      BESSJ = BJ
!      ELSE
!      M = 2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
!      BESSJ = 0.
!      JSUM = 0
!      SUM = 0.
!      BJP = 0.
!      BJ  = 1.
!      DO 12 J = M,1,-1
!      BJM = J*TOX*BJ-BJP
!      BJP = BJ
!      BJ  = BJM
!      IF (ABS(BJ).GT.BIGNO) THEN
!      BJ  = BJ*BIGNI
!      BJP = BJP*BIGNI
!      BESSJ = BESSJ*BIGNI
!      SUM = SUM*BIGNI
!      ENDIF
!      IF (JSUM.NE.0) SUM = SUM+BJ
!      JSUM = 1-JSUM
!      IF (J.EQ.N) BESSJ = BJP
!   12 CONTINUE
!      SUM = 2.*SUM-BJ
!      BESSJ = BESSJ/SUM
!      ENDIF
!      RETURN
!      END
!
!      FUNCTION BESSJ0 (X)
!      REAL *8 X,BESSJ0,AX,FR,FS,Z,FP,FQ,XX
!
!!     This subroutine calculates the First Kind Bessel Function of
!!     order 0, for any real number X. The polynomial approximation by
!!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!!     REFERENCES:
!!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!!     VOL.5, 1962.
!
!      REAL *8 Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
!               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
!      DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
!      -.2073370639D-5,.2093887211D-6 /
!      DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3, &
!      -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
!      DATA R1,R2,R3,R4,R5,R6 /57568490574.D0,-13362590354.D0, &
!      651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0 /
!      DATA S1,S2,S3,S4,S5,S6 /57568490411.D0,1029532985.D0, &
!      9494680.718D0,59272.64853D0,267.8532712D0,1.D0 /
!      IF(X.EQ.0.D0) GO TO 1
!      AX = ABS (X)
!      IF (AX.LT.8.) THEN
!      Y = X*X
!      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
!      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
!      BESSJ0 = FR/FS
!      ELSE
!      Z = 8./AX
!      Y = Z*Z
!      XX = AX-.785398164
!      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
!      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
!      BESSJ0 = SQRT(.636619772/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
!      ENDIF
!      RETURN
!    1 BESSJ0 = 1.D0
!      RETURN
!      END
!! ---------------------------------------------------------------------------
!      FUNCTION BESSJ1 (X)
!      REAL *8 X,BESSJ1,AX,FR,FS,Z,FP,FQ,XX
!!     This subroutine calculates the First Kind Bessel Function of
!!     order 1, for any real number X. The polynomial approximation by
!!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!!     REFERENCES:
!!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!!     VOL.5, 1962.
!      REAL *8 Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6  &
!               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
!      DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4,  &
!      .2457520174D-5,-.240337019D-6 /,P6 /.636619772D0 /
!      DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,   &
!      .8449199096D-5,-.88228987D-6,.105787412D-6 /
!      DATA R1,R2,R3,R4,R5,R6 /72362614232.D0,-7895059235.D0, &
!      242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0 /
!      DATA S1,S2,S3,S4,S5,S6 /144725228442.D0,2300535178.D0, &
!      18583304.74D0,99447.43394D0,376.9991397D0,1.D0 /
!
!      AX = ABS(X)
!      IF (AX.LT.8.) THEN
!      Y = X*X
!      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
!      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
!      BESSJ1 = X*(FR/FS)
!      ELSE
!      Z = 8./AX
!      Y = Z*Z
!      XX = AX-2.35619491
!      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
!      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
!      BESSJ1 = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*SIGN(S6,X)
!      ENDIF
!      RETURN
!      END

end module
