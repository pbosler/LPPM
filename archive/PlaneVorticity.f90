module PlaneVorticityModule
!------------------------------------------------------------------------------
! Lagrangian Particle / Panel Method
!------------------------------------------------------------------------------
!
!> @author
!> Peter Bosler, Department of Mathematics, University of Michigan
!
!> @defgroup PlaneVorticity Planar vorticity distributions
!> defines vorticity distributions for incompressible flow in 2d Cartesian geometry
!
!
! DESCRIPTION:
!> @file
!> defines vorticity distributions for incompressible flow in 2d Cartesian geometry
!
!------------------------------------------------------------------------------
use NumberKindsModule
use LoggerModule
use PlaneGeomModule
use ParticlesModule
use EdgesModule
use PanelsModule
use PlaneMeshModule

implicit none

private
public VorticitySetup
public New, Delete
public InitRankineVortex, SetRankineVortexOnMesh, RankineVortex, RANKINE_N_INT, RANKINE_N_REAL
public NullVorticity
public InitLambDipole, SetLambDipoleOnMesh, LambDipole, LAMB_DIPOLE_N_INT, LAMB_DIPOLE_N_REAL
public InitTwoDipoles, SetTwoDipolesOnMesh, TwoDipoles, TWO_DIPOLES_N_INT, TWO_DIPOLES_N_REAL

!
!----------------
! Types and module constants
!----------------
!
type VorticitySetup
	integer(kint), pointer :: integers(:) => null() 	! integers required by test case definitions
	real(kreal), pointer :: reals(:) => null()			! reals required by test case definitions
end type

integer(kint), parameter :: RANKINE_N_INT = 0, RANKINE_N_REAL = 4, &
							LAMB_DIPOLE_N_INT = 0, LAMB_DIPOLE_N_REAL = 4, &
							TWO_DIPOLES_N_INT = 0, TWO_DIPOLES_N_REAL = 8

real(kreal), parameter :: LAMB_K0 = 3.8317_kreal
!
!----------------
! Logging
!----------------
!
logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logKey = 'VorticitySetup'
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=128) :: logString
character(len=24) :: formatString
!
!----------------
! Interfaces
!----------------
!
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

subroutine NewPrivate(self, nInt, nReal)
	type(VorticitySetup), intent(out) :: self
	integer(kint), intent(in) :: nInt, nReal

	if ( .not. loginit) call InitLogger(log,procRank)

	if (nInt > 0 ) then
		allocate(self%integers(nInt))
		self%integers = 0
	else
		nullify(self%integers)
	endif

	if ( nReal > 0 ) then
		allocate(self%reals(nReal))
		self%reals = 0.0_kreal
	else
		nullify(self%reals)
	endif
end subroutine

subroutine DeletePrivate(self)
	type(VorticitySetup), intent(inout) :: self
	if ( associated(self%integers)) deallocate(self%integers)
	if ( associated(self%reals)) deallocate(self%reals)
	nullify(self%integers)
	nullify(self%reals)
end subroutine

!
!----------------
! Public functions
!----------------
!
subroutine NullVorticity(aMesh, self)
	type(PlaneMesh), intent(inout) :: aMesh
	type(VorticitySetup), intent(in) :: self
end subroutine


subroutine InitRankineVortex(self, xCent, yCent, radius, strength)
	type(VorticitySetup), intent(inout) :: self
	real(kreal), intent(in) :: xCent, yCent, radius, strength

	if ( size(self%reals) /= RANKINE_N_REAL ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'VorticitySetup ERROR : ',' real array size incorrect.')
		return
	endif

	self%reals(1) = xCent
	self%reals(2) = yCent
	self%reals(3) = radius
	self%reals(4) = strength
end subroutine

subroutine SetRankineVortexOnMesh(aMesh, self)
	type(PlaneMesh), intent(inout) :: aMesh
	type(VorticitySetup), intent(in) :: self
	!
	type(Particles), pointer :: aparticles
	type(Panels), pointer :: apanels
	integer(kint) :: j
	real(kreal) :: xyCent(2), radius, strength

	aparticles => aMesh%particles
	apanels => aMesh%panels

	xyCent = [ self%reals(1), self%reals(2)]
	radius = self%reals(3)
	strength = self%reals(4)

	do j=1,aParticles%N
		aParticles%relVort(j) = RankineVortex(aParticles%x0(:,j), xyCent, radius, strength)
	enddo
	do j=1,aPanels%N
		if (aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = 0.0_kreal
		else
			aPanels%relVort(j) = RankineVortex(aPanels%x0(:,j), xyCent, radius, strength)
		endif
	enddo

end subroutine


function RankineVortex(xy, xyCent, radius, strength)
	real(kreal) :: RankineVortex
	real(kreal), intent(in) :: xy(2), xyCent(2), radius, strength
	!
	real(kreal) :: r

	r = sqrt( sum( (xy - xyCent) * (xy - xyCent) ) )

	RankineVortex = 0.0_kreal

	if ( r <= radius) RankineVortex = strength
end function

subroutine InitLambDipole(self, lambR, U0, xCent, yCent)
	type(VorticitySetup), intent(inout) :: self
	real(kreal), intent(in) :: lambR, U0, xCent, ycent

	if ( size(self%reals) /= LAMB_DIPOLE_N_REAL ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'VorticitySetup ERROR : ',' real array size incorrect.')
		return
	endif

	self%reals(1) = lambR
	self%reals(2) = U0
	self%reals(3) = xCent
	self%reals(4) = yCent
end subroutine

subroutine SetLambDipoleOnMesh( aMesh, self)
	type(PlaneMesh), intent(inout) :: aMesh
	type(VorticitySetup), intent(in) :: self
	!
	type(Particles), pointer :: aparticles
	type(Panels), pointer :: apanels
	integer(kint) :: j

	aParticles => aMesh%particles
	aPanels => aMesh%panels

	do j = 1, aParticles%N
		aParticles%relVort(j) = LambDipole( aParticles%x0(:,j), self%reals(1), self%reals(2), self%reals(3), self%reals(4) )
	enddo
	do j = 1, aPanels%N
		if ( aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = 0.0_kreal
		else
			aPanels%relVort(j) = LambDipole(aPanels%x0(:,j), self%reals(1), self%reals(2), self%reals(3), self%reals(4) )
		endif
	enddo
end subroutine

function LambDipole( xy, lambR, U0, xcent, ycent)
	real(kreal) ::LambDipole
	real(kreal), intent(in) :: xy(2), xCent, yCent, lambR, U0
	!
	real(kreal) :: r, k, sintheta, denom
	real(kreal) :: BESSJ0, BESSJ1
	external BESSJ1, BESSJ0

	r = sqrt( (xy(1) - xcent) * (xy(1) - xcent) + (xy(2) - ycent) * (xy(2) - ycent) )

	if ( r > lambR .OR. r < ZERO_TOL) then
		LambDipole = 0.0_kreal
	else
		k = LAMB_K0 / lambR

		sintheta = xy(2) / r

		denom = BESSJ0( LAMB_K0 )

		LambDipole = -2.0_kreal * U0 * k * BESSJ1( k * r) * sintheta / denom
	endif
end function

subroutine InitTwoDipoles(self, xc1, yc1, rad1, u1, xc2, yc2, rad2, u2)
	type(VorticitySetup), intent(inout) :: self
	real(kreal), intent(in) :: xc1, yc1, rad1, u1, xc2, yc2, rad2, u2

	if ( size(self%reals) /= TWO_DIPOLES_N_REAL) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'VorticitySetup ERROR : ',' real array size incorrect.')
		return
	endif

	self%reals(1) = xc1
	self%reals(2) = yc1
	self%reals(3) = rad1
	self%reals(4) = u1
	self%reals(5) = xc2
	self%reals(6) = yc2
	self%reals(7) = rad2
	self%reals(8) = u2
end subroutine


subroutine SetTwoDipolesOnMesh( aMesh, self)
	type(PlaneMesh), intent(inout) :: aMesh
	type(VorticitySetup), intent(in) :: self
	!
	type(Particles), pointer :: aParticles
	type(Panels), pointer :: aPanels
	integer(kint) :: j

	aParticles => amesh%particles
	aPanels => aMesh%panels

	do j = 1, aParticles%N
		aparticles%relVort(j) = TwoDipoles(aParticles%x0(:,j), &
								self%reals(1), self%reals(2), self%reals(3), self%reals(4), &
								self%reals(5), self%reals(6), self%reals(7), self%reals(8) )
	enddo
	do j = 1, aPanels%N
		if ( aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = 0.0_kreal
		else
			aPanels%relVort(j) = TwoDipoles(aPanels%x0(:,j), &
								 self%reals(1), self%reals(2), self%reals(3), self%reals(4), &
								 self%reals(5), self%reals(6), self%reals(7), self%reals(8) )
		endif
	enddo
end subroutine

function TwoDipoles( xy, xc1, yc1, rad1, u1, xc2, yc2, rad2, u2)
	real(kreal) :: TwoDipoles
	real(kreal), intent(in) :: xy(2), xc1, yc1, rad1, u1, xc2, yc2, rad2, u2
	!
	TwoDipoles = LambDipole(xy, rad1, u1, xc1, yc1) + LambDipole(xy, rad2, u2, xc2, yc2)
end function


!
!----------------
! Module methods : module- or type-specific private functions
!----------------
!


subroutine InitLogger(aLog, rank)
	type(Logger), intent(inout) :: aLog
	integer(kint), intent(in) :: rank
	write(logKey,'(A,A,I0.2,A)') trim(logKey),'_',rank,' : '
	if ( rank == 0 ) then
		call New(aLog,logLevel)
	else
		call New(aLog,WARNING_LOGGING_LEVEL)
	endif
	logInit = .TRUE.
end subroutine

end module

FUNCTION BESSJ (N,X)

!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows.
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.

      PARAMETER (IACC = 40,BIGNO = 1.D10, BIGNI = 1.D-10)
      REAL *8 X,BESSJ,BESSJ0,BESSJ1,TOX,BJM,BJ,BJP,SUM
      IF (N.EQ.0) THEN
      BESSJ = BESSJ0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSJ = BESSJ1(X)
      RETURN
      ENDIF
      IF (X.EQ.0.) THEN
      BESSJ = 0.
      RETURN
      ENDIF
      print*,'here'
      TOX = 2./X
      IF (X.GT.FLOAT(N)) THEN
      BJM = BESSJ0(X)
      BJ  = BESSJ1(X)
      DO 11 J = 1,N-1
      BJP = J*TOX*BJ-BJM
      BJM = BJ
      BJ  = BJP
   11 CONTINUE
      BESSJ = BJ
      ELSE
      M = 2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
      BESSJ = 0.
      JSUM = 0
      SUM = 0.
      BJP = 0.
      BJ  = 1.
      DO 12 J = M,1,-1
      BJM = J*TOX*BJ-BJP
      BJP = BJ
      BJ  = BJM
      IF (ABS(BJ).GT.BIGNO) THEN
      BJ  = BJ*BIGNI
      BJP = BJP*BIGNI
      BESSJ = BESSJ*BIGNI
      SUM = SUM*BIGNI
      ENDIF
      IF (JSUM.NE.0) SUM = SUM+BJ
      JSUM = 1-JSUM
      IF (J.EQ.N) BESSJ = BJP
   12 CONTINUE
      SUM = 2.*SUM-BJ
      BESSJ = BESSJ/SUM
      ENDIF
      RETURN
      END

      FUNCTION BESSJ0 (X)
      REAL *8 X,BESSJ0,AX,FR,FS,Z,FP,FQ,XX

!     This subroutine calculates the First Kind Bessel Function of
!     order 0, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.

      REAL *8 Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
      -.2073370639D-5,.2093887211D-6 /
      DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3, &
      -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
      DATA R1,R2,R3,R4,R5,R6 /57568490574.D0,-13362590354.D0, &
      651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0 /
      DATA S1,S2,S3,S4,S5,S6 /57568490411.D0,1029532985.D0, &
      9494680.718D0,59272.64853D0,267.8532712D0,1.D0 /
      IF(X.EQ.0.D0) GO TO 1
      AX = ABS (X)
      IF (AX.LT.8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ0 = FR/FS
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-.785398164
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ0 = SQRT(.636619772/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
      ENDIF
      RETURN
    1 BESSJ0 = 1.D0
      RETURN
      END
! ---------------------------------------------------------------------------
      FUNCTION BESSJ1 (X)
      REAL *8 X,BESSJ1,AX,FR,FS,Z,FP,FQ,XX
!     This subroutine calculates the First Kind Bessel Function of
!     order 1, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
      REAL *8 Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4,  &
      .2457520174D-5,-.240337019D-6 /,P6 /.636619772D0 /
      DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,   &
      .8449199096D-5,-.88228987D-6,.105787412D-6 /
      DATA R1,R2,R3,R4,R5,R6 /72362614232.D0,-7895059235.D0, &
      242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0 /
      DATA S1,S2,S3,S4,S5,S6 /144725228442.D0,2300535178.D0, &
      18583304.74D0,99447.43394D0,376.9991397D0,1.D0 /

      AX = ABS(X)
      IF (AX.LT.8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ1 = X*(FR/FS)
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-2.35619491
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ1 = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*SIGN(S6,X)
      ENDIF
      RETURN
      END

