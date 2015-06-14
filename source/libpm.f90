!
!     file fish.f
!
!  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!  .                                                             .
!  .                  copyright (c) 2004 by UCAR                 .
!  .                                                             .
!  .       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH       .
!  .                                                             .
!  .                      all rights reserved                    .
!  .                                                             .
!  .                                                             .
!  .                      FISHPACK version 5.0                   .
!  .                                                             .
!  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                        F I S H P A C K                        *
!     *                                                               *
!     *                                                               *
!     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
!     *                                                               *
!     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
!     *                                                               *
!     *                  (Version 5.0 , JUNE 2004)                    *
!     *                                                               *
!     *                             BY                                *
!     *                                                               *
!     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
!     *                                                               *
!     *                             OF                                *
!     *                                                               *
!     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
!     *                                                               *
!     *                BOULDER, COLORADO  (80307)  U.S.A.             *
!     *                                                               *
!     *                   WHICH IS SPONSORED BY                       *
!     *                                                               *
!     *              THE NATIONAL SCIENCE FOUNDATION                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!

!     this module is used by all fishpack solvers to allocate
!     real and complex work space
      MODULE fish
	TYPE fishworkspace
	  DOUBLE PRECISION,POINTER,DIMENSION(:) :: rew
	  COMPLEX,POINTER,DIMENSION(:) :: cxw
	END TYPE fishworkspace
	CONTAINS
	SUBROUTINE allocatfish(irwk,icwk,wsave,ierror)
	IMPLICIT NONE
	TYPE (fishworkspace) :: wsave
!       irwk is the required real work space length
!       icwk is the required integer work space length
	INTEGER, INTENT(IN) :: irwk,icwk
!       ierror is set to 20 if the dynamic allocation is unsuccessful
!       (e.g., this would happen if m,n are too large for the computers memory
	INTEGER, INTENT(INOUT) :: ierror
	INTEGER :: istatus
!       first deallocate to avoid memory leakage
!#ifndef G95
 	if(associated(wsave%rew))DEALLOCATE(wsave%rew)
 	if(associated(wsave%cxw))DEALLOCATE(wsave%cxw)
!#endif
!       allocate irwk words of real work space
	if (irwk > 0) then
	     ALLOCATE(wsave%rew(irwk),STAT = istatus)
	end if
!       allocate icwk words of complex work space
	if (icwk > 0) then
	     ALLOCATE(wsave%cxw(icwk),STAT = istatus)
	end if
	ierror = 0
!       flag fatal error if allocation fails
!C       IF (istatus /= 0) THEN
	if (istatus .ne. 0 ) then
	  ierror = 20
	END IF	
	END SUBROUTINE allocatfish

	SUBROUTINE BLK_space(N,M,irwk,icwk)
!       this subroutine computes the real and complex work space
!       requirements (generous estimate) of blktri for N,M values
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: N,M
	INTEGER,INTENT(OUT) :: irwk,icwk
	INTEGER :: L,log2n
!       compute nearest integer greater than or equal to
!       log base 2 of n+1, i.e., log2n is smallest integer
!       such that 2**log2n >= n+1
	log2n = 1
	do
	   log2n = log2n+1
	   if (n+1 <= 2**log2n) EXIT
	end do
	L = 2**(log2n+1)
	irwk = (log2n-2)*L+5+MAX0(2*N,6*M)+log2n+2*n
	icwk = ((log2n-2)*L+5+log2n)/2+3*M+N	
	END SUBROUTINE BLK_space

	SUBROUTINE GEN_space(N,M,irwk)
!       this subroutine computes the real work space
!       requirement (generously) of genbun for the current N,M
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: N,M
	INTEGER,INTENT(OUT) :: irwk
	INTEGER :: log2n
!       compute nearest integer greater than or equal to
!       log base 2 of n+1, i.e., log2n is smallest integer
!       such that 2**log2n >= n+1
	log2n = 1
	do
	   log2n = log2n+1
	   if (n+1 <= 2**log2n) EXIT
	end do
	irwk = 4*N + (10 + log2n)*M	
	END SUBROUTINE GEN_space

	SUBROUTINE fishfin(wsave)
!       this subroutine releases allocated work space
!       fishfin should be called after a fishpack solver has finished
!       TYPE (fishworkspace) variable wsave.
	IMPLICIT NONE
	TYPE (fishworkspace) :: wsave
!!j     INTEGER :: istatus
!#ifndef G95
 	if(associated(wsave%rew))DEALLOCATE(wsave%rew)
 	if(associated(wsave%cxw))DEALLOCATE(wsave%cxw)
!#endif	
	END SUBROUTINE fishfin

      END MODULE fish
!
!
!   file pois3d.f
!
!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!.                                                             .
!.                  copyright (c) 2004 by UCAR                 .
!.                                                             .
!.       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH       .
!.                                                             .
!.                      all rights reserved                    .
!.                                                             .
!.                                                             .
!.                      FISHPACK version 5.0                   .
!.                                                             .
!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!   *                                                               *
!   *                        F I S H P A C K                        *
!   *                                                               *
!   *                                                               *
!   *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
!   *                                                               *
!   *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
!   *                                                               *
!   *                  (Version 5.0 , JUNE 2004)                    *
!   *                                                               *
!   *                             BY                                *
!   *                                                               *
!   *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
!   *                                                               *
!   *                             OF                                *
!   *                                                               *
!   *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
!   *                                                               *
!   *                BOULDER, COLORADO  (80307)  U.S.A.             *
!   *                                                               *
!   *                   WHICH IS SPONSORED BY                       *
!   *                                                               *
!   *              THE NATIONAL SCIENCE FOUNDATION                  *
!   *                                                               *
!   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!   SUBROUTINE POIS3D (LPEROD,L,C1,MPEROD,M,C2,NPEROD,N,A,B,C,LDIMF,
!  +                   MDIMF,F,IERROR)
!
!
! DIMENSION OF           A(N), B(N), C(N), F(LDIMF,MDIMF,N)
! ARGUMENTS
!
! LATEST REVISION        June 2004
!
! PURPOSE                SOLVES THE LINEAR SYSTEM OF EQUATIONS
!                      FOR UNKNOWN X VALUES, WHERE I=1,2,...,L,
!                      J=1,2,...,M, AND K=1,2,...,N
!
!                      C1*(X(I-1,J,K) -2.*X(I,J,K) +X(I+1,J,K)) +
!                      C2*(X(I,J-1,K) -2.*X(I,J,K) +X(I,J+1,K)) +
!                      A(K)*X(I,J,K-1) +B(K)*X(I,J,K)+ C(K)*X(I,J,K+1)
!                      = F(I,J,K)
!
!                      THE INDICES K-1 AND K+1 ARE EVALUATED MODULO N,
!                      I.E. X(I,J,0)=X(I,J,N) AND X(I,J,N+1)=X(I,J,1).
!                      THE UNKNOWNS
!                      X(0,J,K), X(L+1,J,K), X(I,0,K), AND X(I,M+1,K)
!                      ARE ASSUMED TO TAKE ON CERTAIN PRESCRIBED
!                      VALUES DESCRIBED BELOW.
!
! USAGE                  CALL POIS3D (LPEROD,L,C1,MPEROD,M,C2,NPEROD,
!                      N,A,B,C,LDIMF,MDIMF,F,IERROR)
!
! ARGUMENTS
!
! ON INPUT
!                      LPEROD
!                        INDICATES THE VALUES THAT X(0,J,K) AND
!                        X(L+1,J,K) ARE ASSUMED TO HAVE.
!                        = 0  X(0,J,K)=X(L,J,K), X(L+1,J,K)=X(1,J,K)
!                        = 1  X(0,J,K) = 0,      X(L+1,J,K) = 0
!                        = 2  X(0,J,K)=0,        X(L+1,J,K)=X(L-1,J,K)
!                        = 3  X(0,J,K)=X(2,J,K), X(L+1,J,K)=X(L-1,J,K)
!                        = 4  X(0,J,K)=X(2,J,K), X(L+1,J,K) = 0.
!
!                      L
!                        THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
!                        L MUST BE AT LEAST 3.
!
!                      C1
!                        DOUBLE PRECISION CONSTANT IN THE ABOVE LINEAR SYSTEM
!                        OF EQUATIONS TO BE SOLVED.
!
!                      MPEROD
!                        INDICATES THE VALUES THAT X(I,0,K) AND
!                        X(I,M+1,K) ARE ASSUMED TO HAVE.
!                        = 0  X(I,0,K)=X(I,M,K), X(I,M+1,K)=X(I,1,K)
!                        = 1  X(I,0,K)=0,        X(I,M+1,K)=0
!                        = 2  X(I,0,K)=0,        X(I,M+1,K)=X(I,M-1,K)
!                        = 3  X(I,0,K)=X(I,2,K)  X(I,M+1,K)=X(I,M-1,K)
!                        = 4  X(I,0,K)=X(I,2,K)  X(I,M+1,K)=0
!
!                      M
!                        THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
!                        M MUST BE AT LEAST 3.
!
!                      C2
!                        DOUBLE PRECISION CONSTANT IN THE ABOVE LINEAR SYSTEM
!                        OF EQUATIONS TO BE SOLVED.
!
!                      NPEROD
!                        = 0  IF A(1) AND C(N) ARE NOT ZERO.
!                        = 1  IF A(1) = C(N) = 0.
!
!                      N
!                        THE NUMBER OF UNKNOWNS IN THE K-DIRECTION.
!                        N MUST BE AT LEAST 3.
!
!                      A, B, C
!                        ONE-DIMENSIONAL ARRAYS OF LENGTH N THAT
!                        SPECIFY THE COEFFICIENTS IN THE LINEAR
!                        EQUATIONS GIVEN ABOVE.
!
!                        IF NPEROD = 0 THE ARRAY ELEMENTS MUST NOT
!                        DEPEND UPON INDEX K, BUT MUST BE CONSTANT.
!                        SPECIFICALLY,THE SUBROUTINE CHECKS THE
!                        FOLLOWING CONDITION
!                          A(K) = C(1)
!                          C(K) = C(1)
!                          B(K) = B(1)
!                        FOR K=1,2,...,N.
!
!                      LDIMF
!                        THE ROW (OR FIRST) DIMENSION OF THE THREE-
!                        DIMENSIONAL ARRAY F AS IT APPEARS IN THE
!                        PROGRAM CALLING POIS3D.  THIS PARAMETER IS
!                        USED TO SPECIFY THE VARIABLE DIMENSION
!                        OF F.  LDIMF MUST BE AT LEAST L.
!
!                      MDIMF
!                        THE COLUMN (OR SECOND) DIMENSION OF THE THREE
!                        DIMENSIONAL ARRAY F AS IT APPEARS IN THE
!                        PROGRAM CALLING POIS3D.  THIS PARAMETER IS
!                        USED TO SPECIFY THE VARIABLE DIMENSION
!                        OF F.  MDIMF MUST BE AT LEAST M.
!
!                      F
!                        A THREE-DIMENSIONAL ARRAY THAT SPECIFIES THE
!                        VALUES OF THE RIGHT SIDE OF THE LINEAR SYSTEM
!                        OF EQUATIONS GIVEN ABOVE.  F MUST BE
!                        DIMENSIONED AT LEAST L X M X N.
!
! ON OUTPUT
!
!                      F
!                        CONTAINS THE SOLUTION X.
!
!                      IERROR
!                        AN ERROR FLAG THAT INDICATES INVALID INPUT
!                        PARAMETERS.  EXCEPT FOR NUMBER ZERO, A
!                        SOLUTION IS NOT ATTEMPTED.
!                        = 0  NO ERROR
!                        = 1  IF LPEROD .LT. 0 OR .GT. 4
!                        = 2  IF L .LT. 3
!                        = 3  IF MPEROD .LT. 0 OR .GT. 4
!                        = 4  IF M .LT. 3
!                        = 5  IF NPEROD .LT. 0 OR .GT. 1
!                        = 6  IF N .LT. 3
!                        = 7  IF LDIMF .LT. L
!                        = 8  IF MDIMF .LT. M
!                        = 9  IF A(K) .NE. C(1) OR C(K) .NE. C(1)
!                             OR B(I) .NE.B(1) FOR SOME K=1,2,...,N.
!                        = 10 IF NPEROD = 1 AND A(1) .NE. 0
!                             OR C(N) .NE. 0
!                        = 20 If the dynamic allocation of real and
!                             complex work space required for solution
!                             fails (for example if N,M are too large
!                             for your computer)
!
!                        SINCE THIS IS THE ONLY MEANS OF INDICATING A
!                        POSSIBLY INCORRECT CALL TO POIS3D, THE USER
!                        SHOULD TEST IERROR AFTER THE CALL.
!
! SPECIAL CONDITIONS     NONE
!
! I/O                    NONE
!
! PRECISION              SINGLE
!
! REQUIRED files         fish.f,comf.f,fftpack.f
!
! LANGUAGE               FORTRAN 90
!
! HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
!                      1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
!                      LIBRARIES IN JANUARY, 1980.
!                      Revised in June 2004 by John Adams using
!                      Fortran 90 dynamically allocated work space.
!
! PORTABILITY            FORTRAN 90
!
! ALGORITHM              THIS SUBROUTINE SOLVES THREE-DIMENSIONAL BLOCK
!                      TRIDIAGONAL LINEAR SYSTEMS ARISING FROM FINITE
!                      DIFFERENCE APPROXIMATIONS TO THREE-DIMENSIONAL
!                      POISSON EQUATIONS USING THE FFT PACKAGE
!                      FFTPACK WRITTEN BY PAUL SWARZTRAUBER.
!
! TIMING                 FOR LARGE L, M AND N, THE OPERATION COUNT
!                      IS ROUGHLY PROPORTIONAL TO
!                        L*M*N*(LOG2(L)+LOG2(M)+5)
!                      BUT ALSO DEPENDS ON INPUT PARAMETERS LPEROD
!                      AND MPEROD.
!
! ACCURACY               TO MEASURE THE ACCURACY OF THE ALGORITHM A
!                      UNIFORM RANDOM NUMBER GENERATOR WAS USED TO
!                      CREATE A SOLUTION ARRAY X FOR THE SYSTEM GIVEN
!                      IN THE 'PURPOSE' SECTION WITH
!                        A(K) = C(K) = -0.5*B(K) = 1,  K=1,2,...,N
!                      AND, WHEN NPEROD = 1
!                        A(1) = C(N) = 0
!                        A(N) = C(1) = 2.
!
!                      THE SOLUTION X WAS SUBSTITUTED INTO THE GIVEN
!                      SYSTEM AND, USING DOUBLE PRECISION, A RIGHT
!                      SIDE Y WAS COMPUTED.  USING THIS ARRAY Y
!                      SUBROUTINE POIS3D WAS CALLED TO PRODUCE AN
!                      APPROXIMATE SOLUTION Z.  RELATIVE ERROR
!
!                      E = MAX(ABS(Z(I,J,K)-X(I,J,K)))/MAX(ABS(X(I,J,K
!
!                      WAS COMPUTED, WHERE THE TWO MAXIMA ARE TAKEN
!                      OVER I=1,2,...,L, J=1,2,...,M AND K=1,2,...,N.
!                      VALUES OF E ARE GIVEN IN THE TABLE BELOW FOR
!                      SOME TYPICAL VALUES OF L,M AND N.
!
!                      L(=M=N)   LPEROD    MPEROD       E
!                      ------    ------    ------     ------
!
!                        16        0         0        1.E-13
!                        15        1         1        4.E-13
!                        17        3         3        2.E-13
!                        32        0         0        2.E-13
!                        31        1         1        2.E-12
!                        33        3         3        7.E-13
!
! REFERENCES              NONE
! ********************************************************************
      SUBROUTINE POIS3D(LPEROD, L, C1, MPEROD, M, C2, NPEROD, N, A, B, C,&
                        LDIMF, MDIMF, F, IERROR)
      USE fish
      implicit none
      TYPE (fishworkspace) :: w
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: LPEROD
      INTEGER  :: L
      INTEGER  :: MPEROD
      INTEGER  :: M
      INTEGER  :: NPEROD
      INTEGER  :: N
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER  :: IERROR
      DOUBLE PRECISION  :: C1
      DOUBLE PRECISION  :: C2
      DOUBLE PRECISION  :: A(*)
      DOUBLE PRECISION  :: B(*)
      DOUBLE PRECISION  :: C(*)
      DOUBLE PRECISION  :: F(LDIMF,MDIMF,*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LP, MP, NP, K, IRWK, ICWK
!!j      DOUBLE PRECISION, DIMENSION(6) :: SAVE
!-----------------------------------------------
      LP = LPEROD + 1
      MP = MPEROD + 1
      NP = NPEROD + 1
!
!   CHECK FOR INVALID INPUT.
!
      IERROR = 0
      IF (LP<1 .OR. LP>5) IERROR = 1
      IF (L < 3) IERROR = 2
      IF (MP<1 .OR. MP>5) IERROR = 3
      IF (M < 3) IERROR = 4
      IF (NP<1 .OR. NP>2) IERROR = 5
      IF (N < 3) IERROR = 6
      IF (LDIMF < L) IERROR = 7
      IF (MDIMF < M) IERROR = 8
      IF (NP == 1) THEN
         DO K = 1, N
            IF (A(K) /= C(1)) GO TO 102
            IF (C(K) /= C(1)) GO TO 102
            IF (B(K) /= B(1)) GO TO 102
         END DO
         GO TO 104
  102    CONTINUE
         IERROR = 9
      ENDIF
      IF (NPEROD==1 .AND. (A(1)/=0. .OR. C(N)/=0.)) IERROR = 10
! 104 IF (IERROR .NE. 0) GO TO 122
  104 CONTINUE
      IF (IERROR /= 0) goto 999
!     allocate required work space length (generous estimate)
      IRWK=30+L+M+2*N+MAX0(L,M,N)+7*(INT((L+1)/2)+INT((M+1)/2))
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) goto 999
      call pois3dd(LPEROD,L,C1,MPEROD,M,C2,NPEROD,N,A,B,C,LDIMF,&
                   MDIMF,F,IERROR,w%rew)
!     release work space
      CALL FISHFIN (W)      
 999  END SUBROUTINE POIS3D



      SUBROUTINE POIS3DD(LPEROD, L, C1, MPEROD, M, C2, NPEROD, N, A, B,&
                         C, LDIMF, MDIMF, F, IERROR, W)
      implicit none
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: LPEROD
      INTEGER  :: L
      INTEGER , INTENT(IN) :: MPEROD
      INTEGER  :: M
      INTEGER , INTENT(IN) :: NPEROD
      INTEGER  :: N
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER  :: IERROR
      DOUBLE PRECISION  :: C1
      DOUBLE PRECISION  :: C2
      DOUBLE PRECISION  :: A(*)
      DOUBLE PRECISION  :: B(*)
      DOUBLE PRECISION  :: C(*)
      DOUBLE PRECISION  :: F(LDIMF,MDIMF,*)
      DOUBLE PRECISION  :: W(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LP, MP, NP, IWYRT, IWT, IWD, IWBB, IWX, IWY, NH, NHM1,&
                 NODD, I, J, K !!j, NHPK, NHMK
      DOUBLE PRECISION, DIMENSION(6) :: SAVE
!-----------------------------------------------
      LP = LPEROD + 1
      MP = MPEROD + 1
      NP = NPEROD + 1
      IWYRT = L + 1
      IWT = IWYRT + M
      IWD = IWT + MAX0(L,M,N) + 1
      IWBB = IWD + N
      IWX = IWBB + N
      IWY = IWX + 7*((L + 1)/2) + 15
      GO TO (105,114) NP
!
!   REORDER UNKNOWNS WHEN NPEROD = 0.
!
  105 CONTINUE
      NH = (N + 1)/2
      NHM1 = NH - 1
      NODD = 1
      IF (2*NH == N) NODD = 2
      DO I = 1, L
         DO J = 1, M
            DO K = 1, NHM1
               W(K) = F(I,J,NH-K) - F(I,J,K+NH)
               W(K+NH) = F(I,J,NH-K) + F(I,J,K+NH)
            END DO
            W(NH) = 2.*F(I,J,NH)
            GO TO (108,107) NODD
  107       CONTINUE
            W(N) = 2.*F(I,J,N)
  108       CONTINUE
            F(I,J,:N) = W(:N)
         END DO
      END DO
      SAVE(1) = C(NHM1)
      SAVE(2) = A(NH)
      SAVE(3) = C(NH)
      SAVE(4) = B(NHM1)
      SAVE(5) = B(N)
      SAVE(6) = A(N)
      C(NHM1) = 0.
      A(NH) = 0.
      C(NH) = 2.*C(NH)
      SELECT CASE (NODD)
      CASE DEFAULT
         B(NHM1) = B(NHM1) - A(NH-1)
         B(N) = B(N) + A(N)
      CASE (2)
         A(N) = C(NH)
      END SELECT
  114 CONTINUE
      CALL POS3D1 (LP, L, MP, M, N, A, B, C, LDIMF, MDIMF, F, W, W(IWYRT),&
                   W(IWT), W(IWD), W(IWX), W(IWY), C1, C2, W(IWBB))
      GO TO (115,122) NP
  115 CONTINUE
      DO I = 1, L
         DO J = 1, M
            W(NH-1:NH-NHM1:(-1))=0.5*(F(I,J,NH+1:NHM1+NH)+F(I,J,:NHM1))
            W(NH+1:NHM1+NH) = 0.5*(F(I,J,NH+1:NHM1+NH)-F(I,J,:NHM1))
            W(NH) = 0.5*F(I,J,NH)
            GO TO (118,117) NODD
  117       CONTINUE
            W(N) = 0.5*F(I,J,N)
  118       CONTINUE
            F(I,J,:N) = W(:N)
         END DO
      END DO
      C(NHM1) = SAVE(1)
      A(NH) = SAVE(2)
      C(NH) = SAVE(3)
      B(NHM1) = SAVE(4)
      B(N) = SAVE(5)
      A(N) = SAVE(6)
  122 CONTINUE      
      END SUBROUTINE POIS3DD


      SUBROUTINE POS3D1(LP, L, MP, M, N, A, B, C, LDIMF, MDIMF, F, XRT,&
                        YRT, T, D, WX, WY, C1, C2, BB)
      implicit none
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: LP
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: MP
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: LDIMF
      INTEGER , INTENT(IN) :: MDIMF
      DOUBLE PRECISION , INTENT(IN) :: C1
      DOUBLE PRECISION , INTENT(IN) :: C2
      DOUBLE PRECISION  :: A(*)
      DOUBLE PRECISION , INTENT(IN) :: B(*)
      DOUBLE PRECISION  :: C(*)
      DOUBLE PRECISION , INTENT(INOUT) :: F(LDIMF,MDIMF,1)
      DOUBLE PRECISION , INTENT(INOUT) :: XRT(*)
      DOUBLE PRECISION , INTENT(INOUT) :: YRT(*)
      DOUBLE PRECISION  :: T(*)
      DOUBLE PRECISION  :: D(*)
      DOUBLE PRECISION  :: WX(*)
      DOUBLE PRECISION  :: WY(*)
      DOUBLE PRECISION  :: BB(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LR, MR, NR, LRDEL, I, MRDEL, J, IFWRD, IS, K
      DOUBLE PRECISION :: PI, SCALX, DX, DI, SCALY, DY, DJ !!j,DUM
!-----------------------------------------------
      PI = 4.0*ATAN(1.0)
      LR = L
      MR = M
      NR = N
!
!   GENERATE TRANSFORM ROOTS
!
      LRDEL = ((LP - 1)*(LP - 3)*(LP - 5))/3
      SCALX = LR + LRDEL
      DX = PI/(2.*SCALX)
      GO TO (108,103,101,102,101) LP
  101 CONTINUE
      DI = 0.5
      SCALX = 2.*SCALX
      GO TO 104
  102 CONTINUE
      DI = 1.0
      GO TO 104
  103 CONTINUE
      DI = 0.0
  104 CONTINUE
      DO I = 1, LR
         XRT(I) = -4.*C1*SIN((FLOAT(I) - DI)*DX)**2
      END DO
      SCALX = 2.*SCALX
      GO TO (112,106,110,107,111) LP
  106 CONTINUE
      CALL SINTI (LR, WX)
      GO TO 112
  107 CONTINUE
      CALL COSTI (LR, WX)
      GO TO 112
  108 CONTINUE
      XRT(1) = 0.
      XRT(LR) = -4.*C1
      DO I = 3, LR, 2
         XRT(I-1) = -4.*C1*SIN(FLOAT(I - 1)*DX)**2
         XRT(I) = XRT(I-1)
      END DO
      CALL RFFTI (LR, WX)
      GO TO 112
  110 CONTINUE
      CALL SINQI (LR, WX)
      GO TO 112
  111 CONTINUE
      CALL COSQI (LR, WX)
  112 CONTINUE
      MRDEL = ((MP - 1)*(MP - 3)*(MP - 5))/3
      SCALY = MR + MRDEL
      DY = PI/(2.*SCALY)
      GO TO (120,115,113,114,113) MP
  113 CONTINUE
      DJ = 0.5
      SCALY = 2.*SCALY
      GO TO 116
  114 CONTINUE
      DJ = 1.0
      GO TO 116
  115 CONTINUE
      DJ = 0.0
  116 CONTINUE
      DO J = 1, MR
         YRT(J) = -4.*C2*SIN((FLOAT(J) - DJ)*DY)**2
      END DO
      SCALY = 2.*SCALY
      GO TO (124,118,122,119,123) MP
  118 CONTINUE
      CALL SINTI (MR, WY)
      GO TO 124
  119 CONTINUE
      CALL COSTI (MR, WY)
      GO TO 124
  120 CONTINUE
      YRT(1) = 0.
      YRT(MR) = -4.*C2
      DO J = 3, MR, 2
         YRT(J-1) = -4.*C2*SIN(FLOAT(J - 1)*DY)**2
         YRT(J) = YRT(J-1)
      END DO
      CALL RFFTI (MR, WY)
      GO TO 124
  122 CONTINUE
      CALL SINQI (MR, WY)
      GO TO 124
  123 CONTINUE
      CALL COSQI (MR, WY)
  124 CONTINUE
      IFWRD = 1
      IS = 1
  125 CONTINUE
!
!   TRANSFORM X
!
      DO J=1,MR
	 DO K=1,NR
	    DO I=1,LR
               T(I) = F(I,J,K)
	    END DO
            GO TO (127,130,131,134,135),LP
  127       GO TO (128,129),IFWRD
  128       CALL RFFTF (LR,T,WX)
            GO TO 138
  129       CALL RFFTB (LR,T,WX)
            GO TO 138
  130       CALL SINT (LR,T,WX)
            GO TO 138
  131       GO TO (132,133),IFWRD
  132       CALL SINQF (LR,T,WX)
            GO TO 138
  133       CALL SINQB (LR,T,WX)
            GO TO 138
  134       CALL COST (LR,T,WX)
            GO TO 138
  135       GO TO (136,137),IFWRD
  136       CALL COSQF (LR,T,WX)
            GO TO 138
  137       CALL COSQB (LR,T,WX)
  138       CONTINUE
	    DO I=1,LR
               F(I,J,K) = T(I)
	    END DO
	 END DO
      END DO
      GO TO (142,164) IFWRD
!
!   TRANSFORM Y
!
  142 CONTINUE
      DO I=1,LR
	 DO K=1,NR
	    DO J=1,MR
               T(J) = F(I,J,K)
	    END DO
            GO TO (144,147,148,151,152),MP
  144       GO TO (145,146),IFWRD
  145       CALL RFFTF (MR,T,WY)
            GO TO 155
  146       CALL RFFTB (MR,T,WY)
            GO TO 155
  147       CALL SINT (MR,T,WY)
            GO TO 155
  148       GO TO (149,150),IFWRD
  149       CALL SINQF (MR,T,WY)
            GO TO 155
  150       CALL SINQB (MR,T,WY)
            GO TO 155
  151       CALL COST (MR,T,WY)
            GO TO 155
  152       GO TO (153,154),IFWRD
  153       CALL COSQF (MR,T,WY)
            GO TO 155
  154       CALL COSQB (MR,T,WY)
  155       CONTINUE
	    DO J=1,MR
               F(I,J,K) = T(J)
	    END DO
	 END DO
      END DO
      GO TO (159,125) IFWRD
  159 CONTINUE
      DO I = 1, LR
         DO J = 1, MR
            BB(:NR) = B(:NR) + XRT(I) + YRT(J)
            T(:NR) = F(I,J,:NR)
            CALL TRID (NR, A, BB, C, T, D)
            F(I,J,:NR) = T(:NR)
         END DO
      END DO
      IFWRD = 2
      IS = -1
      GO TO 142
  164 CONTINUE
      F(:LR,:MR,:NR) = F(:LR,:MR,:NR)/(SCALX*SCALY)      
      END SUBROUTINE POS3D1


      SUBROUTINE TRID(MR, A, B, C, Y, D)
      implicit none
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: MR
      DOUBLE PRECISION , INTENT(IN) :: A(*)
      DOUBLE PRECISION , INTENT(IN) :: B(*)
      DOUBLE PRECISION , INTENT(IN) :: C(*)
      DOUBLE PRECISION , INTENT(INOUT) :: Y(*)
      DOUBLE PRECISION , INTENT(INOUT) :: D(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: M, MM1, I, IP
      DOUBLE PRECISION :: Z
!-----------------------------------------------
      M = MR
      MM1 = M - 1
      Z = 1./B(1)
      D(1) = C(1)*Z
      Y(1) = Y(1)*Z
      DO I = 2, MM1
         Z = 1./(B(I)-A(I)*D(I-1))
         D(I) = C(I)*Z
         Y(I) = (Y(I)-A(I)*Y(I-1))*Z
      END DO
      Z = B(M) - A(M)*D(MM1)
      IF (Z == 0.) THEN
         Y(M) = 0.
      ELSE
         Y(M) = (Y(M)-A(M)*Y(MM1))/Z
      ENDIF
      DO IP = 1, MM1
         I = M - IP
         Y(I) = Y(I) - D(I)*Y(I+1)
      END DO
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    Version 5.0, Fortran 90 Changes
!-----------------------------------------------------------------------
      END SUBROUTINE TRID

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


      SUBROUTINE COSTI(N, WSAVE)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      DOUBLE PRECISION  :: WSAVE(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NM1, NP1, NS2, K, KC
      DOUBLE PRECISION :: PI, DT, FK !!j , DUM
!-----------------------------------------------
!
      PI = 4.0*ATAN(1.0)
      IF (N <= 3) goto 998
      NM1 = N - 1
      NP1 = N + 1
      NS2 = N/2
      DT = PI/FLOAT(NM1)
      FK = 0.
      DO K = 2, NS2
         KC = NP1 - K
         FK = FK + 1.
         WSAVE(K) = 2.*SIN(FK*DT)
         WSAVE(KC) = 2.*COS(FK*DT)
      END DO
      CALL RFFTI (NM1, WSAVE(N+1))      
 998  END SUBROUTINE COSTI


      SUBROUTINE COST(N, X, WSAVE)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      DOUBLE PRECISION  :: X(*)
      DOUBLE PRECISION  :: WSAVE(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NM1, NP1, NS2, K, KC, MODN, I
      DOUBLE PRECISION :: X1H, X1P3, TX2, C1, T1, T2, XIM2, XI
!-----------------------------------------------
!
      NM1 = N - 1
      NP1 = N + 1
      NS2 = N/2
      IF (N - 2 >= 0) THEN
         IF (N - 2 <= 0) THEN
            X1H = X(1) + X(2)
            X(2) = X(1) - X(2)
            X(1) = X1H
            goto 997
         ENDIF
         IF (N <= 3) THEN
            X1P3 = X(1) + X(3)
            TX2 = X(2) + X(2)
            X(2) = X(1) - X(3)
            X(1) = X1P3 + TX2
            X(3) = X1P3 - TX2
            goto 997
         ENDIF
         C1 = X(1) - X(N)
         X(1) = X(1) + X(N)
         DO K = 2, NS2
            KC = NP1 - K
            T1 = X(K) + X(KC)
            T2 = X(K) - X(KC)
            C1 = C1 + WSAVE(KC)*T2
            T2 = WSAVE(K)*T2
            X(K) = T1 - T2
            X(KC) = T1 + T2
         END DO
         MODN = MOD(N,2)
         IF (MODN /= 0) X(NS2+1) = X(NS2+1) + X(NS2+1)
         CALL RFFTF (NM1, X, WSAVE(N+1))
         XIM2 = X(2)
         X(2) = C1
         DO I = 4, N, 2
            XI = X(I)
            X(I) = X(I-2) - X(I-1)
            X(I-1) = XIM2
            XIM2 = XI
         END DO
         IF (MODN /= 0) X(N) = XIM2
      ENDIF
 997  END SUBROUTINE COST


      SUBROUTINE SINTI(N, WSAVE)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      DOUBLE PRECISION  :: WSAVE(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS2, NP1, K
      DOUBLE PRECISION :: PI, DT !!j , DUM
!-----------------------------------------------
!
      PI = 4.0*ATAN(1.0)
      IF (N <= 1) goto 996
      NS2 = N/2
      NP1 = N + 1
      DT = PI/FLOAT(NP1)
      DO K = 1, NS2
         WSAVE(K) = 2.*SIN(K*DT)
      END DO
      CALL RFFTI (NP1, WSAVE(NS2+1))
 996  END SUBROUTINE SINTI


      SUBROUTINE SINT(N, X, WSAVE)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      DOUBLE PRECISION  :: X(*)
      DOUBLE PRECISION  :: WSAVE(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NP1, IW1, IW2, IW3
!-----------------------------------------------
!
      NP1 = N + 1
      IW1 = N/2 + 1
      IW2 = IW1 + NP1
      IW3 = IW2 + NP1
      CALL SINT1 (N, X, WSAVE, WSAVE(IW1), WSAVE(IW2), WSAVE(IW3))
      END SUBROUTINE SINT


      SUBROUTINE SINT1(N, WAR, WAS, XH, X, IFAC)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER  :: IFAC(*)
      DOUBLE PRECISION  :: WAR(*)
      DOUBLE PRECISION , INTENT(IN) :: WAS(*)
      DOUBLE PRECISION  :: XH(*)
      DOUBLE PRECISION  :: X(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, NP1, NS2, K, KC, MODN
      DOUBLE PRECISION :: SQRT3, XHOLD, T1, T2
!-----------------------------------------------
      DATA SQRT3/ 1.73205080756888/
      XH(:N) = WAR(:N)
      WAR(:N) = X(:N)
      IF (N - 2 <= 0) THEN
         IF (N - 2 /= 0) THEN
            XH(1) = XH(1) + XH(1)
            GO TO 106
         ENDIF
         XHOLD = SQRT3*(XH(1)+XH(2))
         XH(2) = SQRT3*(XH(1)-XH(2))
         XH(1) = XHOLD
         GO TO 106
      ENDIF
      NP1 = N + 1
      NS2 = N/2
      X(1) = 0.
      DO K = 1, NS2
         KC = NP1 - K
         T1 = XH(K) - XH(KC)
         T2 = WAS(K)*(XH(K)+XH(KC))
         X(K+1) = T1 + T2
         X(KC+1) = T2 - T1
      END DO
      MODN = MOD(N,2)
      IF (MODN /= 0) X(NS2+2) = 4.*XH(NS2+1)
      CALL RFFTF1 (NP1, X, XH, WAR, IFAC)
      XH(1) = 0.5*X(1)
      DO I = 3, N, 2
         XH(I-1) = -X(I)
         XH(I) = XH(I-2) + X(I-1)
      END DO
      IF (MODN == 0) XH(N) = -X(N+1)
  106 CONTINUE
      X(:N) = WAR(:N)
      WAR(:N) = XH(:N)
      END SUBROUTINE SINT1


      SUBROUTINE COSQI(N, WSAVE)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      DOUBLE PRECISION  :: WSAVE(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K
      DOUBLE PRECISION :: PIH, DT, FK !!j , DUM
!-----------------------------------------------
!
      PIH = 2.0*ATAN(1.0)
      DT = PIH/FLOAT(N)
      FK = 0.
      DO K = 1, N
         FK = FK + 1.
         WSAVE(K) = COS(FK*DT)
      END DO
      CALL RFFTI (N, WSAVE(N+1))
      END SUBROUTINE COSQI


      SUBROUTINE COSQF(N, X, WSAVE)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      DOUBLE PRECISION  :: X(*)
      DOUBLE PRECISION  :: WSAVE(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      DOUBLE PRECISION :: SQRT2, TSQX
!-----------------------------------------------
      DATA SQRT2/ 1.4142135623731/
!
      IF (N - 2 >= 0) THEN
         IF (N - 2 > 0) GO TO 103
         TSQX = SQRT2*X(2)
         X(2) = X(1) - TSQX
         X(1) = X(1) + TSQX
      ENDIF
      goto 995
  103 CONTINUE
      CALL COSQF1 (N, X, WSAVE, WSAVE(N+1))
 995  END SUBROUTINE COSQF


      SUBROUTINE COSQF1(N, X, W, XH)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      DOUBLE PRECISION  :: X(*)
      DOUBLE PRECISION , INTENT(IN) :: W(*)
      DOUBLE PRECISION  :: XH(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS2, NP2, K, KC, MODN, I
      DOUBLE PRECISION :: XIM1
!-----------------------------------------------
      NS2 = (N + 1)/2
      NP2 = N + 2
      DO K = 2, NS2
         KC = NP2 - K
         XH(K) = X(K) + X(KC)
         XH(KC) = X(K) - X(KC)
      END DO
      MODN = MOD(N,2)
      IF (MODN == 0) XH(NS2+1) = X(NS2+1) + X(NS2+1)
      DO K = 2, NS2
         KC = NP2 - K
         X(K) = W(K-1)*XH(KC) + W(KC-1)*XH(K)
         X(KC) = W(K-1)*XH(K) - W(KC-1)*XH(KC)
      END DO
      IF (MODN == 0) X(NS2+1) = W(NS2)*XH(NS2+1)
      CALL RFFTF (N, X, XH)
      DO I = 3, N, 2
         XIM1 = X(I-1) - X(I)
         X(I) = X(I-1) + X(I)
         X(I-1) = XIM1
      END DO
      END SUBROUTINE COSQF1


      SUBROUTINE COSQB(N, X, WSAVE)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      DOUBLE PRECISION  :: X(*)
      DOUBLE PRECISION  :: WSAVE(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      DOUBLE PRECISION :: TSQRT2, X1
!-----------------------------------------------
      DATA TSQRT2/ 2.82842712474619/
!
      IF (N - 2 <= 0) THEN
         IF (N - 2 /= 0) THEN
            X(1) = 4.*X(1)
            goto 994
         ENDIF
         X1 = 4.*(X(1)+X(2))
         X(2) = TSQRT2*(X(1)-X(2))
         X(1) = X1
         goto 994
      ENDIF
      CALL COSQB1 (N, X, WSAVE, WSAVE(N+1))
 994  END SUBROUTINE COSQB


      SUBROUTINE COSQB1(N, X, W, XH)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      DOUBLE PRECISION  :: X(*)
      DOUBLE PRECISION , INTENT(IN) :: W(*)
      DOUBLE PRECISION  :: XH(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS2, NP2, I, MODN, K, KC
      DOUBLE PRECISION :: XIM1
!-----------------------------------------------
      NS2 = (N + 1)/2
      NP2 = N + 2
      DO I = 3, N, 2
         XIM1 = X(I-1) + X(I)
         X(I) = X(I) - X(I-1)
         X(I-1) = XIM1
      END DO
      X(1) = X(1) + X(1)
      MODN = MOD(N,2)
      IF (MODN == 0) X(N) = X(N) + X(N)
      CALL RFFTB (N, X, XH)
      DO K = 2, NS2
         KC = NP2 - K
         XH(K) = W(K-1)*X(KC) + W(KC-1)*X(K)
         XH(KC) = W(K-1)*X(K) - W(KC-1)*X(KC)
      END DO
      IF (MODN == 0) X(NS2+1) = W(NS2)*(X(NS2+1)+X(NS2+1))
      DO K = 2, NS2
         KC = NP2 - K
         X(K) = XH(K) + XH(KC)
         X(KC) = XH(K) - XH(KC)
      END DO
      X(1) = X(1) + X(1)
      END SUBROUTINE COSQB1


      SUBROUTINE SINQI(N, WSAVE)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      DOUBLE PRECISION  :: WSAVE(*)
!-----------------------------------------------
!
      CALL COSQI (N, WSAVE)
      END SUBROUTINE SINQI


      SUBROUTINE SINQF(N, X, WSAVE)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      DOUBLE PRECISION  :: X(*)
      DOUBLE PRECISION  :: WSAVE(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS2, K, KC
      DOUBLE PRECISION :: XHOLD
!-----------------------------------------------
!
      IF (N == 1) goto 993
      NS2 = N/2
      DO K = 1, NS2
         KC = N - K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
      END DO
      CALL COSQF (N, X, WSAVE)
      X(2:N:2) = -X(2:N:2)
 993  END SUBROUTINE SINQF


      SUBROUTINE SINQB(N, X, WSAVE)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      DOUBLE PRECISION  :: X(*)
      DOUBLE PRECISION  :: WSAVE(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS2, K, KC
      DOUBLE PRECISION :: XHOLD
!-----------------------------------------------
!
      IF (N <= 1) THEN
         X(1) = 4.*X(1)
         goto 992
      ENDIF
      NS2 = N/2
      X(2:N:2) = -X(2:N:2)
      CALL COSQB (N, X, WSAVE)
      DO K = 1, NS2
         KC = N - K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
      END DO
 992  END SUBROUTINE SINQB


      SUBROUTINE RFFTI(N, WSAVE)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      DOUBLE PRECISION  :: WSAVE(*)
!-----------------------------------------------
!
      IF (N == 1) goto 991
      CALL RFFTI1 (N, WSAVE(N+1), WSAVE(2*N+1))
 991  END SUBROUTINE RFFTI


      SUBROUTINE RFFTI1(N, WA, IFAC)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(INOUT) :: IFAC(*)
      DOUBLE PRECISION , INTENT(OUT) :: WA(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(4) :: NTRYH
      INTEGER :: NL, NF, J, NTRY, NQ, NR, I, IS, NFM1, L1, K1, IP,&
                 LD, L2, IDO, IPM, II !!j , IB
      DOUBLE PRECISION :: TPI, ARGH, ARGLD, FI, ARG !!j , DUM
!-----------------------------------------------
      DATA NTRYH(1), NTRYH(2), NTRYH(3), NTRYH(4)/ 4, 2, 3, 5/
      NL = N
      NF = 0
      J = 0
  101 CONTINUE
      J = J + 1
      IF (J - 4 <= 0) THEN
         NTRY = NTRYH(J)
      ELSE
         NTRY = NTRY + 2
      ENDIF
  104 CONTINUE
      NQ = NL/NTRY
      NR = NL - NTRY*NQ
      IF (NR /= 0) GO TO 101
      NF = NF + 1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY == 2) THEN
         IF (NF /= 1) THEN
            IFAC(NF+2:4:(-1)) = IFAC(NF+1:3:(-1))
            IFAC(3) = 2
         ENDIF
      ENDIF
      IF (NL /= 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 8.0*ATAN(1.0)
      ARGH = TPI/FLOAT(N)
      IS = 0
      NFM1 = NF - 1
      L1 = 1
      IF (NFM1 == 0) goto 990
      DO K1 = 1, NFM1
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP - 1
         DO J = 1, IPM
            LD = LD + L1
            I = IS
            ARGLD = FLOAT(LD)*ARGH
            FI = 0.
            DO II = 3, IDO, 2
               I = I + 2
               FI = FI + 1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
            END DO
            IS = IS + IDO
         END DO
         L1 = L2
      END DO
 990  END SUBROUTINE RFFTI1


      SUBROUTINE RFFTB(N, R, WSAVE)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      DOUBLE PRECISION  :: R(*)
      DOUBLE PRECISION  :: WSAVE(*)
!-----------------------------------------------
!
      IF (N == 1) goto 989
      CALL RFFTB1 (N, R, WSAVE, WSAVE(N+1), WSAVE(2*N+1))
 989  END SUBROUTINE RFFTB


      SUBROUTINE RFFTB1(N, C, CH, WA, IFAC)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: IFAC(*)
      DOUBLE PRECISION  :: C(*)
      DOUBLE PRECISION  :: CH(*)
      DOUBLE PRECISION  :: WA(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NF, NA, L1, IW, K1, IP, L2, IDO, IDL1, IX2, IX3, IX4 !!j, I
!-----------------------------------------------
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO K1 = 1, NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP == 4) THEN
            IX2 = IW + IDO
            IX3 = IX2 + IDO
            IF (NA == 0) THEN
               CALL RADB4 (IDO, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
            ELSE
               CALL RADB4 (IDO, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
            ENDIF
            NA = 1 - NA
         ELSE
            IF (IP == 2) THEN
               IF (NA == 0) THEN
                  CALL RADB2 (IDO, L1, C, CH, WA(IW))
               ELSE
                  CALL RADB2 (IDO, L1, CH, C, WA(IW))
               ENDIF
               NA = 1 - NA
            ELSE
               IF (IP == 3) THEN
                  IX2 = IW + IDO
                  IF (NA == 0) THEN
                     CALL RADB3 (IDO, L1, C, CH, WA(IW), WA(IX2))
                  ELSE
                     CALL RADB3 (IDO, L1, CH, C, WA(IW), WA(IX2))
                  ENDIF
                  NA = 1 - NA
               ELSE
                  IF (IP == 5) THEN
                     IX2 = IW + IDO
                     IX3 = IX2 + IDO
                     IX4 = IX3 + IDO
                     IF (NA == 0) THEN
                        CALL RADB5 (IDO, L1, C, CH, WA(IW), WA(IX2),&
                                    WA(IX3), WA(IX4))
                     ELSE
                        CALL RADB5 (IDO, L1, CH, C, WA(IW), WA(IX2),&
                                    WA(IX3), WA(IX4))
                     ENDIF
                     NA = 1 - NA
                  ELSE
                     IF (NA == 0) THEN
                        CALL RADBG(IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
                     ELSE
                        CALL RADBG(IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
                     ENDIF
                     IF (IDO == 1) NA = 1 - NA
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         L1 = L2
         IW = IW + (IP - 1)*IDO
      END DO
      IF (NA == 0) goto 988
      C(:N) = CH(:N)
 988  END SUBROUTINE RFFTB1


      SUBROUTINE RADB2(IDO, L1, CC, CH, WA1)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      DOUBLE PRECISION , INTENT(IN) :: CC(IDO,2,L1)
      DOUBLE PRECISION , INTENT(OUT) :: CH(IDO,L1,2)
      DOUBLE PRECISION , INTENT(IN) :: WA1(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      DOUBLE PRECISION :: TR2, TI2
!-----------------------------------------------
      CH(1,:,1) = CC(1,1,:) + CC(IDO,2,:)
      CH(1,:,2) = CC(1,1,:) - CC(IDO,2,:)
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  CH(I-1,K,1) = CC(I-1,1,K) + CC(IC-1,2,K)
                  TR2 = CC(I-1,1,K) - CC(IC-1,2,K)
                  CH(I,K,1) = CC(I,1,K) - CC(IC,2,K)
                  TI2 = CC(I,1,K) + CC(IC,2,K)
                  CH(I-1,K,2) = WA1(I-2)*TR2 - WA1(I-1)*TI2
                  CH(I,K,2) = WA1(I-2)*TI2 + WA1(I-1)*TR2
               END DO
            END DO
            IF (MOD(IDO,2) == 1) goto 987
         ENDIF
         CH(IDO,:,1) = CC(IDO,1,:) + CC(IDO,1,:)
         CH(IDO,:,2) = -(CC(1,2,:)+CC(1,2,:))
      ENDIF
 987  END SUBROUTINE RADB2


      SUBROUTINE RADB3(IDO, L1, CC, CH, WA1, WA2)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      DOUBLE PRECISION , INTENT(IN) :: CC(IDO,3,L1)
      DOUBLE PRECISION , INTENT(OUT) :: CH(IDO,L1,3)
      DOUBLE PRECISION , INTENT(IN) :: WA1(*)
      DOUBLE PRECISION , INTENT(IN) :: WA2(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      DOUBLE PRECISION::TAUR,TAUI,TR2,CR2,CI3,TI2,CI2,CR3,DR2,DR3,DI2,DI3
!-----------------------------------------------
      DATA TAUR, TAUI/ -.5, 0.866025403784439/
      DO K = 1, L1
         TR2 = CC(IDO,2,K) + CC(IDO,2,K)
         CR2 = CC(1,1,K) + TAUR*TR2
         CH(1,K,1) = CC(1,1,K) + TR2
         CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
         CH(1,K,2) = CR2 - CI3
         CH(1,K,3) = CR2 + CI3
      END DO
      IF (IDO == 1) goto 986
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            TR2 = CC(I-1,3,K) + CC(IC-1,2,K)
            CR2 = CC(I-1,1,K) + TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K) + TR2
            TI2 = CC(I,3,K) - CC(IC,2,K)
            CI2 = CC(I,1,K) + TAUR*TI2
            CH(I,K,1) = CC(I,1,K) + TI2
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
            DR2 = CR2 - CI3
            DR3 = CR2 + CI3
            DI2 = CI2 + CR3
            DI3 = CI2 - CR3
            CH(I-1,K,2) = WA1(I-2)*DR2 - WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2 + WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3 - WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3 + WA2(I-1)*DR3
         END DO
      END DO
 986  END SUBROUTINE RADB3


      SUBROUTINE RADB4(IDO, L1, CC, CH, WA1, WA2, WA3)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      DOUBLE PRECISION , INTENT(IN) :: CC(IDO,4,L1)
      DOUBLE PRECISION , INTENT(OUT) :: CH(IDO,L1,4)
      DOUBLE PRECISION , INTENT(IN) :: WA1(*)
      DOUBLE PRECISION , INTENT(IN) :: WA2(*)
      DOUBLE PRECISION , INTENT(IN) :: WA3(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      DOUBLE PRECISION :: SQRT2, TR1, TR2, TR3, TR4, TI1, TI2, TI3, TI4, CR3, CI3,&
              CR2, CR4, CI2, CI4
!-----------------------------------------------
      DATA SQRT2/ 1.414213562373095/
      DO K = 1, L1
         TR1 = CC(1,1,K) - CC(IDO,4,K)
         TR2 = CC(1,1,K) + CC(IDO,4,K)
         TR3 = CC(IDO,2,K) + CC(IDO,2,K)
         TR4 = CC(1,3,K) + CC(1,3,K)
         CH(1,K,1) = TR2 + TR3
         CH(1,K,2) = TR1 - TR4
         CH(1,K,3) = TR2 - TR3
         CH(1,K,4) = TR1 + TR4
      END DO
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  TI1 = CC(I,1,K) + CC(IC,4,K)
                  TI2 = CC(I,1,K) - CC(IC,4,K)
                  TI3 = CC(I,3,K) - CC(IC,2,K)
                  TR4 = CC(I,3,K) + CC(IC,2,K)
                  TR1 = CC(I-1,1,K) - CC(IC-1,4,K)
                  TR2 = CC(I-1,1,K) + CC(IC-1,4,K)
                  TI4 = CC(I-1,3,K) - CC(IC-1,2,K)
                  TR3 = CC(I-1,3,K) + CC(IC-1,2,K)
                  CH(I-1,K,1) = TR2 + TR3
                  CR3 = TR2 - TR3
                  CH(I,K,1) = TI2 + TI3
                  CI3 = TI2 - TI3
                  CR2 = TR1 - TR4
                  CR4 = TR1 + TR4
                  CI2 = TI1 + TI4
                  CI4 = TI1 - TI4
                  CH(I-1,K,2) = WA1(I-2)*CR2 - WA1(I-1)*CI2
                  CH(I,K,2) = WA1(I-2)*CI2 + WA1(I-1)*CR2
                  CH(I-1,K,3) = WA2(I-2)*CR3 - WA2(I-1)*CI3
                  CH(I,K,3) = WA2(I-2)*CI3 + WA2(I-1)*CR3
                  CH(I-1,K,4) = WA3(I-2)*CR4 - WA3(I-1)*CI4
                  CH(I,K,4) = WA3(I-2)*CI4 + WA3(I-1)*CR4
               END DO
            END DO
            IF (MOD(IDO,2) == 1) goto 985
         ENDIF
         DO K = 1, L1
            TI1 = CC(1,2,K) + CC(1,4,K)
            TI2 = CC(1,4,K) - CC(1,2,K)
            TR1 = CC(IDO,1,K) - CC(IDO,3,K)
            TR2 = CC(IDO,1,K) + CC(IDO,3,K)
            CH(IDO,K,1) = TR2 + TR2
            CH(IDO,K,2) = SQRT2*(TR1 - TI1)
            CH(IDO,K,3) = TI2 + TI2
            CH(IDO,K,4) = -SQRT2*(TR1 + TI1)
         END DO
      ENDIF
 985  END SUBROUTINE RADB4


      SUBROUTINE RADB5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      DOUBLE PRECISION , INTENT(IN) :: CC(IDO,5,L1)
      DOUBLE PRECISION , INTENT(OUT) :: CH(IDO,L1,5)
      DOUBLE PRECISION , INTENT(IN) :: WA1(*)
      DOUBLE PRECISION , INTENT(IN) :: WA2(*)
      DOUBLE PRECISION , INTENT(IN) :: WA3(*)
      DOUBLE PRECISION , INTENT(IN) :: WA4(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      DOUBLE PRECISION :: TR11, TI11, TR12, TI12, TI5, TI4, TR2, TR3, CR2, CR3, CI5,&
              CI4, TI2, TI3, TR5, TR4, CI2, CI3, CR5, CR4, DR3, DR4, DI3,&
              DI4, DR5, DR2, DI5, DI2
!-----------------------------------------------
      DATA TR11, TI11, TR12, TI12/ 0.309016994374947, 0.951056516295154,&
                                 -.809016994374947, 0.587785252292473/
      DO K = 1, L1
         TI5 = CC(1,3,K) + CC(1,3,K)
         TI4 = CC(1,5,K) + CC(1,5,K)
         TR2 = CC(IDO,2,K) + CC(IDO,2,K)
         TR3 = CC(IDO,4,K) + CC(IDO,4,K)
         CH(1,K,1) = CC(1,1,K) + TR2 + TR3
         CR2 = CC(1,1,K) + TR11*TR2 + TR12*TR3
         CR3 = CC(1,1,K) + TR12*TR2 + TR11*TR3
         CI5 = TI11*TI5 + TI12*TI4
         CI4 = TI12*TI5 - TI11*TI4
         CH(1,K,2) = CR2 - CI5
         CH(1,K,3) = CR3 - CI4
         CH(1,K,4) = CR3 + CI4
         CH(1,K,5) = CR2 + CI5
      END DO
      IF (IDO == 1) goto 984
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            TI5 = CC(I,3,K) + CC(IC,2,K)
            TI2 = CC(I,3,K) - CC(IC,2,K)
            TI4 = CC(I,5,K) + CC(IC,4,K)
            TI3 = CC(I,5,K) - CC(IC,4,K)
            TR5 = CC(I-1,3,K) - CC(IC-1,2,K)
            TR2 = CC(I-1,3,K) + CC(IC-1,2,K)
            TR4 = CC(I-1,5,K) - CC(IC-1,4,K)
            TR3 = CC(I-1,5,K) + CC(IC-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K) + TR2 + TR3
            CH(I,K,1) = CC(I,1,K) + TI2 + TI3
            CR2 = CC(I-1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(I,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(I-1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(I,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            DR3 = CR3 - CI4
            DR4 = CR3 + CI4
            DI3 = CI3 + CR4
            DI4 = CI3 - CR4
            DR5 = CR2 + CI5
            DR2 = CR2 - CI5
            DI5 = CI2 - CR5
            DI2 = CI2 + CR5
            CH(I-1,K,2) = WA1(I-2)*DR2 - WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2 + WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3 - WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3 + WA2(I-1)*DR3
            CH(I-1,K,4) = WA3(I-2)*DR4 - WA3(I-1)*DI4
            CH(I,K,4) = WA3(I-2)*DI4 + WA3(I-1)*DR4
            CH(I-1,K,5) = WA4(I-2)*DR5 - WA4(I-1)*DI5
            CH(I,K,5) = WA4(I-2)*DI5 + WA4(I-1)*DR5
         END DO
      END DO
 984  END SUBROUTINE RADB5


      SUBROUTINE RADBG(IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: IP
      INTEGER , INTENT(IN) :: L1
      INTEGER , INTENT(IN) :: IDL1
      DOUBLE PRECISION , INTENT(IN) :: CC(IDO,IP,L1)
      DOUBLE PRECISION , INTENT(INOUT) :: C1(IDO,L1,IP)
      DOUBLE PRECISION , INTENT(INOUT) :: C2(IDL1,IP)
      DOUBLE PRECISION , INTENT(INOUT) :: CH(IDO,L1,IP)
      DOUBLE PRECISION , INTENT(INOUT) :: CH2(IDL1,IP)
      DOUBLE PRECISION , INTENT(IN) :: WA(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER::IDP2,NBD,IPP2,IPPH,K,I,J,JC,J2,L,LC,IS,IDIJ !!j ,IC,IK
      DOUBLE PRECISION::TPI,ARG,DCP,DSP,AR1,AI1,AR1H,DC2,DS2,AR2,AI2,AR2H !!j ,DUM
!-----------------------------------------------
      TPI = 8.0*ATAN(1.0)
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO + 2
      NBD = (IDO - 1)/2
      IPP2 = IP + 2
      IPPH = (IP + 1)/2
      IF (IDO >= L1) THEN
         CH(:,:,1) = CC(:,1,:)
      ELSE
         CH(:,:,1) = CC(:,1,:)
      ENDIF
      DO J = 2, IPPH
         JC = IPP2 - J
         J2 = J + J
         CH(1,:,J) = CC(IDO,J2-2,:) + CC(IDO,J2-2,:)
         CH(1,:,JC) = CC(1,J2-1,:) + CC(1,J2-1,:)
      END DO
      IF (IDO /= 1) THEN
         IF (NBD >= L1) THEN
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = CC(2:IDO-1:2,2*J-1,:) + & 
                                   CC(IDP2-4:IDP2-1-IDO:(-2),2*J-2,:)
               CH(2:IDO-1:2,:,JC) = CC(2:IDO-1:2,2*J-1,:) - & 
                                    CC(IDP2-4:IDP2-1-IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,J) = CC(3:IDO:2,2*J-1,:) - &
                                 CC(IDP2-3:IDP2-IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,JC) = CC(3:IDO:2,2*J-1,:) + &
                                  CC(IDP2-3:IDP2-IDO:(-2),2*J-2,:)
            END DO
         ELSE
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = CC(2:IDO-1:2,2*J-1,:) + &
                                   CC(IDP2-4:IDP2-1-IDO:(-2),2*J-2,:)
               CH(2:IDO-1:2,:,JC) = CC(2:IDO-1:2,2*J-1,:) - &
                                    CC(IDP2-4:IDP2-1-IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,J) = CC(3:IDO:2,2*J-1,:) - &
                                 CC(IDP2-3:IDP2-IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,JC) = CC(3:IDO:2,2*J-1,:) + &
                                  CC(IDP2-3:IDP2-IDO:(-2),2*J-2,:)
            END DO
         ENDIF
      ENDIF
      AR1 = 1.
      AI1 = 0.
      DO L = 2, IPPH
         LC = IPP2 - L
         AR1H = DCP*AR1 - DSP*AI1
         AI1 = DCP*AI1 + DSP*AR1
         AR1 = AR1H
         C2(:,L) = CH2(:,1) + AR1*CH2(:,2)
         C2(:,LC) = AI1*CH2(:,IP)
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO J = 3, IPPH
            JC = IPP2 - J
            AR2H = DC2*AR2 - DS2*AI2
            AI2 = DC2*AI2 + DS2*AR2
            AR2 = AR2H
            C2(:,L) = C2(:,L) + AR2*CH2(:,J)
            C2(:,LC) = C2(:,LC) + AI2*CH2(:,JC)
         END DO
      END DO
      DO J = 2, IPPH
         CH2(:,1) = CH2(:,1) + CH2(:,J)
      END DO
      DO J = 2, IPPH
         JC = IPP2 - J
         CH(1,:,J) = C1(1,:,J) - C1(1,:,JC)
         CH(1,:,JC) = C1(1,:,J) + C1(1,:,JC)
      END DO
      IF (IDO /= 1) THEN
         IF (NBD >= L1) THEN
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = C1(2:IDO-1:2,:,J) - C1(3:IDO:2,:,JC)
               CH(2:IDO-1:2,:,JC) = C1(2:IDO-1:2,:,J) + C1(3:IDO:2,:,JC)
               CH(3:IDO:2,:,J) = C1(3:IDO:2,:,J) + C1(2:IDO-1:2,:,JC)
               CH(3:IDO:2,:,JC) = C1(3:IDO:2,:,J) - C1(2:IDO-1:2,:,JC)
            END DO
         ELSE
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = C1(2:IDO-1:2,:,J) - C1(3:IDO:2,:,JC)
               CH(2:IDO-1:2,:,JC) = C1(2:IDO-1:2,:,J) + C1(3:IDO:2,:,JC)
               CH(3:IDO:2,:,J) = C1(3:IDO:2,:,J) + C1(2:IDO-1:2,:,JC)
               CH(3:IDO:2,:,JC) = C1(3:IDO:2,:,J) - C1(2:IDO-1:2,:,JC)
            END DO
         ENDIF
      ENDIF
      IF (IDO == 1) goto 983
      C2(:,1) = CH2(:,1)
      C1(1,:,2:IP) = CH(1,:,2:IP)
      IF (NBD <= L1) THEN
         IS = -IDO
         DO J = 2, IP
            IS = IS + IDO
            IDIJ = IS
            DO I = 3, IDO, 2
               IDIJ = IDIJ + 2
               C1(I-1,:,J) = WA(IDIJ-1)*CH(I-1,:,J) - WA(IDIJ)*CH(I,:,J)
               C1(I,:,J) = WA(IDIJ-1)*CH(I,:,J) + WA(IDIJ)*CH(I-1,:,J)
            END DO
         END DO
      ELSE
         IS = -IDO
         DO J = 2, IP
            IS = IS + IDO
            DO K = 1, L1
               IDIJ = IS
               C1(2:IDO-1:2,K,J) = WA(IDIJ+1:IDO-2+IDIJ:2)*CH(2:IDO-1:2,K,J)-&
                                   WA(IDIJ+2:IDO-1+IDIJ:2)*CH(3:IDO:2,K,J)
               C1(3:IDO:2,K,J) = WA(IDIJ+1:IDO-2+IDIJ:2)*CH(3:IDO:2,K,J) + &
                                 WA(IDIJ+2:IDO-1+IDIJ:2)*CH(2:IDO-1:2,K,J)
            END DO
         END DO
      ENDIF
 983  END SUBROUTINE RADBG


      SUBROUTINE RFFTF(N, R, WSAVE)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      DOUBLE PRECISION  :: R(*)
      DOUBLE PRECISION  :: WSAVE(*)
!-----------------------------------------------
!
      IF (N == 1) goto 982
      CALL RFFTF1 (N, R, WSAVE, WSAVE(N+1), WSAVE(2*N+1))
 982  END SUBROUTINE RFFTF


      SUBROUTINE RFFTF1(N, C, CH, WA, IFAC)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: IFAC(*)
      DOUBLE PRECISION  :: C(*)
      DOUBLE PRECISION  :: CH(*)
      DOUBLE PRECISION  :: WA(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER::NF,NA,L2,IW,K1,KH,IP,L1,IDO,IDL1,IX2,IX3,IX4 !!j ,I
!-----------------------------------------------
      NF = IFAC(2)
      NA = 1
      L2 = N
      IW = N
      DO K1 = 1, NF
         KH = NF - K1
         IP = IFAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW - (IP - 1)*IDO
         NA = 1 - NA
         IF (IP == 4) THEN
            IX2 = IW + IDO
            IX3 = IX2 + IDO
            IF (NA == 0) THEN
               CALL RADF4 (IDO, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
               GO TO 110
            ENDIF
            CALL RADF4 (IDO, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
            GO TO 110
         ENDIF
         IF (IP == 2) THEN
            IF (NA == 0) THEN
               CALL RADF2 (IDO, L1, C, CH, WA(IW))
               GO TO 110
            ENDIF
            CALL RADF2 (IDO, L1, CH, C, WA(IW))
            GO TO 110
         ENDIF
  104    CONTINUE
         IF (IP == 3) THEN
            IX2 = IW + IDO
            IF (NA == 0) THEN
               CALL RADF3 (IDO, L1, C, CH, WA(IW), WA(IX2))
               GO TO 110
            ENDIF
            CALL RADF3 (IDO, L1, CH, C, WA(IW), WA(IX2))
            GO TO 110
         ENDIF
  106    CONTINUE
         IF (IP == 5) THEN
            IX2 = IW + IDO
            IX3 = IX2 + IDO
            IX4 = IX3 + IDO
            IF (NA == 0) THEN
               CALL RADF5(IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
               GO TO 110
            ENDIF
            CALL RADF5(IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
            GO TO 110
         ENDIF
  108    CONTINUE
         IF (IDO == 1) NA = 1 - NA
         IF (NA == 0) THEN
            CALL RADFG (IDO, IP, L1, IDL1, C, C, C, CH, CH, WA(IW))
            NA = 1
         ELSE
            CALL RADFG (IDO, IP, L1, IDL1, CH, CH, CH, C, C, WA(IW))
            NA = 0
         ENDIF
  110    CONTINUE
         L2 = L1
      END DO
      IF (NA == 1) goto 981
      C(:N) = CH(:N)
 981  END SUBROUTINE RFFTF1


      SUBROUTINE RADF2(IDO, L1, CC, CH, WA1)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      DOUBLE PRECISION , INTENT(IN) :: CC(IDO,L1,2)
      DOUBLE PRECISION , INTENT(OUT) :: CH(IDO,2,L1)
      DOUBLE PRECISION , INTENT(IN) :: WA1(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      DOUBLE PRECISION :: TR2, TI2
!-----------------------------------------------
      CH(1,1,:) = CC(1,:,1) + CC(1,:,2)
      CH(IDO,2,:) = CC(1,:,1) - CC(1,:,2)
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  TR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
                  TI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
                  CH(I,1,K) = CC(I,K,1) + TI2
                  CH(IC,2,K) = TI2 - CC(I,K,1)
                  CH(I-1,1,K) = CC(I-1,K,1) + TR2
                  CH(IC-1,2,K) = CC(I-1,K,1) - TR2
               END DO
            END DO
            IF (MOD(IDO,2) == 1) goto 980
         ENDIF
         CH(1,2,:) = -CC(IDO,:,2)
         CH(IDO,1,:) = CC(IDO,:,1)
      ENDIF
 980  END SUBROUTINE RADF2


      SUBROUTINE RADF3(IDO, L1, CC, CH, WA1, WA2)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      DOUBLE PRECISION , INTENT(IN) :: CC(IDO,L1,3)
      DOUBLE PRECISION , INTENT(OUT) :: CH(IDO,3,L1)
      DOUBLE PRECISION , INTENT(IN) :: WA1(*)
      DOUBLE PRECISION , INTENT(IN) :: WA2(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      DOUBLE PRECISION::TAUR,TAUI,CR2,DR2,DI2,DR3,DI3,CI2,TR2,TI2,TR3,TI3
!-----------------------------------------------
      DATA TAUR, TAUI/ -.5, 0.866025403784439/
      DO K = 1, L1
         CR2 = CC(1,K,2) + CC(1,K,3)
         CH(1,1,K) = CC(1,K,1) + CR2
         CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
         CH(IDO,2,K) = CC(1,K,1) + TAUR*CR2
      END DO
      IF (IDO == 1) goto 979
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            DR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
            CR2 = DR2 + DR3
            CI2 = DI2 + DI3
            CH(I-1,1,K) = CC(I-1,K,1) + CR2
            CH(I,1,K) = CC(I,K,1) + CI2
            TR2 = CC(I-1,K,1) + TAUR*CR2
            TI2 = CC(I,K,1) + TAUR*CI2
            TR3 = TAUI*(DI2 - DI3)
            TI3 = TAUI*(DR3 - DR2)
            CH(I-1,3,K) = TR2 + TR3
            CH(IC-1,2,K) = TR2 - TR3
            CH(I,3,K) = TI2 + TI3
            CH(IC,2,K) = TI3 - TI2
         END DO
      END DO
 979  END SUBROUTINE RADF3


      SUBROUTINE RADF4(IDO, L1, CC, CH, WA1, WA2, WA3)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      DOUBLE PRECISION , INTENT(IN) :: CC(IDO,L1,4)
      DOUBLE PRECISION , INTENT(OUT) :: CH(IDO,4,L1)
      DOUBLE PRECISION , INTENT(IN) :: WA1(*)
      DOUBLE PRECISION , INTENT(IN) :: WA2(*)
      DOUBLE PRECISION , INTENT(IN) :: WA3(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      DOUBLE PRECISION :: HSQT2, TR1, TR2, CR2, CI2, CR3, CI3, CR4, CI4, TR4, TI1,&
              TI4, TI2, TI3, TR3
!-----------------------------------------------
      DATA HSQT2/ 0.7071067811865475/
      DO K = 1, L1
         TR1 = CC(1,K,2) + CC(1,K,4)
         TR2 = CC(1,K,1) + CC(1,K,3)
         CH(1,1,K) = TR1 + TR2
         CH(IDO,4,K) = TR2 - TR1
         CH(IDO,2,K) = CC(1,K,1) - CC(1,K,3)
         CH(1,3,K) = CC(1,K,4) - CC(1,K,2)
      END DO
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  CR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
                  CI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
                  CR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
                  CI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
                  CR4 = WA3(I-2)*CC(I-1,K,4) + WA3(I-1)*CC(I,K,4)
                  CI4 = WA3(I-2)*CC(I,K,4) - WA3(I-1)*CC(I-1,K,4)
                  TR1 = CR2 + CR4
                  TR4 = CR4 - CR2
                  TI1 = CI2 + CI4
                  TI4 = CI2 - CI4
                  TI2 = CC(I,K,1) + CI3
                  TI3 = CC(I,K,1) - CI3
                  TR2 = CC(I-1,K,1) + CR3
                  TR3 = CC(I-1,K,1) - CR3
                  CH(I-1,1,K) = TR1 + TR2
                  CH(IC-1,4,K) = TR2 - TR1
                  CH(I,1,K) = TI1 + TI2
                  CH(IC,4,K) = TI1 - TI2
                  CH(I-1,3,K) = TI4 + TR3
                  CH(IC-1,2,K) = TR3 - TI4
                  CH(I,3,K) = TR4 + TI3
                  CH(IC,2,K) = TR4 - TI3
               END DO
            END DO
            IF (MOD(IDO,2) == 1) goto 978
         ENDIF
         DO K = 1, L1
            TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
            TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
            CH(IDO,1,K) = TR1 + CC(IDO,K,1)
            CH(IDO,3,K) = CC(IDO,K,1) - TR1
            CH(1,2,K) = TI1 - CC(IDO,K,3)
            CH(1,4,K) = TI1 + CC(IDO,K,3)
         END DO
      ENDIF
 978  END SUBROUTINE RADF4


      SUBROUTINE RADF5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      DOUBLE PRECISION , INTENT(IN) :: CC(IDO,L1,5)
      DOUBLE PRECISION , INTENT(OUT) :: CH(IDO,5,L1)
      DOUBLE PRECISION , INTENT(IN) :: WA1(*)
      DOUBLE PRECISION , INTENT(IN) :: WA2(*)
      DOUBLE PRECISION , INTENT(IN) :: WA3(*)
      DOUBLE PRECISION , INTENT(IN) :: WA4(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      DOUBLE PRECISION :: TR11, TI11, TR12, TI12, CR2, CI5, CR3, CI4, DR2, DI2, DR3,&
              DI3, DR4, DI4, DR5, DI5, CR5, CI2, CR4, CI3, TR2, TI2, TR3,&
              TI3, TR5, TI5, TR4, TI4
!-----------------------------------------------
      DATA TR11, TI11, TR12, TI12/ 0.309016994374947, 0.951056516295154,&
                                 -.809016994374947, 0.587785252292473/
      DO K = 1, L1
         CR2 = CC(1,K,5) + CC(1,K,2)
         CI5 = CC(1,K,5) - CC(1,K,2)
         CR3 = CC(1,K,4) + CC(1,K,3)
         CI4 = CC(1,K,4) - CC(1,K,3)
         CH(1,1,K) = CC(1,K,1) + CR2 + CR3
         CH(IDO,2,K) = CC(1,K,1) + TR11*CR2 + TR12*CR3
         CH(1,3,K) = TI11*CI5 + TI12*CI4
         CH(IDO,4,K) = CC(1,K,1) + TR12*CR2 + TR11*CR3
         CH(1,5,K) = TI12*CI5 - TI11*CI4
      END DO
      IF (IDO == 1) goto 977
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            DR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
            DR4 = WA3(I-2)*CC(I-1,K,4) + WA3(I-1)*CC(I,K,4)
            DI4 = WA3(I-2)*CC(I,K,4) - WA3(I-1)*CC(I-1,K,4)
            DR5 = WA4(I-2)*CC(I-1,K,5) + WA4(I-1)*CC(I,K,5)
            DI5 = WA4(I-2)*CC(I,K,5) - WA4(I-1)*CC(I-1,K,5)
            CR2 = DR2 + DR5
            CI5 = DR5 - DR2
            CR5 = DI2 - DI5
            CI2 = DI2 + DI5
            CR3 = DR3 + DR4
            CI4 = DR4 - DR3
            CR4 = DI3 - DI4
            CI3 = DI3 + DI4
            CH(I-1,1,K) = CC(I-1,K,1) + CR2 + CR3
            CH(I,1,K) = CC(I,K,1) + CI2 + CI3
            TR2 = CC(I-1,K,1) + TR11*CR2 + TR12*CR3
            TI2 = CC(I,K,1) + TR11*CI2 + TR12*CI3
            TR3 = CC(I-1,K,1) + TR12*CR2 + TR11*CR3
            TI3 = CC(I,K,1) + TR12*CI2 + TR11*CI3
            TR5 = TI11*CR5 + TI12*CR4
            TI5 = TI11*CI5 + TI12*CI4
            TR4 = TI12*CR5 - TI11*CR4
            TI4 = TI12*CI5 - TI11*CI4
            CH(I-1,3,K) = TR2 + TR5
            CH(IC-1,2,K) = TR2 - TR5
            CH(I,3,K) = TI2 + TI5
            CH(IC,2,K) = TI5 - TI2
            CH(I-1,5,K) = TR3 + TR4
            CH(IC-1,4,K) = TR3 - TR4
            CH(I,5,K) = TI3 + TI4
            CH(IC,4,K) = TI4 - TI3
         END DO
      END DO
 977  END SUBROUTINE RADF5


      SUBROUTINE RADFG(IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
!-----------------------------------------------
! D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: IP
      INTEGER , INTENT(IN) :: L1
      INTEGER , INTENT(IN) :: IDL1
      DOUBLE PRECISION , INTENT(OUT) :: CC(IDO,IP,L1)
      DOUBLE PRECISION , INTENT(INOUT) :: C1(IDO,L1,IP)
      DOUBLE PRECISION , INTENT(INOUT) :: C2(IDL1,IP)
      DOUBLE PRECISION , INTENT(INOUT) :: CH(IDO,L1,IP)
      DOUBLE PRECISION , INTENT(INOUT) :: CH2(IDL1,IP)
      DOUBLE PRECISION , INTENT(IN) :: WA(*)
!-----------------------------------------------
! L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER::IPPH,IPP2,IDP2,NBD,J,K,IS,IDIJ,I,JC,L,LC!!j ,IK,J2,IC 
      DOUBLE PRECISION::TPI,ARG,DCP,DSP,AR1,AI1,AR1H,DC2,DS2,AR2,AI2,AR2H !!j ,DUM
!-----------------------------------------------
      TPI = 8.0*ATAN(1.0)
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP + 1)/2
      IPP2 = IP + 2
      IDP2 = IDO + 2
      NBD = (IDO - 1)/2
      IF (IDO /= 1) THEN
         CH2(:,1) = C2(:,1)
         CH(1,:,2:IP) = C1(1,:,2:IP)
         IF (NBD <= L1) THEN
            IS = -IDO
            DO J = 2, IP
               IS = IS + IDO
               IDIJ = IS
               DO I = 3, IDO, 2
                  IDIJ = IDIJ + 2
                  CH(I-1,:,J)=WA(IDIJ-1)*C1(I-1,:,J)+WA(IDIJ)*C1(I,:,J)
                  CH(I,:,J)=WA(IDIJ-1)*C1(I,:,J)-WA(IDIJ)*C1(I-1,:,J)
               END DO
            END DO
         ELSE
            IS = -IDO
            DO J = 2, IP
               IS = IS + IDO
               DO K = 1, L1
                  IDIJ = IS
                  CH(2:IDO-1:2,K,J) = WA(IDIJ+1:IDO-2+IDIJ:2)*C1(2:IDO-1:2,K,J) + &
                                      WA(IDIJ+2:IDO-1+IDIJ:2)*C1(3:IDO:2,K,J)
                  CH(3:IDO:2,K,J) = WA(IDIJ+1:IDO-2+IDIJ:2)*C1(3:IDO:2,K,J) - &
                                    WA(IDIJ+2:IDO-1+IDIJ:2)*C1(2:IDO-1:2,K,J)
               END DO
            END DO
         ENDIF
         IF (NBD >= L1) THEN
            DO J = 2, IPPH
               JC = IPP2 - J
               C1(2:IDO-1:2,:,J)=CH(2:IDO-1:2,:,J)+CH(2:IDO-1:2,:,JC)
               C1(2:IDO-1:2,:,JC) = CH(3:IDO:2,:,J) - CH(3:IDO:2,:,JC)
               C1(3:IDO:2,:,J) = CH(3:IDO:2,:,J) + CH(3:IDO:2,:,JC)
               C1(3:IDO:2,:,JC) = CH(2:IDO-1:2,:,JC) - CH(2:IDO-1:2,:,J)
            END DO
            GO TO 121
         ENDIF
         DO J = 2, IPPH
            JC = IPP2 - J
            C1(2:IDO-1:2,:,J) = CH(2:IDO-1:2,:,J) + CH(2:IDO-1:2,:,JC)
            C1(2:IDO-1:2,:,JC) = CH(3:IDO:2,:,J) - CH(3:IDO:2,:,JC)
            C1(3:IDO:2,:,J) = CH(3:IDO:2,:,J) + CH(3:IDO:2,:,JC)
            C1(3:IDO:2,:,JC) = CH(2:IDO-1:2,:,JC) - CH(2:IDO-1:2,:,J)
         END DO
         GO TO 121
      ENDIF
      C2(:,1) = CH2(:,1)
  121 CONTINUE
      DO J = 2, IPPH
         JC = IPP2 - J
         C1(1,:,J) = CH(1,:,J) + CH(1,:,JC)
         C1(1,:,JC) = CH(1,:,JC) - CH(1,:,J)
      END DO
!
      AR1 = 1.
      AI1 = 0.
      DO L = 2, IPPH
         LC = IPP2 - L
         AR1H = DCP*AR1 - DSP*AI1
         AI1 = DCP*AI1 + DSP*AR1
         AR1 = AR1H
         CH2(:,L) = C2(:,1) + AR1*C2(:,2)
         CH2(:,LC) = AI1*C2(:,IP)
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO J = 3, IPPH
            JC = IPP2 - J
            AR2H = DC2*AR2 - DS2*AI2
            AI2 = DC2*AI2 + DS2*AR2
            AR2 = AR2H
            CH2(:,L) = CH2(:,L) + AR2*C2(:,J)
            CH2(:,LC) = CH2(:,LC) + AI2*C2(:,JC)
         END DO
      END DO
      DO J = 2, IPPH
         CH2(:,1) = CH2(:,1) + C2(:,J)
      END DO
!
      IF (IDO >= L1) THEN
         CC(:,1,:) = CH(:,:,1)
      ELSE
         CC(:,1,:) = CH(:,:,1)
      ENDIF
      CC(IDO,2:(IPPH-1)*2:2,:) = TRANSPOSE(CH(1,:,2:IPPH))
      CC(1,3:IPPH*2-1:2,:) = TRANSPOSE(CH(1,:,IPP2-2:IPP2-IPPH:(-1)))
      IF (IDO == 1) goto 976
      IF (NBD >= L1) THEN
         CC(2:IDO-1:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(2:IDO-1:2,:,&
            2:IPPH)+CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO&
            -1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
         CC(IDP2-4:IDP2-1-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE =&
            CH(2:IDO-1:2,:,2:IPPH)-CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1))&
            ,SHAPE = (/(IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
         CC(3:IDO:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(3:IDO:2,:,2:&
            IPPH)+CH(3:IDO:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO-1)/&
            2,IPPH-1,L1/),ORDER = (/1,3,2/))
         CC(IDP2-3:IDP2-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE = CH&
            (3:IDO:2,:,IPP2-2:IPP2-IPPH:(-1))-CH(3:IDO:2,:,2:IPPH),SHAPE&
             = (/(IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
         goto 976
      ENDIF
      CC(2:IDO-1:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(2:IDO-1:2,:,2:&
         IPPH)+CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO-1)/2&
         ,IPPH-1,L1/),ORDER = (/1,3,2/))
      CC(IDP2-4:IDP2-1-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE = CH(&
         2:IDO-1:2,:,2:IPPH)-CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE&
          = (/(IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
      CC(3:IDO:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(3:IDO:2,:,2:IPPH)&
         +CH(3:IDO:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO-1)/2,IPPH-1&
         ,L1/),ORDER = (/1,3,2/))
      CC(IDP2-3:IDP2-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE = CH(3:&
         IDO:2,:,IPP2-2:IPP2-IPPH:(-1))-CH(3:IDO:2,:,2:IPPH),SHAPE = (/(&
         IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
 976  END SUBROUTINE RADFG
!
!
!
!=======================================================================
      SUBROUTINE SORTRX_pm(INDEX,N,DATA)!(N,DATA,INDEX)
!=======================================================================
!
!     SORTRX -- SORT, Real input, indeX output
!
!
!     Input:  N     INTEGER
!             DATA  DOUBLE PRECISION
!
!     Output: INDEX INTEGER (DIMENSION N)
!
! This routine performs an in-memory sort of the first N elements of
! array DATA, returning into array INDEX the indices of elements of
! DATA arranged in ascending order.  Thus,
!
!    DATA(INDEX(1)) will be the smallest number in array DATA;
!    DATA(INDEX(N)) will be the largest number in DATA.
!
! The original data is not physically rearranged.  The original order
! of equal input values is not necessarily preserved.
!
!===================================================================
!
! SORTRX uses a hybrid QuickSort algorithm, based on several
! suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
! "pivot key" [my term] for dividing each subsequence is chosen to be
! the median of the first, last, and middle values of the subsequence;
! and the QuickSort is cut off when a subsequence has 9 or fewer
! elements, and a straight insertion sort of the entire array is done
! at the end.  The result is comparable to a pure insertion sort for
! very short arrays, and very fast for very large arrays (of order 12
! micro-sec/element on the 3081K for arrays of 10K elements).  It is
! also not subject to the poor performance of the pure QuickSort on
! partially ordered data.
!
! Created:  15 Jul 1986  Len Moss
!
!===================================================================

!      Dimension INDEX(N),DATA(N)
      IMPLICIT NONE
      INTEGER :: N,INDEX(N)
      DOUBLE PRECISION    :: DATA(N)

      INTEGER  :: LSTK(31),RSTK(31),ISTK
      INTEGER  :: L,R,I,J,P,INDEXP,INDEXT
      DOUBLE PRECISION     :: DATAP

!     QuickSort Cutoff
!
!     Quit QuickSort-ing when a subsequence contains M or fewer
!     elements and finish off at end with straight insertion sort.
!     According to Knuth, V.3, the optimum value of M is around 9.

      INTEGER , PARAMETER ::M=9

!===================================================================
!
!     Make initial guess for INDEX

      DO 50 I=1,N
         INDEX(I)=I
   50    CONTINUE

!     If array is short, skip QuickSort and go directly to
!     the straight insertion sort.

      IF (N.LE.M) GOTO 900

!===================================================================
!
!     QuickSort
!
!     The "Qn:"s correspond roughly to steps in Algorithm Q,
!     Knuth, V.3, PP.116-117, modified to select the median
!     of the first, last, and middle elements as the "pivot
!     key" (in Knuth's notation, "K").  Also modified to leave
!     data in place and produce an INDEX array.  To simplify
!     comments, let DATA[I]=DATA(INDEX(I)).

! Q1: Initialize
      ISTK=0
      L=1
      R=N

  200 CONTINUE

! Q2: Sort the subsequence DATA[L]..DATA[R].
!
!     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
!     r > R, and L <= m <= R.  (First time through, there is no
!     DATA for l < L or r > R.)

      I=L
      J=R

! Q2.5: Select pivot key
!
!     Let the pivot, P, be the midpoint of this subsequence,
!     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
!     so the corresponding DATA values are in increasing order.
!     The pivot key, DATAP, is then DATA[P].

      P=(L+R)/2
      INDEXP=INDEX(P)
      DATAP=DATA(INDEXP)

      IF (DATA(INDEX(L)) .GT. DATAP) THEN
         INDEX(P)=INDEX(L)
         INDEX(L)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF

      IF (DATAP .GT. DATA(INDEX(R))) THEN
         IF (DATA(INDEX(L)) .GT. DATA(INDEX(R))) THEN
            INDEX(P)=INDEX(L)
            INDEX(L)=INDEX(R)
         ELSE
            INDEX(P)=INDEX(R)
         ENDIF
         INDEX(R)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF

!     Now we swap values between the right and left sides and/or
!     move DATAP until all smaller values are on the left and all
!     larger values are on the right.  Neither the left or right
!     side will be internally ordered yet; however, DATAP will be
!     in its final position.

  300 CONTINUE

! Q3: Search for datum on left >= DATAP
!
!     At this point, DATA[L] <= DATAP.  We can therefore start scanning
!     up from L, looking for a value >= DATAP (this scan is guaranteed
!     to terminate since we initially placed DATAP near the middle of
!     the subsequence).

         I=I+1
         IF (DATA(INDEX(I)).LT.DATAP) GOTO 300

  400 CONTINUE

! Q4: Search for datum on right <= DATAP
!
!     At this point, DATA[R] >= DATAP.  We can therefore start scanning
!     down from R, looking for a value <= DATAP (this scan is guaranteed
!     to terminate since we initially placed DATAP near the middle of
!     the subsequence).

         J=J-1
         IF (DATA(INDEX(J)).GT.DATAP) GOTO 400

! Q5: Have the two scans collided?

      IF (I.LT.J) THEN

! Q6: No, interchange DATA[I] <--> DATA[J] and continue

         INDEXT=INDEX(I)
         INDEX(I)=INDEX(J)
         INDEX(J)=INDEXT
         GOTO 300
      ELSE

! Q7: Yes, select next subsequence to sort
!
!     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
!     for all L <= l < I and J < r <= R.  If both subsequences are
!     more than M elements long, push the longer one on the stack and
!     go back to QuickSort the shorter; if only one is more than M
!     elements long, go back and QuickSort it; otherwise, pop a
!     subsequence off the stack and QuickSort it.

         IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=J+1
            RSTK(ISTK)=R
            R=I-1
         ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=L
            RSTK(ISTK)=I-1
            L=J+1
         ELSE IF (R-J .GT. M) THEN
            L=J+1
         ELSE IF (I-L .GT. M) THEN
            R=I-1
         ELSE
! Q8: Pop the stack, or terminate QuickSort if empty
            IF (ISTK.LT.1) GOTO 900
            L=LSTK(ISTK)
            R=RSTK(ISTK)
            ISTK=ISTK-1
         ENDIF
         GOTO 200
      ENDIF

  900 CONTINUE

!===================================================================
!
! Q9: Straight Insertion sort

      DO 950 I=2,N
         IF (DATA(INDEX(I-1)) .GT. DATA(INDEX(I))) THEN
            INDEXP=INDEX(I)
            DATAP=DATA(INDEXP)
            P=I-1
  920       CONTINUE
               INDEX(P+1) = INDEX(P)
               P=P-1
               IF (P.GT.0) THEN
                  IF (DATA(INDEX(P)).GT.DATAP) GOTO 920
               ENDIF
            INDEX(P+1) = INDEXP
         ENDIF
  950    CONTINUE

!===================================================================
!
!     All done

      END
