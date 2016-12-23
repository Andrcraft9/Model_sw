C======================================================================
C MODULE CONTAINS THE SWEEP PROCEDURES NEEDED FOR DIFFERENT CASES
C======================================================================
      SUBROUTINE FACTOR(ND,A,B,C,ETA,RKSI,II,JJ)
      IMPLICIT NONE
C----------------------------------------------------------------------
C UP-DOWN : A(N)*F(N-1) + B(N)*F(N) + C(N)*F(N+1) = ETA(N) ; N=II,JJ
C A(II) AND C(JJ) ARE NOT USED, SO THEY MAY BE SET = 0.
C OPTIMIZED BY DIANSKI N.A.
      INCLUDE '0COM.INC'
      INTEGER II, JJ, ND
      REAL A(ND),B(ND),C(ND), ETA(ND), RKSI(ND), X(MDIM), Y(MDIM)
      INTEGER J1, JJ1
      REAL SCR

      X(II) = - C(II) / B(II)
      Y(II) = ETA(II) / B(II)
      JJ1 = II+1
      DO  J1=JJ1,JJ
       SCR=1.0/(A(J1) * X(J1-1) + B(J1))
       X(J1) = -C(J1) * SCR
       Y(J1) = (ETA(J1) - A(J1)*Y(J1-1)) * SCR
      END DO
      RKSI(JJ)=Y(JJ)
      JJ1 = JJ - 1
      DO J1=JJ1,II,-1
       RKSI(J1) = X(J1) * RKSI(J1+1) + Y(J1)
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE FACTOR8(ND,A,B,C,ETA,RKSI,II,JJ)
      IMPLICIT NONE
C----------------------------------------------------------------------
C UP-DOWN : A(N)*F(N-1) + B(N)*F(N) + C(N)*F(N+1) = ETA(N) ; N=II,JJ
C A(II) AND C(JJ) ARE NOT USED, SO THEY MAY BE SET = 0.
C OPTIMIZED BY DIANSKI N.A.
      INCLUDE '0COM.INC'
      INTEGER II, JJ, ND
      REAL*8 A(ND),B(ND),C(ND), ETA(ND), RKSI(ND), X(MDIM), Y(MDIM)
      INTEGER J1, JJ1
      REAL*8 SCR

      X(II) = - C(II) / B(II)
      Y(II) = ETA(II) / B(II)
      JJ1 = II+1
      DO  J1=JJ1,JJ
       SCR=1.0/(A(J1) * X(J1-1) + B(J1))
       X(J1) = -C(J1) * SCR
       Y(J1) = (ETA(J1) - A(J1)*Y(J1-1)) * SCR
      END DO
      RKSI(JJ)=Y(JJ)
      JJ1 = JJ - 1
      DO J1=JJ1,II,-1
       RKSI(J1) = X(J1) * RKSI(J1+1) + Y(J1)
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE FACTOR2 (ND,A,B,C,ETA,RKSI,II,JJM)
      IMPLICIT NONE
C----------------------------------------------------------------------
C UP-DOWN : A(N)*F(N-1) + B(N)*F(N) + C(N)*F(N+1) = ETA(N) ; N=II,JJ
C A(II) AND C(JJ) ARE NOT USED, SO THEY ARE SET = 0.
      INCLUDE '0COM.INC'
      INTEGER II, JJ, JJ1, JP1, JJM, ND
      REAL    A(ND),B(ND),C(ND),ETA(ND),RKSI(ND),X(MDIM),Y(MDIM)
      INTEGER JLOOP, J1, JM1

      JJ = JJM
      IF (JJ.GT.MM) JJ = JJ - MMD
      A(II) = 0.
      C(JJ) = 0.
      X(II) = - C(II) / B(II)
      Y(II) = ETA(II) / B(II)
      DO JLOOP=II+1,JJM
       J1  = JLOOP
       JM1 = JLOOP-1
       IF (J1 .GT.MM) J1  = J1  - MMD
       IF (JM1.GT.MM) JM1 = JM1 - MMD
       X(J1) = -C(J1) / (A(J1) * X(JM1) + B(J1))
       Y(J1) = (ETA(J1) - A(J1)*Y(JM1)) / (A(J1) * X(JM1) + B(J1))
      ENDDO
      RKSI(JJ)=Y(JJ)
      DO JLOOP=JJM-1,II,-1
       JJ1 = JLOOP
       JP1 = JJ1 + 1
       IF (JJ1.GT.MM) JJ1 = JJ1 - MMD
       IF (JP1.GT.MM) JP1 = JP1 - MMD
       RKSI(JJ1) = X(JJ1) * RKSI(JP1) + Y(JJ1)
      ENDDO
      RETURN
      END
C=====================================================06-23-97 06:25pm
      SUBROUTINE FACTORC(A,B,C,F,X,II,JJ)
      IMPLICIT NONE
C---------------------------------------------------------------------
C "CYCLICHESKAYA PROGONKA", PROGRAMMED BY N.GALKIN.
C     A, B, C :    diagonals of matrix   |                           |
C     F       :    right part            | Bii    Cii     0   ... Aii|
C     X       :    solution              | Aii+1  Bii+1   Cii+1 ...0 |
C     II      :    begin of the vectors  |  0     Aii+2   B_ \   ..0 |
C     JJ      :    end of the vectors    |  0     0..   \   \  \  .0 |
C P,Q,ALPHA,BETA,GAMMA : work arrays     |  ..............\   \  \...|
C     II <= JJ-2 : 3 points.             |  0               \   \ C_ |
C                                        |  Cjj   0 ...       Ajj Bjj|
C---------------------------------------------------------------------
      INCLUDE '0COM.INC'
      INTEGER II, JJ
      REAL A(JJ),B(JJ),C(JJ),F(JJ),X(JJ)
      REAL SCR1, P(MM),Q(MM),ALPHA(MM+1),BETA(MM+1),GAMMA(MM+1)
      INTEGER I
C---------------------------------------------------------------------
      ALPHA(II+1) = -C(II)/B(II)
      BETA (II+1) =  F(II)/B(II)
      GAMMA(II+1) = -A(II)/B(II)
      DO I=II+1,JJ
         SCR1 = -B(I) - A(I)*ALPHA(I)
         ALPHA(I+1) = C(I) / SCR1
         BETA (I+1) = (-F(I) + A(I)*BETA(I)) / SCR1
         GAMMA(I+1) = A(I)*GAMMA(I) / SCR1
      ENDDO
      P(JJ-1) = BETA(JJ)
      Q(JJ-1) = ALPHA(JJ) + GAMMA(JJ)
      DO I=JJ-2,II,-1
         P(I) = ALPHA(I+1)*P(I+1) + BETA (I+1)
         Q(I) = ALPHA(I+1)*Q(I+1) + GAMMA(I+1)
      ENDDO
      X(JJ) = (BETA(JJ+1) + ALPHA(JJ+1)*P(II)) /
     &        (1 - ALPHA(JJ+1)*Q(II) - GAMMA(JJ+1))
      DO I=II,JJ-1
         X(I) = P(I) + X(JJ)*Q(I)
      ENDDO
      RETURN
      END
C======================================================================
C SWEEPS FROM ZALESNY MODEL
C=====================================================================
      SUBROUTINE CFACTOR(ND,A,B,C,F,X,II,JJ)
      IMPLICIT NONE
C---------------------------------------------------------------------
C                   "Cyclicheskaya progonka"
C     MDIM >= MAX(NX,NY)*3-2        : maximul dimension in model
C     A, B, C :    diagonals of matrix   |                           |
C     F       :    right part            | Bii    Cii     0   ... Aii|
C     X       :    solution              | Aii+1  Bii+1   Cii+1 ...0 |
C     II      :    begin of the vectors  |  0     Aii+2   B  \     . |
C     JJ      :    end of the vectors    |  0     0     \   \  \   . |
C                                        |  .             \   \  \   |
C                                        |  0               \   \ C  |
C                                        |  Cjj   0 ...       Ajj Bjj|
C
C A(ii)*X(jj) + B(ii)*X(ii) + C(ii)*X(ii+1) = F(ii)
C A(ii)*X(ii-1) + B(ii)*X(ii) + C(ii)*X(ii+1) = F(ii) , ii=ii+1,jj-1
C A(jj)*X(jj-1) + B(jj)*X(jj) + C(jj)*X(ii) = F(jj)
C
C---------------------------------------------------------------------
      INCLUDE '0COM.INC'
      REAL     SCR1,  P(MDIM),  Q(MDIM)
      REAL     ALPHA(MDIM),  BETA(MDIM),  GAMMA(MDIM)
      INTEGER  II, JJ, I, ND
      REAL     A(ND), B(ND), C(ND), F(ND), X(ND)

C---------------------------------------------------------------------
      ALPHA(II+1) = -C(II)/B(II)
      BETA (II+1) =  F(II)/B(II)
      GAMMA(II+1) = -A(II)/B(II)

      DO 10 I=II+1,JJ
         SCR1 = -B(I) - A(I)*ALPHA(I)
         ALPHA(I+1) = C(I) / SCR1
         BETA (I+1) = (-F(I) + A(I)*BETA(I)) / SCR1
         GAMMA(I+1) = A(I)*GAMMA(I) / SCR1
  10  CONTINUE

      P(JJ-1) = BETA(JJ)
      Q(JJ-1) = ALPHA(JJ) + GAMMA(JJ)
      DO 20 I=JJ-2,II,-1
         P(I) = ALPHA(I+1)*P(I+1) + BETA (I+1)
         Q(I) = ALPHA(I+1)*Q(I+1) + GAMMA(I+1)
  20  CONTINUE

      X(JJ) = (BETA(JJ+1) + ALPHA(JJ+1)*P(II)) /
     &        (1 - ALPHA(JJ+1)*Q(II) - GAMMA(JJ+1))
      DO 30 I=II,JJ-1
         X(I) = P(I) + X(JJ)*Q(I)
  30  CONTINUE

      RETURN
      END
C=====================================================================
      SUBROUTINE FACTORS(ND,A,B,C,F,X,NBLOCS,IIB,JJB)
      IMPLICIT NONE
      INCLUDE '0COM.INC'
      INTEGER  IIB(MLR),   JJB(MLR)
      REAL     KSI(MDIM),  ETA(MDIM), SCR
      INTEGER  II, JJ, I, I1, ILOOP, NB, ND, NBLOCS
      REAL     A(ND), B(ND), C(ND), F(ND), X(ND)
      INTEGER  MSRV
C---------------------------------------------------------------------
      MSRV = ND
      IF (ND.EQ.NX) MSRV = MM
      DO 7 NB=1,NBLOCS
	  II = IIB(NB)
	  JJ = JJB(NB)
         IF (JJ.GT.MSRV) JJ = JJ - MMD
	  KSI(II) = -C(II)/B(II)
	  ETA(II) =  F(II)/B(II)
         DO 10 ILOOP=IIB(NB)+1,JJB(NB),1
             I  = ILOOP
             I1 = ILOOP - 1
             IF ( I.GT.MSRV) I  = I  - MMD
             IF (I1.GT.MSRV) I1 = I1 - MMD
             SCR = B(I) + A(I)*KSI(I1)
             KSI(I) = -C(I) / SCR
             ETA(I) = (F(I) - A(I)*ETA(I1)) / SCR
  10     CONTINUE
	  X(JJ) = ETA(JJ)
         DO 20 ILOOP=JJB(NB)-1,IIB(NB),-1
            I  = ILOOP
            I1 = ILOOP + 1
            IF ( I.GT.MSRV) I  = I  - MMD
            IF (I1.GT.MSRV) I1 = I1 - MMD
            X(I) = KSI(I)*X(I1) + ETA(I)
  20     CONTINUE
  7   CONTINUE

      RETURN
      END
C=====================================================================
      SUBROUTINE FACTORP(ND,A,B,C,F,X)
      IMPLICIT NONE
      INCLUDE '0COM.INC'
      REAL     KSI(MDIM),  ETA(MDIM), SCR
      INTEGER  ND, I 
      REAL     A(ND), B(ND), C(ND), F(ND), X(ND)
C---------------------------------------------------------------------
      KSI(1) = -C(1)/B(1)
      ETA(1) =  F(1)/B(1)
      DO 10 I=2,ND,1
         SCR = B(I) + A(I)*KSI(I-1)
         KSI(I) = -C(I) / SCR
         ETA(I) = (F(I) - A(I)*ETA(I-1)) / SCR
  10  CONTINUE
      X(ND) = ETA(ND)
      DO 20 I=ND-1,1,-1
         X(I) = KSI(I)*X(I+1) + ETA(I)
  20  CONTINUE

      RETURN
      END
C======================================================================
