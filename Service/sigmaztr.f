C======================================================================
      SUBROUTINE Z2S(F1,F2,HHH,MSK,ZSIGMA,ZLVL,
     &                    NX,NY,NZ,NLVL,LEV1,OVER,IERR)
C ü DIANSKY N.A. 18.03.98 15:43
C PROGRAM FOR INTERPOLATION FROM COMMON LEVELS ON SIGMA LEVELS
C F1 - FIRST  3-D FIELD (INPUT DATA ON COMMON LEVELS)
C F2 - SECOND 3-D FIELD (OUTPUT DATA ON SIGMA LEVELS)
C HHH - 2-D FIELD OF BOTTOM TOPOGRAPHY (in centimeters!!!)
C MSK - MASK OF OCEAN GRIDS
C ZSIGMA - VALUES OF SIGMA LEVELS
C ZLVL - VALUES OF COMMON LEVELS (in meters!!!)
C NX,NY - X,Y DIMENSION
C NZ - NUMBER OF SIGMA LEVELS
C NLVL - NUMBER OF COMMON LEVELS
C LEV1 - PARAMETER OF TASK:
C IF LEV1=1 THEN FIRST LEVEL OF F2 IS DEFINITED AS FIST LEVEL OF F1
C IF LEV1=0 THEN FIRST LEVEL OF F2 IS DEFINITED CORRESPONDLY TO ITS
C           LEVEL POSITION
C OVER - UNDEFINITE VALUE
C IERR - ERROR INDICATOR (IF IERR=0 - THERE IS NO ERROR)
      INTEGER NX, NY, NZ, NLVL
      REAL F2(NX,NY,NZ),F1(NX,NY,NLVL),HHH(NX,NY),MSK(NX,NY),
     &     ZSIGMA(NZ),ZLVL(NLVL)
      REAL FZ1,FZ2,Z1,Z2,ZS(150),FZ(150),DS,DEEP
      REAL*4  OVER, OVER5
      INTEGER I,J,K,KDEEP,KUPS,NOCEGRID,LEV1,IERR, KR

C  OVER5 - 5% VICINITY TO OVER
      OVER5=ABS(OVER)-ABS(OVER)/20.0
      IERR=0
      WRITE(*,*)
     &' INTERPOLATION FROM COMMON LEVELS TO SIGMA LEVELS'
      IF (NZ.GT.150) THEN
         WRITE(*,*)' ERROR IN ROUTINE Z2S:'
         WRITE(*,*)' NUMBER OF SIGMA LEVELS IS GREATER THEN 150'
         IERR=1
         RETURN
      END IF
      IF (NLVL.GT.150) THEN
         WRITE(*,*)' ERROR IN ROUTINE Z2S:'
         WRITE(*,*)' NUMBER OF COMMON LEVELS IS GREATER THEN 150'
         IERR=1
         RETURN
      END IF
C  INTERPOLATION
      NOCEGRID=0
      DO J=1,NY
         DO I=1,NX
           IF(MSK(I,J).GT.0.5) THEN
           NOCEGRID=NOCEGRID+1
           DEEP=HHH(I,J)*1.0E-02                  !CM-->M
           DO K=1,NZ
              ZS(K)=ZSIGMA(K)*DEEP
           END DO

              FZ(1)=F1(I,J,1)
           IF(ABS(FZ(1)).GT.OVER5) THEN
              WRITE(*,*)' ERROR IN ROUTINE Z2S:'
              WRITE(*,'(A,2I4,A)')
     &  ' IN POINT ',I,J,' INPUT VALUE OF UPPER LEVEL IS UNDEFINITE!'
              IERR=1
              RETURN
           END IF

C  MAKING PROFILE WITHOUT BOTTOM
           IERR=0
           DO K=2,NLVL
              FZ(K)=F1(I,J,K)
              IF(ABS(FZ(K)).GT.OVER5) THEN
                FZ(K)=FZ(K-1)
                IERR=IERR+1
              END IF
           END DO

           IF(IERR.EQ.NLVL-1) THEN
              WRITE(*,*)' WARNING IN ROUTINE Z2S:'
              WRITE(*,1000) I,J,DEEP,ZLVL(1)
           END IF
           IERR=0

C  FINDING NUMBER OF UPPER SIGMA LEVEL
           KUPS=1
              DO WHILE(ZS(KUPS).LT.ZLVL(1).AND.KUPS.LT.NZ)
               KUPS=KUPS+1
              END DO
           KUPS=KUPS-1
           IF(KUPS.GE.1) THEN
C   FOR UPPER SIGMA LEVELS OF F2
              DO K=1,KUPS
                 F2(I,J,K)=F1(I,J,1)
              ENDDO
           ELSE
              IF (LEV1.NE.0) THEN
                 F2(I,J,1)=F1(I,J,1)
                 KUPS=1
              END IF
           END IF


C  FINDING DEEPEST Z LEVEL
           KDEEP=1
              DO WHILE(ZLVL(KDEEP).LT.DEEP.AND.KDEEP.LT.NLVL)
                 KDEEP=KDEEP+1
              END DO

           KR=1
           FZ1=FZ(1)
           FZ2=FZ1
           Z1=ZLVL(1)-ZLVL(1)/20.0
           Z2=ZLVL(1)

           DO K=KUPS+1,NZ
              DS=ZS(K)
              DO WHILE((DS.LE.Z1.OR.DS.GT.Z2).AND.KR.LT.KDEEP)
                 KR=KR+1
                 FZ1=FZ2
                 Z1=Z2
                 FZ2=FZ(KR)
                 Z2=ZLVL(KR)
              END DO

              IF(DS.GE.ZLVL(KDEEP)) THEN
                F2(I,J,K)=FZ(KDEEP)
              ELSE
                F2(I,J,K)=(FZ1*(Z2-DS)+FZ2*(DS-Z1))/(Z2-Z1)
              END IF
           END DO

           END IF
         ENDDO
      ENDDO
      WRITE(*,*)
     &' FOR CONTROL: NUMBER OF OCEAN GRID IN INTERPOLATION IS'
      WRITE(*,*)  NOCEGRID
      RETURN
1000  FORMAT(' IN POINT I=',I5,',',' J=',I5,' DEEP=',F9.2,'M'/
     &       ' THER IS ONLY ONE LEVEL FOR INTERPOLATION'/
     &       ' ON COMMON LEVEL OF',F9.2,'M')
      END
C======================================================================
      SUBROUTINE S2Z(F1,F2,HHH,MSK,ZSIGMA,ZLVL,
     &                    NX,NY,NZ,NLVL,LEV1,OVER,IERR)
C ü DIANSKY N.A. 18.03.98 15:43
C PROGRAM FOR INTERPOLATION FROM SIGMA LEVELS ON COMMON LEVELS
C F1 - FIRST  3-D FIELD (INPUT DATA ON SIGMA LEVELS)
C F2 - SECOND 3-D FIELD (OUTPUT DATA ON COMMON LEVELS)
C HHH - 2-D FIELD OF BOTTOM TOPOGRAPHY (in centimeters!!!)
C MSK - MASK OF OCEAN GRIDS
C ZSIGMA - VALUES OF SIGMA LEVELS
C ZLVL - VALUES OF COMMON LEVELS (in meters!!!)
C NX,NY - X,Y DIMENSION
C NZ - NUMBER OF SIGMA LEVELS
C NLVL - NUMBER OF COMMON LEVELS
C LEV1 - PARAMETER OF TASK:
C IF LEV1=1 THEN FIRST LEVEL OF F2 IS DEFINITED AS FIST LEVEL OF F1
C IF LEV1=0 THEN FIRST LEVEL OF F2 IS DEFINITED CORRESPONDLY TO ITS
C           LEVEL POSITION
C OVER - UNDEFINITE VALUE
C IERR - ERROR INDICATOR (IF IERR=0 - THERE IS NO ERROR)
      INTEGER NX, NY, NZ, NLVL, LEV1
      INTEGER K, I, J, KR, IERR, KDEEP, KSHAL, NOCEGRID
C
      REAL F1(NX,NY,NZ),F2(NX,NY,NLVL),HHH(NX,NY),MSK(NX,NY),
     &     ZSIGMA(NZ),ZLVL(NLVL)
      REAL FZ1,FZ2,ZS1,ZS2,ZS(150), DLVL, ZDEEP, DEEP
      REAL*4 OVER

      IERR=0
C      WRITE(*,*)
C     &' INTERPOLATION FROM SIGMA LEVELS ON COMMON LEVELS'
      IF (NZ.GT.150) THEN
         WRITE(*,*)' ERROR IN ROUTINE S2Z:'
         WRITE(*,*)' NUMBER OF SIGMA LEVELS IS GREATER THEN 150'
         IERR=1
         RETURN
      END IF
C INTERPOLATION
      NOCEGRID=0
      DO J=1,NY
         DO I=1,NX
           IF(MSK(I,J).GT.0.5) THEN
           NOCEGRID=NOCEGRID+1
           DEEP=HHH(I,J)*1.0E-02                  !CM-->M
           DO K=1,NZ
              ZS(K)=ZSIGMA(K)*DEEP
           END DO
C  FINDING LEVELS NEAREST TO SHALOW AND DEEP SIGMA LEVELS
           ZDEEP=1.0E+10
           KSHAL=0
           DO K=1,NLVL
                 DLVL=ZLVL(K)
              IF(DLVL.LE.ZS(1)) THEN
                 KSHAL=K
              END IF
              IF(ABS(ZS(NZ)-DLVL).LE.ZDEEP) THEN
                 KDEEP=K
                 ZDEEP=ABS(ZS(NZ)-DLVL)
              END IF
           END DO

           IF(KDEEP.EQ.1) THEN
c              WRITE(*,*)' WARNING IN ROUTINE S2Z:'
c              WRITE(*,1000) I,J,DEEP,ZLVL(1)
           END IF

           IF(KSHAL.GE.1) THEN
C   FOR UPPER LAYERS OF F2
              DO K=1,KSHAL
                 F2(I,J,K)=F1(I,J,1)
              ENDDO
           ELSE
              IF (LEV1.NE.0) THEN
                 F2(I,J,1)=F1(I,J,1)
                 KSHAL=1
              END IF
           END IF

           KR=1
           FZ1=F1(I,J,1)
           FZ2=FZ1
           ZS1=ZS(1)-ZS(1)/20.0
           ZS2=ZS(1)

           DO K=KSHAL+1,KDEEP
              DLVL=ZLVL(K)
              DO WHILE((DLVL.LE.ZS1.OR.DLVL.GT.ZS2).AND.KR.LT.NZ)
                 KR=KR+1
                 FZ1=FZ2
                 ZS1=ZS2
                 FZ2=F1(I,J,KR)
                 ZS2=ZS(KR)
              END DO

              IF(DLVL.GE.ZS(NZ)) THEN
                F2(I,J,K)=F1(I,J,NZ)
              ELSE
                F2(I,J,K)=(FZ1*(ZS2-DLVL)+FZ2*(DLVL-ZS1))/(ZS2-ZS1)
              END IF
           END DO

           IF(KDEEP.LT.NLVL) THEN
              DO K=KDEEP+1,NLVL
              F2(I,J,K)=OVER
              ENDDO
           END IF

           END IF
         ENDDO
      ENDDO
      WRITE(*,*)' SIGMA->Z INTERPOLATION: NUMBER OF OCEAN GRID IS',
     &            NOCEGRID
      RETURN
1000  FORMAT(' IN POINT I=',I5,',',' J=',I5,' DEEP=',F9.2,'M'/
     &       ' THER IS ONLY ONE LEVEL FOR INTERPOLATION'/
     &       ' ON COMMON LEVEL OF',F9.2,'M')
      END
C======================================================================
