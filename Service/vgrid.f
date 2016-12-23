C======================================================================
C (C) DIANSKY N.A., 2003
C-----------------------------------------------------09-26-97 11:18am
      SUBROUTINE VGRID()
	IMPLICIT NONE
C FOR SETTING VERTICAL T-,W- GRID LEVELS
C ZW(1) IS SURFACE, ZW(NZ) IS BOTTOM W-LEVELS.
C VERSION WITH CALIBRATION ON LEVITUS LEVELS.
C HH0 IS REFERENCE DEPTH OF OCEAN.
C MEAN WO DEPTH = 3760m(WITH    COS(LAT))
C MEAN WO DEPTH = 3490m(WITHOUT COS(LAT))
C IN 0COM.INC: Z  (NZ)   - T-SIGMA LEVELS
C              HZT(NZ)   - T-SIGMA LEVEL STEPS
C              ZW (NZ+1) - W-SIGMA LEVELS
C              DZ (NZ)   - W-SIGMA LEVEL STEPS

      INCLUDE '0COM.INC'
      INCLUDE '0VGRID.INC'
      REAL    HH0
      INTEGER NLEV
      PARAMETER ( HH0=3500.,  !AVERAGE DEPTH OF WORLD OCEAN
     &           NLEV=33)     !NUMBER OF LEVITUS LEVELS
      REAL DLEV(NLEV)         !LEVITUS HORIZONTS IN METERS FOR ANALITICAL SET
      DATA DLEV/0.,10.,20.,30.,50.0,75.0,100.0,125.0,150.0,200.0,250.0,
     &   300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,1100.0,1200.0,
     &  1300.0,1400.0,1500.0,1750.0,2000.0,2500.0,3000.0,3500.0,4000.0,
     &  4500.0,5000.0,5500.0/
      INTEGER K, NREFL
      REAL    BOTTOM, DEVIH, A, B, UNIDEPTH

C FINDING LEVITUS NUMBER OF REFERENCE DEPTH
      NREFL=1
      DEVIH=DLEV(33)

C     READ(*,*) IFLAGLEV
C     WRITE(*,*) '  IFLAGLEV =', IFLAGLEV

C ANALITICAL SET OF VERTICAL GRID
      DO K=2,33
         IF (ABS(DLEV(K)-HH0).LE.DEVIH) THEN
              DEVIH=ABS(DLEV(K)-HH0)
              NREFL=K
         ENDIF
      ENDDO
      WRITE(*,'(A,I3,A,F8.2)')
     &  ' NUMBER OF LEVITUS HORIZONTS:',NREFL,' FOR H=',HH0

      A=10.0*FLOAT(NREFL-1)/DLEV(NREFL)   ! GRADIENT IN UPPER OCEAN
      B=EXP(1.0)
C MAY USE B1 FOR MORE SLIGHTLY LEVELS IN UPPER OCEAN
C     B1=SQRT((FLOAT(NREFL)/FLOAT(NZ)))

       ZW(1)=0.0                  !SEA SURFACE
       ZW(NZ+1)=1.0               !BOTTOM

      IF (WGR_IN_TGR) THEN
C W-LEVELS ARE ARRANGED IN THE MIDDLES OF T-LAYERS

        IF(ANALITCAL_SET) THEN
C ANALITICAL T-LEVELS SETTING
         DO K=2,NZ-1
          Z(K)=UNIDEPTH((FLOAT(K)-0.5)/(FLOAT(NZ)-0.5),A,B)
         ENDDO
        ELSE
C NON-ANALITICAL T-LEVELS SETTING
C PROVE THE LEVELS
         DO K = 2,NZ
          IF(Z(K).LE.Z(K-1)) THEN
          WRITE(*,'(A,I4,F10.5)')
     &    '  ERROR IN SETING Z-LEVELS IN 0VGRID.INC. HORIZONT #',K,Z(K)
          STOP 1
          END IF
         ENDDO
C CORRECT THE LEVELS FOR 1-DEPTH
         BOTTOM=Z(NZ)+(Z(NZ)-Z(NZ-1))/2.0
         DO K=1,NZ
          Z(K)=Z(K)/BOTTOM
         ENDDO

        END IF
C REGULATING TOP AND BOTTOM T-LEVELS
          Z( 1) =        Z(2)   /3.0
          Z(NZ) =2.0/3.0+Z(NZ-1)/3.0
C W-LEVELS SETTING IN THE MIDDLES OF T-LAYERS
         DO K =2,NZ
           ZW(K  )= (Z (K)  + Z (K-1))/2.0
         ENDDO

        ELSE

C T-LEVELS ARE ARRANGED IN THE MIDDLES OF W-LAYERS
        IF(ANALITCAL_SET) THEN
C ANALITICAL W-LEVEL SETTING
         DO K=3,NZ
         ZW(K)=UNIDEPTH((FLOAT(K-1))/(FLOAT(NZ)-0.5),A,B)
         ENDDO

        ELSE

C NON-ANALITICAL W-LEVELS SETTING
C PROVE THE LEVELS
         DO K = 2,NZ
          IF(Z(K).LE.Z(K-1)) THEN
          WRITE(*,'(A,I4,F10.5)')
     &    '  ERROR IN SETING Z-LEVELS IN 0VGRID.INC. HORIZONT #',K,Z(K)
          STOP 1
          END IF
         ENDDO
C CORRECT THE LEVELS FOR 1-DEPTH
         BOTTOM=Z(NZ)+(Z(NZ)-Z(NZ-1))
         DO K=2,NZ
          ZW(K)=Z(K)/BOTTOM
         ENDDO

        ENDIF
C REGULATING TOP AND BOTTOM LEVELS
          ZW( 2) =  ZW(3)/2.0
          ZW(NZ) = (ZW(NZ-1)+ZW(NZ+1))/2.0

C T-LEVEL SETTING IN THE MIDDLE OF W-LAYER
         DO K =1,NZ
          Z(K )  = (ZW(K+1) + ZW(K))/2.0
         ENDDO
      ENDIF

C T AND W -GRID STEPS:
       HZT(1) = Z (1)
        DZ(1) = ZW(2)
      DO  K=2,NZ
       HZT(K) = Z (K)   - Z (K-1)
        DZ(K) = ZW(K+1) - ZW(K)
      ENDDO

C     BOTTOM=HH0
      BOTTOM=1000.0

      IF (WGR_IN_TGR) THEN
      WRITE(*,*)'  W-LEVELS ARE ARRANGED IN THE MIDDLES OF T-LAYERS.'
      ELSE
      WRITE(*,*)'  T-LEVELS ARE ARRANGED IN THE MIDDLES OF W-LAYERS.'
      END IF

      WRITE(*,110) BOTTOM
  110 FORMAT('  W-LEVELS W-STEPS  T-LEVELS T-STEPS *',F7.2)
      DO K=1,NZ
      WRITE(*,111)ZW(K+1)*BOTTOM,DZ(K)*BOTTOM,Z(K)*BOTTOM,HZT(K)*BOTTOM
      ENDDO
  111 FORMAT(2(2X,2F8.2))
      RETURN
      END
C======================================================================
      REAL FUNCTION UNIDEPTH(X,A,B)
	IMPLICIT NONE
C  UNIVERSAL DIMENSIONLESS FUNCTION OF UNEVEN OCEANOGRAPHYC HORIZONTS
C  CONSTRUCTED ON LEVITUS OCEANOGRAPHYC HORIZONTS
C  X-DIMENSIONLESS LEVEL VALUE FROM [0,1]
      REAL  X, A, B
      UNIDEPTH=(2.0-A)**(X**B)+A*X-1.0
      RETURN
      END
C======================================================================
