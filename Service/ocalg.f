C======================================================================
      SUBROUTINE FF1AR (FF1,CON1,CON2,FF2,LU,NX,NY,NZ)
	IMPLICIT NONE
C----------------------------------------------------------------------
C LINEAR OPERATION: FF2 = CON1 * FF1 + CON2
      INTEGER NX, NY, NZ
      INTEGER M, N, K
      REAL FF1(NX,NY,NZ),FF2(NX,NY,NZ),LU(NX,NY),CON1,CON2
      REAL       OVER
      PARAMETER (OVER=0.0E+00)

!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO N=1,NY
        DO M=1,NX
          IF (LU(M,N).GT.0.5) THEN
            DO K=1,NZ
              FF2(M,N,K) = FF1(M,N,K) * CON1 + CON2
            ENDDO
          ELSE
            DO K=1,NZ
              FF2(M,N,K) = OVER
            ENDDO
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
C======================================================20.02.98 15:28
      SUBROUTINE CYCLIZE (FF,NX,NY,NZ,MMM,MM)
	IMPLICIT NONE
C---------------------------------------------------------------------
C ADDS PERIODICALLY LEFT (M=MMM-1) AND RIGHT (M=MM+1)
C FOR CYCLIC LINES EITHER FOR TT OR UU-ARRAYS.
      INTEGER NX, NY, NZ
      INTEGER MMM, MM, N, K
      REAL FF(NX,NY,NZ)
      
      DO K=1,NZ
!$OMP PARALLEL DO PRIVATE(N)      
          DO N=1,NY
              FF(MMM-1,N,K) = FF(MM ,N,K)
              FF(MM +1,N,K) = FF(MMM,N,K)
          ENDDO
!$OMP END PARALLEL DO      
      ENDDO
      RETURN
      END
C======================================================20.02.98 15:28
      SUBROUTINE CYCLIZE8(FF,NX,NY,NZ,MMM,MM)
	IMPLICIT NONE
C---------------------------------------------------------------------
C ADDS PERIODICALLY LEFT (M=MMM-1) AND RIGHT (M=MM+1)
C FOR CYCLIC LINES EITHER FOR TT OR UU-ARRAYS.
      INTEGER NX, NY, NZ
      INTEGER MMM, MM, N, K
      REAL*8 FF(NX,NY,NZ)

       
      DO K=1,NZ
!$OMP PARALLEL DO PRIVATE(N)  
          DO N=1,NY
              FF(MMM-1,N,K) = FF(MM ,N,K)
              FF(MM +1,N,K) = FF(MMM,N,K)
          ENDDO
!$OMP END PARALLEL DO
      ENDDO
      
      RETURN
      END
C======================================================================
      SUBROUTINE HH2TGR(AIN,LUU,AOUT,LU,OVER,NX,NY,NZ,
     &                  MMM,MM,NNN,NN)
      IMPLICIT NONE
C----------------------------------------------------------------------
C INTERPOLATION: AIN(H-GRID) -> AOUT (T-GRID)
C LUU - OCEAN MASK OF AINP (1.- OCEAN, 0.-LAND)
C LU  - OCEAN MASK OF AOUT (1.- OCEAN, 0.-LAND)
C OVER - SETTING VALUE TO UNDEFINED POINT
C NX,NY - X,Y- UNIVERSAL DIMENSION; NZ - VARIABLE DIMENSION ON Z.


      INTEGER NX, NY, NZ
      INTEGER M, N, K
      REAL AIN(NX,NY,NZ),LUU(NX,NY),AOUT(NX,NY,NZ),LU(NX,NY)
      INTEGER MMM,MM,NNN,NN
      REAL S4,D4,OVER


C 4 POINT INTERPOLATING AINP DEFINED ON LUU-GRID TO AOUT DEFINED ON LU-GRID.

       DO  N=NNN,NN
        DO  M=MMM,MM

         IF (LU(M,N).GT.0.5) THEN
          
          DO  K=1,NZ

           S4=AIN(M-1,N-1,K)*LUU(M-1,N-1)
     &       +AIN(M  ,N-1,K)*LUU(M  ,N-1)
     &       +AIN(M-1,N  ,K)*LUU(M-1,N  )
     &       +AIN(M  ,N  ,K)*LUU(M  ,N  )

           D4=LUU(M-1,N-1)
     &       +LUU(M  ,N-1)
     &       +LUU(M-1,N  )
     &       +LUU(M  ,N  )

           IF (D4.GT.0.5) THEN
            AOUT(M,N,K)=S4/D4
           ELSE
            AOUT(M,N,K)=OVER
           ENDIF
        
          ENDDO
         ENDIF

        ENDDO
       ENDDO



      RETURN
      END

C     FUNCTION OF SENSIBLE HEAT FLUX
      REAL FUNCTION SENS_HEAT(RHO,CP,CD,WND,TATM,TEMP)
      IMPLICIT NONE
      REAL RHO,            !DENSITY OF AIR
     &      CP,            !HEAT CAPACITY OF AIR
     &      CD,            !COEFFICIENT OF EXCHANGE
     &     WND,            !WIND SPEED
     &    TATM,            !TEMPERATURE OF AIR
     &    TEMP             !TEMPERATURE OF SURFACE
      
      SENS_HEAT=
     &    RHO*CP*CD*WND*(TATM-TEMP)

	END FUNCTION

C     FUNCTION OF LATENT HEAT FLUX
      REAL FUNCTION LAT_HEAT(RHO,QL,CD,WND,QATM,QSAT)
      IMPLICIT NONE
      REAL RHO,            !DENSITY OF AIR
     &      QL,            !HEAT OF EVAPORATION/SUBLIMATION
     &      CD,            !COEFFICIENT OF EXCHANGE
     &     WND,            !WIND SPEED
     &    QATM,            !HUMIDITY OF AIR
     &    QSAT             !HUMIDITY OF SATURATED VAPOR
      
      LAT_HEAT=
     &    RHO*QL*CD*WND*(QATM-QSAT)

	END FUNCTION


C     FUNCTION OF SATURATED VAPOR HUMIDITY
	REAL FUNCTION SATUR_HUM(P_SAT,P_A)
	IMPLICIT NONE
      REAL P_SAT,      !PRESSURE OF SATURATED VAPOR(Pa)
     &     P_A         !PRESSURE OF AIR(Pa)

      SATUR_HUM=
     &       0.622*P_SAT/(P_A-0.378*P_SAT)
	END FUNCTION

	
C     FUNCTION OF SATURATED VAPOR PRESSURE OVER OCEAN
      REAL FUNCTION SATUR_PRESS_OCEAN(TEM)
	IMPLICIT NONE
      REAL TEM              !TEMPERATURE OF OCEAN SURFACE(DEGREES C)

      SATUR_PRESS_OCEAN=
     &   10.0**((0.7859+0.03477*TEM)/(1.0+0.00412*TEM)+2.0)
	END FUNCTION


C     FUNCTION OF SATURATED VAPOR PRESSURE OVER ICE/SNOW
      REAL FUNCTION SATUR_PRESS_ICE(TEM)
	IMPLICIT NONE
      REAL TEM              !TEMPERATURE OF ICE/SNOW SURFACE(DEGREES C)

      SATUR_PRESS_ICE=
     &   10.0**((0.7859+0.03477*TEM)/(1.0+0.00412*TEM)+0.00422*TEM+2.0)
	END FUNCTION

      REAL FUNCTION TF(S)
      IMPLICIT NONE
      REAL S
c        TF = -0.054*MAX(S,0.0)
         TF = -0.0575*MAX(S,0.0)
      END  FUNCTION

C======================================================================
      REAL FUNCTION RINTEGR(XX,YY,NN)
C-----------------------------------------------------------------------
C INTEGRATES FUNCTION YY GIVEN ON GRID XX.
      INTEGER NN
      REAL XX(NN),YY(NN)
      INTEGER N
      REAL    R

      R = 0.
      
!$OMP PARALLEL DO REDUCTION(+:R)
      DO N=1,NN-1
       R = R + (YY(N+1)+YY(N))/2.*(XX(N+1)-XX(N))
      ENDDO
!$OMP END PARALLEL DO 
      RINTEGR = R
      RETURN
      END
C======================================================================
      SUBROUTINE EXTU (HH,HHQ,LU,NTO,NX,NY,NZ)
C----------------------------------------------------------------------
C INTERPOLATION: HH -> HHQ. GRID TYPE NTO=1: T->U; NTO=-1: U->T.
C LU  - OCEAN MASK OF HHQ (1.- OCEAN, 0.-LAND)
C NX,NY - X,Y- UNIVERSAL DIMENSION; NZ - VARIABLE DIMENSION ON Z.
      INTEGER NX, NY, NZ
      INTEGER M, N, K
      REAL HH(NX,NY,NZ),HHQ(NX,NY,NZ),LU(NX,NY)
      INTEGER NTO

C 4 POINT INTERPOLATING HH TO HHQ GIVEN ON LU-GRID.
      IF (NTO.EQ.1 ) THEN
      DO  N=1,NY
       DO  M=1,NX
        IF (LU(M,N).GT.0.5) THEN
          DO  K=1,NZ
        HHQ(M,N,K)=(HH(M,N,K)+HH(M,N+1,K)+HH(M+1,N,K)+HH(M+1,N+1,K))/4.
          ENDDO
        ENDIF
       ENDDO
      ENDDO
      ENDIF

      IF (NTO.EQ.-1) THEN
      DO  N=1,NY
       DO  M=1,NX
        IF (LU(M,N).GT.0.5) THEN
          DO  K=1,NZ
        HHQ(M,N,K)=(HH(M,N,K)+HH(M,N-1,K)+HH(M-1,N,K)+HH(M-1,N-1,K))/4.
          ENDDO
        ENDIF
       ENDDO
      ENDDO
      ENDIF
      RETURN
      END
C======================================================================
      SUBROUTINE EXTUD(AIN,LUU,AOUT,LU,OVER,J,NX,NY,NZ)
C----------------------------------------------------------------------
C INTERPOLATION: AIN -> AOUT. GRID TYPE J=1: T->U; J=-1: U->T.
C LUU - OCEAN MASK OF AINP (1.- OCEAN, 0.-LAND)
C LU  - OCEAN MASK OF AOUT (1.- OCEAN, 0.-LAND)
C OVER - SETTING VALUE TO UNDEFINED POINT
C NX,NY - X,Y- UNIVERSAL DIMENSION; NZ - VARIABLE DIMENSION ON Z.
      INTEGER NX, NY, NZ
      INTEGER M, N, K
      REAL AIN(NX,NY,NZ),LUU(NX,NY),AOUT(NX,NY,NZ),LU(NX,NY)
      INTEGER J
      REAL S4,D4,OVER

C 4 POINT INTERPOLATING AINP DEFINED ON LUU-GRID TO AOUT DEFINED ON LU-GRID.
      DO  K=1,NZ
      DO  N=2,NY-1
       DO  M=2,NX-1
        IF (LU(M,N).GT.0.5) THEN
          S4=AIN(M  ,N,K)*LUU(M  ,N)+AIN(M  ,N+J,K)*LUU(M  ,N+J)+
     +       AIN(M+J,N,K)*LUU(M+J,N)+AIN(M+J,N+J,K)*LUU(M+J,N+J)
          D4=             LUU(M  ,N)               +LUU(M  ,N+J)+
     +                    LUU(M+J,N)               +LUU(M+J,N+J)
         IF (D4.GT.0.5) THEN
         AOUT(M,N,K)=S4/D4
         ELSE
         AOUT(M,N,K)=OVER
         ENDIF
        ENDIF
       ENDDO
      ENDDO
      ENDDO

      RETURN
      END

