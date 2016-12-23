C======================================================================
      SUBROUTINE UBARK8(UU,UTR,NOP,DZ,LU,NX,NY,NZ)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C MAKES FULL VELOCITIES OR BAROCLINIC VELOCITIES AND BAROTROPIC
C COMPONENT. TRANSPORT IS NOT STRICKTLY DIVERGENT (WWINT).
C UU = UU' + UTR ; UTR= <UU>, UU' = UU - <UU>.
      INTEGER NX, NY, NZ
      REAL    UU(NX,NY,NZ),LU(NX,NY),DZ(NZ)
      REAL*8 UTR(NX,NY),FF0
      INTEGER NOP
      INTEGER M, N, K

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,FF0)
      DO N=1,NY
       DO M=1,NX
        IF (LU(M,N).GT.0.5) THEN
         IF (NOP.EQ.+1) THEN
           DO K=1,NZ
             UU(M,N,K) = UU(M,N,K) + SNGL(UTR(M,N))
           ENDDO
C BAROCLINIC VELOCITIES AND BAROTROPIC COMPONENT.
         ELSEIF (NOP.EQ.-1) THEN
           FF0 = 0.
           DO K=1,NZ
             FF0 = FF0 + DBLE(UU(M,N,K) * DZ(K))
           ENDDO
           DO K=1,NZ
             UU(M,N,K) = UU(M,N,K) - SNGL(FF0)
           ENDDO
           UTR(M,N) = FF0
         ENDIF
        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
C======================================================================
      SUBROUTINE WWINTC(UU,VV,WW)
      IMPLICIT NONE
C----------------------------------------------------------------------
C COMPUTING WW FROM CONTINUITY EQUATION ON GRID "C".
      INCLUDE '0COM.INC'
C INPUT:
C UU - ZONAL VELOCITY ON UC-GRID
C VV - MERIDIONAL VELOCITY ON VC-GRID
C OUTPUT:
C WW - VERTICAL VELOCITY ON T-GRID
      REAL UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ+1)
      INTEGER M, N, K

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K)
      DO N=NNN,NN
       DO M=MMM,MM
        IF (LU(M,N).GT.0.5) THEN
         WW(M,N,NZ+1)=0.0
         DO K=NZ,1,-1
          WW(M,N,K) = WW(M,N,K+1) + DZ(K)/(DX(M,N)*DY(M,N)*RN)*

     & ( ( UU(M  ,N,K)*HHU(M  ,N)*DYH(M  ,N) -
     &     UU(M-1,N,K)*HHU(M-1,N)*DYH(M-1,N) ) 
     &     + 

     &   ( VV(M,N  ,K)*HHV(M,N  )*DXH(M,N  ) -
     &     VV(M,N-1,K)*HHV(M,N-1)*DXH(M,N-1) ) )
         ENDDO
        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
C======================================================================
      SUBROUTINE UCGR2TGR(UIN,UOUT,KZ)
      IMPLICIT NONE
      INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C INTERPOLATION: UIN ON UC-GRID -> UOUT ON T-GRID
C LU  - T -GRID OCEAN MASK OF AINP (1.- OCEAN, 0.-LAND)
C LCU - UC-GRID OCEAN MASK OF AOUT (1.- OCEAN, 0.-LAND)
C NX,NY - X,Y- UNIVERSAL DIMENSION; NZ - VARIABLE DIMENSION ON Z.
      INTEGER  KZ
      REAL UIN(NX,NY,KZ),UOUT(NX,NY,KZ)
      INTEGER M, N, K

C 2-POINT INTERPOLATING UIN DEFINED ON LCU-GRID TO UOUT DEFINED ON LU-GRID.
              
              DO  K=1,KZ
!$OMP PARALLEL DO PRIVATE(M,N)
            DO  N=NNN,NN
          DO  M=MMM,MM
         UOUT(M,N,K)= LU(M,N)*0.5/DY(M,N)/(HHQ(M,N)+1.0E-07)
     &   *(UIN(M-1,N,K)*HHU(M-1,N)*DYH(M-1,N) +
     &     UIN(M  ,N,K)*HHU(M  ,N)*DYH(M  ,N)  )
          ENDDO
            ENDDO
!$OMP END PARALLEL DO
              ENDDO

      RETURN
      END
C======================================================================
      SUBROUTINE UCGR2TGR_SUBREG(UIN,UOUT,KX1,KX2,KY1,KY2,KZ)
      IMPLICIT NONE
      INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C INTERPOLATION: UIN ON UC-GRID -> UOUT ON T-GRID
C LU  - T -GRID OCEAN MASK OF AINP (1.- OCEAN, 0.-LAND)
C LCU - UC-GRID OCEAN MASK OF AOUT (1.- OCEAN, 0.-LAND)
C NX,NY - X,Y- UNIVERSAL DIMENSION; NZ - VARIABLE DIMENSION ON Z.
C KX1,KX2,KY1,KY2,KZ  SUBREGION OF INTERPOLATION

      INTEGER  KX1,KX2,KY1,KY2,KZ
      REAL UIN(NX,NY,KZ),UOUT(NX,NY,KZ)
      INTEGER M, N, K

C 2-POINT INTERPOLATING UIN DEFINED ON LCU-GRID TO UOUT DEFINED ON LU-GRID.
              
              DO  K=1,KZ
!$OMP PARALLEL DO PRIVATE(M,N)
            DO  N=KY1,KY2
          DO  M=KX1,KX2
         UOUT(M,N,K)= LU(M,N)*0.5/DY(M,N)/(HHQ(M,N)+1.0E-07)
     &   *(UIN(M-1,N,K)*HHU(M-1,N)*DYH(M-1,N) +
     &     UIN(M  ,N,K)*HHU(M  ,N)*DYH(M  ,N)  )
          ENDDO
            ENDDO
!$OMP END PARALLEL DO
              ENDDO

      RETURN
      END
C======================================================================
      SUBROUTINE VCGR2TGR(VIN,VOUT,KZ)
      IMPLICIT NONE
      INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C INTERPOLATION: VIN ON VC-GRID -> VOUT ON T-GRID
C LU  - T -GRID OCEAN MASK OF AINP (1.- OCEAN, 0.-LAND)
C LCV - VC-GRID OCEAN MASK OF AOUT (1.- OCEAN, 0.-LAND)
C NX,NY - X,Y- UNIVERSAL DIMENSION; NZ - VARIABLE DIMENSION ON Z.
      INTEGER  KZ
      REAL     VIN(NX,NY,KZ),VOUT(NX,NY,KZ)
      INTEGER  M,N,K

C 2-POINT INTERPOLATING UIN DEFINED ON LCV-GRID TO UOUT DEFINED ON LU-GRID.

              DO  K=1,KZ
!$OMP PARALLEL DO PRIVATE(M,N)
            DO  N=NNN,NN
          DO  M=MMM,MM
         VOUT(M,N,K)= LU(M,N)*0.5/DX(M,N)/(HHQ(M,N)+1.0E-07) 
     &   *( VIN(M,N-1,K)*HHV(M,N-1)*DXH(M,N-1) +
     &      VIN(M,N  ,K)*HHV(M,N  )*DXH(M,N  )   )
      
          ENDDO
            ENDDO
!$OMP END PARALLEL DO              
              ENDDO

      RETURN
      END
C======================================================================
      SUBROUTINE VCGR2TGR_SUBREG(VIN,VOUT,KX1,KX2,KY1,KY2,KZ)
      IMPLICIT NONE
      INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C INTERPOLATION: VIN ON VC-GRID -> VOUT ON T-GRID
C LU  - T -GRID OCEAN MASK OF AINP (1.- OCEAN, 0.-LAND)
C LCV - VC-GRID OCEAN MASK OF AOUT (1.- OCEAN, 0.-LAND)
C NX,NY - X,Y- UNIVERSAL DIMENSION; NZ - VARIABLE DIMENSION ON Z.
C KX1,KX2,KY1,KY2,KZ  SUBREGION OF INTERPOLATION

      INTEGER  KX1,KX2,KY1,KY2,KZ
      REAL     VIN(NX,NY,KZ),VOUT(NX,NY,KZ)
      INTEGER  M,N,K

C 2-POINT INTERPOLATING UIN DEFINED ON LCV-GRID TO UOUT DEFINED ON LU-GRID.

              DO  K=1,KZ
!$OMP PARALLEL DO PRIVATE(M,N)
            DO  N=KY1,KY2
          DO  M=KX1,KX2

         VOUT(M,N,K)= LU(M,N)*0.5/DX(M,N)/(HHQ(M,N)+1.0E-07) 
     &   *( VIN(M,N-1,K)*HHV(M,N-1)*DXH(M,N-1) +
     &      VIN(M,N  ,K)*HHV(M,N  )*DXH(M,N  )   )
      
          ENDDO
            ENDDO
!$OMP END PARALLEL DO              
              ENDDO

      RETURN
      END
C======================================================================
      SUBROUTINE UCGR2TGR8(UIN,UOUT,KZ)
      IMPLICIT NONE
      INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C INTERPOLATION: UIN ON UC-GRID -> UOUT ON T-GRID
C LU  - T -GRID OCEAN MASK OF AINP (1.- OCEAN, 0.-LAND)
C LCU - UC-GRID OCEAN MASK OF AOUT (1.- OCEAN, 0.-LAND)
C NX,NY - X,Y- UNIVERSAL DIMENSION; NZ - VARIABLE DIMENSION ON Z.
      INTEGER  KZ
      REAL*8 UIN(NX,NY,KZ)
      REAL  UOUT(NX,NY,KZ)
      INTEGER M, N, K

C 2-POINT INTERPOLATING UIN DEFINED ON LCU-GRID TO UOUT DEFINED ON LU-GRID.

              DO  K=1,KZ
!$OMP PARALLEL DO PRIVATE(M,N)
            DO  N=NNN,NN
          DO  M=MMM,MM
         UOUT(M,N,K)= LU(M,N)*0.5/DY(M,N)/(HHQ(M,N)+1.0E-07)
     &   *(UIN(M-1,N,K)*HHU(M-1,N)*DYH(M-1,N) +
     &     UIN(M  ,N,K)*HHU(M  ,N)*DYH(M  ,N)  )
          ENDDO
            ENDDO
!$OMP END PARALLEL DO
              ENDDO

      RETURN
      END
C======================================================================
      SUBROUTINE VCGR2TGR8(VIN,VOUT,KZ)
      IMPLICIT NONE
      INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C INTERPOLATION: VIN ON VC-GRID -> VOUT ON T-GRID
C LU  - T -GRID OCEAN MASK OF AINP (1.- OCEAN, 0.-LAND)
C LCV - VC-GRID OCEAN MASK OF AOUT (1.- OCEAN, 0.-LAND)
C NX,NY - X,Y- UNIVERSAL DIMENSION; NZ - VARIABLE DIMENSION ON Z.
      INTEGER  KZ
      REAL*8 VIN(NX,NY,KZ)
      REAL  VOUT(NX,NY,KZ)
      INTEGER  M,N,K

C 2-POINT INTERPOLATING UIN DEFINED ON LCV-GRID TO UOUT DEFINED ON LU-GRID.
              DO  K=1,KZ
!$OMP PARALLEL DO PRIVATE(M,N)
            DO  N=NNN,NN
          DO  M=MMM,MM
         VOUT(M,N,K)= LU(M,N)*0.5/DX(M,N)/(HHQ(M,N)+1.0E-07) 
     &   *( VIN(M,N-1,K)*HHV(M,N-1)*DXH(M,N-1) +
     &      VIN(M,N  ,K)*HHV(M,N  )*DXH(M,N  )   )
          ENDDO
            ENDDO
!$OMP END PARALLEL DO
              ENDDO

      RETURN
      END
C======================================================================
      SUBROUTINE UCGR2TGR_2D(UIN,UOUT)
      IMPLICIT NONE
      INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C INTERPOLATION: UIN ON UC-GRID -> UOUT ON T-GRID
C LU  - T -GRID OCEAN MASK OF AINP (1.- OCEAN, 0.-LAND)
C LCU - UC-GRID OCEAN MASK OF AOUT (1.- OCEAN, 0.-LAND)
C NX,NY - X,Y- UNIVERSAL DIMENSION; NZ - VARIABLE DIMENSION ON Z.
      INTEGER  KZ
      REAL UIN(NX,NY),UOUT(NX,NY)
      INTEGER M, N

C 2-POINT INTERPOLATING UIN DEFINED ON LCU-GRID TO UOUT DEFINED ON LU-GRID.
              

!$OMP PARALLEL DO PRIVATE(M,N)
            DO  N=NNN,NN
          DO  M=MMM,MM
         UOUT(M,N)= LU(M,N)*0.5/DY(M,N)
     &   *(UIN(M-1,N)*DYH(M-1,N) +
     &     UIN(M  ,N)*DYH(M  ,N)  )
          ENDDO
            ENDDO
!$OMP END PARALLEL DO


      RETURN
      END
C======================================================================
      SUBROUTINE VCGR2TGR_2D(VIN,VOUT)
      IMPLICIT NONE
      INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C INTERPOLATION: VIN ON VC-GRID -> VOUT ON T-GRID
C LU  - T -GRID OCEAN MASK OF AINP (1.- OCEAN, 0.-LAND)
C LCV - VC-GRID OCEAN MASK OF AOUT (1.- OCEAN, 0.-LAND)
C NX,NY - X,Y- UNIVERSAL DIMENSION; NZ - VARIABLE DIMENSION ON Z.
      INTEGER  KZ
      REAL     VIN(NX,NY),VOUT(NX,NY)
      INTEGER  M,N

C 2-POINT INTERPOLATING UIN DEFINED ON LCV-GRID TO UOUT DEFINED ON LU-GRID.


!$OMP PARALLEL DO PRIVATE(M,N)
            DO  N=NNN,NN
          DO  M=MMM,MM
         VOUT(M,N)= LU(M,N)*0.5/DX(M,N)
     &   *( VIN(M,N-1)*DXH(M,N-1) +
     &      VIN(M,N  )*DXH(M,N  )   )
      
          ENDDO
            ENDDO
!$OMP END PARALLEL DO              


      RETURN
      END
