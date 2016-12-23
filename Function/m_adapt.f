C======================================================================
      SUBROUTINE U_BY_PRESSURE_GRADIENT(UU,DEN,TAU)
C  CALCULATION ZONAL VELOCITY FORCED BY PRESSURE GRADIENT
C----------------------------------------------------------------------
      IMPLICIT NONE
C***************************************************************
C  Calculation of Pressure gradients by (Alexeev&Zalesny,1993) - formulas
C                 on grid "C":
C  den = DEN(x,y,s)  - density deviation from ro0
C               s                                          
C  RIX5 = [d/dx($[Z_s*den-Z*den_s]ds)-(den*dZ/dx-d den/dx*Z)] 
C               0 
C  Calculate velocity                                1
C  U(x,y,s)  = -m*(g/ro0)*1/2* RIX5
C                                    0
C                              1
C  DINX(x,y) = -m*(g/ro0)*1/2* $RIX5(s)ds
C                              0

C***************************************************************
C PRESSURE GRADIENTS BY (ALEXEEV&ZALESNY,1993) - FORMULA.

      INCLUDE '0COM.INC'
      REAL  TAU    !TIME STEP

      REAL UU(NX,NY,NZ),   !ZONAL VELOCITY
     &     DEN(NX,NY,NZ)   !WATER DENSITY


C     REAL(8) DINX(NX,NY),RIX5I  !ZONAL DEPTH-INTERGATED PRESSURE GRADSIENT

	INTEGER M, N, K
C HELP VARIABLES: NUMBERING IS IN ACCORDING KEYBOARD
      REAL(8) F6M5,RIX5(NZ)

C     DINX=0.0

C  COMPUTING IN X-DIRECTION

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,F6M5,RIX5)      
      DO N=NNN,NN
       DO M=MMM,MM

        IF (LCU(M,N).GT.0.5) THEN

C        H5 = DBLE(HHQ(M  ,N  ))
C        H6 = DBLE(HHQ(M+1,N  ))

C        F5 = DBLE(Z(1)*DEN(M  ,N  ,1)*H5)
C        F6 = DBLE(Z(1)*DEN(M+1,N  ,1)*H6)

C        F6M5 =DBLE( Z3D(M+1,N,1)*DEN(M+1,N,1)-Z3D(M,N,1)*DEN(M,N,1) )
         F6M5 =DBLE(
     &     ( (DEN(M+1,N,1)* Z3D(M+1,N,2)-DEN(M+1,N,2)*Z3D(M+1,N,1))
     &       *Z3D(M+1,N,1)/(Z3D(M+1,N,2)-Z3D(M+1,N,1)) )
     &    -( (DEN(M  ,N,1)* Z3D(M  ,N,2)-DEN(M  ,N,2)*Z3D(M  ,N,1))
     &       *Z3D(M  ,N,1)/(Z3D(M  ,N,2)-Z3D(M  ,N,1)) )
     &             )

C        RIX5(1)=(F6-F5)-Z(1)*(H6*DEN(M,N,1)-H5*DEN(M+1,N,1))
         RIX5(1)=F6M5 - DBLE(
     &           (Z3D(M+1,N,1)*DEN(M  ,N,1)-Z3D(M,N,1)*DEN(M+1,N,1)) )         
	   
	   DO K=2,NZ

C         F5 = F5 + DBLE( Z(K)*DEN(M  ,N,K-1) - Z(K-1)*DEN(M  ,N,K) )*H5
C         F6 = F6 + DBLE( Z(K)*DEN(M+1,N,K-1) - Z(K-1)*DEN(M+1,N,K) )*H6
          F6M5 = F6M5 + DBLE(
     &   ( Z3D(M+1,N,K)*DEN(M+1,N,K-1) - Z3D(M+1,N,K-1)*DEN(M+1,N,K) )   
     &  -( Z3D(M  ,N,K)*DEN(M  ,N,K-1) - Z3D(M  ,N,K-1)*DEN(M  ,N,K) ) ) 

C         RIX5(K) = (F6-F5)-DBLE( Z(K)*(H6*DEN(M,N,K)-H5*DEN(M+1,N,K)))
          RIX5(K)=  F6M5 - DBLE(
     &              (Z3D(M+1,N,K)*DEN(M,N,K)-Z3D(M,N,K)*DEN(M+1,N,K)) )
         ENDDO

C  IF NEED MAY BE UNCOMMENTED
C	   RIX5I = 0.0
C        DO K=1,NZ
C
C         RIX5I = RIX5I + RIX5(K) * DZ(K)
C
C        ENDDO
C
C COMPUTING DINX 
C        DINX(M,N) =     DBLE(-RM(N)*G/2.0/DXT(M)) *RIX5I

C DEFINE NEW ZONAL VELOCITY 
         DO K=1,NZ

          UU(M,N,K)=SNGL(DBLE(UU(M,N,K))
     &             -DBLE(G)/2.0d0/dble(DXT(M,N)*RN)*dble(TAU)*RIX5(K))

         ENDDO

        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF(MMD.NE.0) THEN
C         CALL CYCLIZE8(DINX,NX,NY, 1,MMM,MM)
          CALL CYCLIZE (UU  ,NX,NY,NZ,MMM,MM)
      END IF

      RETURN
      END
C======================================================================
      SUBROUTINE V_BY_PRESSURE_GRADIENT(VV,DEN,TAU)
C  CALCULATION MERIDIONAL VELOCITY FORCED BY PRESSURE GRADIENT
C----------------------------------------------------------------------
      IMPLICIT NONE
C***************************************************************
C  Calculation of Pressure gradients by (Alexeev&Zalesny,1993) - formulas
C                 on grid "C":
C  den = DEN(x,y,s)  - density deviation from ro0
C               s                                          
C  RIY5 = [d/dy($[Z_s*den-Z*den_s]ds)-(den*dZ/dy-d den/dyZ)] 
C               0 
C  Calculate velocity                               1
C  V(x,y,s)  = -n*(g/ro0)*1/2* RIY5
C                                    0
C                              1
C  DINY(x,y) = -n*(g/ro0)*1/2* $RIX5(s)ds
C                              0


C***************************************************************
C PRESSURE GRADIENTS BY (ALEXEEV&ZALESNY,1993) - FORMULA.

      INCLUDE '0COM.INC'
      REAL  TAU    !TIME STEP

      REAL VV (NX,NY,NZ),    !MERIDIONAL VELOCITY
     &     DEN(NX,NY,NZ)     !WATER DENSITY


C     REAL(8) DINY(NX,NY),RIY5I    !MERIDIONAL DEPTH-INTERGATED PRESSURE GRADSIENTS

	INTEGER M, N, K
C HELP VARIABLES: NUMBERING IS IN ACCORDING KEYBOARD
      REAL(8) F8M5,RIY5(NZ)

C     DINY=0.0

C  COMPUTING IN Y-DIRECTION

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,F8M5,RIY5)  
      DO N=NNN,NN
       DO M=MMM,MM

        IF (LCV(M,N).GT.0.5) THEN

C	   F8M5 = DBLE( Z3D(M,N+1,1)*DEN(M,N+1,1)-Z3D(M,N,1)*DEN(M,N,1) )

         F8M5 =DBLE(
     &     ( (DEN(M,N+1,1)* Z3D(M,N+1,2)-DEN(M,N+1,2)*Z3D(M,N+1,1))
     &       *Z3D(M,N+1,1)/(Z3D(M,N+1,2)-Z3D(M,N+1,1)) )
     &    -( (DEN(M,N  ,1)* Z3D(M,N  ,2)-DEN(M,N  ,2)*Z3D(M,N  ,1))
     &       *Z3D(M,N  ,1)/(Z3D(M,N  ,2)-Z3D(M,N  ,1)) )
     &             )

         RIY5(1)=  F8M5 - DBLE(
     &           (Z3D(M,N+1,1)*DEN(M,N,1)-Z3D(M,N,1)*DEN(M,N+1,1)) )  

         DO K=2,NZ

          F8M5 = F8M5 + DBLE(
     &   ( Z3D(M,N+1,K)*DEN(M,N+1,K-1) - Z3D(M,N+1,K-1)*DEN(M,N+1,K) )   
     &  -( Z3D(M,N  ,K)*DEN(M,N  ,K-1) - Z3D(M,N  ,K-1)*DEN(M,N  ,K) ) ) 

          RIY5(K)=  F8M5 - DBLE(
     &              (Z3D(M,N+1,K)*DEN(M,N,K)-Z3D(M,N,K)*DEN(M,N+1,K)) )

         ENDDO

C  IF NEED MAY BE UNCOMMENTED
C         RIY5I = 0.0
C        DO K=1,NZ
C
C         RIY5I = RIY5I + RIY5(K) * DBLE(DZ(K))
C
C        ENDDO
C
C COMPUTING DINY FOR ISLAND STREAM FUNCTION
C        DINY(M,N) =      DBLE(-RN *G/2.0/DYT(N))*RIY5I

         DO K=1,NZ

          VV(M,N,K)=SNGL(DBLE(VV(M,N,K))
     &             -DBLE(G)/2.0d0/dble(DYT(M,N)*RN)*dble(TAU)*RIY5(K))

         ENDDO

        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO
      IF(MMD.NE.0) THEN
C         CALL CYCLIZE8(DINY,NX,NY, 1,MMM,MM)
          CALL CYCLIZE (VV  ,NX,NY,NZ,MMM,MM)
      END IF

      RETURN
      END
C======================================================================
      SUBROUTINE INTERNAL_INERTIA_OSCILLATION(UU,VV,TAU)
C----------------------------------------------------------------------
      IMPLICIT NONE
C***************************************************************
C  Calculation of UU,VV - baroclinic velocities on grid "C" :
C
C    Implicit time scheme:
C    (UU-UU0)/tau - l*VV = 0
C
C    (VV-VV0)/tau + l*UU = 0
C
C***************************************************************

      INCLUDE '0COM.INC'
      INCLUDE '0ADEXPIMP.INC'

      REAL  TAU

      REAL   UU(NX,NY,NZ), !    VELOCITIES IN X DIRECTION
     &       VV(NX,NY,NZ)  !    VELOCITIES IN Y DIRECTION

      REAL,ALLOCATABLE:: UU0(:,:,:),VV0(:,:,:)

      REAL DENOM1,DENOM2,DENOM3
      INTEGER M, N, K

      ALLOCATE (UU0(NX,NY,NZ),  !OLD VELOCITIES IN X DIRECTION
     &          VV0(NX,NY,NZ))  !OLD VELOCITIES IN Y DIRECTION
      UU0=UU
      VV0=VV

C U CORIOLIS ROTATION
!$OMP PARALLEL DO
!$OMP&PRIVATE(N,M,K,DENOM1,DENOM2,DENOM3)       
      DO N=NNN,NN
       DO M=MMM,MM

        IF (LCU(M,N).GT.0.5) THEN
        
c       DENOM1=LCV(M,N)+LCV(M+1,N)+LCV(M,N-1)+LCV(M+1,N-1)+1E-10
        DENOM1=4.0

	  DENOM2=1.0+TAU**2*(RLH(M,N)**2+RLH(M,N-1)**2)/2.0
     &                   *ADEQ_IMP**2
	  DENOM3=1.0-TAU**2*(RLH(M,N)**2+RLH(M,N-1)**2)/2.0
     &                   *ADEQ_IMP*ADEQ_EXP

        DO K=1,NZ
          UU(M,N,K)=1.0/DENOM2  *  (
     &          UU0(M,N,K)*DENOM3 + TAU/DENOM1

     &        *( RLH(M,N  )*(VV0(M  ,N  ,K)+VV0(M+1,N  ,K))
     &         + RLH(M,N-1)*(VV0(M  ,N-1,K)+VV0(M+1,N-1,K)) )
     &                                              )
        ENDDO

        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

C V CORIOLIS ROTATION
!$OMP PARALLEL DO
!$OMP&PRIVATE(N,M,K,DENOM1,DENOM2,DENOM3)       
      DO N=NNN,NN
       DO M=MMM,MM

        IF (LCV(M,N).GT.0.5) THEN

c        DENOM1=LCU(M,N)+LCU(M,N+1)+LCU(M-1,N)+LCU(M-1,N+1)+1E-10
        DENOM1=4.0
	  
        DENOM2=1.0+TAU**2*(RLH(M,N)**2+RLH(M-1,N)**2)/2.0
     &                   *ADEQ_IMP**2
	  DENOM3=1.0-TAU**2*(RLH(M,N)**2+RLH(M-1,N)**2)/2.0
     &                   *ADEQ_IMP*ADEQ_EXP

        DO K=1,NZ
          VV(M,N,K)=1.0/DENOM2  *  (
     &          VV0(M,N,K)*DENOM3 - TAU/DENOM1

     &        *( RLH(M  ,N)*(UU0(M  ,N  ,K)+UU0(M  ,N+1,K))
     &         + RLH(M-1,N)*(UU0(M-1,N  ,K)+UU0(M-1,N+1,K)) )
     &                                              )


        ENDDO

        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF(MMD.NE.0) THEN
          CALL CYCLIZE(UU,NX,NY,NZ,MMM,MM)
          CALL CYCLIZE(VV,NX,NY,NZ,MMM,MM)
      END IF
      
      DEALLOCATE(VV0,UU0)

      RETURN
      END
C======================================================================

C     CALCULATING ATMPSPHERIC PRESSURE GRADIENTS
      SUBROUTINE ATM_PRESSURE_GRADIENTS(SLPR,DINX,DINY)
      IMPLICIT NONE
	
      INCLUDE '0COM.INC'
      REAL SLPR(NX,NY)
	REAL(8) DINX(NX,NY),DINY(NX,NY)
	INTEGER M,N

!$OMP PARALLEL DO
	DO N=NNN,NN
        DO M=MMM,MM
          IF(LCU(M,N).GT.0.5) THEN
            DINX(M,N)=-DBLE(10.0/RH0
     &               *(SLPR(M+1,N)-SLPR(M,N))/DXT(M,N)/RN)
          ELSE
	      DINX(M,N)=0.0
          END IF

          IF(LCV(M,N).GT.0.5) THEN
            DINY(M,N)=-DBLE(10.0/RH0
     &               *(SLPR(M,N+1)-SLPR(M,N))/DYT(M,N)/RN)
          ELSE
	      DINY(M,N)=0.0
          END IF

        END DO
	END DO
!$OMP END PARALLEL DO      
      IF(MMD.NE.0) THEN
          CALL CYCLIZE8(DINX,NX,NY, 1,MMM,MM)
          CALL CYCLIZE8(DINY,NX,NY, 1,MMM,MM)
      END IF
	
      RETURN
	END
C======================================================================
      SUBROUTINE SPHERIC_ADDV(UU,VV,UUT,VVT,TAU)
C----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE '0COM.INC'
      INCLUDE '0ADEXPIMP.INC'

      REAL  TAU
      REAL UU(NX,NY,NZ),UUT(NX,NY,NZ), !VELOCITIES IN X DIRECTION
     &     VV(NX,NY,NZ),VVT(NX,NY,NZ)  !VELOCITIES IN Y DIRECTION
      
	REAL,ALLOCATABLE:: UU0(:,:,:),VV0(:,:,:), SPA(:,:,:)
      REAL DENOM1,DENOM2,DENOM3
      INTEGER M, N, K

      ALLOCATE (UU0(NX,NY,NZ),  !OLD VELOCITIES IN X DIRECTION
     &          VV0(NX,NY,NZ),  !OLD VELOCITIES IN Y DIRECTION
     &          SPA(NX,NY,NZ) ) !SPHERICAL ADVECTION TERM
      UU0=UU
      VV0=VV

CCCC    SPHERICAL ADVECTION TERM CALCULATION
!$OMP PARALLEL DO
!$OMP&PRIVATE(N,M,K)      
      DO N=NNN,NN
	   DO M=MMM,MM

            DO K=1,NZ	
             SPA(M,N,K)=1.0/(DXB(M,N)*DYB(M,N)*RN)
     &      *  (   (VVT(M+1,N  ,K)*DYT(M+1,N  )-VVT(M,N,K)*DYT(M,N))
     &           - (UUT(M  ,N+1,K)*DXT(M  ,N+1)-UUT(M,N,K)*DXT(M,N))  )

     &      -  (   (VVT(M+1,N  ,K)-VVT(M  ,N  ,K))/DXB(M,N)/RN
     &           - (UUT(M  ,N+1,K)-UUT(M  ,N  ,K))/DYB(M,N)/RN   )  
            END DO
	   END DO
      END DO
!$OMP END PARALLEL DO
     
      IF(MMD.NE.0) THEN
          CALL CYCLIZE(SPA,NX,NY,NZ,MMM,MM)
      END IF

C U A LA CORIOLIS ROTATION

!$OMP PARALLEL DO
!$OMP&PRIVATE(N,M,K,DENOM1,DENOM2,DENOM3)      
      DO N=NNN,NN
       DO M=MMM,MM

        IF (LCU(M,N).GT.0.5) THEN
        
c        DENOM1=LCV(M,N)+LCV(M+1,N)+LCV(M,N-1)+LCV(M+1,N-1)+1E-10    
        DENOM1=4.0

        DO K=1,NZ
	  
         DENOM2=1.0+TAU**2*(SPA(M,N,K)**2+SPA(M,N-1,K)**2)/2.0
     &                    *ADEQ_IMP**2 
         DENOM3=1.0-TAU**2*(SPA(M,N,K)**2+SPA(M,N-1,K)**2)/2.0
     &                    *ADEQ_IMP*ADEQ_EXP 

          UU(M,N,K)=1.0/DENOM2  *  (
     &          UU0(M,N,K)*DENOM3 + TAU/DENOM1
     &      *( SPA(M,N  ,K)*(VV0(M  ,N  ,K)+VV0(M+1,N  ,K))
     &       + SPA(M,N-1,K)*(VV0(M  ,N-1,K)+VV0(M+1,N-1,K)) )
     &                                                                 )
        ENDDO

        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO


C V A LA CORIOLIS ROTATION
!$OMP PARALLEL DO
!$OMP&PRIVATE(N,M,K,DENOM1,DENOM2,DENOM3)       
      DO N=NNN,NN
       DO M=MMM,MM

        IF (LCV(M,N).GT.0.5) THEN

c        DENOM1=LCU(M,N)+LCU(M,N+1)+LCU(M-1,N)+LCU(M-1,N+1)+1E-10 
        DENOM1=4.0

        DO K=1,NZ
	  
        DENOM2=1.0+TAU**2*(SPA(M,N,K)**2+SPA(M-1,N,K)**2)/2.0
     &                    *ADEQ_IMP**2 
        DENOM3=1.0-TAU**2*(SPA(M,N,K)**2+SPA(M-1,N,K)**2)/2.0
     &                    *ADEQ_IMP*ADEQ_EXP 

          VV(M,N,K)=1.0/DENOM2  *  (
     &        VV0(M,N,K)*DENOM3 - TAU/DENOM1
     &      *( SPA(M  ,N,K)*(UU0(M  ,N  ,K)+UU0(M  ,N+1,K))
     &       + SPA(M-1,N,K)*(UU0(M-1,N  ,K)+UU0(M-1,N+1,K)) )
C     &      *( SPA(M  ,N,K)*HH(M  ,N)*(DXT(M  ,N  )*UU0(M  ,N  ,K)
C     &                                +DXT(M  ,N+1)*UU0(M  ,N+1,K))
C     &       + SPA(M-1,N,K)*HH(M-1,N)*(DXT(M-1,N  )*UU0(M-1,N  ,K)
C     &                                +DXT(M-1,N+1)*UU0(M-1,N+1,K)) )
     &                                              )
        ENDDO

        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF(MMD.NE.0) THEN
          CALL CYCLIZE(UU,NX,NY,NZ,MMM,MM)
          CALL CYCLIZE(VV,NX,NY,NZ,MMM,MM)
      END IF

      DEALLOCATE(SPA,VV0,UU0)
      RETURN
      END
C======================================================================
C MODULES FOR SEA LEVEL GRADIENT CALCULATION
      SUBROUTINE SLH_GRADIENTS(SLH,SLHGRX,SLHGRY)
C-------------(C)DIANSKY N.A. dinar@inm.ras.ru,GUSEV A.V---------------
      IMPLICIT NONE
C***************************************************************
C  CALCULATION OF SLH DYNAMIC GRADIENTS 
      INCLUDE '0COM.INC'

      REAL*8   SLH(NX,NY), !SEA LEVEL HEIGHT (INPUT)
     &	  SLHGRX(NX,NY), !SLH DYNAMIC GRADSIENT IN X DIRECTION(CALCULATED)
     &	  SLHGRY(NX,NY)  !SLH DYNAMIC GRADSIENT IN Y DIRECTION(CALCULATED)

	INTEGER M, N

!$OMP PARALLEL DO PRIVATE(M,N)
      DO N=NNN,NN
       DO M=MMM,MM
C  COMPUTING IN X-DIRECTION
        SLHGRX(M,N)= DBLE( LCU(M,N)*GRV/DXT(M,N)/RN ) *
     &                  ( SLH(M+1,N)-SLH(M,N) )

C  COMPUTING IN Y-DIRECTION
	  SLHGRY(M,N)= DBLE( LCV(M,N)*GRV/DYT(M,N)/RN) *
     &                  ( SLH(M,N+1)-SLH(M,N) ) 
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF(MMD.NE.0) THEN 
          CALL CYCLIZE8(SLHGRX,NX,NY, 1,MMM,MM)
          CALL CYCLIZE8(SLHGRY,NX,NY, 1,MMM,MM)
      END IF

      RETURN
      END
C======================================================================

C     CALCULATING ATMPSPHERIC PRESSURE GRADIENTS FOR ICE DYNAMICS
      SUBROUTINE ATM_PRESSURE_GRADIENTS_4ICE(SLPR,DINX,DINY,
     &                                       AISTOT,MISTOT,MITOT)
      IMPLICIT NONE
	
      INCLUDE '0COM.INC'
      REAL SLPR(NX,NY),AISTOT(NX,NY),MISTOT(NX,NY),MITOT(NX,NY)
	REAL(8) DINX(NX,NY),DINY(NX,NY)
	INTEGER M,N

!$OMP PARALLEL DO PRIVATE(M,N)
	DO N=NNN,NN
        DO M=MMM,MM
          IF(LCU(M,N).GT.0.5) THEN

C           IF((MISTOT(M,N)+MISTOT(M+1,N)).GT. 20.0) THEN

            DINX(M,N)=-DBLE(10.0*(MITOT(M,N)+MITOT(M+1,N))
     &               *(SLPR(M+1,N)-SLPR(M  ,N))
     &               /(DXT(M,N)*RN*(MISTOT(M+1,N)+MISTOT(M,N)+20.0)))
C           END IF

          ELSE
	      DINX(M,N)=0.0
          END IF

          IF(LCV(M,N).GT.0.5) THEN

C           IF((MISTOT(M,N)+MISTOT(M,N+1)).GT. 20.0) THEN

            DINY(M,N)=-DBLE(10.0*(MITOT(M,N)+MITOT(M,N+1))
     &               *(SLPR(M,N+1)-SLPR(M  ,N))
     &               /(DYT(M,N)*RN*(MISTOT(M,N+1)+MISTOT(M,N)+20.0)))
C           END IF
          ELSE
	      DINY(M,N)=0.0
          END IF

        END DO
	END DO
!$OMP END PARALLEL DO     

      IF(MMD.NE.0) THEN
          CALL CYCLIZE8(DINX,NX,NY, 1,MMM,MM)
          CALL CYCLIZE8(DINY,NX,NY, 1,MMM,MM)
      END IF
	
      RETURN
	END
C======================================================================
      SUBROUTINE ICE_CORIOLIS_ADAPT(UU,VV,TAU,AISTOT,MISTOT)
C----------------------------------------------------------------------
      IMPLICIT NONE
C***************************************************************
C  Calculation of UU,VV - baroclinic velocities on grid "C" :
C
C    Implicit time scheme:
C    (UU-UU0)/tau - l*VV = 0
C
C    (VV-VV0)/tau + l*UU = 0
C
C***************************************************************

      INCLUDE '0COM.INC'
      REAL  TAU

      REAL   UU(NX,NY), !    VELOCITIES IN X DIRECTION
     &       VV(NX,NY)  !    VELOCITIES IN Y DIRECTION
      REAL AISTOT(NX,NY),MISTOT(NX,NY)
      REAL,ALLOCATABLE:: UU0(:,:),VV0(:,:)

      REAL DENOM2
      INTEGER M, N

      ALLOCATE (UU0(NX,NY),  !OLD VELOCITIES IN X DIRECTION
     &          VV0(NX,NY))  !OLD VELOCITIES IN Y DIRECTION
      UU0=UU
      VV0=VV

C U CORIOLIS ROTATION
!$OMP PARALLEL DO
!$OMP&PRIVATE(N,M,DENOM2)       
      DO N=NNN,NN
       DO M=MMM,MM

        IF (LCU(M,N).GT.0.5) THEN
C           IF((MISTOT(M,N)+MISTOT(M+1,N)).GT.20.0) THEN           

	  DENOM2=1.0+TAU**2*(RLH(M,N)**2+RLH(M,N-1)**2)/2.0

          UU(M,N)=1.0/DENOM2  *  (
     &          UU0(M,N) + TAU/4.0

     &        *( RLH(M,N  )*(VV0(M  ,N  )+VV0(M+1,N  ))
     &         + RLH(M,N-1)*(VV0(M  ,N-1)+VV0(M+1,N-1)) )

     &                                              )
C           END IF

        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

C V CORIOLIS ROTATION
!$OMP PARALLEL DO
!$OMP&PRIVATE(N,M,DENOM2)       
      DO N=NNN,NN
       DO M=MMM,MM

        IF (LCV(M,N).GT.0.5) THEN
C           IF((MISTOT(M,N)+MISTOT(M,N+1)).GT.20.0) THEN
      
	  DENOM2=1.0+TAU**2*(RLH(M,N)**2+RLH(M-1,N)**2)/2.0

          VV(M,N)=1.0/DENOM2  *  (
     &          VV0(M,N) - TAU/4.0

     &        *( RLH(M  ,N)*(UU0(M  ,N  )+UU0(M  ,N+1))
     &         + RLH(M-1,N)*(UU0(M-1,N  )+UU0(M-1,N+1)) )
     &                                              )


C           ENDIF

        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF(MMD.NE.0) THEN
          CALL CYCLIZE(UU,NX,NY,1,MMM,MM)
          CALL CYCLIZE(VV,NX,NY,1,MMM,MM)
      END IF
      
      DEALLOCATE(VV0,UU0)

      RETURN
      END
C======================================================================
      SUBROUTINE ICE_SPHERIC_ADDV(UU,VV,UUT,VVT,TAU)
C----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE '0COM.INC'
      REAL  TAU
      REAL UU(NX,NY),UUT(NX,NY), !VELOCITIES IN X DIRECTION
     &     VV(NX,NY),VVT(NX,NY)  !VELOCITIES IN Y DIRECTION
      
	REAL,ALLOCATABLE:: UU0(:,:),VV0(:,:), SPA(:,:)
      REAL DENOM1,DENOM2
      INTEGER M, N

      ALLOCATE (UU0(NX,NY),  !OLD VELOCITIES IN X DIRECTION
     &          VV0(NX,NY),  !OLD VELOCITIES IN Y DIRECTION
     &          SPA(NX,NY) ) !SPHERICAL ADVECTION TERM
      UU0=UU
      VV0=VV

CCCC    SPHERICAL ADVECTION TERM CALCULATION
!$OMP PARALLEL DO
!$OMP&PRIVATE(N,M)      
      DO N=NNN,NN
	   DO M=MMM,MM
	
             SPA(M,N)=1.0/(DXB(M,N)*DYB(M,N)*RN)
     &      *  (   (VVT(M+1,N  )*DYT(M+1,N  )-VVT(M,N)*DYT(M,N))
     &           - (UUT(M  ,N+1)*DXT(M  ,N+1)-UUT(M,N)*DXT(M,N))  )

     &      -  (   (VVT(M+1,N  )-VVT(M  ,N  ))/DXB(M,N)/RN
     &           - (UUT(M  ,N+1)-UUT(M  ,N  ))/DYB(M,N)/RN   )  
	   END DO
      END DO
!$OMP END PARALLEL DO
     
      IF(MMD.NE.0) THEN
          CALL CYCLIZE(SPA,NX,NY,1,MMM,MM)
      END IF

C U A LA CORIOLIS ROTATION

!$OMP PARALLEL DO
!$OMP&PRIVATE(N,M,DENOM1,DENOM2)      
      DO N=NNN,NN
       DO M=MMM,MM

        IF (LCU(M,N).GT.0.5) THEN
        
C	  DENOM1=DYH(M,N)*HHU(M,N)
        DENOM1=1.0        
	  
         DENOM2=1.0+TAU**2*(SPA(M,N)**2+SPA(M,N-1)**2)/2.0
          UU(M,N)=1.0/DENOM2  *  (
     &          UU0(M,N) + TAU/(4.0*DENOM1)
     &      *( SPA(M,N  )*(VV0(M  ,N  )+VV0(M+1,N  ))
     &       + SPA(M,N-1)*(VV0(M  ,N-1)+VV0(M+1,N-1)) )
C     &      *( SPA(M,N  ,K)*HH(M,N  )*(DYT(M  ,N  )*VV0(M  ,N  ,K)
C     &                                +DYT(M+1,N  )*VV0(M+1,N  ,K))
C     &       + SPA(M,N-1,K)*HH(M,N-1)*(DYT(M  ,N-1)*VV0(M  ,N-1,K)
C     &                                +DYT(M+1,N-1)*VV0(M+1,N-1,K)) )
     &                                                                 )

        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO


C V A LA CORIOLIS ROTATION
!$OMP PARALLEL DO
!$OMP&PRIVATE(N,M,DENOM1,DENOM2)       
      DO N=NNN,NN
       DO M=MMM,MM

        IF (LCV(M,N).GT.0.5) THEN

C	  DENOM1=DXH(M,N)*HHV(M,N)
        DENOM1=1.0
	  
        DENOM2=1.0+TAU**2*(SPA(M,N)**2+SPA(M-1,N)**2)/2.0
          VV(M,N)=1.0/DENOM2  *  (
     &        VV0(M,N) - TAU/(4.0*DENOM1)
     &      *( SPA(M  ,N)*(UU0(M  ,N  )+UU0(M  ,N+1))
     &       + SPA(M-1,N)*(UU0(M-1,N  )+UU0(M-1,N+1)) )
C     &      *( SPA(M  ,N,K)*HH(M  ,N)*(DXT(M  ,N  )*UU0(M  ,N  ,K)
C     &                                +DXT(M  ,N+1)*UU0(M  ,N+1,K))
C     &       + SPA(M-1,N,K)*HH(M-1,N)*(DXT(M-1,N  )*UU0(M-1,N  ,K)
C     &                                +DXT(M-1,N+1)*UU0(M-1,N+1,K)) )
     &                                              )

        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF(MMD.NE.0) THEN
          CALL CYCLIZE(UU,NX,NY,1,MMM,MM)
          CALL CYCLIZE(VV,NX,NY,1,MMM,MM)
      END IF

      DEALLOCATE(SPA,VV0,UU0)
      RETURN
      END
