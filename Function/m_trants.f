C======================================================================
      SUBROUTINE TRXPRP(FF,UU,AX,TAU,SLH,SLH0)
      IMPLICIT NONE
C----------------------------------------------------------------------
C X-TRANSPORT AND DIFFUSION ON T - GRID.
C X-TRANPORT SUPPOSES BOUNDARY VELOCITIES & FI TO BE ZERO.
      INCLUDE '0COM.INC'
C IMPLICIT - EXPLICIT CONTROL TRANSPORT T & S TIME SCHEME:
      INCLUDE '0TRTSS.INC'

      REAL FF(NX,NY,NZ),UU(NX,NY,NZ),AX(NX,NY,NZ)
      REAL A(NX),B(NX),C(NX),ETA(NX),RKSI(NX)
      REAL(8) SLH(NX,NY),SLH0(NX,NY)
      INTEGER M, N, K, IG, II, JJ
      INTEGER M1, M9, MLOOP
      REAL    BP, BP0,DP, DM, PP, PM
      REAL    TAU

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,BP,BP0,DP,DM,PP,PM,A,B,C,ETA,RKSI,
!$OMP&       IG,II,JJ,MLOOP,M1,M9)      
      DO K=1,NZ
       DO N=NNN,NN
        DO IG=1,ABS(LRX(N))
         II=IIX(IG,N)
         JJ=JJX(IG,N)
         IF (JJ.GT.MM) JJ = JJ - MMD

         DO MLOOP=II,JJX(IG,N)
          M  = MLOOP
          M1 = MLOOP - 1
          M9 = MLOOP + 1
          IF ( M.GT.MM) M  = M  - MMD
          IF (M1.GT.MM) M1 = M1 - MMD
          IF (M9.GT.MM) M9 = M9 - MMD

          BP = (HHQ(M,N)-SNGL(SLH (M,N)))*DX(M,N)*DY(M,N)*RN / TAU
          BP0= (HHQ(M,N)-SNGL(SLH0(M,N)))*DX(M,N)*DY(M,N)*RN / TAU

          DP = AX(M ,N,K)*DYH(M ,N)/DXT(M ,N)/RN
     &        *HHU(M ,N)

          DM = AX(M1,N,K)*DYH(M1,N)/DXT(M1,N)/RN
     &        *HHU(M1,N)

          PP = UU(M ,N,K)*DYH(M ,N)
     &        *HHU(M ,N)/2.0

          PM = UU(M1,N,K)*DYH(M1,N)
     &        *HHU(M1,N)/2.0

C LAST "2" ABOVE IS DUE TO APROXIMATION IN SPACE.
          C(M) =  PP * TRTS_IMP - DP
          A(M) = -PM * TRTS_IMP - DM
          B(M) =  BP + DP + DM
          ETA(M) = BP0* FF(M,N,K) - PP * TRTS_EXP * FF(M9,N,K)
     &                            + PM * TRTS_EXP * FF(M1,N,K)
         ENDDO
C IMPLICIT: A(M)=-DM-PM;C(M)=-DP+PP;B(M)=BP+DM+DP;ETA(M)=BP*FF(M,N,K)
C IF (M=II) : PM = 0, A = -DM ; IF (M=JJ) : PP = 0, C = -DP.

         IF (LRX(N).GT.0) THEN     !(2nd BOUNDARY CONDITION ONLY)
           B(II) = B(II) + A(II)
           B(JJ) = B(JJ) + C(JJ)
          CALL FACTOR2(NX,A,B,C,ETA,RKSI,II,JJX(IG,N))
         ELSE
          CALL FACTORC(A,B,C,ETA,RKSI,MMM,MM)
         ENDIF
         DO MLOOP=II,JJX(IG,N)
          M  = MLOOP
          IF (M.GT.MM) M  = M  - MMD
          FF(M,N,K) = RKSI(M)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END
C======================================================================
      SUBROUTINE TRYPR(FF,VV,AY1,TAU,SLH,SLH0)
      IMPLICIT NONE
C----------------------------------------------------------------------
C Y-TRANSPORT AND DIFFUSION.
C Y-TRANPORT SUPPOSES BOUNDARY VELOCITIES & FI TO BE ZERO.
C WATER BOUNDARIES AND SOLID BOUNDARIES ARE DIFFERENT.
      INCLUDE '0COM.INC'
C IMPLICIT - EXPLICIT CONTROL TRANSPORT T & S TIME SCHEME:
      INCLUDE '0TRTSS.INC'
      REAL FF(NX,NY,NZ),VV(NX,NY,NZ),
     &     AY1(NX,NY,NZ)
      REAL    A(NY),B(NY),C(NY),ETA(NY),RKSI(NY)
	REAL(8) SLH(NX,NY),SLH0(NX,NY)
      INTEGER M, N, K, IG, II, JJ
	REAL    BP, BP0, DP, DM, PP, PM
      REAL    TAU

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,BP,BP0,DP,DM,PP,PM,A,B,C,ETA,RKSI,
!$OMP&         IG,II,JJ)
      DO K=1,NZ
       DO M=MMM,MM
        DO IG=1,LRY(M)
         II=IIY(IG,M)
         JJ=JJY(IG,M)
         DO N=II,JJ
          BP = (HHQ(M,N)-SNGL(SLH (M,N)))*DX(M,N)*DY(M,N)*RN / TAU
          BP0= (HHQ(M,N)-SNGL(SLH0(M,N)))*DX(M,N)*DY(M,N)*RN / TAU

          DP = AY1(M,N  ,K)*DXH(M,N  )/DYT(M,N  )/RN
     &       *HHV(M,N  )

          DM = AY1(M,N-1,K)*DXH(M,N-1)/DYT(M,N-1)/RN
     &       *HHV(M,N-1)

          PP =  VV(M,N  ,K)*DXH(M,N  )
     &       *HHV(M,N  )/2.0

          PM =  VV(M,N-1,K)*DXH(M,N-1)
     &       *HHV(M,N-1)/2.0
C LAST "2" ABOVE IS DUE TO APROXIMATION IN SPACE.
          C(N) =  PP * TRTS_IMP - DP
          A(N) = -PM * TRTS_IMP - DM
          B(N) =  BP + DP + DM
          ETA(N) = BP0* FF(M,N,K) - PP * TRTS_EXP * FF(M,N+1,K)
     &                            + PM * TRTS_EXP * FF(M,N-1,K)
         END DO
C IMPLICIT: A(N)=-DM-PM;C(N)=-DP+PP;B(N)=BP+DM+DP;ETA(N)=BP*FF(M,N,K)
C IF (N=II) : PM = 0, A = -DM ; IF (N=JJ) : PP = 0, C = -DP.
         
C        2nd BOUNDARY CONDITION ONLY
           B(II) = B(II) + A(II)
           B(JJ) = B(JJ) + C(JJ)
  
         CALL FACTOR(NY,A,B,C,ETA,RKSI,II,JJ)
         DO N=II,JJ
          FF(M,N,K) = RKSI(N)
         END DO
        END DO
       END DO
      END DO
!$OMP END PARALLEL DO
      RETURN
      END
C======================================================================
      SUBROUTINE TRZ3PR(FF,FF0,WW,A_Z,TAU,IGRZ,SLH,SLH0)
      IMPLICIT NONE
C----------------------------------------------------------------------
C Z-TRANSPORT AND 3D VERTICAL DIFFUSION..
C ADVECTION IS CRANK-NICOLSON, DIFFUSION IS PURE IMPLICIT.
C----------------------------------------------------------------------
      INCLUDE '0COM.INC'
C IMPLICIT - EXPLICIT CONTROL TRANSPORT T & S TIME SCHEME:
      INCLUDE '0TRTSS.INC'

      REAL FF(NX,NY,NZ),FF0(NX,NY),WW(NX,NY,NZ+1),A_Z(NX,NY,NZ+1)
      REAL A(NZ),B(NZ),C(NZ),ETA(NZ),RKSI(NZ)
      INTEGER IGRZ(NX,NY) !TYPE OF SEA SURFACE BOUNDARY CONDITION (1/2)
      INTEGER M, N, K
      REAL    BP,BP0, DP, DM, PP, PM
      REAL    TAU
	REAL(8) SLH(NX,NY),SLH0(NX,NY)
!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,BP,BP0,DP,DM,PP,PM,A,B,C,ETA,RKSI)
      DO N=NNN,NN
       DO M=MMM,MM
        IF (LU(M,N).GT.0.5) THEN
         BP = (HHQ(M,N)-SNGL(SLH (M,N))) / TAU
         BP0= (HHQ(M,N)-SNGL(SLH0(M,N))) / TAU
C INTERNAL POINTS.
         DO K=2,NZ-1
          DP = A_Z(M,N,K+1) / HHQ(M,N) / HZT(K+1) / DZ(K)
          DM = A_Z(M,N,K)   / HHQ(M,N) / HZT(K  ) / DZ(K)
          PP =  WW(M,N,K+1) / DZ(K) / 2.
          PM =  WW(M,N,K  ) / DZ(K) / 2.
C "2." ABOVE IS DUE TO APROXIMATION.
          C(K) =  PP * TRTS_IMP- DP
          A(K) = -PM * TRTS_IMP- DM
          B(K) =  BP + DP + DM
          ETA(K) = BP0*FF(M,N,K) - PP * TRTS_EXP*FF(M,N,K+1)
     +                           + PM * TRTS_EXP*FF(M,N,K-1)
         ENDDO
C IMPLICIT: A(K)=-DM-PM;C(K)=-DP+PP;B(K)=BP+DP+DM;ETA(K)=BP*FF(M,N,K)
C SURFACE POINT: TWO CASES.
         K = 1
         DP = A_Z(M,N,K+1)/ HHQ(M,N) / HZT(K+1) / DZ(K)
         DM = A_Z(M,N,K)  / HHQ(M,N) / HZT(K  ) / DZ(K)
         PP =  WW(M,N,K+1)/ DZ(K) / 2.
         A(K) = 0.0    !NOT USE IN FACTOR      
         C(K) =  PP * TRTS_IMP - DP
C AZ -> AZD, PM IS 0 DUE TO WW(M,N,1) = 0.
         IF(IGRZ(M,N).EQ.1) THEN
          B(K) =  BP + DP + DM
          ETA(K) = BP0*FF(M,N,K) - PP * TRTS_EXP * FF(M,N,K+1)
     +                                    + DM * FF0(M,N)
         ELSEIF(IGRZ(M,N).EQ.2) THEN
          B(K) =  BP + DP
          ETA(K) = BP0*FF(M,N,K) - PP * TRTS_EXP * FF(M,N,K+1)
     +                                         + FF0(M,N) / DZ(K)
         ENDIF
C BOTTOM POINT. FLUX IS ZERO FOR ALL FIELDS.
         K = NZ
         DM = A_Z(M,N,K) / HHQ(M,N) / HZT(K  ) / DZ(K)
         PM =  WW(M,N,K) / DZ(K) / 2.
         A(K) = -PM * TRTS_IMP - DM
         C(K)=0.0     !NOT USE IN FACTOR  
         B(K) =  BP + DM
         ETA(K) = BP0* FF(M,N,K) + PM * TRTS_EXP* FF(M,N,K-1)

         CALL FACTOR(NZ,A,B,C,ETA,RKSI,1,NZ)
         DO K=1,NZ
          FF(M,N,K)=RKSI(K)
         ENDDO
        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END
C========================================================27.11.99 19:43
      SUBROUTINE SWRADPEN(TT,DIVSWRAD,QSWR,LU,TAU,NX,NY,NZ)
	IMPLICIT NONE
C ADD SOLAR SHORTWAVE HEATING PENETRATED BELOW THE OCEAN SURFACE.
      INTEGER NX,NY,NZ
      REAL DIVSWRAD(NZ,NX,NY),TT(NX,NY,NZ),LU(NX,NY),QSWR(NX,NY),TAU
      INTEGER  M,N,K

!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO N = 2, NY-1
      DO M = 2, NX-1

      IF(LU(M,N).GT.0.5) THEN
      DO K=1,NZ
C ADD SOLAR SHORTWAVE HEATING PENETRATED BELOW THE OCEAN SURFACE.
         TT(M,N,K)=TT(M,N,K)+TAU*0.4*QSWR(M,N)*DIVSWRAD(K,M,N)
      ENDDO
      END IF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
C======================================================================
      SUBROUTINE TS_VAR_BY_SSH_CHANGE(FF,SLH,SLH0)
      IMPLICIT NONE
C----------------------------------------------------------------------
C TEMPERATURE AND SALINITY VARIATION DUE TO
C SEA LEVEL CHANGE
C----------------------------------------------------------------------
      INCLUDE '0COM.INC'

      REAL FF(NX,NY,NZ)
      INTEGER M, N, K
      REAL(8) BP,BP0,DV, DENOM
	REAL(8) SLH(NX,NY),SLH0(NX,NY)
!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,BP,BP0,DV,DENOM)
       DO N=NNN,NN
        DO M=MMM,MM
         IF (LU(M,N).GT.0.5) THEN
           BP = DBLE(HHQ(M,N))-SLH (M,N)
           BP0= DBLE(HHQ(M,N))-SLH0(M,N)       
           DV =     (SLH(M,N) -SLH0(M,N))/2.0D0
	   
         DENOM=(BP0- DMIN1(DV,0.0D0))
     &        /(BP + DMAX1(DV,0.0D0))
     
          DO K=1,NZ
            FF(M,N,K)=SNGL(DENOM*DBLE(FF(M,N,K)))
          END DO   
         
         ENDIF
        ENDDO
       ENDDO
!$OMP END PARALLEL DO

      RETURN
      END

C==============================================================
C  TRANSPORT UPWIND MPDATA-SCHEME USING DIFFUSIVITY REDUCTION
      SUBROUTINE TRAN_TS_MPDATA(FF,FF0,TAU,UU,VV,WW,ADDVAL,NITER,
     &                                     A_Z,IGRZ,SLH,SLH0)
	IMPLICIT NONE
      INCLUDE '0COM.INC'
C IMPLICIT - EXPLICIT CONTROL TRANSPORT T & S TIME SCHEME:
      INCLUDE '0TRTSS.INC'
            
      REAL  FF(NX,NY,NZ),FF0(NX,NY),TAU,ADDVAL   !ADDVAL IS VALUE ADDED TO FF TO MAKE IT POSITIVE
	REAL  UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ+1), ! TRANSPORTING VELOCITIES
     &     A_Z(NX,NY,NZ+1) 
      INTEGER NITER       !NITER - NUMBER OF ITERATIONS FOR REDUCIND DIFFUSIVITY


      INTEGER IGRZ(NX,NY)
      REAL FX_P,FX_M,FY_P,FY_M,FZ_P,FZ_M   !FLUXES THROUGH CELL EDGES
	REAL(8) SLH(NX,NY),SLH0(NX,NY)     !SEA SURFACE HEIGTS FOR CONSERVATION

C HELP VARIABLES FOR FACTOR PROCEDURE 
      REAL    BP,BP0, DP, DM, PP, PM
      REAL DISCONT
      REAL A(NZ),B(NZ),C(NZ),ETA(NZ),RKSI(NZ)        

      INTEGER M,N,K 
           
      REAL,allocatable:: FO(:,:,:)               !OLD VALUE OF FF

      allocate(FO(NX,NY,NZ))

C ADDING ADDVAL AND STORING OLD VALUE OF FF 

      FO=FF+ADDVAL

c!$OMP PARALLEL DO      
c      DO K=1,NZ
c         FO(:,:,K)=FF(:,:,K)+ADDVAL
c	END DO
c!$OMP END PARALLEL DO

C   THE MAIN TRANSPORT PROCEDURE

!$OMP PARALLEL DO PRIVATE(M,N,K,BP,BP0,FX_P,FX_M,FY_P,FY_M,
!$OMP&       DP,DM,PP,PM,DISCONT,A,B,C,ETA,RKSI )
	DO N=NNN,NN
         DO M=MMM,MM
          IF(LU(M,N).GT.0.5) THEN
           BP = (HHQ(M,N)-SNGL(SLH (M,N))) / TAU
           BP0= (HHQ(M,N)-SNGL(SLH0(M,N))) / TAU
            
            DO K=2,NZ-1

	       FX_P=(UU(M  ,N  ,K)+ABS(UU(M  ,N  ,K)))
     &           *DYH(M  ,N    ) * HHU (M  ,N    )*FO(M  ,N  ,K)
     &        +   (UU(M  ,N  ,K)-ABS(UU(M  ,N  ,K)))
     &           *DYH(M  ,N    ) * HHU (M  ,N    )*FO(M+1,N  ,K)

	       FX_M=(UU(M-1,N  ,K)+ABS(UU(M-1,N  ,K)))
     &           *DYH(M-1,N    ) * HHU (M-1,N    )*FO(M-1,N  ,K)
     &        +   (UU(M-1,N  ,K)-ABS(UU(M-1,N  ,K)))
     &           *DYH(M-1,N    ) * HHU (M-1,N    )*FO(M  ,N  ,K)

	       FY_P=(VV(M  ,N  ,K)+ABS(VV(M  ,N  ,K)))
     &           *DXH(M  ,N    ) * HHV (M  ,N    )*FO(M  ,N  ,K)
     &        +   (VV(M  ,N  ,K)-ABS(VV(M  ,N  ,K)))
     &           *DXH(M  ,N    ) * HHV (M  ,N    )*FO(M  ,N+1,K)

	       FY_M=(VV(M  ,N-1,K)+ABS(VV(M  ,N-1,K)))
     &           *DXH(M  ,N-1  ) * HHV (M  ,N-1  )*FO(M  ,N-1,K)
     &        +   (VV(M  ,N-1,K)-ABS(VV(M  ,N-1,K)))
     &           *DXH(M  ,N-1  ) * HHV (M  ,N-1  )*FO(M  ,N  ,K)

          DP = A_Z(M,N,K+1) / HHQ(M,N) / HZT(K+1) / DZ(K)
          DM = A_Z(M,N,K)   / HHQ(M,N) / HZT(K  ) / DZ(K)
          PP =  WW(M,N,K+1) / DZ(K) / 2.
          PM =  WW(M,N,K  ) / DZ(K) / 2. 
          DISCONT= PP - PM           
          
          C(K) =  PP * TRTS_IMP- DP
          A(K) = -PM * TRTS_IMP- DM
c          B(K) =  BP + DP + DM + MAX(DISCONT,0.0)   !explicit-implicit scheme
          B(K) =  BP + DP + DM + DISCONT*TRTS_IMP
          ETA(K) = BP0*FF(M,N,K) - PP * TRTS_EXP*FO(M,N,K+1)
     +                           + PM * TRTS_EXP*FO(M,N,K-1)
c     +                        - MIN(DISCONT,0.0)*FO(M,N,K  ) !explicit-implicit scheme
     +                        -  TRTS_EXP*DISCONT*FO(M,N,K  )  !explicit scheme
     +      -(FX_P - FX_M + FY_P - FY_M)/(DX(M,N)*DY(M,N)*RN*2.0)            
            END DO

C IMPLICIT: A(K)=-DM-PM;C(K)=-DP+PP;B(K)=BP+DP+DM;ETA(K)=BP*FF(M,N,K)
C SURFACE POINT: TWO CASES.
         K = 1
	       FX_P=(UU(M  ,N  ,K)+ABS(UU(M  ,N  ,K)))
     &           *DYH(M  ,N    ) * HHU (M  ,N    )*FO(M  ,N  ,K)
     &        +   (UU(M  ,N  ,K)-ABS(UU(M  ,N  ,K)))
     &           *DYH(M  ,N    ) * HHU (M  ,N    )*FO(M+1,N  ,K)

	       FX_M=(UU(M-1,N  ,K)+ABS(UU(M-1,N  ,K)))
     &           *DYH(M-1,N    ) * HHU (M-1,N    )*FO(M-1,N  ,K)
     &        +   (UU(M-1,N  ,K)-ABS(UU(M-1,N  ,K)))
     &           *DYH(M-1,N    ) * HHU (M-1,N    )*FO(M  ,N  ,K)

	       FY_P=(VV(M  ,N  ,K)+ABS(VV(M  ,N  ,K)))
     &           *DXH(M  ,N    ) * HHV (M  ,N    )*FO(M  ,N  ,K)
     &        +   (VV(M  ,N  ,K)-ABS(VV(M  ,N  ,K)))
     &           *DXH(M  ,N    ) * HHV (M  ,N    )*FO(M  ,N+1,K)

	       FY_M=(VV(M  ,N-1,K)+ABS(VV(M  ,N-1,K)))
     &           *DXH(M  ,N-1  ) * HHV (M  ,N-1  )*FO(M  ,N-1,K)
     &        +   (VV(M  ,N-1,K)-ABS(VV(M  ,N-1,K)))
     &           *DXH(M  ,N-1  ) * HHV (M  ,N-1  )*FO(M  ,N  ,K)

         DP = A_Z(M,N,K+1)/ HHQ(M,N) / HZT(K+1) / DZ(K)
         DM = A_Z(M,N,K)  / HHQ(M,N) / HZT(K  ) / DZ(K)
         PP =  WW(M,N,K+1)/ DZ(K) / 2.
         PM = 0.0
         DISCONT= PP - PM
         
         A(K) = 0.0    !NOT USE IN FACTOR      
         C(K) =  PP * TRTS_IMP - DP
C AZ -> AZD, PM IS 0 DUE TO WW(M,N,1) = 0.
         IF(IGRZ(M,N).EQ.1) THEN

c          B(K) =  BP + DP + DM + MAX(DISCONT,0.0)   !explicit-implicit scheme
          B(K) =  BP + DP + DM + DISCONT* TRTS_IMP
          ETA(K) = BP0*FF(M,N,K) - PP * TRTS_EXP * FO(M,N,K+1)
     +                                    + DM * FF0(M,N) 
     +                        - TRTS_EXP *DISCONT*FO(M,N,K  )  !explicit scheme
     +      -(FX_P - FX_M + FY_P - FY_M)/(DX(M,N)*DY(M,N)*RN)

         ELSEIF(IGRZ(M,N).EQ.2) THEN
          
          B(K) =  BP + DP + DISCONT* TRTS_IMP
          ETA(K) = BP0*FF(M,N,K) - PP * TRTS_EXP * FO(M,N,K+1)
     +                                          + FF0(M,N) / DZ(K)
c     +                        - MIN(DISCONT,0.0)*FO(M,N,K  ) !explicit-implicit scheme
     +                        - TRTS_EXP *DISCONT*FO(M,N,K  )  !explicit scheme
     +       -(FX_P - FX_M + FY_P - FY_M)/(DX(M,N)*DY(M,N)*RN*2.0)   
         
         ENDIF

C BOTTOM POINT. FLUX IS ZERO FOR ALL FIELDS.
         K = NZ
	       FX_P=(UU(M  ,N  ,K)+ABS(UU(M  ,N  ,K)))
     &           *DYH(M  ,N    ) * HHU (M  ,N    )*FO(M  ,N  ,K)
     &        +   (UU(M  ,N  ,K)-ABS(UU(M  ,N  ,K)))
     &           *DYH(M  ,N    ) * HHU (M  ,N    )*FO(M+1,N  ,K)

	       FX_M=(UU(M-1,N  ,K)+ABS(UU(M-1,N  ,K)))
     &           *DYH(M-1,N    ) * HHU (M-1,N    )*FO(M-1,N  ,K)
     &        +   (UU(M-1,N  ,K)-ABS(UU(M-1,N  ,K)))
     &           *DYH(M-1,N    ) * HHU (M-1,N    )*FO(M  ,N  ,K)

	       FY_P=(VV(M  ,N  ,K)+ABS(VV(M  ,N  ,K)))
     &           *DXH(M  ,N    ) * HHV (M  ,N    )*FO(M  ,N  ,K)
     &        +   (VV(M  ,N  ,K)-ABS(VV(M  ,N  ,K)))
     &           *DXH(M  ,N    ) * HHV (M  ,N    )*FO(M  ,N+1,K)

	       FY_M=(VV(M  ,N-1,K)+ABS(VV(M  ,N-1,K)))
     &           *DXH(M  ,N-1  ) * HHV (M  ,N-1  )*FO(M  ,N-1,K)
     &        +   (VV(M  ,N-1,K)-ABS(VV(M  ,N-1,K)))
     &           *DXH(M  ,N-1  ) * HHV (M  ,N-1  )*FO(M  ,N  ,K)

         DM = A_Z(M,N,K) / HHQ(M,N) / HZT(K  ) / DZ(K)
         PM =  WW(M,N,K) / DZ(K) / 2.
	   PP = 0.0
         DISCONT= PP - PM
                  
         A(K) = -PM * TRTS_IMP - DM
         C(K)=0.0     !NOT USE IN FACTOR  
          B(K) =  BP + DP + DM + DISCONT* TRTS_IMP  !explicit-implicit scheme
C         B(K) =  BP + DP + DM 
         ETA(K) = BP0* FF(M,N,K) + PM * TRTS_EXP* FO(M,N,K-1)
c     +                        - MIN(DISCONT,0.0)*FO(M,N,K  ) !explicit-implicit scheme
     +                        - TRTS_EXP *DISCONT*FO(M,N,K  )  !explicit scheme
     +      -(FX_P - FX_M + FY_P - FY_M)/(DX(M,N)*DY(M,N)*RN*2.0) 

         CALL FACTOR(NZ,A,B,C,ETA,RKSI,1,NZ)
         DO K=1,NZ
          FF(M,N,K)=RKSI(K)
         ENDDO

           END IF
         END DO
	END DO
!$OMP END PARALLEL DO

C REMOVING ADDVAL
      FF=FF-ADDVAL
c!$OMP PARALLEL DO      
c      DO K=1,NZ
c         FF (:,:,K)=FF(:,:,K)-ADDVAL
c	END DO
c!$OMP END PARALLEL DO

      deallocate(FO )
      RETURN
	END

C==============================================================
C  TRANSPORT EXPLICIT SCHEME USING SIGMA-DIFFUSION
      SUBROUTINE TRANTS_AD_BASH_SDIFF(FF,TAU,UU,VV,WW,AX,AY,SLH,SLH0)
	IMPLICIT NONE
      INCLUDE '0COM.INC'
      INCLUDE '1TRANTS.INC'

      REAL  FF(NX,NY,NZ),TAU,TAU_INNER
	REAL  UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ+1) ! TRANSPORTING VELOCITIE
      REAL  AX(NX,NY,NZ),AY(NX,NY,NZ)
      REAL FX_P,FX_M,FY_P,FY_M,FZ_P,FZ_M   !FLUXES THROUGH CELL EDGES
	REAL(8) SLH(NX,NY),SLH0(NX,NY),SL_P,SL_M     !SEA SURFACE HEIGTS FOR CONSERVATION    
      
	REAL BP,BP0,RHS
      INTEGER M,N,K,ITER
           
      REAL,allocatable:: FO(:,:,:),FM(:,:,:),FP(:,:,:)            !OLD VALUE OF FF

      allocate(FO(NX,NY,NZ),FM(NX,NY,NZ),FP(NX,NY,NZ))

      TAU_INNER=TAU/FLOAT(NITER_TRANS)
      FO=FF
	FP=FF
C   THE MAIN TRANSPORT PROCEDURE

C SOLVING TRANSPORT EQUATIONS WITH SEVERAL ITERATIONS
      DO ITER=1,NITER_TRANS
      
       FM=(1.0+ALPHA_TRANS)*FF-ALPHA_TRANS*FP
!$OMP PARALLEL DO PRIVATE(M,N,K,FX_P,FX_M,FY_P,FY_M,FZ_P,FZ_M,
!$OMP&              BP,BP0,RHS,SL_P,SL_M)      
      DO N=NNN,NN
        DO M=MMM,MM
          IF(LU(M,N).GT.0.5) THEN
           
            SL_P=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER)
     &           +SLH (M,N)*DFLOAT(ITER))/DFLOAT(NITER_TRANS)
            SL_M=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER+1)
     &           +SLH (M,N)*DFLOAT(ITER-1))/DFLOAT(NITER_TRANS)                       
            
            BP = (HHQ(M,N)-SNGL(SL_P)) / TAU_INNER
            BP0= (HHQ(M,N)-SNGL(SL_M)) / TAU_INNER
            
            DO K=1,NZ 
	       FX_P=DYH(M  ,N   ) *HHU(M  ,N  ) * LCU(M  ,N  )
     &          *( UU(M  ,N  ,K)*(FM(M+1,N  ,K) +FM(M  ,N  ,K))/2.0
     &           - AX(M  ,N  ,K)*(FO(M+1,N  ,K) -FO(M  ,N  ,K))
     &           /DXT(M  ,N  )/RN )

	       FX_M=DYH(M-1,N   ) *HHU(M-1,N  ) * LCU(M-1,N  )
     &          *( UU(M-1,N  ,K)*(FM(M  ,N  ,K) +FM(M-1,N  ,K))/2.0
     &           - AX(M-1,N  ,K)*(FO(M  ,N  ,K) -FO(M-1,N  ,K))
     &           /DXT(M-1,N  )/RN )

	       FY_P=DXH(M  ,N   ) *HHV(M  ,N  ) * LCV(M  ,N  )
     &          *( VV(M  ,N  ,K)*(FM(M  ,N+1,K) +FM(M  ,N  ,K))/2.0
     &           - AY(M  ,N  ,K)*(FO(M  ,N+1,K) -FO(M  ,N  ,K))
     &           /DYT(M  ,N  )/RN )

	       FY_M=DXH(M  ,N-1 ) *HHV(M  ,N-1) * LCV(M  ,N-1)
     &          *( VV(M  ,N-1,K)*(FM(M  ,N  ,K) +FM(M  ,N-1,K))/2.0
     &           - AY(M  ,N-1,K)*(FO(M  ,N  ,K) -FO(M  ,N-1,K))
     &           /DYT(M  ,N-1)/RN )

            IF(K.LT.NZ)THEN
	       FZ_P=WW(M  ,N  ,K+1)*(FM(M  ,N  ,K) +FM(M  ,N  ,K+1))/2.0
            ELSE
             FZ_P=0.0
	      END IF
            
            IF(K.GT.1 )THEN
	       FZ_M=WW(M  ,N  ,K  )*(FM(M  ,N  ,K) +FM(M  ,N  ,K-1))/2.0
            ELSE
             FZ_M=0.0
	      END IF
            
	      RHS=( FX_P-FX_M + FY_P-FY_M )/(DX(M,N)*DY(M,N)*RN)
     &                      +(FZ_P-FZ_M )/ DZ(K)

            FF(M,N,K)=(BP0*FO(M,N,K)-RHS)/BP

            END DO

          END IF
        END DO
      END DO
!$OMP END PARALLEL DO

      IF(MMD.NE.0)  THEN
      CALL CYCLIZE (FF,NX,NY,NZ,MMM,MM)
	END IF
      
      FP=FO
      FO=FF

      END DO
C   END OF ITERATIVE LOOP
      
      deallocate(FP,FM,FO )
      RETURN
	END

C==============================================================
C  TRANSPORT EXPLICIT SCHEME USING UNIVERSAL DIFFUSION AND ADAMS-BASHFORT TRANSPORT
      SUBROUTINE TRANTS_AD_BASH_UNIDIF(FF,TAU,UU,VV,WW,SLRX,SLRY,
     &                    AMX,AMY,AFNT,AFNV,SLH,SLH0,
     &                    AS,AZ,AR,BZ,BR,BZR)
                                           
	IMPLICIT NONE
      INCLUDE '0COM.INC'
	INCLUDE '1TRANTS.INC'


      REAL  FF(NX,NY,NZ),TAU,TAU_INNER
	REAL  UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ+1) ! TRANSPORTING VELOCITIES

      REAL FX_P,FX_M,FY_P,FY_M,FZ_P,FZ_M   !FLUXES THROUGH CELL EDGES
	REAL(8) SLH(NX,NY),SLH0(NX,NY),SL_P,SL_M !SEA SURFACE HEIGHTS FOR CONSERVATION    
      
	REAL BP,BP0,RHS_T
      INTEGER M,N,K,KFM1,KFP1,ITER

      REAL  AS,AZ,AR,BZ,BR,BZR
	REAL    SLRX(NX,NY,NZ  ),  AMX(NX,NY,NZ  ), 
     &        SLRY(NX,NY,NZ  ),  AMY(NX,NY,NZ  ),
     &        AFNT(NX,NY,NZ+1), AFNV(NX,NY,NZ+1)
           
      REAL,allocatable:: FO(:,:,:),RHS(:,:,:),FP(:,:,:)               !OLD VALUE OF FF
      REAL AMXX(NX,NZ+1),AMXZ(NX,NZ+1),AMZX(NX,NZ+1),AMZZ(NX,NZ+1),
     &     BMYY(NY,NZ+1),BMYZ(NY,NZ+1),BMZY(NY,NZ+1),BMZZ(NY,NZ+1)

      REAL AFUNT, AFUNV, AFUNTP1, AFUNVP1, DC1, DCX, DCY
      REAL AMUW,ALFR,ALFZ,BLFR,BLFZ,BLFZR,HLAMZ,HLAMP

      allocate(FO(NX,NY,NZ),RHS(NX,NY,NZ),FP(NX,NY,NZ))

      RHS=0.0
      TAU_INNER=TAU/FLOAT(NITER_TRANS)

C   THE MAIN TRANSPORT PROCEDURE

C SOLVING TRANSPORT EQUATIONS WITH SEVERAL ITERATIONS
      
	FP=FF
      
      DO ITER=1,NITER_TRANS
          
       FO=(1.0+ALPHA_TRANS)*FF-ALPHA_TRANS*FP
       FP=FF

!$OMP PARALLEL DO PRIVATE(M,N,K,FX_P,FX_M,FY_P,FY_M,FZ_P,FZ_M,
!$OMP&              BP,BP0,SL_P,SL_M,RHS_T)      
      DO N=NNN,NN
        DO M=MMM,MM
          IF(LU(M,N).GT.0.5) THEN
            
            SL_P=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER)
     &           +SLH (M,N)*DFLOAT(ITER))/DFLOAT(NITER_TRANS)
            SL_M=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER+1)
     &           +SLH (M,N)*DFLOAT(ITER-1))/DFLOAT(NITER_TRANS)                       
            
            BP = (HHQ(M,N)-SNGL(SL_P)) / TAU_INNER
            BP0= (HHQ(M,N)-SNGL(SL_M)) / TAU_INNER
            
            DO K=1,NZ 
	       FX_P=UU(M  ,N  ,K)*DYH(M  ,N  ) *HHU(M  ,N  )
     &          *(FO(M  ,N  ,K) +FO(M+1,N  ,K))/2.0

	       FX_M=UU(M-1,N  ,K)*DYH(M-1,N  ) *HHU(M-1,N  )
     &          *(FO(M  ,N  ,K) +FO(M-1,N  ,K))/2.0

	       FY_P=VV(M  ,N  ,K)*DXH(M  ,N  ) *HHV(M  ,N  )
     &          *(FO(M  ,N  ,K) +FO(M  ,N+1,K))/2.0

	       FY_M=VV(M  ,N-1,K)*DXH(M  ,N-1) *HHV(M  ,N-1)
     &          *(FO(M  ,N  ,K) +FO(M  ,N-1,K))/2.0

            IF(K.LT.NZ)THEN
	       FZ_P=WW(M  ,N  ,K+1)
     &          *(FO(M  ,N  ,K) +FO(M  ,N  ,K+1))/2.0
            ELSE
             FZ_P=0.0
	      END IF
            
            IF(K.GT.1 )THEN
	       FZ_M=WW(M  ,N  ,K  )
     &          *(FO(M  ,N  ,K) +FO(M  ,N  ,K-1))/2.0
            ELSE
             FZ_M=0.0
	      END IF
            
	      RHS_T=( FX_P-FX_M + FY_P-FY_M )/(DX(M,N)*DY(M,N)*RN)
     &                        +(FZ_P-FZ_M )/ DZ(K)

            FF(M,N,K)=(BP0*FF(M,N,K)-RHS_T)/BP

            END DO

          END IF
        END DO
      END DO
!$OMP END PARALLEL DO
      
      IF(MMD.NE.0)  THEN
      CALL CYCLIZE (FF,NX,NY,NZ,MMM,MM)
	END IF

      END DO
      
      FO=FF


C Z-DIFFUZION IN XZ-ZX DIRECTION
!$OMP PARALLEL DO PRIVATE(M,N,K,AMXX,AMXZ,AMZX,AMZZ,AFUNT,AFUNV,
!$OMP&    AFUNTP1,AFUNVP1,ALFR,AMUW,ALFZ,BLFR,BLFZ,BLFZR,
!$OMP&    HLAMZ,HLAMP,DC1,DCX,FX_P,FX_M,FZ_P,FZ_M,KFM1,KFP1)
      DO N=NNN,NN
      AMXX =0.0
      AMXZ =0.0
      AMZX =0.0
      AMZZ =0.0      
      
C CALCULATING EXCHANGE COEFFICIENTS   
         DO M=MMM-1,MM
         
            IF(LCU(M,N).GT.0.5) THEN
         
             DO K=2,NZ
               AFUNT = AFNT(M,N,K)
               AFUNV = AFNV(M,N,K)
             
               AFUNTP1 = AFNT(M+1,N,K)
               AFUNVP1 = AFNV(M+1,N,K)

               ALFR=0.5*(AFUNT+AFUNTP1)
C  COEF. OF DIFFUSION IN W-POINTS AS FUNCTION ON Z:
               AMUW=0.25*(AMX(M,N,K-1)+AMX(M,N,K))*(AFUNV+AFUNVP1)
     &                  *DYH(M,N)/DXT(M,N)/RN
               ALFZ=(1.0-ALFR)*AR+AZ
               BLFR=BR*ALFR               
               BLFZ=(1.0-ALFR)*(BR+BZR)+BZ
              BLFZR=ALFR*BZR               
               ALFR=ALFR*AR

	         HLAMZ=    (Z3D(M+1,N,K  )-Z3D(M,N,K  )
     &                   +Z3D(M+1,N,K-1)-Z3D(M,N,K-1))*0.5
               HLAMP=SLRX(M,N,K)                     !ISOPICNAL

C  DIFFUSION ALONG SIGMA LEVELS:
               AMXX(M,K)=AS*AMUW*HZT(K)*HHU(M,N)
C  ADDITIONAL TERM FOR DIFFUSION ALONG ISOPICNALS AND HORIZONTAL LEVELS:
               AMXZ(M,K)=AMUW*( ALFR*(HLAMP  -ALPHA_GM*(HLAMP-HLAMZ)) 
     &                                               +   ALFZ*HLAMZ  )
               AMZX(M,K)=AMUW*( ALFR*(HLAMP  +ALPHA_GM*(HLAMP-HLAMZ)) 
     &                                               +   ALFZ*HLAMZ  )
               
               AMZZ(M,K)=AMUW*(BLFZ*HLAMZ**2+BLFR*HLAMP**2+
     +                        BLFZR*HLAMZ   *     HLAMP)/
     /                        (HHU(M,N)*HZT(K))         
             END DO
         
            END IF        
         
         END DO

C CALCULATING DIFFUZIONAL FLUXES THROUGH X- ABD Z- EDGES         
         DO K=1,NZ
            
            KFM1=MAX(1,K-1)
            KFP1=MIN(NZ,K+1)            
            
            IF(K.EQ.1.OR.K.EQ.NZ)  THEN
               DC1 = 0.5/DZ(K)
            ELSE
               DC1 = 0.25/DZ(K)
            END IF           
          
          DO M=MMM,MM
            IF(LU(M,N).GT.0.5) THEN

            IF(LCU(M-1,N).LT.0.5.OR.LCU(M,N).LT.0.5) THEN
               DCX = 0.5/DZ(K)
            ELSE
               DCX = 0.25/DZ(K)
            END IF           
C       FLUX X+
	     FX_P= AMXX(M  ,K   )*(FO(M+1,N,K   )-FO(M  ,N,K   )
     &                          +FO(M+1,N,KFM1)-FO(M  ,N,KFM1))
     &          +AMXX(M  ,K+1 )*(FO(M+1,N,K   )-FO(M  ,N,K   )
     &                          +FO(M+1,N,KFP1)-FO(M  ,N,KFP1))
     &          -AMXZ(M,  K   )*(FO(M  ,N,K   )-FO(M  ,N,KFM1)
     &                          +FO(M+1,N,K   )-FO(M+1,N,KFM1))
     &          -AMXZ(M,  K+1 )*(FO(M  ,N,KFP1)-FO(M  ,N,K   )
     &                          +FO(M+1,N,KFP1)-FO(M+1,N,K   ))
C       FLUX X-
           FX_M= AMXX(M-1,K   )*(FO(M  ,N,K   )-FO(M-1,N,K   )
     &                          +FO(M  ,N,KFM1)-FO(M-1,N,KFM1))
     &          +AMXX(M-1,K+1 )*(FO(M  ,N,K   )-FO(M-1,N,K  )
     &                          +FO(M  ,N,KFP1)-FO(M-1,N,KFP1))
     &          -AMXZ(M-1,K   )*(FO(M  ,N,K   )-FO(M  ,N,KFM1)
     &                          +FO(M-1,N,K   )-FO(M-1,N,KFM1))
     &          -AMXZ(M-1,K+1 )*(FO(M  ,N,KFP1)-FO(M  ,N,K   )
     &                          +FO(M-1,N,KFP1)-FO(M-1,N,K   ))
C       FLUX Z+          
           FZ_P=-AMZX(M,  K+1 )*(FO(M+1,N,KFP1)-FO(M  ,N,KFP1)
     &                          +FO(M+1,N,K   )-FO(M  ,N,K   ))
     &          -AMZX(M-1,K+1 )*(FO(M  ,N,KFP1)-FO(M-1,N,KFP1)
     &                          +FO(M  ,N,K   )-FO(M-1,N,K   ))
     &          +AMZZ(M,  K+1 )*(FO(M  ,N,KFP1)-FO(M  ,N,K   )
     &                          +FO(M+1,N,KFP1)-FO(M+1,N,K   ))
     &          +AMZZ(M-1,K+1 )*(FO(M  ,N,KFP1)-FO(M  ,N,K   )
     &                          +FO(M-1,N,KFP1)-FO(M-1,N,K   ))
C       FLUX Z- 
           FZ_M=-AMZX(M,  K   )*(FO(M+1,N,K   )-FO(M  ,N,K   )
     &                          +FO(M+1,N,KFM1)-FO(M  ,N,KFM1))
     &          -AMZX(M-1,K   )*(FO(M  ,N,K   )-FO(M-1,N,K   )
     &                          +FO(M  ,N,KFM1)-FO(M-1,N,KFM1))
     &          +AMZZ(M,  K   )*(FO(M  ,N,K   )-FO(M  ,N,KFM1)
     &                          +FO(M+1,N,K   )-FO(M+1,N,KFM1))
     &          +AMZZ(M-1,K   )*(FO(M  ,N,K   )-FO(M  ,N,KFM1)
     &                          +FO(M-1,N,K   )-FO(M-1,N,KFM1))
C    ADDING FLUXES IN THE EQUATION RIGHT HAND SIDE
            RHS(M,N,K)=RHS(M,N,K)-((FX_P-FX_M)*DC1+(FZ_P-FZ_M)*DCX)
     &                           /(DX(M,N)*DY(M,N)*RN)
            END IF 
          END DO
         END DO
      
      
      END DO
!$OMP END PARALLEL DO


C Z-DIFFUZION IN YZ-ZY DIRECTION
!$OMP PARALLEL DO PRIVATE(M,N,K,BMYY,BMYZ,BMZY,BMZZ,AFUNT,AFUNV,
!$OMP&    AFUNTP1,AFUNVP1,ALFR,AMUW,ALFZ,BLFR,BLFZ,BLFZR,
!$OMP&    HLAMZ,HLAMP,DC1,DCY,FY_P,FY_M,FZ_P,FZ_M,KFM1,KFP1)      
      DO M=MMM,MM
      BMYY =0.0
      BMYZ =0.0
      BMZY =0.0
      BMZZ =0.0      
      
C CALCULATING EXCHANGE COEFFICIENTS    
         DO N=NNN-1,NN
         
            IF(LCV(M,N).GT.0.5) THEN
         
             DO K=2,NZ
               AFUNT = AFNT(M,N,K)
               AFUNV = AFNV(M,N,K)
             
               AFUNTP1 = AFNT(M,N+1,K)
               AFUNVP1 = AFNV(M,N+1,K)

               ALFR=0.5*(AFUNT+AFUNTP1)
C  COEF. OF DIFFUSION IN W-POINTS AS FUNCTION ON Z:
               AMUW=0.25*(AMY(M,N,K-1)+AMY(M,N,K))*(AFUNV+AFUNVP1)
     &                  *DXH(M,N)/DYT(M,N)/RN
               ALFZ=(1.0-ALFR)*AR+AZ
               BLFR=BR*ALFR               
               BLFZ=(1.0-ALFR)*(BR+BZR)+BZ
              BLFZR=ALFR*BZR               
               ALFR=ALFR*AR

	         HLAMZ=    (Z3D(M,N+1,K  )-Z3D(M,N,K  )
     &                   +Z3D(M,N+1,K-1)-Z3D(M,N,K-1))*0.5
               HLAMP=SLRY(M,N,K)                     !ISOPICNAL

C  DIFFUSION ALONG SIGMA LEVELS:
               BMYY(N,K)=AS*AMUW*HZT(K)*HHV(M,N)
C  ADDITIONAL TERM FOR DIFFUSION ALONG ISOPICNALS AND HORIZONTAL LEVELS:
               BMYZ(N,K)=AMUW*( ALFR*(HLAMP  -ALPHA_GM*(HLAMP-HLAMZ)) 
     &                                               +   ALFZ*HLAMZ  )

               BMZY(N,K)=AMUW*( ALFR*(HLAMP  +ALPHA_GM*(HLAMP-HLAMZ)) 
     &                                               +   ALFZ*HLAMZ  )
               
               BMZZ(N,K)=AMUW*(BLFZ*HLAMZ**2+BLFR*HLAMP**2+
     +                        BLFZR*HLAMZ   *     HLAMP)/
     /                        (HHV(M,N)*HZT(K))         
             END DO
         
            END IF        
         
         END DO

C CALCULATING DIFFUZIONAL FLUXES THROUGH Y- ABD Z- EDGES         
         DO K=1,NZ
            
            IF(K.EQ.1.OR.K.EQ.NZ)  THEN
               DC1 = 0.5/DZ(K)
            ELSE
               DC1 = 0.25/DZ(K)
            END IF

            KFM1=MAX(1,K-1)
            KFP1=MIN(NZ,K+1)           
          
          DO N=NNN,NN
            IF(LU(M,N).GT.0.5) THEN
           
            IF(LCV(M,N-1).LT.0.5.OR.LCV(M,N).LT.0.5) THEN
               DCY = 0.5/DZ(K)
            ELSE
               DCY = 0.25/DZ(K)
            END IF
C       FLUX Y+
	     FY_P= BMYY(N  ,K   )*(FO(M,N+1,K   )-FO(M,N  ,K   )
     &                          +FO(M,N+1,KFM1)-FO(M,N  ,KFM1))
     &          +BMYY(N  ,K+1 )*(FO(M,N+1,K   )-FO(M,N  ,K   )
     &                          +FO(M,N+1,KFP1)-FO(M,N  ,KFP1))
     &          -BMYZ(N,  K   )*(FO(M,N  ,K   )-FO(M,N  ,KFM1)
     &                          +FO(M,N+1,K   )-FO(M,N+1,KFM1))
     &          -BMYZ(N,  K+1 )*(FO(M,N  ,KFP1)-FO(M,N  ,K   )
     &                          +FO(M,N+1,KFP1)-FO(M,N+1,K   ))
C       FLUX Y-
           FY_M= BMYY(N-1,K   )*(FO(M,N  ,K   )-FO(M,N-1,K   )
     &                          +FO(M,N  ,KFM1)-FO(M,N-1,KFM1))
     &          +BMYY(N-1,K+1 )*(FO(M,N  ,K   )-FO(M,N-1,K  )
     &                          +FO(M,N  ,KFP1)-FO(M,N-1,KFP1))
     &          -BMYZ(N-1,K   )*(FO(M,N  ,K   )-FO(M,N  ,KFM1)
     &                          +FO(M,N-1,K   )-FO(M,N-1,KFM1))
     &          -BMYZ(N-1,K+1 )*(FO(M,N  ,KFP1)-FO(M,N  ,K   )
     &                          +FO(M,N-1,KFP1)-FO(M,N-1,K   ))
C       FLUX Z+            
           FZ_P=-BMZY(N,  K+1 )*(FO(M,N+1,KFP1)-FO(M,N  ,KFP1)
     &                          +FO(M,N+1,K   )-FO(M,N  ,K   ))
     &          -BMZY(N-1,K+1 )*(FO(M,N  ,KFP1)-FO(M,N-1,KFP1)
     &                          +FO(M,N  ,K   )-FO(M,N-1,K   ))
     &          +BMZZ(N,  K+1 )*(FO(M,N  ,KFP1)-FO(M,N  ,K   )
     &                          +FO(M,N+1,KFP1)-FO(M,N+1,K   ))
     &          +BMZZ(N-1,K+1 )*(FO(M,N  ,KFP1)-FO(M,N  ,K   )
     &                          +FO(M,N-1,KFP1)-FO(M,N-1,K   ))
C       FLUX Z-
           FZ_M=-BMZY(N,  K   )*(FO(M,N+1,K   )-FO(M,N  ,K   )
     &                          +FO(M,N+1,KFM1)-FO(M,N  ,KFM1))
     &          -BMZY(N-1,K   )*(FO(M,N  ,K   )-FO(M,N-1,K   )
     &                          +FO(M,N  ,KFM1)-FO(M,N-1,KFM1))
     &          +BMZZ(N,  K   )*(FO(M,N  ,K   )-FO(M,N  ,KFM1)
     &                          +FO(M,N+1,K   )-FO(M,N+1,KFM1))
     &          +BMZZ(N-1,K   )*(FO(M,N  ,K   )-FO(M,N  ,KFM1)
     &                          +FO(M,N-1,K   )-FO(M,N-1,KFM1))
C    ADDING FLUXES IN THE EQUATION RIGHT HAND SIDE
            RHS(M,N,K)=RHS(M,N,K)-((FY_P-FY_M)*DC1+(FZ_P-FZ_M)*DCY)
     &                           /(DX(M,N)*DY(M,N)*RN)
            END IF 
          END DO
         END DO
      
      
      END DO
!$OMP END PARALLEL DO

C   THE MAIN DIFFUSION PROCEDURE


!$OMP PARALLEL DO PRIVATE(M,N,K,BP,BP0)      
      DO N=NNN,NN
        DO M=MMM,MM
          IF(LU(M,N).GT.0.5) THEN
            
            BP = (HHQ(M,N)-SNGL(SLH(M,N))) / TAU
            BP0= (HHQ(M,N)-SNGL(SLH(M,N))) / TAU
            
            DO K=1,NZ 

            FF(M,N,K)=(BP0*FO(M,N,K)-RHS(M,N,K))/BP

            END DO

          END IF
        END DO
      END DO
!$OMP END PARALLEL DO
      
      IF(MMD.NE.0)  THEN
      CALL CYCLIZE (FF,NX,NY,NZ,MMM,MM)
	END IF
      
c	FP=FO

      deallocate(FP,RHS,FO )
      RETURN
	END
C======================================================================
      SUBROUTINE DIFF_Z_IMPL(FF,FF0,A_Z,TAU,IGRZ,SLH)
      IMPLICIT NONE
C----------------------------------------------------------------------
C Z-TRANSPORT AND 3D VERTICAL DIFFUSION..
C ADVECTION IS CRANK-NICOLSON, DIFFUSION IS PURE IMPLICIT.
C----------------------------------------------------------------------
      INCLUDE '0COM.INC'

      REAL FF(NX,NY,NZ),FF0(NX,NY),A_Z(NX,NY,NZ+1)
      REAL A(NZ),B(NZ),C(NZ),ETA(NZ),RKSI(NZ)
      INTEGER IGRZ(NX,NY) !TYPE OF SEA SURFACE BOUNDARY CONDITION (1/2)
      INTEGER M, N, K
      REAL    BP, DP, DM 
      REAL    TAU
	REAL(8) SLH(NX,NY)
!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,BP,DP,DM,A,B,C,ETA,RKSI)
      DO N=NNN,NN
       DO M=MMM,MM
        IF (LU(M,N).GT.0.5) THEN
         BP = (HHQ(M,N)-SNGL(SLH (M,N))) / TAU
C INTERNAL POINTS.
         DO K=2,NZ-1
          DP = A_Z(M,N,K+1) / HHQ(M,N) / HZT(K+1) / DZ(K)
          DM = A_Z(M,N,K)   / HHQ(M,N) / HZT(K  ) / DZ(K)
C "2." ABOVE IS DUE TO APROXIMATION.
          C(K) =  - DP
          A(K) =  - DM
          B(K) =  BP + DP + DM
          ETA(K) = BP*FF(M,N,K) 
         ENDDO
C IMPLICIT: A(K)=-DM-PM;C(K)=-DP+PP;B(K)=BP+DP+DM;ETA(K)=BP*FF(M,N,K)
C SURFACE POINT: TWO CASES.
         K = 1
         DP = A_Z(M,N,K+1)/ HHQ(M,N) / HZT(K+1) / DZ(K)
         DM = A_Z(M,N,K)  / HHQ(M,N) / HZT(K  ) / DZ(K)

         A(K) = 0.0    !NOT USE IN FACTOR      
         C(K) = - DP
C AZ -> AZD, PM IS 0 DUE TO WW(M,N,1) = 0.
         IF(IGRZ(M,N).EQ.1) THEN
          B(K) =  BP + DP + DM
          ETA(K) = BP*FF(M,N,K) + DM * FF0(M,N)
         ELSEIF(IGRZ(M,N).EQ.2) THEN
          B(K) =  BP + DP
          ETA(K) = BP*FF(M,N,K) + FF0(M,N) / DZ(K)
         ENDIF
C BOTTOM POINT. FLUX IS ZERO FOR ALL FIELDS.
         K = NZ
         DM = A_Z(M,N,K) / HHQ(M,N) / HZT(K  ) / DZ(K)
         A(K) = - DM
         C(K)=0.0     !NOT USE IN FACTOR  
         B(K) =  BP + DM
         ETA(K) = BP* FF(M,N,K) 

         CALL FACTOR(NZ,A,B,C,ETA,RKSI,1,NZ)
         DO K=1,NZ
          FF(M,N,K)=RKSI(K)
         ENDDO
        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END

C-------- CALL ISOPYC(TT,TAU,2.E+6,0)----------------
C-------- CALL ISOPYC(SS,TAU,1.5+7,0)----------------
      SUBROUTINE ISOPYC(FF,TAU,COEF,IC)

	IMPLICIT NONE
      INCLUDE '0COM.INC'
      INCLUDE '0CEAN.INC'


      REAL  FF(NX,NY,NZ),TAU
      INTEGER M,N,K
      INTEGER M1,KKA,KK,N1,IC
      REAL,ALLOCATABLE:: F0(:,:,:),DEN1(:,:,:)
      REAL DD,FFA,AMN,AMN2,AMN3,COEFMAX,COEFMAX2,COEFMAX3,COEF1,
     * FFB,FLXF,FLXFT,CAMN,COEF,ANGLE,ANGLECR,ZN
      INTEGER KKMAX,KKMIN,NIT,NNIT
      REAL COEF0,COEF0MIN,TTMAX,TTMIN
      
      ALLOCATE(F0(NX,NY,NZ),DEN1(NX,NY,NZ))
      NNIT=6
      COEF0MIN=1.5E+6
      TTMIN=1.0
      TTMAX=5.0
      CAMN=0.15
      ANGLECR=0.008
!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO N=1,NY
         DO M=1,NX
            DO K=1,NZ
               F0(M,N,K)=0.0
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO N=NNN,NN
         DO M=MMM,MM
            IF(LU(M,N).GT.0.5) THEN
               DEN1(M,N,1)=DEN_P(M,N,1)
               DO K=2,NZ
                  DEN1(M,N,K)=AMAX1(DEN_P(M,N,K),DEN1(M,N,K-1))
               END DO
            END IF
         END DO
      END DO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(M,N,K,M1,N1,KKA,KK,DD,FFA,AMN,AMN2,AMN3,
!$OMP& COEFMAX,COEFMAX2,COEFMAX3,COEF0,COEF1,FFB,FLXF,FLXFT,ANGLE,ZN,
!$OMP& KKMAX,KKMIN,NIT)
      DO N=NNN,NN
         DO M=MMM,MM
            IF(LU(M,N).GT.0.5) THEN
               DO K=1,NZ
                  DD=DEN1(M,N,K)
                  FFA=FF(M,N,K)
                  COEF0=COEF
                  IF(IC.EQ.1) THEN
                     IF(TT(M,N,K).LT.TTMIN) THEN
                        COEF0=COEF0MIN
                     ELSE
                        IF(TT(M,N,K).GT.TTMAX) THEN
                           COEF0=COEF
                        ELSE
                           COEF0=COEF0MIN+(COEF-COEF0MIN)*
     *                     (TT(M,N,K)-TTMIN)/(TTMAX-TTMIN)
                        END IF
                     END IF
                  END IF
                  COEF0=COEF0/2.0
C-------------------M+1----------------------------
                  M1=M+1
                  IF(M1.GT.MM) THEN
                     M1=MMM
                  END IF
                  IF(LU(M1,N).GT.0.5) THEN
                     IF(DD.LT.DEN1(M1,N,1).OR.DD.GT.DEN1(M1,N,NZ)) THEN
                        KKA=0
                        AMN=1.0
                     ELSE
                        KKMIN=1
                        KKMAX=NZ
                        DO NIT=1,NNIT
                           KKA=(KKMAX+KKMIN)/2
                           IF(DD.GT.DEN1(M1,N,KKA)) THEN
                              KKMIN=KKA
                           ELSE
                              KKMAX=KKA
                           END IF
                        END DO
                        KKA=KKMIN
                        ZN=AMAX1(DEN1(M1,N,KKA+1)-DEN1(M1,N,KKA),3.E-7)
                        AMN=(DD-DEN1(M1,N,KKA))/ZN
C                        DO KK=1,NZ-1
C                           IF(DD.GE.DEN1(M1,N,KK).AND.
C     &                      DD.LE.DEN1(M1,N,KK+1)) THEN  
C                              KKA=KK
C                           ZN=AMAX1(DEN1(M1,N,KK+1)-DEN1(M1,N,KK),3.E-7)
C                              AMN=(DD-DEN1(M1,N,KK))/ZN
C                           END IF
C                        END DO
                     END IF
                     IF(KKA.NE.0) THEN
                        ANGLE=
     =                   ABS(HHQ(M1,N)*(Z(KKA)+AMN*(Z(KKA+1)-Z(KKA)))-
     -                       HHQ(M,N)*Z(K))/(DXT(M,N)*RN)
                        IF(ANGLE.LT.ANGLECR) THEN
                         COEFMAX=CAMN*RN*RN*DXT(M,N)*DX(M,N)/TAU
                         AMN2=(1.0-AMN)*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M1,N)*DY(M1,N)*DZ(KKA)*HHQ(M1,N))
                         AMN3=AMN*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M1,N)*DY(M1,N)*DZ(KKA+1)*HHQ(M1,N))
                         COEFMAX2=CAMN*COEFMAX*AMN2
                         COEFMAX3=CAMN*COEFMAX*AMN3
                         COEF1=AMIN1(COEF0,COEFMAX)
                         COEF1=AMIN1(COEF1,COEFMAX2)
                         COEF1=AMIN1(COEF1,COEFMAX3)
                         FFB=FF(M1,N,KKA)+
     +                    AMN*(FF(M1,N,KKA+1)-FF(M1,N,KKA))
                         FLXF=(FFB-FFA)*COEF1/(RN*DXT(M,N))
                         FLXFT=TAU*FLXF/(RN*DX(M,N))
                         F0(M,N,K)=F0(M,N,K)+FLXFT
                         F0(M1,N,KKA)=F0(M1,N,KKA)-FLXFT*AMN2
                         F0(M1,N,KKA+1)=F0(M1,N,KKA+1)-FLXFT*AMN3
                        END IF
                     END IF
                  END IF
C-------------------M-1----------------------------
                  M1=M-1
                  IF(M1.LT.MMM) THEN
                     M1=MM
                  END IF
                  IF(LU(M1,N).GT.0.5) THEN
                     IF(DD.LT.DEN1(M1,N,1).OR.DD.GT.DEN1(M1,N,NZ)) THEN
                        KKA=0
                        AMN=1.0
                     ELSE
                        KKMIN=1
                        KKMAX=NZ
                        DO NIT=1,NNIT
                           KKA=(KKMAX+KKMIN)/2
                           IF(DD.GT.DEN1(M1,N,KKA)) THEN
                              KKMIN=KKA
                           ELSE
                              KKMAX=KKA
                           END IF
                        END DO
                        KKA=KKMIN
                        ZN=AMAX1(DEN1(M1,N,KKA+1)-DEN1(M1,N,KKA),3.E-7)
                        AMN=(DD-DEN1(M1,N,KKA))/ZN
C                        DO KK=1,NZ-1
C                           IF(DD.GE.DEN1(M1,N,KK).AND.
C     &                      DD.LE.DEN1(M1,N,KK+1)) THEN  
C                              KKA=KK
C                           ZN=AMAX1(DEN1(M1,N,KK+1)-DEN1(M1,N,KK),3.E-7)
C                              AMN=(DD-DEN1(M1,N,KK))/ZN
C                           END IF
C                        END DO
                     END IF
                     IF(KKA.NE.0) THEN
                        ANGLE=
     =                   ABS(HHQ(M1,N)*(Z(KKA)+AMN*(Z(KKA+1)-Z(KKA)))-
     -                       HHQ(M,N)*Z(K))/(DXT(M1,N)*RN)
                        IF(ANGLE.LT.ANGLECR) THEN
                         COEFMAX=CAMN*RN*RN*DXT(M1,N)*DX(M,N)/TAU
                         AMN2=(1.0-AMN)*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M1,N)*DY(M1,N)*DZ(KKA)*HHQ(M1,N))
                         AMN3=AMN*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M1,N)*DY(M1,N)*DZ(KKA+1)*HHQ(M1,N))
                         COEFMAX2=CAMN*COEFMAX*AMN2
                         COEFMAX3=CAMN*COEFMAX*AMN3
                         COEF1=AMIN1(COEF0,COEFMAX)
                         COEF1=AMIN1(COEF1,COEFMAX2)
                         COEF1=AMIN1(COEF1,COEFMAX3)
                         FFB=FF(M1,N,KKA)+
     +                    AMN*(FF(M1,N,KKA+1)-FF(M1,N,KKA))
                         FLXF=(FFB-FFA)*COEF1/(RN*DXT(M1,N))
                         FLXFT=TAU*FLXF/(RN*DX(M,N))
                         F0(M,N,K)=F0(M,N,K)+FLXFT
                         F0(M1,N,KKA)=F0(M1,N,KKA)-FLXFT*AMN2
                         F0(M1,N,KKA+1)=F0(M1,N,KKA+1)-FLXFT*AMN3
                        END IF
                     END IF
                  END IF
C-------------------N+1----------------------------
                  N1=N+1
                  IF(LU(M,N1).GT.0.5) THEN
                     IF(DD.LT.DEN1(M,N1,1).OR.DD.GT.DEN1(M,N1,NZ)) THEN
                        KKA=0
                        AMN=1.0
                     ELSE
                        KKMIN=1
                        KKMAX=NZ
                        DO NIT=1,NNIT
                           KKA=(KKMAX+KKMIN)/2
                           IF(DD.GT.DEN1(M,N1,KKA)) THEN
                              KKMIN=KKA
                           ELSE
                              KKMAX=KKA
                           END IF
                        END DO
                        KKA=KKMIN
                        ZN=AMAX1(DEN1(M,N1,KKA+1)-DEN1(M,N1,KKA),3.E-7)
                        AMN=(DD-DEN1(M,N1,KKA))/ZN
C                        DO KK=1,NZ-1
C                           IF(DD.GE.DEN1(M,N1,KK).AND.
C     &                      DD.LE.DEN1(M,N1,KK+1)) THEN  
C                              KKA=KK
C                           ZN=AMAX1(DEN1(M,N1,KK+1)-DEN1(M,N1,KK),3.E-7)
C                              AMN=(DD-DEN1(M,N1,KK))/ZN
C                           END IF
C                        END DO
                     END IF
                     IF(KKA.NE.0) THEN
                        ANGLE=
     =                   ABS(HHQ(M,N1)*(Z(KKA)+AMN*(Z(KKA+1)-Z(KKA)))-
     -                       HHQ(M,N)*Z(K))/(DYT(M,N)*RN)
                        IF(ANGLE.LT.ANGLECR) THEN
                         COEFMAX=CAMN*RN*RN*DYT(M,N)*DY(M,N)/TAU
                         AMN2=(1.0-AMN)*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M,N1)*DY(M,N1)*DZ(KKA)*HHQ(M,N1))
                         AMN3=AMN*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M,N1)*DY(M,N1)*DZ(KKA+1)*HHQ(M,N1))
                         COEFMAX2=CAMN*COEFMAX*AMN2
                         COEFMAX3=CAMN*COEFMAX*AMN3
                         COEF1=AMIN1(COEF0,COEFMAX)
                         COEF1=AMIN1(COEF1,COEFMAX2)
                         COEF1=AMIN1(COEF1,COEFMAX3)
                         FFB=FF(M,N1,KKA)+
     +                    AMN*(FF(M,N1,KKA+1)-FF(M,N1,KKA))
                         FLXF=(FFB-FFA)*COEF1/(RN*DYT(M,N))
                         FLXFT=FLXF*TAU/(RN*DY(M,N))
                         F0(M,N,K)=F0(M,N,K)+FLXFT
                         F0(M,N1,KKA)=F0(M,N1,KKA)-FLXFT*AMN2
                         F0(M,N1,KKA+1)=F0(M,N1,KKA+1)-FLXFT*AMN3
                        END IF
                     END IF
                  END IF
C-------------------N-1----------------------------
                  N1=N-1
                  IF(LU(M,N1).GT.0.5) THEN
                     IF(DD.LT.DEN1(M,N1,1).OR.DD.GT.DEN1(M,N1,NZ)) THEN
                        KKA=0
                        AMN=1.0
                     ELSE
                        KKMIN=1
                        KKMAX=NZ
                        DO NIT=1,NNIT
                           KKA=(KKMAX+KKMIN)/2
                           IF(DD.GT.DEN1(M,N1,KKA)) THEN
                              KKMIN=KKA
                           ELSE
                              KKMAX=KKA
                           END IF
                        END DO
                        KKA=KKMIN
                        ZN=AMAX1(DEN1(M,N1,KKA+1)-DEN1(M,N1,KKA),3.E-7)
                        AMN=(DD-DEN1(M,N1,KKA))/ZN
C                        DO KK=1,NZ-1
C                           IF(DD.GE.DEN1(M,N1,KK).AND.
C     &                      DD.LE.DEN1(M,N1,KK+1)) THEN  
C                              KKA=KK
C                           ZN=AMAX1(DEN1(M,N1,KK+1)-DEN1(M,N1,KK),3.E-7)
C                              AMN=(DD-DEN1(M,N1,KK))/ZN
C                           END IF
C                        END DO
                     END IF
                     IF(KKA.NE.0) THEN
                        ANGLE=
     =                   ABS(HHQ(M,N1)*(Z(KKA)+AMN*(Z(KKA+1)-Z(KKA)))-
     -                       HHQ(M,N)*Z(K))/(DYT(M,N1)*RN)
                        IF(ANGLE.LT.ANGLECR) THEN
                         COEFMAX=CAMN*RN*RN*DYT(M,N1)*DY(M,N)/TAU
                         AMN2=(1.0-AMN)*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M,N1)*DY(M,N1)*DZ(KKA)*HHQ(M,N1))
                         AMN3=AMN*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M,N1)*DY(M,N1)*DZ(KKA+1)*HHQ(M,N1))
                         COEFMAX2=CAMN*COEFMAX*AMN2
                         COEFMAX3=CAMN*COEFMAX*AMN3
                         COEF1=AMIN1(COEF0,COEFMAX)
                         COEF1=AMIN1(COEF1,COEFMAX2)
                         COEF1=AMIN1(COEF1,COEFMAX3)
                         FFB=FF(M,N1,KKA)+
     +                    AMN*(FF(M,N1,KKA+1)-FF(M,N1,KKA))
                         FLXF=(FFB-FFA)*COEF1/(RN*DYT(M,N1))
                         FLXFT=FLXF*TAU/(RN*DY(M,N))
                         F0(M,N,K)=F0(M,N,K)+FLXFT
                         F0(M,N1,KKA)=F0(M,N1,KKA)-FLXFT*AMN2
                         F0(M,N1,KKA+1)=F0(M,N1,KKA+1)-FLXFT*AMN3
                        END IF
                     END IF
                  END IF
               END DO
            END IF
         END DO
      END DO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO N=NNN,NN
         DO M=MMM,MM
            IF(LU(M,N).GT.0.5) THEN
               DO K=1,NZ
                  FF(M,N,K)=FF(M,N,K)+F0(M,N,K)
               END DO
            END IF
         END DO
      END DO
!$OMP END PARALLEL DO
      CALL CYCLIZE (FF,NX,NY,NZ,MMM,MM)
      DEALLOCATE(F0,DEN1)
      RETURN
      END
C---------------call zdiff(TT,TAU,2.E+6)
      SUBROUTINE ZDIFF(FF,TAU,COEF)

	IMPLICIT NONE
      INCLUDE '0COM.INC'
      INCLUDE '0CEAN.INC'


      REAL  FF(NX,NY,NZ),TAU
      INTEGER M,N,K
      INTEGER M1,KKA,KK,N1
      REAL,ALLOCATABLE:: F0(:,:,:),DEN1(:,:,:)
      REAL DD,FFA,AMN,AMN2,AMN3,COEFMAX,COEFMAX2,COEFMAX3,COEF1,
     * FFB,FLXF,FLXFT,CAMN,COEF,ANGLE,ANGLECR,ZN
      INTEGER KKMAX,KKMIN,NIT,NNIT
      REAL COEF0,COEF0MIN,TTMAX,TTMIN
      
      ALLOCATE(F0(NX,NY,NZ),DEN1(NX,NY,NZ))
      NNIT=6
      CAMN=0.15
!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO N=1,NY
         DO M=1,NX
            DO K=1,NZ
               F0(M,N,K)=0.0
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO N=NNN,NN
         DO M=MMM,MM
            IF(LU(M,N).GT.0.5) THEN
               DO K=1,NZ
                  DEN1(M,N,K)=HHQ(M,N)*Z(K)
               END DO
            END IF
         END DO
      END DO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(M,N,K,M1,N1,KKA,KK,DD,FFA,AMN,AMN2,AMN3,
!$OMP& COEFMAX,COEFMAX2,COEFMAX3,COEF0,COEF1,FFB,FLXF,FLXFT,ZN,
!$OMP& KKMAX,KKMIN,NIT)
      DO N=NNN,NN
         DO M=MMM,MM
            IF(LU(M,N).GT.0.5) THEN
               DO K=1,NZ
                  DD=DEN1(M,N,K)
                  FFA=FF(M,N,K)
                  COEF0=COEF/2.0
C-------------------M+1----------------------------
                  M1=M+1
                  IF(M1.GT.MM) THEN
                     M1=MMM
                  END IF
                  IF(LU(M1,N).GT.0.5) THEN
                     IF(DD.LT.DEN1(M1,N,1).OR.DD.GT.DEN1(M1,N,NZ)) THEN
                        KKA=0
                        AMN=1.0
                     ELSE
                        KKMIN=1
                        KKMAX=NZ
                        DO NIT=1,NNIT
                           KKA=(KKMAX+KKMIN)/2
                           IF(DD.GT.DEN1(M1,N,KKA)) THEN
                              KKMIN=KKA
                           ELSE
                              KKMAX=KKA
                           END IF
                        END DO
                        KKA=KKMIN
                        ZN=AMAX1(DEN1(M1,N,KKA+1)-DEN1(M1,N,KKA),3.E-7)
                        AMN=(DD-DEN1(M1,N,KKA))/ZN
C                        DO KK=1,NZ-1
C                           IF(DD.GE.DEN1(M1,N,KK).AND.
C     &                      DD.LE.DEN1(M1,N,KK+1)) THEN  
C                              KKA=KK
C                           ZN=AMAX1(DEN1(M1,N,KK+1)-DEN1(M1,N,KK),3.E-7)
C                              AMN=(DD-DEN1(M1,N,KK))/ZN
C                           END IF
C                        END DO
                     END IF
                     IF(KKA.NE.0) THEN
                         COEFMAX=CAMN*RN*RN*DXT(M,N)*DX(M,N)/TAU
                         AMN2=(1.0-AMN)*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M1,N)*DY(M1,N)*DZ(KKA)*HHQ(M1,N))
                         AMN3=AMN*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M1,N)*DY(M1,N)*DZ(KKA+1)*HHQ(M1,N))
                         COEFMAX2=CAMN*COEFMAX*AMN2
                         COEFMAX3=CAMN*COEFMAX*AMN3
                         COEF1=AMIN1(COEF0,COEFMAX)
                         COEF1=AMIN1(COEF1,COEFMAX2)
                         COEF1=AMIN1(COEF1,COEFMAX3)
                         FFB=FF(M1,N,KKA)+
     +                    AMN*(FF(M1,N,KKA+1)-FF(M1,N,KKA))
                         FLXF=(FFB-FFA)*COEF1/(RN*DXT(M,N))
                         FLXFT=TAU*FLXF/(RN*DX(M,N))
                         F0(M,N,K)=F0(M,N,K)+FLXFT
                         F0(M1,N,KKA)=F0(M1,N,KKA)-FLXFT*AMN2
                         F0(M1,N,KKA+1)=F0(M1,N,KKA+1)-FLXFT*AMN3
                     END IF
                  END IF
C-------------------M-1----------------------------
                  M1=M-1
                  IF(M1.LT.MMM) THEN
                     M1=MM
                  END IF
                  IF(LU(M1,N).GT.0.5) THEN
                     IF(DD.LT.DEN1(M1,N,1).OR.DD.GT.DEN1(M1,N,NZ)) THEN
                        KKA=0
                        AMN=1.0
                     ELSE
                        KKMIN=1
                        KKMAX=NZ
                        DO NIT=1,NNIT
                           KKA=(KKMAX+KKMIN)/2
                           IF(DD.GT.DEN1(M1,N,KKA)) THEN
                              KKMIN=KKA
                           ELSE
                              KKMAX=KKA
                           END IF
                        END DO
                        KKA=KKMIN
                        ZN=AMAX1(DEN1(M1,N,KKA+1)-DEN1(M1,N,KKA),3.E-7)
                        AMN=(DD-DEN1(M1,N,KKA))/ZN
C                        DO KK=1,NZ-1
C                           IF(DD.GE.DEN1(M1,N,KK).AND.
C     &                      DD.LE.DEN1(M1,N,KK+1)) THEN  
C                              KKA=KK
C                           ZN=AMAX1(DEN1(M1,N,KK+1)-DEN1(M1,N,KK),3.E-7)
C                              AMN=(DD-DEN1(M1,N,KK))/ZN
C                           END IF
C                        END DO
                     END IF
                     IF(KKA.NE.0) THEN
                         COEFMAX=CAMN*RN*RN*DXT(M1,N)*DX(M,N)/TAU
                         AMN2=(1.0-AMN)*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M1,N)*DY(M1,N)*DZ(KKA)*HHQ(M1,N))
                         AMN3=AMN*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M1,N)*DY(M1,N)*DZ(KKA+1)*HHQ(M1,N))
                         COEFMAX2=CAMN*COEFMAX*AMN2
                         COEFMAX3=CAMN*COEFMAX*AMN3
                         COEF1=AMIN1(COEF0,COEFMAX)
                         COEF1=AMIN1(COEF1,COEFMAX2)
                         COEF1=AMIN1(COEF1,COEFMAX3)
                         FFB=FF(M1,N,KKA)+
     +                    AMN*(FF(M1,N,KKA+1)-FF(M1,N,KKA))
                         FLXF=(FFB-FFA)*COEF1/(RN*DXT(M1,N))
                         FLXFT=TAU*FLXF/(RN*DX(M,N))
                         F0(M,N,K)=F0(M,N,K)+FLXFT
                         F0(M1,N,KKA)=F0(M1,N,KKA)-FLXFT*AMN2
                         F0(M1,N,KKA+1)=F0(M1,N,KKA+1)-FLXFT*AMN3
                     END IF
                  END IF
C-------------------N+1----------------------------
                  N1=N+1
                  IF(LU(M,N1).GT.0.5) THEN
                     IF(DD.LT.DEN1(M,N1,1).OR.DD.GT.DEN1(M,N1,NZ)) THEN
                        KKA=0
                        AMN=1.0
                     ELSE
                        KKMIN=1
                        KKMAX=NZ
                        DO NIT=1,NNIT
                           KKA=(KKMAX+KKMIN)/2
                           IF(DD.GT.DEN1(M,N1,KKA)) THEN
                              KKMIN=KKA
                           ELSE
                              KKMAX=KKA
                           END IF
                        END DO
                        KKA=KKMIN
                        ZN=AMAX1(DEN1(M,N1,KKA+1)-DEN1(M,N1,KKA),3.E-7)
                        AMN=(DD-DEN1(M,N1,KKA))/ZN
C                        DO KK=1,NZ-1
C                           IF(DD.GE.DEN1(M,N1,KK).AND.
C     &                      DD.LE.DEN1(M,N1,KK+1)) THEN  
C                              KKA=KK
C                           ZN=AMAX1(DEN1(M,N1,KK+1)-DEN1(M,N1,KK),3.E-7)
C                              AMN=(DD-DEN1(M,N1,KK))/ZN
C                           END IF
C                        END DO
                     END IF
                     IF(KKA.NE.0) THEN
                         COEFMAX=CAMN*RN*RN*DYT(M,N)*DY(M,N)/TAU
                         AMN2=(1.0-AMN)*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M,N1)*DY(M,N1)*DZ(KKA)*HHQ(M,N1))
                         AMN3=AMN*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M,N1)*DY(M,N1)*DZ(KKA+1)*HHQ(M,N1))
                         COEFMAX2=CAMN*COEFMAX*AMN2
                         COEFMAX3=CAMN*COEFMAX*AMN3
                         COEF1=AMIN1(COEF0,COEFMAX)
                         COEF1=AMIN1(COEF1,COEFMAX2)
                         COEF1=AMIN1(COEF1,COEFMAX3)
                         FFB=FF(M,N1,KKA)+
     +                    AMN*(FF(M,N1,KKA+1)-FF(M,N1,KKA))
                         FLXF=(FFB-FFA)*COEF1/(RN*DYT(M,N))
                         FLXFT=FLXF*TAU/(RN*DY(M,N))
                         F0(M,N,K)=F0(M,N,K)+FLXFT
                         F0(M,N1,KKA)=F0(M,N1,KKA)-FLXFT*AMN2
                         F0(M,N1,KKA+1)=F0(M,N1,KKA+1)-FLXFT*AMN3
                     END IF
                  END IF
C-------------------N-1----------------------------
                  N1=N-1
                  IF(LU(M,N1).GT.0.5) THEN
                     IF(DD.LT.DEN1(M,N1,1).OR.DD.GT.DEN1(M,N1,NZ)) THEN
                        KKA=0
                        AMN=1.0
                     ELSE
                        KKMIN=1
                        KKMAX=NZ
                        DO NIT=1,NNIT
                           KKA=(KKMAX+KKMIN)/2
                           IF(DD.GT.DEN1(M,N1,KKA)) THEN
                              KKMIN=KKA
                           ELSE
                              KKMAX=KKA
                           END IF
                        END DO
                        KKA=KKMIN
                        ZN=AMAX1(DEN1(M,N1,KKA+1)-DEN1(M,N1,KKA),3.E-7)
                        AMN=(DD-DEN1(M,N1,KKA))/ZN
C                        DO KK=1,NZ-1
C                           IF(DD.GE.DEN1(M,N1,KK).AND.
C     &                      DD.LE.DEN1(M,N1,KK+1)) THEN  
C                              KKA=KK
C                           ZN=AMAX1(DEN1(M,N1,KK+1)-DEN1(M,N1,KK),3.E-7)
C                              AMN=(DD-DEN1(M,N1,KK))/ZN
C                           END IF
C                        END DO
                     END IF
                     IF(KKA.NE.0) THEN
                         COEFMAX=CAMN*RN*RN*DYT(M,N1)*DY(M,N)/TAU
                         AMN2=(1.0-AMN)*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M,N1)*DY(M,N1)*DZ(KKA)*HHQ(M,N1))
                         AMN3=AMN*DX(M,N)*DY(M,N)*DZ(K)*HHQ(M,N)/
     /                    (DX(M,N1)*DY(M,N1)*DZ(KKA+1)*HHQ(M,N1))
                         COEFMAX2=CAMN*COEFMAX*AMN2
                         COEFMAX3=CAMN*COEFMAX*AMN3
                         COEF1=AMIN1(COEF0,COEFMAX)
                         COEF1=AMIN1(COEF1,COEFMAX2)
                         COEF1=AMIN1(COEF1,COEFMAX3)
                         FFB=FF(M,N1,KKA)+
     +                    AMN*(FF(M,N1,KKA+1)-FF(M,N1,KKA))
                         FLXF=(FFB-FFA)*COEF1/(RN*DYT(M,N1))
                         FLXFT=FLXF*TAU/(RN*DY(M,N))
                         F0(M,N,K)=F0(M,N,K)+FLXFT
                         F0(M,N1,KKA)=F0(M,N1,KKA)-FLXFT*AMN2
                         F0(M,N1,KKA+1)=F0(M,N1,KKA+1)-FLXFT*AMN3
                     END IF
                  END IF
               END DO
            END IF
         END DO
      END DO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO N=NNN,NN
         DO M=MMM,MM
            IF(LU(M,N).GT.0.5) THEN
               DO K=1,NZ
                  FF(M,N,K)=FF(M,N,K)+F0(M,N,K)
               END DO
            END IF
         END DO
      END DO
!$OMP END PARALLEL DO
      CALL CYCLIZE (FF,NX,NY,NZ,MMM,MM)
      DEALLOCATE(F0,DEN1)
      RETURN
      END
C==============================================================
C  TRANSPORT EXPLICIT SCHEME without DIFFUSION
      SUBROUTINE TRANTS_AD_BASH_NODIFF(FF,TAU,UU,VV,WW,SLH,SLH0)
	IMPLICIT NONE
      INCLUDE '0COM.INC'
      INCLUDE '1TRANTS.INC'

      REAL  FF(NX,NY,NZ),TAU,TAU_INNER
	REAL  UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ+1) ! TRANSPORTING VELOCITIE

      REAL FX_P,FX_M,FY_P,FY_M,FZ_P,FZ_M   !FLUXES THROUGH CELL EDGES
	REAL(8) SLH(NX,NY),SLH0(NX,NY),SL_P,SL_M     !SEA SURFACE HEIGTS FOR CONSERVATION    
      
	REAL BP,BP0,RHS
      INTEGER M,N,K,ITER
           
      REAL,allocatable:: FO(:,:,:),FM(:,:,:),FP(:,:,:)            !OLD VALUE OF FF

      allocate(FO(NX,NY,NZ),FM(NX,NY,NZ),FP(NX,NY,NZ))

      TAU_INNER=TAU/FLOAT(NITER_TRANS)
      FO=FF
	FP=FF
C   THE MAIN TRANSPORT PROCEDURE

C SOLVING TRANSPORT EQUATIONS WITH SEVERAL ITERATIONS
      DO ITER=1,NITER_TRANS
      
       FM=(1.0+ALPHA_TRANS)*FO-ALPHA_TRANS*FP
!$OMP PARALLEL DO PRIVATE(M,N,K,FX_P,FX_M,FY_P,FY_M,FZ_P,FZ_M,
!$OMP&              BP,BP0,RHS,SL_P,SL_M)      
      DO N=NNN,NN
        DO M=MMM,MM
          IF(LU(M,N).GT.0.5) THEN
           
            SL_P=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER)
     &           +SLH (M,N)*DFLOAT(ITER))/DFLOAT(NITER_TRANS)
            SL_M=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER+1)
     &           +SLH (M,N)*DFLOAT(ITER-1))/DFLOAT(NITER_TRANS)                       
            
            BP = (HHQ(M,N)-SNGL(SL_P)) / TAU_INNER
            BP0= (HHQ(M,N)-SNGL(SL_M)) / TAU_INNER
            
            DO K=1,NZ 
	       FX_P=DYH(M  ,N   ) *HHU(M  ,N  ) * LCU(M  ,N  )
     &          *( UU(M  ,N  ,K)*(FM(M+1,N  ,K) +FM(M  ,N  ,K))/2.0 )

	       FX_M=DYH(M-1,N   ) *HHU(M-1,N  ) * LCU(M-1,N  )
     &          *( UU(M-1,N  ,K)*(FM(M  ,N  ,K) +FM(M-1,N  ,K))/2.0 )

	       FY_P=DXH(M  ,N   ) *HHV(M  ,N  ) * LCV(M  ,N  )
     &          *( VV(M  ,N  ,K)*(FM(M  ,N+1,K) +FM(M  ,N  ,K))/2.0 )

	       FY_M=DXH(M  ,N-1 ) *HHV(M  ,N-1) * LCV(M  ,N-1)
     &          *( VV(M  ,N-1,K)*(FM(M  ,N  ,K) +FM(M  ,N-1,K))/2.0 )

            IF(K.LT.NZ)THEN
	       FZ_P=WW(M  ,N  ,K+1)*(FM(M  ,N  ,K) +FM(M  ,N  ,K+1))/2.0
            ELSE
             FZ_P=0.0
	      END IF
            
            IF(K.GT.1 )THEN
	       FZ_M=WW(M  ,N  ,K  )*(FM(M  ,N  ,K) +FM(M  ,N  ,K-1))/2.0
            ELSE
             FZ_M=0.0
	      END IF
            
	      RHS=( FX_P-FX_M + FY_P-FY_M )/(DX(M,N)*DY(M,N)*RN)
     &                      +(FZ_P-FZ_M )/ DZ(K)

            FF(M,N,K)=(BP0*FO(M,N,K)-RHS)/BP

            END DO

          END IF
        END DO
      END DO
!$OMP END PARALLEL DO

      IF(MMD.NE.0)  THEN
      CALL CYCLIZE (FF,NX,NY,NZ,MMM,MM)
	END IF
      
      
      
      IF(ITER.EQ.1) THEN
	
	       FM=FF
!$OMP PARALLEL DO PRIVATE(M,N,K,FX_P,FX_M,FY_P,FY_M,FZ_P,FZ_M,
!$OMP&              BP,BP0,RHS,SL_P,SL_M)      
      DO N=NNN,NN
        DO M=MMM,MM
          IF(LU(M,N).GT.0.5) THEN
           
            SL_P=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER)
     &           +SLH (M,N)*DFLOAT(ITER))/DFLOAT(NITER_TRANS)
            SL_M=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER+1)
     &           +SLH (M,N)*DFLOAT(ITER-1))/DFLOAT(NITER_TRANS)                       
            
            BP = (HHQ(M,N)-SNGL(SL_P)) / TAU_INNER
            BP0= (HHQ(M,N)-SNGL(SL_M)) / TAU_INNER
            
            DO K=1,NZ 
	       FX_P=DYH(M  ,N   ) *HHU(M  ,N  ) * LCU(M  ,N  )
     &          *( UU(M  ,N  ,K)*(FM(M+1,N  ,K) +FM(M  ,N  ,K))/2.0 )

	       FX_M=DYH(M-1,N   ) *HHU(M-1,N  ) * LCU(M-1,N  )
     &          *( UU(M-1,N  ,K)*(FM(M  ,N  ,K) +FM(M-1,N  ,K))/2.0 )

	       FY_P=DXH(M  ,N   ) *HHV(M  ,N  ) * LCV(M  ,N  )
     &          *( VV(M  ,N  ,K)*(FM(M  ,N+1,K) +FM(M  ,N  ,K))/2.0 )

	       FY_M=DXH(M  ,N-1 ) *HHV(M  ,N-1) * LCV(M  ,N-1)
     &          *( VV(M  ,N-1,K)*(FM(M  ,N  ,K) +FM(M  ,N-1,K))/2.0 )

            IF(K.LT.NZ)THEN
	       FZ_P=WW(M  ,N  ,K+1)*(FM(M  ,N  ,K) +FM(M  ,N  ,K+1))/2.0
            ELSE
             FZ_P=0.0
	      END IF
            
            IF(K.GT.1 )THEN
	       FZ_M=WW(M  ,N  ,K  )*(FM(M  ,N  ,K) +FM(M  ,N  ,K-1))/2.0
            ELSE
             FZ_M=0.0
	      END IF
            
	      RHS=( FX_P-FX_M + FY_P-FY_M )/(DX(M,N)*DY(M,N)*RN)
     &                      +(FZ_P-FZ_M )/ DZ(K)

            FF(M,N,K)=(BP0*FO(M,N,K)-RHS)/BP

            END DO

          END IF
        END DO
      END DO
!$OMP END PARALLEL DO

      IF(MMD.NE.0)  THEN
      CALL CYCLIZE (FF,NX,NY,NZ,MMM,MM)
	END IF
      
      
      END IF


      FP=FO
      FO=FF

      END DO
C   END OF ITERATIVE LOOP
      
      deallocate(FP,FM,FO )
      RETURN
	END
C=====================================================================
      SUBROUTINE DENSITY_SLOPES(DEN,SLRX,SLRY)
      IMPLICIT NONE
C CALCULATING SLOPE FOR ISOPICNAL DIFFUSION
C DEN0 -- SURFACE DENSITY
C DEN  -- DENSITY
C (SLRX,SLRY) -- HH*(Dx(DENSITY),Dy(DENSITY))/(dDENSITY/dZ)
C              - ISOPICNAL SLOPES ON C-GRIDS
C EPSRHO -- MIN VERTICAL DENSITY GRADIENT (COMMONLY 1.0E-05)
C PERIODICITY IS TAKEN TO ACCOUNT BY MASKS LU, LCU AND LCV
C VERTICAL DENSITY GRADIENT IS AGREEMENT WITH THE OPERATOR OF
c ISOPICNAL DIFFUSION
      INCLUDE '0COM.INC'
C EPSRHO -- MIN VERTICAL DENSITY GRADIENT (COMMONLY 1.0E-12)
      REAL      EPSRHO, SLOPE,ANGLE,ANGLE_MAX,SLOPE_BOT
      PARAMETER(EPSRHO=1.0E-12,ANGLE_MAX=0.008)

      REAL DEN(NX,NY,NZ),SLRX(NX,NY,NZ),SLRY(NX,NY,NZ)
      INTEGER M,N,K
      REAL RHOZ,RHOX,RHOY
      

!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO N=NNN,NN
         DO M=MMM,MM
            IF(LU(M,N).GT.0.5) THEN
               DO K=2,NZ
                  DEN(M,N,K)=MAX(DEN(M,N,K),DEN(M,N,K-1))
               END DO
            END IF
         END DO
      END DO
!$OMP END PARALLEL DO

C DEFINE ARRAYS ON DEPTH LEVELS

!$OMP PARALLEL DO PRIVATE(M,N,K,RHOZ,RHOX,RHOY,SLOPE,ANGLE,SLOPE_BOT)
      DO K=2,NZ
         DO N=NNN,NN
            DO M=MMM,MM
               IF(LCU(M,N).GT.0.5) THEN
C VERTICAL DENSITY GRADIENT ON XC-W-GRID AS MEAN ON T-POINTS
C P-GRADIENT IS ADDED ACCORDING TO BRYDEN FORMULA (VARIANT 2)

      RHOZ=(DEN(M,N,K)-DEN(M,N,K-1)+DEN(M+1,N,K)-DEN(M+1,N,K-1))/
     &                                   (2.0*HHU(M,N)*HZT(K))

            RHOZ=MAX(RHOZ,EPSRHO)

            SLOPE=( DEN(M+1,N,K-1)-DEN(M,N,K-1)
     &                   +DEN(M+1,N,K  )-DEN(M,N,K  ) )/(2.0*RHOZ)

            SLOPE_BOT=( Z3D(M+1,N,K  )-Z3D(M,N,K  )
     &                     +Z3D(M+1,N,K-1)-Z3D(M,N,K-1) )/2.0

C  ISOPICNAL SLOPE ON XC-W-GRID AS MEAN ON T-POINTS
            ANGLE=( SLOPE- SLOPE_BOT )/(DXT(M,N)*RN)

            IF(ABS(ANGLE).LT.ANGLE_MAX) THEN
             SLRX(M,N,K)=SLOPE
            ELSE
             SLRX(M,N,K)=SLOPE_BOT+SIGN(ANGLE_MAX,ANGLE)*DXT(M,N)*RN
	      END IF
               
               END IF

               IF(LCV(M,N).GT.0.5) THEN
C  VERTICAL DENSITY GRADIENT ON YC-W-GRID AS MEAN ON T-POINTS
C P-GRADIENT IS ADDED ACCORDING TO BRYDEN FORMULA (VARIANT 2)

      RHOZ=(DEN(M,N,K)-DEN(M,N,K-1)+DEN(M,N+1,K)-DEN(M,N+1,K-1))/
     &                                   (2.0*HHV(M,N)*HZT(K))

            RHOZ=MAX(RHOZ,EPSRHO)

            SLOPE=( DEN(M,N+1,K-1)-DEN(M,N,K-1)
     &                   +DEN(M,N+1,K  )-DEN(M,N,K  ) )/(2.0*RHOZ)
            
            SLOPE_BOT=( Z3D(M,N+1,K  )-Z3D(M,N,K  )
     &                     +Z3D(M,N+1,K-1)-Z3D(M,N,K-1) )/2.0
            
            ANGLE=( SLOPE- SLOPE_BOT)/(DYT(M,N)*RN)

            IF(ABS(ANGLE).LT.ANGLE_MAX) THEN
             SLRY(M,N,K)=SLOPE
            ELSE
             SLRY(M,N,K)=SLOPE_BOT+SIGN(ANGLE_MAX,ANGLE)*DYT(M,N)*RN
	      END IF

               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
      
	IF(MMD.NE.0) THEN
	 CALL CYCLIZE(SLRX,NX,NY,NZ,MMM,MM)
	 CALL CYCLIZE(SLRY,NX,NY,NZ,MMM,MM)      
      END IF

      RETURN
      END
C==============================================================
C  TRANSPORT EXPLICIT SCHEME USING UNIVERSAL DIFFUSION AND ADAMS-BASHFORT TRANSPORT
      SUBROUTINE UNIVERSAL_DIFFUSION(FF,TAU,SLRX,SLRY,
     &                    AMX,AMY,AFNT,AFNV,
     &                    AS,AZ,AR,BZ,BR,BZR,SLH)
                                           
	IMPLICIT NONE
      INCLUDE '0COM.INC'
	INCLUDE '1TRANTS.INC'

      REAL  FF(NX,NY,NZ),TAU, BP
	REAL(8) SLH(NX,NY)

      REAL FX_P,FX_M,FY_P,FY_M,FZ_P,FZ_M   !FLUXES THROUGH CELL EDGES   
      
      INTEGER M,N,K,KFM1,KFP1,ITER

      REAL  AS,AZ,AR,BZ,BR,BZR
	REAL    SLRX(NX,NY,NZ  ),  AMX(NX,NY,NZ  ), 
     &        SLRY(NX,NY,NZ  ),  AMY(NX,NY,NZ  ),
     &        AFNT(NX,NY,NZ+1), AFNV(NX,NY,NZ+1)
           
      REAL,allocatable:: RHS(:,:,:)
      REAL AMXX(NX,NZ+1),AMXZ(NX,NZ+1),AMZX(NX,NZ+1),AMZZ(NX,NZ+1),
     &     BMYY(NY,NZ+1),BMYZ(NY,NZ+1),BMZY(NY,NZ+1),BMZZ(NY,NZ+1)

      REAL AFUNT, AFUNV, AFUNTP1, AFUNVP1, DC1, DCX, DCY
      REAL AMUW,ALFR,ALFZ,BLFR,BLFZ,BLFZR,HLAMZ,HLAMP

      allocate(RHS(NX,NY,NZ))

      RHS=0.0    
      

C Z-DIFFUZION IN XZ-ZX DIRECTION
!$OMP PARALLEL DO PRIVATE(M,N,K,AMXX,AMXZ,AMZX,AMZZ,AFUNT,AFUNV,
!$OMP&    AFUNTP1,AFUNVP1,ALFR,AMUW,ALFZ,BLFR,BLFZ,BLFZR,
!$OMP&    HLAMZ,HLAMP,DC1,DCX,FX_P,FX_M,FZ_P,FZ_M,KFM1,KFP1)
      DO N=NNN,NN
      AMXX =0.0
      AMXZ =0.0
      AMZX =0.0
      AMZZ =0.0      
      
C CALCULATING EXCHANGE COEFFICIENTS   
         DO M=MMM-1,MM
         
            IF(LCU(M,N).GT.0.5) THEN
         
             DO K=2,NZ
               AFUNT = AFNT(M,N,K)
               AFUNV = AFNV(M,N,K)
             
               AFUNTP1 = AFNT(M+1,N,K)
               AFUNVP1 = AFNV(M+1,N,K)

               ALFR=0.5*(AFUNT+AFUNTP1)
C  COEF. OF DIFFUSION IN W-POINTS AS FUNCTION ON Z:
               AMUW=0.25*(AMX(M,N,K-1)+AMX(M,N,K))*(AFUNV+AFUNVP1)
     &                  *DYH(M,N)/DXT(M,N)/RN
               ALFZ=(1.0-ALFR)*AR+AZ
               BLFR=BR*ALFR               
               BLFZ=(1.0-ALFR)*(BR+BZR)+BZ
              BLFZR=ALFR*BZR               
               ALFR=ALFR*AR

	         HLAMZ=    (Z3D(M+1,N,K  )-Z3D(M,N,K  )
     &                   +Z3D(M+1,N,K-1)-Z3D(M,N,K-1))*0.5
               HLAMP=SLRX(M,N,K)                     !ISOPICNAL

C  DIFFUSION ALONG SIGMA LEVELS:
               AMXX(M,K)=AS*AMUW*HZT(K)*HHU(M,N)
C  ADDITIONAL TERM FOR DIFFUSION ALONG ISOPICNALS AND HORIZONTAL LEVELS:
               AMXZ(M,K)=AMUW*( ALFR*(HLAMP  -ALPHA_GM*(HLAMP-HLAMZ)) 
     &                                               +   ALFZ*HLAMZ  )
               AMZX(M,K)=AMUW*( ALFR*(HLAMP  +ALPHA_GM*(HLAMP-HLAMZ)) 
     &                                               +   ALFZ*HLAMZ  )
               
               AMZZ(M,K)=AMUW*(BLFZ*HLAMZ**2+BLFR*HLAMP**2+
     +                        BLFZR*HLAMZ   *     HLAMP)/
     /                        (HHU(M,N)*HZT(K))         
             END DO
         
            END IF        
         
         END DO

C CALCULATING DIFFUZIONAL FLUXES THROUGH X- ABD Z- EDGES         
         DO K=1,NZ
            
            KFM1=MAX(1,K-1)
            KFP1=MIN(NZ,K+1)            
            
            IF(K.EQ.1.OR.K.EQ.NZ)  THEN
               DC1 = 0.5/DZ(K)
            ELSE
               DC1 = 0.25/DZ(K)
            END IF           
          
          DO M=MMM,MM
            IF(LU(M,N).GT.0.5) THEN

            IF(LCU(M-1,N).LT.0.5.OR.LCU(M,N).LT.0.5) THEN
               DCX = 0.5/DZ(K)
            ELSE
               DCX = 0.25/DZ(K)
            END IF           
C       FLUX X+
	     FX_P= AMXX(M  ,K   )*(FF(M+1,N,K   )-FF(M  ,N,K   )
     &                          +FF(M+1,N,KFM1)-FF(M  ,N,KFM1))
     &          +AMXX(M  ,K+1 )*(FF(M+1,N,K   )-FF(M  ,N,K   )
     &                          +FF(M+1,N,KFP1)-FF(M  ,N,KFP1))
     &          -AMXZ(M,  K   )*(FF(M  ,N,K   )-FF(M  ,N,KFM1)
     &                          +FF(M+1,N,K   )-FF(M+1,N,KFM1))
     &          -AMXZ(M,  K+1 )*(FF(M  ,N,KFP1)-FF(M  ,N,K   )
     &                          +FF(M+1,N,KFP1)-FF(M+1,N,K   ))
C       FLUX X-
           FX_M= AMXX(M-1,K   )*(FF(M  ,N,K   )-FF(M-1,N,K   )
     &                          +FF(M  ,N,KFM1)-FF(M-1,N,KFM1))
     &          +AMXX(M-1,K+1 )*(FF(M  ,N,K   )-FF(M-1,N,K  )
     &                          +FF(M  ,N,KFP1)-FF(M-1,N,KFP1))
     &          -AMXZ(M-1,K   )*(FF(M  ,N,K   )-FF(M  ,N,KFM1)
     &                          +FF(M-1,N,K   )-FF(M-1,N,KFM1))
     &          -AMXZ(M-1,K+1 )*(FF(M  ,N,KFP1)-FF(M  ,N,K   )
     &                          +FF(M-1,N,KFP1)-FF(M-1,N,K   ))
C       FLUX Z+          
           FZ_P=-AMZX(M,  K+1 )*(FF(M+1,N,KFP1)-FF(M  ,N,KFP1)
     &                          +FF(M+1,N,K   )-FF(M  ,N,K   ))
     &          -AMZX(M-1,K+1 )*(FF(M  ,N,KFP1)-FF(M-1,N,KFP1)
     &                          +FF(M  ,N,K   )-FF(M-1,N,K   ))
     &          +AMZZ(M,  K+1 )*(FF(M  ,N,KFP1)-FF(M  ,N,K   )
     &                          +FF(M+1,N,KFP1)-FF(M+1,N,K   ))
     &          +AMZZ(M-1,K+1 )*(FF(M  ,N,KFP1)-FF(M  ,N,K   )
     &                          +FF(M-1,N,KFP1)-FF(M-1,N,K   ))
C       FLUX Z- 
           FZ_M=-AMZX(M,  K   )*(FF(M+1,N,K   )-FF(M  ,N,K   )
     &                          +FF(M+1,N,KFM1)-FF(M  ,N,KFM1))
     &          -AMZX(M-1,K   )*(FF(M  ,N,K   )-FF(M-1,N,K   )
     &                          +FF(M  ,N,KFM1)-FF(M-1,N,KFM1))
     &          +AMZZ(M,  K   )*(FF(M  ,N,K   )-FF(M  ,N,KFM1)
     &                          +FF(M+1,N,K   )-FF(M+1,N,KFM1))
     &          +AMZZ(M-1,K   )*(FF(M  ,N,K   )-FF(M  ,N,KFM1)
     &                          +FF(M-1,N,K   )-FF(M-1,N,KFM1))
C    ADDING FLUXES IN THE EQUATION RIGHT HAND SIDE
            RHS(M,N,K)=RHS(M,N,K)-((FX_P-FX_M)*DC1+(FZ_P-FZ_M)*DCX)
     &                           /(DX(M,N)*DY(M,N)*RN)
            END IF 
          END DO
         END DO
      
      
      END DO
!$OMP END PARALLEL DO


C Z-DIFFUZION IN YZ-ZY DIRECTION
!$OMP PARALLEL DO PRIVATE(M,N,K,BMYY,BMYZ,BMZY,BMZZ,AFUNT,AFUNV,
!$OMP&    AFUNTP1,AFUNVP1,ALFR,AMUW,ALFZ,BLFR,BLFZ,BLFZR,
!$OMP&    HLAMZ,HLAMP,DC1,DCY,FY_P,FY_M,FZ_P,FZ_M,KFM1,KFP1)      
      DO M=MMM,MM
      BMYY =0.0
      BMYZ =0.0
      BMZY =0.0
      BMZZ =0.0      
      
C CALCULATING EXCHANGE COEFFICIENTS    
         DO N=NNN-1,NN
         
            IF(LCV(M,N).GT.0.5) THEN
         
             DO K=2,NZ
               AFUNT = AFNT(M,N,K)
               AFUNV = AFNV(M,N,K)
             
               AFUNTP1 = AFNT(M,N+1,K)
               AFUNVP1 = AFNV(M,N+1,K)

               ALFR=0.5*(AFUNT+AFUNTP1)
C  COEF. OF DIFFUSION IN W-POINTS AS FUNCTION ON Z:
               AMUW=0.25*(AMY(M,N,K-1)+AMY(M,N,K))*(AFUNV+AFUNVP1)
     &                  *DXH(M,N)/DYT(M,N)/RN
               ALFZ=(1.0-ALFR)*AR+AZ
               BLFR=BR*ALFR               
               BLFZ=(1.0-ALFR)*(BR+BZR)+BZ
              BLFZR=ALFR*BZR               
               ALFR=ALFR*AR

	         HLAMZ=    (Z3D(M,N+1,K  )-Z3D(M,N,K  )
     &                   +Z3D(M,N+1,K-1)-Z3D(M,N,K-1))*0.5
               HLAMP=SLRY(M,N,K)                     !ISOPICNAL

C  DIFFUSION ALONG SIGMA LEVELS:
               BMYY(N,K)=AS*AMUW*HZT(K)*HHV(M,N)
C  ADDITIONAL TERM FOR DIFFUSION ALONG ISOPICNALS AND HORIZONTAL LEVELS:
               BMYZ(N,K)=AMUW*( ALFR*(HLAMP  -ALPHA_GM*(HLAMP-HLAMZ)) 
     &                                               +   ALFZ*HLAMZ  )

               BMZY(N,K)=AMUW*( ALFR*(HLAMP  +ALPHA_GM*(HLAMP-HLAMZ)) 
     &                                               +   ALFZ*HLAMZ  )
               
               BMZZ(N,K)=AMUW*(BLFZ*HLAMZ**2+BLFR*HLAMP**2+
     +                        BLFZR*HLAMZ   *     HLAMP)/
     /                        (HHV(M,N)*HZT(K))         
             END DO
         
            END IF        
         
         END DO

C CALCULATING DIFFUZIONAL FLUXES THROUGH Y- ABD Z- EDGES         
         DO K=1,NZ
            
            IF(K.EQ.1.OR.K.EQ.NZ)  THEN
               DC1 = 0.5/DZ(K)
            ELSE
               DC1 = 0.25/DZ(K)
            END IF

            KFM1=MAX(1,K-1)
            KFP1=MIN(NZ,K+1)           
          
          DO N=NNN,NN
            IF(LU(M,N).GT.0.5) THEN
           
            IF(LCV(M,N-1).LT.0.5.OR.LCV(M,N).LT.0.5) THEN
               DCY = 0.5/DZ(K)
            ELSE
               DCY = 0.25/DZ(K)
            END IF
C       FLUX Y+
	     FY_P= BMYY(N  ,K   )*(FF(M,N+1,K   )-FF(M,N  ,K   )
     &                          +FF(M,N+1,KFM1)-FF(M,N  ,KFM1))
     &          +BMYY(N  ,K+1 )*(FF(M,N+1,K   )-FF(M,N  ,K   )
     &                          +FF(M,N+1,KFP1)-FF(M,N  ,KFP1))
     &          -BMYZ(N,  K   )*(FF(M,N  ,K   )-FF(M,N  ,KFM1)
     &                          +FF(M,N+1,K   )-FF(M,N+1,KFM1))
     &          -BMYZ(N,  K+1 )*(FF(M,N  ,KFP1)-FF(M,N  ,K   )
     &                          +FF(M,N+1,KFP1)-FF(M,N+1,K   ))
C       FLUX Y-
           FY_M= BMYY(N-1,K   )*(FF(M,N  ,K   )-FF(M,N-1,K   )
     &                          +FF(M,N  ,KFM1)-FF(M,N-1,KFM1))
     &          +BMYY(N-1,K+1 )*(FF(M,N  ,K   )-FF(M,N-1,K  )
     &                          +FF(M,N  ,KFP1)-FF(M,N-1,KFP1))
     &          -BMYZ(N-1,K   )*(FF(M,N  ,K   )-FF(M,N  ,KFM1)
     &                          +FF(M,N-1,K   )-FF(M,N-1,KFM1))
     &          -BMYZ(N-1,K+1 )*(FF(M,N  ,KFP1)-FF(M,N  ,K   )
     &                          +FF(M,N-1,KFP1)-FF(M,N-1,K   ))
C       FLUX Z+            
           FZ_P=-BMZY(N,  K+1 )*(FF(M,N+1,KFP1)-FF(M,N  ,KFP1)
     &                          +FF(M,N+1,K   )-FF(M,N  ,K   ))
     &          -BMZY(N-1,K+1 )*(FF(M,N  ,KFP1)-FF(M,N-1,KFP1)
     &                          +FF(M,N  ,K   )-FF(M,N-1,K   ))
     &          +BMZZ(N,  K+1 )*(FF(M,N  ,KFP1)-FF(M,N  ,K   )
     &                          +FF(M,N+1,KFP1)-FF(M,N+1,K   ))
     &          +BMZZ(N-1,K+1 )*(FF(M,N  ,KFP1)-FF(M,N  ,K   )
     &                          +FF(M,N-1,KFP1)-FF(M,N-1,K   ))
C       FLUX Z-
           FZ_M=-BMZY(N,  K   )*(FF(M,N+1,K   )-FF(M,N  ,K   )
     &                          +FF(M,N+1,KFM1)-FF(M,N  ,KFM1))
     &          -BMZY(N-1,K   )*(FF(M,N  ,K   )-FF(M,N-1,K   )
     &                          +FF(M,N  ,KFM1)-FF(M,N-1,KFM1))
     &          +BMZZ(N,  K   )*(FF(M,N  ,K   )-FF(M,N  ,KFM1)
     &                          +FF(M,N+1,K   )-FF(M,N+1,KFM1))
     &          +BMZZ(N-1,K   )*(FF(M,N  ,K   )-FF(M,N  ,KFM1)
     &                          +FF(M,N-1,K   )-FF(M,N-1,KFM1))
C    ADDING FLUXES IN THE EQUATION RIGHT HAND SIDE
            RHS(M,N,K)=RHS(M,N,K)-((FY_P-FY_M)*DC1+(FZ_P-FZ_M)*DCY)
     &                           /(DX(M,N)*DY(M,N)*RN)
            END IF 
          END DO
         END DO
      
      
      END DO
!$OMP END PARALLEL DO

C   THE MAIN DIFFUSION PROCEDURE


!$OMP PARALLEL DO PRIVATE(M,N,K,BP)      
      DO N=NNN,NN
        DO M=MMM,MM
          IF(LU(M,N).GT.0.5) THEN
            
            BP = (HHQ(M,N)-SNGL(SLH(M,N))) / TAU
            
            DO K=1,NZ 

            FF(M,N,K)=FF(M,N,K)-RHS(M,N,K)/BP

            END DO

          END IF
        END DO
      END DO
!$OMP END PARALLEL DO
      
      IF(MMD.NE.0)  THEN
      CALL CYCLIZE (FF,NX,NY,NZ,MMM,MM)
	END IF

      deallocate(RHS)
      RETURN
	END
C==============================================================
C  TRANSPORT CABARET SCHEME without DIFFUSION
      SUBROUTINE TRANTS_CABARET(FF,FF_0,TAU,UU,VV,WW,
     &                          FF_FLX,FF_FLY,FF_FLZ,SLH,SLH0)
	IMPLICIT NONE
      INCLUDE '0COM.INC'
      INCLUDE '1TRANTS.INC'

      REAL  FF(NX,NY,NZ),FF_0(NX,NY,NZ),TAU,TAU_INNER
	REAL  UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ+1) ! TRANSPORTING VELOCITIES
	REAL  FF_FLX(NX,NY,NZ),FF_FLY(NX,NY,NZ),FF_FLZ(NX,NY,NZ+1)

      REAL FX_P,FX_M,FY_P,FY_M,FZ_P,FZ_M   !FLUXES THROUGH CELL EDGES
	REAL(8) SLH(NX,NY),SLH0(NX,NY),SL_P,SL_M,SL_H     !SEA SURFACE HEIGTS FOR CONSERVATION    
      
	REAL BP,BP0,RHS, VEC_X,VEC_Y,VEC_Z,VEC_MIN,VEC_MAX
      INTEGER M,N,K,ITER
	INTEGER P,Q,R
           
      REAL,allocatable::  FO_FLX(:,:,:),FO_FLY(:,:,:),FO_FLZ(:,:,:),
     &          FO(:,:,:),F_CORR(:,:,:),G_CORR(:,:,:),H_CORR(:,:,:)

      allocate(FO_FLX(NX,NY,NZ),FO_FLY(NX,NY,NZ),FO_FLZ(NX,NY,NZ+1),
     &            FO(NX,NY,NZ))

      TAU_INNER=TAU/FLOAT(NITER_TRANS)
 
C   THE MAIN TRANSPORT PROCEDURE

C SOLVING TRANSPORT EQUATIONS WITH SEVERAL ITERATIONS
      DO ITER=1,NITER_TRANS

       FO=FF
      IF(ITER.GT.1) THEN
	 FF_0=FF
      END IF

!$OMP PARALLEL DO PRIVATE(M,N,K,FX_P,FX_M,FY_P,FY_M,FZ_P,FZ_M,
!$OMP&              BP,BP0,RHS,SL_P,SL_M,SL_H)      

! PREDICTOR STEP FOR HALF-STEP       
       DO N=NNN,NN
        DO M=MMM,MM
          IF(LU(M,N).GT.0.5) THEN
           
            SL_P=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER)
     &           +SLH (M,N)*DFLOAT(ITER))/DFLOAT(NITER_TRANS)
            SL_M=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER+1)
     &           +SLH (M,N)*DFLOAT(ITER-1))/DFLOAT(NITER_TRANS)
            SL_H=(SL_P+SL_M)/2.0
            
            BP = (HHQ(M,N)-SNGL(SL_H)) / TAU_INNER *2.0
            BP0= (HHQ(M,N)-SNGL(SL_M)) / TAU_INNER *2.0
            
            DO K=1,NZ 
	       FX_P=DYH(M  ,N   ) *    HHU(M  ,N  ) * LCU(M  ,N  )
     &          *( UU(M  ,N  ,K)* FF_FLX(M  ,N  ,K))

	       FX_M=DYH(M-1,N   ) *    HHU(M-1,N  ) * LCU(M-1,N  )
     &          *( UU(M-1,N  ,K)* FF_FLX(M-1,N  ,K))

	       FY_P=DXH(M  ,N   ) *HHV(M  ,N  ) * LCV(M  ,N  )
     &          *( VV(M  ,N  ,K)* FF_FLY(M  ,N  ,K) )

	       FY_M=DXH(M  ,N-1 ) *    HHV(M  ,N-1) * LCV(M  ,N-1)
     &          *( VV(M  ,N-1,K)* FF_FLY(M  ,N-1,K) )

            IF(K.LT.NZ)THEN
	       FZ_P=WW(M  ,N  ,K+1)*FF_FLZ(M  ,N  ,K+1) 
            ELSE
             FZ_P=0.0
	      END IF
            
            IF(K.GT.1 )THEN
	       FZ_M=WW(M  ,N  ,K  )*FF_FLZ(M  ,N  ,K  ) 
            ELSE
             FZ_M=0.0
	      END IF
            
	      RHS=( FX_P-FX_M + FY_P-FY_M )/(DX(M,N)*DY(M,N)*RN)
     &                      +(FZ_P-FZ_M )/ DZ(K)

            FF(M,N,K)=(BP0*FF(M,N,K)-RHS)/BP

            END DO

          END IF
        END DO
       END DO
!$OMP END PARALLEL DO

      IF(MMD.NE.0)  THEN
      CALL CYCLIZE (FF,NX,NY,NZ,MMM,MM)
	END IF
      
!STORING OLD VALUES FOR FLUX VARIABLES
	FO_FLX=FF_FLX
	FO_FLY=FF_FLY
	FO_FLZ=FF_FLZ

! UPDATING FLUX VARIABLES
!$OMP PARALLEL DO PRIVATE(M,N,K,P,Q,R)
       DO N=NNN,NN
        DO M=MMM,MM
          
	     IF(LCU(M,N).GT.0.5) THEN
            DO K=1,NZ
	       P=NINT(0.5*(SIGN(1.0,UU(M,N,K))-1.0))
C             FF_FLX(M,N,K)=2.0*FF(M-P,N,K)-FO_FLX(M-2*P-1,N,K)
             FF_FLX(M,N,K)=FF(M-P,N,K)
     &                  +(FF(M-P,N,K)-FO_FLX(M-2*P-1,N,K))
     &                  *DXT(M,N)/DXT(M-2*P-1,N)
	      END DO
           END IF

           IF(LCV(M,N).GT.0.5) THEN
            DO K=1,NZ
	       Q=NINT(0.5*(SIGN(1.0,VV(M,N,K))-1.0))
C             FF_FLY(M,N,K)=2.0*FF(M,N-Q,K)-FO_FLY(M,N-2*Q-1,K)
             FF_FLY(M,N,K)=FF(M,N-Q,K)
     &                  +(FF(M,N-Q,K)-FO_FLY(M,N-2*Q-1,K))
     &                  *DYT(M,N)/DYT(M,N-2*Q-1)                 
            END DO    
           END IF
          

          IF(LU(M,N).GT.0.5) THEN
           DO K=2,NZ
	      R=NINT(0.5*(SIGN(1.0,WW(M,N,K))-1.0))
C            FF_FLZ(M,N,K)=2.0*FF(M,N,K-R-1)-FO_FLZ(M,N,K-2*R-1)
            FF_FLZ(M,N,K)=FF(M,N,K-R-1)
     &                  +(FF(M,N,K-R-1)-FO_FLZ(M,N,K-2*R-1))
     &                  *HZT(K)/HZT(K-2*R-1)

           END DO
          END IF

        END DO
	 END DO
!$OMP END PARALLEL DO

       IF(MMD.NE.0) THEN
          CALL CYCLIZE(FF_FLX,NX,NY,NZ,MMM,MM)
          CALL CYCLIZE(FF_FLY,NX,NY,NZ,MMM,MM)
          CALL CYCLIZE(FF_FLZ,NX,NY,NZ+1,MMM,MM)       
       END IF

      IF(CABARET_MONOTON.GT.0) THEN  ! MONOTONIZER START
!MONOTONIZER - SOURCE TERM ESTIMATION
       ALLOCATE(F_CORR(NX,NY,NZ),G_CORR(NX,NY,NZ),H_CORR(NX,NY,NZ))

!$OMP PARALLEL DO PRIVATE(M,N,K)
       DO N=NNN,NN
        DO M=MMM,MM
          IF(LU(M,N).GT.0.5) THEN
	      DO K=1,NZ

	       F_CORR(M,N,K)=2.0*FF(M,N,K)-FO(M,N,K)-FF_0(M,N,K)
     &         +0.5*TAU_INNER*(UU(M,N,K)+UU(M-1,N,K))
     &              *(FO_FLX(M,N,K)-FO_FLX(M-1,N,K))/(RN*DX(M,N))
	       
             G_CORR(M,N,K)=2.0*FF(M,N,K)-FO(M,N,K)-FF_0(M,N,K)
     &         +0.5*TAU_INNER*(VV(M,N,K)+VV(M,N-1,K))
     &              *(FO_FLY(M,N,K)-FO_FLY(M,N-1,K))/(RN*DY(M,N))
	       
             H_CORR(M,N,K)=2.0*FF(M,N,K)-FO(M,N,K)-FF_0(M,N,K)
     &         +0.5*TAU_INNER*(WW(M,N,K+1)+WW(M,N,K))
     &              *(FO_FLZ(M,N,K+1)-FO_FLZ(M,N,K))/(HHQ(M,N)*DZ(K))  
                    
            END DO
          ELSE
	     F_CORR(M,N,1:NZ)=0.0
           G_CORR(M,N,1:NZ)=0.0
           H_CORR(M,N,1:NZ)=0.0
          END IF
        END DO
	 END DO
!$OMP END PARALLEL DO

       IF(MMD.NE.0) THEN
          CALL CYCLIZE(F_CORR,NX,NY,NZ,MMM,MM)
          CALL CYCLIZE(G_CORR,NX,NY,NZ,MMM,MM)
          CALL CYCLIZE(H_CORR,NX,NY,NZ,MMM,MM)       
       END IF

!MONOTONIZER - FLUX VARIABLE LIMITATION

!$OMP PARALLEL DO PRIVATE(M,N,K,P,Q,R,VEC_X,VEC_Y,VEC_Z,
!$OMP&                                  VEC_MIN,VEC_MAX)
       DO N=NNN,NN
        DO M=MMM,MM

          IF(LCU(M,N).GT.0.5) THEN
           DO K=1,NZ
	      P=NINT(0.5*(SIGN(1.0,UU(M,N,K))-1.0))
	      
            VEC_X=FO_FLX(M-2*P-1,N,K)+F_CORR(M-P,N,K)
	      VEC_Y=FO_FLX(M      ,N,K)+F_CORR(M-P,N,K)
	      VEC_Z=FF(M-P,N,K)+F_CORR(M-P,N,K)
            VEC_MIN=MIN(VEC_X,VEC_Y,VEC_Z)
            VEC_MAX=MAX(VEC_X,VEC_Y,VEC_Z)

	      FF_FLX(M,N,K)=MIN(MAX(FF_FLX(M,N,K),VEC_MIN),VEC_MAX)

	     END DO
          END IF

          IF(LCV(M,N).GT.0.5) THEN
           DO K=1,NZ
	      Q=NINT(0.5*(SIGN(1.0,VV(M,N,K))-1.0))

	      VEC_X=FO_FLY(M,N-2*Q-1,K)+G_CORR(M,N-Q,K)
	      VEC_Y=FO_FLY(M,N      ,K)+G_CORR(M,N-Q,K)
	      VEC_Z=FF(M,N-Q,K)+G_CORR(M,N-Q,K)
            VEC_MIN=MIN(VEC_X,VEC_Y,VEC_Z)
            VEC_MAX=MAX(VEC_X,VEC_Y,VEC_Z)

	      FF_FLY(M,N,K)=MIN(MAX(FF_FLY(M,N,K),VEC_MIN),VEC_MAX)
                
           END DO
          END IF

          IF(LU(M,N).GT.0.5) THEN
           DO K=2,NZ
	      R=NINT(0.5*(SIGN(1.0,WW(M,N,K))-1.0))
	      
            VEC_X=FO_FLZ(M,N,K-2*R-1)+H_CORR(M,N,K-R-1)
	      VEC_Y=FO_FLZ(M,N,K      )+H_CORR(M,N,K-R-1)
	      VEC_Z=FF(M,N,K-R-1)+H_CORR(M,N,K-R-1)
            VEC_MIN=MIN(VEC_X,VEC_Y,VEC_Z)
            VEC_MAX=MAX(VEC_X,VEC_Y,VEC_Z)

	      FF_FLZ(M,N,K)=MIN(MAX(FF_FLZ(M,N,K),VEC_MIN),VEC_MAX)

           END DO
          END IF

        END DO
	 END DO
!$OMP END PARALLEL DO

       IF(MMD.NE.0) THEN
          CALL CYCLIZE(FF_FLX,NX,NY,NZ,MMM,MM)
          CALL CYCLIZE(FF_FLY,NX,NY,NZ,MMM,MM)
          CALL CYCLIZE(FF_FLZ,NX,NY,NZ+1,MMM,MM)       
       END IF
       
       DEALLOCATE(H_CORR,G_CORR,F_CORR)

      END IF ! MONOTONIZER END

! CORRECTOR STEP

!$OMP PARALLEL DO PRIVATE(M,N,K,FX_P,FX_M,FY_P,FY_M,FZ_P,FZ_M,
!$OMP&              BP,BP0,RHS,SL_P,SL_M,SL_H)      
       DO N=NNN,NN
        DO M=MMM,MM
          IF(LU(M,N).GT.0.5) THEN
           
            SL_P=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER)
     &           +SLH (M,N)*DFLOAT(ITER))/DFLOAT(NITER_TRANS)
            SL_M=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER+1)
     &           +SLH (M,N)*DFLOAT(ITER-1))/DFLOAT(NITER_TRANS)
            SL_H=(SL_P+SL_M)/2.0
            
            BP = (HHQ(M,N)-SNGL(SL_P)) / TAU_INNER *2.0
            BP0= (HHQ(M,N)-SNGL(SL_H)) / TAU_INNER *2.0
            
            DO K=1,NZ 
	       FX_P=DYH(M  ,N   ) *    HHU(M  ,N  ) * LCU(M  ,N  )
     &          *( UU(M  ,N  ,K)* FF_FLX(M  ,N  ,K))

	       FX_M=DYH(M-1,N   ) *    HHU(M-1,N  ) * LCU(M-1,N  )
     &          *( UU(M-1,N  ,K)* FF_FLX(M-1,N  ,K))

	       FY_P=DXH(M  ,N   ) *HHV(M  ,N  ) * LCV(M  ,N  )
     &          *( VV(M  ,N  ,K)* FF_FLY(M  ,N  ,K) )

	       FY_M=DXH(M  ,N-1 ) *    HHV(M  ,N-1) * LCV(M  ,N-1)
     &          *( VV(M  ,N-1,K)* FF_FLY(M  ,N-1,K) )

            IF(K.LT.NZ)THEN
	       FZ_P=WW(M  ,N  ,K+1)*FF_FLZ(M  ,N  ,K+1) 
            ELSE
             FZ_P=0.0
	      END IF
            
            IF(K.GT.1 )THEN
	       FZ_M=WW(M  ,N  ,K  )*FF_FLZ(M  ,N  ,K  ) 
            ELSE
             FZ_M=0.0
	      END IF
            
	      RHS=( FX_P-FX_M + FY_P-FY_M )/(DX(M,N)*DY(M,N)*RN)
     &                      +(FZ_P-FZ_M )/ DZ(K)

            FF(M,N,K)=(BP0*FF(M,N,K)-RHS)/BP

            END DO

          END IF
        END DO
       END DO
!$OMP END PARALLEL DO

       IF(MMD.NE.0)  THEN
        CALL CYCLIZE (FF,NX,NY,NZ,MMM,MM)
	 END IF


       CALL CABARET_BOUNDARY_CORRECT(FF,FF_FLX,FF_FLY,FF_FLZ)


       IF(MMD.NE.0) THEN
          CALL CYCLIZE(FF_FLX,NX,NY,NZ,MMM,MM)
          CALL CYCLIZE(FF_FLY,NX,NY,NZ,MMM,MM)
          CALL CYCLIZE(FF_FLZ,NX,NY,NZ+1,MMM,MM)       
       END IF


      END DO
C   END OF ITERATIVE LOOP
      
      deallocate(FO,FO_FLZ,FO_FLY,FO_FLX)
      RETURN
	END
C=====================================================================================
      SUBROUTINE CABARET_BOUNDARY_CORRECT(FF,FF_FLX,FF_FLY,FF_FLZ)
	IMPLICIT NONE
	INCLUDE '0COM.INC'
	REAL FF(NX,NY,NZ),FF_FLX(NX,NY,NZ),
     &                  FF_FLY(NX,NY,NZ),FF_FLZ(NX,NY,NZ+1)

      INTEGER M,N,K
      
!$OMP PARALLEL DO PRIVATE(M,N,K)
       DO N=NNN,NN
        DO M=MMM,MM
          
          IF(LLU(M,N).GT.0.5.AND.LCU(M,N).LT.0.5) THEN
          
            DO K=1,NZ
             FF_FLX(M,N,K)=(FF(M,N,K)*LU(M,N)+FF(M+1,N,K)*LU(M+1,N))
     &                             /(LU(M,N)+LU(M+1,N))
	      END DO
	    
          END IF

          IF(LLV(M,N).GT.0.5.AND.LCV(M,N).LT.0.5) THEN
            
            DO K=1,NZ
             FF_FLY(M,N,K)=(FF(M,N,K)*LU(M,N)+FF(M,N+1,K)*LU(M,N+1))
     &                             /(LU(M,N)+LU(M,N+1))             
            END DO	     
          
          END IF

          IF(LU(M,N).GT.0.5) THEN
           FF_FLZ(M,N,1)=FF(M,N,1)
           FF_FLZ(M,N,NZ+1)=FF(M,N,NZ)
          END IF

        END DO
	 END DO
!$OMP END PARALLEL DO

       IF(MMD.NE.0) THEN
          CALL CYCLIZE(FF_FLX,NX,NY,NZ,MMM,MM)
          CALL CYCLIZE(FF_FLY,NX,NY,NZ,MMM,MM)
          CALL CYCLIZE(FF_FLZ,NX,NY,NZ+1,MMM,MM)       
       END IF


      RETURN
	END
C=====================================================================================
      SUBROUTINE CABARET_FLUX_VAR_CORRECT(FF,FF_0,FF_FLX,FF_FLY,FF_FLZ)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
      INTEGER M,N,K
      REAL FF(NX,NY,NZ),FF_0(NX,NY,NZ),
     & FF_FLX(NX,NY,NZ),FF_FLY(NX,NY,NZ),FF_FLZ(NX,NY,NZ+1)
      REAL DELTA_FM,DELTA_FP

!$OMP PARALLEL DO PRIVATE(M,N,K,DELTA_FM,DELTA_FP)	 
       DO N=NNN,NN
	  DO M=MMM,MM
	  
          IF(LCU(M,N).GT.0.5) THEN
	     DO K=1,NZ
	      DELTA_FM=FF(M  ,N,K)-FF_0(M  ,N,K)
	      DELTA_FP=FF(M+1,N,K)-FF_0(M+1,N,K)
            FF_FLX(M,N,K)=FF_FLX(M,N,K)
     &      +(DELTA_FM+DELTA_FP)/2.0
	     END DO
          END IF

          IF(LCV(M,N).GT.0.5) THEN
	     DO K=1,NZ
	      DELTA_FM=FF(M,N  ,K)-FF_0(M,N  ,K)
	      DELTA_FP=FF(M,N+1,K)-FF_0(M,N+1,K)
            FF_FLY(M,N,K)=FF_FLY(M,N,K)
     &      +(DELTA_FM+DELTA_FP)/2.0
	     END DO         
          END IF	   

         
         IF(LU(M,N).GT.0.5) THEN
	     DO K=2,NZ
	      DELTA_FM=FF(M,N,K-1)-FF_0(M,N,K-1)
	      DELTA_FP=FF(M,N,K  )-FF_0(M,N,K  )
            FF_FLZ(M,N,K)=FF_FLZ(M,N,K)
     &      +(DELTA_FM+DELTA_FP)/2.0
	     END DO        	   
         END IF                 
        
        END DO
       END DO
!$OMP END PARALLEL DO

	RETURN
	END
C===================================================================================
C  TRANSPORT EXPLICIT SCHEME without DIFFUSION
      SUBROUTINE TRANTS_MATSUNO_NODIFF(FF,TAU,UU,VV,WW,SLH,SLH0)
	IMPLICIT NONE
      INCLUDE '0COM.INC'
      INCLUDE '1TRANTS.INC'

      REAL  FF(NX,NY,NZ),TAU,TAU_INNER
	REAL  UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ+1) ! TRANSPORTING VELOCITIE

      REAL FX_P,FX_M,FY_P,FY_M,FZ_P,FZ_M   !FLUXES THROUGH CELL EDGES
	REAL(8) SLH(NX,NY),SLH0(NX,NY),SL_P,SL_M     !SEA SURFACE HEIGTS FOR CONSERVATION    
      
	REAL BP,BP0,RHS
      INTEGER M,N,K,ITER
           
      REAL,allocatable:: FO(:,:,:),FM(:,:,:)            !OLD VALUE OF FF

      allocate(FO(NX,NY,NZ),FM(NX,NY,NZ))

      TAU_INNER=TAU/FLOAT(NITER_TRANS)

C   THE MAIN TRANSPORT PROCEDURE

C SOLVING TRANSPORT EQUATIONS WITH SEVERAL ITERATIONS
      DO ITER=1,NITER_TRANS
       
       FO=FF      

!$OMP PARALLEL DO PRIVATE(M,N,K,FX_P,FX_M,FY_P,FY_M,FZ_P,FZ_M,
!$OMP&              BP,BP0,RHS,SL_P,SL_M)      
      DO N=NNN,NN
        DO M=MMM,MM
          IF(LU(M,N).GT.0.5) THEN
           
            SL_P=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER)
     &           +SLH (M,N)*DFLOAT(ITER))/DFLOAT(NITER_TRANS)
            SL_M=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER+1)
     &           +SLH (M,N)*DFLOAT(ITER-1))/DFLOAT(NITER_TRANS)                       
            
            BP = (HHQ(M,N)-SNGL(SL_P)) / TAU_INNER
            BP0= (HHQ(M,N)-SNGL(SL_M)) / TAU_INNER
            
            DO K=1,NZ 
	       FX_P=DYH(M  ,N   ) *HHU(M  ,N  ) * LCU(M  ,N  )
     &          *( UU(M  ,N  ,K)*(FO(M+1,N  ,K) +FO(M  ,N  ,K))/2.0 )

	       FX_M=DYH(M-1,N   ) *HHU(M-1,N  ) * LCU(M-1,N  )
     &          *( UU(M-1,N  ,K)*(FO(M  ,N  ,K) +FO(M-1,N  ,K))/2.0 )

	       FY_P=DXH(M  ,N   ) *HHV(M  ,N  ) * LCV(M  ,N  )
     &          *( VV(M  ,N  ,K)*(FO(M  ,N+1,K) +FO(M  ,N  ,K))/2.0 )

	       FY_M=DXH(M  ,N-1 ) *HHV(M  ,N-1) * LCV(M  ,N-1)
     &          *( VV(M  ,N-1,K)*(FO(M  ,N  ,K) +FO(M  ,N-1,K))/2.0 )

            IF(K.LT.NZ)THEN
	       FZ_P=WW(M  ,N  ,K+1)*(FO(M  ,N  ,K) +FO(M  ,N  ,K+1))/2.0
            ELSE
             FZ_P=0.0
	      END IF
            
            IF(K.GT.1 )THEN
	       FZ_M=WW(M  ,N  ,K  )*(FO(M  ,N  ,K) +FO(M  ,N  ,K-1))/2.0
            ELSE
             FZ_M=0.0
	      END IF
            
	      RHS=( FX_P-FX_M + FY_P-FY_M )/(DX(M,N)*DY(M,N)*RN)
     &                      +(FZ_P-FZ_M )/ DZ(K)

            FF(M,N,K)=(BP0*FO(M,N,K)-RHS)/BP

            END DO

          END IF
        END DO
      END DO
!$OMP END PARALLEL DO

      IF(MMD.NE.0)  THEN
      CALL CYCLIZE (FF,NX,NY,NZ,MMM,MM)
	END IF
      
       FM=FF
!$OMP PARALLEL DO PRIVATE(M,N,K,FX_P,FX_M,FY_P,FY_M,FZ_P,FZ_M,
!$OMP&              BP,BP0,RHS,SL_P,SL_M)      
      DO N=NNN,NN
        DO M=MMM,MM
          IF(LU(M,N).GT.0.5) THEN
           
            SL_P=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER)
     &           +SLH (M,N)*DFLOAT(ITER))/DFLOAT(NITER_TRANS)
            SL_M=(SLH0(M,N)*DFLOAT(NITER_TRANS-ITER+1)
     &           +SLH (M,N)*DFLOAT(ITER-1))/DFLOAT(NITER_TRANS)                       
            
            BP = (HHQ(M,N)-SNGL(SL_P)) / TAU_INNER
            BP0= (HHQ(M,N)-SNGL(SL_M)) / TAU_INNER
            
            DO K=1,NZ 
	       FX_P=DYH(M  ,N   ) *HHU(M  ,N  ) * LCU(M  ,N  )
     &          *( UU(M  ,N  ,K)*(FM(M+1,N  ,K) +FM(M  ,N  ,K))/2.0 )

	       FX_M=DYH(M-1,N   ) *HHU(M-1,N  ) * LCU(M-1,N  )
     &          *( UU(M-1,N  ,K)*(FM(M  ,N  ,K) +FM(M-1,N  ,K))/2.0 )

	       FY_P=DXH(M  ,N   ) *HHV(M  ,N  ) * LCV(M  ,N  )
     &          *( VV(M  ,N  ,K)*(FM(M  ,N+1,K) +FM(M  ,N  ,K))/2.0 )

	       FY_M=DXH(M  ,N-1 ) *HHV(M  ,N-1) * LCV(M  ,N-1)
     &          *( VV(M  ,N-1,K)*(FM(M  ,N  ,K) +FM(M  ,N-1,K))/2.0 )

            IF(K.LT.NZ)THEN
	       FZ_P=WW(M  ,N  ,K+1)*(FM(M  ,N  ,K) +FM(M  ,N  ,K+1))/2.0
            ELSE
             FZ_P=0.0
	      END IF
            
            IF(K.GT.1 )THEN
	       FZ_M=WW(M  ,N  ,K  )*(FM(M  ,N  ,K) +FM(M  ,N  ,K-1))/2.0
            ELSE
             FZ_M=0.0
	      END IF
            
	      RHS=( FX_P-FX_M + FY_P-FY_M )/(DX(M,N)*DY(M,N)*RN)
     &                      +(FZ_P-FZ_M )/ DZ(K)

            FF(M,N,K)=(BP0*FO(M,N,K)-RHS)/BP

            END DO

          END IF
        END DO
      END DO
!$OMP END PARALLEL DO

      IF(MMD.NE.0)  THEN
      CALL CYCLIZE (FF,NX,NY,NZ,MMM,MM)
	END IF
      

      END DO
C   END OF ITERATIVE LOOP
      
      deallocate(FM,FO)
      RETURN
	END