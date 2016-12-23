C======================================================================
C OCEAN GENERAL CIRCULATION MODEL. FOR PERIODIC CASE.
C (C) BAGNO A.V., 1997, 1998.
C------------------------------------------------------06-23-97 04:03pm
C ADVECTION IS CRANK-NICOLSON, DIFFUSION IS PURE IMPLICIT.
C======================================================================
      SUBROUTINE TRXU(FF,UU,AX,TAU)
      IMPLICIT NONE
C----------------------------------------------------------------------
C X-TRANSPORT AND DIFFUSION OF FF ON U-GRID.
C X-TRANPORT SUPPOSES BOUNDARY VELOCITIES & FI TO BE ZERO.
      INCLUDE '0COM.INC'
      INCLUDE '0TRUVS.INC'
      REAL TAU
      REAL FF(NX,NY,NZ),UU(NX,NY,NZ),AX(NX,NY,NZ)
      REAL A(NX),B(NX),C(NX),ETA(NX),RKSI(NX)
      INTEGER M, N, K, IG, IGRX
      INTEGER II, JJ, JJP, MLOOP
      INTEGER M1, M9
      REAL    BP, DP, DM, PM, PP

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,IG,II,JJ,JJP,MLOOP,M1,M9,
!$OMP&        BP,DP,DM,PP,PM,A,B,C,ETA,RKSI)
      DO K=1,NZ
       DO N=NNN,NN
        DO IG=1,ABS(LRX(N))
         II=IIX(IG,N)

         IF (LRX(N).GT.0) THEN
		JJ=JJX(IG,N)-1
		JJP=JJ
         IF (JJ.GT.MM) JJ = JJ - MMD
         ELSE
		JJP=MM
         END IF

         DO MLOOP=II,JJP
          M  = MLOOP
          M1 = MLOOP - 1
          M9 = MLOOP + 1

          IF  ( M.GT.MM) M  = M  - MMD
          IF  (M1.GT.MM) M1 = M1 - MMD
          IF  (M9.GT.MM) M9 = M9 - MMD

          BP = HHU(M,N)*DXT(M,N)*DYH(M,N)*RN / TAU

          DP = AX(M9,N,K)*DY(M9,N)/DX(M9,N)*HHQ(M9,N)/RN
          DM = AX(M ,N,K)*DY(M ,N)/DX(M ,N)*HHQ(M ,N)/RN

          PP = ( UU(M ,N,K)*HHU(M ,N)*DYH(M ,N) +
     &           UU(M9,N,K)*HHU(M9,N)*DYH(M9,N) ) /4.

          PM = ( UU(M1,N,K)*HHU(M1,N)*DYH(M1,N) +
     &           UU(M ,N,K)*HHU(M ,N)*DYH(M ,N) ) /4.

C LAST "2" ABOVE IS DUE TO APROXIMATION IN SPACE.
          C(M) =  PP * TRUV_IMP  - DP
          A(M) = -PM * TRUV_IMP  - DM
          B(M) =  BP + DP + DM
          ETA(M) = BP * FF(M,N,K) - PP * FF(M9,N,K) * TRUV_EXP
     &                            + PM * FF(M1,N,K) * TRUV_EXP
         ENDDO
C IMPLICIT: A(M)=-DM-PM;C(M)=-DP+PP;B(M)=BP+DM+DP;ETA(M)=BP*FF(M,N,K)
C IF (M=II) : PM = 0, A = -DM ; IF (M=JJ) : PP = 0, C = -DP.

         IF (LRX(N).GT.0) THEN

          CALL FACTOR2(NX,A,B,C,ETA,RKSI,II,JJP)

         ELSE

          CALL FACTORC(A,B,C,ETA,RKSI,MMM,MM)

         ENDIF

          DO MLOOP=II,JJP
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
      SUBROUTINE TRYU(FF,VV,AY,TAU,IGRY)
      IMPLICIT NONE
C------------------------------------------------------09-13-96 06:08pm
C Y-TRANSPORT AND DIFFUSION OF FF ON U-GRID.
C Y-TRANPORT SUPPOSES BOUNDARY VELOCITIES & FI TO BE ZERO.
      INCLUDE '0COM.INC'
      INCLUDE '0TRUVS.INC'
      REAL TAU
	REAL FF(NX,NY,NZ),VV(NX,NY,NZ),AY(NX,NY,NZ)
      REAL A(NY),B(NY),C(NY),ETA(NY),RKSI(NY)
	INTEGER M, N, K, IG, IGRY
      INTEGER II, JJ
      REAL    BP, DP, DM, PM, PP

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,IG,II,JJ,
!$OMP&        BP,DP,DM,PP,PM,A,B,C,ETA,RKSI)
      DO  K=1,NZ
       DO  M=MMM,MM
        DO IG=1,LRYU(M)
         II=IIYU(IG,M)
         JJ=JJYU(IG,M)+1
         DO  N=II,JJ

          BP =  HHU(M,N)*DXT(M,N)*DYH(M,N)*RN/ TAU

          DP = ( AY(M  ,N  ,K) + AY(M+1,N  ,K)
     &         + AY(M  ,N+1,K) + AY(M+1,N+1,K) )*HH(M,N  )
     &           *DXB(M,N  )/DYB(M,N  )/4./RN
             
          DM = ( AY(M  ,N  ,K) + AY(M+1,N  ,K)
     &         + AY(M  ,N-1,K) + AY(M+1,N-1,K) )*HH(M,N-1)
     &           *DXB(M,N-1)/DYB(M,N-1)/4./RN 
          

          PP =( VV(M  ,N  ,K)*DXH(M  ,N  )*HHV(M  ,N  )
     &         +VV(M+1,N  ,K)*DXH(M+1,N  )*HHV(M+1,N  )  )/4.0

          PM =( VV(M  ,N-1,K)*DXH(M  ,N-1)*HHV(M  ,N-1)
     &         +VV(M+1,N-1,K)*DXH(M+1,N-1)*HHV(M+1,N-1)  )/4.0


C LAST "2" ABOVE IS DUE TO APROXIMATION IN SPACE.
          C(N) =  PP * TRUV_IMP - DP
          A(N) = -PM * TRUV_IMP - DM
          B(N) =  BP + DP + DM
          ETA(N) = BP * FF(M,N,K) - PP * TRUV_EXP * FF(M,N+1,K)
     &                            + PM * TRUV_EXP * FF(M,N-1,K)

	   END DO
C IMPLICIT: A(N)=-DM-PM;C(N)=-DP+PP;B(N)=BP+DM+DP;ETA(N)=BP*FF(M,N,K)
C IF (N=II) : PM = 0, A = -DM ; IF (N=JJ) : PP = 0, C = -DP.
C IF(IGRY.EQ.1) : POINTS OUTSIDE U-GRID ARE TAKEN, SO THEY MUST BE = 0,
C OR U', FI = 0. SOURCE IN ETA APPEAR.
         IF (IGRY.EQ.2) THEN
          B(II) = B(II) + A(II)
          B(JJ) = B(JJ) + C(JJ)
         ELSE
          B(II) = B(II) - A(II)
          B(JJ) = B(JJ) - C(JJ)
 	   ENDIF

	   CALL FACTOR(NY,A,B,C,ETA,RKSI,II,JJ)

	   DO  N=II,JJ
          FF(M,N,K) = RKSI(N)
	   END DO

	  END DO
	 END DO
	END DO
!$OMP END PARALLEL DO

      RETURN
      END

C======================================================================
      SUBROUTINE TRXV(FF,UU,AX,TAU,IGRX)
      IMPLICIT NONE
C----------------------------------------------------------------------
C X-TRANSPORT AND DIFFUSION OF FF ON U-GRID.
C X-TRANPORT SUPPOSES BOUNDARY VELOCITIES & FI TO BE ZERO.
      INCLUDE '0COM.INC'
      INCLUDE '0TRUVS.INC'
      REAL TAU
	REAL FF(NX,NY,NZ),UU(NX,NY,NZ),AX(NX,NY,NZ)
      REAL A(NX),B(NX),C(NX),ETA(NX),RKSI(NX)
      INTEGER M, N, K, IG, IGRX
      INTEGER II, JJ, MLOOP
      REAL    BP, DP, DM, PM, PP
      INTEGER M1, M9

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,IG,II,JJ,MLOOP,M1,M9,
!$OMP&        BP,DP,DM,PP,PM,A,B,C,ETA,RKSI)
      DO K=1,NZ
       DO N=NNN,NN
        DO IG=1,ABS(LRXU(N))
         II=IIXU(IG,N)
         JJ=JJXU(IG,N)+1
         IF (JJ.GT.MM) JJ = JJ - MMD

         DO MLOOP=II,JJXU(IG,N)
          M  = MLOOP
          M1 = MLOOP - 1
          M9 = MLOOP + 1
          IF ( M.GT.MM) M  = M  - MMD
          IF (M1.GT.MM) M1 = M1 - MMD
          IF (M9.GT.MM) M9 = M9 - MMD

          BP = HHV(M,N)*DXH(M,N)*DYT(M,N)*RN / TAU


          DP = (AX(M ,N  ,K) + AX(M ,N+1,K)
     &        + AX(M9,N  ,K) + AX(M9,N+1,K) )*HH(M ,N)
     &          *DYB(M ,N) /DXB(M ,N) /4./RN
             
          DM = (AX(M ,N  ,K) + AX(M ,N+1,K)
     &        + AX(M1,N  ,K) + AX(M1,N+1,K) )*HH(M1,N)
     &          *DYB(M1,N) /DXB(M1,N) /4./RN 
          
          PP =( UU(M  ,N  ,K)*DYH(M  ,N  )*HHU(M,N  )
     &         +UU(M  ,N+1,K)*DYH(M  ,N+1)*HHU(M,N+1) )/4.0

          PM =( UU(M1,N  ,K)*DYH(M1,N  )*HHU(M1,N  )
     &         +UU(M1,N+1,K)*DYH(M1,N+1)*HHU(M1,N+1)  )/4.0


C LAST "2" ABOVE IS DUE TO APROXIMATION IN SPACE.
          C(M) =  PP * TRUV_IMP - DP
          A(M) = -PM * TRUV_IMP - DM
          B(M) =  BP + DP + DM
          ETA(M) = BP * FF(M,N,K) - PP * FF(M9,N,K) * TRUV_EXP
     &                            + PM * FF(M1,N,K) * TRUV_EXP
         ENDDO
C IMPLICIT: A(M)=-DM-PM;C(M)=-DP+PP;B(M)=BP+DM+DP;ETA(M)=BP*FF(M,N,K)
C IF (M=II) : PM = 0, A = -DM ; IF (M=JJ) : PP = 0, C = -DP.
C IF(IGRX.EQ.1) : POINTS OUTSIDE T-GRID IS TAKEN, SO THEY MUST BE = 0,
C OR U', FI = 0.
         IF (LRXU(N).GT.0) THEN
          IF (IGRX.EQ.2) THEN
           B(II) = B(II) + A(II)
           B(JJ) = B(JJ) + C(JJ)
          ELSE
           B(II) = B(II) - A(II)
           B(JJ) = B(JJ) - C(JJ)          
          ENDIF
          CALL FACTOR2(NX,A,B,C,ETA,RKSI,II,JJXU(IG,N))
         ELSE
          CALL FACTORC(A,B,C,ETA,RKSI,MMM,MM)
         ENDIF
         DO MLOOP=II,JJXU(IG,N)
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
      SUBROUTINE TRYV(FF,VV,AY,TAU)
      IMPLICIT NONE
C------------------------------------------------------09-13-96 06:08pm
C Y-TRANSPORT AND DIFFUSION OF FF ON V-GRID.
C Y-TRANPORT SUPPOSES BOUNDARY VELOCITIES & FI TO BE ZERO.
      INCLUDE '0COM.INC'
	INCLUDE '0TRUVS.INC'
      REAL TAU
      REAL FF(NX,NY,NZ),VV(NX,NY,NZ),AY(NX,NY,NZ)
      REAL A(NY),B(NY),C(NY),ETA(NY),RKSI(NY)
	INTEGER M, N, K, IG, IGRY
      INTEGER II, JJ
      REAL    BP, DP, DM, PM, PP

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,IG,II,JJ,
!$OMP&        BP,DP,DM,PP,PM,A,B,C,ETA,RKSI)
      DO K=1,NZ
       DO M=MMM,MM
        DO IG=1,LRY(M)
         II=IIY(IG,M)
         JJ=JJY(IG,M)-1
         DO N=II,JJ
          BP = HHV(M,N)*DXH(M,N)*DYT(M,N)*RN/TAU


          DP =  AY(M,N+1,K)*HHQ(M,N+1)*DX(M,N+1)/DY(M,N+1)/RN
          DM =  AY(M,N  ,K)*HHQ(M,N  )*DX(M,N  )/DY(M,N  )/RN          	

          PP =( VV(M,N+1,K)*HHV(M,N+1)*DXH(M,N+1)
     &         +VV(M,N  ,K)*HHV(M,N  )*DXH(M,N  ) )/4.0
          
          PM =( VV(M,N-1,K)*HHV(M,N-1)*DXH(M,N-1)
     &         +VV(M,N  ,K)*HHV(M,N  )*DXH(M,N  ) )/4.0


C LAST "2" ABOVE IS DUE TO APROXIMATION IN SPACE.
          C(N) =  PP * TRUV_IMP - DP
          A(N) = -PM * TRUV_IMP - DM
          B(N) =  BP + DP + DM
          ETA(N) = BP * FF(M,N,K) - PP * FF(M,N+1,K) * TRUV_EXP
     &                            + PM * FF(M,N-1,K) * TRUV_EXP

	   END DO
C IMPLICIT: A(N)=-DM-PM;C(N)=-DP+PP;B(N)=BP+DM+DP;ETA(N)=BP*FF(M,N,K)
C IF (N=II) : PM = 0, A = -DM ; IF (N=JJ) : PP = 0, C = -DP.
C IF(IGRY.EQ.1) : POINTS OUTSIDE U-GRID ARE TAKEN, SO THEY MUST BE = 0,
C OR U', FI = 0. SOURCE IN ETA APPEAR.


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
c     Bottom friction parametrization is added
C======================================================================
      SUBROUTINE TRZU(UU,FF0,WW,A_Z,TAU,IGRZ)
      IMPLICIT NONE
C----------------------------------------------------------------------
C VERTICAL TRANSPORT AND 3D DIFFUSION OF FF ON U-GRID.
      INCLUDE '0COM.INC'
      INCLUDE '0TRUVS.INC'
      INCLUDE '0BOTTFR.INC'
      REAL TAU
      REAL UU(NX,NY,NZ),FF0(NX,NY),WW(NX,NY,NZ+1),A_Z(NX,NY,NZ+1)
      REAL A(NZ),B(NZ),C(NZ),ETA(NZ),RKSI(NZ)
      INTEGER M, N, K
      INTEGER IGRZ
      REAL    BP, DP, DM, PM, PP
      REAL    FORCE

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,FORCE,
!$OMP&        BP,DP,DM,PP,PM,A,B,C,ETA,RKSI)
      DO N=NNN,NN
       DO M=MMM,MM
        IF (LCU(M,N).GT.0.5) THEN
         BP = HHU(M,N) / TAU
	   FORCE =(FF0(M+1,N)*DY(M+1,N)+FF0(M,N)*DY(M,N))   !WIND FROM T-GRID
     &                    /(2*DYH(M,N))   
C INTERNAL POINTS.
         DO K=2,NZ-1
          DP =(A_Z(M,N,K+1)+A_Z(M+1,N,K+1))/2.0/HHU(M,N)/HZT(K+1)/DZ(K)
          DM =(A_Z(M,N,K  )+A_Z(M+1,N,K  ))/2.0/HHU(M,N)/HZT(K  )/DZ(K)
          
          PP =( WW(M  ,N,K+1)*DY(M  ,N)*HHQ(M  ,N)
     &        + WW(M+1,N,K+1)*DY(M+1,N)*HHQ(M+1,N))
     &         /(DYH(M,N)*HHU(M,N)) / DZ(K) / 4.0
          PM =( WW(M  ,N,K  )*DY(M  ,N)*HHQ(M  ,N)
     &        + WW(M+1,N,K  )*DY(M+1,N)*HHQ(M+1,N)) 
     &         /(DYH(M,N)*HHU(M,N)) / DZ(K) / 4.0

          C(K) =  PP * TRUV_IMP - DP
          A(K) = -PM * TRUV_IMP - DM
          B(K) =  BP + DP + DM
          ETA(K) = BP*UU(M,N,K) - PP * UU(M,N,K+1) * TRUV_EXP
     &                          + PM * UU(M,N,K-1) * TRUV_EXP
         ENDDO
C IMPLICIT: A(K)=-DM-PM;C(K)=-DP+PP;B(K)=BP+DP+DM;ETA(K)=BP*UU(M,N,K)
C SURFACE POINT: TWO CASES.
         K = 1
         DP =(A_Z(M,N,K+1)+A_Z(M+1,N,K+1))/2.0/HHU(M,N)/HZT(K+1)/DZ(K)
         DM =(A_Z(M,N,K  )+A_Z(M+1,N,K  ))/2.0/HHU(M,N)/HZT(K  )/DZ(K)
         
         PP =( WW(M  ,N,K+1)*DY(M  ,N)*HHQ(M  ,N)
     &       + WW(M+1,N,K+1)*DY(M+1,N)*HHQ(M+1,N)) 
     &        /(DYH(M,N)*HHU(M,N))/ DZ(K) / 4.0

         C(K) =  PP * TRUV_IMP - DP
C AZD -> AZ1 IN DM, PM IS 0 DUE TO WW(M,N,1) = 0.
         IF(IGRZ.EQ.1) THEN
          B(K) =  BP + DP + DM
          ETA(K) = BP*UU(M,N,K) - PP * UU(M,N,K+1) * TRUV_EXP
     &                          + DM * FORCE
         ELSEIF(IGRZ.EQ.2) THEN
          B(K) =  BP + DP
          ETA(K) = BP*UU(M,N,K) - PP * UU(M,N,K+1) * TRUV_EXP
     +                          + FORCE / DZ(K)
         ENDIF
C BOTTOM POINT.
         K = NZ
         DM =(A_Z(M,N,K  )+A_Z(M+1,N,K  ))/2.0/HHU(M,N)/HZT(K  )/DZ(K)
         PM =( WW(M  ,N,K  )*DY(M  ,N)*HHQ(M  ,N)
     &       + WW(M+1,N,K  )*DY(M+1,N)*HHQ(M+1,N))
     &        /(DYH(M,N)*HHU(M,N)) / DZ(K) / 4.0
         
         A(K) = -PM * TRUV_IMP - DM
         B(K) =  BP + DM
c----------------------------------------------------------------------
c        Bottom Friction: set it on diagonal (implicit friction)
     +        + (BottomFriction(M,N)+ BottomFriction(M+1,N))/2.0
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c        Set rhs bottom point
c----------------------------------------------------------------------
         ETA(K) =  BP*UU(M,N,K) + PM * UU(M,N,K-1) * TRUV_EXP
c----------------------------------------------------------------------
c        Bottom Friction: set it on diagonal (explicit friction)
c    -     - (BottomFriction(M,N)+ BottomFriction(M+1,N))/2.0*UU(M,N,K)
c----------------------------------------------------------------------
         CALL FACTOR(NZ,A,B,C,ETA,RKSI,1,NZ)
         DO K=1,NZ
          UU(M,N,K)=RKSI(K)
         ENDDO
        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
C======================================================================
c     Bottom friction parametrization is added
C======================================================================
      SUBROUTINE TRZV(VV,FF0,WW,A_Z,TAU,IGRZ)
      IMPLICIT NONE
C----------------------------------------------------------------------
C VERTICAL TRANSPORT AND 3D DIFFUSION OF VV ON V-GRID.
      INCLUDE '0COM.INC'
      INCLUDE '0TRUVS.INC'
      INCLUDE '0BOTTFR.INC'
      REAL TAU
      REAL VV(NX,NY,NZ),FF0(NX,NY),WW(NX,NY,NZ+1),A_Z(NX,NY,NZ+1)
      REAL A(NZ),B(NZ),C(NZ),ETA(NZ),RKSI(NZ)
      INTEGER M, N, K
      INTEGER IGRZ
      REAL    BP, DP, DM, PM, PP
      REAL    FORCE

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,FORCE,
!$OMP&        BP,DP,DM,PP,PM,A,B,C,ETA,RKSI)
      DO N=NNN,NN
       DO M=MMM,MM
        IF (LCV(M,N).GT.0.5) THEN
         BP = HHV(M,N) / TAU
C  WIND STRESS INTERPOLATION FROM T-GRID
C        SIMPLE:
C	   FORCE =(FF0(M,N+1) + FF0(M,N))/2.0   
C        WITH COSINE OF LATITUDE:
	   FORCE =(FF0(M,N+1)*DX(M,N+1)+FF0(M,N)*DX(M,N))
     &          /(2.0*DXH(M,N)) 

C INTERNAL POINTS.
         DO K=2,NZ-1
          DP =(A_Z(M,N,K+1)+A_Z(M,N+1,K+1))/2.0/HHV(M,N)/HZT(K+1)/DZ(K)
          DM =(A_Z(M,N,K  )+A_Z(M,N+1,K  ))/2.0/HHV(M,N)/HZT(K  )/DZ(K)
          
          PP =( WW(M,N  ,K+1)*DX(M,N  )*HHQ(M,N  )
     &         +WW(M,N+1,K+1)*DX(M,N+1)*HHQ(M,N+1))
     &         /(DXH(M,N)*HHV(M,N)) / DZ(K)/4.0
          PM =( WW(M,N  ,K  )*DX(M,N  )*HHQ(M,N  )
     &         +WW(M,N+1,K  )*DX(M,N+1)*HHQ(M,N+1))
     &         /(DXH(M,N)*HHV(M,N)) / DZ(K)/4.0

          C(K) =  PP * TRUV_IMP - DP
          A(K) = -PM * TRUV_IMP - DM
          B(K) =  BP + DP + DM
          ETA(K) = BP*VV(M,N,K) - PP * VV(M,N,K+1) * TRUV_EXP
     &                          + PM * VV(M,N,K-1) * TRUV_EXP
         ENDDO
C IMPLICIT: A(K)=-DM-PM;C(K)=-DP+PP;B(K)=BP+DP+DM;ETA(K)=BP*VV(M,N,K)
C SURFACE POINT: TWO CASES.
         K = 1
         DP =(A_Z(M,N,K+1)+A_Z(M,N+1,K+1))/2.0/HHV(M,N)/HZT(K+1)/DZ(K)
         DM =(A_Z(M,N,K  )+A_Z(M,N+1,K  ))/2.0/HHV(M,N)/HZT(K  )/DZ(K)
         
         PP =( WW(M,N  ,K+1)*DX(M,N  )*HHQ(M,N  )
     &        +WW(M,N+1,K+1)*DX(M,N+1)*HHQ(M,N+1))
     &         /(DXH(M,N)*HHV(M,N)) / DZ(K)/4.0
         C(K) =  PP * TRUV_IMP - DP
C AZD -> AZ1 IN DM, PM IS 0 DUE TO WW(M,N,1) = 0.
         IF(IGRZ.EQ.1) THEN
          B(K)   = BP + DP + DM
          ETA(K) = BP*VV(M,N,K) - PP * VV(M,N,K+1) * TRUV_EXP
     &                          + DM * FORCE
         ELSEIF(IGRZ.EQ.2) THEN
          B(K)   = BP + DP
          ETA(K) = BP*VV(M,N,K) - PP * VV(M,N,K+1) * TRUV_EXP
     &                          + FORCE / DZ(K)
         ENDIF
c
C BOTTOM POINT.
         K = NZ
         DM =(A_Z(M,N,K  )+A_Z(M,N+1,K  ))/2.0/HHV(M,N)/HZT(K  )/DZ(K)
         PM =( WW(M,N  ,K  )*DX(M,N  )*HHQ(M,N  )
     &        +WW(M,N+1,K  )*DX(M,N+1)*HHQ(M,N+1))
     &         /(DXH(M,N)*HHV(M,N))/DZ(K)/4.0
         PP = 0.0

         A(K) = -PM * TRUV_IMP - DM
         B(K) =  BP + DM
c----------------------------------------------------------------------
c        Bottom Friction: set it on diagonal (implicit friction)
     +        + (BottomFriction(M,N)+ BottomFriction(M,N+1))/2.0
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c        Set the rhs on bottom point
c----------------------------------------------------------------------
          ETA(K) = BP*VV(M,N,K) + PM * VV(M,N,K-1) * TRUV_EXP
c----------------------------------------------------------------------
c        Bottom Friction: set it on diagonal (explicit friction)
c    -     - (BottomFriction(M,N)+ BottomFriction(M,N+1))/2.0*VV(M,N,K)
c----------------------------------------------------------------------

          CALL FACTOR(NZ,A,B,C,ETA,RKSI,1,NZ)
         DO K=1,NZ
          VV(M,N,K)=RKSI(K)
         ENDDO
        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
