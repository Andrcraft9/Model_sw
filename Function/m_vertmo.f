C======================================================================
C OCEAN GENERAL CIRCULATION MODEL PARAMETRIZATIONS OF VERTICAL MIXING
C (C) DIANSKY N.A., 2008.   
C======================================================================
      SUBROUTINE MONIN_OBUHOV(DEN,UU,VV,RIT,
     &                        ANZT,ANUMAXT,ANUBGRT,
     &                        ANZU,ANUMAXU,ANUBGRU,TAUX,TAUY,AFNT)
      IMPLICIT NONE
C----------------------------------------------------------------------
*     Vertical Mixing Coefficients Calculation
*     according to Monin-Obukhov (or Kochergin's) parametrizations .
C ANU = ((D U/D Z)^2+(D V/D Z)^2-CONST*G/R0*(D DEN/D Z) )
C HZT(K) = Z(K)-Z(K-1), K=2,NZ

      INCLUDE '0COM.INC'
      INCLUDE '1SEAICE.INC'

      REAL  DEN(NX,NY,NZ  ),    !POTENTIAL DENSYTY
     &       UU(NX,NY,NZ  ),	!ZONAL VELOCITY
     &       VV(NX,NY,NZ  ),    !MERIDIONAL VELOCITY
     &      RIT(NX,NY,NZ  ),    !ARRAY FOR RICHRDSON IS USED FOR STORE 
                              !              VAISJALIJA BRENDTA FREQ
     &     TAUX(NX,NY     ),	!ZONAL WIND STRESS
     &     TAUY(NX,NY     ),  !MERIDIONAL WIND STRESS
     &     ANZT(NX,NY,NZ+1),  !DIFFUSION COEFFICIENT
     &     ANZU(NX,NY,NZ+1),  !VISCOSITY COEFFICIENT
     &     AFNT(NX,NY,NZ+1),
     &     ANUMAXT,ANUBGRT,   !MAX AND MIN DIFFUSION COEFFICIENT
     &     ANUMAXU,ANUBGRU    !MAX AND MIN VISCOSITY COEFFICIENT

      REAL R1,R2,U1,U2,V1,V2,CKro,Q05,VBF,C20,C21,ZZZ,A2,Uz,Vz,DZ21,
     &     ANUPU,ANUPT !ADDITIONAL MIXING IN UPPER LAYER 
      INTEGER M, N, K

	INTEGER, ALLOCATABLE:: KTOP(:,:),KBOT(:,:)   !HELP ARRAYS FOR UPPER
                                                   !AND BOTTOM LEVELS WHERE 
	                                             !MIXING IS NEED
* IAKOVLEV PARAMETRIZATION
*     PARAMETERS FOR COEFFICIENTS CUTOFFS.
      REAL D,SL,DIMDEPTH,UNDIMDEPTH,NU_UNSTABLE
      PARAMETER (DIMDEPTH=300.0,ANUPU=500.0,ANUPT=500.0) !UPPER MIXED LAYER COEFFICIENTS ARE BIGGER
      PARAMETER (NU_UNSTABLE=500.0)

*     DEPTH SCALE PARAMETER [CM]
*     PARAMETER (D=2500.0,SL= 0.05*D)  ! MARHUK ET. AL. 1976 0.05
*     PARAMETER (D=2500.0, C20=0.4/D, C21=1.2/D/D)  !Iakovlev, 2008
*     PARAMETER (D=1250.0, C20=0.4/D, C21=2.0/D/D)  !Diansky,  2008
*     PARAMETER (D=1000.0,C20=0.4/D,C21=3.0/D/D)    !Diansky,  2008
	
      REAL PIP180,DPIP180              !FOR DEGREES TO RADIANS CONVERS
	PARAMETER(PIP180=3.141592653/180.0,
     &         DPIP180=3.1415926535897/180.0D00)

C--------ADDITIONAL UPPER MIXING PARAMETERS FROM VOLDIN----------------------
	REAL ANUPU2, ANUPT2
      PARAMETER (ANUPU2=500.0,ANUPT2=500.0) !UPPER MIXED LAYER COEFFICIENTS ARE BIGGER

C VOLODIN'S PARAMETRIZATION IN UPPER LAYER	
	REAL STRMAX1,STRMIN1,STRMAX2,STRMIN2,
     &     HMMAX1,HMMAX2,HMMIN1,HMMIN2,DCRT

      PARAMETER(
     &         STRMIN1=0.4,	
     &		  STRMAX1=1.0,	
     &		  STRMIN2=1.0,	
     &		  STRMAX2=8.0,	
     &		  HMMIN1=300.0,	
     &		  HMMAX1=3000.0,	
     &		  HMMIN2=3000.0,
     &		  HMMAX2=10000.0,	
     &		  DCRT=300.0)						
c      PARAMETER(
c     &          STRMIN1=0.4,	
c     &		  STRMAX1=1.0,	
c     &		  STRMIN2=1.0,	
c     &		  STRMAX2=8.0,	
c     &		  HMMIN1=200.0,	
c     &		  HMMAX1=2000.0,	
c     &		  HMMIN2=2000.0,
c     &		  HMMAX2=8000.0,	
c     &		  DCRT=500.0)
	
	REAL STR,HM,ARG,SANZT,ANZT1(NZ)
	INTEGER KK
C----------------------------------------------------------------------


C FOR MOMENTUM:
!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,R2,U2,V2,DZ21,R1,U1,V1,D,SL,VBF,CKro,UNDIMDEPTH,
!$OMP&        Uz,Vz,Q05,STR,HM,ANZT1,SANZT,KK,ARG)
      DO N=NNN,NN
       DO M=MMM,MM
        IF (LU(M,N).GT.0.5) THEN

          UNDIMDEPTH = DIMDEPTH/HHQ(M,N)

          R2 = DEN(M,N,1)

          U2 = LU(M,N)*0.5/DY(M,N)*                   !U ON T-GRID
     &        (UU(M-1,N,1)*HHU(M-1,N)*DYH(M-1,N) +
     &         UU(M  ,N,1)*HHU(M  ,N)*DYH(M  ,N) )/HHQ(M,N)

          V2 = LU(M,N)*0.5/DX(M,N)*                   !V ON T-GRID
     &        (VV(M,N-1,1)*HHV(M,N-1)*DXH(M,N-1) +
     &         VV(M,N  ,1)*HHV(M,N  )*DXH(M,N  ) )/HHQ(M,N)

	   DO K=2,NZ

          DZ21= HZT(K)*HHQ(M,N)
           R1 = R2
           R2 = DEN(M,N,K)
           U1 = U2
           V1 = V2

          U2 = LU(M,N)*0.5/DY(M,N)*                   !U ON T-GRID
     &        (UU(M-1,N,K)*HHU(M-1,N)*DYH(M-1,N) +
     &         UU(M  ,N,K)*HHU(M  ,N)*DYH(M  ,N) )/HHQ(M,N)

          V2 = LU(M,N)*0.5/DX(M,N)*                   !V ON T-GRID
     &        (VV(M,N-1,K)*HHV(M,N-1)*DXH(M,N-1) +
     &         VV(M,N  ,K)*HHV(M,N  )*DXH(M,N  ) )/HHQ(M,N)

***************Iakovlev procerure**************************
*C         zzz=0.5*(z(k)+z(k+1))       ! Rough wall? Waves?
*
*          IF( zzz .LE. 0.5*D) then
*            A2= (D-zzz)*zzz
*            SL=C20*A2*(1.-C21*A2)
*          ELSE
*            SL= 0.07*D                   !Polyakov&Kowalik 0.07
*cc          SL= 0.05*hz(k)               !Old Marchuk paper
*cc          SL= 0.05*MAX(D,hz(k))        !Combination
*          END IF
*
*		SL=SL**2
*
**     Simplified Vasala-Brendt **2 Frequency  
*      ro1=ro(i,j,k  )                                 
*      ro2=ro(i,j,k+1)
*      VB= 2.*g*(ro2-ro1)/hz(k)/(2.+ro1+ro2)
*     Vertical shift
*      uz= (u(i,j,k+1)-u(i,j,k))/hz(k)
*      vz= (v(i,j,k+1)-v(i,j,k))/hz(k)
*     Formulas for momentum
*      Q= MAX( 0., uz*uz +vz*vz -CKro*VB)
*      az(i,j,k)= MIN(azmax, azmin+SL*SQRT(Q))
*
**     Heat ans Salt
*
*c      IF(VB.GE.0.) then
*c      azt(i,j,k)= aztbk      ! mixing is due to advection only -AWI.
*c      azs(i,j,k)= aztbk
*c      else
*c      azt(i,j,k)= amixer
*c      azs(i,j,k)= amixer
*c      end if
*
*
**     Standart Kochergin&Polyakov
*      azt(i,j,k)= MIN(amixer,aztbk+0.01*SL*SQRT(Q))
*      azs(i,j,k)= MIN(amixer,azsbk+0.01*SL*SQRT(Q))
*
*
*      IF(k.LE.kmixed) then             ! Upper mixed layer
*      azt(i,j,k)= azt(i,j,k)+amixer
*      az (i,j,k)= az (i,j,k)+amixer
*      azs(i,j,k)= azs(i,j,k)+amixer
*      end if
*
*
************************************************************

***********Analogous procedure in sigma coordinate**********************
*     PARAMETER (D=1000.0, C20=0.4/D, C21=3.0/D/D)  !Diansky,  2008
	
 	D=MAX(200.0,1000.0*HHQ(M,N)/350000.0)
c	C20=0.4/D, C21=3.0/D/D

c		ZZZ=ZW(K)*HHQ(M,N)
c          IF( ZZZ .LE. 0.5*D) then
c            A2= (D-ZZZ)*ZZZ
c            SL=C20*A2*(1.-C21*A2)
c          ELSE
c            SL= 0.025*D                  !Diansky, 2008
ccc          SL= 0.07*D                   !Polyakov&Kowalik 0.07
ccc          SL= 0.05*hz(k)               !Old Marchuk paper
ccc          SL= 0.05*MAX(D,hz(k))        !Combination
c          END IF

          SL= 0.025*D                  !Diansky, 2008
		SL=SL**2

*     Simplified Vasala-Brendt **2 Frequency  
          VBF= G*(R2-R1)/DZ21  !VAISJALIJA BRENDTA FREQUENCY

	    RIT(M,N,K)=VBF       !ARRAY FOR RICHRDSON IS USED FOR STORING OF
		                     !VAISJALIJA BRENDTA FREQUENCY


*     Vertical shift
          Uz = (U2-U1)/DZ21
          Vz = (V2-V1)/DZ21


		IF(VBF.GE.0.) THEN	 ! Stable Stratification.
			CKro=1E-02       ! 0.01 - Kochergin, 1987

		Q05= SQRT( AMAX1( 0.0, Uz**2 +Vz**2 - CKro*VBF) )

C FOR MOMENTUM:
c		       ANZU(M,N,K)= AMIN1(ANUMAXU*(1.0-AFNT(M,N,K))
c     &                           +ANUBGRU*AFNT(M,N,K),ANUBGRU + SL*Q05)
		       ANZU(M,N,K)= AMIN1(ANUMAXU,ANUBGRU + SL*Q05)

C FOR TEMPERATURE AND SALINITY:
		       ANZT(M,N,K)= AMIN1(ANUMAXT,ANUBGRT + 0.01*SL*Q05)

            ELSE                 ! Unstable Stratification.

                  CKro=1E+03       ! Kochergin, 1987
		Q05= SQRT( AMAX1( 0.0, Uz**2 +Vz**2 - CKro*VBF) )
                  
C FOR MOMENTUM:
c		      ANZU(M,N,K)= AMIN1(ANUMAXU*(1.0-AFNT(M,N,K))
c     &                           +ANUBGRU*AFNT(M,N,K),ANUBGRU + SL*Q05)
		      ANZU(M,N,K)= AMIN1(ANUMAXU,ANUBGRU + SL*Q05)
C FOR TEMPERATURE AND SALINITY:
c		      ANZT(M,N,K)= NU_UNSTABLE*(1.0-AFNT(M,N,K))
c     &                                 +ANUMAXT*AFNT(M,N,K)
		      ANZT(M,N,K)= NU_UNSTABLE
		END IF

C FOR UPPER LAYER
          IF(ZW(K).LE.UNDIMDEPTH) THEN

C FOR MOMENTUM:
		ANZU(M,N,K)= AMAX1(ANUPU,ANZU(M,N,K))

C FOR TEMPERATURE AND SALINITY:
		ANZT(M,N,K)= AMAX1(ANUPT*AICE0(M,N),ANZT(M,N,K))
       
	    END IF

         ENDDO

C---- ADDITIONAL DIFFUSION DUE TO PROPAGATION OF WAVE FLUCTUATIONS--
         STR=RH0*SQRT(TAUX(M,N)*TAUX(M,N)+TAUY(M,N)*TAUY(M,N))*
     *    AICE0(M,N)
c         strout(m,n)=str
         STR=AMAX1(STR,STRMIN1)
         STR=AMIN1(STR,STRMAX2)
         IF(STR.LE.STRMAX1) THEN
            HM=HMMIN1+(HMMAX1-HMMIN1)*(STR-STRMIN1)/(STRMAX1-STRMIN1)
         ELSE
            HM=HMMIN2+(HMMAX2-HMMIN2)*(STR-STRMIN2)/(STRMAX2-STRMIN2)
         END IF
         DO K=2,NZ
            IF(ZW(K)*HHQ(M,N).LE.HM) THEN
               ANZU(M,N,K)=AMAX1(ANZU(M,N,K),ANUPU2)
               ANZT(M,N,K)=AMAX1(ANZT(M,N,K),ANUPT2)
            END IF
         END DO
         DO K=2,NZ
            ANZT1(K)=ANZT(M,N,K)
         END DO

C---- ADDITIONAL DIFFUSION DUE TO PROPAGATION OF INSTABILITY FLUCTUATIONS 
         DO K=2,NZ

            IF(ANZT1(K).LT.100.0.AND.HHQ(M,N)*ZW(K).LE.30000.) THEN
               SANZT=0.0
               DO KK=2,NZ
                  IF(K.NE.KK) THEN
                     ARG=-ABS(ZW(K)-ZW(KK))*HHQ(M,N)/DCRT
                     SANZT=AMAX1(SANZT,(ANZT(M,N,KK)-ANUBGRT)*EXP(ARG))
                  END IF
               END DO
               ANZT(M,N,K)=ANZT(M,N,K)+SANZT
            END IF
         END DO


        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

C---ADDITION CORRECTION OF UPPER MIXING FOR TEMPERATURE AND SALINITY-----------
C       FOR UPPER CONVECTIVE INSTABILITY

c	ALLOCATE (KTOP(NX,NY),KBOT(NX,NY))
c
c      CALL DEF_TS_UPPER_MIX_LEV(DEN,KTOP,KBOT,LU,NX,NY,NZ) 
c
c!$OMP PARALLEL DO PRIVATE(M,N,K)
c
c      DO N=1,NY
c         DO M=1,NX
c            
c            IF(LU(M,N).GT.0.5.AND.YT(N).LE.LAT_CRIT_4D) THEN
c
c                K=KTOP(M,N)+1
c	          
c                DO WHILE(K.LE.KBOT(M,N))
c		    
c                ANZT(M,N,K)=ANUMAXT  
c                K=K+1
c
c	          END DO
c
c            END IF
c         END DO
c      END DO
c!$OMP END PARALLEL DO
c
c	DEALLOCATE (KBOT,KTOP)
C------------------------------------------------------------------------------


      IF(MMD.NE.0) CALL CYCLIZE(ANZU,NX,NY,NZ,MMM,MM) !ENOUGH FOR U,V ONLY



      RETURN
      END
C======================================================26.11.99 19:17
      SUBROUTINE RICHNUM(DEN,UU,VV,TAUX,TAUY,RIT)
      IMPLICIT NONE
C----------------------------------------------------------------------
C RICHARDSON' NUMBER ON U-GRID. (D UV/D Z .NE. 0). K=2,NZ.
C RICHNU = G/R0*(D DEN/D Z) / ((D U/D Z)^2+(D V/D Z)^2), D Z = HZT*HH
C HZT(K) = Z(K)-Z(K-1), K=2,NZ
C RIT - RICHARDSON NUMBER ON T-GRID

      INCLUDE '0COM.INC'
      REAL DEN(NX,NY,NZ),UU(NX,NY,NZ),VV(NX,NY,NZ),
     &     RIT(NX,NY,NZ),TAUX(NX,NY),TAUY(NX,NY)

      REAL R1,R2,U1,U2,V1,V2,U_STAR, RLT
      INTEGER M, N, K
      REAL    CORIOLIS

C FOR MOMENTUM:
!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,RLT,CORIOLIS,U_STAR,R2,U2,V2,R1,U1,V1)
      DO N=NNN,NN
       DO M=MMM,MM
        IF (LU(M,N).GT.0.5) THEN

          RLT=(RLH(M,N  )+RLH(M-1,N  )
     &        +RLH(M,N-1)+RLH(M-1,N-1))/4.0

          CORIOLIS=MAX(ABS(RLT),2.5E-05)   !FREEZING ON 10N,S

          U_STAR = SQRT(SQRT(TAUX(M,N)**2+TAUY(M,N)**2))
          RIT(M,N,1)=0.5*U_STAR/CORIOLIS

          R2 = DEN(M,N,1)

          U2 = LU(M,N)*0.5/DY(M,N)*                   !U ON T-GRID
     &        (UU(M-1,N,1)*HHU(M-1,N)*DYH(M-1,N) +
     &         UU(M  ,N,1)*HHU(M  ,N)*DYH(M  ,N) )/HHQ(M,N)

          V2 = LU(M,N)*0.5/DX(M,N)*                   !V ON T-GRID
     &        (VV(M,N-1,1)*HHV(M,N-1)*DXH(M,N-1) +
     &         VV(M,N  ,1)*HHV(M,N  )*DXH(M,N  ) )/HHQ(M,N)

	   DO K=2,NZ
          R1 = R2
          R2 = DEN(M,N,K)
          U1 = U2
          V1 = V2
          U2 = LU(M,N)*0.5/DY(M,N)*                   !U ON T-GRID
     &        (UU(M-1,N,K)*HHU(M-1,N)*DYH(M-1,N) +
     &         UU(M  ,N,K)*HHU(M  ,N)*DYH(M  ,N) )/HHQ(M,N)

          V2 = LU(M,N)*0.5/DX(M,N)*                   !V ON T-GRID
     &        (VV(M,N-1,K)*HHV(M,N-1)*DXH(M,N-1) +
     &         VV(M,N  ,K)*HHV(M,N  )*DXH(M,N  ) )/HHQ(M,N)

C DIRECT DISCRETE FORM OF Ri
          RIT(M,N,K)=MIN(G*(R2-R1)/2.0*(ZW(K+1)-ZW(K-1))*HHQ(M,N)/
     &                    ((U2-U1)**2+(V2-V1)**2+0.0025),1E+05)
         ENDDO

        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO
C     IF(MMD.NE.0) CALL CYCLIZE(RIT,NX,NY,NZ,MMM,MM)

      RETURN
      END
C======================================================02-09-99 08:41pm
      SUBROUTINE PPMIXT(RIT,AZH,AZ0,AZB,TAUX,TAUY)
      IMPLICIT NONE
      INCLUDE '0COM.INC'
      INCLUDE '1SEAICE.INC'
C----------------------------------------------------------------------
C VERT. MIXING, DEPENDING ON RICHARDSON' NUMBER.
C DIFF = AZ0 / (1+5*RI) + AZB,  K=2,NZ. (< NZ+1).
C NPOW := PP: T=1, U=2; MOM2: T=3, U=2.
      REAL      UNSTABLE, DIMDEPTH
      PARAMETER(UNSTABLE = 500.0, DIMDEPTH=300.0)
      REAL RIT(NX,NY,NZ),AZH(NX,NY,NZ+1)
      REAL AZ0,AZB,AZBADD,UNDIMDEPTH
      INTEGER M, N, K
C VOLODIN'S PARAMETRIZATION IN UPPER LAYER
	REAL STRMAX1,STRMIN1,STRMAX2,STRMIN2,
     &     HMMAX1,HMMAX2,HMMIN1,HMMIN2,DCRT
       REAL TTMIN,TTMAX,STR,HM
       REAL TAUX(NX,NY),TAUY(NX,NY)

      PARAMETER(
     &         STRMIN1=0.4,
     &		  STRMAX1=1.0,
     &		  STRMIN2=1.0,
     &		  STRMAX2=8.0,
     &		  HMMIN1=300.0,
     &		  HMMAX1=3000.0,
     &		  HMMIN2=3000.0,
     &		  HMMAX2=10000.0,
     &		  DCRT=300.0)

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,AZBADD,UNDIMDEPTH,STR,HM)
      DO N=1,NY
       DO M=1,NX

        IF (LU(M,N).GT.0.5) THEN
          STR=RH0*SQRT(TAUX(M,N)*TAUX(M,N)+TAUY(M,N)*TAUY(M,N))
          STR=AMAX1(STR,STRMIN1)
          STR=AMIN1(STR,STRMAX2)
          IF(STR.LE.STRMAX1) THEN
             HM=HMMIN1+(HMMAX1-HMMIN1)*(STR-STRMIN1)/(STRMAX1-STRMIN1)
          ELSE
             HM=HMMIN2+(HMMAX2-HMMIN2)*(STR-STRMIN2)/(STRMAX2-STRMIN2)
          END IF
          IF(AICE0(M,N).LT.0.65) THEN
             HM=DIMDEPTH*1.5
          END IF
          UNDIMDEPTH = AMAX1(DIMDEPTH,HM)/HHQ(M,N)
          DO K=2,NZ

           IF(ZW(K).LE.UNDIMDEPTH) THEN
             AZBADD=UNSTABLE              !ADDITION MIXING IN UPPER LAYER FOR T&S
           ELSE
             AZBADD=AZB
           END IF

           IF(RIT(M,N,K).GE.0.0) THEN
             AZH(M,N,K) = AZ0/(1.0+5.0*RIT(M,N,K))**3 + AZBADD
           ELSE
             AZH(M,N,K  ) = UNSTABLE
           END IF

          ENDDO

C  FOR UNSTABLE STRATIFICATION
C          DO K=2,NZ
C           IF(RIT(M,N,K).LT.0.0) THEN
C MIXING IS DECREASED IN VICINITY OF EQUATOR IN UNSTABLE SITUATION
C          AZH(M,N,K  ) = UNSTABLE*(1.0-0.7*(RN/RM(N))**10)
C           AZH(M,N,K  ) = UNSTABLE
C          AZH(M,N,K+1) = UNSTABLE
C CONVECTIVE DIFFUSION AS dZ*dZ*VAJSJALA-BRENTA FREQUENCY
c           AZH(M,N,K) = AZ0+
c    +     HHQ(M,N)*HZT(K)*
c    *     SQRT( ABS(DEN(M,N,K)-DEN(M,N,K-1))*G*HHQ(M,N)*HZT(K) )
c    +     0.3*HHQ(M,N)**2*HZT(K-1)*HZT(K)/(12.0*3600.0)

C           END IF
C          ENDDO

        ENDIF

       ENDDO
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
C======================================================02-09-99 08:41pm
      SUBROUTINE PPMIXU(RIT,AZH,AZ0,AZB,TAUX,TAUY)
      IMPLICIT NONE
      INCLUDE '0COM.INC'
      INCLUDE '1SEAICE.INC'
C----------------------------------------------------------------------
C MOMENTUM VERT. MIXING, DEPENDING ON RICHARDSON' NUMBER.
C DIFF = AZ0 / (1+5*RI)**2 + AZB, K=2,NZ. (< NZ+1).
C NPOW := PP: T=1, U=2; MOM2: T=3, U=2.
      REAL UNSTABLE,DIMDEPTH,UNDIMDEPTH,AZBADD
      PARAMETER(UNSTABLE = 5.0E+02)
      REAL RIT(NX,NY,NZ),AZH(NX,NY,NZ+1)
      REAL AZ0, AZB
      PARAMETER(DIMDEPTH=300.0)
      INTEGER M, N, K
C VOLODIN'S PARAMETRIZATION IN UPPER LAYER
	REAL STRMAX1,STRMIN1,STRMAX2,STRMIN2,
     &     HMMAX1,HMMAX2,HMMIN1,HMMIN2,DCRT
       REAL TAUX(NX,NY),TAUY(NX,NY)
       REAL HM,STR

      PARAMETER(
     &         STRMIN1=0.4,
     &		  STRMAX1=1.0,
     &		  STRMIN2=1.0,
     &		  STRMAX2=8.0,
     &		  HMMIN1=300.0,
     &		  HMMAX1=3000.0,
     &		  HMMIN2=3000.0,
     &		  HMMAX2=10000.0,
     &		  DCRT=300.0)

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,UNDIMDEPTH,AZBADD,STR,HM)
      DO N=1,NY
       DO M=1,NX

        IF (LU(M,N).GT.0.5) THEN
          STR=RH0*SQRT(TAUX(M,N)*TAUX(M,N)+TAUY(M,N)*TAUY(M,N))
          STR=AMAX1(STR,STRMIN1)
          STR=AMIN1(STR,STRMAX2)
          IF(STR.LE.STRMAX1) THEN
             HM=HMMIN1+(HMMAX1-HMMIN1)*(STR-STRMIN1)/(STRMAX1-STRMIN1)
          ELSE
             HM=HMMIN2+(HMMAX2-HMMIN2)*(STR-STRMIN2)/(STRMAX2-STRMIN2)
          END IF
          IF(AICE0(M,N).LT.0.65) THEN
             HM=DIMDEPTH*1.5
          END IF
          UNDIMDEPTH = AMAX1(DIMDEPTH,HM)/HHQ(M,N)
          DO K=2,NZ

           IF(ZW(K).LE.UNDIMDEPTH) THEN
             AZBADD=UNSTABLE              !ADDITION MIXING IN UPPER LAYER FOR T&S
           ELSE
             AZBADD=AZB
           END IF

           IF(RIT(M,N,K).GT.0.0) THEN

           AZH(M,N,K) = AZ0 /(1.0 + 5.0*RIT(M,N,K))**2 + AZB

           ELSE

           AZH(M,N,K  ) = AZ0  + AZB
C          AZH(M,N,K+1) = AZ0   ! + AZB

           END IF
          ENDDO
        ENDIF

       ENDDO
      ENDDO
!$OMP END PARALLEL DO
      IF(MMD.NE.0) CALL CYCLIZE(AZH,NX,NY,NZ,MMM,MM)

      RETURN
      END
C======================================================================
      SUBROUTINE TSMIX(TT,SS,DEN,DEEPMIX,LATMIN,LATMAX)
      IMPLICIT NONE
C----------------------------------------------------------------------
C COMPUTING DENSITY AND MIXING T,S CONSERVATIVELY.
      INCLUDE '0COM.INC'
      REAL TT(NX,NY,NZ),SS(NX,NY,NZ),DEN(NX,NY,NZ),
     &              DEEPMIX,       !DEEP OF MIXING[CM]
     &              LATMIN,LATMAX  !HORIZONTAL BOUNDARIES OF MIXING
      INTEGER NLATMIN,NLATMAX,N,M,K,K1,K2,KDEEPMIX
      REAL UNDIMDEEP
      INCLUDE '0DENP.INC'
      REAL    S1, S2, S3, TU, SU, DEN0


C DEFINITION LATITUDE BOUNDARY FOR TSMIX USING
      NLATMIN=MAX(NINT((LATMIN-RLAT)/DYST)+NNN,NNN)
      NLATMAX=MIN(NINT((LATMAX-RLAT)/DYST)+NNN,NN )
      IF(NLATMIN.GT.NLATMAX) THEN
       WRITE(*,*)'   ERORR IN SUBROUTINE TSMIX:'
       WRITE(*,*)'   NLATMIN > NLATMAX'
       STOP
      END IF

      DO N=NLATMIN,NLATMAX
      DO M=MMM,MM
      IF (LU(M,N).GT.0.5) THEN

       UNDIMDEEP=DEEPMIX/HHQ(M,N)
C DEFINITION LAST OF LEVEL VERTICAL INDEX FOR MIXING
       KDEEPMIX=2
       DO WHILE(ZW(KDEEPMIX+1).LE.UNDIMDEEP.AND.KDEEPMIX.LT.NZ)
       KDEEPMIX=KDEEPMIX+1
       ENDDO

C MIXING IF UNSTABLE.

       DEN(M,N,1)=DENP(TT(M,N,1),SS(M,N,1)+35.0,0.0)
       DO K=2,KDEEPMIX
       DEN(M,N,K)=DENP(TT(M,N,K),SS(M,N,K)+35.0,0.0)
        IF( DEN(M,N,K-1).GT.DEN(M,N,K) ) THEN
         DO K1=K,2,-1
          IF( DEN(M,N,K1-1).LT.DEN(M,N,K1) ) GOTO 10
           S1 = 0.0
           S2 = 0.0
           S3 = 0.0

           DO K2=K1-1,K
            S1 = S1 + TT(M,N,K2) * DZ(K2)
            S2 = S2 +              DZ(K2)
            S3 = S3 + SS(M,N,K2) * DZ(K2)
           ENDDO

           TU = S1 / S2
           SU = S3 / S2
           DEN0=DENP(TU,SU+35.0,0.0)
           DO K2=K1-1,K
             DEN(M,N,K2) = DEN0
              TT(M,N,K2) = TU
              SS(M,N,K2) = SU
           ENDDO
        ENDDO
   10   CONTINUE

        ENDIF

       ENDDO

      ENDIF

      ENDDO
      ENDDO

      RETURN
      END
C======================================================================
      SUBROUTINE TSMIXV
C------ UPPER MIXING PROCEDURE, WHICH MIX AS:
C                            Density
C    +--+---+--+--+--+--+--+--+--
C    |                *
C    |                *
C  3 +                 *
C    |                  *--------  KTOP
C    | Unstable        *|      /|\
C  6 + stratification * |       |
C    |               *  |       |
C    |              *   |       |
C  9 +             *    |       |
C    |            *     |       | MIXIN DEPTH
C    |            *     |       |
C 12 +            *     |       |
C    |             *    |       |
C    |             *    |       |
C 15 +              *   |       |
C    |               *  |       |
C    | Stable         * |       |
C 18 + stratification  *|      \|/
C    |                  *--------  KBOT
C    |                   *
C .. +                    *
C    |                     *
C    |                      *
C NZ + Depth
C
      IMPLICIT NONE
      INCLUDE '0COM.INC'
      INCLUDE '0CEAN.INC'
      INTEGER M,N,K,KTOP,KBOT
      REAL SUMT,SUMS,SUMPT,SUM0,ALP
      ALP=1.0
!$OMP PARALLEL DO PRIVATE(M,N,K,KTOP,KBOT,SUMT,SUMS,SUMPT,SUM0)

      DO N=NNN,NN
         DO M=MMM,MM
            IF(LU(M,N).GT.0.5.AND.ABS(YT(N)).LE.LAT_CRIT_4D) THEN
               KTOP=0
               KBOT=0
               DO K=2,NZ
                  IF(DEN(M,N,K).LE.DEN(M,N,K-1)) THEN
                     KTOP=K-1
                     GO TO 10
                  END IF
               END DO
   10          CONTINUE
               IF(KTOP.NE.0) THEN
                  DO K=KTOP+1,NZ
                     IF(DEN(M,N,K).GT.DEN(M,N,KTOP)) THEN
                        GO TO 20
                     ELSE
                        KBOT=K
                     END IF
                  END DO
               END IF
   20          CONTINUE
               IF(KTOP.NE.0.AND.KBOT.GT.KTOP) THEN
                  SUMS=0.0
                  SUMT=0.0
                  SUMPT=0.0
                  SUM0=0.0
                  DO K=KTOP,KBOT
                     SUMS=SUMS+SS(M,N,K)*DZ(K) 
                     SUMT=SUMT+TT(M,N,K)*DZ(K)
                     SUMPT=SUMPT+PASS_TRACER(M,N,K)*DZ(K)
                     SUM0=SUM0+DZ(K)
                  END DO
                  SUMT=SUMT/SUM0
                  SUMS=SUMS/SUM0
                  SUMPT=SUMPT/SUM0
                  DO K=KTOP,KBOT
                     TT(M,N,K)=ALP*SUMT+(1.0-ALP)*TT(M,N,K)
                     SS(M,N,K)=ALP*SUMS+(1.0-ALP)*SS(M,N,K)
                     PASS_TRACER(M,N,K)=
     =                ALP*SUMPT+(1.0-ALP)*PASS_TRACER(M,N,K)
                  END DO
               END IF
            END IF
         END DO
      END DO
!$OMP END PARALLEL DO
      RETURN
      END  
C======================================================================
      SUBROUTINE DEF_TS_UPPER_MIX_LEV(DEN,KTOP,KBOT,LU,NX,NY,NZ)
C------NEW UPPER MIXING PROCEDURE, WHICH MIX AS:
C
C                   Density DEN(NX,NY,NZ)
C    +--+---+--+--+--+--+--+--+--
C    |                *
C    |                *
C  3 +                 *
C    |                  *---------- KTOP(NX,NY)
C    | Unstable        *|      /|\
C  6 + stratification * |       |
C    |               *  |       |
C    |              *   |       |
C  9 +             *    |       |
C    |            *     |       | MIXIN DEPTH
C    |            *     |       |
C 12 +            *     |       |
C    |             *    |       |
C    |             *    |       |
C 15 +              *   |       |
C    |               *  |       |
C    | Stable         * |       |
C 18 + stratification  *|      \|/
C    |                  *---------- KBOT(NX,NY)
C    |                   *
C... +                    *
C    |                     *
C    |                      *
C NZ + Depth
C
      IMPLICIT NONE

	INTEGER NX,NY,NZ

      REAL DEN(NX,NY,NZ),LU(NX,NY)

	INTEGER KTOP(NX,NY),KBOT(NX,NY)

	INTEGER M,N,K

	KTOP=0
	KBOT=0
 
!$OMP PARALLEL DO PRIVATE(M,N,K)

      DO N=1,NY
         DO M=1,NX

            IF(LU(M,N).GT.0.5) THEN

                K=1
	          DO WHILE(K.LE.NZ-1.AND.DEN(M,N,K).LT.DEN(M,N,K+1))
                         K=K+1
	          END DO
              
                KTOP(M,N)=K

	          K=K+1
         
                DO WHILE(K.LE.NZ  .AND.DEN(M,N,K).LE.DEN(M,N,KTOP(M,N)))
                         K=K+1
                END DO

                KBOT(M,N)=MIN(K,NZ)

            END IF
         END DO
      END DO
!$OMP END PARALLEL DO

      RETURN
      END          

C======================================================07.06.07 12:20
      SUBROUTINE RICHNUM_MOD(DEN,UU,VV,TAUX,TAUY,RIT)
      IMPLICIT NONE
C----------------------------------------------------------------------
C RICHARDSON' NUMBER ON U-GRID. (D UV/D Z .NE. 0). K=2,NZ.
C RICHNU = G/R0*(D DEN/D Z) / ((D U/D Z)^2+(D V/D Z)^2), D Z = HZT*HH
C HZT(K) = Z(K)-Z(K-1), K=2,NZ
C RIT - RICHARDSON NUMBER ON T-GRID

      INCLUDE '0COM.INC'
      REAL DEN(NX,NY,NZ),UU(NX,NY,NZ),VV(NX,NY,NZ),
     &     RIT(NX,NY,NZ),TAUX(NX,NY),TAUY(NX,NY)

      REAL R1,R2,U1,U2,V1,V2,U_STAR,RLT
      INTEGER M, N, K
      REAL    CORIOLIS

C FOR MOMENTUM:
      DO N=NNN,NN
       DO M=MMM,MM
        IF (LU(M,N).GT.0.5) THEN
         
          RLT=(RLH(M,N  )+RLH(M-1,N  )
     &        +RLH(M,N-1)+RLH(M-1,N-1))/4.0

          CORIOLIS=MAX(ABS(RLT),2.5E-05)   !FREEZING ON 10N,S
C STORAGE ECMAN DEPTH IN FIST LEVEL OF RIU
          U_STAR = SQRT(SQRT(0.25*(TAUX(M,N)+TAUX(M-1,N))**2
     &                      +0.25*(TAUY(M,N)+TAUY(M,N-1))**2))
          RIT(M,N,1)=0.5*U_STAR/CORIOLIS

          R2 = DEN(M,N,1)

          U2 = LU(M,N)*0.5/DY(M,N)*                   !U ON T-GRID
     &        (UU(M-1,N,1)*HHU(M-1,N)*DYH(M-1,N) +
     &         UU(M  ,N,1)*HHU(M  ,N)*DYH(M  ,N) )/HHQ(M,N)

          V2 = LU(M,N)*0.5/DX(M,N)*                   !V ON T-GRID
     &        (VV(M,N-1,1)*HHV(M,N-1)*DXH(M,N-1) +
     &         VV(M,N  ,1)*HHV(M,N  )*DXH(M,N  ) )/HHQ(M,N)

	   DO K=2,NZ
          R1 = R2
          R2 = DEN(M,N,K)
          U1 = U2
          V1 = V2
          U2 = LU(M,N)*0.5/DY(M,N)*                   !U ON T-GRID
     &        (UU(M-1,N,K)**HHU(M-1,N)*DYH(M-1,N) +
     &         UU(M  ,N,K)**HHU(M  ,N)*DYH(M  ,N) )/HHQ(M,N)

          V2 = LU(M,N)*0.5/DX(M,N)*                   !V ON T-GRID
     &        (VV(M,N-1,K)*HHV(M,N-1)*DXH(M,N-1) +
     &         VV(M,N  ,K)*HHV(M,N  )*DXH(M,N  ) )/HHQ(M,N)

C DIRECT DISCRETE FORM OF Ri
          RIT(M,N,K)=MIN(G*(R2-R1)/2.0*(ZW(K+1)-ZW(K-1))*HHQ(M,N)/
     &                    ((U2-U1)**2+(V2-V1)**2+0.0025),1E+05)
         ENDDO

        ENDIF
       ENDDO
      ENDDO

C     IF(MMD.NE.0) CALL CYCLIZE(RIT,NX,NY,NZ,MMM,MM)

      RETURN
      END
C======================================================================
C OCEAN GENERAL CIRCULATION MODEL PARAMETRIZATIONS OF VERTICAL MIXING
C (C) DIANSKY N.A., 2008.   
C======================================================================
      SUBROUTINE MO_MIX_MODIFIED(DEN,UU,VV,RIT,
     &                        ANZT,ANUMAXT,ANUBGRT,
     &                        ANZU,ANUMAXU,ANUBGRU,TAUX,TAUY,AFNT)
      IMPLICIT NONE
C----------------------------------------------------------------------
*     Vertical Mixing Coefficients Calculation
*     according to Monin-Obukhov (or Kochergin's) parametrizations .
C ANU = ((D U/D Z)^2+(D V/D Z)^2-CONST*G/R0*(D DEN/D Z) )
C HZT(K) = Z(K)-Z(K-1), K=2,NZ

      INCLUDE '0COM.INC'
      INCLUDE '1SEAICE.INC'

      REAL  DEN(NX,NY,NZ  ),    !POTENTIAL DENSYTY
     &       UU(NX,NY,NZ  ),	!ZONAL VELOCITY
     &       VV(NX,NY,NZ  ),    !MERIDIONAL VELOCITY
     &      RIT(NX,NY,NZ  ),    !ARRAY FOR RICHRDSON IS USED FOR STORE 
                              !              VAISJALIJA BRENDTA FREQ
     &     TAUX(NX,NY     ),	!ZONAL WIND STRESS
     &     TAUY(NX,NY     ),  !MERIDIONAL WIND STRESS
     &     ANZT(NX,NY,NZ+1),  !DIFFUSION COEFFICIENT
     &     ANZU(NX,NY,NZ+1),  !VISCOSITY COEFFICIENT
     &     AFNT(NX,NY,NZ+1),
     &     ANUMAXT,ANUBGRT,   !MAX AND MIN DIFFUSION COEFFICIENT
     &     ANUMAXU,ANUBGRU    !MAX AND MIN VISCOSITY COEFFICIENT

      REAL R1,R2,U1,U2,V1,V2,CKro,Q05,VBF,C20,C21,ZZZ,A2,Uz,Vz,DZ21,
     &     ANUPU,ANUPT !ADDITIONAL MIXING IN UPPER LAYER 
      INTEGER M, N, K

	INTEGER, ALLOCATABLE:: KTOP(:,:),KBOT(:,:)   !HELP ARRAYS FOR UPPER
                                                   !AND BOTTOM LEVELS WHERE 
	                                             !MIXING IS NEED
* IAKOVLEV PARAMETRIZATION
*     PARAMETERS FOR COEFFICIENTS CUTOFFS.
      REAL D,SL,DIMDEPTH,UNDIMDEPTH,NU_UNSTABLE
      PARAMETER (DIMDEPTH=300.0,ANUPU=500.0,ANUPT=500.0) !UPPER MIXED LAYER COEFFICIENTS ARE BIGGER
      PARAMETER (NU_UNSTABLE=500.0)


      REAL PIP180,DPIP180              !FOR DEGREES TO RADIANS CONVERS
	PARAMETER(PIP180=3.141592653/180.0,
     &         DPIP180=3.1415926535897/180.0D00)

C--------ADDITIONAL UPPER MIXING PARAMETERS FROM VOLDIN----------------------
	REAL ANUPU2, ANUPT2
      PARAMETER (ANUPU2=500.0,ANUPT2=500.0) !UPPER MIXED LAYER COEFFICIENTS ARE BIGGER

C VOLODIN'S PARAMETRIZATION IN UPPER LAYER	
	REAL STRMAX1,STRMIN1,STRMAX2,STRMIN2,
     &     HMMAX1,HMMAX2,HMMIN1,HMMIN2,DCRT

      PARAMETER(
     &         STRMIN1=0.4,	
     &		  STRMAX1=1.0,	
     &		  STRMIN2=1.0,	
     &		  STRMAX2=8.0,	
     &		  HMMIN1=300.0,	
     &		  HMMAX1=3000.0,	
     &		  HMMIN2=3000.0,
     &		  HMMAX2=10000.0,	
     &		  DCRT=300.0)						
	
	REAL STR,HM,ARG,SANZT,ANZT1(NZ)
	INTEGER KK
      REAL      DPTHMIX
      PARAMETER(DPTHMIX=1.5E-04) !DENSITY CRITERIUM APPROPRIATE HALF
	                           !DEGREE TEMPERATURE JUMP
C----------------------------------------------------------------------


C FOR MOMENTUM:
!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,K,R2,U2,V2,DZ21,R1,U1,V1,VBF,CKro,UNDIMDEPTH,
!$OMP&        Uz,Vz,Q05,STR,HM,ANZT1,SANZT,KK,ARG,D,SL)
      DO N=NNN,NN
       DO M=MMM,MM
        IF (LU(M,N).GT.0.5) THEN

C MIXED LAYER DEPTH ESTIMATION
        K=2

        DO WHILE (DEN(M,N,K)-DEN(M,N,1).LT.DPTHMIX.AND.K.LT.NZ)
         K=K+1
        END DO
        
        D=HHQ(M,N)*(Z(K-1)+Z(K))/6.0
	  
        SL=(D*0.025)**2

          UNDIMDEPTH = DIMDEPTH/HHQ(M,N)

          R2 = DEN(M,N,1)

          U2 = LU(M,N)*0.5/DY(M,N)*                   !U ON T-GRID
     &        (UU(M-1,N,1)*HHU(M-1,N)*DYH(M-1,N) +
     &         UU(M  ,N,1)*HHU(M  ,N)*DYH(M  ,N) )/HHQ(M,N)

          V2 = LU(M,N)*0.5/DX(M,N)*                   !V ON T-GRID
     &        (VV(M,N-1,1)*HHV(M,N-1)*DXH(M,N-1) +
     &         VV(M,N  ,1)*HHV(M,N  )*DXH(M,N  ) )/HHQ(M,N)

	   DO K=2,NZ

          DZ21= HZT(K)*HHQ(M,N)
           R1 = R2
           R2 = DEN(M,N,K)
           U1 = U2
           V1 = V2

          U2 = LU(M,N)*0.5/DY(M,N)*                   !U ON T-GRID
     &        (UU(M-1,N,K)*HHU(M-1,N)*DYH(M-1,N) +
     &         UU(M  ,N,K)*HHU(M  ,N)*DYH(M  ,N) )/HHQ(M,N)

          V2 = LU(M,N)*0.5/DX(M,N)*                   !V ON T-GRID
     &        (VV(M,N-1,K)*HHV(M,N-1)*DXH(M,N-1) +
     &         VV(M,N  ,K)*HHV(M,N  )*DXH(M,N  ) )/HHQ(M,N)


***********Analogous procedure in sigma coordinate**********************

*     Simplified Vasala-Brendt **2 Frequency  
          VBF= G*(R2-R1)/DZ21  !VAISJALIJA BRENDTA FREQUENCY

*     Vertical shift
          Uz = (U2-U1)/DZ21
          Vz = (V2-V1)/DZ21
	    
          RIT(M,N,K)=VBF/(Uz**2+Vz**2+1E-10)       !RICHARDSON NUMBER

          IF(RIT(M,N,K).NE.0.0) THEN
		CKro=0.725*(RIT(M,N,K)+0.186
     &       -SQRT(RIT(M,N,K)**2-0.316*RIT(M,N,K)+0.034596))/RIT(M,N,K)
          ELSE
            CKro=1.34
	    END IF

		Q05= SQRT( AMAX1( 0.0, Uz**2 +Vz**2 - CKro*VBF) )

C FOR MOMENTUM:
		       ANZU(M,N,K)= AMIN1(ANUMAXU,ANUBGRU + SL*Q05)

C FOR TEMPERATURE AND SALINITY:
              IF(RIT(M,N,K).GE.0.0) THEN
		       ANZT(M,N,K)= AMIN1(ANUMAXT,ANUBGRT + CKro*SL*Q05)
              ELSE
		       ANZT(M,N,K)= NU_UNSTABLE
	        END IF

C FOR UPPER LAYER
          IF(ZW(K).LE.UNDIMDEPTH) THEN

C FOR MOMENTUM:
		ANZU(M,N,K)= AMAX1(ANUPU,ANZU(M,N,K))

C FOR TEMPERATURE AND SALINITY:
		ANZT(M,N,K)= AMAX1(ANUPT*AICE0(M,N),ANZT(M,N,K))
       
	    END IF

         ENDDO

C---- ADDITIONAL DIFFUSION DUE TO PROPAGATION OF WAVE FLUCTUATIONS--
C         STR=RH0*SQRT(TAUX(M,N)*TAUX(M,N)+TAUY(M,N)*TAUY(M,N))*
C     *    AICE0(M,N)
c         strout(m,n)=str
C         STR=AMAX1(STR,STRMIN1)
C         STR=AMIN1(STR,STRMAX2)
C         IF(STR.LE.STRMAX1) THEN
C            HM=HMMIN1+(HMMAX1-HMMIN1)*(STR-STRMIN1)/(STRMAX1-STRMIN1)
C         ELSE
C            HM=HMMIN2+(HMMAX2-HMMIN2)*(STR-STRMIN2)/(STRMAX2-STRMIN2)
C         END IF
C         DO K=2,NZ
C            IF(ZW(K)*HHQ(M,N).LE.HM) THEN
C               ANZU(M,N,K)=AMAX1(ANZU(M,N,K),ANUPU2)
C               ANZT(M,N,K)=AMAX1(ANZT(M,N,K),ANUPT2)
C            END IF
C         END DO
C         DO K=2,NZ
C            ANZT1(K)=ANZT(M,N,K)
C         END DO

C---- ADDITIONAL DIFFUSION DUE TO PROPAGATION OF INSTABILITY FLUCTUATIONS 
C         DO K=2,NZ

C            IF(ANZT1(K).LT.100.0.AND.HHQ(M,N)*ZW(K).LE.30000.) THEN
C               SANZT=0.0
C               DO KK=2,NZ
C                  IF(K.NE.KK) THEN
C                     ARG=-ABS(ZW(K)-ZW(KK))*HHQ(M,N)/DCRT
C                     SANZT=AMAX1(SANZT,(ANZT(M,N,KK)-ANUBGRT)*EXP(ARG))
C                  END IF
C               END DO
C               ANZT(M,N,K)=ANZT(M,N,K)+SANZT
C            END IF
C         END DO


        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF(MMD.NE.0) CALL CYCLIZE(ANZU,NX,NY,NZ,MMM,MM) !ENOUGH FOR U,V ONLY



      RETURN
      END
