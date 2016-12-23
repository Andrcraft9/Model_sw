C======================================================================
C (C) ARCTIC@INM.RAS.RU-----------TSICF4N.FOR--------04.05.05
C LH = LH(Q2M,TA,PA)
C SEA ICE MODEL FOR WORLD OCEAN. THERMAL. (NJ'S NOV'01 VERS.)
C---------------------------------------------------------------
C-----------------------------------------------
      SUBROUTINE ICETH(TAU)
C-----------------------------------------------
C     VERSION 07.11.2001.
      IMPLICIT NONE
      INCLUDE '0COM.INC'
      INCLUDE '1SEAICE.INC'
      INCLUDE '0CEAN.INC'
	INCLUDE '2ICE4FLUX.INC'

      REAL TAU,SH,LH,AREA
      INTEGER  K,M,N, MITER, MARK
      REAL TFC,COEFW,Q,CTB,QS,
     &     ACCUR,HI,C1,C2,HS,T1,T2, F2,QSI,QAS,X,QIW,QI,
     &     QDELTA,HOLD,AIOLD,DAICE, STOW,DHSNOW,STOICE,DHICEN,
     &     VOL,DHICE,QAI,TW,SWI
      INCLUDE '0TF.INC'
      INCLUDE '0ICEPAR.INC'

      REAL T00,SALIN,TSRF,EPSs
	REAL HHWK,HICEF
	REAL SENS_HEAT,LAT_HEAT,SATUR_PRESS_ICE,SATUR_HUM

c      P=1.01

      ACCUR= 1E-3
      MITER= 20

	TICE=0.0
	TSNOW=0.0

!$OMP PARALLEL DO PRIVATE(M,N,K,SALIN,TFC,DHSNOW,DHICE,SWI,
!$OMP&     HI,C2,HS,C1,T00,T1,T2,TSRF,EPSs,QS,SH,LH,F2,QSI,QAS,
!$OMP&      QAI,COEFW,X,CTB,TW,QIW,QI,QDELTA,HOLD,AIOLD,DAICE,
!$OMP&      STOW,STOICE,DHICEN,MARK,VOL,HHWK,HICEF,Q,AREA)              
      DO N=NNN,NN
	  DO M=MMM,MM
	          
           IF (LU(M,N).GT.0.5) THEN
C           1. OLD SNOW AND ICE.
             SALIN=SS(M,N,1)+35.0
             TFC= TF(SALIN)

              DO K=1,MGRAD

               DHSNOW =0.
               DHICE  =0.
	         SWI=0.0

                IF( AICE(M,N,K) .GT. AIMIN) THEN  !CASE IF ICE IS AVAILABLE

                   HI= HICE(M,N,K)/AICE(M,N,K)
                   IF(HI.LT.1.0E-6) WRITE(*,*) 'HI',M,N,HI
                   C2= CKICE /HI

C                 SNOWFALL IF AIR TEMPERATURE LESS THEN 0C.
                   IF( TA(M,N) .LT. 0.0) THEN
	             AREA=(TA(M,N)+5.0)/10.0
                   AREA=MAX(AREA,0.0)
                   AREA=MIN(AREA,1.0)
                   HSNOW(M,N,K)=HSNOW(M,N,K)
     &                    +ROW*PR(M,N)*TAU*AICE(M,N,K)*(1.0-AREA)/ROSDRY
                   END IF

C                  1.1. UPPER SURFACE.

C                   1.1.1. THE CASE WITH SNOW COVER.
C                       VOLUMETRIC FUSION HEAT.
            
                   IF( HSNOW(M,N,K) .GT. HSMIN) THEN !CASE IF SNOW IS AVAILABLE
                      HS= HSNOW(M,N,K)/AICE(M,N,K)
                      C1= CKSNOW/HS
                      IF(C1+C2.LT.1.E-8) 
     &                   WRITE(*,*) 'C1+C2',M,N,C1,C2,C1+C2


C                     BISECTION METHOD.
                   T00= -50.
                   T1= 1.0
      
	            DO WHILE((T1-T00) .GT. ACCUR)
                   T2=0.5*(T00+T1)
                   TSRF=T2+273.15

                   IF(TSRF.GT.373.0.OR.TSRF.LT.173.0) THEN
                     WRITE(*,*) 'IN QMAX T =',TSRF
                     STOP
                   END IF

c                  Saturation vapor pressure EPSs:
                   EPSs=SATUR_PRESS_ICE(T2)
c                  Specific humidity at the surface Qs:
                   QS=SATUR_HUM(EPSs,SLPR(M,N))
                   
                   SH=SENS_HEAT(ROA,CPA,CDHI,E0+WIND(M,N),TA(M,N),T2)
                   LH=LAT_HEAT(ROA,QLI,CDHI,E0+WIND(M,N),QA(M,N),QS)

                   F2= SH + LH
     *               +  C1*( (C2*TFC+C1*T2)/(C1+C2) - T2)
     &             + (1.-asdry)*SW(M,N) +Es*LW(M,N) -Es*Sigma*(TSRF**4)
           
                   IF( F2 .LT. 0.) THEN
                    T1=T2
                   ELSE
                    T00=T2
                   END IF
                  END DO 
           
                  TSNOW(M,N,K)= T2

                  IF (TSNOW(M,N,K) .LT. TMELT) THEN
C                   NO SNOW MELTING.
                    TICE(M,N,K)= (C1*TSNOW(M,N,K) +C2*TFC)/(C1+C2)
                  ELSE

C                   SNOW MELTING.
                    TSNOW(M,N,K)= TMELT

C                   NEW ICE TEMPERATURE.
                    TICE(M,N,K)= (C1*TMELT +C2*TFC)/(C1+C2)

C                   HEAT FLUX IN ICE.
                    QSI= C1*(TMELT-TICE(M,N,K))

C            TOTAL HEAT FLUX AT SNOW SURFACE.
C            SPECIFIC HUMIDITY AT THE SURFACE QS:
                    TSRF=TMELT+273.15      
      
                   IF(TSRF.GT.373.0.OR.TSRF.LT.173.0) THEN
                     WRITE(*,*) 'IN QMAX T =',TSRF
                     STOP
                   END IF

c                 Saturation vapor pressure EPSs:
                  EPSs=SATUR_PRESS_ICE(TMELT)
c                 Specific humidity at the surface Qs:
                  Qs=SATUR_HUM(EPSs,SLPR(M,N))
c     LH= (E0+roa*QLi*CDHi*wind)*(Q2m(i,j)-Qs)
                  
                  SH=SENS_HEAT(ROA,CPA,CDHI,E0+WIND(M,N),TA(M,N),TMELT)
                  LH= LAT_HEAT(ROA,QLI,CDHI,E0+WIND(M,N),QA(M,N),QS)

                  QAS= SH + LH - QSI                
     &            + (1.-aswet)*SW(M,N) +Es*LW(M,N) -Es*Sigma*(TSRF**4)

C                 SNOW MASS CHANGE.
                  DHSNOW= AICE(M,N,K)*MIN(-TAU*QAS/QSNOW, 0.)

C     IF THE MASS OF MELTED SNOW IS TOO LARGE, THE REMAINING HEAT
C     QDELTA MELTS UNDERLYING ICE.

C     VOLUMETRIC SPECIFIC HEAT OF FUSION.
C     ICE TEMPERATURE IS SET TO TMELT.
                   IF( HSNOW(M,N,K)+DHSNOW .LT. 0.) THEN
                    DHICE= DHICE + QSNOW*(HSNOW(M,N,K)+DHSNOW)/QICE
                    DHSNOW = -HSNOW(M,N,K)
                    TICE(M,N,K)= TMELT
                   END IF

                   HSNOW(M,N,K)= HSNOW(M,N,K) + DHSNOW

                  END IF ! END OF SNOW MELTING.

           
           SWI=0.0
            
            
                  ELSE  !END OF CASE IF SNOW IS AVAILABLE

C           1.1.2. NO SNOW.

C           BISECTION METHOD.
                   T00=-50.
                   T1=1.0
	    
          DO WHILE((T1-T00) .GT. ACCUR)
                   T2=0.5*(T00+T1)
                   TSRF=T2+273.15
      
              IF(TSRF.GT.373.0.OR.TSRF.LT.173.0) THEN
                 WRITE(*,*) 'IN QMAX T =',TSRF
                 STOP
              END IF
C         QS=QMAX(TSRF,P,0.0E0)
c              Saturation vapor pressure EPSs:
              EPSs=SATUR_PRESS_ICE(T2)
c              Specific humidity at the surface Qs:
              Qs=SATUR_HUM(EPSs,SLPR(M,N))
              
              SH=SENS_HEAT(ROA,CPA,CDHI,E0+WIND(M,N),TA(M,N),T2)
              LH= LAT_HEAT(ROA,QLI,CDHI,E0+WIND(M,N),QA(M,N),QS)

              F2= SH + LH
     *         +C2*( TFC - T2) + (1.-aI0)*(1.-aidry)*SW(M,N)
     &         +Ei*LW(M,N) -Ei*Sigma*(TSRF**4)
          
              IF( F2 .LT. 0.) THEN
               T1=T2
              ELSE
               T00=T2
              END IF
	    END DO

140       TICE(M,N,K)= T2

          SWI = aI0 * (1.-aidry)*SW(M,N)

          IF( TICE(M,N,K) .GE. TMELT) THEN         !ICE MELTING
               
               TICE(M,N,K)= TMELT
C               NEW ENERGY BALANCE AT ICE SURFACE.
C               SATURATION VAPOR PRESSURE EPSS:
               TSRF=TMELT+273.15
	
              IF(TSRF.GT.373.0.OR.TSRF.LT.173.0) THEN
                WRITE(*,*) 'IN QMAX T =',TSRF
                STOP
              END IF
c     Saturation vapor pressure EPSs:
             EPSs=SATUR_PRESS_ICE(TMELT)
c     Specific humidity at the surface Qs:
             Qs=SATUR_HUM(EPSs,SLPR(M,N))
C     SPECIFIC HUMIDITY AT THE SURFACE QS:
              SH=SENS_HEAT(ROA,CPA,CDHI,E0+WIND(M,N),TA(M,N),TMELT)
              LH= LAT_HEAT(ROA,QLI,CDHI,E0+WIND(M,N),QA(M,N),QS)

             QAI= SH + LH 
     &       + C2*(TFC-TMELT)
     &       + (1.-aI0)*(1.-aiwet)*SW(M,N)
     &       + Ei*LW(M,N) -Ei*Sigma*(TSRF**4)

             SWI = aI0 * (1.-aiwet)*SW(M,N)
C             ICE MASS CHANGE.
             DHICE= DHICE+AICE(M,N,K)*MIN(-TAU*QAI/QICE, 0.)
           END IF   !END OF ICE MELTING

         END IF     !  END OF CASE IF SNOW IS UNAVAILABLE

C     1.2. HEAT FLUXES AT THE BASEMENT OF THE ICE.
C          ICE TEMPERATURE IS NOT LOWER TFC.

C       PARAMETRIZATION BY EBERT & CURRY, 1993.
         COEFW=2.

      
         IF(AICE(M,N,K).LT.1.E-8) WRITE(*,*) 'AICE',M,N,AICE(M,N,K)
         
         X= 0.01*HICE(M,N,K)/AICE(M,N,K)
         X= MAX(X,0.01)
	   X= MIN(X,3.00)
         CTB= 1.26E-2*COEFW /SQRT(X)


       TW= MAX(TFC,TT(M,N,1))
      QIW= ROW*CPW8*CTB*(TW-TFC)
      QI = C2*(TICE(M,N,K)-TFC)

C     LOWER ICE SURFACE MASS CHANGE
      DHICE= DHICE -AICE(M,N,K)*TAU*(QIW+QI)/QICE

C     IF TOO MUCH ICE WAS MELTED, THE REMAINING HEAT IS GOING TO
C     THE VERY UPPER OCEAN LAYER.
C     WARNING: THIS IS ESSENTIALLY EXPLICIT TIME SCHEME.

        IF( (HICE(M,N,K)+DHICE) .LT. 0.)THEN
           QDELTA= QICE*(-DHICE-HICE(M,N,K))
C           TT(I,J,1)= TT(I,J,1)+QDELTA/(HHQ(I,J)*(ZW(2)-ZW(1))*ROW*CPW8)
           HEATICE2OC(M,N)=HEATICE2OC(M,N)+QDELTA
           DHICE= -HICE(M,N,K)
           AICE(M,N,K)=0.

C             COMPLETE SNOW MELTING
              DHSNOW= DHSNOW-HSNOW(M,N,K)
              HSNOW(M,N,K)=0.

        END IF

        HOLD= HICE(M,N,K)
        HICE(M,N,K)= HOLD +DHICE

      END IF !END OF CASE IF ICE IS AVAILABLE

C     CHANGE DUE TO SIDE MELTING
       IF(DHICE.LT.0..AND. HICE(M,N,K).GT.1.E-3) THEN
         AIOLD=AICE(M,N,K)
         if(hold.lt.1.e-8) write(*,*) 'hold',M,N,hold
         DAICE=CMELT*AIOLD*DHICE/HOLD
         AICE(M,N,K)= AICE(M,N,K)+DAICE

C             SNOW MASS CHANGE
              IF(HSNOW(M,N,K).GT.1.E-6) THEN
       if(aiold.lt.1.e-8) write(*,*) 'aiold',M,N,aiold
              STOW=DAICE*HSNOW(M,N,K)/AIOLD
              DHSNOW= DHSNOW+STOW
              HSNOW(M,N,K)=HSNOW(M,N,K)+STOW
              END IF

      END IF

      DHSNOWT(M,N) = DHSNOWT(M,N) +DHSNOW
      DHICET(M,N) = DHICET(M,N) +DHICE
      SWICE(M,N)=SWICE(M,N)+AICE(M,N,K)*SWI*TAU
      
      
        IF(HICE(M,N,K)  .LT. HIMIN) THEN
          HICE(M,N,K) = 0.
          AICE(M,N,K) = 0.
          HSNOW(M,N,K)= 0.
        END IF


         IF(HSNOW(M,N,K).LT.HSMIN) HSNOW(M,N,K)=0.

C     SNOW TO ICE CONVERSATION.
      IF(AICE(M,N,K). GT. AIMIN)THEN
        HS=HSNOW(M,N,K)/AICE(M,N,K)
        IF(HS.GT.HSMAX) THEN
          STOICE=GAMMA*TAU*HSNOW(M,N,K)/(1.0+GAMMA*TAU)
          HSNOW(M,N,K)= HSNOW(M,N,K)-STOICE
          HICE(M,N,K) = HICE(M,N,K) +STOICE*ROSDRY/ROI
        END IF
      END IF

      END DO  ! LOOP OVER ICE THICKNESS GRADATIONS

C  2. NEW ICE FORMATION.

      DHICEN= 0.
      MARK= 0
      DO K=1,NZ
         VOL=HHQ(M,N)*DZ(K)
         TFC= TF(SS(M,N,K)+35.0)
         IF( TT(M,N,K) .LT. TFC) THEN
            HHWK=ZW(K+1)*HHQ(M,N)
C     Restriction on Sea Ice Growth Depth
C           HICEF=4.E+8/(HICE(M,N,1)**2) 
C           IF(HHWK.LT.HICEF)  THEN            
            IF(HHWK.LT.5000.0) THEN
               DHICEN=DHICEN + CPW8*ROW*VOL*(TFC-TT(M,N,K))/QICE
               MARK= 1
            END IF
            TT(M,N,K)=TFC
         END IF
      END DO

      IF (MARK .EQ. 1) THEN
C       NEW HEAT CONTENT
         Q= (TICE(M,N,1)+TFC)*HICE(M,N,1)/2. +TFC*DHICEN
C          NEW ICE MASS
           HICE(M,N,1)= HICE(M,N,1) + DHICEN
           DHICET(M,N)= DHICET(M,N) + DHICEN
C            NEW ICE TEMPERATURE
             IF( HICE(M,N,1) .GT. 1.) THEN
               TICE(M,N,1)= Q/HICE(M,N,1)
             ELSE
               TICE(M,N,1)= TFC
             END IF

C     ICE COMPACTNESS CHANGE DUE TO NEW ICE (HIBLER,1979)
      AICE(M,N,1)= AICE(M,N,1)+DHICEN*AICE0(M,N)/HREF
      END IF
           
C           SAICE=0.0
C           DO K=1,MGRAD
C              AICE(M,N,K)=MIN(AICE(M,N,K),1.0E0)
C              AICE(M,N,K)=MAX(AICE(M,N,K),0.0E0)
C              SAICE=SAICE+AICE(M,N,K)
C           END DO
C           AICE0(M,N)=1.0-SAICE
            
            END IF    !     (LU>0.5)
        END DO
	END DO
!$OMP END PARALLEL DO        
 

      IF(MMD.NE.0) THEN
       CALL CYCLIZE(AICE,NX,NY,MGRAD,MMM,MM)
       CALL CYCLIZE(HICE,NX,NY,MGRAD,MMM,MM)
	 CALL CYCLIZE(HSNOW,NX,NY,MGRAD,MMM,MM)
C      CALL CYCLIZE(AICE0,NX,NY,1,MMM,MM)
      END IF

C     AICE0=MAX(AICE0,0.0)
C     AICE0=MIN(AICE0,1.0)

      
      RETURN
      END
C======================================================================
      SUBROUTINE TRX2D_DD(FF,UU,TAU,MGRAD)
      IMPLICIT NONE
C----------------------------------------------------------------------
C X-TRANSPORT AND DIFFUSION ON T - GRID.
C X-TRANPORT SUPPOSES BOUNDARY VELOCITIES & FI TO BE ZERO.
      INCLUDE '0COM.INC'
C IMPLICIT - EXPLICIT CONTROL TRANSPORT T & S TIME SCHEME:
      INCLUDE '0TRSIM.INC'
      
      INTEGER MGRAD
      REAL FF(NX,NY,MGRAD),UU(NX,NY)
      REAL A(NX),B(NX),C(NX),ETA(NX),RKSI(NX)
      INTEGER M, N, K, IG, II, JJ
      INTEGER M1, M9, MLOOP
      REAL    BP, DP, DM, PP, PM
      REAL    TAU,DISCONT

      DO K=1,MGRAD

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,BP,DP,DM,PP,PM,A,B,C,ETA,RKSI,
!$OMP&       IG,II,JJ,MLOOP,M1,M9,DISCONT)      
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

          BP = DX(M,N)*DY(M,N)*RN / TAU

          PP = UU(M ,N)*DYH(M ,N)/2.0
          PM = UU(M1,N)*DYH(M1,N)/2.0

          DP = ABS(PP)
          DM = ABS(PM)
          DISCONT= PP - PM

C LAST "2" ABOVE IS DUE TO APROXIMATION IN SPACE.
          C(M) =  PP * TRTS_IMP - DP * TRTS_IMP
          A(M) = -PM * TRTS_IMP - DM * TRTS_IMP
          B(M) =  BP + DP * TRTS_IMP + DM * TRTS_IMP
     &          + DISCONT * TRTS_IMP
          ETA(M) = BP * FF(M,N,K) - PP * TRTS_EXP * FF(M9,N,K)
     &                            + PM * TRTS_EXP * FF(M1,N,K)
     &                 + DP * TRTS_EXP *(FF(M9,N,K)-FF(M ,N,K))
     &                 - DM * TRTS_EXP *(FF(M ,N,K)-FF(M1,N,K))
     &            - DISCONT * TRTS_EXP * FF(M ,N,K)
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
!$OMP END PARALLEL DO      
      
      ENDDO


      RETURN
      END
C======================================================================
      SUBROUTINE TRY2D_DD(FF,VV,TAU,MGRAD)
      IMPLICIT NONE
C----------------------------------------------------------------------
C Y-TRANSPORT AND DIFFUSION.
C Y-TRANPORT SUPPOSES BOUNDARY VELOCITIES & FI TO BE ZERO.
C WATER BOUNDARIES AND SOLID BOUNDARIES ARE DIFFERENT.
      INCLUDE '0COM.INC'
C IMPLICIT - EXPLICIT CONTROL TRANSPORT T & S TIME SCHEME:
      INCLUDE '0TRSIM.INC'
      INTEGER MGRAD
      REAL FF(NX,NY,MGRAD),VV(NX,NY)
      REAL    A(NY),B(NY),C(NY),ETA(NY),RKSI(NY)
      INTEGER M, N, K, IG, II, JJ
	REAL    BP, DP, DM, PP, PM
      REAL    TAU,DISCONT

      DO K=1,MGRAD

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,BP,DP,DM,PP,PM,A,B,C,ETA,RKSI,
!$OMP&         IG,II,JJ,DISCONT)
       DO M=MMM,MM
        DO IG=1,LRY(M)
         II=IIY(IG,M)
         JJ=JJY(IG,M)
         DO N=II,JJ
          BP = DX(M,N)*DY(M,N)*RN /TAU

          PP =  VV(M,N  )*DXH(M,N  )/2.0
          PM =  VV(M,N-1)*DXH(M,N-1)/2.0

          DP = ABS(PP)
          DM = ABS(PM)
          DISCONT= PP - PM

C LAST "2" ABOVE IS DUE TO APROXIMATION IN SPACE.
          C(N) =  PP * TRTS_IMP - DP * TRTS_IMP
          A(N) = -PM * TRTS_IMP - DM * TRTS_IMP
          B(N) =  BP + DP * TRTS_IMP + DM * TRTS_IMP
     &          + DISCONT * TRTS_IMP
          ETA(N) = BP * FF(M,N,K)  - PP * TRTS_EXP * FF(M,N+1,K)
     &                             + PM * TRTS_EXP * FF(M,N-1,K)
     &                 + DP * TRTS_EXP *(FF(M,N+1,K)-FF(M,N  ,K))
     &                 - DM * TRTS_EXP *(FF(M,N  ,K)-FF(M,N-1,K))
     &            - DISCONT * TRTS_EXP * FF(M,N  ,K)
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
!$OMP END PARALLEL DO      
      
      END DO

      RETURN
      END
C---------------------------------------------------------------
C IAKOVLEV NG. VER. 1.2 02.05.02. IDEA BY I. POLYAKOV.
C SEA ICE THICKNESS REDISTRIBUTION WITH T,S RECALCULATIONS.
      SUBROUTINE REDIS()

      IMPLICIT NONE
      INCLUDE '0COM.INC'
      INCLUDE '1SEAICE.INC'
      INCLUDE '0ICEPAR.INC'
      INTEGER M,N,GRAD

!$OMP PARALLEL DO PRIVATE(M,N,GRAD)      
      DO N=NNN,NN
        DO M=MMM,MM
          IF(LU(M,N).GT.0.5) THEN
c     First gradation water is eliminated by converging
c     in a case of negative open water concentration.
C            IF( AICE0(M,N).LT. 0.) THEN
C             AICE(M,N,1)= AICE0(M,N)+AICE(M,N,1)
C             AICE0(M,N)= 0.
C            END IF

c     Concentration deficit is passed to next thicker category.
            DO GRAD=1, MGRAD-1
             IF( AICE(M,N,GRAD).LT.0.) THEN

              AICE(M,N,GRAD+1)=AICE(M,N,GRAD+1)+AICE(M,N,GRAD)
              AICE(M,N,GRAD)= 0.

              HICE(M,N,GRAD+1)=HICE(M,N,GRAD+1)+HICE(M,N,GRAD)
              HICE(M,N,GRAD)= 0.
             
              HSNOW(M,N,GRAD+1)=HSNOW(M,N,GRAD+1)+HSNOW(M,N,GRAD)
              HSNOW(M,N,GRAD)= 0.
            
             END IF
            END DO

c     Reestablishing thickness distribution.

c     To next thicker category.
            DO GRAD=1,MGRAD-1
             IF(HICE(M,N,GRAD).GT.(HMAX(GRAD)*AICE(M,N,GRAD))) THEN
              
              AICE(M,N,GRAD+1)=AICE(M,N,GRAD+1)+AICE(M,N,GRAD)
              AICE(M,N,GRAD)= 0.

              HICE(M,N,GRAD+1)=HICE(M,N,GRAD+1)+HICE(M,N,GRAD)
              HICE(M,N,GRAD)= 0.

              HSNOW(M,N,GRAD+1)=HSNOW(M,N,GRAD+1)+HSNOW(M,N,GRAD)
              HSNOW(M,N,GRAD)= 0.

             END IF
            END DO

c     To next thinner category.
            DO GRAD=MGRAD,2,-1
             IF( HICE(M,N,GRAD) .LE. HMAX(GRAD-1)*AICE(M,N,GRAD)) THEN
              
              AICE(M,N,GRAD-1)=AICE(M,N,GRAD-1)+AICE(M,N,GRAD)
              AICE(M,N,GRAD)= 0.
              
              HICE(M,N,GRAD-1)=HICE(M,N,GRAD-1)+HICE(M,N,GRAD)
              HICE(M,N,GRAD)= 0.
              
              HSNOW(M,N,GRAD-1)=HSNOW(M,N,GRAD-1)+HSNOW(M,N,GRAD)
              HSNOW(M,N,GRAD)= 0.

             END IF
            END DO


c     Small floes cut off.
C            DO GRAD=1,MGRAD
C             IF(HICE(M,N,GRAD).LT.1E-6.OR.AICE(M,N,GRAD).LT.1E-3) THEN
C              AICE(M,N,GRAD)=0.
C              HICE(M,N,GRAD)=0.
C              HSNOW(M,N,GRAD)=0.
C             END IF
C            END DO

            END IF
          END DO
        END DO
!$OMP END PARALLEL DO
      RETURN
      END
C======================================================================
      SUBROUTINE TRXU2D_DD(FF,UU,TAU)
      IMPLICIT NONE
C----------------------------------------------------------------------
C X-TRANSPORT AND DIFFUSION OF FF ON U-GRID.
C X-TRANPORT SUPPOSES BOUNDARY VELOCITIES & FI TO BE ZERO.
      INCLUDE '0COM.INC'
      INCLUDE '0TRSIV.INC'
      REAL TAU,DISCONT
      REAL FF(NX,NY),UU(NX,NY)
      REAL A(NX),B(NX),C(NX),ETA(NX),RKSI(NX)
      INTEGER M, N, IG, IGRX
      INTEGER II, JJ, JJP, MLOOP
      INTEGER M1, M9
      REAL    BP, DP, DM, PM, PP

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,IG,II,JJ,JJP,MLOOP,M1,M9,
!$OMP&        BP,DP,DM,PP,PM,A,B,C,ETA,RKSI,DISCONT)
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

          BP = DXT(M,N)*DYH(M,N)*RN / TAU

          PP = ( UU(M ,N)*DYH(M ,N) +
     &           UU(M9,N)*DYH(M9,N) ) /4.

          PM = ( UU(M1,N)*DYH(M1,N) +
     &           UU(M ,N)*DYH(M ,N) ) /4.

          DP = ABS(PP)
          DM = ABS(DM)

          DISCONT=- (PP - PM)

C LAST "2" ABOVE IS DUE TO APROXIMATION IN SPACE.
          C(M) =  PP * TRUV_IMP - DP * TRUV_IMP
          A(M) = -PM * TRUV_IMP - DM * TRUV_IMP
          B(M) =  BP + DP * TRUV_IMP + DM * TRUV_IMP
     &          + DISCONT * TRUV_IMP
          ETA(M) = BP * FF(M,N) - PP * TRUV_EXP * FF(M9,N)
     &                          + PM * TRUV_EXP * FF(M1,N)
     &                 + DP * TRUV_EXP *(FF(M9,N)-FF(M ,N))
     &                 - DM * TRUV_EXP *(FF(M ,N)-FF(M1,N))
     &            - DISCONT * TRUV_EXP * FF(M ,N)
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
           FF(M,N) = RKSI(M)
          ENDDO

	  ENDDO
       ENDDO
!$OMP END PARALLEL DO
      
      RETURN
      END
C======================================================================
      SUBROUTINE TRYU2D_DD(FF,VV,TAU,IGRY)
      IMPLICIT NONE
C------------------------------------------------------09-13-96 06:08pm
C Y-TRANSPORT AND DIFFUSION OF FF ON U-GRID.
C Y-TRANPORT SUPPOSES BOUNDARY VELOCITIES & FI TO BE ZERO.
      INCLUDE '0COM.INC'
      INCLUDE '0TRSIV.INC'
      REAL TAU,DISCONT
	REAL FF(NX,NY),VV(NX,NY)
      REAL A(NY),B(NY),C(NY),ETA(NY),RKSI(NY)
	INTEGER M, N, IG, IGRY
      INTEGER II, JJ
      REAL    BP, DP, DM, PM, PP

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,IG,II,JJ,
!$OMP&        BP,DP,DM,PP,PM,A,B,C,ETA,RKSI,DISCONT)
       DO  M=MMM,MM
        DO IG=1,LRYU(M)
         II=IIYU(IG,M)
         JJ=JJYU(IG,M)+1
         DO  N=II,JJ

          BP =  DXT(M,N)*DYH(M,N)*RN/ TAU

          PP =( VV(M  ,N  )*DXH(M  ,N  )
     &         +VV(M+1,N  )*DXH(M+1,N  ) )/4.0

          PM =( VV(M  ,N-1)*DXH(M  ,N-1)
     &         +VV(M+1,N-1)*DXH(M+1,N-1) )/4.0

          DP = ABS(PP)
          DM = ABS(PM)
          DISCONT= - (PP - PM)          

C LAST "2" ABOVE IS DUE TO APROXIMATION IN SPACE.
          C(N) =  PP * TRUV_IMP - DP * TRUV_IMP
          A(N) = -PM * TRUV_IMP - DM * TRUV_IMP
          B(N) =  BP + DP * TRUV_IMP + DM * TRUV_IMP
     &          + DISCONT * TRUV_IMP
          ETA(N) = BP * FF(M,N)  - PP * TRUV_EXP * FF(M,N+1)
     &                           + PM * TRUV_EXP * FF(M,N-1)
     &                 + DP * TRUV_EXP *(FF(M,N+1)-FF(M,N  ))
     &                 - DM * TRUV_EXP *(FF(M,N  )-FF(M,N-1))
     &            - DISCONT * TRUV_EXP * FF(M,N  )

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

	   CALL FACTOR2(NY,A,B,C,ETA,RKSI,II,JJ)

	   DO  N=II,JJ
          FF(M,N) = RKSI(N)
	   END DO

	  END DO
	 END DO
!$OMP END PARALLEL DO

      RETURN
      END
C======================================================================
      SUBROUTINE TRXV2D_DD(FF,UU,TAU,IGRX)
      IMPLICIT NONE
C----------------------------------------------------------------------
C X-TRANSPORT AND DIFFUSION OF FF ON U-GRID.
C X-TRANPORT SUPPOSES BOUNDARY VELOCITIES & FI TO BE ZERO.
      INCLUDE '0COM.INC'
      INCLUDE '0TRSIV.INC'
      REAL TAU,DISCONT
	REAL FF(NX,NY),UU(NX,NY)
      REAL A(NX),B(NX),C(NX),ETA(NX),RKSI(NX)
      INTEGER M, N, IG, IGRX
      INTEGER II, JJ, MLOOP
      REAL    BP, DP, DM, PM, PP
      INTEGER M1, M9

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,IG,II,JJ,MLOOP,M1,M9,
!$OMP&        BP,DP,DM,PP,PM,A,B,C,ETA,RKSI)
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

          BP = DXH(M,N)*DYT(M,N)*RN / TAU

          PP =( UU(M ,N  )*DYH(M  ,N  )
     &         +UU(M ,N+1)*DYH(M  ,N+1))/4.0

          PM =( UU(M1,N  )*DYH(M1,N  )
     &         +UU(M1,N+1)*DYH(M1,N+1) )/4.0
          
          DP=ABS(PP)
          DM=ABS(PM)
          DISCONT= - (PP - PM)

C LAST "2" ABOVE IS DUE TO APROXIMATION IN SPACE.
          C(M) =  PP * TRUV_IMP - DP * TRUV_IMP
          A(M) = -PM * TRUV_IMP - DM * TRUV_IMP
          B(M) =  BP + DP * TRUV_IMP + DM * TRUV_IMP
     &          + DISCONT * TRUV_IMP
          ETA(M) = BP * FF(M,N) - PP * TRUV_EXP * FF(M9,N)
     &                          + PM * TRUV_EXP * FF(M1,N)
     &                 + DP * TRUV_EXP *(FF(M9,N)-FF(M ,N))
     &                 - DM * TRUV_EXP *(FF(M ,N)-FF(M1,N))
     &            - DISCONT * TRUV_EXP * FF(M ,N)
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
          FF(M,N) = RKSI(M)
         ENDDO

	  ENDDO
       ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
C======================================================================
      SUBROUTINE TRYV2D_DD(FF,VV,TAU)
      IMPLICIT NONE
C------------------------------------------------------09-13-96 06:08pm
C Y-TRANSPORT AND DIFFUSION OF FF ON V-GRID.
C Y-TRANPORT SUPPOSES BOUNDARY VELOCITIES & FI TO BE ZERO.
      INCLUDE '0COM.INC'
	INCLUDE '0TRSIV.INC'
      REAL TAU,DISCONT
      REAL FF(NX,NY),VV(NX,NY)
      REAL A(NY),B(NY),C(NY),ETA(NY),RKSI(NY)
	INTEGER M, N, IG, IGRY
      INTEGER II, JJ
      REAL    BP, DP, DM, PM, PP

!$OMP PARALLEL DO
!$OMP&PRIVATE(M,N,IG,II,JJ,
!$OMP&        BP,DP,DM,PP,PM,A,B,C,ETA,RKSI)
       DO M=MMM,MM
        DO IG=1,LRY(M)
         II=IIY(IG,M)
         JJ=JJY(IG,M)-1
         DO N=II,JJ
          BP = DXH(M,N)*DYT(M,N)*RN/TAU       	

          PP =( VV(M,N+1)*DXH(M,N+1)
     &         +VV(M,N  )*DXH(M,N  ) )/4.0
          
          PM =( VV(M,N-1)*DXH(M,N-1)
     &         +VV(M,N  )*DXH(M,N  ) )/4.0

          DP=ABS(PP)
          DM=ABS(PM)

          DISCONT= - (PP - PM)

C LAST "2" ABOVE IS DUE TO APROXIMATION IN SPACE.
          C(N) =  PP * TRUV_IMP - DP * TRUV_IMP
          A(N) = -PM * TRUV_IMP - DM * TRUV_IMP
          B(N) =  BP + DP * TRUV_IMP + DM * TRUV_IMP
     &          + DISCONT * TRUV_IMP
          ETA(N) = BP * FF(M,N)  - PP * TRUV_EXP * FF(M,N+1)
     &                           + PM * TRUV_EXP * FF(M,N-1)
     &                 + DP * TRUV_EXP *(FF(M,N+1)-FF(M,N  ))
     &                 - DM * TRUV_EXP *(FF(M,N  )-FF(M,N-1))
     &            - DISCONT * TRUV_EXP * FF(M,N  )

	   END DO
C IMPLICIT: A(N)=-DM-PM;C(N)=-DP+PP;B(N)=BP+DM+DP;ETA(N)=BP*FF(M,N,K)
C IF (N=II) : PM = 0, A = -DM ; IF (N=JJ) : PP = 0, C = -DP.
C IF(IGRY.EQ.1) : POINTS OUTSIDE U-GRID ARE TAKEN, SO THEY MUST BE = 0,
C OR U', FI = 0. SOURCE IN ETA APPEAR.


         CALL FACTOR2(NY,A,B,C,ETA,RKSI,II,JJ)
	   DO N=II,JJ
          FF(M,N) = RKSI(N)
         END DO
        END DO
       END DO
!$OMP END PARALLEL DO
      RETURN
      END
C==============================================================
CCCC  SUBROUTINE FOR COMPUTING ICE VELOCITY DUE TO WIND AND WATER STRESS
      SUBROUTINE ICE_WIND_WATER_DRIFT(UICE,VICE,CAI,CWI,TAU,
     &           UWND,VWND,UU1,VV1,AISTOT,MISTOT)
	IMPLICIT NONE
	INCLUDE '0COM.INC'
      REAL UICE(NX,NY),VICE(NX,NY),CAI(NX,NY),CWI(NX,NY),TAU
	REAL UWND(NX,NY),VWND(NX,NY),UU1(NX,NY),VV1(NX,NY)
	REAL AISTOT(NX,NY),MISTOT(NX,NY)
	INTEGER M,N

!$OMP PARALLEL DO
	DO N=NNN,NN
        DO M=MMM,MM
          
          IF(LCU(M,N).GT.0.5) THEN 
C           IF((MISTOT(M,N)+MISTOT(M+1,N)).GT.20.0) THEN

            UICE(M,N)=(UICE(M,N)
     &      +TAU*((CAI(M,N)*UWND(M,N)+CAI(M+1,N)*UWND(M+1,N))/2.0
     &           +(CWI(M,N)+CWI(M+1,N))/2.0*UU1(M,N)))
     &        /(1.0+TAU*((CAI(M,N)+CAI(M+1,N))/2.0
     &                  +(CWI(M,N)+CWI(M+1,N))/2.0))
            
C            END IF
          END IF

          IF(LCV(M,N).GT.0.5) THEN
C           IF((MISTOT(M,N)+MISTOT(M,N+1)).GT.20.0) THEN
      
            VICE(M,N)=(VICE(M,N)
     &      +TAU*((CAI(M,N)*VWND(M,N)+CAI(M,N+1)*VWND(M,N+1))/2.0
     &           +(CWI(M,N)+CWI(M,N+1))/2.0*VV1(M,N)))
     &        /(1.0+TAU*((CAI(M,N)+CAI(M,N+1))/2.0
     &                  +(CWI(M,N)+CWI(M,N+1))/2.0))
            
C            END IF
          END IF

	  END DO
	END DO
!$OMP END PARALLEL DO
        IF(MMD.NE.0) THEN
         CALL CYCLIZE (UICE,NX,NY,1,MMM,MM)
         CALL CYCLIZE (VICE,NX,NY,1,MMM,MM)
	  END IF

	RETURN
	END
C==============================================================
CCCC  SUBROUTINE FOR COMPUTING ICE VELOCITY DUE TO CAVITATING LIQUID RHEOLOGY
	SUBROUTINE ICE_CAV_LIQ_RHEOL(UICE,VICE,PICE,TAU,AISTOT,MISTOT,CX)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
      REAL UICE(NX,NY),VICE(NX,NY),PICE(NX,NY),TAU
	REAL AISTOT(NX,NY),MISTOT(NX,NY)
	REAL PITOT(NX,NY),DIVVEL,HEV,CX
      INTEGER M,N
      
	DIVVEL=0.0
	PITOT=0.0

!$OMP PARALLEL DO PRIVATE(M,N,DIVVEL,HEV)
	DO N=NNN,NN
        DO M=MMM,MM
          IF(LU(M,N).GT.0.5) THEN
                DIVVEL= UICE(M,N)*DYH(M,N)-UICE(M-1,N)*DYH(M-1,N)
     &          +       VICE(M,N)*DXH(M,N)-VICE(M,N-1)*DXH(M,N-1)

                HEV=-SIGN(0.5,DIVVEL)+0.5

	          PITOT(M,N)=PICE(M,N)*HEV
          END IF
        END DO
	END DO
!$OMP END PARALLEL DO       

       IF(MMD.NE.0) THEN
        CALL CYCLIZE(PITOT,NX,NY,1,MMM,MM)
	 END IF
      
c	 CALL FILTER4_TGR_2D(PITOT,CX,1)

C     CALCULATING NEW UICE 
     
!$OMP PARALLEL DO PRIVATE(M,N)      
      DO N=NNN,NN
	  DO M=MMM,MM
          IF(LCU(M,N).GT.0.5) THEN
	     
C           IF((MISTOT(M,N)+MISTOT(M+1,N)).GT.20.0) THEN
            UICE(M,N)=UICE(M,N)
     &           -2.0/(MISTOT(M,N)+MISTOT(M+1,N)+20.0)
     &                *(PITOT(M+1,N)-PITOT(M,N))/(DXT(M,N)*RN)*TAU
C           END IF

          END IF

C     CALCULATING NEW VICE      

          IF(LCV(M,N).GT.0.5) THEN
	     
C           IF((MISTOT(M,N)+MISTOT(M,N+1)).GT.20.0) THEN
            VICE(M,N)=VICE(M,N)
     &           -2.0/(MISTOT(M,N)+MISTOT(M,N+1)+20.0)
     &                *(PITOT(M,N+1)-PITOT(M,N))/(DYT(M,N)*RN)*TAU
C           END IF

          END IF
        END DO
	END DO
!$OMP END PARALLEL DO 
        IF(MMD.NE.0) THEN
         CALL CYCLIZE (UICE,NX,NY,1,MMM,MM)
         CALL CYCLIZE (VICE,NX,NY,1,MMM,MM)
	  END IF

	RETURN 
	END
C==============================================================
CCCC  SUBROUTINE FOR COMPUTING ICE VELOCITY DUE TO SEA LEVEL & PRESSURE GRADIENTS
      SUBROUTINE ICE_SLPRES_GRAD(UICE,VICE,SLHGRX,SLHGRY,DINX,DINY,
     &                           AISTOT,MISTOT,TAU)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
      REAL UICE(NX,NY),VICE(NX,NY),TAU
	REAL AISTOT(NX,NY),MISTOT(NX,NY)
      REAL(8) SLHGRX(NX,NY),SLHGRY(NX,NY),DINX(NX,NY),DINY(NX,NY)
      INTEGER M,N

!$OMP PARALLEL DO PRIVATE(M,N)
      DO N=NNN,NN
	  DO M=MMM,MM
          IF(LCU(M,N).GT.0.5) THEN
C          IF((MISTOT(M,N)+MISTOT(M+1,N)).GT.20.0) THEN
            
            UICE(M,N)=UICE(M,N)+SNGL(SLHGRX(M,N)+DINX(M,N))*TAU
           
C          END IF
          END IF

          IF(LCV(M,N).GT.0.5) THEN
C          IF((MISTOT(M,N)+MISTOT(M,N+1)).GT.20.0) THEN
            
            VICE(M,N)=VICE(M,N)+SNGL(SLHGRY(M,N)+DINY(M,N))*TAU
           
C          END IF
          END IF
         END DO
	END DO
!$OMP END PARALLEL DO

        IF(MMD.NE.0) THEN
         CALL CYCLIZE (UICE,NX,NY,1,MMM,MM)
         CALL CYCLIZE (VICE,NX,NY,1,MMM,MM)
	  END IF	

      RETURN
	END
C--arctica@inm.ras.ru---------Started   04.01.08-----------M_RHEVP.FOR--
C EVP: idea by E. Hunke ea., prog. by N. Iakovlev,             12.03.08
C modif. Bagno A. on C-grid, spheric coord.sys.             
      SUBROUTINE RHEOLOGY_EVP (PICE, dte,Tdamp, AISTOT,MISTOT)
c dte - internal t-steps < Tdamp - elas. wave param. < TAU
c  Parameters for EVP integration :  NICESTEP = 60  !  NJ: dt = 1 Hr, Nst = 60. 
c  Tdamp = TAU / 3.    !  elast. waves  ||  dte   = TAU / REAL(NICEstep)  ++++
c  PARAMETER (C0=20.0, Pcr= 0.5e5, extr2=4., dmin=1.e-11)  ; warn: Pcr .NE. CFR.
c-----------------------------------------------------------------------
c Force caused by rheology stresses. Ver. 4.0 21.05.2004 + Later?.
c EVP sea ice rheology based on CICE v.3.0 time scheme.
c Regularization Nreg={1,2} = by {Hunke, Harder}.  Here Nreg = 1.
c delta_min ~ 2-5 e-9 sec. - just inform (=~div).
c
      IMPLICIT NONE
      INCLUDE '0COM.INC'
      INCLUDE '1SEAICE.INC'

      REAL AISTOT(NX,NY), MISTOT(NX,NY),PICE(NX,NY),
C     &     S11(NX,NY),S22(NX,NY),   !!!   UICE(NX,NY),VICE(NX,NY),
     &     dte, Tdamp, Pmean_XX, Pmean_XY, Fu,Fv, extr2, dmin,
     &     delta_XX, delta_XY, det1,det2, C_Damp_XX,C_Damp_XY,SLU
      INTEGER M,N,Nreg
	REAL e11(NX,NY),e22(NX,NY),e12(NX,NY)
      PARAMETER (Nreg=1,  extr2=4., dmin=1.e-11 )

      E11  =0.
      E22  =0.
      E12  =0.
CC        !!!!!   10:29 12.03.08   ++++++++++++
c Damping constant 615 kg/m2 ~ 61 cm.

c AB: e11 + e22 =e1=Div, e2=Tens(a-div), sigma_(11,22,1,2) in (Pice,T) pts.
c     e12, sigma_12 - (Shear, a-Rotor) - in (P) - 9 pts big shablon.
c Warn: outside LCU and LCV must be u_i=0 and v_i=0.
c Note: boundary pts can be done with 2nd b.c.
c

	det1 = 1.+ 0.5*dte/Tdamp
	det2 = 1.+ 0.5*dte*extr2/Tdamp


!$OMP PARALLEL DO PRIVATE(M,N)
      DO N=NNN,NN                       
        DO M=MMM,MM
C  CALCULATING DEFORMATION RATES E11,E22 ON LU-GRID
         IF(LU(M,N).GT.0.5) THEN     
C           IF(AISTOT(M,N).GE.0.01) THEN
           
           e11(M,N)= 1.0/DX(M,N)/RN*( uice(M,N)-uice(M-1,N) +
     &     (vice(M,N)+vice(M,N-1)) /2.0 * (DXH(M,N)-DXH(M,N-1))/DY(M,N))

           e22(M,N)= 1.0/DY(M,N)/RN*( Vice(M,N)-Vice(M,N-1) +
     &     (Uice(M,N)+Uice(M-1,N)) /2.0 * (DYH(M,N)-DYH(M-1,N))/DX(M,N))              
             
C	     END IF
         END IF

C  CALCULATING DEFORMATION RATES E12 ON LUU-GRID
         IF(LUH(M,N).GT.0.5) THEN     
C          IF((AISTOT(M  ,N  )+AISTOT(M+1,N  )
C     &        +AISTOT(M  ,N+1)+AISTOT(M+1,N+1)).GE.0.04) THEN
           
           e12(M,N)= (                   DXB(M  ,N  )/DYB(M  ,N  )
     &     *(UICE(M  ,N+1)/DXT(M  ,N+1)-UICE(M  ,N  )/DXT(M  ,N  ))    
     &                                  +DYB(M  ,N  )/DXB(M  ,N  )
     &     *(VICE(M+1,N  )/DYT(M+1,N  )-VICE(M  ,N  )/DYT(M  ,N  ))  
     &                                                     )/(RN*2.0)           
C          END IF
         END IF

        END DO
	END DO
!$OMP END PARALLEL DO
      
	IF(MMD.NE.0) THEN
          CALL CYCLIZE(E11 ,NX,NY, 1,  MMM,MM)
          CALL CYCLIZE(E22 ,NX,NY, 1,  MMM,MM)
	    CALL CYCLIZE(E12 ,NX,NY, 1,  MMM,MM)
      END IF

!$OMP PARALLEL DO PRIVATE(M,N,C_Damp_XX,delta_XX,Pmean_XX,
!$OMP&                        C_Damp_XY,delta_XY,Pmean_XY,SLU)
      
   
      DO N=NNN,NN                     
        DO M=MMM,MM
C   SOLVING EVOLUTIONAL EQUATION FOR STRESS COMPONENTS SIGMA1,SIGMA2 ON LU-GRID           
         IF(LU(M,N).GT.0.5) THEN     
C           IF(AISTOT(M,N).GE.0.01) THEN

	       C_Damp_XX = 6150.*Tdamp* RN **2 *DX(M,N)*DY(M,N) /dte**2
              delta_XX =(e11(M,N)+e22(M,N))**2
     &            +(    (e11(M,N)-e22(M,N))**2 
     &                 +(e12(M  ,N  )**2+e12(M  ,N-1)**2
     &                  +e12(M-1,N  )**2+e12(M-1,N-1)**2) )/extr2
              delta_XX = SQRT(delta_XX)
              Pmean_XX = Pice(M,N)
c
c Regularization: new Pice, delta.--------------------------------
              if(nreg .EQ. 1) then
                 Pmean_XX = MIN (Pmean_XX,C_Damp_XX*delta_XX)       ! El. Hunke
                 delta_XX = MAX (dmin,delta_XX)              
C               else
C                Pmean = Pmean*delta/(delta+dmin)       ! M. Harder.
C               delta=delta+dmin                       
              end if
c Stress tensor components. sig_1 = sig_D, sig_2 = sig_T.
c Implicit dte - time steps for sigma_ij.
c
              sigma1(M,N) = (  sigma1(M,N) + 0.5*dte*Pmean_XX
     &         *((e11(M,N)+e22(M,N))/delta_XX-1.0)/Tdamp    )/det1
c
              sigma2(M,N) = (  sigma2(M,N) + 0.5*dte*Pmean_XX
     &          *(e11(M,N)-e22(M,N))/delta_XX/Tdamp        )/det2	

C           END IF
         END IF

C   SOLVING EVOLUTIONAL EQUATION FOR STRESS COMPONENT SIGMA12,SIGMA2 ON LUU-GRID   
         IF(LUH(M,N).GT.0.5) THEN     
C           IF((AISTOT(M  ,N  )+AISTOT(M+1,N  )
C     &        +AISTOT(M  ,N+1)+AISTOT(M+1,N+1)).GE.0.04) THEN
              SLU=LU(M  ,N  )+LU(M+1,N  )+LU(M  ,N+1)+LU(M+1,N+1)+1E-6
	        C_Damp_XY = 6150.*Tdamp* RN **2 *DXB(M,N)*DYB(M,N) /dte**2
               
               delta_XY =  (  (e11(M  ,N  )+e22(M  ,N  ))**2*LU(M  ,N  )
     &                       +(e11(M+1,N  )+e22(M+1,N  ))**2*LU(M+1,N  )
     &                       +(e11(M  ,N+1)+e22(M  ,N+1))**2*LU(M  ,N+1)
     &                       +(e11(M+1,N+1)+e22(M+1,N+1))**2*LU(M+1,N+1) 
     &                                                        )/SLU
     &                   +(
     &                      ( (e11(M  ,N  )-e22(M  ,N  ))**2*LU(M  ,N  )
     &                       +(e11(M+1,N  )-e22(M+1,N  ))**2*LU(M+1,N  )
     &                       +(e11(M  ,N+1)-e22(M  ,N+1))**2*LU(M  ,N+1)
     &                       +(e11(M+1,N+1)-e22(M+1,N+1))**2*LU(M+1,N+1)
     &                                                        )/SLU
     &                   +4.0*e12(M  ,N  )**2      )/extr2
              
              delta_XY = SQRT(delta_XY)
              Pmean_XY = (Pice(M  ,N  )*LU(M  ,N  )
     &                   +Pice(M+1,N  )*LU(M+1,N  )
     &                   +Pice(M  ,N+1)*LU(M  ,N+1)
     &                   +Pice(M+1,N+1)*LU(M+1,N+1))/SLU
c
c Regularization: new Pice, delta.--------------------------------
              if(nreg .EQ. 1) then
                 Pmean_XY = MIN (Pmean_XY,C_Damp_XY*delta_XY)       ! El. Hunke
                 delta_XY = MAX (dmin,delta_XY)              
C               else
C                Pmean = Pmean*delta/(delta+dmin)       ! M. Harder.
C               delta=delta+dmin                       
              end if
c Stress tensor components. sig_1 = sig_D, sig_2 = sig_T.
c Implicit dte - time steps for sigma_ij.

       s12(M,N)    = (s12(M,N) + 0.5*dte*Pmean_XY
     &                       * e12(M,N) /delta_XY/Tdamp )/det2

C           END IF
         END IF

        END DO
	END DO
!$OMP END PARALLEL DO 


	IF(MMD.NE.0) THEN
          CALL CYCLIZE(SIGMA1 ,NX,NY, 1,  MMM,MM)
          CALL CYCLIZE(SIGMA2 ,NX,NY, 1,  MMM,MM)
	    CALL CYCLIZE(S12    ,NX,NY, 1,  MMM,MM)
      END IF
c Rheol. Force Fu, Fv = D sig_ij / D x_j.  -----------------------
c

!$OMP PARALLEL DO PRIVATE(FU,FV)      

C CALCULATING RHEOLOGY FORCES AND NEW VELOCITIES      
      DO N=NNN,NN
        DO M=MMM,MM
          IF(LCU(M,N).GT.0.5) THEN   
C            IF((MISTOT(M,N)+MISTOT(M+1,N)).GT.20) THEN
c
       FU= (  (SIGMA1(M+1,N)      -SIGMA1(M,N))/DXT(M,N)
     &       +(SIGMA2(M+1,N)*DY(M+1,N)**2
     &        -SIGMA2(M  ,N)*DY(M  ,N)**2)/(DXT(M,N)*DYH(M,N)**2)
     &   +2.0*(S12(M,N  )*DXB(M,N  )**2  
     &     -   S12(M,N-1)*DXB(M,N-1)**2 ) /(DXT(M,N)**2*DYH(M,N))
     &                                    )/(2.0*RN)

	 UICE(M,N)= UICE(M,N)+ 2.0/(MISTOT(M,N)+MISTOT(M+1,N)+20.0)*Fu*dte
	  
C            END IF
          END IF
        
          IF(LCV(M,N).GT.0.5) THEN
C            IF((MISTOT(M,N)+MISTOT(M,N+1)).GT.20) THEN
       
       
       FV= (  (SIGMA1(M,N+1)      -SIGMA1(M,N))/DYT(M,N)
     &       -(SIGMA2(M,N+1)*DX(M,N+1)**2
     &        -SIGMA2(M,N  )*DX(M,N  )**2)/(DYT(M,N)*DXH(M,N)**2)
     &   +2.0*(S12(M  ,N)*DYB(M  ,N)**2
     &       - S12(M-1,N)*DYB(M-1,N)**2 ) /(DYT(M,N)**2*DXH(M,N))
     &                                    )/(2.0*RN)

	 VICE(M,N)=VICE(M,N)+ 2.0/(MISTOT(M,N)+MISTOT(M,N+1)+20.0)*Fv*dte
C            END IF
          END IF
        END DO
      END DO
!$OMP END PARALLEL DO
c znak inside Fu, Fv is correct. F = ~Lapl (u; v) .
c NJ: Delta and divergence in nodes for ice ridging calculations.
c if(itime.EQ.nstep) delta_ice(i,j)=d1, div_ice(i,j)= e11+e22 .

C--------NORTH POLE FILTRATION------------------------------------------
	   CALL NPS_FILTER_T_GRD_NW(UICE,1 ) !3D ARRAY FILTERED & ITS DIM ON Z 
	   CALL NPS_FILTER_T_GRD_NW(VICE,1 )           
C----------------------------------------------------------------------   

	IF(MMD.NE.0) THEN
          CALL CYCLIZE(UICE,NX,NY, 1,  MMM,MM)
          CALL CYCLIZE(VICE,NX,NY, 1,  MMM,MM)
      END IF
      RETURN
      END

C==============================================================
C  TRANSPORT UPWIND MPDATA-SCHEME USING DIFFUSIVITY REDUCTION
      SUBROUTINE TRAN_ICE_MPDATA(FF,MGRAD,TAU,UU,VV,NITER)
	IMPLICIT NONE
      INCLUDE '0COM.INC'

      INTEGER MGRAD,NITER       !NITER - NUMBER OF ITERATIONS FOR REDUCIND DIFFUSIVITY            
      REAL  FF(NX,NY,MGRAD),TAU
	REAL  UU(NX,NY),VV(NX,NY) ! TRANSPORTING VELOCITIES


      REAL FX_P,FX_M,FY_P,FY_M   !FLUXES THROUGH CELL EDGES
    

      INTEGER M,N,K 
           
      REAL FO(NX,NY,MGRAD)               !OLD VALUE OF FF
      
      FO=FF

C   THE MAIN TRANSPORT PROCEDURE

!$OMP PARALLEL DO PRIVATE(M,N,K,FX_P,FX_M,FY_P,FY_M)
	DO N=NNN,NN
         DO M=MMM,MM
          IF(LU(M,N).GT.0.5) THEN
            
            DO K=1,MGRAD

	       FX_P=(UU(M  ,N  )+ABS(UU(M  ,N  )))
     &           *DYH(M  ,N  )    *FO(M  ,N  ,K)
     &        +   (UU(M  ,N  )-ABS(UU(M  ,N  )))
     &           *DYH(M  ,N  )    *FO(M+1,N  ,K)

	       FX_M=(UU(M-1,N  )+ABS(UU(M-1,N  )))
     &           *DYH(M-1,N  )    *FO(M-1,N  ,K)
     &        +   (UU(M-1,N  )-ABS(UU(M-1,N  )))
     &           *DYH(M-1,N  )    *FO(M  ,N  ,K)

	       FY_P=(VV(M  ,N  )+ABS(VV(M  ,N  )))
     &           *DXH(M  ,N  )    *FO(M  ,N  ,K)
     &        +   (VV(M  ,N  )-ABS(VV(M  ,N  )))
     &           *DXH(M  ,N  )    *FO(M  ,N+1,K)

	       FY_M=(VV(M  ,N-1)+ABS(VV(M  ,N-1)))
     &           *DXH(M  ,N-1)    *FO(M  ,N-1,K)
     &        +   (VV(M  ,N-1)-ABS(VV(M  ,N-1)))
     &           *DXH(M  ,N-1)    *FO(M  ,N  ,K)


            FF(M,N,K)=FO(M,N,K)-(FX_P-FX_M+FY_P-FY_M)*TAU
     &                              /(DX(M,N)*DY(M,N)*RN*2.0)
           
            END DO

           END IF
         END DO
	END DO
!$OMP END PARALLEL DO

      RETURN
	END

