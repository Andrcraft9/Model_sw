C======================================================================
      REAL FUNCTION AAAVEH(F_F,NZZ)
	IMPLICIT NONE
C----------------------------------------------------------------------
C GIVES AVERAGE OF F_F ON H-GRID WITH SPHERICITY AND VERTICAL 
C                                                    LEVEL DISTRIBUTION 
C F_F - INPUT ARRAY
C LUH - OCEAN T-GRID MASK
C RQ1 = COS(LALITUDE)
C DZ  - VERTICAL LAYER THIKNESS
C---------------------------------------------------------------------- 
      INCLUDE '0COM.INC'
      INTEGER NZZ, M, N, K
      REAL    F_F(NX,NY,NZZ)
      REAL    FAV, SWEIGHT,  DZK, WEIGHT
      
      FAV     =0.0
      SWEIGHT =0.0
      
      DO K=1,NZZ
       DZK=DZ(K)
       
       DO N=NNN-1,NN
        DO M=MMM-1,MM
	
           WEIGHT = LUH(M,N) * DXB(M,N)
     &                       * DYB(M,N)
          SWEIGHT = SWEIGHT + WEIGHT
          FAV     = FAV + WEIGHT * F_F(M,N,K)
	 
	   END DO
       END DO
      END DO
      AAAVEH = FAV / SWEIGHT
      RETURN
      END      
C======================================================================
      REAL FUNCTION WAVER(F_F)
	IMPLICIT NONE
C----------------------------------------------------------------------
C GIVES AVERAGE WITH SPHERICITY. USED FOR A SURFACE FIELD.
      INCLUDE '0COM.INC'
      REAL    F_F(NX,NY)
      INTEGER M, N
      REAL    FAV, WIE
      FAV=0.0
      WIE=0.0
!$OMP PARALLEL DO REDUCTION(+:WIE,FAV) PRIVATE(M,N)      
      DO N=NNN,NN
      DO M=MMM,MM
        IF (LU(M,N).GT.0.5) THEN
          WIE = WIE + DX(M,N)*DY(M,N)
          FAV = FAV + DX(M,N)*DY(M,N)*F_F(M,N)
        ENDIF
	END DO
	END DO
!$OMP END PARALLEL DO
      WAVER = FAV / WIE
      RETURN
      END
C======================================================02-21-96 04:14pm
      REAL FUNCTION FAVER(F_F)
C----------------------------------------------------------------------
C GIVES AVERAGE WITH TOPOGRAPHY AND SPHERICITY.
C USED FOR A LAYER OF A 3D-FIELD.
      IMPLICIT NONE      
      INCLUDE '0COM.INC'
      REAL F_F(NX,NY)
	REAL FAV, WIE
      INTEGER M, N

      FAV=0.0
      WIE=0.0
!$OMP PARALLEL DO REDUCTION(+:WIE,FAV) PRIVATE(M,N)
      DO N=NNN,NN
      DO M=MMM,MM
        IF (LU(M,N).GT.0.5) THEN
          WIE = WIE  + HHQ(M,N)*DX(M,N)*DY(M,N)
          FAV = FAV  + HHQ(M,N)*DX(M,N)*DY(M,N)* F_F(M,N)
        ENDIF
	END DO
	END DO
!$OMP END PARALLEL DO
      FAVER = FAV / WIE
      RETURN
      END
C======================================================================
      REAL FUNCTION TAVER(P_P)
C----------------------------------------------------------------------
C GIVES LINEAR AVERAGE OF A 3D - FIELD.
      INCLUDE '0COM.INC'
      INTEGER K
      REAL P_P(NX,NY,NZ)
      REAL P
      REAL     FAVER
      EXTERNAL FAVER
      P = 0.0
      DO K=1,NZ
       P = P + FAVER(P_P(1,1,K)) * DZ(K)
      ENDDO
      TAVER = P
      RETURN
      END 
C======================================================02-21-96 04:15pm
      SUBROUTINE EKINET(TT,SS,UU,VV,WW,UBRTR,VBRTR,SLH,
     &                  FF0,ECL,ECL1,ETR2,EWW2,EDN)
	IMPLICIT NONE
C----------------------------------------------------------------------
C GIVES DENSITY OF KINETIC ENERGY: BAROCLINIC = RO0*<(UU**2+VV**2)/2.>,
C VERTICAL (W): RO0*<WW**2/2.>, BAROTROPIC = RO0*(<UTR**2>+<VTR**2>)/2.
C POTENTIAL = <(-DEN + DEN(4.,35))*G*HHQ*Z>.
C WITHOUT SMOTHING.
      INCLUDE '0COM.INC'
c      INCLUDE '0DENP.INC'
	INCLUDE '0FUNCDEF.INC'
C INPUT BASE ARRAYS:
      REAL TT(NX,NY,NZ),SS(NX,NY,NZ),
     &     UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ+1)
      REAL*8 UBRTR(NX,NY),VBRTR(NX,NY),SLH(NX,NY),SPHERIC_NORM2
C Auxiliary ARRAYS
      REAL ECL,ECL1,ETR2,EWW2,EDN   !ENERGIES
      REAL FFX(NZ+1),FFY(NZ+1),FF0(NX,NY)
      INTEGER M, N, K
      REAL    ETRV, ETRU, RPV,  RPU, ECLZ, WEIZFI
      REAL    WEI

C FULL AND BAROTROPIC KINETIC ENERGY DEFINITION

      RPU  = 0.0
      RPV  = 0.0
      ETRU = 0.0
      ETRV = 0.0

      ECL  = 0.0
      WEI  = 0.0
      ECL1 = 0.0
	ETR2 = 0.0

!$OMP PARALLEL DO REDUCTION(+:WEI,ECL,ECL1,ETR2)
!$OMP&PRIVATE(M,N,K,WEIZFI,ECLZ)		
      DO N=NNN,NN
      DO M=MMM,MM

       IF (LU(M,N).GT.0.5) THEN

	 WEIZFI=HHQ(M,N)*DX(M,N)*DY(M,N)
	 ECLZ= ((UU(M,N,1)+UU(M-1,N,1))**2+
     +	    (VV(M,N,1)+VV(M,N-1,1))**2 )/8.0*DZ(1)
       
C      ENEGRY AT FIRST LEVEL	 
	 ECL1=ECL1 + ECLZ * WEIZFI  

C      BARORTOPIC ENEGRY 
	 ETR2=ETR2+ ((UBRTR(M,N)+UBRTR(M-1,N))**2+
     +	         (VBRTR(M,N)+VBRTR(M,N-1))**2 )/8.0* WEIZFI  

	 DO K=2,NZ
        ECLZ = ECLZ+((UU(M,N,K)+UU(M-1,N,K))**2+
     +               (VV(M,N,K)+VV(M,N-1,K))**2)/8.0*DZ(K)
	 END DO

	 ECL = ECL + ECLZ * WEIZFI
       WEI = WEI +        WEIZFI

       ENDIF
	END DO
	END DO
!$OMP END PARALLEL DO
      ECL = ECL / WEI         !FULL (BARORTOPIC+BAROCLINIC) ENEGRY 
      ECL1= ECL1/ WEI/DZ(1)   !ENEGRY AT FIRST LEVEL	 
      ETR2= ETR2/ WEI         !BARORTOPIC ENEGRY 
      ECL = ECL - ETR2        !BAROCLINIC ENEGRY(MAST BE > 0) 

C W-KINETIC ENERGY, ACCURATE PROCEDURE WITHOUT USING ZERO BOUNDARY
C POINTS. VERTICAL RMS WW-VELOCITY WIL BE IN FF0.

!$OMP PARALLEL DO PRIVATE(M,N,K,FFX,FFY)
      DO N=NNN,NN
       DO M=MMM,MM
        IF (LU(M,N).GT.0.5) THEN
C STARTING FROM Z(1), WHICH IS DEPTH OF WW(.,2),  UNTILL Z(NZ-1) AND
C WW(M,N,NZ). NEXT POINT CORRESPOND TO ZERO WW AND BOTTOM.
         DO K = 1,NZ-1
          FFX(K) = Z(K)
          FFY(K) = WW(M,N,K+1)**2 / 2.
         ENDDO
C LENGTH IN FFX IS NOT 1, SO
         FF0(M,N) = RINTEGR(FFX,FFY,NZ-1) / (FFX(NZ-1) - FFX(1))
        ELSE
         FF0(M,N) = 0.
        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      EWW2 = FAVER (FF0)

C------------------------------------------------------------------------------C
C POTENTIAL ENERGY = - <(DEN(TT,SS)-DEN(4.,35.))*G*Z(K)*HHQ(M,N))>
         
	   CALL NORM2_ON_TTGR8(SLH,SPHERIC_NORM2)
	   EDN=SNGL(SPHERIC_NORM2)/2.0*GRV*RH0

	
	RETURN

      END
C======================================================02-21-96 04:14pm
      SUBROUTINE NORM2_ON_TTGR8(FF,SPHERIC_NORM2)
	IMPLICIT NONE
C----------------------------------------------------------------------
C GIVES QUADRATIC NORM ON LU-GRID(TEMPERATURE) (WITH SPHERICITY).
      INCLUDE '0COM.INC'
      REAL*8 FF(NX,NY)    !ARRAY ON H-GRID
      REAL*8 FF2,WIE, WEIGHT, SPHERIC_NORM2
      INTEGER M, N

      FF2 =0.0
      WIE =0.0

!$OMP PARALLEL DO REDUCTION(+:WIE,FF2) PRIVATE(M,N,WEIGHT)
	DO N=NNN,NN
        DO M=MMM,MM
        IF (LU(M,N).GT.0.5) THEN
C	    WEIGHT = RQ1(N)*HX*HY
	    WEIGHT = DBLE(DX(M,N)*DY(M,N))
          WIE = WIE  +  WEIGHT
          FF2 = FF2  +  FF(M,N)**2*WEIGHT
        ENDIF
	END DO
	END DO
!$OMP END PARALLEL DO
      SPHERIC_NORM2 = DSQRT(FF2 / WIE)

      RETURN
      END
C======================================================================
      SUBROUTINE FZONAL(FF,FFZON)
C---------------------------------------------------------------------
C GIVES ZONAL AVERAGE FFZON OF SURFACE FIELD FF - LIKE TEMP.
      INCLUDE '0COM.INC'
      REAL FF(NX,NY),FFZON(NY)
      INTEGER M, N
      REAL RZ, FZ

      DO N=NNN,NN
       RZ = 0.
       FZ = 0.
       DO M=MMM,MM
        IF (LU(M,N).GT.0.5) THEN
         RZ = RZ + DX(M,N)*DY(M,N)
         FZ = FZ + FF(M,N)*DX(M,N)*DY(M,N)
        ENDIF
       ENDDO
       FFZON(N) = FZ / RZ
      ENDDO
      RETURN
      END

C======================================================================
      SUBROUTINE STREAM_FUNCTION_CALC(UB,SF)
	IMPLICIT NONE
      INCLUDE '0COM.INC'
      REAL(8) UB(NX,NY)
	REAL SF(NX,NY)
      INTEGER M,N
      
	SF=0.0
!$OMP PARALLEL DO PRIVATE(M,N)
	DO M=MMM,MM
	   DO N=NNN,NN
           SF(M,N)=SF(M,N-1)
     &            +SNGL(UB(M,N))*HHU(M,N)
     &                          * DYH(M,N)*RN* LCU(M,N)*1E-12
         END DO
	END DO
!$OMP END PARALLEL DO      
	IF(MMD.NE.0) CALL CYCLIZE (SF,NX,NY,1,MMM,MM)

	RETURN
	END
C======================================================31.10.99 17:23
C OCEAN GENERAL CIRCULATION MODEL. ZONAL STREAM FUNCTION.
C (C) A.V. BAGNO@IMN.RAS.RU.
      SUBROUTINE ZONFI(VV,XPT4,VVZON4,ZFI4,NZY4,EZFI4,ZLEV,NLEV,OVR)
	IMPLICIT NONE
C---------------------------------------------------------------------
C GIVES: VVZON - ZONAL AVERAGE OF ARRAY VV; XPT - NO. OF ZONAL POINTS;
C ZFI - MERIOD.STREAM FUNCTION; EZFI - DENSITY OF VVZON KINETIC ENERGY;
C NZY - NZY(NY) <= NLEV, LAST POINT OF VVZON(N,K).
C DZ1 - HEIGHT OF VV - BOXES (IN CM).
C LU - T-MASK (LIKE LU).
      INCLUDE '0COM.INC'
      INTEGER NLEV
      REAL VV(NX,NY,NLEV),
     &     ZLEV(NLEV),OVR
      REAL XPT4(NBASINS,NY,NLEV),VVZON4(NBASINS,NY,NLEV),
     & ZFI4(NBASINS,NY,NLEV+1),EZFI4(NBASINS)
C HELP ARRAYS:
      INTEGER   NHELP
      PARAMETER(NHELP=150)
      INTEGER NZY(NY), NZY4(NBASINS,NY)
      REAL EZFI, DZ1(NHELP),  ZFI(NY,NHELP+1),
     &        XPT(NY,NHELP),VVZON(NY,NHELP)
	INTEGER M, N, K, LBAS
      REAL    S1, S2, ZZ
      INTEGER NZY0
      REAL    VZ, XM, FF1, FF2

c       IF (NLEV.GT.NHELP) THEN
c       ENDIF
       DO K = 2,NLEV-1
        DZ1(K)=(ZLEV(K+1) - ZLEV(K-1))/2. * 100.
       ENDDO
       DZ1(1)= (ZLEV(2) - ZLEV(1))/2.    * 100.
       DZ1(NLEV)=(ZLEV(NLEV) - ZLEV(NLEV-1)) * 100.

      DO LBAS = 1,NBASINS    ! BASIN LOOP

      DO N=NNN,NN-1
       NZY (N) = 0
       ZFI (N,NLEV+1) = 0.
C OVR REQUIRES CAREFULL USE OF LATERAL BOUNDS OF ZFI.
       DO K=1,NLEV
C XPT AND VV ZONAL INTEGRAL.
        XPT  (N,K) = 0.
        VVZON(N,K) = OVR
         ZFI (N,K) = OVR
        XM = 0.
        VZ = 0.
        DO M=MMM,MM
         IF (LU(M,N).GT.0.5.AND.
     &           (LBAS.EQ.1.OR.LBAS.EQ.LBASINS(M,N)) ) THEN
c          IF (VV(M,N,K).GT.OVR) THEN
           IF (VV(M,N,K).GT.OVR/2.0) THEN      !Gusev's variant
           XM = XM + 1.  
           VZ = VZ + VV(M,N,K) 
          ENDIF
         ENDIF
        ENDDO
        IF (XM.GT.0.5) THEN
         XPT  (N,K) = XM
         VVZON(N,K) = VZ
        ENDIF
       ENDDO
      ENDDO
C MAKING NZY(NY) <= NLEV, LAST POINT OF VVZON(N,K).
C NZY = 0: NO POINT; NZY = 1: ONLY 1 POINT, VVZON=0; NZY =2,_: O.K.
      DO N=NNN,NN-1
       IF (VVZON(N,1).GT.OVR/2.0) THEN
        DO K=1,NLEV
         NZY0 = K
c        IF (K.EQ.NLEV.OR.VVZON(N,K+1).EQ.OVR) EXIT
         IF (K.EQ.NLEV.OR.VVZON(N,K+1).LE.OVR/2.0) EXIT !Gusev's variant
        ENDDO
        NZY(N) = NZY0
       ENDIF
      ENDDO

C MAKING VVZON BAROCLINIC - IT MAY BE NECESSARY.
      DO N=NNN,NN-1
       IF (NZY(N).GE.1) THEN
        S1 = 0.
        S2 = 0.
        DO K=1,NZY(N)
         S1 = S1 + VVZON(N,K) * DZ1(K)
         S2 = S2 +              DZ1(K)
        ENDDO
        DO K=1,NZY(N)
         VVZON(N,K) = VVZON(N,K) - S1/S2
        ENDDO
       ENDIF
      ENDDO
C 2 WAYS FOR MERIOD.STRFUN VERT.INTEGR.CAN BE USED. THIS IS FROM
C BOTTOM - LIKE IN WWINT. AT BOTTOM ZFI IS SET 0.
C SIGN '-VV' MEANS POSITIVE TRANSPORT WHEN IT IS CLOCKWISE.
      DO N=NNN,NN-1
       ZFI(N,NZY(N)+1) = 0.
       IF (NZY(N).GE.1) THEN
        ZZ = 0.
        DO K=NZY(N),1,-1
         ZZ = ZZ - VVZON(N,K) * DX(M,N)*RN *DZ1(K)
         ZFI(N,K) = ZZ
        ENDDO
       ENDIF
      ENDDO
C COMPUTING DENSITY OF VVZON KINETIC ENERGY. BOUNDARY ZERO POINTS ON
C NORTH AND SOUTH ARE NOT TAKEN.
C MAKING VV ZONAL AVERAGE.
      FF1 = 0.
      FF2 = 0.
      DO N=NNN,NN-1
       IF (NZY(N).GE.1) THEN
        DO K=1,NZY(N)
         VVZON(N,K) = VVZON(N,K) / XPT(N,K)
         FF1 = FF1 + VVZON(N,K)**2 /2. * DZ1(K)
         FF2 = FF2 +                     DZ1(K)
        ENDDO
       ENDIF
      ENDDO
      EZFI = FF1 / FF2
C MAKING  FF1, FF2:  ___  / RM(N) IS NOT NECESSARY.
C SAVING RESULTS FOR ALL BASINS.
       EZFI4 (LBAS) = EZFI
       DO N=NNN,NN-1
        NZY4 (LBAS,N) = NZY (N)
        ZFI4 (LBAS,N,NLEV+1) = 0.
        DO K=1,NLEV
         XPT4  (LBAS,N,K) =  XPT   (N,K)
         VVZON4(LBAS,N,K) =  VVZON (N,K)
         ZFI4  (LBAS,N,K) =  ZFI   (N,K)
        ENDDO
       ENDDO

      ENDDO
      RETURN
      END
C======================================================14.05.99 21:19
      SUBROUTINE MTRANSB(FF,SLRY,AFNT, AFNV,VV,VBRT,AMY,MHT,
     &	AS,AZ,AR,CONVCON)
      IMPLICIT NONE
C @N.A.DIANSKY 14.05.99 19:55 (dinar@inm.ras.ru)
C GIVES MERIDIONAL ADVECTIVE TRANSPORT OF TRACER FF AND
C MERIDIONAL COMPLEX DIFFUSIVE FLUXES OF TRACER FF.
C PURE DIFFUSIVE FLUX ALONG SIGMA LEVELS - WITH LIQUID BOUNDARY FBS, FBN.
C TRANSPORT ON LIQUID WALLS IS SUPPOSED TO BE ABSENT.
C FF --  ADVECTIVE TRACER
C SLRY - HH*Dy(DENSITY)/(dDENSITY/dZ) FOR ISOPICNAL DIFFUSIVITY
C VV --  MERIDIONAL VELOCITY IN Y - DIRECTION
C VBRT - MERIDIONAL BARTROPIC VELOCITY 
C AMY -- COEFFICIENT OF DIFFUSION IN Y-DIRECTION
C FBS -- SOUTH BOUNDARY VALUES
C FBN -- NORTH BOUNDARY VALUES
C MHT -- MERIDONAL ADVECTIVE TRANSPORTS AND DIFFUSIVE FLUXES OF FF
C        FOR "NBASINS" BASINS
C     HERE THIRD INDEX DENOTES:
C        1 - TOTAL (SUM OF 2 AND 3)
C        2 - ADVECTIVE TRANSPORT (SUM OF 4,5 AND 6)
C        3 - DIFFUSIVE FLUX (SUM OF 8 AND 9)
C        4 - TRANSPORT BY ZONALLY AVERAGE  VV OF ZONALLY AVERAGE FF (YTATAV);
C        5 - AZONAL TRANSPORT VV OF AZONAL FF (YTATV);
C        6 - TRANSPORT BY STREAMFUNCTION (YTATVB);
C        7 - PURE DIFFUSION ALONG SIGMA LEVELS;
C        8 - DIFFUSION ALONG SIGMA LEVELS AS PART OF COMPLEX DIFFUSION;
C        9 - ADDITIONAL TERM FOR DIFFUSION ALONG ISOPICNALS AND
C            HORIZONTAL LEVELS;
C GRYW -- DIRIHLET OR NEIMANN BOUNDARY CONDITION DIFFUSION
C AS -- WEIGHT COEFFICIENT FOR DIFFUSION ALONG SIGMA-LEVELS
C AZ -- WEIGHT COEFFICIENT FOR DIFFUSION ALONG HORIZONTAL LEVELS
C AR -- WEIGHT COEFFICIENT FOR DIFFUSION ALONG ISOPICNALS
C CONVCON - CONSTANT TO CONVERSE TO AN OTHER UNIT
      INCLUDE '0COM.INC'
      REAL FF(NX,NY,NZ),SLRY(NX,NY,NZ),VV(NX,NY,NZ),
     &     MHT(NBASINS,NY,9),AMY(NX,NY,NZ),
     &     AS,AZ,AR,CONVCON,ALFR,ALFZ
	REAL    AFNT(NX,NY,NZ+1), AFNV(NX,NY,NZ+1)
	REAL*8 VBRT(NX,NY)


      DOUBLE PRECISION YTATAV,YTATV,YTATVB, !ADVECTION TRANSPORT
     &                                 YTD, !S-DIFFUSION TRANSPORT
     &                     FLYY,FLYZ,       !Z-DIFFUSION TRANSPORT
     &       	FFVH,FFVB,FFT,FFM             !HELP VARIABLE

      INTEGER M, N, K, LBAS
      REAL    DC1, AMUW, AFUNTP1, AFUNVP1, AFUNT, AFUNV

      DO LBAS = 1,NBASINS     !BASIN LOOP

      DO N=NNN,NN-1
       YTD	=0.0D+00
       YTATAV =0.0D+00
       YTATV	=0.0D+00
       YTATVB	=0.0D+00
       FLYY	=0.0D+00
       FLYZ	=0.0D+00

       DO K=1,NZ

C OLD:INNER INTERVAL IS TAKEN.
C OLD:THERE WAS NO IF HERE BECAUSE IT WAS SUPPOSED THAT ALL LINES IN Y
C BELONG TO REGION. WHEN LRX >= 2 PROBLEMS COULD BE HERE AND IN MERFI.
C TRANSPORT ON X-BOUNDARY IS SUPPOSED TO BE 0 DUE TO ZERO VV, FI AND
C DIFFUSION - DUE TO INSULATION.
        FFVH = 0.0D+00
        FFVB = 0.0D+00
        FFT  = 0.0D+00
        FFM  = 0.0D+00
C ZONAL AVERAGE FOR YTATAV, ZONAL MHT.
        DO M=MMM,MM
         IF(LCV(M,N).GT.0.5.AND.
     &     (LBAS.EQ.1.OR.LBASINS(M,N).EQ.LBAS)) THEN
          FFVH = FFVH +   VV(M,N,K)*HHV(M,N)
          FFVB = FFVB + VBRT(M,N  )*HHV(M,N)
          FFT  = FFT  + (FF(M,N+1,K)+FF(M,N,K))/2.0
          FFM  = FFM  + 1.0
         ENDIF
        ENDDO

        IF(FFM.GT.0.5) THEN
        FFVH = FFVH / FFM
        FFVB = FFVB / FFM
        FFT  = FFT  / FFM

C Y-TRANSPORT AS <V><T>
        YTATAV  = YTATAV + FFVH * FFT * FFM *DXH(M,N)*RN* DZ(K)
        END IF

C AZONAL B-FULL, B-TROPIC & DIFFUSIVE MERID.TRANSPORTS.
C IN YTATVB FFT IS SUBTRACTED FOR GREATER ACCURACY.
        
	  DO M=MMM,MM

         IF(LCV(M,N).GT.0.5.AND.
     &     (LBAS.EQ.1.OR.LBASINS(M,N).EQ.LBAS)) THEN
C Y-TRANSPORT AS <V'T'>
          YTATV = YTATV +
     &    (VV(M,N,K)*HHV(M,N) - FFVH) *
     &          ((FF(M,N+1,K)+FF(M,N,K))/2. - FFT ) *DXH(M,N)*RN*DZ(K)

C Y-TRANSPORT AS <VB'T'>
          YTATVB = YTATVB +
     &    (VBRT(M,N)*HHV(M,N) - FFVB) *
     &          ((FF(M,N+1,K)+FF(M,N,K))/2. - FFT ) *DXH(M,N)*RN*DZ(K)

C Y-TRANSPORT BY DIFFUSION ON SIMA LEVELS
          YTD  = YTD - AMY(M,N,K)*HHV(M,N)*
     &           (FF(M,N+1,K)-FF(M,N,K))*DXH(M,N)/DYT(M,N)*DZ(K)
         ENDIF
        ENDDO
       ENDDO

C COMPLEX DIFFUSIVE MERIDIONAL FLUXES
       DO K=2,NZ

            IF(K.EQ.2)  THEN
               DC1 = 1.5
            ELSEIF(K.LT.NZ) THEN
               DC1 = 1.0
            ELSE
               DC1 = 1.5
            END IF

        DO M=MMM,MM
         IF(LCV(M,N).GT.0.5.AND.
     &     (LBAS.EQ.1.OR.LBASINS(M,N).EQ.LBAS)) THEN
C  DIFFUSION ALONG SIGMA LEVELS AS PART OF COMPLEX DIFFUSION:
C  WEIGHT COEFFICIENTS FOR Z-AND-ISOPICNAL DIFFUSION
             
             AFUNT = AFNT(M,N,K)
             AFUNV = AFNV(M,N,K)
             
             AFUNTP1 = AFNT(M,N+1,K)
             AFUNVP1 = AFNV(M,N+1,K)
             
             ALFR=0.5*(AFUNT+AFUNTP1)
             ALFZ=(1.0-ALFR)*AR+AZ
             ALFR=ALFR*AR
C  COEF. OF DIFFUSION IN W-POINTS AS FUNCTION ON DEPTH:
          AMUW=DC1*0.25*(AMY(M,N,K-1)+AMY(M,N,K))*DXH(M,N)/DYT(M,N)
     &                      *(AFUNV+AFUNVP1)
          FLYY = FLYY -
     -                AS*AMUW*HZT(K)*HHV(M,N)*
     *               ((FF(M,N+1,K  )-FF(M,N,K  ))+
     +                (FF(M,N+1,K-1)-FF(M,N,K-1)) )/2.

C  ADDITIONAL TERM FOR DIFFUSION ALONG ISOPICNALS AND HORIZONTAL LEVELS:
          FLYZ = FLYZ +
     +              AMUW*(ALFR*SLRY(M,N,K)+
     +                    ALFZ*ZW(K)*(HHQ(M,N+1)-HHQ(M,N)))*
     *            ((FF(M,N+1,K)-FF(M,N+1,K-1))+
     +             (FF(M,N  ,K)-FF(M,N  ,K-1)) )/2.0
         ENDIF
        ENDDO
       ENDDO

       MHT(LBAS,N,1)=(YTATAV+YTATV+FLYY+FLYZ)*CONVCON
       MHT(LBAS,N,2)=(YTATAV+YTATV)*CONVCON
       MHT(LBAS,N,3)=(FLYY+FLYZ)*CONVCON
       MHT(LBAS,N,4)= YTATAV *CONVCON
       MHT(LBAS,N,5)= YTATV*CONVCON
       MHT(LBAS,N,6)= YTATVB*CONVCON
       MHT(LBAS,N,7)= YTD *CONVCON
       MHT(LBAS,N,8)= FLYY*CONVCON
       MHT(LBAS,N,9)= FLYZ*CONVCON

      ENDDO
C LIQUID BOUNDARY POINTS : TRANSPORT IS ZERO, DIFFUSION IS ZERO ONLY
C WHEN IGRYW = 2 , IF IGRYW= 1 FBN, FBS ARE USED.
C N=NNN-1 & N=NN - BOUNDARY FOR TRANSPORT AND IT IS SUPPOSED TO BE 0.
C      DO K = 1,9
C      MHT(LBAS,NNN-1,K) = 0.0E+00
C      MHT(LBAS,NN   ,K) = 0.0E+00
C      ENDDO

      ENDDO       !END OF BASIN LOOP
      RETURN
      END
C======================================================================
      SUBROUTINE MHTBAL(FLUX,BMHT,CONVCF)
C @N.A.DIANSKY 31.10.99 17:23 (dinar@inm.ras.ru)
C GIVES MERIDIONAL ADVECTIVE TRANSPORT AS BALANCE WITH SEA SURFACE FLUX
C FLUX - SEA SURFACE FLUX
C BMHT(L,N,1) - ZONAL MEAN FLUX
C BMHT(L,N,2) - MERIDONAL ADVECTIVE TRANSPORTS AS BALANCE WITH
C               SEA SURFACE FLUX  (ALL 1&2 GIVEN ON H(U)-GRID)
C               WHERE L - NUMBER OF BASIN, N- NUMBER OF LATITUDE
C CONVCF - COEFFICIENT OF CONVERSATION
      INCLUDE '0COM.INC'
      REAL    FLUX(NX,NY),BMHT(NBASINS,NY,2),CONVCF
	INTEGER M, N, LBAS
      REAL    FFM, TRMEAN, TRINT

      DO LBAS = 1,NBASINS
      DO N = 1,NY
       BMHT(LBAS,N,1)=0.0
       BMHT(LBAS,N,2)=0.0
      ENDDO

      TRINT=0.0

      DO N=NNN,NN-1
        TRMEAN=0.0
        FFM=0.0
        DO M=MMM,MM
         IF(LCV(M,N).GT.0.5.AND.
     &     (LBAS          .EQ.1   .OR.
     &      LBASINS(M,N)  .EQ.LBAS.OR.
     &      LBASINS(M,N+1).EQ.LBAS)) THEN
          TRMEAN = TRMEAN+(FLUX(M,N+1)+FLUX(M,N))*DXH(M,N)*DYT(M,N)/2.0
          FFM=FFM+DXH(M,N)*DYT(M,N)
         ENDIF
        ENDDO

        IF(FFM.GT.0.5) THEN
        BMHT(LBAS,N,1)= TRMEAN/FFM
        TRINT = TRINT + TRMEAN*RN**2*CONVCF
        BMHT(LBAS,N,2)= TRINT
        END IF

       ENDDO

      ENDDO

      RETURN
      END
C=======================================================================================
CCCC SUBROUTINE OF MERIDIONAL OVERTURNING CIRCULATION COMPUTATION IN SIGMA-COORDINATES      
      SUBROUTINE MOC_CALC(MOC,NX,NY,NZ,LU,VV,HHQ,MMM,MM,NNN,NN,DZ)
      IMPLICIT NONE
      
      INTEGER NX,NY,NZ,MMM,MM,NNN,NN
	REAL LU(NX,NY),VV(NX,NY,NZ), HHQ(NX,NY), ! TEMPERATURE MASK, MERIDIONAL VELOCITY AND TOPOGRAPHY
     &    MOC(NX,NY,NZ+1),DZ(NZ)               ! OVERTURNING STREAMFUNCTION
      
      REAL V_AVE
      INTEGER M,N,K

	MOC=0.0

!$OMP PARALLEL DO PRIVATE(M,N,K,V_AVE)
	DO N=NNN,NN
	  DO M=MMM,MM
        
         IF(LU(M,N).GT.0.5) THEN
	    
          V_AVE=0.0          
C  CALCULATING BAROTROPIC VELOCITY          
          DO K=1,NZ
          V_AVE=V_AVE+VV(M,N,K)*DZ(K)
          END DO
C  COMPUTING MOC BY USING BAROCLINIC VELOCITY
          DO K=1,NZ
           MOC(M,N,K+1)=MOC(M,N,K)+(VV(M,N,K)-V_AVE)*DZ(K)*HHQ(M,N)
          END DO	  
	   
         END IF
        
        END DO
      END DO
!$OMP END PARALLEL DO

	RETURN
	END
C=======================================================================================
CCCC SUBROUTINE OF MERIDIONAL OVERTURNING INTEGRATION IN Z-COORDINATES
      SUBROUTINE MOC_INTEGRATION(MOC,NX,NY,NZ,NBASINS,MMM,MM,NNN,NN,
     &                           MOC_INT,LU,LBASINS,DX,FACTOR,OVER)
      IMPLICIT NONE
	INTEGER NX,NY,NZ,NBASINS,LBAS,MMM,MM,NNN,NN,M,N,K
	REAL MOC(NX,NY,NZ),MOC_INT(NBASINS,NY,NZ),
     &      LU(NX,NY),DX(NX,NY),FACTOR,OVER
      INTEGER LBASINS(NX,NY)
	INTEGER,ALLOCATABLE::  NPT(:,:,:)
      
      ALLOCATE (NPT(NBASINS,NY,NZ))
      
	NPT=0
      MOC_INT=0.0

C  INTEGRATING OVERTURNING SF OVER LON-LINES
      DO LBAS=1,NBASINS
!$OMP PARALLEL DO PRIVATE(M,N,K)	 
       DO K=1,NZ       
        DO N=NNN,NN
	   DO M=MMM,MM
          
          IF(LU(M,N).GT.0.5.AND.MOC(M,N,K).GT.OVER/2.0.AND.
     &     (LBAS.EQ.1.OR.LBASINS(M,N).EQ.LBAS)) THEN

             MOC_INT(LBAS,N,K)=MOC_INT(LBAS,N,K)
     &                        +MOC(M,N,K)*DX(M,N)*FACTOR
	       NPT(LBAS,N,K)=NPT(LBAS,N,K)+1
          END IF
	   
         END DO
        END DO
       END DO
!$OMP END PARALLEL DO
	END DO
      
C  FILLINF MISSING VALIES BY UNDEF VALUES      
      DO LBAS=1,NBASINS
!$OMP PARALLEL DO PRIVATE(N,K)	 
       DO K=1,NZ       
        DO N=NNN,NN
          IF(NPT(LBAS,N,K).EQ.0) THEN
             MOC_INT(LBAS,N,K)=OVER            
          END IF
        END DO
       END DO
!$OMP END PARALLEL DO
      END DO


      DEALLOCATE(NPT)
	RETURN
	END
C======================================================
      SUBROUTINE TS_FLUXES_4_MTRANS(FF,ADV_FL_X,ADV_FL_Y,
     &                               DIFF_FL_X,DIFF_FL_Y)
	IMPLICIT NONE
      INCLUDE '0COM.INC'
	INCLUDE '0CEAN.INC'
      REAL FF(NX,NY,NZ)
	REAL ADV_FL_X(NX,NY), ADV_FL_Y(NX,NY),
     &    DIFF_FL_X(NX,NY),DIFF_FL_Y(NX,NY)
      INTEGER M,N,K
      REAL AS,AZ,AR,ALFR,ALFZ
      REAL    DC1, AMUW, AFUNTP1, AFUNVP1, AFUNT, AFUNV
      
      ADV_FL_X=0.0
      ADV_FL_Y=0.0
      
      DIFF_FL_X=0.0
      DIFF_FL_Y=0.0

      AR=WCLRD
	AS=WCLSD
	AZ=WCLZD

!$OMP PARALLEL DO PRIVATE(M,N,K)      
      DO N=NNN,NN
       DO M=MMM,MM
	  
        IF(LCU(M,N).GT.0.5) THEN
         DO K=1,NZ
	    ADV_FL_X(M,N)=ADV_FL_X(M,N)
     &     +UU(M,N,K)*(FF(M,N,K)+FF(M+1,N,K))/2.0
     &     *HHU(M,N)*DZ(K)
         END DO
	  END IF
        
        IF(LCV(M,N).GT.0.5) THEN
         DO K=1,NZ
	    ADV_FL_Y(M,N)=ADV_FL_Y(M,N)
     &     +VV(M,N,K)*(FF(M,N,K)+FF(M,N+1,K))/2.0
     &     *HHV(M,N)*DZ(K)
         END DO
	  END IF
       
       END DO
	END DO
!OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(M,N,K,DC1,AFUNT,AFUNV,AFUNTP1,AFUNVP1,
!$OMP&         ALFR,ALFZ,AMUW)
      DO N=NNN,NN
       DO K=2,NZ        
            
          IF(K.EQ.2)  THEN
            DC1 = 1.5
          ELSEIF(K.LT.NZ) THEN
            DC1 = 1.0
          ELSE
            DC1 = 1.5
          END IF
        
        DO M=MMM,MM
         IF(LCU(M,N).GT.0.5) THEN
             
             AFUNT = AFNT(M,N,K)
             AFUNV = AFNV(M,N,K)
             
             AFUNTP1 = AFNT(M+1,N,K)
             AFUNVP1 = AFNV(M+1,N,K)
             
             ALFR=0.5*(AFUNT+AFUNTP1)
             ALFZ=(1.0-ALFR)*AR+AZ
             ALFR=ALFR*AR
C  COEF. OF DIFFUSION IN W-POINTS AS FUNCTION ON DEPTH:
          AMUW=DC1*0.25*(AMXT(M,N,K-1)+AMXT(M,N,K))/(DXT(M,N)*RN)
     &                      *(AFUNV+AFUNVP1)
          
          DIFF_FL_X(M,N) = DIFF_FL_X(M,N) -
     -                AS*AMUW*HZT(K)*HHU(M,N)*
     *               ((FF(M+1,N,K  )-FF(M,N,K  ))+
     +                (FF(M+1,N,K-1)-FF(M,N,K-1)) )/2.   +
     
     +              AMUW*(ALFR*SLRX(M,N,K)+
     +                    ALFZ*ZW(K)*(HHQ(M+1,N)-HHQ(M,N)))*
     *            ((FF(M+1,N,K)-FF(M+1,N,K-1))+
     +             (FF(M,N  ,K)-FF(M,N  ,K-1)) )/2.0       
         
         END IF
        END DO
	 
       END DO
	END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(M,N,K,DC1,AFUNT,AFUNV,AFUNTP1,AFUNVP1,
!$OMP&         ALFR,ALFZ,AMUW)
      DO M=MMM,MM
       DO K=2,NZ        
            
          IF(K.EQ.2)  THEN
            DC1 = 1.5
          ELSEIF(K.LT.NZ) THEN
            DC1 = 1.0
          ELSE
            DC1 = 1.5
          END IF
        
        DO N=NNN,NN
         IF(LCV(M,N).GT.0.5) THEN
             
             AFUNT = AFNT(M,N,K)
             AFUNV = AFNV(M,N,K)
             
             AFUNTP1 = AFNT(M,N+1,K)
             AFUNVP1 = AFNV(M,N+1,K)
             
             ALFR=0.5*(AFUNT+AFUNTP1)
             ALFZ=(1.0-ALFR)*AR+AZ
             ALFR=ALFR*AR
C  COEF. OF DIFFUSION IN W-POINTS AS FUNCTION ON DEPTH:
          AMUW=DC1*0.25*(AMYT(M,N,K-1)+AMYT(M,N,K))/(DYT(M,N)*RN)
     &                      *(AFUNV+AFUNVP1)
          
          DIFF_FL_Y(M,N) = DIFF_FL_Y(M,N) -
     -                AS*AMUW*HZT(K)*HHV(M,N)*
     *               ((FF(M,N+1,K  )-FF(M,N,K  ))+
     +                (FF(M,N+1,K-1)-FF(M,N,K-1)) )/2.   +
     
     +              AMUW*(ALFR*SLRY(M,N,K)+
     +                    ALFZ*ZW(K)*(HHQ(M,N+1)-HHQ(M,N)))*
     *            ((FF(M,N+1,K)-FF(M,N+1,K-1))+
     +             (FF(M,N  ,K)-FF(M,N  ,K-1)) )/2.0       
         
         END IF
        END DO
	 
       END DO
	END DO
!$OMP END PARALLEL DO

      IF(MMD.NE.0) THEN
	 CALL CYCLIZE ( ADV_FL_X,NX,NY,1,MMM,MM)
	 CALL CYCLIZE (DIFF_FL_X,NX,NY,1,MMM,MM)
      END IF

      RETURN
	END
C========================================================
      SUBROUTINE MTRANS_ZONAVE(ADV_Y,DIFF_Y,NX,NY,MMM,MM,NNN,NN,
     &                         LU,DX,FACTOR,CONVCON,
     &                         MHT,NBASINS,LBASINS)
      IMPLICIT NONE
      INTEGER NX,NY,MMM,MM,NNN,NN,NBASINS
      REAL ADV_Y(NX,NY),DIFF_Y(NX,NY),MHT(NBASINS,NY,3),LU(NX,NY),
     &                      DX(NX,NY),FACTOR,CONVCON
	INTEGER LBASINS(NX,NY)
      INTEGER M,N,LBAS

      MHT=0.0

      DO LBAS=1,NBASINS
!$OMP PARALLEL DO PRIVATE(M,N)
       DO N=NNN,NN
        DO M=MMM,MM
         IF(LU(M,N).GT.0.5.AND.
     &     (LBAS.EQ.1.OR.LBASINS(M,N).EQ.LBAS)) THEN	 
           MHT(LBAS,N,1)=MHT(LBAS,N,1)+
     &                   (ADV_Y(M,N)+DIFF_Y(M,N))*DX(M,N)*FACTOR*CONVCON
           MHT(LBAS,N,2)=MHT(LBAS,N,2)+
     &                   (ADV_Y(M,N)            )*DX(M,N)*FACTOR*CONVCON
           MHT(LBAS,N,3)=MHT(LBAS,N,3)+
     &                   (           DIFF_Y(M,N))*DX(M,N)*FACTOR*CONVCON              
         END IF
        END DO
	 END DO
!$OMP END PARALLEL DO      
      END DO

	RETURN
	END

C=================================================================
      SUBROUTINE MOC_CALC_INT_XZ(VV_Z,NX,NY,NLEV,LU,MOC_INT,
     &                NBASINS,LBASINS,MMM,MM,NNN,NN,ZLEV,DX,FACTOR,OVER)
      IMPLICIT NONE
	INTEGER NX,NY,NLEV,NBASINS,MMM,MM,NNN,NN
	INTEGER LBASINS(NX,NY)
      REAL VV_Z(NX,NY,NLEV),LU(NX,NY),DX(NX,NY),FACTOR,OVER,
     &       ZLEV(NLEV),MOC_INT(NBASINS,NY,NLEV)
      
	REAL CALC_V,CALC_M,V_ZAVE, CALC_MOC
      INTEGER NPT
      REAL,   ALLOCATABLE:: DZ1(:)
      INTEGER M,N,K,KM1,LBAS
      
      ALLOCATE(DZ1(NLEV))
      
	DO K=2,NLEV-1
	DZ1(K)=(ZLEV(K+1)-ZLEV(K-1))/2.0*100.0
	END DO
	
      DZ1(1)   =(ZLEV(1)+ZLEV(2))/2.0*100.0
	DZ1(NLEV)=(ZLEV(NLEV)-ZLEV(NLEV-1))*100.0

C  CORRECTING V-VELOCITY TO PROVIDE ZONAL AND DEPTH INTEGRAL TO BE EQUAL 0
!$OMP PARALLEL DO PRIVATE(M,N,K,CALC_M,CALC_V,V_ZAVE)	
      DO N=NNN,NN
	 CALC_M=0.0
	 CALC_V=0.0

	 DO K=1,NLEV
        DO M=MMM,MM
          IF(LU(M,N).GT.0.5.AND.ABS(VV_Z(M,N,K)).GT.0.0) THEN
             CALC_M=CALC_M+DZ1(K)*DX(M,N)
             CALC_V=CALC_V+DZ1(K)*DX(M,N)*VV_Z(M,N,K)
          END IF
        END DO
       END DO
	 
	 IF(CALC_M.GT.0.0) THEN
	 
       V_ZAVE=CALC_V/CALC_M

	  DO K=1,NLEV       
         DO M=MMM,MM
          IF(LU(M,N).GT.0.5.AND.ABS(VV_Z(M,N,K)).GT.0.0) THEN
           VV_Z(M,N,K)=VV_Z(M,N,K)-V_ZAVE
          END IF
         END DO
        END DO
       END IF

	END DO
!$OMP END PARALLEL DO

      MOC_INT=0.0

C  INTEGRATING OVERTURNING SF OVER LON-LINES
      DO LBAS=1,NBASINS
!$OMP PARALLEL DO PRIVATE(M,N,K,CALC_MOC,NPT,KM1)	 
       DO N=NNN,NN
        
        DO K=1,NLEV 
	   CALC_MOC=0.0
	   NPT=0
 
         DO M=MMM,MM
          IF(LU(M,N).GT.0.5.AND.ABS(VV_Z(M,N,K)).GT.0.0.AND.
     &     (LBAS.EQ.1.OR.LBASINS(M,N).EQ.LBAS)) THEN
             CALC_MOC=CALC_MOC
     &           +VV_Z(M,N,K)*DX(M,N)*FACTOR
	       NPT=NPT+1
          END IF
         END DO
         
         
         IF(NPT.GT.0) THEN
           KM1=MAX(K-1,1)
           MOC_INT(LBAS,N,K)=MOC_INT(LBAS,N,KM1)
     &                   +CALC_MOC*DZ1(K)
         ELSE
           MOC_INT(LBAS,N,K)=OVER
	   END IF
        END DO

       END DO
!$OMP END PARALLEL DO
	END DO
      

	DEALLOCATE(DZ1)
	RETURN
	END
C======================================================
      SUBROUTINE TS_BALANCE_4_MTRANS(FF,ADV_FLUX,DIFF_FLUX)
	IMPLICIT NONE
      INCLUDE '0COM.INC'
	INCLUDE '0CEAN.INC'
      REAL FF(NX,NY,NZ)
	REAL ADV_FLUX(NX,NY), DIFF_FLUX(NX,NY)
      INTEGER M,N,K
      REAL AS,AZ,AR,ALFR,ALFZ
      REAL    DC1, AMUW, AFUNTP1, AFUNVP1, AFUNT, AFUNV
      REAL ADVFX_P,ADVFX_M,ADVFY_P,ADVFY_M,DIFF_P,DIFF_M

      ADV_FLUX=0.0
      
      DIFF_FLUX=0.0

      AR=WCLRD
	AS=WCLSD
	AZ=WCLZD

!$OMP PARALLEL DO PRIVATE(M,N,K,ADVFX_P,ADVFX_M,ADVFY_P,ADVFY_M)      
      DO N=NNN,NN
       DO M=MMM,MM
	  
        IF(LU(M,N).GT.0.5) THEN
         DO K=1,NZ
	    ADVFX_P=LCU(M  ,N  ) * UU(M  ,N  ,K)
     &           *(FF(M  ,N  ,K)+FF(M+1,N  ,K))
     &           *HHU(M  ,N  ) *DYH(M  ,N  )
          
          ADVFX_M=LCU(M-1,N  ) * UU(M-1,N  ,K)
     &           *(FF(M  ,N  ,K)+FF(M-1,N  ,K))
     &           *HHU(M-1,N  ) *DYH(M-1,N  )
          
          ADVFY_P=LCV(M  ,N  ) * VV(M  ,N  ,K)
     &           *(FF(M  ,N  ,K)+FF(M  ,N+1,K))
     &           *HHV(M  ,N  ) *DXH(M  ,N  )
          
          ADVFY_M=LCV(M  ,N-1) * VV(M  ,N-1,K)
     &           *(FF(M  ,N  ,K)+FF(M  ,N-1,K))
     &           *HHV(M  ,N-1) *DXH(M  ,N-1)

	    ADV_FLUX(M,N)=ADV_FLUX(M,N)
     &    -(ADVFX_P-ADVFX_M+ADVFY_P-ADVFY_M)*DZ(K)
     &           /(2.0*DX(M,N)*DY(M,N)*RN)
         END DO
	  END IF
       
       END DO
	END DO
!OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(M,N,K,DC1,AFUNT,AFUNV,AFUNTP1,AFUNVP1,
!$OMP&         ALFR,ALFZ,AMUW,DIFF_P,DIFF_M)
      DO N=NNN,NN
       DO K=2,NZ        
            
          IF(K.EQ.2)  THEN
            DC1 = 1.5
          ELSEIF(K.LT.NZ) THEN
            DC1 = 1.0
          ELSE
            DC1 = 1.5
          END IF
        
        DO M=MMM,MM
         IF(LU(M,N).GT.0.5) THEN
             
             AFUNT = AFNT(M,N,K)
             AFUNV = AFNV(M,N,K)
             
             AFUNTP1 = AFNT(M+1,N,K)
             AFUNVP1 = AFNV(M+1,N,K)
             
             ALFR=0.5*(AFUNT+AFUNTP1)
             ALFZ=(1.0-ALFR)*AR+AZ
             ALFR=ALFR*AR
C  COEF. OF DIFFUSION IN W-POINTS AS FUNCTION ON DEPTH:
          AMUW=LCU(M,N)*DC1*0.25*(AMXT(M,N,K-1)+AMXT(M,N,K))
     *           /(DXT(M,N)*RN)*(AFUNV+AFUNVP1)*DYH(M,N)
          
          DIFF_P =  -  AS*AMUW*HZT(K)*HHU(M,N)*
     *               ((FF(M+1,N,K  )-FF(M,N,K  ))+
     +                (FF(M+1,N,K-1)-FF(M,N,K-1)) )/2.   +
     
     +              AMUW*(ALFR*SLRX(M,N,K)+
     +                    ALFZ*ZW(K)*(HHQ(M+1,N)-HHQ(M,N)))*
     *            ((FF(M+1,N,K)-FF(M+1,N,K-1))+
     +             (FF(M,N  ,K)-FF(M,N  ,K-1)) )/2.0       

             AFUNT = AFNT(M-1,N,K)
             AFUNV = AFNV(M-1,N,K)
             
             AFUNTP1 = AFNT(M,N,K)
             AFUNVP1 = AFNV(M,N,K)
             
             ALFR=0.5*(AFUNT+AFUNTP1)
             ALFZ=(1.0-ALFR)*AR+AZ
             ALFR=ALFR*AR
C  COEF. OF DIFFUSION IN W-POINTS AS FUNCTION ON DEPTH:
          AMUW=LCU(M-1,N)*DC1*0.25*(AMXT(M-1,N,K-1)+AMXT(M-1,N,K))
     &                      /(DXT(M-1,N)*RN)
     &                      *(AFUNV+AFUNVP1)*DYH(M-1,N)
          
          DIFF_M =  -  AS*AMUW*HZT(K)*HHU(M-1,N)*
     *               ((FF(M,N,K  )-FF(M-1,N,K  ))+
     +                (FF(M,N,K-1)-FF(M-1,N,K-1)) )/2.   +
     
     +              AMUW*(ALFR*SLRX(M-1,N,K)+
     +                    ALFZ*ZW(K)*(HHQ(M,N)-HHQ(M-1,N)))*
     *            ((FF(M  ,N,K)-FF(M  ,N,K-1))+
     +             (FF(M-1,N,K)-FF(M-1,N,K-1)) )/2.0  
         
	    DIFF_FLUX(M,N)=DIFF_FLUX(M,N)
     &    -(DIFF_P-DIFF_M)/(DX(M,N)*DY(M,N)*RN)
         END IF
        END DO
	 
       END DO
	END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(M,N,K,DC1,AFUNT,AFUNV,AFUNTP1,AFUNVP1,
!$OMP&         ALFR,ALFZ,AMUW,DIFF_P,DIFF_M)
      DO M=MMM,MM
       DO K=2,NZ        
            
          IF(K.EQ.2)  THEN
            DC1 = 1.5
          ELSEIF(K.LT.NZ) THEN
            DC1 = 1.0
          ELSE
            DC1 = 1.5
          END IF
        
        DO N=NNN,NN
         IF(LU(M,N).GT.0.5) THEN
             
             AFUNT = AFNT(M,N,K)
             AFUNV = AFNV(M,N,K)
             
             AFUNTP1 = AFNT(M,N+1,K)
             AFUNVP1 = AFNV(M,N+1,K)
             
             ALFR=0.5*(AFUNT+AFUNTP1)
             ALFZ=(1.0-ALFR)*AR+AZ
             ALFR=ALFR*AR
C  COEF. OF DIFFUSION IN W-POINTS AS FUNCTION ON DEPTH:
          AMUW=LCV(M,N)*DC1*0.25*(AMYT(M,N,K-1)+AMYT(M,N,K))
     &                      /(DYT(M,N)*RN)
     &                      *(AFUNV+AFUNVP1)*DXH(M,N)
          
          DIFF_P =  -  AS*AMUW*HZT(K)*HHV(M,N)*
     *               ((FF(M,N+1,K  )-FF(M,N,K  ))+
     +                (FF(M,N+1,K-1)-FF(M,N,K-1)) )/2.   +
     
     +              AMUW*(ALFR*SLRY(M,N,K)+
     +                    ALFZ*ZW(K)*(HHQ(M,N+1)-HHQ(M,N)))*
     *            ((FF(M,N+1,K)-FF(M,N+1,K-1))+
     +             (FF(M,N  ,K)-FF(M,N  ,K-1)) )/2.0       
             
             AFUNT = AFNT(M,N-1,K)
             AFUNV = AFNV(M,N-1,K)
             
             AFUNTP1 = AFNT(M,N,K)
             AFUNVP1 = AFNV(M,N,K)
             
             ALFR=0.5*(AFUNT+AFUNTP1)
             ALFZ=(1.0-ALFR)*AR+AZ
             ALFR=ALFR*AR
C  COEF. OF DIFFUSION IN W-POINTS AS FUNCTION ON DEPTH:
          AMUW=LCV(M,N-1)*DC1*0.25*(AMYT(M,N-1,K-1)+AMYT(M,N-1,K))
     &                      /(DYT(M,N-1)*RN)
     &                      *(AFUNV+AFUNVP1)*DXH(M,N-1)
          
          DIFF_M =  -  AS*AMUW*HZT(K)*HHV(M,N-1)*
     *               ((FF(M,N,K  )-FF(M,N-1,K  ))+
     +                (FF(M,N,K-1)-FF(M,N-1,K-1)) )/2.   +
     
     +              AMUW*(ALFR*SLRY(M,N-1,K)+
     +                    ALFZ*ZW(K)*(HHQ(M,N)-HHQ(M,N-1)))*
     *            ((FF(M,N  ,K)-FF(M,N  ,K-1))+
     +             (FF(M,N-1,K)-FF(M,N-1,K-1)) )/2.0            

	    DIFF_FLUX(M,N)=DIFF_FLUX(M,N)
     &    -(DIFF_P-DIFF_M)/(DX(M,N)*DY(M,N)*RN)

         END IF
        END DO
	 
       END DO
	END DO
!$OMP END PARALLEL DO

      IF(MMD.NE.0) THEN
	 CALL CYCLIZE ( ADV_FLUX,NX,NY,1,MMM,MM)
	 CALL CYCLIZE (DIFF_FLUX,NX,NY,1,MMM,MM)
      END IF

      RETURN
	END
C========================================================
      SUBROUTINE MTRANS_INTEGRAL(ADV_FLUX,DIFF_FLUX,
     &                         NX,NY,MMM,MM,NNN,NN,
     &                         LU,DX,DY,FACTOR,CONVCON,
     &                         MHT,NBASINS,LBASINS)
      IMPLICIT NONE
      INTEGER NX,NY,MMM,MM,NNN,NN,NBASINS
      REAL ADV_FLUX(NX,NY),DIFF_FLUX(NX,NY),
     &          MHT(NBASINS,NY,3),LU(NX,NY),
     &          DX(NX,NY),DY(NX,NY),FACTOR,CONVCON

      REAL MHT_ADV,MHT_DIF
	INTEGER LBASINS(NX,NY)
      INTEGER M,N,LBAS

      MHT=0.0

      DO LBAS=1,NBASINS
       DO N=NN,NNN,-1
        
        MHT_ADV=0.0
        MHT_DIF=0.0

        DO M=MMM,MM
         IF(LU(M,N).GT.0.5.AND.
     &     (LBAS.EQ.1.OR.LBASINS(M,N).EQ.LBAS)) THEN	 
         
         MHT_ADV= MHT_ADV+ ADV_FLUX(M,N)*DX(M,N)*FACTOR  
         MHT_DIF= MHT_DIF+DIFF_FLUX(M,N)*DX(M,N)*FACTOR
                   
         END IF
        END DO
        
        MHT(LBAS,N-1,2)=MHT(LBAS,N,2)+
     &                   MHT_ADV*DY(M,N)*FACTOR*CONVCON
        MHT(LBAS,N-1,3)=MHT(LBAS,N,3)+
     &                   MHT_DIF*DY(M,N)*FACTOR*CONVCON
        MHT(LBAS,N-1,1)=MHT(LBAS,N-1,2)+MHT(LBAS,N-1,3)

	 END DO
      END DO

	RETURN
	END
