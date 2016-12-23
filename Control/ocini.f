C===============================================================================
      SUBROUTINE  ATM_AR_2_ZERO()
      IMPLICIT NONE
      INCLUDE '2FUNDCONST.INC'
      INCLUDE '1ATMFORCING.INC'

      QBL = 0.0
      PME = 0.0
      QSW = 0.0
      SIC = 0.0
      TSA = 0.0
      TXA = 0.0
      TYA = 0.0
      SLP = 0.0
      LWAT= 0.0
      SWAT= 0.0
      PRAT= 0.0
      TAT = 0.0
      QAT = 0.0
      UAT = 0.0
      VAT = 0.0

      RETURN
      END

C======================================================================
C(C) DIANSKY N.A. (dinar@inm.ras.ru),
C PROGRAM MODULES FOR INPUT PARAMETERS FOR THE OGCM
C=======================================================15.10.03
      SUBROUTINE OCZEROSET()
      use mod_shallowater
      IMPLICIT NONE
      INTEGER MEMORY_INDEX
C---------------------------------------------------------------------
      INCLUDE '0COM.INC'
      INCLUDE '0CEAN.INC'
      INCLUDE '1LEVS.INC'
	  INCLUDE '1SEAICE.INC'

C SET ZERO TO STREAM FUNCTION ARRAYS
       MEMORY_INDEX=0.0

c----------------------------------------------------------------------
       call allocate_shallowater()
       slh     = 0.0 !SEA LEVEL 1
       slh0    = 0.0 !SEA LEVEL 0
       rbottom = 0.0 !BOTTOM STRESS
       slhgrx  = 0.0
       slhgry  = 0.0
       ubrtr   = 0.0
       vbrtr   = 0.0 !BAROTROPIC VELOCITIES
       dinx    = 0.0
       diny    = 0.0
c----------------------------------------------------------------------

       MEMORY_INDEX=MEMORY_INDEX+9*NX*NY*2
C SET ZERO TO HYDRODYNAMICS  ARRAYS
       TT     =0.0
       SS     =0.0  !TEMPERATURE, SALINITY AND
       DEN    =0.0  !DENSITY
       DEN_P  =0.0  !DENSITY
       AGE    =0.0  !DENSITY
       PASS_TRACER  =0.0
       UU     =0.0
       VV     =0.0  !BAROCLINIC VELOCITIES
       XXT    =0.0
       YYT    =0.0
       WW     =0.0  !VERTICAL VELOCITY


       MEMORY_INDEX=MEMORY_INDEX+10*NX*NY*NZ+ NX*NY*(NZ+1)

C SET ZERO TO VISCOSITY AND DIFFUSION  ARRAYS

C    LATERAL VISCOUS AND DIFFUSION COEFFICIENTS FOR T, S, U, V:
C    HORIZONTAL:
       AMXT   =0.0
       AMYT   =0.0
       AMXU   =0.0
       AMYU   =0.0
       PT_DIFF_X=0.0
       PT_DIFF_Y=0.0
       MEMORY_INDEX=MEMORY_INDEX+6*NX*NY*NZ
C    LATERAL DIFFUSION FUNCTION FOR T, S:
       AFNT   =0.0
       AFNV   =0.0
       SLRX   =0.0
       SLRY   =0.0
       MEMORY_INDEX=MEMORY_INDEX+2*NX*NY*(NZ+1)+2*NX*NY*NZ

C    VERTICAL VISCOUS AND DIFFUSION FUNCTIONS
       RIT    =0.0            !RICHARDSON NUMBERs ON T-GRID
       ANZT   =0.0
       ANZU   =0.0                 !VERT. DIFF.&VISC.COEF.

       MEMORY_INDEX=MEMORY_INDEX+NX*NY*NZ+2*NX*NY*(NZ+1)


C SET ZERO TO SURFACE  ARRAYS

       T0        =0.0
       S0        =0.0        !SST, SSS (USED ARRAYS)
       AG0       =0.0
	 PT_FORC   =0.0
       T0B       =0.0
       S0B       =0.0        !SST, SSS (REFERENCE ARRAYS)
       TAUX      =0.0
       TAUY      =0.0        !SEA SURFACE WIND STRESS
       QBAL      =0.0
       QSAL      =0.0        !BALC HEAT, SALINITY FLUX
       QSWR      =0.0        !SW RADIATION
       DIVSWRAD  =0.0
       SICE      =0.0        !SWRAD PENETRATION,SEA ICE
       DKFT      =0.0
       DKFS      =0.0        !T,S RELAXING DIFF.COEFFICIENTS

       TA        =0.0
       QA        =0.0
       PR        =0.0
       WIND      =0.0
       LW        =0.0
       SW        =0.0
       SLPR      =0.0
       UWND      =0.0
       VWND      =0.0
       RUNOFF    =0.0

       MEMORY_INDEX=MEMORY_INDEX+23*NX*NY+NX*NY*NZ

C SET ZERO TO WORK ARRAYS
       WOR    =0.0
       ACCUM  =0.0

       MEMORY_INDEX=MEMORY_INDEX+10*NX*NY

C SET ZERO TO Z- ARRAYS

C       ALEV    =0.0     !ARRAY ON LEVITUS HORIZONTS

       MEMORY_INDEX=MEMORY_INDEX+NX*NY*NLEV+NLEV
       MEMORY_INDEX=MEMORY_INDEX+3*NX*NY    !FOR STORAGE

       HICE   =0.0
       TICE   =0.0
       AICE   =0.0
       HSNOW  =0.0
       TSNOW  =0.0

       AICE0  =0.0

       UICE   =0.0
       VICE   =0.0

       SIGMA1 =0.0
       SIGMA2 =0.0
       S12    =0.0

       MEMORY_INDEX=MEMORY_INDEX+NX*NY*MGRAD*5+NX*NY*6

       WRITE(*,'(1X,A,I10)') ' NUMBER OF ALL GRID POINTS =',MEMORY_INDEX
       MEMORY_INDEX=MEMORY_INDEX*4
       WRITE(*,'(1X,A,I10,A,F10.3,A,F10.6,A)')
     &  ' TOTAL MEMORY =',MEMORY_INDEX,' Bytes,',
     &    FLOAT(MEMORY_INDEX)*2.0**(-20),' Mbytes,',
     &    FLOAT(MEMORY_INDEX)*2.0**(-30),' Gbytes.'

	 WRITE(*,*)

      RETURN
      END

C=======================================================28.04.00 21:56
      SUBROUTINE OCPAR(IGT,IGS)
	  IMPLICIT NONE
C---------------------------------------------------------------------
C DEFINE MODEL PARAMETERS
      INCLUDE '0COM.INC'
      INCLUDE '0CEAN.INC'
      INCLUDE '0FUNCDEF.INC'
C ANUD  - INPUT VALUE OF VERTICAL DIFFUSION COEFFICIENT IN UPPER LAEYR
C ANUZT - INPUT VALUE OF BACKGRAUND VERTICAL DIFFUSION COEFFICIENT FOR T&S
C ANUZU - INPUT VALUE OF BACKGRAUND VERTICAL DIFFUSION COEFFICIENT FOR U&V
C DSWT  - SWITCHING DEPTH FOR TYPE OF LATERAL DIFFUSION (Z ON TOP) [METER]
C DSWV  - SWITCHING DEPTH FOR VALUE OF LATERAL DIFFUSION (E-FOLD)  [METER]
C ZFRAC - RESIDUAL FRACTION IN DEPTH OF VALUE OF LATERAL DIFFUSION [0;1]
      REAL AMUXT,AMUYT,AMUXU,AMUYU,DSWT,DSWV,ZFRAC
	  REAL RLU,RLV,RLT
C+++++++++++++++++++GRID PARAMETER +++++++++++++++++++++++++++++++++
C MLRDEF -- DEFINED SECTIONITY OF OCEAN REGION
C NGRIDT,NGRIDH,NGRIDU,NGRIDXC,NGRIDYC -- DEFINED NUMBER OF OCEAN POINTS
C                                         IN APPROPRIATE MASK
      INTEGER IGT,IGS   ! TYPE OF SEA SURFACE BOUNDARY CONDITION FOR T&S
      INTEGER MLRDEF,NGRIDT,NGRIDH,NGRIDU,NGRIDXC,NGRIDYC
      CHARACTER(1024)
     &  T_MASK_FILE,    !NAME OF FILE WITH TEMPERATURE POINT SEA-LAND MASK
     &  BOTTOM_TOPOGRAPHY_FILE, !NAME OF FILE WITH BOTTOM TOPOGRAPHY
     &  HELP_STRING

	INTEGER  M, N, K, IERR

C--------PARAMETERS FOR SOUTH AND NORTH BOUNDARY Y-DIFFUSION SET-----------
      REAL     RLATMAX,   !MAX LATITUDE (RLAT IS MIN LATITUDE)
     &         RLATLENTH, !LENTH OF INFLUENCE IN DEGRIDS
     &         POINTLAT,  !LATITUDE AT GRID
     &         DIFFUSION_ON_BONDARIES
	PARAMETER(RLATMAX=RLAT+DYST*(NN-NNN),RLATLENTH=2.0,
     &          DIFFUSION_ON_BONDARIES=1E+06)
C--------------------------------------------------------------------------

C ZERO SETTING MOST OF ALL OF BASE OCEAN ARRAYS
      CALL OCZEROSET()

C DEFINE PARAMETERS OF TASK
C DESCRIPTION OF PARAMETERS SEE IN FILE WITH MAME FILEPAR

      OPEN (90,FILE='oceanmodel.par',STATUS='OLD')
      READ (90,*) KSW_INP !SWITCHS FOR TASK:INP.COND
      READ (90,*) KSW_FI  !STREAM FUNCTION
      READ (90,*) KSW_TSL !T&S LATERAL  TRANSPORT&MIX
      READ (90,*) KSW_TSV !T&S VERTICAL TRANSP&MIX
      READ (90,*) KSW_UVL !U&V LATERAL TRANSPORT&MIX
      READ (90,*) KSW_UVV !U&V VERTICAL TRANSPORT&MIX
      READ (90,*) KSW_PTL !PASSIVE TRACER LATERAL TRANSPORT&MIX
      READ (90,*) KSW_PTV !PASSIVE TRACER VERTICAL TRANSPORT&MIX
      READ (90,*) KSW_LQB !NUMBER OF LIQUID BOUNDARIES INCLUDED(0-NO INCL)
      READ (90,*)  RDR,RMU   !BOTTOM FRICTION & ITERATION PARAMETER FOR SF
      READ (90,*)  AMUXT,AMUYT !LATERAL DIFFUSION & VICOSITY
      READ (90,*)  AMUXU,AMUYU !IN X AND Y DIRECTION (2ND ORDER)
      READ (90,*)  CX4T,CY4T   !LATERAL DIFFUSION & VICOSITY
      READ (90,*)  CX4U,CY4U   !IN X AND Y DIRECTION (4TH ORDER)
      READ (90,*)  ANUMAXT,ANUBGRT  !MAX(TOP) & BACKGROUND DIFFUS.VERT.COEF
      READ (90,*)  ANUMAXU,ANUBGRU  !MAX & BACKGROUND VISKOS.VERT.COEF
      READ (90,*) WCLSD, WCLZD, WCLRD  !Weight Coefficient Lateral S,Z,R-Diffusion
      READ (90,*)       WCVZD1,WCVRD1 !Weight Coefficient Vertical Z,R-Diffusion(1stPart)
      READ (90,*) WCVZD2,WCVRD2,WCVZRD !Weight Coefficient Vertical Z,R,ZR-Diffusion(2nd Part)
      READ (90,*) DSWT,DSWV,ZFRAC !SWITCHING DEPTH FOR TYPE&VALUE OF LATERAL DIFFUSION
      HELP_STRING =' '
      READ (90,'(A)') HELP_STRING   ! FILE WITH T-MASK'
	CALL GET_FIRST_LEXEME(HELP_STRING ,T_MASK_FILE   )

      HELP_STRING =' '
	READ (90,'(A)') HELP_STRING  ! FILE WITH BOTTOM TOPOGRAPHY'
	CALL GET_FIRST_LEXEME(HELP_STRING ,BOTTOM_TOPOGRAPHY_FILE  )

      CLOSE(90)

      WRITE(*,'(I7,A)') KSW_INP,' -SWITCHING KEYS FOR TASK:INP.COND'
      WRITE(*,'(I7,A)') KSW_FI ,' -SEA LEVEL(+)/STREAM FUNCTION(-)'
      WRITE(*,'(I7,A)') KSW_TSL,' -T&S LATERAL  TRANSPORT&MIX'
      WRITE(*,'(I7,A)') KSW_TSV,' -T&S VERTICAL TRANSP&MIX'
      WRITE(*,'(I7,A)') KSW_UVL,' -U&V LATERAL TRANSPORT&MIX'
      WRITE(*,'(I7,A)') KSW_UVV,' -U&V VERTICAL TRANSPORT&MIX'
      WRITE(*,'(I7,A)') KSW_PTL,' -PASSIVETRACER LATERAL TRANSPORT&MIX'
      WRITE(*,'(I7,A)') KSW_PTV,' -PASSIVETRACER VERTICAL TRANSPORT&MIX'
      WRITE(*,'(I7,A)') KSW_LQB,' -NUMBER OF LIQUID BOUNDARIES INCLUDED'
      WRITE(*,'(E12.4,F8.4,A)')  RDR,RMU,
     &' BOTTOM FRICTION AND ITERATION PARAMETER FOR STREAM FUNCTION'
      WRITE(*,'(2E12.4,A)')AMUXT,AMUYT,
     &'[CM**2/S] -LATERAL DIFFUSION OF 2nd ORDER IN X & Y DIRECTION'
      WRITE(*,'(2E12.4,A)')AMUXU,AMUYU,
     &'[CM**2/S] -LATERAL VISCOSITY OF 2nd ORDER IN X & Y DIRECTION'
      WRITE(*,'(2E12.4,A)')CX4T,CY4T,
     &'[UNDIM,CM**4/S] -LATERAL DIFFUSION OF 4th ORDER IN X&Y DIRECTION'
      WRITE(*,'(2E12.4,A)')CX4U,CY4U,
     &'[UNDIM,CM**4/S] -LATERAL VISCOSITY OF 4th ORDER IN X&Y DIRECTION'
      WRITE(*,'(2F9.4,A)')ANUMAXT,ANUBGRT,'-MAX(TOP)&BACKGROUND DIFFUS.'
      WRITE(*,'(2F9.4,A)')ANUMAXU,ANUBGRU,'-MAX(TOP)&BACKGROUND VISKOS.'
      WRITE(*,'(3F7.1,A)')DSWT,DSWV,ZFRAC,
     &   '-SWITCHING DEPTHS AND RESIDUAL FRACTION OF LATERAL DIFFUSION'
      WRITE(*,'(A,A)')
     &' FILE WITH TEM. POINT SEA-LAND MASK: ',
     &      T_MASK_FILE(1:LEN_TRIM (T_MASK_FILE))
      WRITE(*,'(A,A)')
     &'        FILE WITH BOTTOM TOPOGRAPHY: ',
     &      BOTTOM_TOPOGRAPHY_FILE(1:LEN_TRIM (BOTTOM_TOPOGRAPHY_FILE))

C----TYPES OF BOUNDARY CONDITIONS FOR U AND V--------------------------

c      IGRYU = 1    ! NON-SLIP CONDITION FOR U ON Y-WALLS
       IGRYU = 2    !     SLIP CONDITION FOR U ON Y-WALLS
       IGRZU = 2    ! 2-ND CONDITION FOR U ON SEA SURFACE

c              IGRXV = 1    ! NON-SLIP CONDITION FOR V ON X-WALLS
               IGRXV = 2    !     SLIP CONDITION FOR V ON X-WALLS
               IGRZV = 2    ! 2-ND CONDITION FOR V ON SEA SURFACE
C----------------------------------------------------------------------

       IGRZT  = MIN(iabs(IGT),2) ! TYPE OF CONDITION FOR T ON SEA SURFACE.
       IGRZS  = MIN(iabs(IGS),2) ! TYPE OF CONDITION FOR S ON SEA SURFACE.
       IGRZAGE= 1

C----------------------------------------------------------------------

C SETTING VERTICAL T-,W- GRID LEVELS
C IF INPUT INTEGER NUMBER EQUALS
C +-1 - THEN W-LEVELS ARE ARRANGED IN THE MIDDLES OF T-LAYERS
C +-2 - THEN T-LEVELS ARE ARRANGED IN THE MIDDLES OF W-LAYERS
C > 0 - THEN LEVELS ARE ANALITICALLY SETTED
C < 0 - THEN LEVELS ARE DIRECT TAKEN FROM LEVITUS LEVELS
      CALL VGRID()

C DEFINE GRID STEPS AND CORIOLIS PARAMETER
      CALL BASINPAR()
C GRID PARAMETER SETTING
      CALL GRIDCON(T_MASK_FILE,
     &             MLRDEF,NGRIDT,NGRIDH,NGRIDU,NGRIDXC,NGRIDYC)
	  INDEXTV=0
      NVECTD=0

      DO N=NNN,NN
        DO M=MMM,MM
           IF(LU(M,N).GT.0.5) THEN
              NVECTD=NVECTD+1
              INDEXTV(M,N)=NVECTD
           END IF
        END DO
      END DO
      WRITE(*,'(A,I7)')'   TOTAL VECTOR DIMENSION ON T-GRID =',NVECTD

C CALCULATION OF H-GRID VECTOR NUMBERING AND ITS TOTAL VECTOR DIMENSION

      INDEXHV=0
	NVECHD=0

C      DO M=MMM,MM
C        DO IG=1,LRYU(M)
C         II = IIYU(IG,M)
C         JJ = JJYU(IG,M)
C         DO N=II,JJ
C           NVECHD=NVECHD+1
C           INDEXHV(M,N)=NVECHD
C	   END DO
C        END DO
C      END DO

      DO N=NNN,NN
        DO M=MMM,MM
           IF(LUU(M,N).GT.0.5) THEN
              NVECHD=NVECHD+1
              INDEXHV(M,N)=NVECHD
           END IF
        END DO
      END DO

      WRITE(*,'(A,I7)')'   TOTAL VECTOR DIMENSION ON H-GRID =',NVECHD

C DEFINE BOTTOM TOPOGRAPHY IN METERS
C READ BOTTON TOPOGRAPHY ON T-GRID
c      CALL RDSTD(' ',BOTTOM_TOPOGRAPHY_FILE,1,HH,LU,NX,NY,1,
c     &                      MMM,MM,NNN,NN,1,1,IERR)
c      WRITE(*,'(A,F9.2)')'   MEAN DEPTH  [M]:',WAVER(HH)

C DEFINE BOTTOM TOPOGRAPHY IN METERS
c      CALL RDSTD(' ',BOTTOM_TOPOGRAPHY_FILE,1,HH,LUH,NX,NY,1,
c     &                      MMM-1,MM,NNN-1,NN,1,1,IERR)
c      WRITE(*,'(A,F9.2)')'   MEAN DEPTH  [M]:',AAAVEH(HH,1)

      HH = 1000.0

C  CONVERS BOTTOM UNIT FROM METER TO CENTIMETER
      CALL FF1AR (HH,100.0,0.0,HH,LUH,NX,NY,1)
      IF(MMD.NE.0) CALL CYCLIZE(HH ,NX,NY,1,MMM,MM)
      WRITE(*,'(A,F9.2)')'   MEAN DEPTH [CM]:',AAAVEH(HH,1)

C INTERPOLATING HH GIVEN ON H-GRID(LUH) TO HHQ GIVEN ON T-GRID(LU).
C INTERPOLATION: HH -> HHQ. GRID TYPE NTO=1: T->U; NTO=-1: U->T.
      CALL EXTUD(HH,LUH,HHQ,LU,0.0,-1,NX,NY,1)
      IF(MMD.NE.0) CALL CYCLIZE(HHQ,NX,NY,1,MMM,MM)


      DO N=2,NY-1
       DO M=2,NX-1
C INTERPOLATING HH GIVEN ON H-GRID(LUH) TO HHU GIVEN ON U-GRID(LCU).
         HHU(M,N)=(HH(M,N)+HH(M,N-1))/2.0
C INTERPOLATING HH GIVEN ON H-GRID(LUH) TO HHV GIVEN ON V-GRID(LCV).
         HHV(M,N)=(HH(M,N)+HH(M-1,N))/2.0
       END DO
	END DO

      IF(MMD.NE.0) THEN
        CALL CYCLIZE(HHU,NX,NY,1,MMM,MM)
        CALL CYCLIZE(HHV,NX,NY,1,MMM,MM)
      END IF

      INCLUDE '1LATDIFF.INC'


C COEFFICIENTS ON EXTERNAL BOUNDARIES
      DO N = 1,NNN-1
       DO M = 1,NX
        DO K=1,NZ
         AMXT(M,N,K)=AMXT(M,NNN,K)
         AMYT(M,N,K)=AMYT(M,NNN,K)
         AMXU(M,N,K)=AMXU(M,NNN,K)
         AMYU(M,N,K)=AMYU(M,NNN,K)
        ENDDO
       END DO
      ENDDO

      DO N = NN+1,NY
       DO M = 1,NX
        DO K=1,NZ
         AMXT(M,N,K)=AMXT(M,NN,K)
         AMYT(M,N,K)=AMYT(M,NN,K)
         AMXU(M,N,K)=AMXU(M,NN,K)
         AMYU(M,N,K)=AMYU(M,NN,K)
        ENDDO
       END DO
      ENDDO

      DO N = 1,NY
       DO M = 1,MMM-1
        DO K=1,NZ
         AMXT(M,N,K)=AMXT(MMM,N,K)
         AMYT(M,N,K)=AMYT(MMM,N,K)
         AMXU(M,N,K)=AMXU(MMM,N,K)
         AMYU(M,N,K)=AMYU(MMM,N,K)
        ENDDO
       END DO
      ENDDO

      DO N = 1,NY
       DO M = MM+1,NX
        DO K=1,NZ
         AMXT(M,N,K)=AMXT(MM,N,K)
         AMYT(M,N,K)=AMYT(MM,N,K)
         AMXU(M,N,K)=AMXU(MM,N,K)
         AMYU(M,N,K)=AMYU(MM,N,K)
        ENDDO
       END DO
      ENDDO

      IF(IABS(KSW_TSL).EQ.4) THEN
	   AMXT=AMUXT
	   AMYT=AMUYT
	END IF

	IF(IABS(KSW_UVL).EQ.4) THEN
	   AMXU=AMUXU
	   AMYU=AMUYU
	END IF

      IF(MMD.NE.0) THEN
        CALL CYCLIZE (AMXT,NX,NY,NZ,MMM,MM)
        CALL CYCLIZE (AMYT,NX,NY,NZ,MMM,MM)
        CALL CYCLIZE (AMXU,NX,NY,NZ,MMM,MM)
        CALL CYCLIZE (AMYU,NX,NY,NZ,MMM,MM)
	END IF


      IF(KSW_TSV.LT.0)  THEN
c solar shortwave heating penetrate divergence function.
C => Shortwave penetration is a double exponential as follows:
         CALL SWRABSOR(DIVSWRAD,HHQ,LU,ZW,NX,NY,NZ)
                write(*,*)
     &'  Shortwave penetration divergent coefficients have been defined'
	ELSE
                       DIVSWRAD = 0.0
      END IF

      IF(IABS(KSW_TSL).GT.2)  THEN
         CALL VERTFUN( AFNT, AFNV,DSWT,DSWV,ZFRAC,1)
               WRITE(*,*)
     &'  Setting vertical switching function for lateral diffusion.'
      END IF

C---------------------------------------------------------------------72
C DEFINITION OF 3D FIXED COEFFICIENTS FOR VERTICAL DIFFUSION AND VISKOSITY
      CALL DIFCFTU3D(ANZT, ANZU,LU,
     &               ANUMAXT,ANUBGRT,ANUMAXU,ANUBGRU,NX,NY,NZ)

C---------------------------------------------------------------------72
C DEFINITION OF COEFFICIENTS FOR SIMMETRIC FILTRATION OVER NORTH POLE
C      ALONG LATITUDE CIRCLE
      CALL NPSF_COEF_DEF()

C-----------LIST OF SOME PARAMETERS OF OGCM---------------------------72
C CASE WITH ISLAND(S) -- DEFINITION OF ISLAND STREAM FUNCTION
      WRITE(*,'('' RDR;RMU;EPS='',5E10.3)') RDR,RMU,EPS
      WRITE(*,*) '   LATERAL DIFFUSION TYPE:'
      WRITE(*,'(''   WCLSD, WCLZD, WCLRD  ='',3F8.5)')
     &               WCLSD, WCLZD, WCLRD
      WRITE(*,'(''         WCVZD1,WCVRD1  ='',8X,2F8.5)')
     &                     WCVZD1,WCVRD1
      WRITE(*,'(''  WCVZD2,WCVRD2,WCVZRD  ='',8X,3F8.5)')
     &              WCVZD2,WCVRD2,WCVZRD

      WRITE(*,'('' TT & SS: AMXT,AMYT ='',2G12.4)')
     &               AMXT(MMM,NNN,1),AMYT(MMM,NNN,1)
      WRITE(*,'('' UU & VV: AMXU,AMYU ='',2G12.4)')
     &               AMXU(MMM,NNN,1),AMYU(MMM,NNN,1)
      WRITE(*,'('' ANUMAXT,ANUBGRT  ='',2F12.4,A)')
     &             ANUMAXT,ANUBGRT,' -MAX(TOP)&BACKGR.TS VERT.DIFFUS.'
      WRITE(*,'('' ANUMAXU,ANUBGRU  ='',2F12.4,A)')
     &             ANUMAXU,ANUBGRU,' -MAX(TOP)&BACKGR. VERT.VISKOSITY'


      RETURN
      END
C======================================================================
      SUBROUTINE VERTFUN(AFNT, AFNV, DSW_TYPE,DSW_VALUE,ZFRACTION,IFIG)
	IMPLICIT NONE
      INCLUDE '0COM.INC'
	INCLUDE '0FUNCDEF.INC'
      REAL     AFNT(NX,NY,NZ+1), AFNV(NX,NY,NZ+1)
	REAL  XPOINT,YPOINT
c	PARAMETER(XPOINT=-37.0,YPOINT=53.5)  !COORDINATES of "C"point FOR VERT.FUN DEMO
	PARAMETER(XPOINT=+323.0,YPOINT=53.5) !COORDINATES of "C"point FOR VERT.FUN DEMO
C DSW_TYPE - swithing depth of type  of lateral diffusion [METERS]
C DSW_VALUE - swithing depth of value of its coefficient [METERS]
C ZFRACTION - fraction of value of its coefficient in depth
C IFIG - PLOTING PARAMETER(0 - NOT,1 - PLOT)
      REAL DSW_TYPE,DSW_VALUE, ZFRACTION
      INTEGER IFIG
      CHARACTER*80 TXT_FIG         !FOR OUTPUT GRAPHICS AS TEXT
      INTEGER   IWIDTH_FIG
      PARAMETER(IWIDTH_FIG=20)
      INTEGER   M, N, K, L, IC, JC
	REAL      AV, AS, AT
C Setting vertical function for automatic switching lateral diffusion
c from horizontal to isopicnal type
C and decreasing with depth function of coefficient of lateral diffusion
c
c    1|_______  ______________________________
c     |       \/                              |
c  0.5|       /\                              |
c     |      /| \_____________________________| ZFRACTION
c    0|_____/ ||dsw_VALUE/H                   |
c     |_______||______________________________|
c     0    dsw_type/H    zw(k)                1

      DO K = 1,NZ+1
	    DO N = NNN-1,NN+1
              DO M = MMM-1,MM+1
          		IF(LU(M,N).GT.0.5) THEN

		        AS=(0.025-0.005)*EXP(-ZW(K)*HHQ(M,N)*0.0001)+0.0045
                  AT=0.5*(TANH(AS*(ZW(K)*HHQ(M,N)*0.01-DSW_TYPE))+1.0)
		        AV=ZFRACTION+(1.0-ZFRACTION)*
     *                 EXP(-ZW(K)*HHQ(M,N)/(100.0*DSW_VALUE))
                  AFNT(M,N,K) = AT
                  AFNV(M,N,K) = AV
	            END IF

	        ENDDO
	    ENDDO
      ENDDO

      IF (IFIG.NE.0) THEN

C TEXT GRAPHIC OF VERTICAL SWITCH FUNCTION ON OWS CHARLY
          WRITE(*,*)
          WRITE(*,*)'   VERTICAL SWITCH FUNCTIONS IN VICINITY OF POINT'
          WRITE(*,'(2(A,F10.5))')'      LON =',XPOINT,'   LAT =',YPOINT
         IC=INT((XPOINT-RLON)/DXST)+1
	 IC=MIN(IC,MM)      !CORRECTION IF NEED
	 IC=MAX(IC,MMM)     !CORRECTION IF NEED

         JC=INT((YPOINT-RLAT)/DYST)+1
	 JC=MIN(JC,NN)      !CORRECTION IF NEED
	 JC=MAX(JC,NNN)     !CORRECTION IF NEED

C     &     '   VERTICAL SWITCH FUNCTION IN MIDDLE POINT OF INDIAN OCEAN'
C          IC=INT(( 80.0-RLON)/DXST)+1
C          JC=INT((-20.0-RLAT)/DYST)+1
          WRITE(*,'(1X,A,F7.2,A,F7.2,A)')
     & '    WITH COORDINATES OF ',RLON+FLOAT(IC-1)*DXST,' W, ',
     &      RLAT+FLOAT(JC-1)*DYST,' N.'
          WRITE(*,'(1X,A,F7.1)')
     &'                       SWITCH DEPTH FOR TYPE OF DIFFUSION:',
     &                 DSW_TYPE
          WRITE(*,'(1X,A,2F7.1)')
     &'    SWITCH DEPTH&RESIDUAL FRACTION FOR VALUE OF DIFFUSION:',
     &                 DSW_VALUE,ZFRACTION
          DO L = 1,IWIDTH_FIG+15
             TXT_FIG(L:L)=' '
          ENDDO

          TXT_FIG(1 : 8)='Depth[m]'
          TXT_FIG(10:10)='0'
          TXT_FIG(IWIDTH_FIG/4  +10-1:IWIDTH_FIG/4  +10+2)='0.25'
          TXT_FIG(IWIDTH_FIG/2  +10-1:IWIDTH_FIG/2  +10+2)='0.50'
          TXT_FIG(IWIDTH_FIG/4*3+10-1:IWIDTH_FIG/4*3+10+2)='0.75'
          TXT_FIG(IWIDTH_FIG    +10  :IWIDTH_FIG    +10  )='1'
          WRITE(*,'(1X,A)') TXT_FIG(1:IWIDTH_FIG+15)
          DO L = 1,IWIDTH_FIG+15
              TXT_FIG(L:L)='-'
          ENDDO
          TXT_FIG(10:10)='+'
          TXT_FIG(IWIDTH_FIG/4  +10:IWIDTH_FIG/4  +10)='+'
          TXT_FIG(IWIDTH_FIG/2  +10:IWIDTH_FIG/2  +10)='+'
          TXT_FIG(IWIDTH_FIG/4*3+10:IWIDTH_FIG/4*3+10)='+'
          TXT_FIG(IWIDTH_FIG    +10:IWIDTH_FIG    +10)='+'
          WRITE(*,'(1X,A)') TXT_FIG(1:IWIDTH_FIG+15)
          DO K = 1,NZ+1,2
              DO L = 1,IWIDTH_FIG+15
                  TXT_FIG(L:L)=' '
              ENDDO
              WRITE(TXT_FIG(1:9),'(F7.1,2H |)') ZW(K)*HH(IC,JC)*1E-02
              AT = AFNT(IC,JC,K)
              AV = AFNV(IC,JC,K)
              L=NINT(AT*IWIDTH_FIG)+10
              WRITE(TXT_FIG(L:L),'(A1)') '*'
              L=NINT(AV*IWIDTH_FIG)+10
              WRITE(TXT_FIG(L:L),'(A1)') '+'
              WRITE(*,'(1X,A)') TXT_FIG(1:IWIDTH_FIG+15)
          ENDDO
          WRITE(*,*)
      END IF

      RETURN
      END
C=========================================================27.11.99 17:38
      SUBROUTINE SWRABSOR(DIVSWRAD,HHQ,LU,ZW,NX,NY,NZ)
	  IMPLICIT NONE
      INTEGER NX, NY, NZ,M,N,K
      REAL DIVSWRAD(NZ,NX,NY),HHQ(NX,NY),LU(NX,NY),ZW(NZ+1)
c-----------------------------------------------------------------------
c     Solar Shortwave energy penetrates below the ocean surface. Clear
c     water assumes energy partitions between two exponentials as
c     follows:
c
c     58% of the energy decays with a 35 cm e-folding scale
c     42% of the energy decays with a 23 m  e-folding scale
c
c     Paulson and Simpson (1977 Irradiance measurements in the upper
c                               ocean JPO 7, 952-956)
c     Also see ... Jerlov (1968 Optical oceanography. Elsevier)
c                  A General Circulation Model for Upper Ocean
c                  Simulaton (Rosati and Miyakoda JPO vol 18,Nov 1988)
c-----------------------------------------------------------------------
c
C => Shortwave penetration is a double exponential as follows:
C     PARAMETER(RPART = 0.58,EFOLD1 = 35.0,EFOLD2 = 23.0E+02)
C NEW APPROUCH OF USING SHORT WAVE RADIATION: UPPER PART
C OF ABOUT 60% ADDED TO HEAT FLUX AND RESIDUAL PART
C OF 40% OF THE ENERGY DECAYS WITH A 20 M E-FOLDING SCALE
      REAL      RPART, EFOLD1, EFOLD2
      PARAMETER(RPART = 0.0,EFOLD1 = 35.0,EFOLD2 = 20.0E+02)
      REAL SWARG1,SWARG2,PEN1,PEN2,SPACESUM,SUMK,POINTSUM

      SPACESUM=0.0
      POINTSUM=0.0

      DO N = 2,NY-1
      DO M = 2,NX-1

      IF(LU(M,N).GT.0.5) THEN

      POINTSUM=POINTSUM+1.0

      PEN2=1.0
      SUMK=0.0
      DO K=1,NZ
        SWARG1 = -MIN(ZW(K+1)*HHQ(M,N)/EFOLD1,70.0)
        SWARG2 = -MIN(ZW(K+1)*HHQ(M,N)/EFOLD2,70.0)
        PEN1 = PEN2
        PEN2 = RPART*EXP(SWARG1) + (1.0-RPART)*EXP(SWARG2)
        DIVSWRAD(K,M,N) =(PEN1 - PEN2)/(ZW(K+1)-ZW(K))/HHQ(M,N)
        SUMK=SUMK+PEN1-PEN2
      ENDDO

      DO K=1,NZ
        DIVSWRAD(K,M,N) = DIVSWRAD(K,M,N)/SUMK
        SPACESUM=SPACESUM+DIVSWRAD(K,M,N)*(ZW(K+1)-ZW(K))*HHQ(M,N)
      ENDDO

      END IF
      ENDDO
      ENDDO

      WRITE(*,'(2X,A,F10.2)')
     &'SUM OF SW DIVERGENCE COEFFICIENT =',SPACESUM
      WRITE(*,'(2X,A,F10.2,A)')'FOR ',POINTSUM,' T-GRID POINTS'

      RETURN
      END
C======================================================================
      SUBROUTINE DIFCFTU3D(ANZT,ANZU,LU,
     &                ANUMAXT,ANUBGRT,ANUMAXU,ANUBGRU,
     &                         NX,NY,NZ)
      IMPLICIT NONE
C DEFINITION OF 3D COEFFICIENTS FOR VERTICAL DIFFUSION
      INTEGER NX,NY,NZ
      REAL ANZT(NX,NY,NZ+1), !COEFICIENT OF DIFFUSION FOR T AND S
     &     ANZU(NX,NY,NZ+1)  !COEFICIENT OF VISKOSITY FOR U AND V
      REAL ANUMAXT,ANUBGRT,ANUMAXU,ANUBGRU  !MAX(TOP)&BACKGR. TS&U VERT.COEF
      REAL LU(NX,NY)
	INTEGER M, N, K

      DO N=1,NY
         DO M=1,NX
C SET UPPER COEFFICIENTS FOR TEMPERATURE VERTICAL DIFFUSION
            ANZT(M,N,1) = ANUMAXT*LU(M,N)
C SET UPPER COEFFICIENTS FOR MOMENTUM VERTICAL VISKOSIITY
            ANZU(M,N,1) = ANUMAXU*LU(M,N)
         ENDDO
      ENDDO

      DO K=2,NZ+1
      DO N=1,NY
         DO M=1,NX
C SET COEFFICIENTS FOR TEMPERATURE VERTICAL DIFFUSION
            ANZT(M,N,K) = ANUBGRT*LU(M,N)
C SET COEFFICIENTS FOR VELOCITY VERTICAL VISKOSIITY
            ANZU(M,N,K) = ANUBGRU*LU(M,N)
         ENDDO
      ENDDO
      ENDDO

      RETURN
      END
C=================================================================
      SUBROUTINE ICEPAR()
      IMPLICIT NONE

      INCLUDE '0COM.INC'
      INCLUDE '0CEAN.INC'
      INCLUDE '0FUNCDEF.INC'

      OPEN (90,FILE='seaicemodel.par',STATUS='OLD')
      READ (90,*) KSW_IML   ! LATERAL  CHANGE OF ICE-SNOW MASS AND COMPACTNESS (TRANSPORT)
      READ (90,*) KSW_IMV   ! VERTICAL CHANGE OF ICE-SNOW MASS AND COMPACTNESS (THERMODYNAMICS)
      READ (90,*) KSW_IMI   ! INTERNAL CHANGE OF ICE-SNOW MASS AND COMPACTNESS (REDISTRIBUTION)
      READ (90,*) KSW_IVL   ! LATERAL  CHANGE OF ICE VELOCITY (TRANSPORT-DIFFUSION)
      READ (90,*) KSW_IVV   ! VERTICAL CHANGE OF ICE VELOCITY (WIND AND WATER STRESS)
      READ (90,*) KSW_IVA_SLG, KSW_IVA_COR  ! ADAPTATION      OF ICE VELOCITY (ATMPRESSURE GRADIENT + SEALEVEL & CORIOLIS)
      READ (90,*) KSW_IVR   ! RHEOLOGY        OF ICE VELOCITY
	CLOSE(90)

      WRITE(*,'(I7 ,A)') KSW_IML, ' LATERAL  CHANGE OF ICE-SNOW MASS'
      WRITE(*,'(I7 ,A)') KSW_IMV, ' VERTICAL CHANGE OF ICE-SNOW MASS'
	WRITE(*,'(I7 ,A)') KSW_IMI, ' INTERNAL CHANGE OF ICE-SNOW MASS'
	WRITE(*,'(I7 ,A)') KSW_IVL, ' LATERAL  CHANGE OF ICE VELOCITY'
	WRITE(*,'(I7 ,A)') KSW_IVV, ' VERTICAL CHANGE OF ICE VELOCITY'
	WRITE(*,'(2I3,A)') KSW_IVA_SLG,
     &                  KSW_IVA_COR, ' ADAPTATION OF ICE VELOCITY'
	WRITE(*,'(I7,A)') KSW_IVR, ' RHEOLOGY        OF ICE VELOCITY'

      RETURN
	END
