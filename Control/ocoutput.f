
      SUBROUTINE XY_OUTPUT(PATH2DATA,NUMOFREC)
      use mpi_parallel_tools
      use mod_shallowater
      implicit none

C WRITING 3D XY-FILD DATA TO NUMOFREC
      CHARACTER*(*) PATH2DATA    !PATH TO OUTPUT FILES
      INTEGER NUMOFREC           !NUMBER OF RECORD
      INCLUDE '0COM.INC'
c      INCLUDE '0CEAN.INC'
      INCLUDE '1LREC.INC'        !SET LONG OF UNIQUE RECL AND UNDEF VALUE
      INCLUDE '1STDOUT.INC'
      INCLUDE '2STDOUT.INC'

      REAL, ALLOCATABLE::  UTR(:,:),VTR(:,:)

      INTEGER  IERR,M,N,K

c	  allocate(UTR(NX,NY), VTR(NX,NY))

c--------- ! OUTPUT ON MODEL GRID ----------------------------------------------
            allocate(vtr(nx_start:nx_end, ny_start:ny_end))
            vtr(nx_start:nx_end, ny_start:ny_end) =
     &            sngl(slh(nx_start:nx_end, ny_start:ny_end))

c            print *, vtr(nx_start:nx_end, ny_start:ny_end)

c            call wdstd_parallel(path2data,'XY/sl.dat',NUMOFREC,
c     &                          vtr,LU,NX,NY,
c     &                          MMM,MM,NNN,NN)
            call PWDSTD(path2data,'XY/sl.dat',NUMOFREC,
     &              vtr, LU,
     &              NX,nx_start,nx_end,
     &              NY,ny_start,ny_end,1,
     &              MMM,MM,NNN,NN,1,1,IERR)

c           if (rank .eq. 0) then
c               vtr = 0
c               do n = ny_start, ny_end
c                   do m = nx_start, nx_end
c                       vtr(m, n) = sngl(slh(m,n))
c                   enddo
c               enddo
c
c               call WDSTD(PATH2DATA,'XY/sl.dat',NUMOFREC,VTR,LU,NX,NY,1,
c     &                MMM,MM,NNN,NN,1,1,IERR)
c           endif

c        CALL STREAM_FUNCTION_CALC(UBRTR,STRFUN)


c        CALL WDSTD(PATH2DATA,'XY/sf.dat',NUMOFREC,STRFUN,LUH,NX,NY,1,
c     &      MMM-1,MM,NNN-1,NN,1,1,IERR)

c       IF(IABS(KSW_TSL).GT.1.OR.IABS(KSW_TSV).GT.1) THEN
C USE ARRAY VTR FOR WRITE MIXED LAYER DEPTH
c        CALL MLDEF(TT,SS,DEN_P,VTR)
c        CALL WDSTD(PATH2DATA,'XY/mxld.dat',NUMOFREC,VTR,
c     &            LU,NX,NY,1,MMM,MM,NNN,NN,1,1,IERR)
c       END IF

!STERICAL LEVEL CALCULATE
c       UTR =0.0

!$OMP PARALLEL DO
c          DO N = NNN,NN
c            DO M = MMM,MM
c      	  IF(lu(M,N).GT.0.5) THEN
c               DO K = 1,NZ
c                UTR(M,N) =UTR(M,N) - DZ(K)*(DEN(M,N,K)/RH0)*HHQ(M,N)       !USE DEN AS HELP ARRAY
c               ENDDO
c              END IF
c	     ENDDO
c         ENDDO
!$OMP END PARALLEL DO

c        CALL WDSTD(PATH2DATA,'XY/sls.dat',NUMOFREC,UTR,LU,NX,NY,1,
c     &      MMM,MM,NNN,NN,1,1,IERR)

      deallocate(vtr)

      RETURN
      END
C=======================================================================

C======================================================================
C (C) DIANSKY N.A.,16.10.99 12:32. (dinar@inm.ras.ru)
C WRITING CONTROL POINTS(PREFIX "cp") AND OTHERS ARRAYS
      SUBROUTINE OCPWRITE(PATH2OCP,IWRITE,IGT,IGS,KALCSTMON)
	  use mpi_parallel_tools
	  use mod_shallowater
	  IMPLICIT NONE
      INCLUDE '0COM.INC'
      INCLUDE '0CEAN.INC'
	INCLUDE '1SEAICE.INC'
      INCLUDE '1LREC.INC'       !SET LONG OF UNIQUE RECL
      INCLUDE '0FUNCDEF.INC'
      INCLUDE '0ICEPAR.INC'
      CHARACTER*(*) PATH2OCP      !PATH TO CONTROL POINT
      INTEGER IWRITE,   !=1 THEN BRIEF WRITE;=2 THEN FULL WRITE
     &        KALCSTMON ! CALCULATOR FOR MONTHLY MEANS
      INTEGER IGT,IGS   ! TYPE OF SEA SURFACE BOUNDARY CONDITION FOR T&S
C-----------------------------------------------------------------------
      REAL       CMPS_MMPDAY, GCMPS_WPM2
      PARAMETER (CMPS_MMPDAY=8.64E+05, !CONVERSION FACTOR [CM/S]=>[MM/DAY]
     &      GCMPS_WPM2=CPWRH0*1.0E+04) !CONVERSION FACTOR [GRAD*CM/S]=>[W/M**2]

      REAL ZNZP1(NZ+1)              !OUTPUT Z LEVELS
C ARRAYS FOR MERIDIONAL STREAM FUNCTION DESIGN
      REAL FFZON(NY)
      REAL TTA(NZ),SSA(NZ),DDA(NZ)   !AQUATORY MEAN TEM.,SAL.,DENSITYT
      CHARACTER(2048) FNAME
	CHARACTER*256 FORCING_TITLE
      INTEGER      M, N, K, IERR, LBAS, NDISPC
      REAL    STARTIME,TSTEP   !START-TIME AND INCREMENT IN SECONDS
      REAL SENS_HEAT,LAT_HEAT,SATUR_HUM,SATUR_PRESS_OCEAN
	REAL MSK(NX,NY)
	MSK=1.0               !IMAGINE REGION MASK IS EQUAL TO 1.0 EVERYWHERE

      STARTIME=0.0
      TSTEP=3600.0

C SET LEVELS IN CP FILES AS HORIZONT NUMBERS
      DO K = 1,NZ+1
       ZNZP1(K)=FLOAT(K)
      ENDDO

      NDISPC=(NN-NNN)/79
      IF(NDISPC.LT.1) NDISPC=1
C DEFINE VERTICAL VELOCITY FOR DIFFUSIVE OF T&S

c      CALL WDSTD(PATH2OCP,'cptt.dat', 1,TT,  LU ,NX,NY,NZ,
c     &                       MMM,MM,NNN,NN,1,NZ,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'cptt.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,NZ,1,
c     &     XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &                    ZNZP1,'Model Potential Temperature[C]','tt')

c      CALL WDSTD(PATH2OCP,'cpss.dat', 1,SS,  LU ,NX,NY,NZ,
c     &                       MMM,MM,NNN,NN,1,NZ,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'cpss.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,NZ,1,
c     &     XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &                    ZNZP1,'Model Salinity -35PPt','ss')

c      CALL WDSTD(PATH2OCP,'cpage.dat', 1,AGE,  LU ,NX,NY,NZ,
c     &                       MMM,MM,NNN,NN,1,NZ,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'cpage.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,NZ,1,
c     &     XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &                    ZNZP1,'Water ideal age, days','age')

c      CALL WDSTD(PATH2OCP,'cpuu.dat', 1,UU, LCU,NX,NY,NZ,
c     &      MMM-1,MM,NNN  ,NN,1,NZ,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'cpuu.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+2,NN-NNN+1,NZ,1,
c     &          XU(MMM-1:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &               ZNZP1,'Zonal velocity[cm/s]','u')

c      CALL WDSTD(PATH2OCP,'cpvv.dat', 1,VV, LCV,NX,NY,NZ,
c     &      MMM  ,MM,NNN-1,NN,1,NZ,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'cpvv.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+2,NZ,1,
c     &          XT(MMM:MM),0.0,YV(NNN-1:NN),0.0,STARTIME,TSTEP,
c     &               ZNZP1,'Meridional velocity[cm/s]','v')


      IF (KSW_FI.GE.1) THEN   !FOR SEA LEVEL ALG

      CALL WDSTD8(PATH2OCP,'cpslh8.dat', 1,SLH , LU ,NX,NY,1,
     &                       MMM,MM,NNN,NN,1,1,IERR)
      CALL WDSTD8(PATH2OCP,'cpslh8.dat', 2,SLH0, LU ,NX,NY,1,
     &                       MMM,MM,NNN,NN,1,1,IERR)

		XXT(:,:,1)=SNGL(SLH (:,:))
 		XXT(:,:,2)=SNGL(SLH0(:,:))

      CALL WDSTD(PATH2OCP,'cpslh.dat', 1,XXT, LU ,NX,NY,2,
     &                       MMM,MM,NNN,NN,1,2,IERR)
      CALL FULFNAME(FNAME,PATH2OCP,'cpslh.dat',IERR)
      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,2,1,
     &      XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
     &               ZNZP1,'Sea Level Height[cm]','sl')

      END IF

c	XXT(:,:,1)=SNGL(UBRTR(:,:))
c      CALL WDSTD(PATH2OCP,'ubc.dat', 1,XXT(:,:,1),LCU,NX,NY,1,
c     &                       MMM-1,MM,NNN,NN,1,1,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'ubc.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+2,NN-NNN+1,1,1,
c     &      XU(MMM-1:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &               ZNZP1,'Barotropic veloscity[cm/c]','u')

c	YYT(:,:,1)=SNGL(VBRTR(:,:))
c      CALL WDSTD(PATH2OCP,'vbc.dat', 1,YYT(:,:,1),LCV,NX,NY,1,
c     &                       MMM,MM,NNN-1,NN,1,1,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'vbc.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+2,1,1,
c     &      XT(MMM:MM),0.0,YV(NNN-1:NN),0.0,STARTIME,TSTEP,
c     &               ZNZP1,'Barotropic veloscity[cm/c]','v')


C	WRITE MONTHLY MEAN COUNTER
c	    CALL FULFNAME(FNAME,PATH2OCP,'cpmmc.txt',IERR)
c          OPEN (30,FILE=FNAME,ERR=125)
c	    WRITE(30,*,ERR=126) KALCSTMON
c	    WRITE(*,*) 'WRITE MONTHLY MEAN COUNTER IN cpmmc.txt'
c	    CLOSE(30)

C     WRITE ACCUMULATED MONTHLY MEANS
c	    CALL WDSTD(PATH2OCP,'cpmmean.dat', 1,ACCUM,MSK,NX,NY,9,
c     &             1,NX, 1,NY, 1,9,IERR)
c          CALL FULFNAME(FNAME,PATH2OCP,'cpmmean.dat',IERR)
c          CALL CREATE_CTLFILE(FNAME,UNDEF,NX,NY,9,1,
c     &        XT,0.0,YT,0.0,STARTIME,TSTEP,
c     &                    ZNZP1,'MONTHLY MEANS ACCUMULATED','mm')

c      IF (IWRITE.EQ.1) RETURN   !EXIT FOR BRIEF WRITE


C      IF (IGT.GT.1) THEN
C DEFINITION TYPE OF SEA SURFACE CONDITION FOR TEMPERATURE

c	IF(NZ.LT.20) THEN
c	  WRITE(*,'(A)')' ERROR IN WRITE sfl.dat IN SBR. OCPWRITE:'
c        WRITE(*,'(A,I3,A)')' NZ=',NZ,' < 20(NUMBER OF OUTPUT VARIABLES)'
c	END IF

!$OMP PARALLEL DO PRIVATE(M,N)
c      DO N = NNN,NN
c         DO M = MMM,MM
c         IF(LU(M,N).GE.0.5) THEN
C             IF(IGRZT(M,N).EQ.2) THEN
C FLUXES FOR TEMPERATURE:
C DEFINITION DOWNWARD LW-RADIATION W/M**2
c         XXT(M,N,1)=0.001*LW(M,N)
C DEFINITION DOWNWARD SW-RADIATION W/M**2
c         XXT(M,N,2)=0.001*SW(M,N)
C DEFINITION TEMPERATURE OF ATMOSPHERE, GRAD C
c         XXT(M,N,3)=TA(M,N)
C DEFINITION HUMIDITY OF ATMOSPHERE, KG/KG
c         XXT(M,N,4)=QA(M,N)
C DEFINITION SEA LEVEL PRESSURE, PA
c         XXT(M,N,5)=SLPR(M,N)
C DEFINITION BALANCE HEAT FLUX FOR W/M**2
c         XXT(M,N,6)= GCMPS_WPM2*QBAL(M,N)
C DEFINITION SHORT WAVE RADIATION FOR W/M**2
c         XXT(M,N,7)= GCMPS_WPM2*QSWR(M,N)
C DEFINITION NET HEAT FLUX FOR W/M**2
c         XXT(M,N,8)= GCMPS_WPM2*T0(M,N)
C DEFINITION RELAX HEAT FLUX FOR W/M**2
c         XXT(M,N,9)= DKFT(M,N)*GCMPS_WPM2*(T0B(M,N)-TT(M,N,1))
C              END IF
C              IF(IGRZS(M,N).EQ.2) THEN
C FLUXES FOR SALINITY:
C DEFINITION NET FRESH WATER FLUX FOR MM/DAY
c         XXT(M,N,10)= CMPS_MMPDAY*QSAL(M,N)
C DEFINITION NET RIVER RUNOFF, MM/DAY
c         XXT(M,N,11)= CMPS_MMPDAY*RUNOFF(M,N)
C DEFINITION NET FRESH WATER FLUX FOR MM/DAY
c         XXT(M,N,12)=-CMPS_MMPDAY*S0(M,N)
c     &                /max(SS(M,N,1)+35.0,0.01)
C DEFINITION RELAXING FRESH WATER FLUX FOR MM/DAY
c         XXT(M,N,13)= DKFS(M,N)*CMPS_MMPDAY*(S0B(M,N)-SS(M,N,1))
c     &                                  /max(SS(M,N,1)+35.0,0.01)
C              ENDIF
c         XXT(M,N,14)=SICE(M,N)

c         XXT(M,N,15)=0.001*SENS_HEAT(ROA,CPA,CDHW,E0+WIND(M,N)
c     &                            ,TA(M,N),TT(M,N,1))
c	   XXT(M,N,16)=0.001*LAT_HEAT(ROA,QLW,CDHW,E0+WIND(M,N),QA(M,N),
c     &               SATUR_HUM(SATUR_PRESS_OCEAN(TT(M,N,1)),SLPR(M,N)))

c         XXT(M,N,17)=0.001*(EW*(LW(M,N)-SIGMA*(TT(M,N,1)+273.15)**4))
c         ENDIF
c         ENDDO
c      ENDDO
!$OMP END PARALLEL DO

c      FORCING_TITLE=
c     &'1-DWLWrad,2-DWSWrad,3-Tatm,4-Qatm,5-Patm,6-NetHFlux,7-SW-Balance,
c     &8-ForcHeat,9-RlxHF,10-PmE,11-Runoff,12-NetWF,13-RlxWF,14-IceConcen
c     &tr,15-SensHeat,16-LatHeat,17-LW-Balance'
c      CALL WDSTD(PATH2OCP,'sfl.dat', 1,XXT,LU,NX,NY,17,
c     &                          MMM,MM,NNN,NN,1,17,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'sfl.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,17,1,
c     &            XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,ZNZP1,
c     &                    FORCING_TITLE,'sfl')


c      CALL WDSTD(PATH2OCP,'wws.dat', 1,WW,  LU ,NX,NY,NZ+1,
c     &                       MMM,MM,NNN,NN,1,NZ+1,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'wws.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,NZ+1,1,
c     &            XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &                    ZNZP1,'Vertical velocity[cm/s]','w')

c      CALL WDSTD(PATH2OCP,'dens.dat',1,DEN, LU ,NX,NY,NZ,
c     &                       MMM,MM,NNN,NN,1,NZ,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'dens.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,NZ,1,
c     &            XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &                    ZNZP1,'Density[g/cm**3]','r')

c      CALL WDSTD(PATH2OCP,'txo.dat', 1,TAUX*RH0, LU,NX,NY,1,
c     &                       MMM,MM,NNN,NN,1,1,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'txo.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,1,1,
c     &            XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &      ZNZP1,'Meridional surface wind stress[din/cm**2]','tx')

c      CALL WDSTD(PATH2OCP,'tyo.dat', 1,TAUY*RH0, LU,NX,NY,1,
c     &                       MMM,MM,NNN,NN,1,1,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'tyo.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,1,1,
c     &            XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &      ZNZP1,'Zonal surface wind stress[din/cm**2]','ty')

c      CALL WDSTD(PATH2OCP,'uwnd.dat', 1,UWND*0.01, LU,NX,NY,1,
c     &                       MMM,MM,NNN,NN,1,1,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'uwnd.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,1,1,
c     &            XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &      ZNZP1,'Zonal surface wind speed[m/s]','u')

c      CALL WDSTD(PATH2OCP,'vwnd.dat', 1,VWND*0.01, LU,NX,NY,1,
c     &                       MMM,MM,NNN,NN,1,1,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'vwnd.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,1,1,
c     &            XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &      ZNZP1,'Meridional surface wind speed[m/s]','v')


c       IF(KSW_TSV.GE.2.OR.KSW_UVV.GE.2) THEN
C STORE RICHARDSON NUMBER AS ITS ln
c       DO K=1,NZ
c          DO N=NNN,NN
c           DO M=MMM,MM
c              IF(LU(M,N).GT.0.5) THEN
c                 RIT(M,N,K)=ALOG(RIT(M,N,K)+1.0E-8)
c              END IF
c           END DO
c          END DO
c       END DO
c      CALL WDSTD(PATH2OCP,'ri.dat', 1,RIT, LU ,NX,NY,NZ+1,
c     &                       MMM,MM,NNN,NN,1,NZ+1,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'ri.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,NZ+1,1,
c     &             XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &                    ZNZP1,'Richardson Number','Ri')

c      CALL WDSTD(PATH2OCP,'nut.dat',1,ANZT,LU ,NX,NY,NZ+1,
c     &                       MMM,MM,NNN,NN,1,NZ+1,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'nut.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,NZ+1,1,
c     &             XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &   ZNZP1,'Vertical diffusion coefficient[cm**2/c]','nut')

c      CALL WDSTD(PATH2OCP,'nuu.dat',1,ANZU,LU,NX,NY,NZ+1,
c     &                       MMM,MM,NNN,NN,1,NZ+1,IERR)
c      CALL FULFNAME(FNAME,PATH2OCP,'nuu.dat',IERR)
c      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,NZ+1,1,
c     &             XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &   ZNZP1,'Vertical viscosity coefficient[cm**2/c]','nut')
c       END IF



C HORIZONTALLY AVERAGE T, S AND DENSITY PROFILES.
c      DO K = 1,NZ
c       TTA(K) = FAVER (TT(1,1,K))
c       SSA(K) = FAVER (SS(1,1,K))
c       DDA(K) =(FAVER(DEN(1,1,K))+RH0-1.0)*1.0E+03
c      ENDDO
c      WRITE(*,*)'  MEAN VERTICAL PROFILES:'
c      WRITE(*,*)'  Tem     Sal     Den'
c      DO K=1,NZ
c      WRITE(*,'(1X,3F8.3)') TTA(K),SSA(K),DDA(K)
c      ENDDO
c      CALL FZONAL(HHQ,FFZON)
c      WRITE(*,*)'   HHQ ZONAL AVERAGE '
c      WRITE(*,'(1X,16I5)')(NINT(FFZON(N)*1.0E-2),N=NNN,NN,NDISPC*6)
c	WRITE(*,*) (NINT(FFZON(N)*1.0E-2),N=NNN,NN,NDISPC*6)

      RETURN
125   WRITE(*,*) 'ERROR IN OPENING FILE cpmmc.txt !!!'
      STOP

126   WRITE(*,*) 'ERROR IN WRITING IN FILE cpmmc.txt !!!'
      STOP
      END
