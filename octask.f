C 26.09.09 14.04
C MAIN MODULE FOR OGCM (C) GUSEV A.V, DIANSKY N.A., STARTED AT 18.10.99 20:59

      PROGRAM MAIN
      use mpi_parallel_tools
      IMPLICIT NONE

C+++++++++++++++++PETSc+++++++++++++++++++++++++++++++++++++++++++++++++
c#include <petsc/finclude/petscsys.h>
c#include <petsc/finclude/petscvec.h>
c#include <petsc/finclude/petscmat.h>
c#include <petsc/finclude/petscpc.h>
c#include <petsc/finclude/petscksp.h>
c#include <petsc/finclude/petscviewer.h>

      REAL RELAXT,RELAXS       !VALUE OF RELAXING COEFFICIENT FOR T&S[CM/S]
C+++++++++++++++++TASK CONTROL PARAMETERS+++++++++++++++++++++++++++++++
      INCLUDE '2TASKCNTRL.INC'
      INCLUDE '1BASINPAR.INC'
c      include 'mpif.h'

      real*8 :: w_time_global, w_time

      integer :: locn
C+++++++++++++++++TIME CONTROL PARAMETERS+++++++++++++++++++++++++++++++
      REAL
     & WRXYZE,WRXYZI, !EXTERNAL&INTERNAL PERIODS IN DAYS TO WRITE 3D ARRAYS
     & WRPTE, WRPTI,  !EXTERNAL&INTERNAL PERIODS IN DAYS TO WRITE 3D ARRAYS
     & WRXYE, WRXYI,  !EXTERNAL&INTERNAL PERIODS IN DAYS TO WRITE 2D ARRAYS
     & WRYZE, WRYZI   !EXTERNAL&INTERNAL PERIODS IN DAYS TO WRITE 2D ARRAYS

      INTEGER IST,  !TIME STEP COUNTER (CALCULATED)
     &        IST0  !INITIAL TIME STEP (DOES NOT CHANGE)

      INTEGER INITYEAR, !INITIAL YEAR OF  TIME INTEGRATION
     &        LWRINT,   !PERIOD IN STEPS TO WRITE INTEGRAL VALUES (INPUT)
     &        LWRCP,    !PERIOD IN STEPS TO WRITE  TO WRITE CTRL POINT
     &        LWRLOC,   !PERIOD IN STEPS TO WRITE  TO WRITE LOCAL INFORMATION
     &        NSLOR,    !NUMBER OF ITERATIONS (CALCULATED)
     &        NOFCOM,   !NUMBER OF LINES IN "oct.par" (CALCULATED)
     & LWRXYZE,LWRXYZI, !EXTERNAL&INTERNAL PERIODS IN STEPS TO WRITE 3D ARRAYS
     & LWRPTE, LWRPTI,  !EXTERNAL&INTERNAL PERIODS IN STEPS TO WRITE 3D PASSIVE TRACER ARRAYS
     & LWRXYE, LWRXYI,  !EXTERNAL&INTERNAL PERIODS IN STEPS TO WRITE 2D ARRAYS
     & LWRYZE, LWRYZI,  !EXTERNAL&INTERNAL PERIODS IN STEPS TO WRITE 2D ARRAYS
     & NRECXYZ,     !No.OF RECORD TO WRITE XYZ ARARAYS
     & NRECPT,      !No.OF RECORD TO WRITE PASSIVE TRACER XYZ ARARAYS
     & NRECXY,      !No.OF RECORD TO WRITE XY  ARARAYS
     & NRECYZ       !No.OF RECORD TO WRITE YZ  ARARAYS

C        CONDITIONS FOR MONTHLY OUTPUT
C        TRUE - OUTPUT IN EVERY CALENDAR MONTH
C       FALSE - OUTPUT WITH GIVEN TIME INTERVAL

       LOGICAL MO_OUTPUT_XYZ,  ! FOR MODEL VARIABLES XYZ-ARRAYS
     &         MO_OUTPUT_PT,   ! FOR PASSIVE TRACER  XYZ-ARRAYS
     &         MO_OUTPUT_XY,   ! FOR MODEL VARIABLES XY -ARRAYS
     &         MO_OUTPUT_YZ,   ! FOR MODEL VARIABLES YZ -ARRAYS
     &         MO_LWRINT       ! FOR INTEGRAL PARAMETERS

C+++++++++++++++++SEA SURFACE CONTROL PARAMETERS++++++++++++++++++++++++
      INTEGER IGT,IGS,IGWS,IGICE,IAVT,IAVS,NICESTEP   !,ICESTEP
C         TYPES OF SS CONDITION FOR T,S:
C   IGT,IGS=
C         1 - FIRST CONDITION
C         (ONLY T,S AT OCEAN SURFACE ARE USED)
C         2 - SECOND CONDITION
C         (HEAT AND FRESH WATER BALANCES ARE USED)
C         3 - THIRD CONDITION
C         (HEAT AND SALT FLUXES ARE CALCULATED BY BULK-FORMULAE)
C
C         TYPES OF SS CONDITION FOR U,V:
C   IGWS=
C         1 - WIND STRESS COMPONENTS ARE READ FROM FILES,
C             ATM PRESSURE DOES NOT AFFECT DYNAMICS
C         2 - WIND STRESS COMPONENTS ARE READ FROM FILES,
C             ATM PRESSURE AFFECTS DYNAMICS
C         3 - WIND STRESS COMPONENTS ARE CALCULATED BY BULK FORMULAE,
C             ATM PRESSURE DOES NOT AFFECT DYNAMICS
C         4 - WIND STRESS COMPONENTS ARE CALCULATED BY BULK FORMULAE,
C             ATM PRESSURE AFFECTS DYNAMICS
C
C         USAGE OF ICE BLOCK IN OCEAN MODEL:
C   IGICE=
C         0 - ICE BLOCK IS NOT USED
C        >0 - ICE BLOCK IS     USED
C
C IAVT=  1-DO NOT REMOVE SPACE AVERAGING HEAT FLUX,
C        2-       REMOVE SPACE AVERAGING HEAT FLUX; (NOT USED)
C IAVS=  1-DO NOT REMOVE SPACE AVERAGING WRESH WATER FLUX,
C        2-       REMOVE SPACE AVERAGING WRESH WATER FLUX.
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C INPUT PARAMETERS OF TASK AND NAMES OF FILES WITH INPUT DATA
      CHARACTER(1024) COMMENTS(1024)
      CHARACTER(1024) FILEPAR,             !PARAMETER FOLE NAME
     &               PATH2OCP,             !PATH TO CONTROL POINTS(RESULTS)
     &               PATH2OCEANDATA,       !PATH TO SURFACE DATA ON OCEAN GRID
     &               PATH2ATMFORCING,      !PATH TO SURFACE DATA ON ATM   GRID
     &               BLANK
      CHARACTER(2048) FNAME
      CHARACTER(1024) SS_FILES(20),          !FILES WITH SEA SURFACE DATA
     &                ATMASK,                !FILE WITH ATMOSPHERIC SEA-LAND MASK (1-LAND,0-OCEAN)
     &                SS_FILES_FULLNAME(20)  !FILES WITH SEA SURFACE DATA (FULLNAMES)

      CHARACTER(1024) FLQBW(2) !FILES WITH LIQUID BOUNDARY VALUES FOR TEM & SAL

C----------------------------------------------------------------------
C PARAMETER FOR DATA INTERPOLATION FROM MONTH STEP ON THE MODEL STEP:
C WEIGHT MONTH COEFFICIENTS FOR IN-MONTH INTERPOLATION:
      REAL TAU,     !MODEL TIME STEP IN SECONDS
     &     TAUH,    !MODEL TIME STEP IN HOURS
     &     TAUD,    !MODEL TIME STEP IN DAYS
     &     DAYF

      INTEGER MONTH_OF_4YR,    !THE PRESENT MONTH IN 4-YEAR
     &        NEWMONTH,        !COMMAND TO READING DATA(>1-Y,0-NO)
     &        MONTH,           !THE PRESENT MONTH
     &        MONTH_OF_YEAR,   !MONTH OF YEAR
     &        NSTEP_PER_DAY,   !number of step per day
     &        INIT_DAY_OF_YEAR_FOR_LOCAL_OUTPUT

C----------------LOCAL COUNTERS----------------------------------------
      INTEGER M, K, IERR, MONTHWR, KALCSTMON, ISTREC, NHELP
C----------------------------------------------------------------------

      REAL TIME_COUNTER,       !TIME COUNTER IN HOURS
     &     DAY_COUNTER,        !DAY COUNTER
     &     DAY_OF_4YR,         !THE PRESENT DAY OF 4-YEAR PERIODICITY
     &     DAY_OF_MONTH,       !THE PRESENT DAY OF THE PRESENT MONTH
     &     TIME_START,         !START TIME FOR OUTPUT FILES
     &     TIME_INCREMENT      !INCREMENT TIME FOR OUTPUT FILES


      REAL SECNDS_OF_DAY        !CURRENT SECONDS IN DAY

      INTEGER NDAYS_IN_4YR(0:48),!INTEGER DAY DISTRIBUTIONS IN 4-YEARS
     &        M_DAY,             !MODEL ELAPSED DAY COUNTER STARTING FROM ZERO
     &        M_SEC_OF_MIN,      !SECOND COUNTER IN MINUTE
     &        M_MIN_OF_HOUR,     !MINUTE COUNTER IN HOUR
     &        M_HOUR_OF_DAY,     !HOUR COUNTER IN DAY
     &        M_DAY_OF_MONTH,    !DAY COUNTER IN MONTH
     &        M_DAY_OF_YEAR,     !DAY COUNTER IN YEAR
     &        M_DAY_OF_4YR,      !DAY COUNTER IN 4-YEARS
     &        M_MON_OF_YEAR,     !MON COUNTER IN YEAR
     &        M_MON_OF_4YR,      !MON COUNTER IN 4-YEARS
     &        M_YEAR_OF_4YR,     !YEAR COUNTER IN 4YRS
     &        M_YEAR,            !YEAR COUNTER
     &        M_4YR,             !COUNTER OF 4-YR GROUPS
     &        M_TIME_CHANGED(7), !CHANGE INDICATOR OF TIME
     &        KEY_TIME_PRINT     !KEY OF PRINTING TIME:0-NOT,1-PRINT

      INCLUDE '1DAYDIST.INC'  !DAY DISTRIBUTION IN 4-YEAR


C Init MPI
	  call MPI_INIT(ierr)
      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

      period = (/1,1/)
      p_size = (/0,0/)
      ierr = 0
!      CART_COMM = 0

      call MPI_COMM_RANK(PETSC_COMM_WORLD,rank,ierr)
      call MPI_COMM_SIZE(PETSC_COMM_WORLD,procs,ierr)
      call MPI_DIMS_CREATE(procs,2,p_size,IERR)
      call MPI_CART_CREATE(PETSC_COMM_WORLD,2,p_size,period,
     &                     0,CART_COMM,ierr)
      call MPI_CART_COORDS(CART_COMM,rank,2,p_coord,ierr)

c-----------------------------------NX-------------------------------------------
      locn = floor(real(NX)/real(p_size(1)))
      nx_start = locn*p_coord(1) + 1
      if ( p_coord(1) .EQ. p_size(1)-1 ) then
          locn = NX - nx_start + 1
      endif
      nx_end = nx_start + locn - 1
      nx_start = nx_start
c     border area
      bnd_x1 = nx_start - 1
      if (bnd_x1 < 1) bnd_x1 = 1
      bnd_x2 = nx_end + 1
      if (bnd_x2 > NX) bnd_x2 = NX

c-----------------------------------NY-------------------------------------------
      locn = floor(real(NY)/real(p_size(2)))
      ny_start = locn*p_coord(2) + 1
      if ( p_coord(2) .EQ. p_size(2)-1 ) then
          locn = NY - ny_start + 1
      endif
      ny_end = ny_start + locn - 1
      ny_start = ny_start
c     border area
      bnd_y1 = ny_start - 1
      if (bnd_y1 < 1) bnd_y1 = 1
      bnd_y2 = ny_end + 1
      if (bnd_y2 > NY) bnd_y2 = NY

!      print *, nx_start, nx_end, ny_start, ny_end
!      print *, p_coord, p_size, nx_start, nx_end, ny_start, ny_end


C ZERO SETTING OF:
      BLANK=' '
      CALL ATM_AR_2_ZERO()

      M_DAY=0          !MODEL ELAPSED DAY COUNTER STARTING FROM ZERO
      SECNDS_OF_DAY =0 !CURRENT SECONDS IN DAY  ,OUTPUT
      M_SEC_OF_MIN  =0 !SECOND COUNTER IN MINUTE,OUTPUT
      M_MIN_OF_HOUR =0 !MINUTE COUNTER IN HOUR  ,OUTPUT
      M_HOUR_OF_DAY =0 !HOUR COUNTER IN DAY     ,OUTPUT
      M_DAY_OF_MONTH=0 !DAY COUNTER IN MONTH    ,OUTPUT
      M_DAY_OF_YEAR =0 !DAY COUNTER IN YEAR     ,OUTPUT
      M_DAY_OF_4YR  =0 !DAY COUNTER IN 4-YEARS  ,OUTPUT
      M_MON_OF_YEAR =0 !MON COUNTER IN YEAR     ,OUTPUT
      M_MON_OF_4YR  =0 !MON COUNTER IN 4-YEARS  ,OUTPUT
      M_YEAR_OF_4YR =0 !YEAR COUNTER IN 4YRS    ,OUTPUT
      M_YEAR        =0 !YEAR COUNTER            ,OUTPUT
      M_4YR         =0 !COUNTER OF 4-YR GROUPS  ,OUTPUT
      M_TIME_CHANGED=0 !CHANGE INDICATOR OF TIME,OUTPUT
      KEY_TIME_PRINT=0 !KEY OF PRINTING TIME:0-NOT,1-PRINT

      DO K=1,49
       NDAYS_IN_4YR(K-1)=INT(DAYS_IN_4YR(K))
	END DO

C-----------------------------------------------------------------------
      FILEPAR='octask.par'
C     READING PARAMETERS FROM FILE

      CALL READPAR(FILEPAR,COMMENTS,NOFCOM)
      READ(COMMENTS( 1),*) TAU      !TIME STEP IN SECONDS
      READ(COMMENTS( 2),*) DAYF     !DURATION OF RUN IN DAYS
      READ(COMMENTS( 3),*) IST      !TIME STEP
      READ(COMMENTS( 4),*) INITYEAR !INITIAL YEAR OF TIME INTEGRATION
      READ(COMMENTS( 5),*) LWRCP  !PERIOD IN STEPS TO WRITE CP
      READ(COMMENTS( 6),*) LWRINT !PERIOD IN STEPS TO WRITE INTEGRALS
      READ(COMMENTS( 7),*) LWRLOC !PERIOD IN STEPS TO WRITE LOCAL INFORMATION
      READ(COMMENTS( 8),*) WRXYZE,WRXYZI !EXT&INT PERIODS IN DAYS
      READ(COMMENTS( 9),*) WRPTE, WRPTI  !EXT&INT PERIODS IN DAYS
      READ(COMMENTS(10),*) WRXYE, WRXYI  !EXT&INT PERIODS IN DAYS
      READ(COMMENTS(11),*) WRYZE, WRYZI  !EXT&INT PERIODS IN DAYS
      READ(COMMENTS(12),*) IGT,IGS,IGWS   !TEM,SAL,WINDSTRES TYPES
      READ(COMMENTS(13),*) IGICE,NICESTEP !ICE USAGE(1-ON,2-OFF)/NUM OF ICE STEPS ON 1 OCSTEP
      READ(COMMENTS(14),*) IAVT,IAVS      !TEM.&SAL. SPACE AVERAGING
      READ(COMMENTS(15),*) RELAXT,RELAXS  !RELAXATION COEFFICIENTS[CM/S]
      CALL GET_FIRST_LEXEME(COMMENTS(16), PATH2OCP    )    !PATH TO DIR. WITH CP
      CALL GET_FIRST_LEXEME(COMMENTS(17), PATH2OCEANDATA)  !PATH TO OCEAN DATA
      CALL GET_FIRST_LEXEME(COMMENTS(18), SS_FILES(1)   )  !FILE WITH SST
      CALL GET_FIRST_LEXEME(COMMENTS(19), SS_FILES(2)   )  !FILE WITH SSS
      CALL GET_FIRST_LEXEME(COMMENTS(20), FLQBW(1)      )  !FILE WITH LIQ.WAL OF T
      CALL GET_FIRST_LEXEME(COMMENTS(21), FLQBW(2)      )  !FILE WITH LIQ.WAL OF S
      CALL GET_FIRST_LEXEME(COMMENTS(22), PATH2ATMFORCING )!PATH TO DIR. WITH FLUXES
      CALL GET_FIRST_LEXEME(COMMENTS(23), SS_FILES(3)   )  !FILE WITH TAUX
      CALL GET_FIRST_LEXEME(COMMENTS(24), SS_FILES(4)   )  !FILE WITH TAUY
      CALL GET_FIRST_LEXEME(COMMENTS(25), SS_FILES(5)   )  !FILE WITH HEAT BALANCE
      CALL GET_FIRST_LEXEME(COMMENTS(26), SS_FILES(6)   )  !FILE WITH SW-RAD BALANCE
      CALL GET_FIRST_LEXEME(COMMENTS(27), SS_FILES(7)   )  !FILE WITH PRECIP-EVAP
      CALL GET_FIRST_LEXEME(COMMENTS(28), SS_FILES(8)   )  !FILE WITH ICE MASK
      CALL GET_FIRST_LEXEME(COMMENTS(29), SS_FILES(9)   )  !FILE WITH SST ON ATMGRID
	  CALL GET_FIRST_LEXEME(COMMENTS(30), SS_FILES(10)  )  !FILE WITH PRESSURE
      CALL GET_FIRST_LEXEME(COMMENTS(31), SS_FILES(11)  )  !FILE WITH RIVER RUNOFF
      CALL GET_FIRST_LEXEME(COMMENTS(32), SS_FILES(12)  )  !FILE WITH DW-LW-RAD
      CALL GET_FIRST_LEXEME(COMMENTS(33), SS_FILES(13)  )  !FILE WITH DW-SW-RAD
      CALL GET_FIRST_LEXEME(COMMENTS(34), SS_FILES(14)  )  !FILE WITH PRECIPITATION
      CALL GET_FIRST_LEXEME(COMMENTS(35), SS_FILES(15)  )  !FILE WITH TEMP OF SAT
      CALL GET_FIRST_LEXEME(COMMENTS(36), SS_FILES(16)  )  !FILE WITH HUMIDITY
      CALL GET_FIRST_LEXEME(COMMENTS(37), SS_FILES(17)  )  !FILE WITH U-WIND SPEED
      CALL GET_FIRST_LEXEME(COMMENTS(38), SS_FILES(18)  )  !FILE WITH V-WIND SPEED
      CALL GET_FIRST_LEXEME(COMMENTS(39), ATMASK)          !FILE WITH ATM MASK

	  INITYEAR=MAX(1,INITYEAR)  !THERE IS NO ZERO YEAR (SEE CALENDAR PLEASE)

      IST0=IST
      TAUH = TAU / 3600.00
      TAUD = TAU /86400.00
      NSTEP_PER_DAY=NINT(86400.0/TAU)       !number of step per day


C SET THE OUTPUT CONTROLS IN TIME
      MO_OUTPUT_XYZ=WRXYZI.GE.30.0
      MO_OUTPUT_PT =WRPTI .GE.30.0
      MO_OUTPUT_XY =WRXYI .GE.30.0
      MO_OUTPUT_YZ =WRYZI .GE.30.0
      MO_LWRINT    =LWRINT*TAUD.GE.30.0

      IF(MO_LWRINT) THEN
      LWRINT=NINT(FLOAT(LWRINT)*TAUH/(24.0*30.5))
      END IF

      IF(MO_OUTPUT_XYZ) THEN
      LWRXYZE=NINT(WRXYZE/30.5)      !EXT. PERIOD IN MONTHS TO WRITE 3D ARRAYS
      LWRXYZI=NINT(WRXYZI/30.5)      !INT. PERIOD IN MONTHS TO WRITE 3D ARRAYS
      ELSE
      LWRXYZE=NINT(WRXYZE*24.0/TAUH) !EXT. PERIOD IN STEPS  TO WRITE 3D ARRAYS
      LWRXYZI=NINT(WRXYZI*24.0/TAUH) !INT. PERIOD IN STEPS  TO WRITE 3D ARRAYS
      END IF

      IF(MO_OUTPUT_PT) THEN
      LWRPTE =NINT(WRPTE/30.5)      !EXT. PERIOD IN MONTHS TO WRITE 3D ARRAYS
      LWRPTI =NINT(WRPTI/30.5)      !INT. PERIOD IN MONTHS TO WRITE 3D ARRAYS
      ELSE
      LWRPTE =NINT(WRPTE*24.0/TAUH) !EXT. PERIOD IN STEPS  TO WRITE 3D ARRAYS
      LWRPTI =NINT(WRPTI*24.0/TAUH) !INT. PERIOD IN STEPS  TO WRITE 3D ARRAYS
      END IF

      IF(MO_OUTPUT_XY) THEN
      LWRXYE =NINT(WRXYE /30.5)      !EXT. PERIOD IN MONTHS TO WRITE 2D ARRAYS
      LWRXYI =NINT(WRXYI /30.5)      !INT. PERIOD IN MONTHS TO WRITE 2D ARRAYS
      ELSE
      LWRXYE =NINT(WRXYE *24.0/TAUH) !EXT. PERIOD IN STEPS  TO WRITE 2D ARRAYS
      LWRXYI =NINT(WRXYI *24.0/TAUH) !INT. PERIOD IN STEPS  TO WRITE 2D ARRAYS
      END IF

      IF(MO_OUTPUT_YZ) THEN
      LWRYZE =NINT(WRYZE /30.5)      !EXT. PERIOD IN MONTHS TO WRITE 2D ARRAYS
      LWRYZI =NINT(WRYZI /30.5)      !INT. PERIOD IN MONTHS TO WRITE 2D ARRAYS
      ELSE
      LWRYZE =NINT(WRYZE *24.0/TAUH) !EXT. PERIOD IN STEPS  TO WRITE 2D ARRAYS
      LWRYZI =NINT(WRYZI *24.0/TAUH) !INT. PERIOD IN STEPS  TO WRITE 2D ARRAYS
      END IF

      LWRXYZE= MAX(LWRXYZE,LWRXYZI)  !CORRECT EXTERNAL PERIOD IF NEED
      LWRPTE = MAX(LWRPTE ,LWRPTI )  !CORRECT EXTERNAL PERIOD IF NEED
      LWRXYE = MAX(LWRXYE ,LWRXYI )  !CORRECT EXTERNAL PERIOD IF NEED
      LWRYZE = MAX(LWRYZE ,LWRYZI )  !CORRECT EXTERNAL PERIOD IF NEED

C OCEAN PARAMETERS SETTING
      CALL OCPAR(IGT,IGS)

      KALCSTMON=0


C-----------------------------------------------------------------------
      if (rank .eq. 0) then
      WRITE(*,'(2X,5HSTEP:,F7.2,4HHRS;,F9.2,4HSEC;,F9.5,4HDAY.)')
     &                     TAUH,      TAU,       TAUD
      WRITE(*,'(2X,15HDURATION OF RUN:,F9.2,6H DAYS.)') DAYF
      endif

C-----------------GETTING TIME FOR TIME CYCLE---------------------------
      CALL TIMESERV(0)

      NEWMONTH     =2
      TIME_COUNTER =FLOAT(IST)*TAUH
      DAY_COUNTER  =FLOAT(IST)*TAUH/24.0
      DAY_OF_4YR   =MOD(DAY_COUNTER,DAYS_IN_4YR(49))

      MONTH_OF_4YR=2
      DO WHILE (DAY_OF_4YR.GE.DAYS_IN_4YR(MONTH_OF_4YR) )
         MONTH_OF_4YR = MONTH_OF_4YR+1
      END DO
      MONTH_OF_4YR = MONTH_OF_4YR-1

      MONTH_OF_YEAR = MOD(MONTH_OF_4YR,12)
      IF (MONTH_OF_YEAR.EQ.0) THEN
         MONTH_OF_YEAR = 12
      END IF
      MONTHWR = MONTH_OF_YEAR    !FOR WRITING MONTHLY MEAN
C            MONTHWR = MONTH  !FOR WRITING MONTHLY MEAN VALUE FOR PRESENT MON

      MONTH=48*INT(DAY_COUNTER/DAYS_IN_4YR(49))+MONTH_OF_4YR

C-------------------------------------------------------------------------------
C START OF TIME INTEGRATION OF OCEAN MODEL
C-------------------------------------------------------------------------------
      if (rank .eq. 0) then
          WRITE(*,'(A,1I7)') ' '
          WRITE(*,'(A,1I7)') 'START OF TIME INTEGRATION OF OCEAN MODEL'
          WRITE(*,'(A,1I7)') ' '
      endif

      call start_timer(w_time_global)
      DO WHILE(TIME_COUNTER.LT.DAYF*24.0)

        IST = IST + 1
        TIME_COUNTER= TIME_COUNTER + TAUH

        CALL MODEL_TIME_DEF(
     &        IST,            !STEP COUNTER,            INPUT
     &        TAU,            !TIME STEP IN SECONDS,    INPUT
     &        NDAYS_IN_4YR,   !INTEGER DAY DISTRIBUTION IN 4-YEARS (49 MONTHS)
     &        M_DAY,          !MODEL ELAPSED DAY COUNTER STARTING FROM ZERO
     &        SECNDS_OF_DAY,  !CURRENT SECONDS IN DAY  ,OUTPUT
     &        M_SEC_OF_MIN,   !SECOND COUNTER IN MINUTE,OUTPUT
     &        M_MIN_OF_HOUR,  !MINUTE COUNTER IN HOUR  ,OUTPUT
     &        M_HOUR_OF_DAY,  !HOUR COUNTER IN DAY     ,OUTPUT
     &        M_DAY_OF_MONTH, !DAY COUNTER IN MONTH    ,OUTPUT
     &        M_DAY_OF_YEAR,  !DAY COUNTER IN YEAR     ,OUTPUT
     &        M_DAY_OF_4YR,   !DAY COUNTER IN 4-YEARS  ,OUTPUT
     &        M_MON_OF_YEAR,  !MON COUNTER IN YEAR     ,OUTPUT
     &        M_MON_OF_4YR,   !MON COUNTER IN 4-YEARS  ,OUTPUT
     &        M_YEAR_OF_4YR,  !YEAR COUNTER IN 4YRS    ,OUTPUT
     &        M_YEAR,         !YEAR COUNTER            ,OUTPUT
     &        M_4YR,          !COUNTER OF 4-YR GROUPS  ,OUTPUT
     &        M_TIME_CHANGED, !CHANGE INDICATOR OF TIME,OUTPUT
     &        KEY_TIME_PRINT, !KEY OF PRINTING TIME:0-NOT,1-PRINT
     &        INITYEAR)       !INITIAL REAL-TIME YEAR

        DAY_COUNTER = TIME_COUNTER/24.0
        DAY_OF_4YR  = (DAY_OF_4YR*24.0  + TAUH)/24.0

C TIME POSITIONING FOR CALCULATION OF
C WEIGT COEFFICIENTS FOR INTERPOLATION FROM MONTH STEP ON DAY STEP

        IF(DAY_OF_4YR.GE.DAYS_IN_4YR(MONTH_OF_4YR+1) ) THEN
C CHANGE PRESENT MONTH:
          NEWMONTH=NEWMONTH+1  !COMMAND TO READ DATA IF MONTH CHANGED
             MONTH_OF_4YR = MONTH_OF_4YR+1
             MONTH        = MONTH       +1
          IF (MONTH_OF_4YR.EQ.49) THEN
             DAY_OF_4YR = DAY_OF_4YR - DAYS_IN_4YR(49)
             MONTH_OF_4YR =   1
          ENDIF
             MONTH_OF_YEAR = MOD(MONTH_OF_4YR,12)
          IF (MONTH_OF_YEAR.EQ.0) THEN
             MONTH_OF_YEAR = 12
          END IF
             MONTHWR = MONTH_OF_YEAR-1  !FOR WRITING MONTHLY MEAN
          IF (MONTHWR.EQ.0) THEN         !VALUES FOR CLIMATIC YEAR
             MONTHWR = 12               !
          END IF                         !
C         MONTHWR = MONTH-1 !FOR WRITING MONTHLY MEAN VALUE FOR PRESENT MON
        ENDIF
        DAY_OF_MONTH= DAY_OF_4YR-DAYS_IN_4YR(MONTH_OF_4YR)

C===================OGCM INTEGRATION OVER ONE STEP========================
C------------------- SWALLOW WATER ----------------------------------------------
        CALL ADAPTATION_MODULE(TAU,NSLOR,IGT,IGS,IGWS)
C------------------- END SWALLOW WATER ------------------------------------------
C===================END OF OGCM INTEGRATION OVER ONE STEP========================

C=================== OUTPUT XY ==================================================
        IF(MO_OUTPUT_XY) THEN
C WRITE DATA FOR 15-TH DAY OF MONTH
           IF (15.0-0.5*TAUD.LE.DAY_OF_MONTH.AND.
     &          DAY_OF_MONTH.LE.15.0+0.5*TAUD.AND.
     &                    MOD(MONTH,LWRXYI).EQ.0) THEN
              NRECXY = MOD(MONTH,LWRXYE)/LWRXYI
              IF (NRECXY.EQ.0) NRECXY=(LWRXYE/LWRXYI)

              if (rank .eq. 0) then
              WRITE(*,'(A,1I7)')'  WRITE XY-DATA: REC= ',NRECXY
              WRITE(*,'(A,I5,A,F6.2,A,I5))')
     &                '  MONTH=',MONTH_OF_YEAR,
     &                '; DAY OF MONTH=',DAY_OF_MONTH,
     &                '; YEAR=',MONTH/12+1
              CALL MODEL_TIME_PRINT(M_HOUR_OF_DAY,
     &                            M_MIN_OF_HOUR,
     &                            M_SEC_OF_MIN,
     &                            M_DAY_OF_MONTH,
     &                            M_MON_OF_YEAR,
     &                            M_YEAR,
     &                            M_DAY_OF_YEAR,
     &                            M_DAY_OF_4YR)
              endif

              call start_timer(w_time)
              CALL XY_OUTPUT(PATH2OCP,NRECXY)
              call end_timer(w_time)
              if (rank.eq.0) print*,'XY_OUT:',w_time

           END IF
        ELSE
           IF (MOD(IST,LWRXYI).EQ.0) THEN
               NRECXY = MOD(IST,LWRXYE)/LWRXYI
               IF (NRECXY.EQ.0) NRECXY=(LWRXYE/LWRXYI)

               if (rank .eq. 0) then
               WRITE(*,'(A,1I7)')'  WRITE XY-DATA: REC= ',NRECXY
               CALL MODEL_TIME_PRINT(M_HOUR_OF_DAY,
     &                            M_MIN_OF_HOUR,
     &                            M_SEC_OF_MIN,
     &                            M_DAY_OF_MONTH,
     &                            M_MON_OF_YEAR,
     &                            M_YEAR,
     &                            M_DAY_OF_YEAR,
     &                            M_DAY_OF_4YR)
              endif

              call start_timer(w_time)
              CALL XY_OUTPUT(PATH2OCP,NRECXY)
              call end_timer(w_time)
              if (rank.eq.0) print*,'XY_OUT:',w_time

           END IF
        END IF
C=================== END OUTPUT XY ==============================================

        NEWMONTH =0 ! FOR NOT CHANGING MONTH
      END DO      ! END OF INTEGRATION

C-------------------------------------------------------------------------------
C END OF TIME INTEGRATION OF OCEAN MODEL
C-------------------------------------------------------------------------------
      call end_timer(w_time_global)
      if (rank .eq. 0) then
          WRITE(*,'(A,1I7)') ' '
          WRITE(*,'(A,1I7)') 'END OF TIME INTEGRATION OF OCEAN MODEL'
          WRITE(*,'(A,1I7)') ' '
          print *, "TIME GLOBAL:", w_time_global
      endif

c      call SIMPLE_OUT()

C WRITING CONTROL POINT FOR END OF RUN
c      CALL OCPWRITE(PATH2OCP,2,IGT,IGS,KALCSTMON)  !COMPLETE WRITE

C=================== CREATE CTL XY ==============================================
      IF(MO_OUTPUT_XY) THEN
         NRECXY=MIN(MONTH/LWRXYI,LWRXYE/LWRXYI)
         TIME_START    = MAX(0.5,FLOAT(MONTH-LWRXYE)+0.5)*30.0*86400.0
         TIME_INCREMENT= FLOAT(LWRXYI)*30.0*86400.0
      ELSE
         NRECXY=MIN(IST/LWRXYI,LWRXYE/LWRXYI)
         TIME_START    =MAX(0.0,FLOAT(IST-LWRXYE)*TAU)
         TIME_INCREMENT=FLOAT(LWRXYI)*TAU
      END IF

      if (rank .eq. 0) then
          CALL XY_CTL(PATH2OCP,NRECXY,TIME_START,TIME_INCREMENT)
      endif
C=================== END CREATE CTL XY ===========================================

C--------STOP TIME --------------------------------------
      CALL TIMESERV(1)

      CALL MPI_FINALIZE(IERR)

      STOP
      END


!      SUBROUTINE SIMPLE_OUT()
!      IMPLICIT NONE
!
!      INCLUDE '0COM.INC'
!      INCLUDE '0CEAN.INC'
!
!      INTEGER M, N
!
!      do M = MMM, MM, 10
!          do N = NNN, NN, 5
!              print *, SLH(M, N)
!          enddo
!      enddo
!      ENDSUBROUTINE
