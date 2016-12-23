C============TIME SERVICE FOR MODEL: GETTING START TIME====================
      SUBROUTINE TIMESERV(FLAGTIME)
	IMPLICIT NONE
      INTEGER FLAGTIME     !SWITCH FOR GET TIME:
C                          !IF FLAGTIME=0 => GET START TIME
C                          !IF FLAGTIME>0 => DEFINE ELAPSED TIME
      COMMON/ELAPSTIME/KDAY0,NHOUR0,NMINUTE0,NSEC0,NDECS0,  !START TIME
     &                 NDAY1,NHOUR1,NMINUTE1,NSEC1,NDECS1   !ELAPSED TIME
      INTEGER KDAY0,NHOUR0,NMINUTE0,NSEC0,NDECS0,    !START TIME
     &        NDAY1,NHOUR1,NMINUTE1,NSEC1,NDECS1     !ELAPSED TIME

      IF(FLAGTIME.EQ.0) THEN
C---------------- START TIME ------------------------------------------
       KDAY0=0
       NDAY1=1
c      CALL GETTIM(NHOUR0,NMINUTE0,NSEC0,NDECS0) !IF INSTRING PROCEDURE USED
       CALL GETIME(NHOUR0,NMINUTE0,NSEC0,NDECS0)
       WRITE(*,'(2X,A,I3,A,I2,A,I2,A,I2,A,I2)')'START DAY: ',
     &       NDAY1,'; TIME: ',NHOUR0,':',NMINUTE0,':',NSEC0,',',NDECS0
      RETURN
      END IF
      IF(FLAGTIME.NE.0) THEN
C----------DEFINING ELAPSED TIME --------------------------------------
c      CALL GETTIM(NHOUR1,NMINUTE1,NSEC1,NDECS1)
       CALL GETIME(NHOUR1,NMINUTE1,NSEC1,NDECS1)
       NDECS1  =NDECS1-NDECS0
       NSEC1   =NSEC1-   NSEC0
       NMINUTE1 =NMINUTE1- NMINUTE0
       NHOUR1  =NHOUR1-  NHOUR0
       IF(NDECS1.LT.0) THEN
         NDECS1=100+NDECS1
         NSEC1   =NSEC1-1
       END IF
       IF(NSEC1.LT.0) THEN
         NSEC1  =60+NSEC1
         NMINUTE1=NMINUTE1-1
       END IF
       IF(NMINUTE1.LT.0) THEN
         NMINUTE1=60+NMINUTE1
         NHOUR1 =NHOUR1-1
       END IF
       IF(NHOUR1.LT.0) THEN
         NHOUR1=24+NHOUR1
         KDAY0=KDAY0+1
       ELSE
         KDAY0=0
       END IF

       IF(KDAY0.EQ.1) THEN
          NDAY1=NDAY1+1
       END IF

       WRITE(*,'(2X,A,I3,A,I2,A,I2,A,I2,A,I2)')'ELAPSED DAY: ',
     &       NDAY1,'; TIME: ',NHOUR1,':',NMINUTE1,':',NSEC1,',',NDECS1
      RETURN
      END IF

      RETURN
      END
C======================================================================
      SUBROUTINE GETIME(NHOUR,NMINUTE,NSEC,NDECS)
	IMPLICIT NONE
      INTEGER NHOUR,NMINUTE,NSEC,NDECS   !NUMB.OF HRS,MIN,SEC,SEC101
      REAL TIME_IN_SEC,SEC

      CALL CPU_TIME(TIME_IN_SEC)

      NHOUR  = INT(TIME_IN_SEC/3600.0)

      NMINUTE = INT((TIME_IN_SEC-FLOAT(NHOUR)*3600.0)/60.0)

      SEC    = TIME_IN_SEC-FLOAT(NHOUR)*3600.0-FLOAT(NMINUTE)*60.0

      NSEC   = INT(SEC)
      NDECS  = INT(SEC*100.0)-NSEC*100

      RETURN
      END
C======================================================================
      SUBROUTINE MODEL_TIME_DEF(
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

      IMPLICIT NONE 
C+++++++++++++++++TIME CONTROL PARAMETERS+++++++++++++++++++++++++++++++
      REAL   TAU,                 !TIME STEP IN SECONDS
     &       SECNDS_OF_DAY        !CURRENT SECONDS IN DAY

      INTEGER IST             !TIME STEP COUNTER

      INTEGER NSTEP_PER_DAY,     !NUMBER OF STEP PER DAY
     &        NDAYS_IN_4YR(0:48),!INTEGER DAY DISTRIBUTIONS IN 4-YEARS
     &        M_DAY,             !MODEL ELAPSED DAY COUNTER STARTING FROM ZERO
     &        M_SEC_OF_MIN,      !SECOND COUNTER IN MINUTE
     &        M_MIN_OF_HOUR,     !MINUTE COUNTER IN HOUR
     &        M_HOUR_OF_DAY,     !HOUR COUNTER IN DAY
     &        M_DAY_OF_MONTH,    !DAY  COUNTER IN MONTH
     &        M_DAY_OF_YEAR,     !DAY  COUNTER IN YEAR
     &        M_DAY_OF_4YR,      !DAY  COUNTER IN 4-YEARS
     &        M_MON_OF_YEAR,     !MON  COUNTER IN YEAR
     &        M_MON_OF_4YR,      !MON COUNTER IN 4-YEARS
     &        M_YEAR_OF_4YR,     !YEAR COUNTER IN 4YRS
     &        M_YEAR,            !YEAR COUNTER 
     &        M_4YR,             !COUNTER OF 4-YR GROUPS
     &        M_TIME_CHANGED(7),!INDICATOR OF TIME CHANGED (0-NOT,1-CHANGED) FOR
                                 !1-SEC,2-MIN,3-HOUR,4-DAY,5-MONTH,6-YEAR,7-4YRS
     &        KEY_TIME_PRINT,    !KEY OF PRINTING TIME:0-NOT,1-PRINT     
     &        INITYEAR           !INITIAL REAL-TIME YEAR



C------------------ INTERNAL VARIABLES: --------------------------------
C HELP VARIABLE PREFFIX I DENOTES INSTANT
      INTEGER I_DAY_OF_YEAR,    !DAY COUNTER IN YEAR
     &        I_DAY_OF_4YR,     !DAY COUNTER IN 4-YEAR
     &        I_4YR,            !COUNTER OF 4-YR GROUPS
     &        I_SEC_OF_MIN,     !SECOND COUNTER IN MINUTE
     &        I_MIN_OF_HOUR,    !MINUTE COUNTER IN HOUR
     &        I_HOUR_OF_DAY,    !HOUR COUNTER IN DAY
     &        I_DAY_OF_MONTH,   !DAY COUNTER IN MONTH
     &        I_MON_OF_YEAR,    !MON COUNTER IN YEAR
     &        I_YEAR,           !YEAR COUNTER 
     &        I_MON_OF_4YR,     !MON COUNTER IN 4-YEARS
     &        I_YEAR_OF_4YR     !YEAR COUNTER IN 4YRS



      CHARACTER  MONTH_NAME(12)*3
      DATA MONTH_NAME/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG',
     &           'SEP','OCT','NOV','DEC'/
      INTEGER I,K
	INTEGER NSHIFT        ! NUMBER OF YEARS (0-3) FOR COUNTER SHIFT  
	                      ! IF THE REAL-TIME INITIAL YEAR DOES NOT FOLLOW
	                      ! THE LEAP YEAR
      INTEGER IHELP8, KHELP8
C-----------------------------------------------------------------------
      NSHIFT=MOD(INITYEAR-1,4)

      M_TIME_CHANGED=0
      
      NSTEP_PER_DAY=NINT(86400.0/TAU)
      
      IHELP8=INT(NSTEP_PER_DAY)
      
      M_DAY=INT((IST-1)/IHELP8)+1

      
      KHELP8=MOD(IST-1,IHELP8)+1
      SECNDS_OF_DAY = TAU*FLOAT(KHELP8)
      

      I_4YR=(M_DAY-1+NDAYS_IN_4YR(NSHIFT*12))/NDAYS_IN_4YR(48)+1

      I_DAY_OF_4YR=MOD(M_DAY-1+NDAYS_IN_4YR(NSHIFT*12),
     &                 NDAYS_IN_4YR(48))+1

          I_MON_OF_4YR=1
          I=49
1000      K=(I_MON_OF_4YR+I)/2
          IF (I_DAY_OF_4YR.LE.NDAYS_IN_4YR(K-1)) I=K
          IF (I_DAY_OF_4YR.GT.NDAYS_IN_4YR(K-1)) I_MON_OF_4YR=K
          IF (I.GT.I_MON_OF_4YR+1) GO TO 1000


      I_YEAR_OF_4YR = (I_MON_OF_4YR-1)/12 +1
      
      I_YEAR        = (I_4YR-1)*4 + I_YEAR_OF_4YR +INITYEAR-1-NSHIFT

      I_MON_OF_YEAR = MOD(I_MON_OF_4YR-1,12)+1
      
      I_DAY_OF_MONTH= I_DAY_OF_4YR-NDAYS_IN_4YR(I_MON_OF_4YR-1)
      
      I_DAY_OF_YEAR = I_DAY_OF_4YR-NDAYS_IN_4YR((I_YEAR_OF_4YR-1)*12)
                
      I_HOUR_OF_DAY = INT((SECNDS_OF_DAY-1.0)/3600.0)+1

      I_MIN_OF_HOUR = INT((SECNDS_OF_DAY-1.0)/60.0)+1   ! MIN IN DAY YET

      I_SEC_OF_MIN  = INT(SECNDS_OF_DAY) -  I_MIN_OF_HOUR * 60

      I_MIN_OF_HOUR = I_MIN_OF_HOUR      -  I_HOUR_OF_DAY * 60
                                        
      IF(M_SEC_OF_MIN  .NE.I_SEC_OF_MIN  ) THEN
         M_SEC_OF_MIN    = I_SEC_OF_MIN
         M_TIME_CHANGED(1)=1  
      END IF  

      IF(M_MIN_OF_HOUR .NE.I_MIN_OF_HOUR ) THEN
         M_MIN_OF_HOUR   = I_MIN_OF_HOUR
         M_TIME_CHANGED(2)=1  
      END IF 

      IF(M_HOUR_OF_DAY .NE.I_HOUR_OF_DAY ) THEN
         M_HOUR_OF_DAY   = I_HOUR_OF_DAY
         M_TIME_CHANGED(3)=1
      END IF  
       
      IF(M_DAY_OF_MONTH.NE.I_DAY_OF_MONTH) THEN
         M_DAY_OF_4YR    = I_DAY_OF_4YR
         M_DAY_OF_YEAR   = I_DAY_OF_YEAR
         M_DAY_OF_MONTH  = I_DAY_OF_MONTH
         M_TIME_CHANGED(4)=1  
      END IF
 
      IF(M_MON_OF_YEAR .NE.I_MON_OF_YEAR ) THEN
         M_MON_OF_YEAR   = I_MON_OF_YEAR
	   M_MON_OF_4YR    = I_MON_OF_4YR
         M_TIME_CHANGED(5)=1  
      END IF 
 
      IF(M_YEAR .NE.I_YEAR ) THEN
         M_YEAR        = I_YEAR
         M_YEAR_OF_4YR = I_YEAR_OF_4YR
         M_TIME_CHANGED(6)=1  
      END IF 
      
      IF(M_4YR .NE.I_4YR ) THEN
         M_4YR   = I_4YR
         M_TIME_CHANGED(7)=1  
      END IF 

      IF(KEY_TIME_PRINT.NE.0) THEN

      WRITE(*,'(4(A,I2.2), A,I4.4, A,I3.3, A,I4.4)')
     &    '   MODEL TIME: ',
     &                   M_HOUR_OF_DAY,':',
     &                   M_MIN_OF_HOUR,':',
     &                   M_SEC_OF_MIN, '  ',
     &                   M_DAY_OF_MONTH,
     &        MONTH_NAME(M_MON_OF_YEAR),M_YEAR,
     &        ';  DAY IN YEAR:',M_DAY_OF_YEAR,
     &        ',  DAY IN 4YRS:',M_DAY_OF_4YR      
      
      END IF 
      END SUBROUTINE MODEL_TIME_DEF
C======================================================================
      SUBROUTINE MODEL_TIME_PRINT(M_HOUR_OF_DAY,
     &                            M_MIN_OF_HOUR,
     &                            M_SEC_OF_MIN,
     &                            M_DAY_OF_MONTH,
     &                            M_MON_OF_YEAR,
     &                            M_YEAR,
     &                            M_DAY_OF_YEAR,
     &                            M_DAY_OF_4YR)
      IMPLICIT NONE


      INTEGER                     M_HOUR_OF_DAY,
     &                            M_MIN_OF_HOUR,
     &                            M_SEC_OF_MIN,
     &                            M_DAY_OF_MONTH,
     &                            M_MON_OF_YEAR,
     &                            M_YEAR,
     &                            M_DAY_OF_YEAR,
     &                            M_DAY_OF_4YR

      CHARACTER  MONTH_NAME(12)*3
      DATA MONTH_NAME/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG',
     &           'SEP','OCT','NOV','DEC'/  
     
     
         
      WRITE(*,'(4(A,I2.2), A,I4.4, A,I3.3, A,I4.4)')
     &    '   MODEL TIME: ',
     &                   M_HOUR_OF_DAY,':',
     &                   M_MIN_OF_HOUR,':',
     &                   M_SEC_OF_MIN, '  ',
     &                   M_DAY_OF_MONTH,
     &        MONTH_NAME(M_MON_OF_YEAR),M_YEAR,
     &        ';  DAY IN YEAR:',M_DAY_OF_YEAR,
     &        ',  DAY IN 4YRS:',M_DAY_OF_4YR       
      RETURN
	END
