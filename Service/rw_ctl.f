      SUBROUTINE CREATE_CTLFILE(FNAME,UNDEF,NX,NY,NZ,NT,
     &                  X0,HX,Y0,HY,T0,HT,Z,TITLE,VARNAME)
      use mpi_parallel_tools
	  use mod_shallowater
      IMPLICIT NONE 
      CHARACTER*(*) FNAME
      CHARACTER*(*) TITLE
      INTEGER       NX, NY, NZ, NT   !DIMENSION OF DATA
      REAL  X0(1),HX,Y0(1),HY, !INITIAL POINTS AND STEPS[DEG] OF HORIZONTAL GRID
                               !IF HX OR HY <0, THEN LEVELS ARE USED FOR APPROPRIATE VARIABLE
     &      T0(1),HT,          !INITIAL TIME AND STEP[SEC] IN TIME
     &      Z(NZ),             !VERTICAL LEVELS
     &      UNDEF              !UNDEFINIT VALUE
	INTEGER I
      CHARACTER*(*) VARNAME
      INTEGER NYRHT,NMOHT,NDYHT,NHRHT,NMNHT  !YEARS,MONTHS,DAYS,HOURS,MINUTES
      INTEGER NYRT0,NMOT0,NDYT0,NHRT0,NMNT0  !YEARS,MONTHS,DAYS,HOURS,MINUTES
      CHARACTER(128)  NAMECTL,NAMEDAT,NAMEDAT2
      CHARACTER  TFSTEP*2,MON*3,MONTH(12)*3
      DATA MONTH/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG',
     &           'SEP','OCT','NOV','DEC'/
C     CHARACTER(256) ZSTRING
      INTEGER K

      MON=MONTH(1)

!-------- create extensions ------------
      NAMEDAT = FNAME
      NAMECTL = FNAME
      CALL EXTNAME(NAMECTL,'.ctl')
      CALL EXTNAME(NAMEDAT,'.dat')

!--------  remove path from NAMEDAT ( .../.../xxxx.dat -> xxxx.dat )  ----
      K  = LEN(NAMEDAT)
C     DO WHILE(NAMEDAT(K:K).ne.'/'.and.NAMEDAT(K:K).ne.'\')
      DO WHILE(NAMEDAT(K:K).ne.'/')
          K = K - 1
          IF (K.LE.0) EXIT
      ENDDO
      WRITE(NAMEDAT2, '(A)') NAMEDAT( K+1: )

!-------- write to file     ------------
      OPEN (40,FILE=NAMECTL,STATUS='UNKNOWN')
       WRITE(40,'(A,A)') 'DSET    ^', NAMEDAT2
       WRITE(40,'(A,A)') 'TITLE    ', TITLE

       WRITE(40,'(A,E12.5,A)'  ) 'UNDEF   ',UNDEF,'  ! gap value'
	  
	 IF(HX.GT.0.0) THEN
        WRITE(40,'(A,I6,A,2(F14.7,5X))') 
     &          'XDEF  ',NX,'  LINEAR   ', X0(1),HX
	 ELSE
       WRITE(40,'(A,I6,A,15(F14.7,1X)/(22X,15(F14.7,1X)))')
     &          'XDEF  ',NX,'  LEVELS  ',(X0(I),I=1,NX)	 
	 END IF

	 IF(HY.GT.0.0) THEN
        WRITE(40,'(A,I6,A,2(F14.7,5X))') 
     &          'YDEF  ',NY,'  LINEAR   ', Y0(1),HY
	 ELSE
        WRITE(40,'(A,I6,A,15(F14.7,1X)/(22X,15(F14.7,1X)))')
     &          'YDEF  ',NY,'  LEVELS  ',(Y0(I),I=1,NY)	 
	 END IF

C TIME GRID:
      CALL SEC_TO_YR_MO_HR_MN(T0,NYRT0,NMOT0,NDYT0,NHRT0,NMNT0)
              IF(NMOT0.NE.0) THEN
              MON=MONTH(NMOT0)
              END IF

      CALL SEC_TO_YR_MO_HR_MN(HT,NYRHT,NMOHT,NDYHT,NHRHT,NMNHT)

      IF(NYRHT.NE.0) THEN
       TFSTEP='yr'
       K=NYRHT
      ELSEIF(NMOHT.NE.0) THEN
       TFSTEP='mo'
       K=NMOHT
      ELSEIF(NDYHT.NE.0) THEN
       TFSTEP='dy'
       K=NDYHT
      ELSEIF(NHRHT.NE.0) THEN
       TFSTEP='hr'
       K=NHRHT
      ELSEIF(NMNHT.NE.0) THEN
       TFSTEP='mn'
       K=NMNHT
      ELSE
      WRITE(*,*)' WARNING! TIME STEP ERROR IN CREATING CTL FILE',NAMECTL
      END IF

        WRITE(40,'(A,I6,A,I2.2,A,I2.2,A,I2.2,A,I4.4,5X,I4,A)')
     &      'TDEF  ',NT,'  LINEAR  ',NHRT0,':',NMNT0,'Z',
     &   NDYT0,MON,NYRT0,K,TFSTEP

      WRITE(40,'(A,I6,A,5(F14.7,1X)/(22X,5(F14.7,1X)))')
     &           'ZDEF  ',NZ,'  LEVELS  ',Z
c      WRITE(40,*)
       WRITE(40,'(A)')  'VARS 1  ! number of variables'
       WRITE(40,'(A,A,I5,A)')
     &            VARNAME, '        ', NZ,  '  1 variables '
       WRITE(40,'(A)') 'ENDVARS'
       WRITE(40,*)
      CLOSE (40,STATUS='KEEP')
      END
C====================================================================
      SUBROUTINE EXTNAME(NAME,EXT)
C--------------------------------------------------------------------
C     Construction of extension  for the  NAME
c      INPUT : NAME
c              EXT  (must consist of for symbols '.dat' for example)
c         OUTPUT: NAME = NAME + EXT
C--------------------------------------------------------------------
      CHARACTER*(*)  NAME
      CHARACTER*(*)  EXT
      INTEGER        N
C     INTEGER        NCH, NCHB
      DO  N=LEN(NAME)-3,2,-1
      IF (Name(N:N).EQ.'.') GO TO 5
      END DO
      N=N+1
  5   NAME(N:N+4) = EXT

      RETURN
      END
C======================================================================
      SUBROUTINE SEC_TO_YR_MO_HR_MN(SEC,NYR,NMO,NDY,NHR,NMN)
      INTEGER NYR,NMO,NDY,NHR,NMN  !NUMBERS OF YEARS,MONTHS,DAYS,HOURS,MINUTES
      REAL SEC, !TIME IN SECONDS
     &     YR_SEC,MO_SEC,SEC1
      PARAMETER (YR_SEC=360.0*86400.0,MO_SEC=30.0*86400.0)

      SEC1=SEC+0.10
      NYR  = INT(SEC1/YR_SEC)
      SEC1=(SEC1-FLOAT(NYR)*YR_SEC)

      NMO  = INT(SEC1/MO_SEC)
      SEC1=(SEC1-FLOAT(NMO)*MO_SEC)

      NDY  = INT(SEC1/86400.0)
      SEC1=(SEC1-FLOAT(NDY)*86400.0)

      NHR  = INT(SEC1/3600.0)
      SEC1=(SEC1-FLOAT(NHR)*3600.0)

      NMN  = INT(SEC1/60.0)

      RETURN
      END
C======================================================================
      SUBROUTINE READ_CTLFILE(CTLNAME,DATNAME,UNDEF,NX,NY,NZ,NT,
     &               X0,HX,Y0,HY,T0,HT,Z0,HZ,NVAR,TITLE,VARNAME,INFORM)
	IMPLICIT NONE 
      CHARACTER*(*) CTLNAME,DATNAME
      CHARACTER*(*) TITLE
      INTEGER       NX, NY, NZ, NT,   !DIMENSION OF DATA
     &              NVAR              !NUMBER OF VARIABLES
      REAL  X0(1),HX,Y0(1),HY, !INITIAL POINTS AND STEPS[DEG] OF HORIZONTAL GRID
                               !IF HX OR HY <0, THEN LEVELS ARE USED FOR APPROPRIATE VARIABLE
     &      T0(1),HT,          !INITIAL TIME AND STEP[SEC] IN TIME
     &      Z0(1),HZ,          !INITIAL POINTS AND STEPS[DEG] OF VERTICAL GRID
                               !IF HZ <0, THEN LEVELS ARE USED 
     &      UNDEF              !UNDEFINIT VALUE

      CHARACTER*(*) VARNAME(1)
	CHARACTER(LEN=32) VARTYPE, GRIDTYPE, TIMEINIT, TIMESTEP

      CHARACTER(256) STRING
      INTEGER I,L,L1,L2,NB,LINE,INFORM,LPRINT

C  EXTERNAL FUNCTIONS
	INTEGER NONBLANK

	IF(INFORM.EQ.0) THEN

		LPRINT=0

	ELSE

		LPRINT=1

	END IF

		INFORM=0

      OPEN (40,FILE=CTLNAME,STATUS='OLD',ACTION='READ')
      
	DATNAME=''
		
	DO LINE=1,512
       
	 READ(40,'(A)',ERR=21) STRING
	 NB=NONBLANK(STRING,1)
	 IF(NB.EQ.0) GO TO 7
	 IF(STRING(NB:NB+6).EQ.'ENDVARS'.OR.
     &    STRING(NB:NB+6).EQ.'endvars'.OR.
     &    STRING(NB:NB+6).EQ.'Endvars'    ) THEN
		
		CLOSE (40)

c Printing information about CTL file if need
	IF(LPRINT.NE.0) THEN

       WRITE(*,'(A,A)') 'DSET    ^', DATNAME
       WRITE(*,'(A,A)') 'TITLE    ', TITLE

       WRITE(*,'(A,E12.5,A)'  ) 'UNDEF   ',UNDEF,'  ! gap value'
	  
	 IF(HX.GT.0.0) THEN
        WRITE(*,'(A,I6,A,2(G14.7,5X))') 
     &          'XDEF  ',NX,'  LINEAR   ', X0(1),HX
	 ELSE
       WRITE(*,'(A,I6,A,15(G14.7,1X)/(22X,15(G14.7,1X)))')
     &          'XDEF  ',NX,'  LEVELS  ',(X0(I),I=1,NX)	 
	 END IF

	 IF(HY.GT.0.0) THEN
        WRITE(*,'(A,I6,A,2(G14.7,5X))') 
     &          'YDEF  ',NY,'  LINEAR   ', Y0(1),HY
	 ELSE
        WRITE(*,'(A,I6,A,15(G14.7,1X)/(22X,15(G14.7,1X)))')
     &          'YDEF  ',NY,'  LEVELS  ',(Y0(I),I=1,NY)	 
	 END IF
		
	  
	  WRITE(*,*) 'TDEF  ',NT,'  LINEAR  ',TIMEINIT,TIMESTEP


      WRITE(*,'(A,I6,A,5(G14.7,1X)/(22X,5(G14.7,1X)))')
     &           'ZDEF  ',NZ,'  LEVELS  ',(Z0(I),I=1,NZ)
c      WRITE(40,*)
       WRITE(*,'(A)')  'VARS 1  ! number of variables'
       WRITE(*,'(A,A,I5,A)')
     &            VARNAME, '        ', NZ,  '  1 variables '
       WRITE(*,'(A)') 'ENDVARS'
       WRITE(*,*)

	END IF


	    RETURN
	 END IF
C DEFINE NUMBER OF VARIABLES
	 IF(STRING(NB:NB+3).EQ.'VARS'.OR.
     &    STRING(NB:NB+3).EQ.'Vars'.OR.
     &    STRING(NB:NB+3).EQ.'vars'    ) THEN
         
	    READ(STRING,*,ERR=22) VARTYPE,NVAR 
		
		DO L=1,NVAR
		READ(40,*,ERR=22) VARNAME(L),I
		END DO
	    INFORM=INFORM+1

	 END IF 

C DEFINE NAME OF FILE WITH DATA
	 IF(STRING(NB:NB+3).EQ.'DSET'.OR.
     &    STRING(NB:NB+3).EQ.'dset'.OR.
     &    STRING(NB:NB+3).EQ.'Dset'    ) THEN
          
		DO  L1=NB+4,LEN(STRING)
          IF (STRING(L1:L1).NE.' '.AND.STRING(L1:L1).NE.'^') GO TO 5
          END DO
5       CONTINUE
		DO  L2=MIN(LEN(STRING),LEN(DATNAME)+L1-1),NB+4,-1
          IF (STRING(L2:L2).NE.' ') GO TO 6
          END DO 
	    L2=L2-1
6       CONTINUE
        
	  IF(L1.GT.L2) THEN		       
	     WRITE(*,'(A,A)') ' ERROR IN ',CTLNAME
	     WRITE(*,*) ' ERROR IN DEFINITION NAME OF FILE WITH DATA!'
	     RETURN
	  ELSE
           DATNAME(1:L2-L1+1)=STRING(L1:L2)
	     INFORM=INFORM+1
	  END IF

	 END IF

C DEFINE TITLE OF DATA
	 IF(STRING(NB:NB+4).EQ.'TITLE'.OR.
     &    STRING(NB:NB+4).EQ.'Title'.OR.
     &    STRING(NB:NB+4).EQ.'title'    ) THEN

		READ(STRING,*,ERR=14) VARTYPE,TITLE
          INFORM=INFORM+1

	 END IF

C DEFINE UNDEFINIT VALUE
	 IF(STRING(NB:NB+4).EQ.'UNDEF'.OR.
     &    STRING(NB:NB+4).EQ.'undef'.OR.
     &    STRING(NB:NB+4).EQ.'Undef'    ) THEN
         
	    READ(STRING,*,ERR=8) VARTYPE,UNDEF 
	    INFORM=INFORM+1
	 END IF   
     
C DEFINE X GRID
	 IF(STRING(NB:NB+3).EQ.'XDEF'.OR.
     &    STRING(NB:NB+3).EQ.'Xdef'.OR.
     &    STRING(NB:NB+3).EQ.'xdef'    ) THEN

		READ(STRING,*,ERR=10) VARTYPE,NX,GRIDTYPE

          IF(GRIDTYPE.EQ.'LINEAR'.OR.
     &       GRIDTYPE.EQ.'linear'.OR.
     &       GRIDTYPE.EQ.'Linear'    ) THEN
		
		READ(STRING,*,ERR=10) VARTYPE,NX,GRIDTYPE,X0(1),HX         
	    
		   DO I=2,NX
		   X0(I)=X0(I-1)+HX
		   END DO
	       INFORM=INFORM+1

	    ELSE IF(GRIDTYPE.EQ.'LEVELS'.OR.
     &            GRIDTYPE.EQ.'levels'.OR.
     &            GRIDTYPE.EQ.'Levels'    ) THEN

		        BACKSPACE (40)
	            READ(40,*,ERR=10) VARTYPE,L,GRIDTYPE,(X0(I),I=1,NX)
	            HX=0.0
	            INFORM=INFORM+1
		
		ELSE 
  

	    END IF

	 END IF
	    

C DEFINE Y GRID
	 IF(STRING(NB:NB+3).EQ.'YDEF'.OR.
     &    STRING(NB:NB+3).EQ.'Ydef'.OR.
     &    STRING(NB:NB+3).EQ.'ydef'    ) THEN

		READ(STRING,*,ERR=11) VARTYPE,NY,GRIDTYPE

          IF(GRIDTYPE.EQ.'LINEAR'.OR.
     &       GRIDTYPE.EQ.'linear'.OR.
     &       GRIDTYPE.EQ.'Linear'    ) THEN
		
		READ(STRING,*,ERR=11) VARTYPE,NY,GRIDTYPE,Y0(1),HY  
		       
	       DO I=2,NY
		   Y0(I)=Y0(I-1)+HY
		   END DO
	       INFORM=INFORM+1

	    ELSE IF(GRIDTYPE.EQ.'LEVELS'.OR.
     &            GRIDTYPE.EQ.'levels'.OR.
     &            GRIDTYPE.EQ.'Levels'    ) THEN

		        BACKSPACE (40)
	            READ(40,*,ERR=11) VARTYPE,L,GRIDTYPE,(Y0(I),I=1,NY)
                  HY=0.0
	            INFORM=INFORM+1
		
		ELSE 
		

	    END IF

	 END IF

C DEFINE Z GRID
	 IF(STRING(NB:NB+3).EQ.'ZDEF'.OR.
     &    STRING(NB:NB+3).EQ.'Zdef'.OR.
     &    STRING(NB:NB+3).EQ.'zdef'    ) THEN

		READ(STRING,*,ERR=12) VARTYPE,NZ,GRIDTYPE

          IF(GRIDTYPE.EQ.'LINEAR'.OR.
     &       GRIDTYPE.EQ.'linear'.OR.
     &       GRIDTYPE.EQ.'Linear'    ) THEN

 		READ(STRING,*,ERR=12) VARTYPE,NZ,GRIDTYPE,Z0(1),HZ        

	       DO I=2,NZ
		   Z0(I)=Z0(I-1)+HZ
		   END DO
	       INFORM=INFORM+1

	    ELSE IF(GRIDTYPE.EQ.'LEVELS'.OR.
     &            GRIDTYPE.EQ.'levels'.OR.
     &            GRIDTYPE.EQ.'Levels'    ) THEN

		        BACKSPACE (40)
	            READ(40,*,ERR=12) VARTYPE,L,GRIDTYPE,(Z0(I),I=1,NZ)
	            HZ=0.0
	            INFORM=INFORM+1
		
		ELSE 
		

	    END IF

	 END IF
	    

                
C DEFINE T GRID
	 IF(STRING(NB:NB+3).EQ.'TDEF'.OR.
     &    STRING(NB:NB+3).EQ.'tdef'.OR.
     &    STRING(NB:NB+3).EQ.'Tdef'    ) THEN
	    
		READ(STRING,*,ERR=13) VARTYPE,NT,GRIDTYPE,TIMEINIT,TIMESTEP

          IF(GRIDTYPE.EQ.'LINEAR'.OR.
     &       GRIDTYPE.EQ.'linear'.OR.
     &       GRIDTYPE.EQ.'Linear'    ) THEN
		
		DO  L2=LEN(TIMESTEP),2,-1
          IF (TIMESTEP(L2:L2).NE.' ') GO TO 26
          END DO 

		   		   		
26        CONTINUE			

		READ(TIMESTEP(1:L2-2),*,ERR=13) HT

          IF    (TIMESTEP(L2-1:L2).EQ.'yr') THEN
                 HT=HT*365.0*24.0*3600.0

          ELSEIF(TIMESTEP(L2-1:L2).EQ.'mo') THEN
                 HT=HT*30.0*24.0*3600.0

          ELSEIF(TIMESTEP(L2-1:L2).EQ.'dy') THEN
                 HT=HT     *24.0*3600.0

          ELSEIF(TIMESTEP(L2-1:L2).EQ.'hr') THEN
                 HT=HT          *3600.0

          ELSEIF(TIMESTEP(L2-1:L2).EQ.'mn') THEN
                 HT=HT          *60.0

          ELSEIF(TIMESTEP(L2-1:L2).EQ.'mo') THEN
                 HT=30.0*24.0*3600.0

          ELSE
	       WRITE(*,'(A,A)') ' ERROR IN ',CTLNAME
	       WRITE(*,*) ' IN DEFINITION TIME STEP FOR T-GRID!'
             CLOSE (40)
             RETURN
          END IF

          
c         READ(STRING(NB+4:L1-1),*,ERR=11) NT
c	    READ(STRING(L1+7:LEN(STRING)),*,ERR=11) T0(1),HT
c	         DO I=2,NT
c			 T0(I)=T0(I-1)+HT
c			 END DO
c	    INFORM=INFORM+1

	    ELSE

	    WRITE(*,'(A,A)') ' WARNING IN ',CTLNAME
	    WRITE(*,*) ' GRID TIME IS NOT LINEAR!'

	    END IF

	 END IF   


7	CONTINUE
      END DO


8       CONTINUE
	  WRITE(*,'(A,A)') ' ERROR IN ',CTLNAME
	  WRITE(*,*) ' IN DEFINITION UNDEFINIT VALUE!'
        CLOSE (40)
      RETURN

10      CONTINUE
	  WRITE(*,'(A,A)') ' ERROR IN ',CTLNAME
	  WRITE(*,*) ' IN DEFINITION X-GRID!'
        CLOSE (40)
      RETURN	

11      CONTINUE
	  WRITE(*,'(A,A)') ' ERROR IN ',CTLNAME
	  WRITE(*,*) ' IN DEFINITION Y-GRID!'
        CLOSE (40)
      RETURN
						
12      CONTINUE
	  WRITE(*,'(A,A)') ' ERROR IN ',CTLNAME
	  WRITE(*,*) ' IN DEFINITION Z-GRID!'
        CLOSE (40)
      RETURN
						
13      CONTINUE
	  WRITE(*,'(A,A)') ' ERROR IN ',CTLNAME
	  WRITE(*,*) ' IN DEFINITION T-GRID!'
        CLOSE (40)
      RETURN

14      CONTINUE
	  WRITE(*,'(A,A)') ' ERROR IN ',CTLNAME
	  WRITE(*,*) ' IN DEFINITION TITLE OF DATA!'
        CLOSE (40)
      RETURN

21      CONTINUE
	  WRITE(*,*) ' WARNING IN ',CTLNAME
	  WRITE(*,*) ' EXIT BY EOF, NOT BY ENDVARS!'
        CLOSE (40)	  	
	RETURN

22      CONTINUE
	  WRITE(*,'(A,A)') ' ERROR IN ',CTLNAME
	  WRITE(*,*) ' IN DEFINITION VARS OR VARIABLE NAME!'
        CLOSE (40)
      RETURN

      END
C====================================================================
      INTEGER FUNCTION NONBLANK(STRING,NB1)
	IMPLICIT NONE
C--------------------------------------------------------------------
C     Calculate non blank symbol position
c      INPUT : STRING - string
c              NB1    - BEGINNIG POSITION
c      OUTPUT: NONBLANK  - non blank position
C--------------------------------------------------------------------
      CHARACTER*(*)  STRING
	INTEGER NB,NB1

      DO  NB=NB1,LEN(STRING)
      IF (STRING(NB:NB).NE.' ') GO TO 5
      END DO

	NONBLANK=0
	RETURN

5     CONTINUE

      NONBLANK=NB

      IF (NB.GT.LEN(STRING)-5) THEN 
	 WRITE(*,*) ' WARNING FOR NON BLANK POSITION IN STRING:'
	 WRITE(*,'(A)') STRING
 	 NONBLANK=0
      END IF

      RETURN
      END
C====================================================================
      INTEGER FUNCTION CHARPOSITION(STRING,NB1,CHAR)
	IMPLICIT NONE
C--------------------------------------------------------------------
C     Calculate non blank symbol position
c      INPUT : STRING - string
c              NB1    - BEGINNIG POSITION
c      OUTPUT: CHARPOSITION - position of character CHAR
C--------------------------------------------------------------------
      CHARACTER*(*)  STRING
	CHARACTER CHAR*1
	INTEGER NB,NB1

      DO  NB=NB1,LEN(STRING)
      IF (STRING(NB:NB).EQ.CHAR) GO TO 5
      END DO

	 WRITE(*,*) ' WARNING IN STRING:'
	 WRITE(*,'(A)') STRING
 	 WRITE(*,'(A,A)') ' THERE IS NO CHARACTER ', CHAR
	 CHARPOSITION=0
	RETURN

5     CONTINUE

      CHARPOSITION=NB
      RETURN
      END
C======================================================================
	SUBROUTINE COMPARE_GRIDS(MAXGRIDLEN,IERR,
     &	         UNDEFA,NXA,NYA,NZA,NTA,XA,HXA,YA,HYA,TA,HTA,ZA,HZA,
     &             UNDEFB,NXB,NYB,NZB,NTB,XB,HXB,YB,HYB,TB,HTB,ZB,HZB)
	IMPLICIT NONE
	INTEGER MAXGRIDLEN,IERR,I
C PARAMETERS FOR TEMPERATURE
	REAL  UNDEFA,XA(MAXGRIDLEN),HXA,
     &             YA(MAXGRIDLEN),HYA,TA,HTA,
     &             ZA(MAXGRIDLEN),HZA 
	INTEGER NXA,NYA,NZA,NTA
 
C PARAMETERS FOR SALINITY
	REAL  UNDEFB,XB(MAXGRIDLEN),HXB,
     &             YB(MAXGRIDLEN),HYB,TB,HTB,
     &             ZB(MAXGRIDLEN),HZB 
	INTEGER NXB,NYB,NZB,NTB

	IERR=0

	IF(NXA.GT.MAXGRIDLEN) THEN
	IERR=1
	WRITE(*,'(A,I7,A,I7)')
     &      ' NUBERS OF X-GRID POINTS ',NXA,' > ',' MAX BE',MAXGRIDLEN
	RETURN
	END IF

	IF(NYA.GT.MAXGRIDLEN) THEN
	IERR=1
	WRITE(*,'(A,I7,A,I7)')
     &      ' NUBERS OF Y-GRID POINTS ',NYA,' > ',' MAX BE',MAXGRIDLEN
	RETURN
	END IF
		
	IF(NZA.GT.MAXGRIDLEN) THEN
	IERR=1
	WRITE(*,'(A,I7,A,I7)')
     &      ' NUBERS OF Z-GRID POINTS ',NZA,' > ',' MAX BE',MAXGRIDLEN
	RETURN
	END IF

	IF(NXA.NE.NXB) THEN
	IERR=1
	WRITE(*,'(A,I7,A,A,I7)')
     &      ' NUBERS OF A X-GRID POINTS ',NXA,' DOES NOT EQUAL TO',
     &      ' NUBERS OF B X-GRID POINTS ',NXB
	RETURN
	END IF
	
	IF(NYA.NE.NYB) THEN
	IERR=1
	WRITE(*,'(A,I7,A,A,I7)')
     &      ' NUBERS OF A Y-GRID POINTS ',NYA,' DOES NOT EQUAL TO',
     &      ' NUBERS OF B Y-GRID POINTS ',NYB
	RETURN
	END IF

	IF(NZA.NE.NZB) THEN
	IERR=1
	WRITE(*,'(A,I7,A,A,I7)')
     &      ' NUBERS OF A Z-GRID POINTS ',NZA,' DOES NOT EQUAL TO',
     &      ' NUBERS OF B Z-GRID POINTS ',NZB
	RETURN
	END IF

	IF(NTA.NE.NTB) THEN
	IERR=1
	WRITE(*,'(A,I7,A,A,I7)')
     &      ' NUBERS OF A Y-GRID POINTS ',NYA,' DOES NOT EQUAL TO',
     &      ' NUBERS OF B Y-GRID POINTS ',NYB
	RETURN
	END IF

	DO I=1,NXA
	 IF(XA(I).NE.XB(I)) THEN
	 WRITE(*,*)' X-GRIDS (A) AND (B) ARE DIFFERENT!'
	 IERR=2
	 END IF
	END DO

	DO I=1,NYA
	 IF(YA(I).NE.YB(I)) THEN
	 WRITE(*,*)' Y-GRIDS (A) AND (B) ARE DIFFERENT!'
	 IERR=2
	 END IF
	END DO

	DO I=1,NZA
	 IF(ZA(I).NE.ZB(I)) THEN
	 WRITE(*,*)' Z-GRIDS (A) AND (B) ARE DIFFERENT!'
	 IERR=2
	 END IF
	END DO
			
	RETURN
	END

