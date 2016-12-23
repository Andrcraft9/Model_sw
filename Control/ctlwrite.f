C=======================================================================
      SUBROUTINE XY_CTL(PATH2DATA,NUMOFREC,STARTIME,TSTEP)
	  use mpi_parallel_tools
	  use mod_shallowater
C WRITING 3D XY-FILD DATA TO NUMOFREC
      CHARACTER*(*) PATH2DATA    !PATH TO OUTPUT FILES
      INTEGER NUMOFREC           !NUMBER OF RECORD
      REAL    STARTIME,TSTEP   !START-TIME AND INCREMENT IN SECONDS
      INCLUDE '0COM.INC'
      INCLUDE '0CEAN.INC'
      INCLUDE '1LREC.INC'
	  INCLUDE '1STDOUT.INC'
	  INCLUDE '2STDOUT.INC'

      CHARACTER FNAME*128
      INTEGER  IERR

C======================= GEO_OUTPUT.EQ.0 ========================================
C WRITING BAROTROPIC COMPONETS

C STORE BAROTROPIC VELOCITIES:
c        CALL FULFNAME(FNAME,PATH2DATA,'XY/ub.dat',IERR)
c        CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,1,NUMOFREC,
c     &            XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &                0.0,'Barotropic Zonal Velocity[CM/S]','ub')

c        CALL FULFNAME(FNAME,PATH2DATA,'XY/vb.dat',IERR)
c        CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,1,NUMOFREC,
c     &            XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &           0.0,'Barotropic Meridional Velocity[CM/S]','vb')

C STORE BAROTROPIC VELOCITIES:
c        CALL FULFNAME(FNAME,PATH2DATA,'XY/taux.dat',IERR)
c        CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,1,NUMOFREC,
c     &            XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &                0.0,'Zonal wind stress [CM**2/S**2]','tx')

c        CALL FULFNAME(FNAME,PATH2DATA,'XY/tauy.dat',IERR)
c        CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,1,NUMOFREC,
c     &            XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &           0.0,'Meridional wind stress[CM**2/S**2]','ty')

C      STORE SEA LEVEL

      CALL FULFNAME(FNAME,PATH2DATA,'XY/sl.dat',IERR)
      CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,1,
     &           NUMOFREC,
     &           XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
     &           0.0,'Model Sea Level[CM]','sl')

C      STORE STREAMFUNCTION
c         CALL FULFNAME(FNAME,PATH2DATA,'XY/sf.dat',IERR)
c         CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+2,NN-NNN+2,1,NUMOFREC,
c     &        XU(MMM-1:MM),0.0,YV(NNN-1:NN),0.0,STARTIME,TSTEP,
c     &                     0.0,'Stream Function[SV]','sf')

c       IF(IABS(KSW_TSL).GT.1.OR.IABS(KSW_TSV).GT.1) THEN
C USE ARRAY VTR FOR WRITE MIXED LAYER DEPTH
c        CALL FULFNAME(FNAME,PATH2DATA,'XY/mxld.dat',IERR)
c        CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,1,NUMOFREC,
c     &          XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &                    0.0,'Model mixed layer depth[m]','mld')
c       END IF

C      STORE SEA LEVEL
c         CALL FULFNAME(FNAME,PATH2DATA,'XY/sls.dat',IERR)
c         CALL CREATE_CTLFILE(FNAME,UNDEF,MM-MMM+1,NN-NNN+1,1,NUMOFREC,
c     &           XT(MMM:MM),0.0,YT(NNN:NN),0.0,STARTIME,TSTEP,
c     &                    0.0,'Sterical Sea Level[CM]','sl')

      RETURN
      END
