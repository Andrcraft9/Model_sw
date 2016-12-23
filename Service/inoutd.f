C======================================================================
C  VERSION#1 OF INPUT-OUTPUT SUBROUTINES
C======================================================================
      SUBROUTINE FULFNAME(NAME,PATH,FILEN,IERR)
	IMPLICIT NONE
C  DEFINIT SUM: NAME = PATH + FILEN REMOVING BLANKS FROM PATH
      CHARACTER*(*) PATH, FILEN, NAME
      CHARACTER*1 FNDEVIDER
      INTEGER N, LP, IERR, NP1, NP2, LF, NF1, NF2
C  FNDEVIDER - DEVIDER IN FULL FILE NAME WITH PATH
C  CHANGE IT ACCORDING TO OPERATION SISTEM
C  FOR DOS:
C     FNDEVIDER='\'
C  FOR UNIX(DOS ALSO):
      FNDEVIDER='/'

      IERR=0
C     NAME=FILEN
C     RETURN
      DO N=1,LEN(NAME)
         NAME(N:N)=' '
      END DO
      LP=LEN(PATH)
      IF(LP.LT.1) THEN
         IERR=1
         WRITE(*,'(2X,A)')'ERROR IN SUBROUTINE FULFNAME:'
         WRITE(*,'(2X,A,A)')' ERROR IN PATH TO FILE ',
     &                        NAME(1:LEN_TRIM(NAME))
         RETURN
      END IF
C  FINED REAL INITIAL AND LAST DEFINIT POSITION IN PATH WITHOUT BLANK
          NP1=0
      IF(LP.GT.1) THEN
          N=1
          DO WHILE(N.LT.LP)
             IF(PATH(N:N).EQ.' '.AND.PATH(N+1:N+1).NE.' ') NP1=N
             IF(PATH(N:N).NE.' '.AND.PATH(N+1:N+1).EQ.' ') EXIT
             N=N+1
          END DO
          IF(N.EQ.LP.AND.PATH(LP:LP).EQ.' ') THEN
             NP2 = 0
          ELSE
             NP2 = N
          END IF
      ELSE
          IF(PATH(1:1).EQ.' ') THEN
              NP2=0
          ELSE
              NP2=1
          END IF
      END IF
      LP=NP2-NP1                !LONG OF REAL PATH WITOUT BLANKS

      LF=LEN(FILEN)
      IF(LF.LT.1) THEN
         IERR=1
         WRITE(*,'(2X,A)') 'ERROR IN SUBROUTINE FULFNAME:'
         WRITE(*,'(2X,A)') 'ERROR IN FILE NAME: '
         WRITE(*,*)'PATH: ', PATH(1:LEN_TRIM(PATH)),
     &          ';  NAME:',  NAME(1:LEN_TRIM(PATH))
         RETURN
      END IF
C  FINED REAL INITIAL AND LAST DEFINIT POSITION IN FILENAME WITHOUT BLANK
          NF1=0
      IF(LF.GT.1) THEN
          N=1
          DO WHILE(N.LT.LF)
             IF(FILEN(N:N).EQ.' '.AND.FILEN(N+1:N+1).NE.' ') NF1=N
             IF(FILEN(N:N).NE.' '.AND.FILEN(N+1:N+1).EQ.' ') EXIT
             N=N+1
          END DO
          IF(N.EQ.LF.AND.FILEN(LF:LF).EQ.' ') THEN
             NF2=0
          ELSE
             NF2 = N
          END IF
      ELSE
          IF(FILEN(1:1).EQ.' ') THEN
              NF2=0
          ELSE
              NF2=1
          END IF
      END IF
      LF=NF2-NF1                !LONG OF REAL FILENAME WITOUT BLANKS

      IF(LF.LE.0) THEN
      IERR=1
      WRITE(*,'(2X,A)')'ERROR IN SUBROUTINE FULFNAME:'
      WRITE(*,'(2X,A)')'THERE IS NO FILE NAME!'
      RETURN
      END IF

      IF(LP+LF.GT.LEN(NAME)) THEN
      IERR=1
      WRITE(*,'(2X,A)')'ERROR IN SUBROUTINE FULFNAME:'
      WRITE(*,'(2X,A)') 'ERROR IN FILE NAME: '
      WRITE(*,*)'PATH: ', PATH(1:LEN_TRIM(PATH)),
     &       ';  NAME:',  NAME(1:LEN_TRIM(NAME))
      WRITE(*,'(2X,A)')'LEN OF FULNAME < PATH+FILENAME:'
      RETURN
      END IF

      IF(LP.GT.0) THEN
          NAME(1:LP) = PATH(NP1+1:NP2)
          IF(NAME(LP:LP)       .EQ.FNDEVIDER.AND.
     &       FILEN(NF1+1:NF1+1).EQ.FNDEVIDER) THEN
            LP=LP-1
          END IF
          IF(NAME(LP:LP)       .NE.FNDEVIDER.AND.
     &       FILEN(NF1+1:NF1+1).NE.FNDEVIDER) THEN
            LP=LP+1
            NAME(LP:LP)=FNDEVIDER
          END IF
          NAME(LP+1:LP+LF) = FILEN(NF1+1:NF2)
      ELSE
          NAME(1:LF) = FILEN(NF1+1:NF2)
      END IF

      RETURN
      END

C======================================================================

      SUBROUTINE RDSTD(PATH,FNAME,NFILD,FILD,LU,NX,NY,NZ,
     &   NXB,NXE,NYB,NYE,NZB,NZE,IERR)
	IMPLICIT NONE
C  NX,NY,NZ - GENERAL DIMESION OF FILD
C  NXB,NXE,NYB,NYE,NZB,NZE - GRID COORDINATES OF TREAT ARRAY SUBDOMAIN
C      WHERE INDEX B DENOTES BEGIN, AND E - END

C     This subroutine fills (read) array FILD from unformatted deirect
C     file of DIOGIN standard
C     PATH           - path to file (i.g. 'F:\ARAB')
C     FNAME          - name of file (i.g.: 'taux.std')
C     NFILD          - number of field in file (on t)
C     FILD(NX,NY,NZ) - field array
C     LU(NX,NY)    - OCEAN MASK
C---------------------------------------------------------------------
      INCLUDE '1LREC.INC'       !SET LONG OF UNIQUE RECL
      CHARACTER*(*) PATH, FNAME
      INTEGER       NFILD,NRECF,NX,NY,NZ,I,J,K,IERR
      REAL          FILD(NX,NY,NZ), LU(NX,NY)
      CHARACTER(4096) NAMOFILE
      INTEGER  NXE, NXB, NYE, NYB, NZB, NZE, L, KR

C  DEFINITION FULL FILE NAME
      CALL FULFNAME(NAMOFILE,PATH,FNAME,IERR)
      IF (IERR.NE.0) GO TO 100

C  CHECK OF CORRECTNESS OF GRID COORDINATES OF TREAT ARRAY PART
      IF(NXE.GT.NX.OR.NXB.LT.1.OR.NXB.GT.NXE) THEN
          IERR=1
          GOTO 103
      END IF
      IF(NYE.GT.NY.OR.NYB.LT.1.OR.NYB.GT.NYE) THEN
          IERR=2
          GOTO 103
      END IF
      IF(NZE.GT.NZ.OR.NZB.LT.1.OR.NZB.GT.NZE) THEN
          IERR=3
          GOTO 103
      END IF
C      WRITE(*,'(2X,A)')'OPEN DIRECT FILE FOR READING:'
C      WRITE(*,'(2X,A)')  NAMOFILE
C      WRITE(*,'(2X,A,I4,A,I4)')'RECL LENTH IN WORDS:', NX,' X',NY
      OPEN(40,FILE=NAMOFILE,STATUS='OLD',ACCESS='DIRECT',
     &   FORM='UNFORMATTED',RECL=(NXE-NXB+1)*(NYE-NYB+1)*LRECL,ERR=101)

      NRECF=(NFILD-1)*(NZE-NZB+1)     !INITIAL NUMBER OF RECORD

      L=0
      KR=0
      DO K=NZB,NZE
      KR=KR+1
      READ(40,REC=NRECF+KR,ERR=102) ((FILD(I,J,K),I=NXB,NXE),J=NYB,NYE)

C  FULLING UNDEFINITE POINTS BY ZERO INSTED UNDEF
          DO J=NYB,NYE
            DO I=NXB,NXE

            IF(ABS(LU(I,J)).LT.0.5) THEN
               L=L+1
               FILD(I,J,K)=0.0
            END IF

            ENDDO
         ENDDO

      END DO

      CLOSE(40)

      L=(NXE-NXB+1)*(NYE-NYB+1)-L/KR
      WRITE(*,'(1X,A,A)')
     &         'INPUT DATA from ',NAMOFILE(1:LEN_TRIM (NAMOFILE))
      WRITE(*,'(7X,A,I7,A,I7,A,I8,A)')
     &         'DIMENSION OF FIELD =',NXE-NXB+1,' *',NYE-NYB+1,
     &        ' (',L,'-OCEAN POINTS)'
C     WRITE(*,'(2X,A)')'CLOSE DIRECT FILE:'
C     WRITE(*,'(2X,A)')  NAMOFILE
      RETURN
100   WRITE(*,'(2X,A)')'ERROR IN FULL NAME OF FILE FOR READING: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      STOP
101   WRITE(*,'(2X,A)')'ERROR IN OPEN FILE FOR READING: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      STOP
102   WRITE(*,'(2X,A)')'ERROR IN READING FROM FILE: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      WRITE(*,'(18H ERROR IN READING ,I3,6H LEVEL)') K
      STOP
103   WRITE(*,'(2X,A)')'ERROR IN READING FROM FILE: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      WRITE(*,'(2X,A,I3,A)')'ERROR IN GRID DIAPASON OF ',
     &                       IERR,' - COORDINATE'
      STOP
      END

C======================================================================

      subroutine  wdstd_parallel(path,fname,nfild,fild,LU,NX,NY,
     &   NXB,NXE,NYB,NYE)
      use mpi_parallel_tools
      implicit NONE
      include '1LREC.INC'

      character*(*) :: path, fname
      integer :: NX, NY, NXB, NXE, NYB, NYE
      integer :: nx_s, nx_e, ny_s, ny_e
      integer :: nfild, nrecf, offset, buffsize
      integer :: thefile, ierr, len_ierr, temp
      real :: LU(NX, NY)
c      real,allocatable :: fild(:, :)
      real :: fild(nx_start:nx_end, ny_start:ny_end)
      character(4096) :: namofile
      character(MPI_MAX_ERROR_STRING) :: s_ierr

      nx_s = max(nx_start, NXB)
      nx_e = min(nx_end, NXE)
      ny_s = max(ny_start, NYB)
      ny_e = min(ny_end, NYE)

      call FULFNAME(namofile, path, fname, ierr)

      call MPI_FILE_OPEN(CART_COMM, namofile,
     &                   MPI_MODE_WRONLY + MPI_MODE_CREATE,
     &                   MPI_INFO_NULL, thefile, ierr)


      nrecf = (nfild-1) * (NXE-NXB+1) * (NYE-NYB+1) * LRECL

      buffsize = (nx_e - nx_s + 1) * (ny_e - ny_s + 1)

      offset = ( ((nx_s - NXB) * (ny_e - ny_s + 1)) +
     &             ((ny_s - NYB) * (NXE - NXB + 1)) ) * LRECL

c      print *, rank, nx_s, nx_e, ny_s, ny_e
c      print *, rank, offset, nrecf+offset, buffsize

      call MPI_FILE_SET_VIEW(thefile,
     &                    int(nrecf + offset, KIND=MPI_OFFSET_KIND),
     &                    MPI_REAL4, MPI_REAL4, 'native',
     &                    MPI_INFO_NULL, ierr)

      call MPI_FILE_WRITE_ALL(thefile, fild(nx_s:nx_e, ny_s:ny_e),
     &                    buffsize, MPI_REAL4,
     &                    MPI_STATUS_IGNORE, ierr)

c      call MPI_Error_string(ierr, s_ierr, len_ierr, temp)
c      print *, rank, s_ierr(1:len_ierr)

      call MPI_FILE_CLOSE(thefile, ierr)

      call MPI_FINALIZE(ierr)

      stop

      return
      end subroutine

C======================================================================

      SUBROUTINE PWDSTD(PATH,FNAME,NFILD,FILD,LU,
     &   NX,OFFNXIN,ENDNXIN,NY,OFFNYIN,ENDNYIN,NZ,
     &   NXB,NXE,NYB,NYE,NZB,NZE,IERRR)
      use mpi_parallel_tools
      IMPLICIT NONE
C  NX,NY,NZ - GENERAL DIMESION OF FILD
C  NXB,NXE,NYB,NYE,NZB,NZE - GRID COORDINATES OF TREAT ARRAY SUBDOMAIN
C      WHERE INDEX B DENOTES BEGIN, AND E - END
C     This subroutine fills (WRITE) array FILD to unformatted deirect
C     file of DIOGIN standard
C     PATH           - path to file (i.g. 'F:\ARAB')
C     FNAME          - name of file (i.g.: 'taux.std')
C     NFILD          - number of field in file (on t)
C     FILD(NX,NY,NZ) - field array
C     LU(NX,NY)    - OCEAN MASK
C---------------------------------------------------------------------
      INCLUDE '1LREC.INC'       !SET LONG OF UNIQUE RECL

      CHARACTER CHUNKED*2048,CHUNKED_SIZE*2048,STRIPING_UNIT*2048,
     &          BUFFER_SIZE*2048,CHUNKED_ITEM*2048
      CHARACTER*(*) PATH, FNAME
      INTEGER HFILE,FI
      INTEGER(KIND=MPI_OFFSET_KIND) DISP
      INTEGER       OFFNXIN,ENDNXIN,OFFNYIN,ENDNYIN
      INTEGER       OFFNXINNB,ENDNXINNB,OFFNYINNB,ENDNYINNB
      INTEGER       NFILD,NRECF,NX,NY,NZ,I,J,K,IERR,IERRR
      REAL          FILD(OFFNXIN:ENDNXIN,OFFNYIN:ENDNYIN,NZ)
      REAL          LU(NX,NY)
      REAL,ALLOCATABLE :: FILDALL(:,:),LUALL(:,:)
      CHARACTER*4096 NAMOFILE
      INTEGER      NXE, NXB, NYE, NYB, NZB, NZE, L, KR
      INTEGER TSUBARR, SIZES3(3),LOCSIZES3(3),OFFSET3(3),TOTSIZE
      INTEGER SIZES2(2), LOCSIZES2(2), OFFSET2(2)
c      CALL GETPARAREANB(NX,OFFNXINNB,ENDNXINNB,1)
c      CALL GETPARAREANB(NY,OFFNYINNB,ENDNYINNB,2)
      OFFNXINNB = nx_start
      ENDNXINNB = nx_end
      OFFNYINNB = ny_start
      ENDNYINNB = ny_end
C  DEFINITION FULL FILE NAME
      CALL FULFNAME(NAMOFILE,PATH,FNAME,IERRR)
!      WRITE(*,*) "WRITING TO ", TRIM(NAMOFILE)
      IF (IERRR.NE.0) GO TO 100
C  CHECK OF CORRECTNESS OF GRID COORDINATES OF TREAT ARRAY PART
      IF(NXE.GT.NX.OR.NXB.LT.1.OR.NXB.GT.NXE) THEN
          IERRR=1
          GOTO 103
      END IF
      IF(NYE.GT.NY.OR.NYB.LT.1.OR.NYB.GT.NYE) THEN
          IERRR=2
          GOTO 103
      END IF
      IF(NZE.GT.NZ.OR.NZB.LT.1.OR.NZB.GT.NZE) THEN
          IERRR=3
          GOTO 103
      END IF
      DISP = (NXE-NXB+1)*(NYE-NYB+1)*LRECL*(NFILD-1)*(NZE-NZB+1)
      OFFSET3 = (/MAX(OFFNXINNB-NXB,0),MAX(OFFNYINNB-NYB,0),0/)
      OFFSET2 = (/MAX(OFFNXINNB-NXB,0),MAX(OFFNYINNB-NYB,0)/)
      LOCSIZES3 = (/ENDNXINNB-OFFNXINNB+1,ENDNYINNB-OFFNYINNB+1,
     &                                                NZE-NZB+1/)
      LOCSIZES2 = (/ENDNXINNB-OFFNXINNB+1,ENDNYINNB-OFFNYINNB+1/)
      IF( p_coord(1).EQ.0 ) THEN
        LOCSIZES3(1) = ENDNXINNB-NXB+1
        LOCSIZES2(1) = ENDNXINNB-NXB+1
      ELSE IF( p_coord(1).EQ.p_size(1)-1 ) THEN
        LOCSIZES3(1) = NXE-OFFNXINNB+1
        LOCSIZES2(1) = NXE-OFFNXINNB+1
      END IF
      IF( p_coord(2).EQ.0 ) THEN
        LOCSIZES3(2) = ENDNYINNB-NYB+1
        LOCSIZES2(2) = ENDNYINNB-NYB+1
      ELSE IF( p_coord(2).EQ.p_size(2)-1 ) THEN
        LOCSIZES3(2) = NYE-OFFNYINNB+1
        LOCSIZES2(2) = NYE-OFFNYINNB+1
      END IF
      IF( p_size(1).EQ.1 ) THEN
        LOCSIZES3(1) = NXE-NXB+1
        LOCSIZES2(1) = NXE-NXB+1
      END IF
      IF( p_size(2).EQ.1 ) THEN
        LOCSIZES3(2) = NYE-NYB+1
        LOCSIZES2(2) = NYE-NYB+1
      END IF
      SIZES3 = (/NXE-NXB+1,NYE-NYB+1,NZE-NZB+1/)
      SIZES2 = (/NXE-NXB+1,NYE-NYB+1/)
      TOTSIZE = LOCSIZES3(1)*LOCSIZES3(2)*LOCSIZES3(3)
      CALL MPI_INFO_CREATE(FI,IERR)
      CALL MPI_INFO_SET(FI,"IBM_largeblock_io","true",IERR)
      CALL MPI_INFO_SET(FI,"access_style",
     &                      "write_mostly",IERR)
      CALL MPI_INFO_SET(FI,"collective_buffering",
     &                     "true",IERR)
      WRITE(CHUNKED,1000) NY,NX,NZ
      CALL MPI_INFO_SET(FI,"chunked",CHUNKED,IERR)
      WRITE(CHUNKED_ITEM,1001) LRECL
      CALL MPI_INFO_SET(FI,"chunked_item",
     &                     CHUNKED_ITEM,IERR)
      WRITE(CHUNKED_SIZE,1000) ENDNYINNB-OFFNYINNB+1,
     &                         ENDNXINNB-OFFNXINNB+1,NZE-NZB+1
      CALL MPI_INFO_SET(FI,"chunked_size",
     &                         CHUNKED_SIZE,IERR)
      WRITE(STRIPING_UNIT,1002) (NX*NY*NZ)*LRECL
      CALL MPI_INFO_SET(FI,"striping_unit",
     &                     STRIPING_UNIT,IERR)
!      CALL MPI_INFO_SET(FI,"striping_factor",
!     &                    "16",IERR)
!      CALL MPI_INFO_SET(FI,"cb_nodes",
!     &                     "4",IERR)
      CALL MPI_INFO_SET(FI,"cb_buffer_size",
     &                     STRIPING_UNIT,IERR)
      WRITE(BUFFER_SIZE,1002)
     &(ENDNXINNB-OFFNXINNB+1)*(ENDNYINNB-OFFNYINNB+1)*(NZB-NZE+1)*LRECL
      CALL MPI_INFO_SET(FI,"ind_wr_buffer_size",
     &                  BUFFER_SIZE,IERR)
      CALL MPI_INFO_SET(FI,"cb_block_size",
     &                  BUFFER_SIZE,IERR)
      CALL MPI_FILE_OPEN(CART_COMM,NAMOFILE,
     &                   IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),
     &                   FI,HFILE,IERR)
      IF( IERR.NE.MPI_SUCCESS ) GOTO 101
      IF(SIZES3(3).EQ.1 ) THEN
!        print *, sizes2, locsizes2, offset2
      CALL MPI_TYPE_CREATE_SUBARRAY(2,SIZES2,LOCSIZES2,OFFSET2,
     &            MPI_ORDER_FORTRAN,MPI_REAL,TSUBARR,IERR)
      ELSE
!        print *, sizes3, locsizes3, offset3
      CALL MPI_TYPE_CREATE_SUBARRAY(3,SIZES3,LOCSIZES3,OFFSET3,
     &            MPI_ORDER_FORTRAN,MPI_REAL,TSUBARR,IERR)
      END IF
      CALL MPI_TYPE_COMMIT(TSUBARR,IERR)
      CALL MPI_FILE_SET_VIEW(HFILE,DISP,MPI_REAL,
     &                       TSUBARR,"native",FI,IERR)
!$OMP PARALLEL DO PRIVATE(I,J,K)
      DO J = MAX(NYB,OFFNYINNB),MIN(NYE,ENDNYINNB)
        DO I = MAX(NXB,OFFNXINNB),MIN(NXE,ENDNXINNB)
          IF( ABS(LU(I,J)).LT.0.5 ) THEN
            FILD(I,J,NZB:NZE) = UNDEF
          END IF
        END DO
      END DO
!$OMP END PARALLEL DO
      CALL MPI_FILE_WRITE_ALL(HFILE,
     & FILD(MAX(NXB,OFFNXINNB):MIN(NXE,ENDNXINNB),
     &      MAX(NYB,OFFNYINNB):MIN(NYE,ENDNYINNB),NZB:NZE),
     & TOTSIZE,MPI_REAL,MPI_STATUS_IGNORE,IERR)
      IF( IERR.NE.MPI_SUCCESS ) GOTO 102

      CALL MPI_TYPE_FREE(TSUBARR,IERR)
      CALL MPI_INFO_FREE(FI,IERR)
      CALL MPI_FILE_CLOSE(HFILE,IERR)
      L = 0
      IERRR = 0
      DO I = MAX(NXB,OFFNXINNB),MIN(NXE,ENDNXINNB)
        DO J = MAX(NYB,OFFNYINNB),MIN(NYE,ENDNYINNB)
          IF( ABS(LU(I,J)).LT.0.5 ) THEN
            L = L + NZE-NZB+1
            IERRR = IERRR + NZE-NZB+1
            FILD(I,J,NZB:NZE) = 0
          END IF
        END DO
      END DO
      L=(MIN(NXE,ENDNXINNB)-MAX(NXB,OFFNXINNB)+1)*
     &  (MIN(NYE,ENDNYINNB)-MAX(NYB,OFFNYINNB)+1)-L/(NZE-NZB+1)
      IF(rank.EQ.0) WRITE(*,'(1X,A,A)')
     &         'OUTPUT DATA to ',NAMOFILE(1:LEN_TRIM (NAMOFILE))
C      WRITE(*,'(8X,A,I7,A,I7,A,I8,A,I3)')
C     &         'DIMENSION OF FIELD =',NXE-NXB+1,' *',NYE-NYB+1,
C     &        ' (',L,'-OCEAN POINTS) ON RANK ', RANK
      CLOSE(40)
C      WRITE(*,'(2X,A)')'CLOSE DIRECT FILE:'
C      WRITE(*,'(2X,A)') NAMOFILE
      IERRR=(MIN(NXE,ENDNXINNB)-MAX(NXB,OFFNXINNB)+1)*
     &  (MIN(NYE,ENDNYINNB)-MAX(NYB,OFFNYINNB)+1)-IERRR/(NZE-NZB+1)-L
      IF (IERR.NE.0) THEN
C            WRITE(*,'(2X,A)')  NAMOFILE
C            WRITE(*,'(I7,A,A,I3)') IERRR,
C     &      'ERRORS IN NUMBER OF LOCAL OCEAN HORIZONTAL GRID POINTS.',
C     &      ' ON RANK ', RANK
      END IF
      RETURN
100   WRITE(*,'(2X,A)')'ERROR IN FULL NAME OF FILE FOR WRITING: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      CALL MPI_ABORT(CART_COMM,-1,IERR)
      STOP
101   WRITE(*,'(2X,A)')'ERROR IN OPEN FILE FOR WRITING: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      CALL MPI_ABORT(CART_COMM,-1,IERR)
      STOP
102   WRITE(*,'(2X,A)')'ERROR IN WRITING ON FILE: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      WRITE(*,'(18H ERROR IN WRITING ,I3,6H LEVEL)') K
      CALL MPI_ABORT(CART_COMM,-1,IERR)
      STOP
103   WRITE(*,'(2X,A)')'ERROR IN WRITING TO FILE: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      WRITE(*,'(2X,A,I3,A)')'ERROR IN GRID DIAPASON OF ',
     &                       IERRR,' - COORDINATE'
      CALL MPI_ABORT(CART_COMM,-1,IERR)
      STOP
1000  FORMAT(I6,',',I6,',',I6)
1001  FORMAT(I6)
1002  FORMAT(I18)
      END

C======================================================================

      SUBROUTINE WDSTD(PATH,FNAME,NFILD,FILD,LU,NX,NY,NZ,
     &   NXB,NXE,NYB,NYE,NZB,NZE,IERR)
	  IMPLICIT NONE
C  NX,NY,NZ - GENERAL DIMESION OF FILD
C  NXB,NXE,NYB,NYE,NZB,NZE - GRID COORDINATES OF TREAT ARRAY SUBDOMAIN
C      WHERE INDEX B DENOTES BEGIN, AND E - END

C     This subroutine fills (WRITE) array FILD to unformatted deirect
C     file of DIOGIN standard
C     PATH           - path to file (i.g. 'F:\ARAB')
C     FNAME          - name of file (i.g.: 'taux.std')
C     NFILD          - number of field in file (on t)
C     FILD(NX,NY,NZ) - field array
C     LU(NX,NY)    - OCEAN MASK
C     IERR - ERROR INFORMATION. IF IERR NE 0 IN INPUT, THEN NO PRINTING
C---------------------------------------------------------------------
      INCLUDE '1LREC.INC'       !SET LONG OF UNIQUE RECL
      CHARACTER*(*) PATH, FNAME
      INTEGER       NFILD,NRECF,NX,NY,NZ,I,J,K,IERR
      REAL          FILD(NX,NY,NZ), LU(NX,NY)
      CHARACTER(4096) NAMOFILE
      INTEGER      NXE, NXB, NYE, NYB, NZB, NZE, L, KR, LPRINT

	IF(IERR.EQ.0) THEN
	  LPRINT=1
	ELSE
	  LPRINT=0
	END IF

	IERR=0

C  DEFINITION FULL FILE NAME
      CALL FULFNAME(NAMOFILE,PATH,FNAME,IERR)
      IF (IERR.NE.0) GO TO 100

C  CHECK OF CORRECTNESS OF GRID COORDINATES OF TREAT ARRAY PART
      IF(NXE.GT.NX.OR.NXB.LT.1.OR.NXB.GT.NXE) THEN
          IERR=1
          GOTO 103
      END IF
      IF(NYE.GT.NY.OR.NYB.LT.1.OR.NYB.GT.NYE) THEN
          IERR=2
          GOTO 103
      END IF
      IF(NZE.GT.NZ.OR.NZB.LT.1.OR.NZB.GT.NZE) THEN
          IERR=3
          GOTO 103
      END IF
C      WRITE(*,'(2X,A)')'OPEN DIRECT FILE FOR WRITING:'
C      WRITE(*,'(2X,A)')  NAMOFILE
C      WRITE(*,'(2X,A,I4,A,I4)')'RECL LENTH IN WORDS:', NX,' X',NY
      OPEN(40,FILE=NAMOFILE,STATUS='UNKNOWN',ACCESS='DIRECT',
     &   FORM='UNFORMATTED',RECL=(NXE-NXB+1)*(NYE-NYB+1)*LRECL,ERR=101)

      NRECF=(NFILD-1)*(NZE-NZB+1)     !INITIAL NUMBER OF RECORD

      L   =0
      IERR=0
      KR=0
      DO K=NZB,NZE
         KR=KR+1
C  FULLING UNDEFINITE POINTS BY 0VER INSTED ZERO
          DO J=NYB,NYE
            DO I=NXB,NXE

            IF(ABS(LU(I,J)).LT.0.5) THEN
               L=L+1
               FILD(I,J,K)=UNDEF
            END IF

            ENDDO
         ENDDO

C  WRITING ON THE FILE
       WRITE(40,REC=NRECF+KR,ERR=102)((FILD(I,J,K),I=NXB,NXE),J=NYB,NYE)

C  FULLING UNDEFINITE POINTS BY ZERO INSTED UNDEF
          DO J=NYB,NYE
            DO I=NXB,NXE

            IF(ABS(LU(I,J)).LT.0.5) THEN
               IERR=IERR+1
               FILD(I,J,K)=0.0
            END IF

            ENDDO
         ENDDO
      END DO

      L=(NXE-NXB+1)*(NYE-NYB+1)-L/KR

	IF(LPRINT.EQ.1) THEN
C PRINT INFORMATION ON TERMINAL
      WRITE(*,'(1X,A,A)')
     &         'OUTPUT DATA to ',NAMOFILE(1:LEN_TRIM (NAMOFILE))
      WRITE(*,'(8X,A,I7,A,I7,A,I8,A)')
     &         'DIMENSION OF FIELD =',NXE-NXB+1,' *',NYE-NYB+1,
     &        ' (',L,'-OCEAN POINTS)'
      END IF

      CLOSE(40)
C      WRITE(*,'(2X,A)')'CLOSE DIRECT FILE:'
C      WRITE(*,'(2X,A)') NAMOFILE

      IERR=(NXE-NXB+1)*(NYE-NYB+1)-IERR/KR-L
      IF (IERR.NE.0) THEN
            WRITE(*,'(2X,A)')  NAMOFILE
            WRITE(*,'(I7,A)') IERR,
     &      'ERRORS IN NUMBER OF OCEAN HORIZONTAL GRID POINTS.'
      ENDIF

      RETURN
100   WRITE(*,'(2X,A)')'ERROR IN FULL NAME OF FILE FOR WRITING: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      STOP
101   WRITE(*,'(2X,A)')'ERROR IN OPEN FILE FOR WRITING: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      STOP
102   WRITE(*,'(2X,A)')'ERROR IN WRITING ON FILE: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      WRITE(*,'(18H ERROR IN WRITING ,I3,6H LEVEL)') K
      STOP
103   WRITE(*,'(2X,A)')'ERROR IN WRITING TO FILE: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      WRITE(*,'(2X,A,I3,A)')'ERROR IN GRID DIAPASON OF ',
     &                       IERR,' - COORDINATE'
      STOP
      END

C======================================================================

      SUBROUTINE RDSTD8(PATH,FNAME,NFILD,FILD,LU,NX,NY,NZ,
     &   NXB,NXE,NYB,NYE,NZB,NZE,IERR)
	IMPLICIT NONE
C  NX,NY,NZ - GENERAL DIMESION OF FILD
C  NXB,NXE,NYB,NYE,NZB,NZE - GRID COORDINATES OF TREAT ARRAY SUBDOMAIN
C      WHERE INDEX B DENOTES BEGIN, AND E - END

C     This subroutine fills (read) array FILD from unformatted deirect
C     file of DIOGIN standard
C     PATH           - path to file (i.g. 'F:\ARAB')
C     FNAME          - name of file (i.g.: 'taux.std')
C     NFILD          - number of field in file (on t)
C     FILD(NX,NY,NZ) - field array
C     LU(NX,NY)    - OCEAN MASK
C---------------------------------------------------------------------
      INCLUDE '1LREC.INC'       !SET LONG OF UNIQUE RECL
      CHARACTER*(*) PATH, FNAME
      INTEGER       NFILD,NRECF,NX,NY,NZ,I,J,K,IERR
      REAL*8        FILD(NX,NY,NZ)
      REAL            LU(NX,NY)
      CHARACTER(4096) NAMOFILE
      INTEGER  NXE, NXB, NYE, NYB, NZB, NZE, L, KR

C  DEFINITION FULL FILE NAME
      CALL FULFNAME(NAMOFILE,PATH,FNAME,IERR)
      IF (IERR.NE.0) GO TO 100

C  CHECK OF CORRECTNESS OF GRID COORDINATES OF TREAT ARRAY PART
      IF(NXE.GT.NX.OR.NXB.LT.1.OR.NXB.GT.NXE) THEN
          IERR=1
          GOTO 103
      END IF
      IF(NYE.GT.NY.OR.NYB.LT.1.OR.NYB.GT.NYE) THEN
          IERR=2
          GOTO 103
      END IF
      IF(NZE.GT.NZ.OR.NZB.LT.1.OR.NZB.GT.NZE) THEN
          IERR=3
          GOTO 103
      END IF
C      WRITE(*,'(2X,A)')'OPEN DIRECT FILE FOR READING:'
C      WRITE(*,'(2X,A)')  NAMOFILE
C      WRITE(*,'(2X,A,I4,A,I4)')'RECL LENTH IN WORDS:', NX,' X',NY
      OPEN(40,FILE=NAMOFILE,STATUS='OLD',ACCESS='DIRECT',
     &  FORM='UNFORMATTED',RECL=(NXE-NXB+1)*(NYE-NYB+1)*LRECL*2,ERR=101)

      NRECF=(NFILD-1)*(NZE-NZB+1)     !INITIAL NUMBER OF RECORD

      L=0
      KR=0
      DO K=NZB,NZE
      KR=KR+1
      READ(40,REC=NRECF+KR,ERR=102) ((FILD(I,J,K),I=NXB,NXE),J=NYB,NYE)

C  FULLING UNDEFINITE POINTS BY ZERO INSTED UNDEF
          DO J=NYB,NYE
            DO I=NXB,NXE

            IF(ABS(LU(I,J)).LT.0.5) THEN
               L=L+1
               FILD(I,J,K)=0.0
            END IF

            ENDDO
         ENDDO

      END DO

      CLOSE(40)

      L=(NXE-NXB+1)*(NYE-NYB+1)-L/KR
      WRITE(*,'(1X,A,A)')
     &         'INPUT DATA from ',NAMOFILE(1:LEN_TRIM (NAMOFILE))
      WRITE(*,'(7X,A,I7,A,I7,A,I8,A)')
     &         'DIMENSION OF FIELD =',NXE-NXB+1,' *',NYE-NYB+1,
     &        ' (',L,'-OCEAN POINTS)'
C     WRITE(*,'(2X,A)')'CLOSE DIRECT FILE:'
C     WRITE(*,'(2X,A)')  NAMOFILE
      RETURN
100   WRITE(*,'(2X,A)')'ERROR IN FULL NAME OF FILE FOR READING: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      STOP
101   WRITE(*,'(2X,A)')'ERROR IN OPEN FILE FOR READING: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      STOP
102   WRITE(*,'(2X,A)')'ERROR IN READING FROM FILE: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      WRITE(*,'(18H ERROR IN READING ,I3,6H LEVEL)') K
      STOP
103   WRITE(*,'(2X,A)')'ERROR IN READING FROM FILE: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      WRITE(*,'(2X,A,I3,A)')'ERROR IN GRID DIAPASON OF ',
     &                       IERR,' - COORDINATE'
      STOP
      END

C======================================================================

      SUBROUTINE WDSTD8(PATH,FNAME,NFILD,FILD,LU,NX,NY,NZ,
     &   NXB,NXE,NYB,NYE,NZB,NZE,IERR)
	IMPLICIT NONE
C  NX,NY,NZ - GENERAL DIMESION OF FILD
C  NXB,NXE,NYB,NYE,NZB,NZE - GRID COORDINATES OF TREAT ARRAY SUBDOMAIN
C      WHERE INDEX B DENOTES BEGIN, AND E - END

C     This subroutine fills (WRITE) array FILD to unformatted deirect
C     file of DIOGIN standard
C     PATH           - path to file (i.g. 'F:\ARAB')
C     FNAME          - name of file (i.g.: 'taux.std')
C     NFILD          - number of field in file (on t)
C     FILD(NX,NY,NZ) - field array
C     LU(NX,NY)    - OCEAN MASK
C---------------------------------------------------------------------
      INCLUDE '1LREC.INC'       !SET LONG OF UNIQUE RECL
      CHARACTER*(*) PATH, FNAME
      INTEGER       NFILD,NRECF,NX,NY,NZ,I,J,K,IERR
      REAL*8        FILD(NX,NY,NZ)
      REAL            LU(NX,NY)
      CHARACTER(4096) NAMOFILE
      INTEGER      NXE, NXB, NYE, NYB, NZB, NZE, L, KR

C  DEFINITION FULL FILE NAME
      CALL FULFNAME(NAMOFILE,PATH,FNAME,IERR)
      IF (IERR.NE.0) GO TO 100

C  CHECK OF CORRECTNESS OF GRID COORDINATES OF TREAT ARRAY PART
      IF(NXE.GT.NX.OR.NXB.LT.1.OR.NXB.GT.NXE) THEN
          IERR=1
          GOTO 103
      END IF
      IF(NYE.GT.NY.OR.NYB.LT.1.OR.NYB.GT.NYE) THEN
          IERR=2
          GOTO 103
      END IF
      IF(NZE.GT.NZ.OR.NZB.LT.1.OR.NZB.GT.NZE) THEN
          IERR=3
          GOTO 103
      END IF
C      WRITE(*,'(2X,A)')'OPEN DIRECT FILE FOR WRITING:'
C      WRITE(*,'(2X,A)')  NAMOFILE
C      WRITE(*,'(2X,A,I4,A,I4)')'RECL LENTH IN WORDS:', NX,' X',NY
      OPEN(40,FILE=NAMOFILE,STATUS='UNKNOWN',ACCESS='DIRECT',
     &  FORM='UNFORMATTED',RECL=(NXE-NXB+1)*(NYE-NYB+1)*LRECL*2,ERR=101)

      NRECF=(NFILD-1)*(NZE-NZB+1)     !INITIAL NUMBER OF RECORD

      L   =0
      IERR=0
      KR=0
      DO K=NZB,NZE
         KR=KR+1
C  FULLING UNDEFINITE POINTS BY 0VER INSTED ZERO
          DO J=NYB,NYE
            DO I=NXB,NXE

            IF(ABS(LU(I,J)).LT.0.5) THEN
               L=L+1
               FILD(I,J,K)=DBLE(UNDEF)
            END IF

            ENDDO
         ENDDO

C  WRITING ON THE FILE
       WRITE(40,REC=NRECF+KR,ERR=102)((FILD(I,J,K),I=NXB,NXE),J=NYB,NYE)

C  FULLING UNDEFINITE POINTS BY ZERO INSTED UNDEF
          DO J=NYB,NYE
            DO I=NXB,NXE

            IF(ABS(LU(I,J)).LT.0.5) THEN
               IERR=IERR+1
               FILD(I,J,K)=0.0
            END IF

            ENDDO
         ENDDO
      END DO


      L=(NXE-NXB+1)*(NYE-NYB+1)-L/KR
      WRITE(*,'(1X,A,A)')
     &         'OUTPUT DATA to ',NAMOFILE(1:LEN_TRIM (NAMOFILE))
      WRITE(*,'(8X,A,I7,A,I7,A,I8,A)')
     &         'DIMENSION OF FIELD =',NXE-NXB+1,' *',NYE-NYB+1,
     &        ' (',L,'-OCEAN POINTS)'
      CLOSE(40)
C      WRITE(*,'(2X,A)')'CLOSE DIRECT FILE:'
C      WRITE(*,'(2X,A)') NAMOFILE

      IERR=(NXE-NXB+1)*(NYE-NYB+1)-IERR/KR-L
      IF (IERR.NE.0) THEN
            WRITE(*,'(2X,A)')  NAMOFILE
            WRITE(*,'(I7,A)') IERR,
     &      'ERRORS IN NUMBER OF OCEAN HORIZONTAL GRID POINTS.'
      ENDIF

      RETURN
100   WRITE(*,'(2X,A)')'ERROR IN FULL NAME OF FILE FOR WRITING: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      STOP
101   WRITE(*,'(2X,A)')'ERROR IN OPEN FILE FOR WRITING: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      STOP
102   WRITE(*,'(2X,A)')'ERROR IN WRITING ON FILE: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      WRITE(*,'(18H ERROR IN WRITING ,I3,6H LEVEL)') K
      STOP
103   WRITE(*,'(2X,A)')'ERROR IN WRITING TO FILE: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      WRITE(*,'(2X,A,I3,A)')'ERROR IN GRID DIAPASON OF ',
     &                       IERR,' - COORDINATE'
      STOP
      END
C======================================================================
      SUBROUTINE RDLUSTD(PATH,FNAME,NFILD,FILD,LU,NX,NY,NZ,
     &   NXB,NXE,NYB,NYE,NZB,NZE,IERR)
	IMPLICIT NONE
C  NX,NY,NZ - GENERAL DIMESION OF FILD
C  NXB,NXE,NYB,NYE,NZB,NZE - GRID COORDINATES OF TREAT ARRAY SUBDOMAIN
C      WHERE INDEX B DENOTES BEGIN, AND E - END
C     This subroutine fills (READ) array FILD from unformatted deirect
C     file of DIOGIN standard
C     PATH           - path to file (i.g. 'F:\ARAB')
C     FNAME          - name of file (i.g.: 'taux.std')
C     NFILD          - number of field in file (on t)
C     FILD(NX,NY,NZ) - field array
C     LU(NX,NY)    - OCEAN MASK
C---------------------------------------------------------------------
      INCLUDE '1LREC.INC'       !SET LONG OF UNIQUE RECL
      INTEGER   LRECMAX
      PARAMETER(LRECMAX=131072)
      CHARACTER*(*) PATH, FNAME
      INTEGER       NFILD,NRECF,NX,NY,NZ,I,J,K,L,LN,IERR
      REAL          FILD(NX,NY,NZ), LU(NX,NY)
      REAL*4    BUFFER(LRECMAX)             !BUFFER FOR WRITING
      CHARACTER(4096) NAMOFILE
      INTEGER  NXE, NXB, NYE, NYB, NZB, NZE, KR


C  DEFINITION FULL FILE NAME
      CALL FULFNAME(NAMOFILE,PATH,FNAME,IERR)
      IF (IERR.NE.0) GO TO 100

C  CHECK OF CORRECTNESS OF GRID COORDINATES OF TREAT ARRAY PART
      IF(NXE.GT.NX.OR.NXB.LT.1.OR.NXB.GT.NXE) THEN
          IERR=1
          GOTO 103
      END IF
      IF(NYE.GT.NY.OR.NYB.LT.1.OR.NYB.GT.NYE) THEN
          IERR=2
          GOTO 103
      END IF
      IF(NZE.GT.NZ.OR.NZB.LT.1.OR.NZB.GT.NZE) THEN
          IERR=3
          GOTO 103
      END IF

C  DEFINITION NUMBER OF DEFINED POINTS
      LN=0
      DO J=NYB,NYE
         DO I=NXB,NXE
         IF(ABS(LU(I,J)).GT.0.5) LN=LN+1
         ENDDO
      ENDDO

      IF(LN.EQ.0) THEN
        IERR=1
        WRITE(*,'(2X,A)')'ERROR IN RDLUSTD!'
        WRITE(*,'(2X,A)')'ZERO OCEAN HORIZONTAL GRID POINTS!'
        STOP 1
      END IF
      IF(LN.GT.LRECMAX) THEN
        IERR=1
        WRITE(*,'(2X,A)')'ERROR IN RDLUSTD!'
        WRITE(*,'(2X,A,I7)')'OCEAN HORIZONTAL GRID POINTS',LN
        WRITE(*,'(2X,A,I7)')'IS GREATER THEN LONG OF BUFFER!',LRECMAX
        WRITE(*,'(2X,A)')
     &          'ENLARGE THIS PARAMETER AND RETRANSLATE PROGRAMM'
        STOP 1
      END IF

C      WRITE(*,'(2X,A)') 'OPEN DIRECT FILE FOR READING:'
C      WRITE(*,'(2X,A)')  NAMOFILE
C      WRITE(*,'(2X,A,I7)')'RECL LENTH IN WORDS:', LN
C      WRITE(*,'(2X,A)') '(CALCULATABLE HORIZONTAL GRID POINTS)'
      WRITE(*,'(2X,A,A32,A,I6)')
     &     'INPUT from:',NAMOFILE(1:LEN_TRIM (NAMOFILE)),'; RLW=', LN

      OPEN(40,FILE=NAMOFILE,STATUS='OLD',ACCESS='DIRECT',
     &        FORM='UNFORMATTED',RECL=LN*LRECL,ERR=101)

      NRECF=(NFILD-1)*(NZE-NZB+1)     !INITIAL NUMBER OF RECORD

      KR=0
      DO K=NZB,NZE
         KR=KR+1
         READ(40,REC=NRECF+KR,ERR=102) (BUFFER(L),L=1,LN)
         L=0
         DO J=NYB,NYE
            DO I=NXB,NXE

            IF(ABS(LU(I,J)).GT.0.5) THEN
               L=L+1
               FILD(I,J,K)=BUFFER(L)
            ELSE
               FILD(I,J,K)=0.0
            END IF

            ENDDO
         ENDDO

         IF(L.NE.LN) THEN
            IERR=1
            WRITE(*,'(2X,A)')   'ERROR IN READING FROM FILE:'
            WRITE(*,'(2X,A)')    NAMOFILE
            WRITE(*,'(2X,A,I6)')'LEVEL:',K
            STOP 1
         END IF

      END DO

      CLOSE(40)
C      WRITE(*,'(2X,A)') 'CLOSE DIRECT FILE:'
C      WRITE(*,'(2X,A)')  NAMOFILE

      RETURN
100   WRITE(*,'(2X,A)')'ERROR IN FULL NAME OF FILE FOR READING: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      STOP
101   WRITE(*,'(2X,A)')'ERROR IN OPEN FILE FOR READING: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      STOP
102   WRITE(*,'(2X,A)')'ERROR IN READING FROM FILE: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      WRITE(*,'(18H ERROR IN READING ,I3,6H LEVEL)') K
      STOP
103   WRITE(*,'(2X,A)')'ERROR IN READING FROM FILE: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      WRITE(*,'(2X,A,I3,A)')'ERROR IN GRID DIAPASON OF ',
     &                       IERR,' - COORDINATE'
      STOP
      END

C======================================================================

      SUBROUTINE WDLUSTD(PATH,FNAME,NFILD,FILD,LU,NX,NY,NZ,
     &   NXB,NXE,NYB,NYE,NZB,NZE,IERR)
	IMPLICIT NONE
C  NX,NY,NZ - GENERAL DIMESION OF FILD
C  NXB,NXE,NYB,NYE,NZB,NZE - GRID COORDINATES OF TREAT ARRAY SUBDOMAIN
C      WHERE INDEX B DENOTES BEGIN, AND E - END
C     This subroutine fills (WRITE) array FILD to unformatted deirect
C     file of DIOGIN standard
C     PATH           - path to file (i.g. 'F:\ARAB')
C     FNAME          - name of file (i.g.: 'taux.std')
C     NFILD          - number of field in file (on t)
C     FILD(NX,NY,NZ) - field array
C     LU(NX,NY)    - OCEAN MASK
C---------------------------------------------------------------------
      INCLUDE '1LREC.INC'       !SET LONG OF UNIQUE RECL
      INTEGER   LRECMAX
      PARAMETER(LRECMAX=131072)
      CHARACTER*(*) PATH, FNAME
      INTEGER   NFILD,NRECF,NX,NY,NZ,I,J,K,L,LN,IERR
      REAL      FILD(NX,NY,NZ), LU(NX,NY)
      REAL*4    BUFFER(LRECMAX)  !BUFFER FOR WRITING
      CHARACTER(4096) NAMOFILE
      INTEGER  NXE, NXB, NYE, NYB, NZB, NZE, KR

C  DEFINITION FULL FILE NAME
      CALL FULFNAME(NAMOFILE,PATH,FNAME,IERR)
      IF (IERR.NE.0) GO TO 100

C  CHECK OF CORRECTNESS OF GRID COORDINATES OF TREAT ARRAY PART
      IF(NXE.GT.NX.OR.NXB.LT.1.OR.NXB.GT.NXE) THEN
          IERR=1
          GOTO 103
      END IF
      IF(NYE.GT.NY.OR.NYB.LT.1.OR.NYB.GT.NYE) THEN
          IERR=2
          GOTO 103
      END IF
      IF(NZE.GT.NZ.OR.NZB.LT.1.OR.NZB.GT.NZE) THEN
          IERR=3
          GOTO 103
      END IF

C  DEFINITION NUMBER OF MARCETABLE POINTS
      LN=0
      DO J=NYB,NYE
         DO I=NXB,NXE
         IF(ABS(LU(I,J)).GT.0.5) LN=LN+1
         ENDDO
      ENDDO

      IF(LN.EQ.0) THEN
         IERR=1
         WRITE(*,'(2X,A)') 'ERROR IN WDLUSTD!'
         WRITE(*,'(2X,A)') 'NUMBER OF OCEAN POINTS IS ZERO!'
         STOP 1
      END IF
      IF(LN.GT.LRECMAX) THEN
        IERR=1
        WRITE(*,'(2X,A)')  'ERROR IN WDLUSTD!'
        WRITE(*,'(2X,A,I7)')'OCEAN HORIZONTAL GRID POINTS',LN
        WRITE(*,'(2X,A,I7)')'IS GREATER THEN LONG OF BUFFER!',LRECMAX
        WRITE(*,'(2X,A)')
     &          'ENLARGE THIS PARAMETER AND RETRANSLATE PROGRAMM'
        STOP 1
      END IF

C      WRITE(*,'(2X,A)')  'OPEN DIRECT FILE FOR WRITING:'
C      WRITE(*,'(2X,A)')   NAMOFILE
C      WRITE(*,'(2X,A,I4)')'RECL LENTH IN WORDS:', LN
C      WRITE(*,'(2X,A)')'(CALCULATABLE HORIZONTAL GRID POINTS)'
      WRITE(*,'(2X,A,A32,A,I6)')
     &     'OUTPUT to:',NAMOFILE(1:LEN_TRIM (NAMOFILE)),'; RLW=', LN

      OPEN(40,FILE=NAMOFILE,STATUS='UNKNOWN',ACCESS='DIRECT',
     &        FORM='UNFORMATTED',RECL=LN*LRECL,ERR=101)

      NRECF=(NFILD-1)*(NZE-NZB+1)     !INITIAL NUMBER OF RECORD

      KR=0
      DO K=NZB,NZE
         KR=KR+1
         L=0
         DO J=NYB,NYE
            DO I=NXB,NXE

            IF(ABS(LU(I,J)).GT.0.5) THEN
               L=L+1
               BUFFER(L)=FILD(I,J,K)
            END IF

            ENDDO
         ENDDO

         IF(L.NE.LN) THEN
            IERR=1
            WRITE(*,'(2X,A)')'ERROR IN WRITING TO FILE:'
            WRITE(*,'(2X,A)')  NAMOFILE(1:LEN_TRIM(NAMOFILE))
            WRITE(*,'(2X,A,I6,A,I6)')'LEVEL:',K,'; OCEAN GRID POINTS:',L
            STOP 1
         END IF

         WRITE(40,REC=NRECF+KR,ERR=102) (BUFFER(L),L=1,LN)
      END DO

      CLOSE(40)
C      WRITE(*,'(2X,A)')'CLOSE DIRECT FILE:'
C      WRITE(*,'(2X,A)')  NAMOFILE

      RETURN
100   WRITE(*,'(2X,A)')'ERROR IN FULL NAME OF FILE FOR WRITING: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      STOP
101   WRITE(*,'(2X,A)')'ERROR IN OPEN FILE FOR WRITING: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      STOP
102   WRITE(*,'(2X,A)')'ERROR IN WRITING ON FILE: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      WRITE(*,'(18H ERROR IN WRITING ,I3,6H LEVEL)') K
      STOP
103   WRITE(*,'(2X,A)')'ERROR IN WRITING TO FILE: '
      WRITE(*,'(2X,A)') NAMOFILE(1:LEN_TRIM(NAMOFILE))
      WRITE(*,'(2X,A,I3,A)')'ERROR IN GRID DIAPASON OF ',
     &                       IERR,' - COORDINATE'
      STOP
      END
c=======================================================================
c      This function finds the first lexeme in the string.
c      We call lexeme the character set separated by space or tab.
c      Memory for OUT_STRING must be allocated by caller.
c      Rusakov Noida.
c
      SUBROUTINE GET_FIRST_LEXEME(IN_STRING, OUT_STRING)
      IMPLICIT NONE
      CHARACTER(*) IN_STRING
      CHARACTER(*) OUT_STRING

      !remove leading and trim blanks
      OUT_STRING = ADJUSTL(IN_STRING)
      OUT_STRING = TRIM   (OUT_STRING)

      OUT_STRING = OUT_STRING(1 : INDEX(OUT_STRING, ' '))
      END SUBROUTINE
c======================================================================
