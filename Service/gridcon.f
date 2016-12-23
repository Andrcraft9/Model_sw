C======================================================================
C GRID CONSRUCTION MODULE BY TEMPERATURE MASK.
C DIANSKY N.A.,(dinar@inm.ras.ru)
C------------------------------------------------------20.03.02 19:11
      SUBROUTINE GRIDCON(FTEMASK,MLRDEF,
     &                    NGRIDT,NGRIDH,NGRIDU,NGRIDXC,NGRIDYC)
      IMPLICIT NONE
C SUBROUTIN FOR CONSTRUCTION PASS BOUNDARY, VELOSITY AND BOTTOM MASKS
C USING TEMPERATURE MASK IN DIOGIN STANDART
C ARRAYS HH(NX,NY) AND HHQ(NX,NY) ARE USED AS AUXILIARY ARRAYS
      INCLUDE '0COM.INC'
      INCLUDE '1LREC.INC'       !SET LONG OF UNIQUE RECL
C--------------------------------------------------------------------
      CHARACTER*(*) FTEMASK
      CHARACTER FRMT*16,COMMENT*80
C MLRDEF -- DEFINED SECTIONITY OF OCEAN REGION
C NGRIDT,NGRIDH,NGRIDU,NGRIDXC,NGRIDYC -- DEFINED NUMBER OF OCEAN POINTS
C                                         IN APPROPRIATE MASK
      INTEGER MLRDEF,NGRIDT,NGRIDH,NGRIDU,NGRIDXC,NGRIDYC
C TEMPORARY INTEGER INDEXES
      INTEGER M, N, MLRX, MLRXU, MLRY, MLRYU, IERR

      DATA LRX/NY*0/LRY/NX*0/LRXU/NY*0/LRYU/NX*0/

      WRITE(FRMT,1000) NX
1000  FORMAT('(',I9,'I1)')
C READING DIOGIN MASK FROM:
      OPEN (11,FILE=FTEMASK,STATUS='OLD',RECL=NX*LRECL)
      READ (11,  '(A)') COMMENT(1:MIN(80,NX))
      WRITE(*,'(1X,A)') COMMENT
      DO N=NY,1,-1
         READ(11,FRMT,END=99) (LBASINS(M,N),M=1,NX)
      ENDDO

C CONVERSION INTEGER DIOGIN MASK TO REAL MODEL MASK
      CALL IDM2RMM(LU,LBASINS,NX,NY)

      IF (NBASINS.GT.1) THEN
      READ (11,  '(A)') COMMENT(1:MIN(80,NX))
      WRITE(*,'(1X,A)') COMMENT
      DO N=NY,1,-1
         READ(11,FRMT,END=99) (LBASINS(M,N),M=1,NX)
      ENDDO
      END IF

      CLOSE(11)
C  FORMING MASK FOR DEPTH GRID POINTS
C  FORMING LUH FROM LU, WHICH HAVE LAND NEIBOURS IN LUH.
C CONSTRUCTING ARRAY LUH FOR RELIEF HH.

	   LUH = 0.

         DO N=NNN-1,NN
            DO M=MMM-1,MM
            IF(LU(M  ,N  ).GE.0.5.OR.
     &         LU(M  ,N+1).GE.0.5.OR.
     &         LU(M+1,N  ).GE.0.5.OR.
     &         LU(M+1,N+1).GE.0.5)
     &        LUH(M,N) = 1.0
            ENDDO
         ENDDO
C MODIFICATION LUH FOR FOR PERIODIC CASE.
      IF(MMD.NE.0) THEN
         DO N=NNN-1,NN
            IF(LU(MM   ,N  ).GE.0.5.OR.
     &         LU(MM   ,N+1).GE.0.5.OR.
     &         LU(MMM  ,N  ).GE.0.5.OR.
     &         LU(MMM  ,N+1).GE.0.5)    THEN

               LUH(MMM-1,N) = 1.0
               LUH(MM   ,N) = 1.0
            END IF
         ENDDO
      END IF


C CALCULATION OF NUMBER OF REMARKABLE GRIDS
      NGRIDH=0
      DO N=1,NY
       DO M=1,NX
        IF(LUH(M,N).GT.0.5) NGRIDH=NGRIDH+1
       ENDDO
      ENDDO
      WRITE(*,'(A)')'  BOTTOM MASK BASED ON TEMPERATURE MASK HAS'
      WRITE(*,'(A,I9)')'  NUMBER OF GRIDS =',NGRIDH

C ARRAY LUU FOR VELOSITY GRID POINTS.
C FILLING: LUU=LU.
      DO N=1,NY
       DO M=1,NX
        LUU(M,N) = LU(M,N)
       ENDDO
      ENDDO
C DELETING POINTS IN LUU, WHICH HAVE LAND NEIBOURS IN LU: RIGHT-UP.
         DO N=NNN-1,NN
            DO M=MMM-1,MM
               IF(LU(M  ,N+1).LT.0.5.OR.
     &            LU(M+1,N  ).LT.0.5.OR.
     &            LU(M+1,N+1).LT.0.5)
     &         LUU(M,N) =0.0
            ENDDO
         ENDDO

C CYCLIC CONDITIONS IN X DIRECTION: LUU CORRECTION ONLY
      IF(MMD.NE.0) THEN
      WRITE(*,*)'  CYCLIC CORRECTION OF U-GRID MASK(LUU)'
         CALL LUUCYC(LUU,NX,NY,MMM,MM)
      END IF

C COMPUTING IIX,JJX BY ANALYZING PAIR OF POINTS IN LU LINE.
C WITH EXCLUDING PASS WITH LENGTH < 2
      WRITE(*,*)'  COMPUTING IIX, JJX, LRX  BY LU'
      CALL CONIIX(LU,NX,NY,IIX,JJX,MLR,LRX,MMM,MM,MMD,2)
C COMPUTING IIXU,JJXU BY ANALYZING PAIR OF POINTS IN LUU LINE.
      WRITE(*,*)'  COMPUTING IIXU,JJXU,LRXU BY LUU'
      CALL CONIIX(LUU,NX,NY,IIXU,JJXU,MLR,LRXU,MMM,MM,MMD,1)

C COMPUTING IIY,JJY BY ANALYZING PAIR OF POINTS IN LU LINE.
C EXCLUDING PASS WITH LENGTH < 2
      WRITE(*,*)'  COMPUTING IIY, JJY, LRY  BY LU'
      CALL CONIIY(LU,NX,NY,IIY,JJY,MLR,LRY,2)
C COMPUTING IIYU,JJYU BY ANALYZING PAIR OF POINTS IN LUU LINE.
      WRITE(*,*)'  COMPUTING IIYU,JJYU,LRYU BY LUU'
      CALL CONIIY(LUU,NX,NY,IIYU,JJYU,MLR,LRYU,1)

C FIND MAXIMAL SECTIONITY OF AQUATORY-----------------------------------
       MLRX =0
       MLRXU=0
       MLRY =0
       MLRYU=0
       DO  N=1,NY
         IF(LRX(N) .GT.MLRX ) MLRX = LRX(N)
         IF(LRXU(N).GT.MLRXU) MLRXU=LRXU(N)
       END DO
       DO  M=1,NX
         IF(LRY(M) .GT.MLRY ) MLRY = LRY(M)
         IF(LRYU(M).GT.MLRYU) MLRYU=LRYU(M)
       END DO
       MLRDEF=MAX(MLRX,MLRXU,MLRY,MLRYU)

      WRITE(*,'(A,I7,A,I7)')'   MLR=',MLR,';  MLRDEF =',MLRDEF
      WRITE(*,'(A,I7,A,I7)')'   LRX  =',MLRX, ';  LRY =',MLRY
      WRITE(*,'(A,I7,A,I7)')'   LRXU =',MLRXU,';  LRYU=',MLRYU

      IF(MLR.LT.MLRDEF) THEN
      WRITE(*,'(A)')
     & '   ERROR! PARAMETER OF SECTION IN COMMON BLOCK 0COM.INC'
      WRITE(*,'(A,I7,A,I7)')
     &  '   IS LESS THEN DEFINED: MLR=',MLR,' < MLRDEF=',MLRDEF
      WRITE(*,'(A)')'   REMOVE THIS ERROR!'
      STOP 1
      END IF

      IF(MLR.GT.MLRDEF) THEN
      WRITE(*,'(A)')
     &  '   WARNING! PARAMETER OF SECTION IN COMMON BLOCK 0COM.INC'
      WRITE(*,'(A,I7,A,I7)')
     &  '   IS GREATER THEN DEFINED: MLR=',MLR,' > MLRDEF=',MLRDEF
C      WRITE(*,'(A)')'   DO YOU REMOVE THIS WARNING?(1-CONTINUE,0-STOP)'
C      READ (*,*) NCONT
C      IF(NCONT.EQ.0) STOP 2
      END IF

C TESTING ARRAYS OF BOUNDARIES FOR LENGTH.
      WRITE(*,*)
      CALL IITEST(IIX ,JJX ,MLR,NY,LRX ,2,'IIX-JJX ' )
      WRITE(*,*)
      CALL IITEST(IIY ,JJY ,MLR,NX,LRY ,2,'IIY-JJY ' )
      WRITE(*,*)
      CALL IITEST(IIX ,JJX ,MLR,NY,LRX ,3,'IIX-JJX ' )
      WRITE(*,*)
      CALL IITEST(IIY ,JJY ,MLR,NX,LRY ,3,'IIY-JJY ' )
      WRITE(*,*)
      CALL IITEST(IIXU,JJXU,MLR,NY,LRXU,2,'IIXU-JJXU')
      WRITE(*,*)
      CALL IITEST(IIYU,JJYU,MLR,NX,LRYU,2,'IIYU-JJYU')

C TESTING IIX, JJX BY FULLING LU.
      WRITE(*,*)
      WRITE(*,*)'  RECONSTRUCTION OF MASK BASED ON IIX AND JJX'
      CALL RECONX(HH,LU,IIX,JJX,NX,NY,MLR,LRX,NGRIDT,
     &            MM,MMD,IERR)
      WRITE(*,'(A,I9)')
     &  '   MASK BASED ON IIX, JJX; NUMBER OF GRIDS =',NGRIDT

             IF(IERR.NE.0) THEN
      WRITE(*,'(A,I7)')'   NUMBER OF NOT COINCIEDED GRIDS =',IERR
             ENDIF

C TESTING IIY, JJY BY FULLING LU.
      WRITE(*,*)
      WRITE(*,*)'  RECONSTRUCTION OF MASK BASED ON IIY AND JJY'
      CALL RECONY(HHQ,LU,IIY,JJY,NX,NY,MLR,LRY,NGRIDT,IERR)
      WRITE(*,'(A,I9)')
     &     '   MASK BASED ON IIY, JJY; NUMBER OF GRIDS =',NGRIDT

             IF(IERR.NE.0) THEN
      WRITE(*,'(A,I7)')'   NUMBER OF NOT COINCIEDED GRIDS =',IERR
             ENDIF
C TESTING COMMON TEMPERATURE MASK ARRAY FROM X AND Y DIRECTION
      WRITE(*,*)
      WRITE(*,*)'  RECONSTRUCTION OF TEMPERATURE MASK BASED ON'
      WRITE(*,*)'  (IIX,JJX) AND (IIY,JJY) SIMULTANEOUSLY:'
      NGRIDT=0
      IERR=0
      DO N=1,NY
        DO M=1,NX
           IF(HH(M,N).GT.0.5.OR.HHQ(M,N).GT.0.5) HH(M,N)=1.0
           IF(HH(M,N).NE.LU(M,N)) IERR = IERR+1
           IF(HH(M,N).GT.0.5) NGRIDT=NGRIDT+1
        ENDDO
      ENDDO
           IF(IERR.NE.0) THEN
      WRITE(*,*)
     &     '  FOR TEMPERATURE MASK BASED ON (IIX,JJX) AND (IIY,JJY)'
      WRITE(*,*)'  SIMULTANEOUSLY:'
      WRITE(*,'(A,I7)')'   NUMBER OF NOT COINCIEDED GRIDS =',IERR
           ELSE
      WRITE(*,*)'  TEMPERATURE MASK BASED ON (IIX,JJX) AND (IIY,JJY)'
      WRITE(*,*)'  IS EQUAL TO ORIGINAL MASK'
      WRITE(*,'(A,I9)')'   NUMBER OF GRIDS =',NGRIDT
           ENDIF

C TESTING IIXU, JJXU BY FULLING LU.
      WRITE(*,*)
      WRITE(*,*)'  RECONSTRUCTION OF MASK BASED ON IIXU AND JJXU'
      CALL RECONX(HH,LUU,IIXU,JJXU,NX,NY,MLR,LRXU,NGRIDU,
     &            MM,MMD,IERR)
      WRITE(*,'(A,I9)')
     &   '   MASK BASED ON IIXU, JJXU; NUMBER OF GRIDS=',NGRIDU

             IF(IERR.NE.0) THEN
      WRITE(*,'(A,I7)')'   NUMBER OF NOT COINCIEDED GRIDS =',IERR
      WRITE(*,*)'  CHECK MASK ARRAY AND PASS BOUNDARIES!'
             ENDIF

C TESTING IIYU, JJYU BY FULLING LU.
      WRITE(*,*)
      WRITE(*,*)'  RECONSTRUCTION OF MASK BASED ON IIYU AND JJYU'
      CALL RECONY(HHQ,LUU,IIYU,JJYU,NX,NY,MLR,LRYU,NGRIDU,IERR)
      WRITE(*,'(A,I9)')
     &    '   MASK BASED ON IIYU, JJYU; NUMBER OF GRIDS=',NGRIDU

             IF(IERR.NE.0) THEN
      WRITE(*,'(A,I7)')'   NUMBER OF NOT COINCIEDED GRIDS =',IERR
      WRITE(*,*)'  YOU MUST CHECK MASK ARRAY AND PASS BOUNDARIES!'
             ENDIF
             IF (MMD.GT.0) THEN
      WRITE(*,*)'  SET PERIODICITY TO T-GRID MASK(LU ).'
      CALL CYCLIZE(LU ,NX,NY,1,MMM,MM)
      WRITE(*,*)'  SET PERIODICITY TO U-GRID MASK(LUU).'
      CALL CYCLIZE(LUU,NX,NY,1,MMM,MM)
             ENDIF
C C-GRID MASK CONSRUCTION
C LXC,LYC -- C-GRID MASK (ARAKAWA CLASSIFICATION);
C LU      -- B-GRID TEMPERATURE MASK.
      WRITE(*,*)'  CONSTRUCTION OF C-X-Y-GRID MASKS BASED ON T-GRID:'
      CALL CGRIDMASK(LCU,LCV,LLU,LLV,LU,NX,NY,NGRIDXC,NGRIDYC)
             IF (MMD.GT.0) THEN
      WRITE(*,*)'  SET PERIODICITY TO U-GRID MASK(LCU).'
      CALL CYCLIZE(LCU,NX,NY,1,MMM,MM)
      CALL CYCLIZE(LLU,NX,NY,1,MMM,MM)
      WRITE(*,*)'  SET PERIODICITY TO V-GRID MASK(LCV).'
      CALL CYCLIZE(LCV,NX,NY,1,MMM,MM)
      CALL CYCLIZE(LLV,NX,NY,1,MMM,MM)
             ENDIF

      WRITE(*,'(A,I9)')
     &    '   C-GRID MASK IN X DIRECTION HAS NUMBER OF OCEAN POINTS=',
     &                    NGRIDXC
      WRITE(*,'(A,I9)')
     &    '   C-GRID MASK IN Y DIRECTION HAS NUMBER OF OCEAN POINTS=',
     &                    NGRIDYC

      RETURN
99    WRITE(*,*)'  ERROR IN READING FILE ',FTEMASK(1:LEN_TRIM(FTEMASK))
      STOP 1
      END
C======================================================16.03.99 17:19
      SUBROUTINE CGRIDMASK(LCU,LCV,LLU,LLV,LU,NX,NY,NGRIDXC,NGRIDYC)
      IMPLICIT NONE
C C-GRID MASK CONSRUCTION

C LU      -- C-GRID TEMPERATURE MASK.
      INTEGER NX, NY
      REAL LCU(NX,NY),LCV(NX,NY), !MASKS OF U- V- GRIDS BASED ON T-GRID WITH 0 ON BOUNDARY
     &     LLU(NX,NY),LLV(NX,NY), !MASKS OF U- V- GRIDS BASED ON T-GRID WITH 1 ON BOUNDARY
     &     LU(NX,NY)              !MASK ON T-GRID


C NUMBER OF CALCULATABLE POINTES IN APPROPRIATE MASKS IN X AND Y DIRECTION
      INTEGER NGRIDXC,NGRIDYC
      INTEGER M, N

	LCU=0.0
	LCV=0.0
	LLU=0.0
	LLV=0.0

       DO N=1,NY-1
          DO M=1,NX-1
          LCU(M,N)=ANINT(LU(M+1,N)*LU(M,N))
          LCV(M,N)=ANINT(LU(M,N+1)*LU(M,N))
          LLU(M,N)=MIN(1.0,LU(M+1,N)+LU(M,N))
          LLV(M,N)=MIN(1.0,LU(M,N+1)+LU(M,N))         
          ENDDO
       ENDDO

       NGRIDXC=0
       NGRIDYC=0
       DO N=1,NY-1
          DO M=1,NX-1
          IF(LCU(M,N).GT.0.5) NGRIDXC=NGRIDXC+1
          IF(LCV(M,N).GT.0.5) NGRIDYC=NGRIDYC+1
          ENDDO
       ENDDO
      RETURN
      END
C======================================================18.10.98 19:13
      SUBROUTINE IDM2RMM(RMASK,IMASK,NX,NY)
      IMPLICIT NONE
C CONVERSION INTEGER DIOGIN MASK TO REAL MODEL MASK
      INTEGER NX, NY
      INTEGER M, N
      REAL RMASK(NX,NY)
      INTEGER IMASK(NX,NY)
      DO N=1,NY
         DO M=1,NX
                             RMASK(M,N)=0.0
         IF(IMASK(M,N).EQ.0) RMASK(M,N)=1.0
         END DO
      END DO
      RETURN
      END
C======================================================18.10.98 19:17
      SUBROUTINE RMM2IDM(RMASK,IMASK,NX,NY)
      IMPLICIT NONE
C CONVERSION REAL MODEL MASK TO INTEGER DIOGIN MASK
      INTEGER NX, NY
      INTEGER M, N
      REAL    RMASK(NX,NY)
      INTEGER IMASK(NX,NY)
      DO N=1,NY
         DO M=1,NX
                               IMASK(M,N)=1
         IF(RMASK(M,N).GT.0.5) IMASK(M,N)=0
         END DO
      END DO
      RETURN
      END
C======================================================18.10.98 19:27
      SUBROUTINE RECONX(LU,LU0,IIX,JJX,NX,NY,MLR,LRX,NGRID,
     &                  MM,MMD,IERR)
      IMPLICIT NONE
C TESTING IIX, JJX BY FULLING LU.
C LU  - MASK ARRAY CONSTRUCTED BY IIX AND JJY
C LU0 - MASK ARRAY CONSTRUCTED IN MAIN PROGRAM(REFERENCE ARRAY)
      INTEGER NX, NY, MLR, MM, MMD
      INTEGER M, N, II, JJ, M1, IGAP, NGRID
      REAL LU(NX,NY),LU0(NX,NY)
      INTEGER LRX(NY),IIX(MLR,NY),JJX(MLR,NY),IERR

      IERR=0
      NGRID=0
      DO N=1,NY
       DO M=1,NX
        LU (M,N) =0.0
       ENDDO
      ENDDO
      DO N=1,NY
       DO IGAP=1,IABS(LRX(N))
        II=IIX(IGAP,N)
        JJ=JJX(IGAP,N)
C CYCLIC CONDITION
        IF(II.NE.0)THEN
         DO M=II,JJ
            M1=M
C CYCLIC CONDITION
            IF(M1.GT.MM) M1=M1-MMD
            LU(M1,N) =1.0
            NGRID=NGRID+1
         ENDDO
        ENDIF
       ENDDO
      ENDDO
      DO N=1,NY
       DO M=1,NX
        IF(LU(M,N).NE.LU0(M,N)) IERR = IERR+1
       ENDDO
      ENDDO

       IF(IERR.NE.0) THEN
        WRITE(*,*)
     &    '!!!ERROR IN CODE ARRAY RECONSTRACTING IN X-DIRECTION!!!'
       ENDIF

      RETURN
      END
C======================================================18.10.98 19:28
      SUBROUTINE RECONY(LU,LU0,IIY,JJY,NX,NY,MLR,LRY,NGRID,IERR)
      IMPLICIT NONE
C TESTING IIY, JJY BY FULLING LU.
C LU  - MASK ARRAY CONSTRUCTED BY IIY AND JJY
C LU0 - MASK ARRAY CONSTRUCTED IN MAIN PROGRAM(REFERENCE ARRAY)
      INTEGER NX, NY, MLR
      INTEGER M, N, IGAP, II, JJ, NGRID
      REAL LU(NX,NY),LU0(NX,NY)
      INTEGER LRY(NX),IIY(MLR,NX),JJY(MLR,NX),IERR
      IERR=0
      NGRID=0
      DO N=1,NY
       DO M=1,NX
        LU (M,N) =0.0
       ENDDO
      ENDDO
      DO M=1,NX
       DO IGAP=1,LRY(M)
        II=IIY(IGAP,M)
        JJ=JJY(IGAP,M)
        IF(II.NE.0)THEN
         DO N=II,JJ
          LU(M,N) =1.0
          NGRID=NGRID+1
         ENDDO
        ENDIF
       ENDDO
      ENDDO
      DO N=1,NY
       DO M=1,NX
        IF(LU(M,N).NE.LU0(M,N)) IERR = IERR+1
       ENDDO
      ENDDO

       IF(IERR.NE.0) THEN
        WRITE(*,*)
     &   '!!!ERROR IN CODE ARRAY RECONSTRACTING IN Y-DIRECTION!!!'
       ENDIF

      RETURN
      END
C======================================================21.02.07 16:30
      SUBROUTINE CONIIX(LU,NX,NY,IIX,JJX,MLR,LRX,MMM,MM,MMD,MINPLEN)
	IMPLICIT NONE
C--------------------------------------------------------------------
C CONSTRUCTING ARRAYS IIX,JJX OF BEGINNINGS AND ENDS OF X-LINES
C USING ARRAY LU OF 0,1.
C MINPLEN - MINIMAL PASS LENGTH
      INTEGER   NHSECT
      PARAMETER(NHSECT=1024)
      INTEGER NX, NY, MLR, MMM, MM, MMD
      INTEGER M, N, L, LP, MBEGIN, NLRT, LENOFPASS, MINPLEN
      REAL LU(NX,NY)
      INTEGER IIX(MLR,NY),JJX(MLR,NY),LRX(NY),IB(NHSECT,2)

      IF(MLR.GT.NHSECT-1) THEN
      WRITE(*,'(A,I7,A,I7)')
     &'  ENLARGE PARAMETER NHSECT IN SUBROUTINE CONIIX BECAUSE NHSECT='
     &   ,NHSECT,' <  MLR+1=',MLR+1
      STOP 1
      END IF

C----INITIAL VALUES--------------------------------------------------

          LRX =0
          IIX=0
          JJX=0
          LENOFPASS=0          
       DO N=1,NY

C--CALCULATE BOUNDARIES WITHOUT RESTRICTION IN NON PERIODIC CASE FOR N-LINE-----
	  IB=0
        NLRT=0
        DO M=MMM-1,MM+1
C---------FOR BEGIN OF PASS---------------------------
          IF(LU(M-1,N).LT.0.5.AND.LU(M,N).GT.0.5) THEN
              NLRT=NLRT+1
              MBEGIN=M
              LENOFPASS=0
          END IF

          LENOFPASS=LENOFPASS+1

C---------FOR END OF PASS-----------------------------
          IF(LU(M,N).GT.0.5.AND.LU(M+1,N).LT.0.5) THEN
              IB(NLRT,1)=MBEGIN
              IB(NLRT,2)=M
          END IF

        END DO

C  CHECKING FOR CORRECTION OF PASS BOUNDARIES IN CYCLIC CASE
          IF(MMD.NE.0.AND.NLRT.GT.0) THEN

            IF(LU(MMM,N).GT.0.5.AND.LU(MM,N).GT.0.5) THEN
               IF(NLRT.EQ.1) THEN
C  CYCLE LOOP WITHOUT BOUNDARIES
                  NLRT=-1
               ELSE
C  CYCLE LOOP WITH BOUNDARIES
                  IB(1,1)=IB(NLRT,1)
                  IB(1,2)=IB(1,2)+MMD

			    IB(NLRT,:)=0
	            NLRT=NLRT-1

               END IF
            END IF
          END IF


C---SET RESTRICTION FOR N-LINE--------------------------   

        IF(NLRT.GT.0) THEN 
      
		LP=1
				
	    DO WHILE(LP.LE.NLRT)   !BEGIN OPERATION
				
		LENOFPASS=IB(LP,2)-IB(LP,1)+1

            IF(LENOFPASS.LT.MINPLEN) THEN

              WRITE(*,'(A,I7,A,I7,A,I7,A)')
     &        '  LINE N=',N,' MBEGIN=',IB(LP,1),' MEND=',IB(LP,2),
     &                   '  HAS BEEN EXCLUDED!'
              WRITE(*,'(A,I7,A,I7)')
     &        '  BECAUSE ITS LENGTH=',LENOFPASS,' < ',MINPLEN

			NLRT=NLRT-1 
              
			  IF(LP.LE.NLRT) THEN    

C----------SHIFT PASS BOUNDARIES-----------------------              
				DO L=LP,NLRT
				IB(L,1)=IB(L+1,1)
				IB(L,2)=IB(L+1,2)
				END DO

			  END IF 
			  
			  IB(NLRT+1,:)=0					         

	      ELSE 

		    LP=LP+1
		  
		  END IF
		
		END DO                  !END OPERATION 
	   
	   END IF
			 
              IF (NLRT.GT.MLR) THEN
         WRITE(*,*)'  ERROR IN CONIIX:'
         WRITE(*,'(A,I7,A,I7,A,I7)')
     &           '   MLR=',MLR,' <  LRX(',N,')=',NLRT
         WRITE(*,*)'  SET NEW VALUE OF MLR IN 1BASINPAR.INC FOR'
         WRITE(*,*)'  MAIN PROGRAM, THEN RECOMPILE IT!!!'
              STOP 1
              ENDIF


          DO L=1,IABS(NLRT)
            IIX(L,N)=IB(L,1)
            JJX(L,N)=IB(L,2)
          END DO

          LRX(N) = NLRT
       END DO
      RETURN
      END

C========================================================18.10.98 19:29
      SUBROUTINE CONIIY(LU,NX,NY,IIY,JJY,MLR,LRY,MINPLEN)
      IMPLICIT NONE
C----------------------------------------------------------------------
C CONSTRUCTING ARRAYS IIY,JJY OF BEGINNINGS AND ENDS OF Y-LINES
C USING ARRAY LU OF 0,1.
C MINPLEN - MINIMAL PASS LENGTH
      INTEGER NX, NY, MLR
      INTEGER M, N, K, NBEGIN, NLRT
      REAL LU(NX,NY)
      INTEGER IIY(MLR,NX),JJY(MLR,NX),LRY(NX)
	INTEGER LENOFPASS, MINPLEN
C----INITIAL VALUES----------------------------------------------------
       DO M=1,NX
          LRY(M) =0
          DO K=1,MLR
            IIY(K,M)=0
            JJY(K,M)=0
          END DO
       END DO

       LENOFPASS=0
       
       DO M=2,NX-1
C----CALCULATE BOUNDARIES --------------------------------
          NLRT=0
          DO N=2,NY-1

C---------FOR BEGIN OF PASS---------------------------
          IF(LU(M,N-1).LT.0.5.AND.LU(M,N).GT.0.5) THEN
              NBEGIN=N
              NLRT=NLRT+1
              LENOFPASS=0
            IF (NLRT.GT.MLR) THEN
         WRITE(*,*)'  ERROR IN CONIIY:'
         WRITE(*,'(A,I7,A,I7,A,I7)')
     &           '   MLR=',MLR,' <  LRY(',N,')=',NLRT
         WRITE(*,*)'  SET NEW VALUE OF MLR IN 0COM.INC FOR'
         WRITE(*,*)'  MAIN PROGRAM, THEN RECOMPILE IT!!!'
         STOP 1
            ENDIF
          END IF

          LENOFPASS=LENOFPASS+1

C---------FOR END OF PASS-----------------------------
          IF(LU(M,N).GT.0.5.AND.LU(M,N+1).LT.0.5) THEN
            IF(LENOFPASS.LT.MINPLEN) THEN
              WRITE(*,'(A,I7,A,I7,A,I7,A)')
     &        '  LINE M=',M,' NBEGIN=',NBEGIN,' NEND=',N,
     &                   '  HAS BEEN EXCLUDED!'
              WRITE(*,'(A,I7,A,I7)')
     &        '  BECOSE ITS LENGTH=',LENOFPASS,' < ',MINPLEN
              NLRT=NLRT-1
            ELSE
              IIY(NLRT,M)=NBEGIN
              JJY(NLRT,M)=N
            END IF
          END IF

          END DO
          LRY(M) = NLRT
       END DO
      RETURN
      END
C========================================================18.10.98 19:31
      SUBROUTINE IITEST(IIX,JJX,MLR,NY,LRX,NLEN,TEXT)
      IMPLICIT NONE
C----------------------------------------------------------------------
C TESTS IIX,JJX FOR PRESENCE OF ENDS OF LINES AND LENGTHS OF LINES.
      INTEGER NY, MLR
      INTEGER IG, N, II, JJ, IERR
      INTEGER IIX(MLR,NY),JJX(MLR,NY),LRX(NY)
      CHARACTER *(*) TEXT
	INTEGER   LENTEXT, NLEN

      LENTEXT=LEN(TEXT)
      WRITE(*,'(2A)')  '   BEGIN TESTING ',TEXT
      WRITE(*,'(A,I7)')'   SEARCHING FOR LENGTH <',NLEN
C IF (II.NOT EQUAL.0) COMPUTING DIFFERENCE JJ-II.

      IERR=0
      DO N=1,NY
       DO IG=1,IABS(LRX(N))
        II=IIX(IG,N)
        JJ=JJX(IG,N)
        IF (II.NE.0) THEN
         IF ((JJ-II+1).LT.NLEN) THEN
          IERR=IERR+1
          WRITE(*,100) II,JJ,IG,JJ-II+1,TEXT(LENTEXT-1:LENTEXT),N,LRX(N)
          END IF
        ENDIF
       ENDDO
      ENDDO
      WRITE(*,'(A,I7)')'   END TESTING; ERRORS=',IERR
      RETURN
100   FORMAT('  II,JJ=',2I4,'; IG=',I7,'; LENGTH=',I7,
     &           '; LR',A2,'(',I7,')=',I7)
      END
C========================================================18.10.98 19:33
      SUBROUTINE LUUCYC(LUU,NX,NY,MMM,MM)
      IMPLICIT NONE
      INTEGER NX, NY, MMM, MM
      INTEGER N
      REAL LUU(NX,NY)
C  IT IS ASSUMED THAT PERIODIC LENGTH = MMD
      DO N=1,NY
C  CORRECTION RIGHT BOUNDARY OF VELOSITY GRID COD ARRAY LUU
            IF(LUU(MMM,N).GT.0.5.AND.LUU(MM-1,N).GT.0.5) THEN
               LUU(MM,N)=1.0
            END IF
      END DO

      RETURN
      END
C======================================================================
