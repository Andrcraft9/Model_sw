      SUBROUTINE NPSF_COEF_DEF() 
C DEFINITION OF COEFFICIENTS FOR SIMMETRIC FILTRATION OVER NORTH POLE 
C      ALONG LATITUDE CIRCLE

	IMPLICIT NONE
	INCLUDE '0COM.INC'
C      COMMON /NORTH_POLE_FILTER/
C     &        NUM_LAT_NPSF(NY),   !NUMBER OF SIM FILTER COEFFICIENTS
C	                             !ON GIVEN LATITUDE CIRCLE
C     &        NUM_NPSF_LAT(NY),   !NUMBER OF SIM FILTER 
C	                             !ON GIVEN LATITUDE CIRCLE
C     &        CF_OF_FNPS(MAX_LONG_SFLT,NUM_NPSF)    !COEFFICIENTS OF NORTH POLE FILTER


	REAL  CUT_OFF_PERIOD,REAL_STEP

	INTEGER N,NXA_NEW, 
     &        N_NPSF !NUMBER OF LAT CIRCLE

	REAL PIP180
	PARAMETER(PIP180=PI/180.0)

	NXA_NEW=MM-MMM+1

c	CUT_OFF_PERIOD = 
c     &   4.0*DXST * COS(YT(NUM_LAT_BEG_FILTER)*PIP180 )
	CUT_OFF_PERIOD = 
     &   4.0*DXST * COS((RLAT+DYST*(NUM_LAT_BEG_FILTER-NNN))*PIP180 )
     	
	   N_NPSF=0
      
	CF_OF_FNPS=0.0
	CF_OF_FNPS(1,:)=1.0
	
	NUM_LAT_NPSF = 0 
	
	IF(NUM_LAT_BEG_FILTER.GE.NNN.AND.NUM_LAT_BEG_FILTER.LE.NN) THEN

	DO N= NUM_LAT_BEG_FILTER,NN	
         
	   N_NPSF=N_NPSF+1

	   IF(N_NPSF.GT.NUM_NPSF) THEN

	      WRITE(*,*)' ERROR IN SUBROUTINE NPSF_COEF_DEF()!'
	      WRITE(*,*)
     &      ' NUMBERING OF LATITUDE CIRCLE FOR NPSF IS GREATER THEN'
            WRITE(*,*)
     &      ' ITS MAX NUMBER!'
	      WRITE(*,*)' CORRECT 1BASINPAR.INC FOR NORTH POLE FILTER!'
	      STOP ' ERROR IN SUBROUTINE NPSF_COEF_DEF()!'
	   	
         END IF

         NUM_NPSF_LAT(N)  = N_NPSF

C	   REAL_STEP=DXST * COS( YT(N)*PIP180 )
	   REAL_STEP=DXST * COS( (RLAT+DYST*(N-NNN))*PIP180 )

         NUM_LAT_NPSF(N) =    
     &             MIN( MAX( NINT(8.0*CUT_OFF_PERIOD/REAL_STEP) ,32 ), 
     &                                                   NXA_NEW/2-1 )

	   IF(NUM_LAT_NPSF(N).GT.MAX_LONG_SFLT) THEN

	      WRITE(*,*)' ERROR IN SUBROUTINE NPSF_COEF_DEF()!'
	      WRITE(*,*) NUM_LAT_NPSF(N),' > ',MAX_LONG_SFLT
	      WRITE(*,*)
     &      ' NUMBER OF COEFFICIENTS ON LATITUDE =',N,RLAT+DYST*(N-NNN)
	      WRITE(*,*)' IS GREATER THEN MAX_LONG_SFLT !'
	      WRITE(*,*)' CORRECT 1BASINPAR.INC FOR NORTH POLE FILTER!'
	      STOP ' ERROR IN SUBROUTINE NPSF_COEF_DEF()!'
	   	
	   END IF
	
		CALL SIMFILTER_COEF_DEFINITION(CF_OF_FNPS(:,NUM_NPSF_LAT(N)),
     &                                                NUM_LAT_NPSF(N),
     &                                       CUT_OFF_PERIOD,REAL_STEP  )  
      END DO

          IF(N_NPSF.GT.NUM_NPSF) THEN

            WRITE(*,*)' ERROR IN SUBROUTINE NPSF_COEF_DEF()!'
            WRITE(*,*)
     &      ' COMPLETE NUMBER OF LATITUDE CIRCLE FOR NPSF IS NOT EQUAL'
            WRITE(*,*)
     &      ' TO ITS MAX NUMBER!'
            WRITE(*,*)' CORRECT 1BASINPAR.INC FOR NORTH POLE FILTER!'
            STOP ' ERROR IN SUBROUTINE NPSF_COEF_DEF()!'
	   	
           END IF

      END IF

      END	SUBROUTINE
C======================================================================
      SUBROUTINE NPS_FILTER_T_GRD(VARIN,   !3D INPUT ARRAY FOR FILTRATION
     &                              NNZ )  !ITS DIMENSION ON Z   

C SIMMETRIC FILTER REALIZATION ON T-GRID  
	IMPLICIT NONE
	INCLUDE '0COM.INC'

	INTEGER  NNZ   !DIMENSION ON Z-GRID 
 
      REAL  VARIN(NX,NY,NNZ) !INPUT DATA ARRAY
	
	REAL, ALLOCATABLE:: VAR(:)         !HELP ARRAYS FOR IN- AND OUT-PUT FILTER

	INTEGER K,N,NXA_NEW

	NXA_NEW=MM-MMM+1
   
	ALLOCATE (VAR(NXA_NEW))

C!$OMP PARALLEL DO PRIVATE(N,K,VAR)
	  DO K=1,NNZ		
		DO N= NNN,NN	
	
		IF(NUM_LAT_NPSF(N).GT.0) THEN
				
          VAR(1:NXA_NEW)=VARIN(MMM:MM,N,K)
					
	    CALL FILTER_ALONG_LON(VAR,VARIN(MMM,N,K),HHQ(MMM,N),
     &        CF_OF_FNPS(:,NUM_NPSF_LAT(N)),NXA_NEW,NUM_LAT_NPSF(N))

	    END IF        
		
	    END DO      
        END DO
C!$OMP END PARALLEL DO	
	DEALLOCATE (VAR)

      END	SUBROUTINE
C======================================================================
      SUBROUTINE NPS_FILTER_T_GRD_NW(VARIN,   !3D INPUT ARRAY FOR FILTRATION
     &                                NNZ )  !ITS DIMENSION ON Z   

C SIMMETRIC FILTER REALIZATION ON T-GRID  
	IMPLICIT NONE
	INCLUDE '0COM.INC'

	INTEGER  NNZ   !DIMENSION ON Z-GRID 
 
      REAL  VARIN(NX,NY,NNZ) !INPUT DATA ARRAY
	
	REAL, ALLOCATABLE:: VAR(:)         !HELP ARRAYS FOR IN- AND OUT-PUT FILTER

	INTEGER K,N,NXA_NEW

	NXA_NEW=MM-MMM+1
   
	ALLOCATE (VAR(NXA_NEW))

C!$OMP PARALLEL DO PRIVATE(N,K,VAR)
	  DO K=1,NNZ		
		DO N= NNN,NN	
	
		IF(NUM_LAT_NPSF(N).GT.0) THEN
				
          VAR(1:NXA_NEW)=VARIN(MMM:MM,N,K)
					
	    CALL FILTER_ALONG_LON_NW(VAR,VARIN(MMM,N,K),
     &        CF_OF_FNPS(:,NUM_NPSF_LAT(N)),NXA_NEW,NUM_LAT_NPSF(N))

	    END IF        
		
	    END DO      
        END DO
C!$OMP END PARALLEL DO	
	DEALLOCATE (VAR)

      END SUBROUTINE
C======================================================================
      SUBROUTINE SIMFILTER_COEF_DEFINITION(CFF,LOT,TPER,DELT)
	IMPLICIT NONE
C  DEFINITION COEFFICIENTS OF FILTER
      INTEGER M,K,LOTM1,LOT
      REAL CFF(LOT),DA,D0,PI,D(3),SUM,SUMG,TPER,DELT

      DATA D0/0.355770/D/0.2436983,0.07211497,0.00630165/
     *     PI/3.14159265/
C  THIS PROGRAM CALCULATE WEIGHTS OF LOW-FREQUENCES FILTER
C  METHOD SUGEST BY POTTER,BIGFORD & GLAIZ
C  LOT*2-1 NUMBER OF WEIGHTS
C  LOT - NUMBER OF COEFFICIENTS OF WEIGHTS
C  TPER - PERIOD OF BRAKE
C  DELT - DISCRET OF TIME (THE SAME UNITS)
C  CFF - COEFFICIENTS OF WEIGHTS
      LOTM1=LOT-1
      DA=DELT/TPER*2
      CFF(1) = DA
      DA=DA*PI

      DO  K=1,LOTM1
      CFF(K+1)=(SIN(DA*K)/(PI*K))
C	WRITE(*,*) CFF(K+1)
      END DO

C  TRAPECIAY SMOOS ON END
      CFF(LOT)=CFF(LOT)*0.5
C  POTTER'S WINDOW P310
      SUMG=CFF(1)

      DO  K=2,LOT
      SUM=D0
      DA=PI*FLOAT(K-1)/FLOAT(LOTM1)
            DO M=1,3
            SUM=SUM+2*D(M)*COS(DA*M)
            END DO
      CFF(K)=CFF(K)*SUM
      SUMG=SUMG+2*CFF(K)
      END DO
C      WRITE(*,*) 'SUM =',SUMG

            DO K=1,LOT
            CFF(K)=CFF(K)/SUMG
            END DO

      END SUBROUTINE
C======================================================================
      SUBROUTINE FILTER_ALONG_LON(X,Y,W,CFF,NF,LOT)
	IMPLICIT NONE
C X-INPUT,Y-OUTPUT,W-WEIGHTS,CFF-COEFFICIENTS OF FILTER,
C 2*LOT-1=NUMBER OF COEFFICIENTS OF FILTER, 
C NF-NUMBER OF POINTS OF INPUT & OUTPUT SETS
	INTEGER NF,LOT
      REAL X(NF),Y(NF),W(NF)
      REAL CFF(LOT)
      INTEGER I,L,IPE
C
C FILTRATION AT THE BEGINNING OF THE SET WITH THE ASSUMPTION OF ITS PERIODICITY 
C
      DO I=1,LOT-1
            Y(I)=CFF(1)*X(I)*W(I)
            DO L=1,LOT-1
			IPE=I-L
			IF(IPE.LT.1) THEN 
				
				IPE=IPE+NF

			END IF
			Y(I)=Y(I)+CFF(L+1)*(X(IPE)*W(IPE)+X(I+L)*W(I+L))
            END DO
      END DO
C
C FILTRATION IN THE MIDDLE OF THE SET
C
            DO I=LOT,NF-LOT+1
               Y(I)=CFF(1)*X(I)*W(I)
                  DO L=1,LOT-1
               Y(I)=Y(I)+CFF(L+1)*(X(I-L)*W(I-L)+X(I+L)*W(I+L))
                  END DO
            END DO
C
C FILTRATION AT THE END OF THE SET WITH THE ASSUMPTION OF ITS PERIODICITY 
C
      DO I=NF-LOT+2,NF

            Y(I)=CFF(1)*X(I)*W(I)
            DO L=1,LOT-1
			IPE=I+L
              IF(IPE.GT.NF) THEN 
				IPE=IPE-NF
			END IF
              Y(I)=Y(I)+CFF(L+1)*(X(I-L)*W(I-L)+X(IPE)*W(IPE))
            END DO
      END DO

      DO I=1,NF
      Y(I)=Y(I)/W(I)
	END DO

      END	SUBROUTINE
C======================================================================
      SUBROUTINE FILTER_ALONG_LON_NW(X,Y,CFF,NF,LOT)
	IMPLICIT NONE
C X-INPUT,Y-OUTPUT,W-WEIGHTS,CFF-COEFFICIENTS OF FILTER,
C 2*LOT-1=NUMBER OF COEFFICIENTS OF FILTER, 
C NF-NUMBER OF POINTS OF INPUT & OUTPUT SETS
	INTEGER NF,LOT
      REAL X(NF),Y(NF)
      REAL CFF(LOT)
      INTEGER I,L,IPE
C
C FILTRATION AT THE BEGINNING OF THE SET WITH THE ASSUMPTION OF ITS PERIODICITY 
C
      DO I=1,LOT-1
            Y(I)=CFF(1)*X(I)
            DO L=1,LOT-1
			IPE=I-L
			IF(IPE.LT.1) THEN 
				
				IPE=IPE+NF

			END IF
			Y(I)=Y(I)+CFF(L+1)*(X(IPE)+X(I+L))
            END DO
      END DO
C
C FILTRATION IN THE MIDDLE OF THE SET
C
            DO I=LOT,NF-LOT+1
               Y(I)=CFF(1)*X(I)
                  DO L=1,LOT-1
               Y(I)=Y(I)+CFF(L+1)*(X(I-L)+X(I+L))
                  END DO
            END DO
C
C FILTRATION AT THE END OF THE SET WITH THE ASSUMPTION OF ITS PERIODICITY 
C
      DO I=NF-LOT+2,NF

            Y(I)=CFF(1)*X(I)
            DO L=1,LOT-1
			IPE=I+L
              IF(IPE.GT.NF) THEN 
				IPE=IPE-NF
			END IF
              Y(I)=Y(I)+CFF(L+1)*(X(I-L)+X(IPE))
            END DO
      END DO

      END	SUBROUTINE
C======================================================================
      SUBROUTINE NPFILTER_TGR_4(FF,VISCOEF)
      IMPLICIT NONE
      INCLUDE '0COM.INC'
      
C ONE STEP TIME INTEGRATION VISCOSITY OF 4-TH ORDER DIFFUSION ALONG X-DIRECTION
C EXPLICIT PROCEDURE.
C FF   - VARIABLE ON LU-MASK.
C FFL1 - LAPLAS 1ST ORDER OF FF ON LU MASK WITH ZERO FLUX ON BOUNDARY.
C FFL2 - LAPLAS 2ND ORDER OF FF ON LU-MASK WITH ZERO FLUX ON BOUNDARY.
C VISKOEF - UNDIMENSIONAL DIFFUSION  COEFFICIENT
      INTEGER K,N,M
      REAL FF(NX,NY,NZ),VISCOEF

      REAL FFL1(NX,NY,NZ),FFL2(NX,NY,NZ)

      CALL GRIDLAPLAS0_X1D(FF,FFL1)
      IF(MMD.NE.0) CALL CYCLIZE(FFL1,NX,NY,NZ,MMM,MM)
      CALL GRIDLAPLAS0_X1D(FFL1,FFL2)

!$OMP PARALLEL DO PRIVATE (M,N,K)
      DO N=NNN,NN
	 IF(ABS(YT(N)).GT.LAT_CRIT_4D) THEN
        DO M=MMM,MM
         DO K = 1,NZ
          FF(M,N,K)=FF(M,N,K)-VISCOEF*LU(M,N)*FFL2(M,N,K)
         ENDDO
	  END DO
	 END IF
	END DO
!$OMP END PARALLEL DO
      
      IF(MMD.NE.0) CALL CYCLIZE(FF,NX,NY,NZ,MMM,MM)

      CALL GRIDLAPLAS0_X1D(FF,FFL1)
      IF(MMD.NE.0) CALL CYCLIZE(FFL1,NX,NY,NZ,MMM,MM)

!$OMP PARALLEL DO PRIVATE (M,N,K)
      DO N=NNN,NN
	 IF(ABS(YT(N)).GT.LAT_CRIT_4D) THEN
        DO M=MMM,MM
         DO K = 1,NZ
          FF(M,N,K)=FF(M,N,K)+VISCOEF*LU(M,N)*FFL1(M,N,K)
         ENDDO
	  END DO
	 END IF
	END DO
!$OMP END PARALLEL DO

      IF(MMD.NE.0) CALL CYCLIZE(FF,NX,NY,NZ,MMM,MM)

      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDLAPLAS0_X1D(FF,FFL)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C GIVES LAPLAS ON UNIT GRID.
C FF  - INPUT VARIABLE ON LU-MASK.
C FFL - OUTPUT LAPLAS OF FF ON LU-MASK WITH ZERO FLUX ON BOUNDARY.
      REAL    FF(NX,NY,NZ),FFL(NX,NY,NZ)
      INTEGER M, N, K

!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO K = 1,NZ
      DO N=1,NY
       DO M=1,NX
        IF (LU(M,N).GT.0.5.AND.ABS(YT(N)).GT.LAT_CRIT_4D) THEN
        FFL(M,N,K)=1.0/(DX(M,N)*DY(M,N)*HHQ(M,N))
     &   * (
     &      ( FF(M+1,N,K)- FF(M  ,N,K))*DXT(M  ,N)*DYH(M  ,N)
     &     * HHU(M  ,N)* LU(M+1,N)
     &   -
     &      ( FF(M  ,N,K)- FF(M-1,N,K))*DXT(M-1,N)*DYH(M-1,N)
     &     * HHU(M-1,N)* LU(M-1,N)
     &                              )
        ELSE
	  FFL(M,N,K)= 0.0
        ENDIF
       ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
