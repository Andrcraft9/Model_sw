C======================================================================
C INITIALIZATING SHPHERICAL PARAMETERS OF BASIN.
C DIANSKY N.A.,(dinar@inm.ras.ru)
C------------------------------------------------------20.03.02 19:11
      SUBROUTINE BASINPAR()
C----------------------------------------------------------------------
C  INITIALIZING BASIN PARAMETERS FOR BASIN
      IMPLICIT NONE
      INCLUDE '0COM.INC'
      INTEGER M,N,ierr
	
      REAL     PIP180
      REAL(8) DPIP180              !FOR DEGREES TO RADIANS CONVERS
	PARAMETER(PIP180=3.141592653/180.0,
     &         DPIP180=3.1415926535897/180.0D00)
      REAL BCOS89P9           !MINIMAL COS(LATITUDE) VALUE OF 89.9 DEGREES

      REAL S0,T0,A0,B0,S,T,A,B, A_BIG
	REAL COS_LON,COS_LAT,SIN_LON,SIN_LAT, M_X,M_Y
      
      real mask(nx,ny)
      
C-----INITIALIZATION OF GRID TYPE AND TEMPERATURE GRID
      INCLUDE '1BASPAR.INC'

      
      mask=1.0

C PARAMETERS:
      WRITE(*,'(2X,A)')' BASIN PARAMETERS FROM 0COM.INC:'
      WRITE(*,'(2X,2(A,F10.3))')
C     &' INITIAL LATITUDE =',RLAT,'; INITIAL LONGITUDE =',RLON
C      WRITE(*,'(2X,A,F10.3,A,F10.7,A)')
C     &' STEP ON LONGITUDE =',DXST,'[DRG] =',HX,'[rad]',
C     &' STEP ON LATITUDE  =',DYST,'[DGR] =',HY,'[rad]'
      WRITE(*,'(2X,A,G14.7,A,G14.7,A)')
     &' ERTH RADIUS =',RE,'(cm); rn =',RN,'[cm, or undim]'
      WRITE(*,'(2X,A,G14.7,A)')
     &'EARTH ANGULAR VELOCITY(OMEGA) =',OMEGA,'[RAD/SEC]'
      WRITE(*,'(2X,A,G14.7,A)')
     &'       HEAT CAPACITY OF WATER =',CPW,'[J/GR/GRAD] FOR 35%. SAL'
      WRITE(*,'(2X,A,G14.7,A)')
     &'            REFERENCE DENSITY =',RH0,'[GR/CM**3]'
      WRITE(*,'(2X,A,F10.3,A)')
     &'   FREE FALL ACCELERATION(GRV)=',GRV,'[CM/S**2]'
      WRITE(*,'(2X,A,F10.3,A)')
     &'FREE FALL ACCELERATION/RH0 (G)=',G,'[CM/S**2/(GR/CM**3)]'
      WRITE(*,'(2X,A,G14.7,A)')
     &'                      CPW*RH0 =',CPWRH0,'[J/GRAD/CM**3]'

      WRITE(*,'(2X,5(A,I4))')
     &'  NX=',NX, ';  NY=',NY,';  NZ=',NZ,';  MMM=',MMM, '; NNN=',NNN
      WRITE(*,'(2X,5(A,I4))')
     &' MLR=',MLR,';  MM=',MM,';  NN=',NN,'; MDIM=',MDIM,'; MMD=',MMD
C--------------------------------------------------------------------

C VELOCITY GRID INITIALIZATION

C X-COORDINATE (IN DEGREES)

      DO M=1,NX-1
	  XU(M)=(XT(M)+XT(M+1))/2.0
      END DO

C Y-COORDINATE (IN DEGREES)

      DO N=1,NY-1
	  YV(N)=(YT(N)+YT(N+1))/2.0
      END DO


C MINIMAL COS(LATITUDE) VALUE OF 89.9 DEGREES
      BCOS89P9=COS(89.9*PIP180)

C-----METRIC INITIALIZATION IN CASE OF SPHERICAL GRID---------------- 
      IF(CURVE_GRID.EQ.0) THEN 

C-----INITIALIZATION OF T-GRID X-STEPS IN CENTIMETERS
	
       DO N=1,NY
	  DO M=1,NX-1
	    DXT(M,N)=(XT(M+1)-XT(M))*PIP180*(RE/RN) 
     &           *MAX(COS(YT(N)*PIP180),BCOS89P9)
        END DO
       END DO

C-----INITIALIZATION OF T-GRID Y-STEPS IN CENTIMETERS
	 DO N=1,NY-1
	  DO M=1,NX
	    DYT(M,N)=(YT(N+1)-YT(N))*PIP180*(RE/RN)
        END DO
	 END DO

C-----INITIALIZATION OF U-GRID X-STEPS IN CENTIMETERS
       DO N=1,NY
	  DO M=2,NX-1
	    DX(M,N)=(XU(M)-XU(M-1))*PIP180*(RE/RN) 
     &           *MAX(COS(YT(N)*PIP180),BCOS89P9)
        END DO
       END DO

C-----INITIALIZATION OF V-GRID Y-STEPS IN CENTIMETERS
	 DO N=2,NY-1
	  DO M=1,NX
	    DY(M,N)=(YV(N)-YV(N-1))*PIP180*(RE/RN)
        END DO
	 END DO


C-----INITIALIZATION OF H-GRID X-STEPS IN CENTIMETERS
       DO N=1,NY-1
	  DO M=2,NX-1
	    DXH(M,N)=(XU(M)-XU(M-1))*PIP180*(RE/RN) 
     &           *MAX(COS(YV(N)*PIP180),BCOS89P9)
        END DO
       END DO

C-----INITIALIZATION OF H-GRID Y-STEPS IN CENTIMETERS
	 DO N=2,NY-1
	  DO M=1,NX-1
	    DYH(M,N)=(YV(N)-YV(N-1))*PIP180*(RE/RN)
        END DO
	 END DO

C-----INITIALIZATION OF V-GRID X-STEPS IN CENTIMETERS
       DO N=2,NY-1
	  DO M=2,NX-1
	    DXB(M,N)=(XT(M+1)-XT(M))*PIP180*(RE/RN) 
     &           *MAX(COS(YV(N)*PIP180),BCOS89P9)
        END DO
       END DO

C-----INITIALIZATION OF U-GRID Y-STEPS IN CENTIMETERS
	 DO N=2,NY-1
	  DO M=2,NX-1
	    DYB(M,N)=(YT(N+1)-YT(N))*PIP180*(RE/RN)
        END DO
	 END DO


C CORIOLIS TERMS ON H-GRID
       DO N=1,NY-1
        DO M=1,NX-1
          
          RLH(M,N)=2.*OMEGA*(SIN(PIP180*YV(N))*
     &            COS(PIP180*ROTATION_ON_LAT)+
     &            COS(PIP180*YV(N))*COS(PIP180*XU(M))*
     &            SIN(PIP180*ROTATION_ON_LAT))      
        END DO
       END DO

      END IF

C-----END OF METRIC INITIALIZATION IN CASE OF SPHERICAL GRID----------------

C-----METRIC INITIALIZATION IN CASE OF CURVILINEAR GRID WITH DISPLACED NORTH POLE--- 
      IF(CURVE_GRID.EQ.1) THEN

	  S0 = 2.0*TAN((45.0 + Y_POLE/2.0)*DPIP180)
     &        *COS(X_POLE*DPIP180)

        T0 = 2.0*TAN((45.0 + Y_POLE/2.0)*DPIP180)
     &        *SIN(X_POLE*DPIP180)

        A0 = 2.0*TAN((45.0 + Q_POLE/2.0)*DPIP180)
     &        *COS(P_POLE*DPIP180)

        B0 = 2.0*TAN((45.0 + Q_POLE/2.0)*DPIP180)
     &        *SIN(P_POLE*DPIP180)

C-----INITIALIZATION OF T-GRID X-STEPS AND H-GRID Y-STEPS
       DO N=2,NY-1
        DO M=1,NX-1
          CALL SPEC_METR(XU(M),YT(N),X_POLE,Y_POLE,
     &                     P_POLE,Q_POLE,M_X,M_Y)
          
	    DXT(M,N)=(XT(M+1)-XT(M  ))*PIP180*(RE/RN)*M_X
          DYH(M,N)=(YV(N  )-YV(N-1))*PIP180*(RE/RN)*M_Y

	  END DO
	 END DO
C-----INITIALIZATION OF T-GRID Y-STEPS AND H-GRID X-STEPS
       DO N=1,NY-1
        DO M=2,NX-1
          CALL SPEC_METR(XT(M),YV(N),X_POLE,Y_POLE,
     &                     P_POLE,Q_POLE,M_X,M_Y)
          
	    DYT(M,N)=(YT(N+1)-YT(N  ))*PIP180*(RE/RN)*M_Y
	    DXH(M,N)=(XU(M  )-XU(M-1))*PIP180*(RE/RN)*M_X          

	  END DO
	 END DO

C-----INITIALIZATION OF U-GRID X-STEPS AND V-GRID Y-STEPS
       DO N=2,NY-1
        DO M=2,NX-1
          CALL SPEC_METR(XT(M),YT(N),X_POLE,Y_POLE,
     &                     P_POLE,Q_POLE,M_X,M_Y)
	    DX(M,N)=(XU(M)-XU(M-1))*PIP180*(RE/RN)*M_X 
	    DY(M,N)=(YV(N)-YV(N-1))*PIP180*(RE/RN)*M_Y
	  END DO
	 END DO

C-----INITIALIZATION OF V-GRID X-STEPS AND U-GRID Y-STEPS
       DO N=1,NY-1
        DO M=1,NX-1
          CALL SPEC_METR(XU(M),YV(N),X_POLE,Y_POLE,
     &                     P_POLE,Q_POLE,M_X,M_Y)
          
	    DXB(M,N)=(XT(M+1)-XT(M))*PIP180*(RE/RN)*M_X
	    DYB(M,N)=(YT(N+1)-YT(N))*PIP180*(RE/RN)*M_Y

	  END DO
	 END DO

C METRIC AND  CORIOLIS TERMS ON H-GRID
        DO N=1,NY-1
         DO M=1,NX-1    
           
           S = 2.0*TAN((45.0 + YV(N)/2.0)*DPIP180)
     &       *COS(XU(M)*DPIP180)
           T = 2.0*TAN((45.0 + YV(N)/2.0)*DPIP180)
     &       *SIN(XU(M)*DPIP180)
	     A = ((S*S0-T*T0)*(S-A0) + (T*S0+S*T0)*(T-B0)) / 
     &              ((S - A0)**2 + (T - B0)**2)
	     B = ((T*S0+S*T0)*(S-A0) - (S*S0-T*T0)*(T-B0)) / 
     &              ((S - A0)**2 + (T - B0)**2)
           SIN_LAT=(A**2+B**2-4.0)/(A**2+B**2+4.0)
	     RLH(M,N)=2.*OMEGA*SIN_LAT        
         END DO
        END DO 

	END IF
C----- END OF METRIC INITIALIZATION IN CASE OF SPECIAL CURVILINEAR GRID-----


C-----METRIC INITIALIZATION IN CASE OF CURVILINEAR GRID WITH BOTH DISPLACED POLES--- 
      IF(CURVE_GRID.EQ.2) THEN

      A_BIG =  DTAN((45.0 + DBLE(Y_POLE)/2.0)*DPIP180)

C-----INITIALIZATION OF T-GRID X-STEPS AND H-GRID Y-STEPS
       DO N=2,NY-1
        DO M=1,NX-1
          CALL SPEC_METR_2(XU(M),YT(N),X_POLE,Y_POLE,
     &                     P_POLE,M_X,M_Y)
          
	    DXT(M,N)=(XT(M+1)-XT(M  ))*PIP180*(RE/RN)*M_X
          DYH(M,N)=(YV(N  )-YV(N-1))*PIP180*(RE/RN)*M_Y

	  END DO
	 END DO
C-----INITIALIZATION OF T-GRID Y-STEPS AND H-GRID X-STEPS
       DO N=1,NY-1
        DO M=2,NX-1
          CALL SPEC_METR_2(XT(M),YV(N),X_POLE,Y_POLE,
     &                     P_POLE,M_X,M_Y)
          
	    DYT(M,N)=(YT(N+1)-YT(N  ))*PIP180*(RE/RN)*M_Y
	    DXH(M,N)=(XU(M  )-XU(M-1))*PIP180*(RE/RN)*M_X

	  END DO
	 END DO

C-----INITIALIZATION OF U-GRID X-STEPS AND V-GRID Y-STEPS
       DO N=2,NY-1
        DO M=2,NX-1
          CALL SPEC_METR_2(XT(M),YT(N),X_POLE,Y_POLE,
     &                     P_POLE,M_X,M_Y)
          
	    DX(M,N)=(XU(M)-XU(M-1))*PIP180*(RE/RN)*M_X
	    DY(M,N)=(YV(N)-YV(N-1))*PIP180*(RE/RN)*M_Y

	  END DO
	 END DO

C-----INITIALIZATION OF V-GRID X-STEPS AND U-GRID Y-STEPS
       DO N=1,NY-1
        DO M=1,NX-1
          CALL SPEC_METR_2(XU(M),YV(N),X_POLE,Y_POLE,
     &                     P_POLE,M_X,M_Y)
          
	    DXB(M,N)=(XT(M+1)-XT(M))*PIP180*(RE/RN)*M_X
	    DYB(M,N)=(YT(N+1)-YT(N))*PIP180*(RE/RN)*M_Y

	  END DO
	 END DO

C METRIC AND  CORIOLIS TERMS ON H-GRID
        DO N=1,NY-1
         DO M=1,NX-1

            A = DTAN((45.0 + DBLE(YV(N))/2.0)*DPIP180)
     &       *DCOS(DBLE(XU(M)-P_POLE)*DPIP180)

            B = DTAN((45.0 + DBLE(YV(N))/2.0)*DPIP180)
     &       *DSIN(DBLE(XU(M)-P_POLE)*DPIP180)

	      S = ((A*A_BIG+1.0)*(A+A_BIG)+B**2*A_BIG) / 
     &              ((A+A_BIG)**2 + B**2)
	      T = (B*A_BIG*(A+A_BIG)-B*(A*A_BIG+1.0))  / 
     &              ((A+A_BIG)**2 + B**2)
           SIN_LAT=(S**2+T**2-1.0)/(S**2+T**2+1.0)
	     RLH(M,N)=2.*OMEGA*SIN_LAT        
         END DO
        END DO 

	END IF
C----- END OF METRIC INITIALIZATION IN CASE OF SPECIAL CURVILINEAR GRID-----



      IF(MMD.NE.0) THEN
        CALL CYCLIZE(DXT,NX,NY,1,MMM,MM)
        CALL CYCLIZE(DYT,NX,NY,1,MMM,MM)
        CALL CYCLIZE(DX ,NX,NY,1,MMM,MM)
        CALL CYCLIZE(DY ,NX,NY,1,MMM,MM)
        CALL CYCLIZE(DXH,NX,NY,1,MMM,MM)
        CALL CYCLIZE(DYH,NX,NY,1,MMM,MM)
        CALL CYCLIZE(DXB,NX,NY,1,MMM,MM)
        CALL CYCLIZE(DYB,NX,NY,1,MMM,MM)
        CALL CYCLIZE(RLH,NX,NY,1,MMM,MM)
      END IF

      
      
c      CALL WDSTD(' ','dx.dat',1,dx,mask,NX,NY,1,
c     &                      MMM,MM,NNN,NN,1,1,IERR)
c      CALL CREATE_CTLFILE('dx.dat',-1e32,MM-MMM+1,NN-NNN+1,1,1,
c     &         XT(MMM:MM),0.0,YT(NNN:NN),0.0,0.0,3600.0,0.0,'dx','dx')

c      CALL WDSTD(' ','dy.dat',1,dy,mask,NX,NY,1,
c     &                      MMM,MM,NNN,NN,1,1,IERR)
c      CALL CREATE_CTLFILE('dy.dat',-1e32,MM-MMM+1,NN-NNN+1,1,1,
c     &         XT(MMM:MM),0.0,YT(NNN:NN),0.0,0.0,3600.0,0.0,'dy','dy')

c      CALL WDSTD(' ','rlh.dat',1,rlh,mask,NX,NY,1,
c     &                      MMM-1,MM,NNN-1,NN,1,1,IERR)
c      CALL CREATE_CTLFILE('rlh.dat',-1e32,MM-MMM+2,NN-NNN+2,1,1,
c     &        XU(MMM-1:MM),0.0,YV(NNN-1:NN),0.0,0.0,
c     &            3600.0,0.0,'rlh','rlh')


      RETURN
      END

C========METRIC CALCULATION ON CURVILINEAR GRID WITH NORTH POLE DISPLACED=
      SUBROUTINE SPEC_METR(XOUT,YOUT,X_POLE,Y_POLE,
     &                     P_POLE,Q_POLE,M_X,M_Y)
      IMPLICIT NONE
      REAL     PIP180
      REAL(8) DPIP180              !FOR DEGREES TO RADIANS CONVERS
	PARAMETER(PIP180=3.141592653/180.0,
     &         DPIP180=3.1415926535897/180.0D00)
           
      REAL XOUT,YOUT,X_POLE,Y_POLE,P_POLE,Q_POLE,M_X,M_Y
	REAL(8) SIN_LAT,COS_LAT       	!REVERSE ROTATION

      REAL(8) A,B,S,T,A0,B0,S0,T0                   !AUXILARY VARIABLES
	
	REAL(8) DX_DA,DX_DB,DY_DA,DY_DB,
     &        DA_DS,DA_DT,DB_DS,DB_DT,
     &        DS_DP,DS_DQ,DT_DP,DT_DQ,

     &        DA_DP,DA_DQ,DB_DP,DB_DQ,
     &        DX_DP,DX_DQ,DY_DP,DY_DQ,
     &        DET, HP_DIVIDE_R,HQ_DIVIDE_R

      REAL(8) DFM1(2,2),DF(2,2)
	   	   !TRASFORMATION FROM NEW TO OLD GRID
	  S0 = 2.0*TAN((45.0 + Y_POLE/2.0)*DPIP180)
     &        *COS(X_POLE*DPIP180)

        T0 = 2.0*TAN((45.0 + Y_POLE/2.0)*DPIP180)
     &        *SIN(X_POLE*DPIP180)

        A0 = 2.0*TAN((45.0 + Q_POLE/2.0)*DPIP180)
     &        *COS(P_POLE*DPIP180)

        B0 = 2.0*TAN((45.0 + Q_POLE/2.0)*DPIP180)
     &        *SIN(P_POLE*DPIP180)

      S = 2.0*DTAN((45.0 + DBLE(YOUT)/2.0)*DPIP180)
     &       *DCOS(DBLE(XOUT)*DPIP180)

      T = 2.0*DTAN((45.0 + DBLE(YOUT)/2.0)*DPIP180)
     &       *DSIN(DBLE(XOUT)*DPIP180)

	A = ((S*S0-T*T0)*(S-A0) + (T*S0+S*T0)*(T-B0)) / 
     &              ((S - A0)**2 + (T - B0)**2)
	B = ((T*S0+S*T0)*(S-A0) - (S*S0-T*T0)*(T-B0)) / 
     &              ((S - A0)**2 + (T - B0)**2)


C         NESSESARY LATITUDE
      SIN_LAT=(A**2+B**2-4.0)/(A**2+B**2+4.0)
	SIN_LAT=DMIN1(SIN_LAT, DSIN(90.0*DPIP180))
	SIN_LAT=DMAX1(SIN_LAT,-DSIN(90.0*DPIP180))

	COS_LAT=DSQRT(1.0-SIN_LAT**2)


C         DIFFERENTIAL OF TRANSFORMATION

	DX_DA = -B / (A**2 + B**2)
	DX_DB =  A / (A**2 + B**2)

	DY_DA = A / ( DSQRT(A**2 + B**2) * (1.0+ (A**2 + B**2)/4.0))
	DY_DB = B / ( DSQRT(A**2 + B**2) * (1.0+ (A**2 + B**2)/4.0))

	DA_DS = (S0*(S-A0) + S*S0 - T*T0 + T0*(T-B0)) / 
     &           ( (S-A0)**2 + (T-B0)**2 ) -
     &	2.0*( (S*S0 - T*T0)*(S - A0) + (T*S0 + S*T0)*(T - B0) ) * 
     &        (S -A0) / ( (S-A0)**2 + (T-B0)**2 )**2

	DA_DT = (-T0*(S-A0) + S0*(T-B0) + T*S0 + S*T0) /
     &            ( (S-A0)**2 + (T-B0)**2 ) -
     &	2.0*( (S*S0 - T*T0)*(S - A0) + (T*S0 + S*T0)*(T - B0) ) * 
     &        (T - B0) / ( (S-A0)**2 + (T-B0)**2 )**2

	DB_DS = (T0*(S-A0) + T*S0 + S*T0 - S0*(T-B0)) /
     &              ( (S-A0)**2 + (T-B0)**2 ) -
     &        2.0*( (T*S0 + S*T0)*(S-A0) - (S*S0 - T*T0)*(T-B0))*(S-A0)/
     &             ( (S-A0)**2 + (T-B0)**2 )**2
	DB_DT = (S0*(S-A0) + T0*(T-B0) - S*S0 + T*T0) /
     &              ( (S-A0)**2 + (T-B0)**2 ) - 
     &        2.0*( (T*S0 + S*T0)*(S-A0) - (S*S0 - T*T0)*(T-B0))*(T-B0)/
     &             ( (S-A0)**2 + (T-B0)**2 )**2


      DS_DP = -2.0*DTAN((45.0 + DBLE(YOUT)/2.0)*DPIP180)
     &            *DSIN(DBLE(XOUT)*DPIP180)

      DS_DQ = DCOS(DBLE(XOUT)*DPIP180)
     &      /(DCOS((45.0 + DBLE(YOUT)/2.0)*DPIP180))**2

      DT_DP = 2.0*DTAN((45.0 + DBLE(YOUT)/2.0)*DPIP180)
     /           *DCOS(DBLE(XOUT)*DPIP180)

      DT_DQ = DSIN(DBLE(XOUT)*DPIP180) 
     &      /(DCOS((45.0 + DBLE(YOUT)/2.0)*DPIP180))**2

	DA_DP = DA_DS*DS_DP + DA_DT*DT_DP
	DA_DQ = DA_DS*DS_DQ + DA_DT*DT_DQ
	DB_DP = DB_DS*DS_DP + DB_DT*DT_DP
	DB_DQ = DB_DS*DS_DQ + DB_DT*DT_DQ
	
	DX_DP = DX_DA*DA_dP + DX_DB*DB_DP
	DX_DQ = DX_DA*DA_dQ + DX_DB*DB_DQ
	DY_DP = DY_DA*DA_dP + DY_DB*DB_DP
	DY_DQ = DY_DA*DA_dQ + DY_DB*DB_DQ

	HP_DIVIDE_R = DSQRT((DX_DP*COS_LAT)**2 + (DY_DP)**2)
	HQ_DIVIDE_R = DSQRT((DX_DQ*COS_LAT)**2 + (DY_DQ)**2)

      M_X=SNGL(HP_DIVIDE_R)
      M_Y=SNGL(HQ_DIVIDE_R)

      RETURN
	END

C========METRIC CALCULATION ON CURVILINEAR GRID WITH BOTH DISPLACED POLES=
      SUBROUTINE SPEC_METR_2(XOUT,YOUT,X_POLE,Y_POLE,
     &                     L0_POLE,M_X,M_Y)
      IMPLICIT NONE
      REAL     PIP180
      REAL(8) DPIP180              !FOR DEGREES TO RADIANS CONVERS
	PARAMETER(PIP180=3.141592653/180.0,
     &         DPIP180=3.1415926535897/180.0D00)
           
      REAL XOUT,YOUT,X_POLE,Y_POLE,L0_POLE,M_X,M_Y
	REAL(8) SIN_LAT,COS_LAT       	!REVERSE ROTATION

      REAL(8) A,B,S,T,A_BIG                    !AUXILARY VARIABLES
	
	REAL(8) DS_DA,DS_DB,DT_DA,DT_DB,
     &        DX_DS,DX_DT,DY_DS,DY_DT,
     &        DS_DP,DS_DQ,DT_DP,DT_DQ,

     &        DA_DP,DA_DQ,DB_DP,DB_DQ,
     &        DX_DP,DX_DQ,DY_DP,DY_DQ,
     &        DET, HP_DIVIDE_R,HQ_DIVIDE_R


      REAL(8) DFM1(2,2),DF(2,2)
	   	   !TRASFORMATION FROM NEW TO OLD GRID
      A_BIG =  DTAN((45.0 + DBLE(Y_POLE)/2.0)*DPIP180)


      A = DTAN((45.0 + DBLE(YOUT)/2.0)*DPIP180)
     &       *DCOS(DBLE(XOUT-L0_POLE)*DPIP180)

      B = DTAN((45.0 + DBLE(YOUT)/2.0)*DPIP180)
     &       *DSIN(DBLE(XOUT-L0_POLE)*DPIP180)

	S = ((A*A_BIG+1.0)*(A+A_BIG)+B**2*A_BIG) / 
     &              ((A+A_BIG)**2 + B**2)
	T = (B*A_BIG*(A+A_BIG)-B*(A*A_BIG+1.0))  / 
     &              ((A+A_BIG)**2 + B**2)


C         NESSESARY LATITUDE
      SIN_LAT=(S**2+T**2-1.0)/(S**2+T**2+1.0)
	SIN_LAT=DMIN1(SIN_LAT, DSIN(90.0*DPIP180))
	SIN_LAT=DMAX1(SIN_LAT,-DSIN(90.0*DPIP180))

	COS_LAT=DSQRT(1.0-SIN_LAT**2)


C         DIFFERENTIAL OF TRANSFORMATION

      DA_DP = -DTAN((45.0 + DBLE(YOUT)/2.0)*DPIP180)
     &            *DSIN(DBLE(XOUT-L0_POLE)*DPIP180)

      DA_DQ = DCOS(DBLE(XOUT-L0_POLE)*DPIP180)/2.0
     &      /(DCOS((45.0 + DBLE(YOUT)/2.0)*DPIP180))**2

      DB_DP =  DTAN((45.0 + DBLE(YOUT)/2.0)*DPIP180)
     /           *DCOS(DBLE(XOUT-L0_POLE)*DPIP180)

      DB_DQ = DSIN(DBLE(XOUT-L0_POLE)*DPIP180)/2.0 
     &      /(DCOS((45.0 + DBLE(YOUT)/2.0)*DPIP180))**2

      DS_DA=(1.0+A*A_BIG+A_BIG*(A+A_BIG))
     &     /((A+A_BIG)**2+B**2)
     &  -
     &      2.0*(A+A_BIG)*((A+A_BIG)*(1.0+A*A_BIG)+A_BIG*B**2)
     &     /((A+A_BIG)**2+B**2)**2

      
      DS_DB=2.0*A_BIG*B/((A+A_BIG)**2+B**2)
     &  -
     &      2.0*B*((A+A_BIG)*(1.0+A*A_BIG)+A_BIG*B**2)
     &     /((A+A_BIG)**2+B**2)**2
      
      
	DT_DA=-(2.0*(A+A_BIG)*(A_BIG*(A+A_BIG)*B
     &           -(1.0+A*A_BIG)*B))
     &     /((A+A_BIG)**2+B**2)**2

      
      DT_DB=(A_BIG*(A+A_BIG)-1.0-A*A_BIG)
     &     /((A+A_BIG)**2+B**2)
     &  -    2.0*B*(A_BIG*(A+A_BIG)*B-(1.0+A*A_BIG)*B)
     &     /((A+A_BIG)**2+B**2)**2

      
      DX_DS=-T/(S**2+T**2)
	DX_DT= S/(S**2+T**2)
	DY_DS=2.0*S/(DSQRT(S**2+T**2)*(1.0+S**2+T**2))
	DY_DT=2.0*T/(DSQRT(S**2+T**2)*(1.0+S**2+T**2))

      DS_DP=DS_DA*DA_DP+DS_DB*DB_DP
      DS_DQ=DS_DA*DA_DQ+DS_DB*DB_DQ
      DT_DP=DT_DA*DA_DP+DT_DB*DB_DP
      DT_DQ=DT_DA*DA_DQ+DT_DB*DB_DQ

	DX_DP=DX_DS*DS_DP+DX_DT*DT_DP
	DY_DP=DY_DS*DS_DP+DY_DT*DT_DP
	DX_DQ=DX_DS*DS_DQ+DX_DT*DT_DQ
	DY_DQ=DY_DS*DS_DQ+DY_DT*DT_DQ


	HP_DIVIDE_R = DSQRT((DX_DP*COS_LAT)**2 + (DY_DP)**2)
	HQ_DIVIDE_R = DSQRT((DX_DQ*COS_LAT)**2 + (DY_DQ)**2)

      M_X=SNGL(HP_DIVIDE_R)
      M_Y=SNGL(HQ_DIVIDE_R)

      RETURN
	END

C========TRANSFORMATION OF MODEL COORDINATES TO GEOGRAPHICAL ONES ON CURVILINEAR GRID WITH BOTH DISPLACED POLES=
      SUBROUTINE XMYM_2_XGYG_2(XM,YM,X_POLE,Y_POLE,
     &                     L0_POLE,XG,YG)
      IMPLICIT NONE
      REAL     PIP180
      REAL(8) DPIP180              !FOR DEGREES TO RADIANS CONVERS
	PARAMETER(PIP180=3.141592653/180.0,
     &         DPIP180=3.1415926535897/180.0D00)
           
      REAL XM,YM,X_POLE,Y_POLE,L0_POLE,XG,YG
	REAL(8) SIN_LAT,COS_LAT,SIN_LON,COS_LON,RET_LON,RET_LAT   !REVERSE ROTATION

      REAL(8) A,B,S,T,A_BIG                    !AUXILARY VARIABLES
	
	   	   !TRASFORMATION FROM NEW TO OLD GRID
      A_BIG =  DTAN((45.0 + DBLE(Y_POLE)/2.0)*DPIP180)


      A = DTAN((45.0 + DBLE(YM)/2.0)*DPIP180)
     &       *DCOS(DBLE(XM-L0_POLE)*DPIP180)

      B = DTAN((45.0 + DBLE(YM)/2.0)*DPIP180)
     &       *DSIN(DBLE(XM-L0_POLE)*DPIP180)

	S = ((A*A_BIG+1.0)*(A+A_BIG)+B**2*A_BIG) / 
     &              ((A+A_BIG)**2 + B**2)
	T = (B*A_BIG*(A+A_BIG)-B*(A*A_BIG+1.0))  / 
     &              ((A+A_BIG)**2 + B**2)


C         NESSESARY LONGITUDE

      COS_LON=S/DSQRT(S**2+T**2)
	SIN_LON=T/DSQRT(S**2+T**2)

	    SIN_LON=DMIN1(SIN_LON, 1.0D0)
	    SIN_LON=DMAX1(SIN_LON,-1.0D0)
	    COS_LON=DMIN1(COS_LON, 1.0D0)
	    COS_LON=DMAX1(COS_LON,-1.0D0)

          RET_LON=DSIGN(DACOS(COS_LON)/DPIP180,SIN_LON)+DBLE(X_POLE)

C         NESSESARY LATITUDE
      SIN_LAT=(S**2+T**2-1.0)/(S**2+T**2+1.0)
	SIN_LAT=DMIN1(SIN_LAT, 1D0)
	SIN_LAT=DMAX1(SIN_LAT,-1D0)

	RET_LAT=DASIN(SIN_LAT)/DPIP180


      XG=SNGL(RET_LON)
      YG=SNGL(RET_LAT)

      RETURN
	END
C========TRANSFORMATION OF MODEL COORDINATES TO GEOGRAPHICAL ONES ON SPHERICAL GRID WITH DISPLACED POLES=
      SUBROUTINE XMYM_2_XGYG_SPH(XM,YM,LAMBDA_ROT,PHI_ROT,XG,YG,
     &           ROTVEC_COEFF)
      IMPLICIT NONE
      REAL     PIP180
      REAL(8) DPIP180              !FOR DEGREES TO RADIANS CONVERS
	PARAMETER(PIP180=3.141592653/180.0,
     &         DPIP180=3.1415926535897/180.0D00)
           
      REAL XM,YM,	LAMBDA_ROT,PHI_ROT, XG,YG,ROTVEC_COEFF(4)
	REAL(8) SIN_LAT,COS_LAT,SIN_LON,COS_LON,RET_LON,RET_LAT,   !REVERSE ROTATION
     &        FREE_MEMBER_COSLON,FREE_MEMBER_SINLON

	
	   !REVERSE ROTATION
	    SIN_LAT = DSIN(DPIP180*DBLE(YM)) * 
     &              DCOS(DPIP180*DBLE(PHI_ROT)) +
     &              DCOS(DPIP180*DBLE(XM)) *
     &              DCOS(DPIP180*DBLE(YM)) *
     &              DSIN(DPIP180*DBLE(PHI_ROT))

	    COS_LAT = DSQRT(1D0-SIN_LAT**2)
        
	    RET_LAT = DASIN(SIN_LAT)/DPIP180

C         NESSESARY LATITUDE
	    FREE_MEMBER_COSLON =( DCOS(DPIP180*DBLE(XM)) *
     &                          DCOS(DPIP180*DBLE(YM)) *
     &                          DCOS(DPIP180*DBLE(PHI_ROT)) -
     &                          DSIN(DPIP180*DBLE(YM)) *
     &                          DSIN(DPIP180*DBLE(PHI_ROT))   )/COS_LAT

	    FREE_MEMBER_SINLON =( DSIN(DPIP180*DBLE(XM)) *
     &                          DCOS(DPIP180*DBLE(YM)) )/COS_LAT
	    
		COS_LON=FREE_MEMBER_COSLON*DCOS(DPIP180*DBLE(LAMBDA_ROT))
     &           -FREE_MEMBER_SINLON*DSIN(DPIP180*DBLE(LAMBDA_ROT))

		SIN_LON=FREE_MEMBER_SINLON*DCOS(DPIP180*DBLE(LAMBDA_ROT))
     &           +FREE_MEMBER_COSLON*DSIN(DPIP180*DBLE(LAMBDA_ROT))
	   
	    SIN_LON=DMIN1(SIN_LON, 1.0D0)
	    SIN_LON=DMAX1(SIN_LON,-1.0D0)
	    COS_LON=DMIN1(COS_LON, 1.0D0)
	    COS_LON=DMAX1(COS_LON,-1.0D0)

C         NESSESARY LONGITUDE
          RET_LON=DSIGN(DACOS(COS_LON)/DPIP180,SIN_LON)


      XG=SNGL(RET_LON)
      YG=SNGL(RET_LAT)

C--------DEFINITION OF ANGLES BETWEEN PARALLELS-----------------------
	  ROTVEC_COEFF(1)
     &               = SNGL (( COS_LAT*DCOS(DPIP180*DBLE(PHI_ROT))+
     &                         SIN_LAT*DSIN(DPIP180*DBLE(PHI_ROT)) *
     &                       ( COS_LON*DCOS(DPIP180*DBLE(LAMBDA_ROT))+
     &                         SIN_LON*DSIN(DPIP180*DBLE(LAMBDA_ROT))))
     &                                /DCOS(DPIP180*DBLE(YM)))

        ROTVEC_COEFF(2)
     &               = SNGL ((-DSIN(DPIP180*DBLE(PHI_ROT)) *
     &                       ( SIN_LON*DCOS(DPIP180*DBLE(LAMBDA_ROT))-
     &                         COS_LON*DSIN(DPIP180*DBLE(LAMBDA_ROT))))
     &                                /DCOS(DPIP180*DBLE(YM)))

	  ROTVEC_COEFF(3)=-ROTVEC_COEFF(2)
        ROTVEC_COEFF(4)= ROTVEC_COEFF(1)

      RETURN
	END
C=========CALCULATING NUMBER OF COORDINATE IN GRID ARRAY=========================
      SUBROUTINE GRID_NUM_FROM_COORD(GRD,NGRD,COORD,NCOORD,MMM,MM)
      IMPLICIT NONE
      INTEGER NGRD,      !GRID DIMENSION(INPUT)
     &        NCOORD     !NUM OF COORDINATE TO BE FOUND (OUTPUT)

      REAL GRD(NGRD),    !GRID ARRAY(INPUT)
     &     COORD         !COORDINATE (INPUT)

      INTEGER MMM,MM
	INTEGER M,N,I_NES

        I_NES = MMM
	 
	  N     = MM

1000    M=(I_NES+N)/2	 
	     
           IF (COORD.LT.GRD(M)) N     = M
	     IF (COORD.GE.GRD(M)) I_NES = M
	     IF (N.GT.I_NES+1) GOTO 1000

         IF(ABS(GRD(I_NES)-COORD).LT.ABS(GRD(I_NES+1)-COORD)) THEN
         	 NCOORD=I_NES
         ELSE
         	 NCOORD=I_NES+1
	   END IF
 
            
	RETURN
	END
C=============VECTOR FIELD TRANSFORMATION BY COORDINATE ROTATION
      SUBROUTINE VECTOR_TRANSFORM(FIELD_X_IN, FIELD_Y_IN,
     &                            FIELD_X_OUT,FIELD_Y_OUT,ROTVEC_COEFF)
      IMPLICIT NONE
	REAL FIELD_X_IN, FIELD_Y_IN,
     &     FIELD_X_OUT,FIELD_Y_OUT,ROTVEC_COEFF(4)

      FIELD_X_OUT=FIELD_X_IN*ROTVEC_COEFF(1)+FIELD_Y_IN*ROTVEC_COEFF(2)
      FIELD_Y_OUT=FIELD_X_IN*ROTVEC_COEFF(3)+FIELD_Y_IN*ROTVEC_COEFF(4)

	RETURN
	END
