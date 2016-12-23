C<<<<<<<< DEFINE DIFFUSIVE AND VISKOSITY COEFFICIENTS OF SECOND ORDER>>>>>
      DO N=NNN,NN
       DO M=MMM,MM
C DEFINE ANALITIC FUNCTIONS OF LATITUDE
	POINTLAT=RLAT+DYST*FLOAT(N-NNN)
         
         RLU=(RLH(M,N)+RLH(M  ,N-1))/2.0
         
         RLV=(RLH(M,N)+RLH(M-1,N  ))/2.0
         
         RLT=(RLH(M,N  )+RLH(M-1,N  )
     &       +RLH(M,N-1)+RLH(M-1,N-1))/4.0

         DO K=1,NZ

C VARIANT &1 WITH CORRECTION ON LATITUDE FOR EQUATOR
            
C        FOR TEMPERATURE AND SALINITY

c          AMYT(M,N,K)=AMUYT*DYT(M,N)/(RE/RN)/PI*180.0
c     &    *(1.0+0.0*(1.0-(RLV/2.0/OMEGA)**2)**50)  !Y-DIRECTION ON V-GRID

c         AMXT(M,N,K)=AMUYT*DXT(M,N)/(RE/RN)/PI*180.0 
c     &    *(1.0+0.0*(1.0-(RLU/2.0/OMEGA)**2)**50)  !X-DIRECTION ON U-GRID
c     &               *AMUXT 
c         AMXT(M,N,K)=MAX(AMXT(M,N,K),AMXT(M,NUM_LAT_BEG_FILTER,K))         
         
         
         AMYT(M,N,K)=AMUYT
         AMXT(M,N,K)=AMUYT*AMUXT
c        AMXT(M,N,K)=AMUYT*AMUXT*COS(YT(N)*PI/180.0)  

C         FOR VELOCITY AT T-GRID

c          AMYU(M,N,K)=AMUYU *DY(M,N)/(RE/RN)/PI*180.0
c     &    *(0.7+0.3*(1.0-(RLT/2.0/OMEGA)**2)**50)   !Y-DIRECTION
c          AMXU(M,N,K)=AMUYU *DX(M,N)/(RE/RN)/PI*180.0
c     &    *(0.7+0.3*(1.0-(RLT/2.0/OMEGA)**2)**50)  !X-DIRECTION
c     &              *AMUXU
c	    AMXU(M,N,K)=MAX(AMXU(M,N,K),AMXU(M,NUM_LAT_BEG_FILTER,K))

         AMYU(M,N,K)=AMUYU
         AMXU(M,N,K)=AMUYU*AMUXU
c        AMXU(M,N,K)=AMUYU*AMUXU*COS(YT(N)*PI/180.0)
          
          END DO
        ENDDO
      ENDDO