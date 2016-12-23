C======================================================================
C(C) DIANSKY N.A.(dinar@inm.ras.ru),GUSEV A.V(glavpip@pochta.ru),18.01.2006
C PROGRAM MODULE OF THE OGCM INTEGRATION OVER ONE STEP
C----------------------------------------------------------------------
      SUBROUTINE ADAPTATION_MODULE(TAU,NSLOR,IGT,IGS,IGWS)
      use mpi_parallel_tools
      use mod_shallowater
      IMPLICIT NONE
	  INCLUDE '0COM.INC'
      INCLUDE '0CEAN.INC'
      INCLUDE '0FUNCDEF.INC'

      REAL    TAU,WFLUX_AVER       !TIME STEP IN SECONDS(INPUT)
      INTEGER NSLOR ! AND NBR. OF ITERATIONS(OUTPUT)
      INTEGER IGICE,IGT,IGS,IGWS,NICESTEP
      INTEGER M, N, K
      REAL*8 FREESURFACE

C**********************************************************************C
      FREESURFACE =1D0

C*********************ADAPTATION MODULE********************************

C CALCULATION ZONAL VELOCITY FORCED BY PRESSURE GRADIENT
c        CALL U_BY_PRESSURE_GRADIENT(UU,DEN,TAU)
C CALCULATION MERIDIONAL VELOCITY FORCED BY PRESSURE GRADIENT
c        CALL V_BY_PRESSURE_GRADIENT(VV,DEN,TAU)

C COMPUTE BAROTROPIC COMPONENTS THEN REMOVE THEIR FROM 3-D VELOCITIES
C          PERODICITY INCLUDED AUTOMATICALLY

c        CALL UBARK8 (UU,UBRTR,-1,DZ,LCU,NX,NY,NZ) !REMOVE BAROTROP. COMP
c        CALL UBARK8 (VV,VBRTR,-1,DZ,LCV,NX,NY,NZ) !REMOVE BAROTROP. COMP

c        IF (IABS(KSW_UVL).GE.2.OR.KSW_FI.GT.1) THEN  ! U,V ARE CHANGED
C CALCULATE CORIOLIS ROTATION AND INTERNAL GRAVITY VAWES (PERIODICITY INCLUDED)
c          CALL INTERNAL_INERTIA_OSCILLATION(UU,VV,TAU)
c        END IF

        SLPR = 0.0
        do n = max(ny_start, NNN), min(ny_end, NN)
          do m = max(nx_start, MMM), min(nx_end, MM)
              if (LCU(M,N).GT.0.5) then
                DINX(M,N)=-0.0001
              else
    	        DINX(M,N)=0.0
              endif

              if (LCV(M,N).GT.0.5) then
                DINY(M,N)=-0.0005
              else
    	        DINY(M,N)=0.0
              endif
          enddo
        enddo
C CALCULATING GRADIENTS OF ATMOSPHERIC PRESSURE
c        IF(IABS(IGWS).EQ.2.OR.IABS(IGWS).EQ.4) THEN
c            CALL ATM_PRESSURE_GRADIENTS(SLPR,DINX,DINY)
c        ELSE
c		    DINX=0.0
c		    DINY=0.0
c        END IF


C################## SEA LEVEL MODULE ##################################
      IF(KSW_FI.GT.1) THEN
        SLH0=SLH
	  END IF

C      IF(KSW_FI.GT.1.AND.IGS.GT.1) THEN
C        WFLUX_AVER=WAVER(QSAL+RUNOFF)
C!$OMP PARALLEL DO PRIVATE(M,N)
C        DO N=NNN,NN
C	    DO M=MMM,MM
C	     SLH(M,N)=SLH(M,N)
C     &             -DBLE((QSAL(M,N)+RUNOFF(M,N)-WFLUX_AVER)*TAU)
C          END DO
C	  END DO
C!$OMP END PARALLEL DO
C        IF(MMD.NE.0) CALL CYCLIZE8(SLH,NX,NY,1,MMM,MM)
C	END IF

c        print *, slh(1:NX:10, 1:NY:5)

      IF (KSW_FI.GE.5) THEN
C##############  SEA LEVEL ZALESNY - RUSAKOV #####################
        CALL SEALEVEL_ZALRUS(TAU, FREESURFACE, NSLOR)
      END IF

C############ DEFINE 3-D VELOCITY ON NEW TIME STEP:#####################
C DEFINE VERTICAL VELOCITY ON T-GRID

c       CALL WWINTC(UU,VV,WW)
C--------NORTH POLE FILTRATION------------------------------------------
c	   CALL NPS_FILTER_T_GRD(WW,NZ ) !3D ARRAY FILTERED & ITS DIM ON Z
C----------------------------------------------------------------------
       IF(MMD.NE.0) CALL CYCLIZE(WW,NX,NY,NZ+1,MMM,MM)

C ADD BAROTROPIC COMPONENTS TO BAROCLINIC VELOCITY (PERIODICITY INCLUDED)
!$OMP PARALLEL DO PRIVATE(K)
c       DO K=1,NZ
c	     UU(:,:,K) = UU(:,:,K) + SNGL(UBRTR(:,:))
c	     VV(:,:,K) = VV(:,:,K) + SNGL(VBRTR(:,:))
c	   END DO
!$OMP END PARALLEL DO
C############ END OF DEFINITION 3-D VELOCITY ON NEW TIME STEP:############
c restriction on SLH
c	SLH =MAX(SLH ,-500.0D0)
c	SLH =MIN(SLH ,+500.0D0,HHQ-10.0)

C*********************END OF ADAPTATION MODULE*************************
      RETURN
      END

C##############  SEA LEVEL ZALESNY - RUSAKOV #####################

      SUBROUTINE SEALEVEL_ZALRUS(TAU, FREESURFACE,NSLOR)
      use mpi_parallel_tools
      use mod_shallowater
	  use sea_level_no_split
	  IMPLICIT NONE
      INCLUDE '0CEAN.INC'
      INCLUDE '0FUNCDEF.INC'
	  real :: tau
      real*8 :: FREESURFACE
	  integer :: M, N, ierr, NSLOR, iii
      integer, save :: init = 0
      real*8 :: w_time

      RBOTTOM = RDR

      IF (MMD.NE.0) THEN
          CALL CYCLIZE8(RBOTTOM,NX,NY,1,MMM,MM)
      END IF

      if (init.eq.0) then
          call start_timer(w_time)
          call estimate_matrix_size()
          call FORM_MATRICE_SLH(tau, freesurface)
          init = init + 1
          call end_timer(w_time)
          if (rank.eq.0) print*,'MatForm:',w_time
      endif

      call start_timer(w_time)
      call FORM_RHS(tau, freesurface)
      call end_timer(w_time)
      if (rank.eq.0) print*,'RHSForm:',w_time

      call start_timer(w_time)
      call solve_system()
      call end_timer(w_time)
      if (rank.eq.0) print*,'SolveSys:',w_time

c      call PetscFinalize(ierr)
c      stop

c	    CALL SOLVE_MATRICE_ITER(UBRTR,VBRTR,SLH,NSLOR) !iterative solver ! should be comment off is direct solver is used

      IF(MMD.NE.0) THEN
          CALL CYCLIZE8(UBRTR, NX,NY, 1,MMM,MM)
          CALL CYCLIZE8(VBRTR, NX,NY, 1,MMM,MM)
          CALL CYCLIZE8(SLH,   NX,NY, 1,MMM,MM)
      END IF


      RETURN
      END
C############## END SEA LEVEL ZALESNY - RUSAKOV #####################
