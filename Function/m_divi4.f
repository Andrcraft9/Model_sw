C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDVISCOSITYT4(FF,VISCOEF)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C ONE STEP TIME INTEGRATION VISCOSITY OF 4-TH ORDER DIFFUSION
C EXPLICIT PROCEDURE.
C FF   - VARIABLE ON LU-MASK.
C FFL1 - LAPLAS 1ST ORDER OF FF ON LU MASK WITH ZERO FLUX ON BOUNDARY.
C FFL2 - LAPLAS 2ND ORDER OF FF ON LU-MASK WITH ZERO FLUX ON BOUNDARY.
C VISKOEF - UNDIMENSIONAL DIFFUSION  COEFFICIENT
      INTEGER K
      REAL FF(NX,NY,NZ),VISCOEF

      REAL FFL1(NX,NY,NZ),FFL2(NX,NY,NZ)

      CALL GRIDLAPLAS0(FF,FFL1)
      IF(MMD.NE.0) CALL CYCLIZE(FFL1,NX,NY,NZ,MMM,MM)
      CALL GRIDLAPLAS0(FFL1,FFL2)

!$OMP PARALLEL DO
      DO K = 1,NZ
      FF(:,:,K)=FF(:,:,K)-VISCOEF*LU(:,:)*FFL2(:,:,K)
      ENDDO
!$OMP END PARALLEL DO
      
      IF(MMD.NE.0) CALL CYCLIZE(FF,NX,NY,NZ,MMM,MM)

      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDLAPLAS0(FF,FFL)
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
        IF (LU(M,N).GT.0.5) THEN
        FFL(M,N,K)=1.0/(DX(M,N)*DY(M,N)*HHQ(M,N))
     &   * (
     &      ( FF(M+1,N,K)- FF(M  ,N,K))*DXT(M  ,N)*DYH(M  ,N)
     &     * HHU(M  ,N)* LU(M+1,N)
     &   -
     &      ( FF(M  ,N,K)- FF(M-1,N,K))*DXT(M-1,N)*DYH(M-1,N)
     &     * HHU(M-1,N)* LU(M-1,N)
     &   +
     &      ( FF(M,N+1,K)- FF(M,N  ,K))*DYT(M,N  )*DXH(M,N  )
     &     * HHV(M,N  )* LU(M,N+1)
     &   -
     &      ( FF(M,N  ,K)- FF(M,N-1,K))*DYT(M,N-1)*DXH(M,N-1)
     &     * HHV(M,N-1)* LU(M,N-1)
     &                                              )
        ELSE
	  FFL(M,N,K)= 0.0
        ENDIF
       ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
C======================================19.07.00 14:50(dinar@inm.ras.ru)
      SUBROUTINE NPSFILTER(FF,NZZ,BDGR,IGRID)
      IMPLICIT NONE
C----------------------------------------------------------------------
C TO FILTER FF OVER NP IN X-DIRECTION BY NANSEN FILTER TO REMOVE 2-STEP
C OSCILATIONS
C FF - ARRAY TO BE FILTERED
C NZZ - Z-DIMENSION
C BDGR - BOUNDARY LATITUDE FOR FILTRATION
C IGRID = 1 THEN T-GRID
C IGRID = 2 THEN U-GRID
C IGRID = 3 THEN V-GRID

      INCLUDE '0COM.INC'
      INTEGER NZZ, NF
      INTEGER IGRID
      REAL    FF(NX,NY,NZZ),BDGR,FLTR(NX)
      INTEGER M, N, K
      REAL    FWEITS

C     DGR(N)= (RLAT+DYST/2.+(N-NNN)*DYST)

      IF(IGRID.EQ.3) THEN
C BOUNDARY INDEX FOR V-GRID
      NF=NINT((BDGR-RLAT-DYST/2.)/DYST)+NNN
      ELSE
C BOUNDARY INDEX FOR T-GRID AND U-GRID
      NF=NINT((BDGR-RLAT        )/DYST)+NNN
      END IF

      DO K = 1,NZZ
      DO N=NF,NN
C INITIAL VALUES
       DO M = MMM-1,MM+1
       FLTR(M)=FF(M,N,K)
       ENDDO

C NANSEN FILTRATION
      IF(IGRID.EQ.1) THEN

C FILTRATION ON T-GRID
       DO M = MMM, MM
       IF(LU(M,N).GT.0.5) THEN
       FWEITS=0.25*LU(M-1,N)+0.5+0.25*LU(M+1,N)
       FF(M,N,K)=( 0.25*FLTR(M-1)*LU(M-1,N)
     &            +0.50*FLTR(M  )
     &            +0.25*FLTR(M+1)*LU(M+1,N) )/FWEITS
       END IF
       ENDDO

      ELSE IF(IGRID.EQ.2) THEN

C FILTRATION ON U-GRID
       DO M = MMM, MM
       IF(LCU(M,N).GT.0.5) THEN
       FWEITS=0.25*LCU(M-1,N)+0.5+0.25*LCU(M+1,N)
       FF(M,N,K)=( 0.25*FLTR(M-1)*LCU(M-1,N)
     &            +0.50*FLTR(M  )
     &            +0.25*FLTR(M+1)*LCU(M+1,N) )/FWEITS
       END IF
       ENDDO

      ELSE IF(IGRID.EQ.3) THEN
C FILTRATION ON V-GRID
       DO M = MMM,MM
       IF(LCV(M,N).GT.0.5) THEN
       FWEITS=0.25*LCV(M-1,N)+0.5+0.25*LCV(M+1,N)
       FF(M,N,K)=( 0.25*FLTR(M-1)*LCV(M-1,N)
     &            +0.50*FLTR(M  )
     &            +0.25*FLTR(M+1)*LCV(M+1,N) )/FWEITS
       END IF
       ENDDO

      END IF

      ENDDO
      ENDDO

      IF(MMD.NE.0) THEN
         CALL CYCLIZE(FF,NX,NY,NZ,MMM,MM)
      END IF

      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDVISCOSITYU4(FF,VISCOEF)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C ONE STEP TIME INTEGRATION VISCOSITY OF 4-TH ORDER DIFFUSION
C EXPLICIT PROCEDURE.
C FF   - VARIABLE ON LCU-MASK.
C FFL1 - LAPLAS 1ST ORDER OF FF ON LCU MASK WITH ZERO ON U-BOUNDARY.
C FFL2 - LAPLAS 2ND ORDER OF FF ON LCU-MASK WITH ZERO FLUX ON BOUNDARY.
C VISKOEF - UNDIMENSIONAL DIFFUSION  COEFFICIENT
      INTEGER K
      REAL FF(NX,NY,NZ),VISCOEF


      REAL FFL1(NX,NY,NZ),FFL2(NX,NY,NZ)

	CALL GRIDLAPLASU (FF,FFL1)
      IF(MMD.NE.0) CALL CYCLIZE(FFL1,NX,NY,NZ,MMM,MM)
	CALL GRIDLAPLASU0(FFL1,FFL2)

!$OMP PARALLEL DO
      DO K = 1,NZ
      FF(:,:,K)=FF(:,:,K)-VISCOEF*LCU(:,:)*FFL2(:,:,K)
      ENDDO
!$OMP END PARALLEL DO
      IF(MMD.NE.0) CALL CYCLIZE(FF,NX,NY,NZ,MMM,MM)

      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDVISCOSITYV4(FF,VISCOEF)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C ONE STEP TIME INTEGRATION VISCOSITY OF 4-TH ORDER DIFFUSION
C EXPLICIT PROCEDURE.
C FF   - VARIABLE ON LCV-MASK.
C FFL1 - LAPLAS 1ST ORDER OF FF ON LCV MASK WITH ZERO ON V-BOUNDARY.
C FFL2 - LAPLAS 2ND ORDER OF FF ON LCV-MASK WITH ZERO FLUX ON BOUNDARY.
C VISKOEF - UNDIMENSIONAL DIFFUSION  COEFFICIENT
      INTEGER K
      REAL FF(NX,NY,NZ),VISCOEF


      REAL FFL1(NX,NY,NZ),FFL2(NX,NY,NZ)

      CALL GRIDLAPLASV(FF,FFL1)
      IF(MMD.NE.0) CALL CYCLIZE(FFL1,NX,NY,NZ,MMM,MM)
      CALL GRIDLAPLASV0(FFL1,FFL2)

!$OMP PARALLEL DO
      DO K = 1,NZ
      FF(:,:,K)=FF(:,:,K)-VISCOEF*LCV(:,:)*FFL2(:,:,K)
      ENDDO
!$OMP END PARALLEL DO
      IF(MMD.NE.0) CALL CYCLIZE(FF,NX,NY,NZ,MMM,MM)
      
      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDLAPLASU(FF,FFL)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C GIVES LAPLAS ON UNIT GRID.
C FF  - INPUT VARIABLE ON LCU-MASK.
C FFL - OUTPUT LAPLAS OF FF ON LCU-MASK WITH ZERO ON U-BOUNDARY.
      REAL    FF(NX,NY,NZ),FFL(NX,NY,NZ)
      INTEGER M, N, K

!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO K = 1,NZ
      DO N=1,NY
       DO M=1,NX
        IF (LCU(M,N).GT.0.5) THEN

         FFL(M,N,K)=1.0/(DXT(M,N)*DYH(M,N)*HHU(M,N))
     &   * (
     &      ( FF(M+1,N,K)- FF(M  ,N,K))*DX(M+1,N)*DY(M+1,N)*HHQ(M+1,N  )
     &   -
     &      ( FF(M  ,N,K)- FF(M-1,N,K))*DX(M  ,N)*DY(M  ,N)*HHQ(M  ,N  )
     &   +
     &      ( FF(M,N+1,K)- FF(M,N  ,K))*DXB(M,N  )*DYB(M,N  )
     &     *  HH(M,N  )*LCU(M,N+1)
     &   -
     &      ( FF(M,N  ,K)- FF(M,N-1,K))*DXB(M,N-1)*DYB(M,N-1)
     &     *  HH(M,N-1)*LCU(M,N-1)
     &                                              )
        ELSE
	  FFL(M,N,K)= 0.0
        ENDIF
       ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDLAPLASU0(FF,FFL)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C GIVES LAPLAS ON UNIT GRID.
C FF  - INPUT VARIABLE ON LCU-MASK.
C FFL - OUTPUT LAPLAS OF FF ON LCU-MASK WITH ZERO FLUX ON BOUNDARY.
      REAL    FF(NX,NY,NZ),FFL(NX,NY,NZ)
      INTEGER M, N, K

!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO K = 1,NZ
      DO N=1,NY
       DO M=1,NX
        IF (LCU(M,N).GT.0.5) THEN

         FFL(M,N,K)=1.0/(DXT(M,N)*DYH(M,N)*HHU(M,N))
     &   * (
     &      ( FF(M+1,N,K)- FF(M  ,N,K))*DX(M+1,N)*DY(M+1,N)
     &     *HHQ(M+1,N  )*LCU(M+1,N)
     &   -
     &      ( FF(M  ,N,K)- FF(M-1,N,K))*DX(M  ,N)*DY(M  ,N)
     &     *HHQ(M  ,N  )*LCU(M-1,N)
     &   +
     &      ( FF(M,N+1,K)- FF(M,N  ,K))*DXB(M,N  )*DYB(M,N  )
     &     *  HH(M,N  )*LCU(M,N+1)
     &   -
     &      ( FF(M,N  ,K)- FF(M,N-1,K))*DXB(M,N-1)*DYB(M,N-1)
     &     *  HH(M,N-1)*LCU(M,N-1)
     &                                              )
        ELSE
	  FFL(M,N,K)= 0.0
        ENDIF
       ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDLAPLASV(FF,FFL)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C GIVES LAPLAS ON UNIT GRID.
C FF  - INPUT VARIABLE ON LU-MASK.
C FFL - OUTPUT LAPLAS OF FF ON LCV-MASK WITH ZERO ON V-BOUNDARY.
      REAL    FF(NX,NY,NZ),FFL(NX,NY,NZ)
      INTEGER M, N, K

!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO K = 1,NZ
      DO N=1,NY
       DO M=1,NX
        IF (LCV(M,N).GT.0.5) THEN

         FFL(M,N,K)=1.0/(DXH(M,N)*DYT(M,N)*HHV(M,N))
     &   * (
     &      ( FF(M,N+1,K)- FF(M,N  ,K))*DX(M,N+1)*DY(M,N+1)*HHQ(M,N+1)
     &   -
     &      ( FF(M,N  ,K)- FF(M,N-1,K))*DX(M  ,N)*DY(M  ,N)*HHQ(M,N  )
     &   +
     &      ( FF(M+1,N,K)- FF(M  ,N,K))*DXB(M  ,N)*DYB(M  ,N)
     &     * HH(M  ,N)*LCV(M+1,N)
     &   -
     &      ( FF(M  ,N,K)- FF(M-1,N,K))*DXB(M-1,N)*DYB(M-1,N)
     &     * HH(M-1,N)*LCV(M-1,N)
     &                                              )
        ELSE
	  FFL(M,N,K)= 0.0
        ENDIF
       ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDLAPLASV0(FF,FFL)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C GIVES LAPLAS ON UNIT GRID.
C FF  - INPUT VARIABLE ON LCV-MASK.
C FFL - OUTPUT LAPLAS OF FF ON LCV-MASK WITH ZERO FLUX ON BOUNDARY.
      REAL    FF(NX,NY,NZ),FFL(NX,NY,NZ)
      INTEGER M, N, K

!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO K = 1,NZ
      DO N=1,NY
       DO M=1,NX
        IF (LCV(M,N).GT.0.5) THEN

          FFL(M,N,K)=1.0/(DXH(M,N)*DYT(M,N)*HHV(M,N))
     &   * (
     &      ( FF(M,N+1,K)- FF(M,N  ,K))*DX(M,N+1)*DY(M,N+1)
     &       *HHQ(M,N+1)*LCV(M,N+1)
     &   -
     &      ( FF(M,N  ,K)- FF(M,N-1,K))*DX(M  ,N)*DY(M  ,N)
     &       *HHQ(M,N  )*LCV(M,N-1)
     &   +
     &      ( FF(M+1,N,K)- FF(M  ,N,K))*DXB(M  ,N)*DYB(M  ,N)
     &     * HH(M  ,N)*LCV(M+1,N)
     &   -
     &      ( FF(M  ,N,K)- FF(M-1,N,K))*DXB(M-1,N)*DYB(M-1,N)
     &     * HH(M-1,N)*LCV(M-1,N)
     &                                              )

        ELSE
	  FFL(M,N,K)= 0.0
        ENDIF
       ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDLAPLAS_T_FILTER(FF,VISCOEF)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C ONE STEP TIME INTEGRATION VISCOSITY OF 2-TH ORDER DIFFUSION
C EXPLICIT PROCEDURE.
C FF   - VARIABLE ON LU-MASK.
C FFL1 - LAPLAS 1ST ORDER OF FF ON LU MASK WITH ZERO FLUX ON BOUNDARY.
C FFL2 - LAPLAS 2ND ORDER OF FF ON LU-MASK WITH ZERO FLUX ON BOUNDARY.
C VISKOEF - UNDIMENSIONAL DIFFUSION  COEFFICIENT
      INTEGER K
      REAL FF(NX,NY,NZ),VISCOEF

      REAL FFL1(NX,NY,NZ)

      CALL GRIDLAPLAS0(FF,FFL1)
!$OMP PARALLEL DO
      DO K = 1,NZ
      FF(:,:,K)=FF(:,:,K)+VISCOEF*LU(:,:)*FFL1(:,:,K)
      ENDDO
!$OMP END PARALLEL DO
      IF(MMD.NE.0) CALL CYCLIZE(FF,NX,NY,NZ,MMM,MM)

      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDLAPLAS_U_FILTER(FF,VISCOEF)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C ONE STEP TIME INTEGRATION VISCOSITY OF 4-TH ORDER DIFFUSION
C EXPLICIT PROCEDURE.
C FF   - VARIABLE ON LCU-MASK.
C FFL1 - LAPLAS 1ST ORDER OF FF ON LCU MASK WITH ZERO ON U-BOUNDARY.
C VISKOEF - UNDIMENSIONAL DIFFUSION  COEFFICIENT
      INTEGER K
      REAL FF(NX,NY,NZ),VISCOEF


      REAL FFL1(NX,NY,NZ)

	CALL GRIDLAPLASU (FF,FFL1)

!$OMP PARALLEL DO
      DO K = 1,NZ
      FF(:,:,K)=FF(:,:,K)+VISCOEF*LCU(:,:)*FFL1(:,:,K)
      ENDDO
!$OMP END PARALLEL DO
      IF(MMD.NE.0) CALL CYCLIZE(FF,NX,NY,NZ,MMM,MM)

      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDLAPLAS_V_FILTER(FF,VISCOEF)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C ONE STEP TIME INTEGRATION VISCOSITY OF 4-TH ORDER DIFFUSION
C EXPLICIT PROCEDURE.
C FF   - VARIABLE ON LCV-MASK.
C FFL1 - LAPLAS 1ST ORDER OF FF ON LCV MASK WITH ZERO ON V-BOUNDARY.
C VISKOEF - UNDIMENSIONAL DIFFUSION  COEFFICIENT
      INTEGER K
      REAL FF(NX,NY,NZ),VISCOEF
      REAL FFL1(NX,NY,NZ)

      CALL GRIDLAPLASV(FF,FFL1)

!$OMP PARALLEL DO
      DO K = 1,NZ
      FF(:,:,K)=FF(:,:,K)+VISCOEF*LCV(:,:)*FFL1(:,:,K)
      ENDDO
!$OMP END PARALLEL DO
      IF(MMD.NE.0) CALL CYCLIZE(FF,NX,NY,NZ,MMM,MM)
      
      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE FILTER4_TGR_2D(FF,VISCOEF,MGRAD)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C ONE STEP TIME INTEGRATION VISCOSITY OF 4-TH ORDER DIFFUSION
C EXPLICIT PROCEDURE.
C FF   - VARIABLE ON LU-MASK.
C FFL1 - LAPLAS 1ST ORDER OF FF ON LU MASK WITH ZERO FLUX ON BOUNDARY.
C FFL2 - LAPLAS 2ND ORDER OF FF ON LU-MASK WITH ZERO FLUX ON BOUNDARY.
C VISKOEF - UNDIMENSIONAL DIFFUSION  COEFFICIENT
      INTEGER MGRAD,K
      REAL FF(NX,NY,MGRAD),VISCOEF

      REAL FFL1(NX,NY,MGRAD),FFL2(NX,NY,MGRAD)

      CALL GRDFLT0(FF,FFL1,MGRAD)
      IF(MMD.NE.0) CALL CYCLIZE(FFL1,NX,NY,MGRAD,MMM,MM)
      CALL GRDFLT0(FFL1,FFL2,MGRAD)

!$OMP PARALLEL DO
      DO K=1,MGRAD
      FF(:,:,K)=FF(:,:,K)-VISCOEF*LU(:,:)*FFL2(:,:,K)
      END DO
!$OMP END PARALLEL DO
      
      IF(MMD.NE.0) CALL CYCLIZE(FF,NX,NY,MGRAD,MMM,MM)

      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRDFLT0(FF,FFL,MGRAD)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C GIVES LAPLAS ON UNIT GRID.
C FF  - INPUT VARIABLE ON LU-MASK.
C FFL - OUTPUT LAPLAS OF FF ON LU-MASK WITH ZERO FLUX ON BOUNDARY.
      INTEGER MGRAD
      REAL    FF(NX,NY,MGRAD),FFL(NX,NY,MGRAD)
      INTEGER M, N ,K

      DO K=1,MGRAD
!$OMP PARALLEL DO PRIVATE(M,N)
      DO N=1,NY
       DO M=1,NX
        IF (LU(M,N).GT.0.5) THEN
        FFL(M,N,K)=1.0/(DX(M,N)*DY(M,N))
     &   * (
     &      ( FF(M+1,N,K)- FF(M  ,N,K))*DXT(M  ,N)*DYH(M  ,N)
     &     * LU(M+1,N)
     &   -
     &      ( FF(M  ,N,K)- FF(M-1,N,K))*DXT(M-1,N)*DYH(M-1,N)
     &     * LU(M-1,N)
     &   +
     &      ( FF(M,N+1,K)- FF(M,N  ,K))*DYT(M,N  )*DXH(M,N  )
     &     * LU(M,N+1)
     &   -
     &      ( FF(M,N  ,K)- FF(M,N-1,K))*DYT(M,N-1)*DXH(M,N-1)
     &     * LU(M,N-1)
     &                                              )
        ELSE
	  FFL(M,N,K)= 0.0
        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO
      ENDDO

      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE FILTER4_UGR_2D(FF,VISCOEF)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C ONE STEP TIME INTEGRATION VISCOSITY OF 4-TH ORDER DIFFUSION
C EXPLICIT PROCEDURE.
C FF   - VARIABLE ON LCU-MASK.
C FFL1 - LAPLAS 1ST ORDER OF FF ON LCU MASK WITH ZERO ON U-BOUNDARY.
C FFL2 - LAPLAS 2ND ORDER OF FF ON LCU-MASK WITH ZERO FLUX ON BOUNDARY.
C VISKOEF - UNDIMENSIONAL DIFFUSION  COEFFICIENT
      REAL FF(NX,NY),VISCOEF

      REAL FFL1(NX,NY),FFL2(NX,NY)

	CALL GRDFLTU (FF,FFL1)
      IF(MMD.NE.0) CALL CYCLIZE(FFL1,NX,NY,1,MMM,MM)
	CALL GRDFLTU0(FFL1,FFL2)

      FF(:,:)=FF(:,:)-VISCOEF*LCU(:,:)*FFL2(:,:)

      IF(MMD.NE.0) CALL CYCLIZE(FF,NX,NY,1,MMM,MM)

      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRDFLTU(FF,FFL)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C GIVES LAPLAS ON UNIT GRID.
C FF  - INPUT VARIABLE ON LCU-MASK.
C FFL - OUTPUT LAPLAS OF FF ON LCU-MASK WITH ZERO FLUX ON BOUNDARY.
      REAL    FF(NX,NY),FFL(NX,NY)
      INTEGER M, N

!$OMP PARALLEL DO PRIVATE(M,N)
      DO N=1,NY
       DO M=1,NX
        IF (LCU(M,N).GT.0.5) THEN

         FFL(M,N)=1.0/(DXT(M,N)*DYH(M,N))
     &   * (
     &      ( FF(M+1,N)- FF(M  ,N))*DX(M+1,N)*DY(M+1,N)
     &   -
     &      ( FF(M  ,N)- FF(M-1,N))*DX(M  ,N)*DY(M  ,N)
     &   +
     &      ( FF(M,N+1)- FF(M,N  ))*DXB(M,N  )*DYB(M,N  )
     &   -
     &      ( FF(M,N  )- FF(M,N-1))*DXB(M,N-1)*DYB(M,N-1)
     &                                              )
        ELSE
	  FFL(M,N)= 0.0
        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRDFLTU0(FF,FFL)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C GIVES LAPLAS ON UNIT GRID.
C FF  - INPUT VARIABLE ON LCU-MASK.
C FFL - OUTPUT LAPLAS OF FF ON LCU-MASK WITH ZERO FLUX ON BOUNDARY.
      REAL    FF(NX,NY),FFL(NX,NY)
      INTEGER M, N

!$OMP PARALLEL DO PRIVATE(M,N)
      DO N=1,NY
       DO M=1,NX
        IF (LCU(M,N).GT.0.5) THEN

         FFL(M,N)=1.0/(DXT(M,N)*DYH(M,N))
     &   * (
     &      ( FF(M+1,N)- FF(M  ,N))*DX(M+1,N)*DY(M+1,N)
     &     *LCU(M+1,N)
     &   -
     &      ( FF(M  ,N)- FF(M-1,N))*DX(M  ,N)*DY(M  ,N)
     &     *LCU(M-1,N)
     &   +
     &      ( FF(M,N+1)- FF(M,N  ))*DXB(M,N  )*DYB(M,N  )
     &     *LCU(M,N+1)
     &   -
     &      ( FF(M,N  )- FF(M,N-1))*DXB(M,N-1)*DYB(M,N-1)
     &     *LCU(M,N-1)
     &                                              )
        ELSE
	  FFL(M,N)= 0.0
        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE FILTER4_VGR_2D(FF,VISCOEF)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C ONE STEP TIME INTEGRATION VISCOSITY OF 4-TH ORDER DIFFUSION
C EXPLICIT PROCEDURE.
C FF   - VARIABLE ON LCV-MASK.
C FFL1 - LAPLAS 1ST ORDER OF FF ON LCV MASK WITH ZERO ON V-BOUNDARY.
C FFL2 - LAPLAS 2ND ORDER OF FF ON LCV-MASK WITH ZERO FLUX ON BOUNDARY.
C VISKOEF - UNDIMENSIONAL DIFFUSION  COEFFICIENT
      REAL FF(NX,NY),VISCOEF

      REAL FFL1(NX,NY),FFL2(NX,NY)

      CALL GRDFLTV(FF,FFL1)
      IF(MMD.NE.0) CALL CYCLIZE(FFL1,NX,NY,1,MMM,MM)
      CALL GRDFLTV0(FFL1,FFL2)

      FF(:,:)=FF(:,:)-VISCOEF*LCV(:,:)*FFL2(:,:)

      IF(MMD.NE.0) CALL CYCLIZE(FF,NX,NY,1,MMM,MM)
      
      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRDFLTV(FF,FFL)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C GIVES LAPLAS ON UNIT GRID.
C FF  - INPUT VARIABLE ON LCV-MASK.
C FFL - OUTPUT LAPLAS OF FF ON LCV-MASK WITH ZERO FLUX ON BOUNDARY.
      REAL    FF(NX,NY),FFL(NX,NY)
      INTEGER M, N, K

!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO N=1,NY
       DO M=1,NX
        IF (LCV(M,N).GT.0.5) THEN

          FFL(M,N)=1.0/(DXH(M,N)*DYT(M,N))
     &   * (
     &      ( FF(M,N+1)- FF(M,N  ))*DX(M,N+1)*DY(M,N+1)
     &   -
     &      ( FF(M,N  )- FF(M,N-1))*DX(M  ,N)*DY(M  ,N)
     &   +
     &      ( FF(M+1,N)- FF(M  ,N))*DXB(M  ,N)*DYB(M  ,N)
     &   -
     &      ( FF(M  ,N)- FF(M-1,N))*DXB(M-1,N)*DYB(M-1,N)
     &                                              )

        ELSE
	  FFL(M,N)= 0.0
        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRDFLTV0(FF,FFL)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C----------------------------------------------------------------------
C GIVES LAPLAS ON UNIT GRID.
C FF  - INPUT VARIABLE ON LCV-MASK.
C FFL - OUTPUT LAPLAS OF FF ON LCV-MASK WITH ZERO FLUX ON BOUNDARY.
      REAL    FF(NX,NY),FFL(NX,NY)
      INTEGER M, N, K

!$OMP PARALLEL DO PRIVATE(M,N,K)
      DO N=1,NY
       DO M=1,NX
        IF (LCV(M,N).GT.0.5) THEN

          FFL(M,N)=1.0/(DXH(M,N)*DYT(M,N))
     &   * (
     &      ( FF(M,N+1)- FF(M,N  ))*DX(M,N+1)*DY(M,N+1)
     &       *LCV(M,N+1)
     &   -
     &      ( FF(M,N  )- FF(M,N-1))*DX(M  ,N)*DY(M  ,N)
     &       *LCV(M,N-1)
     &   +
     &      ( FF(M+1,N)- FF(M  ,N))*DXB(M  ,N)*DYB(M  ,N)
     &       *LCV(M+1,N)
     &   -
     &      ( FF(M  ,N)- FF(M-1,N))*DXB(M-1,N)*DYB(M-1,N)
     &       *LCV(M-1,N)
     &                                              )

        ELSE
	  FFL(M,N)= 0.0
        ENDIF
       ENDDO
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDVISCOSITYU8(FF,VISCOEF)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C ONE STEP TIME INTEGRATION VISCOSITY OF 4-TH ORDER DIFFUSION
C EXPLICIT PROCEDURE.
C FF   - VARIABLE ON LCU-MASK.
C FFL1 - LAPLAS 1ST ORDER OF FF ON LCU MASK WITH ZERO ON U-BOUNDARY.
C FFL2 - LAPLAS 2ND ORDER OF FF ON LCU-MASK WITH ZERO FLUX ON BOUNDARY.
C VISKOEF - UNDIMENSIONAL DIFFUSION  COEFFICIENT
      INTEGER K
      REAL FF(NX,NY,NZ),VISCOEF


      REAL FFL1(NX,NY,NZ),FFL2(NX,NY,NZ)

	CALL GRIDLAPLASU (FF,FFL1)
      IF(MMD.NE.0) CALL CYCLIZE(FFL1,NX,NY,NZ,MMM,MM)
	
      CALL GRIDLAPLASU0(FFL1,FFL2)
      IF(MMD.NE.0) CALL CYCLIZE(FFL2,NX,NY,NZ,MMM,MM)

	CALL GRIDLAPLASU0(FFL2,FFL1)
      IF(MMD.NE.0) CALL CYCLIZE(FFL1,NX,NY,NZ,MMM,MM)
	
      CALL GRIDLAPLASU0(FFL1,FFL2)

!$OMP PARALLEL DO
      DO K = 1,NZ
      FF(:,:,K)=FF(:,:,K)-VISCOEF*LCU(:,:)*FFL2(:,:,K)
      ENDDO
!$OMP END PARALLEL DO
      IF(MMD.NE.0) CALL CYCLIZE(FF,NX,NY,NZ,MMM,MM)

      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDVISCOSITYV8(FF,VISCOEF)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C ONE STEP TIME INTEGRATION VISCOSITY OF 4-TH ORDER DIFFUSION
C EXPLICIT PROCEDURE.
C FF   - VARIABLE ON LCV-MASK.
C FFL1 - LAPLAS 1ST ORDER OF FF ON LCV MASK WITH ZERO ON V-BOUNDARY.
C FFL2 - LAPLAS 2ND ORDER OF FF ON LCV-MASK WITH ZERO FLUX ON BOUNDARY.
C VISKOEF - UNDIMENSIONAL DIFFUSION  COEFFICIENT
      INTEGER K
      REAL FF(NX,NY,NZ),VISCOEF


      REAL FFL1(NX,NY,NZ),FFL2(NX,NY,NZ)

      CALL GRIDLAPLASV(FF,FFL1)
      IF(MMD.NE.0) CALL CYCLIZE(FFL1,NX,NY,NZ,MMM,MM)
      
      CALL GRIDLAPLASV0(FFL1,FFL2)
      IF(MMD.NE.0) CALL CYCLIZE(FFL2,NX,NY,NZ,MMM,MM)
      
      CALL GRIDLAPLASV0(FFL2,FFL1)
      IF(MMD.NE.0) CALL CYCLIZE(FFL1,NX,NY,NZ,MMM,MM)
      
      CALL GRIDLAPLASV0(FFL1,FFL2)

!$OMP PARALLEL DO
      DO K = 1,NZ
      FF(:,:,K)=FF(:,:,K)-VISCOEF*LCV(:,:)*FFL2(:,:,K)
      ENDDO
!$OMP END PARALLEL DO
      IF(MMD.NE.0) CALL CYCLIZE(FF,NX,NY,NZ,MMM,MM)
      
      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDVISCOSITYU48(FF,VISCOEF)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C ONE STEP TIME INTEGRATION VISCOSITY OF 4-TH ORDER DIFFUSION
C EXPLICIT PROCEDURE.
C FF   - VARIABLE ON LCU-MASK.
C FFL1 - LAPLAS 1ST ORDER OF FF ON LCU MASK WITH ZERO ON U-BOUNDARY.
C FFL2 - LAPLAS 2ND ORDER OF FF ON LCU-MASK WITH ZERO FLUX ON BOUNDARY.
C VISKOEF - UNDIMENSIONAL DIFFUSION  COEFFICIENT
      INTEGER K
      REAL FF(NX,NY,NZ),VISCOEF


      REAL FFL1(NX,NY,NZ),FFL2(NX,NY,NZ),FFL3(NX,NY,NZ)

	CALL GRIDLAPLASU (FF,FFL1)
      IF(MMD.NE.0) CALL CYCLIZE(FFL1,NX,NY,NZ,MMM,MM)
	
      CALL GRIDLAPLASU0(FFL1,FFL2)
      IF(MMD.NE.0) CALL CYCLIZE(FFL2,NX,NY,NZ,MMM,MM)

	CALL GRIDLAPLASU0(FFL2,FFL3)
      IF(MMD.NE.0) CALL CYCLIZE(FFL3,NX,NY,NZ,MMM,MM)
	
      CALL GRIDLAPLASU0(FFL3,FFL1)

!$OMP PARALLEL DO
      DO K = 1,NZ
      FF(:,:,K)=FF(:,:,K)-VISCOEF*LCU(:,:)*(FFL1(:,:,K)+FFL2(:,:,K))
      ENDDO
!$OMP END PARALLEL DO
      IF(MMD.NE.0) CALL CYCLIZE(FF,NX,NY,NZ,MMM,MM)

      RETURN
      END
C======================================01.04.08 20:25(glavpip@pochta.ru)
      SUBROUTINE GRIDVISCOSITYV48(FF,VISCOEF)
      IMPLICIT NONE
	INCLUDE '0COM.INC'
C ONE STEP TIME INTEGRATION VISCOSITY OF 4-TH ORDER DIFFUSION
C EXPLICIT PROCEDURE.
C FF   - VARIABLE ON LCV-MASK.
C FFL1 - LAPLAS 1ST ORDER OF FF ON LCV MASK WITH ZERO ON V-BOUNDARY.
C FFL2 - LAPLAS 2ND ORDER OF FF ON LCV-MASK WITH ZERO FLUX ON BOUNDARY.
C VISKOEF - UNDIMENSIONAL DIFFUSION  COEFFICIENT
      INTEGER K
      REAL FF(NX,NY,NZ),VISCOEF


      REAL FFL1(NX,NY,NZ),FFL2(NX,NY,NZ),FFL3(NX,NY,NZ)

      CALL GRIDLAPLASV(FF,FFL1)
      IF(MMD.NE.0) CALL CYCLIZE(FFL1,NX,NY,NZ,MMM,MM)
      
      CALL GRIDLAPLASV0(FFL1,FFL2)
      IF(MMD.NE.0) CALL CYCLIZE(FFL2,NX,NY,NZ,MMM,MM)
      
      CALL GRIDLAPLASV0(FFL2,FFL3)
      IF(MMD.NE.0) CALL CYCLIZE(FFL3,NX,NY,NZ,MMM,MM)
      
      CALL GRIDLAPLASV0(FFL3,FFL1)

!$OMP PARALLEL DO
      DO K = 1,NZ
      FF(:,:,K)=FF(:,:,K)-VISCOEF*LCV(:,:)*(FFL1(:,:,K)+FFL2(:,:,K))
      ENDDO
!$OMP END PARALLEL DO
      IF(MMD.NE.0) CALL CYCLIZE(FF,NX,NY,NZ,MMM,MM)
      
      RETURN
      END
