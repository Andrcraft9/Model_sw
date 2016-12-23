C======================================================================
ccc  File m_bottomfriction.f
ccc  Contains:
ccc        BOTTOM_FRICTION. - bottom friction parametrization.
ccc
ccc  Description: 
ccc        Bottom friction parametrization. Update bottom_friction array
ccc
ccc  Arguments:
ccc        UU, VV       - gorizontal velocities
ccc  
ccc  Other Parameter are stored in 0BOTTFR.INC file
ccc See: OPA 8.1 Ocean General Circulation Model Reference Manual
C======================================================================
      SUBROUTINE BOTTOM_FRICTION(UU,VV)
      IMPLICIT NONE
C----------------------------------------------------------------------
C VERTICAL TRANSPORT AND 3D DIFFUSION OF FF ON U-GRID.
      INCLUDE '0COM.INC'
      INCLUDE '0BOTTFR.INC'
      REAL UU(NX,NY,NZ),VV(NX,NY,NZ)
      INTEGER M, N
      REAL    bu, bv, dz1
      REAL, PARAMETER :: little = 1e-6

      dz1 = 1.0/DZ(NZ)
C
C 0. bottom friction is zero
C ------------
C
      IF ( type_fric.EQ.0 ) THEN
          BottomFriction = 0.0  
       
C 1. linear bottom friction
C -------------------------
C
      ELSEIF ( type_fric.EQ.1 ) THEN
          BottomFriction = Cb_l * dz1
         
C
C 2. nonlinear bottom friction
C ----------------------------
C
      ELSEIF ( type_fric.EQ.2 ) THEN
!$OMP PARALLEL DO PRIVATE(M,N,BU,BV)          
          DO N=NNN,NN
              DO M=MMM,MM
c                  bv  = ( LCV(M,N)   * VV(M,N,NZ)   ! on temperature mask
c     &                  + LCV(M,N-1) * VV(M,N-1,NZ)) /
c     &	                (LCV(M,N) + LCV(M,N-1) + little)
c                  bu  = ( LCU(M,N)   * UU(M,N,NZ) 
c     &                  + LCU(M-1,N) * UU(M-1,N,NZ)) /
c     &	                (LCU(M,N) + LCU(M-1,N) + little)
C  Another means of avaraging on temperature mask
                  bv  = (VV(M,N,NZ) + VV(M,N-1,NZ)) /2.0
                  bu  = (UU(M,N,NZ) + UU(M-1,N,NZ)) /2.0
                  BottomFriction(M,N) = Cb_nl
     &                            * sqrt(bu*bu + bv*bv+Ebottom) * dz1
              ENDDO
          ENDDO
!$OMP END PARALLEL DO
C
C 3. E r r o r
C ------------    
C
      ELSE
          STOP 'bottomfriction'
      ENDIF !type_fric

      IF(MMD.NE.0) CALL CYCLIZE(BottomFriction,NX,NY,1,MMM,MM)
	
      RETURN
      END
C======================================================================
