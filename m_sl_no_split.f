C================================================
C MODULES FOR SEA LEVEL CALCULATIONS PREPERED BY Rusakov.
c No splitting here. We form full matrice and solve it
C======================================================================
      module sea_level_no_split
#include <petsc/finclude/petscdef.h>

      use mpi_parallel_tools
      use mod_shallowater
      implicit none

c======================================================================
c      DECLARATIONS
c======================================================================
      INCLUDE '0COM.INC'
      INCLUDE '0SWEXPIMP.INC'

      PetscInt  Istart, Iend, its
      PetscReal norm
      Mat       matrix
      Vec       RHS
      Vec       SOL
      KSP       ksp

      integer nvar, locnozero
      integer, allocatable :: UINDX(:,:),VINDX(:,:),SLHINDX(:,:)

      real*8  :: A(7)
      integer :: posI(1), JA(7)

      real*8  :: valsol(3), valrhs(3), val(7)
      integer :: pos(3)

      real*8, pointer :: retarray(:)

      PetscInt nnz
c      integer nnz    ! number of nonzeroes

      INCLUDE '1SWEQPRECOND.INC'  !input parameters of solver and number for precondition

      integer maxnnz   ! size of A and JA

!
!     Arrays of indexes. connect two dimensional problem with one dimensional
!      vector excluding the ground.
!
	  INTEGER  IS_FORMED, MATRICE_PERMUTED,
     &                      PRECONDITIONER_FORMED

	  INTEGER, PARAMETER::DEBUG = 0, !set debug to 1 to see some useful information
     &                      DO_REORDERING_FOR_SPARSITY = 1

      INTEGER, SAVE :: factors(8), ! array of pointer for Super LU
     &                   isfactored  = 0

c======================================================================
c     END OF DECLARATIONS
c======================================================================

      contains

c======================================================================
c     SUBROUTINES
c======================================================================
      subroutine estimate_matrix_size()
      implicit none
      integer m, n, k, k1, k2, k3
      integer, allocatable :: nozero(:), offset(:)
      integer ierr
c----------------------------------------------------------------------
c
c      Estimate matrix size, Form numbering of variables.
c
c----------------------------------------------------------------------
      allocate(nozero(procs), offset(procs))
      locnozero = 0
      do n = max(ny_start,NNN),min(ny_end,NN)
          do m = max(nx_start,MMM),min(nx_end,MM)
              if ( LCU(M,N).GT.0.5 ) locnozero = locnozero + 1
              if ( LCV(M,N).GT.0.5 ) locnozero = locnozero + 1
              if ( LU (M,N).GT.0.5 ) locnozero = locnozero + 1
          enddo
      enddo
      call MPI_ALLGATHER(locnozero,1,MPI_INTEGER,
     &                        nozero,1,MPI_INTEGER,
     &                          CART_COMM, ierr)
      nvar = sum(nozero)
      offset(1) = 0
      do M = 1, procs-1
          offset(M+1) = offset(M) + nozero(M)
      enddo

c      print *, rank, locnozero, offset(rank+1), nvar

c     construct indexes
c
c      allocate(UINDX(NX, NY))
c      allocate(VINDX(NX, NY))
c      allocate(SLHINDX(NX, NY))
c      K = 0
c      do n = NNN, NN
c          do m = MMM, MM
      allocate(UINDX(bnd_x1:bnd_x2, bnd_y1:bnd_y2))
      allocate(VINDX(bnd_x1:bnd_x2, bnd_y1:bnd_y2))
      allocate(SLHINDX(bnd_x1:bnd_x2, bnd_y1:bnd_y2))
      K = offset(rank + 1)
      do n = max(ny_start,NNN),min(ny_end,NN)
          do m = max(nx_start,MMM),min(nx_end,MM)
              if ( LCU(M,N).GT.0.5 ) then
                  UINDX(M,N) = K
                  K = K + 1
              else
                  UINDX(M,N) = 0
              endif
              if ( LCV(M,N).GT.0.5 ) then
                  VINDX(M,N) = K
                  K = K + 1
              else
                  VINDX(M,N) = 0
              endif
              if ( LU(M,N).GT.0.5 ) then
                  SLHINDX(M,N) = K
                  K = K + 1
              else
                  SLHINDX(M,N) = 0
              endif
          enddo
      enddo

c      write(*,*) '|pcoord:',p_coord,'|  ', '|rank:',rank,'|  ',
c     &           '|nx:',nx_start,nx_end, '|  ',
c     &           '|ny:',ny_start,ny_end, '|  '

c      write(*,*) '|pcoord:',p_coord,'|  ', '|rank:',rank,'|  ',
c     &           '|b_x:',bnd_x1,bnd_x2, '|  ',
c     &           '|b_y:',bnd_y1,bnd_y2, '|  '

c      call MPI_FINALIZE(ierr)
c      stop

c
c      write(*,*)rank,'nx_str:',uindx(nx_start,ny_start:ny_end)
c      write(*,*)rank,'nx_end:',uindx(nx_end,ny_start:ny_end)
c      write(*,*)rank,'ny_str:',uindx(nx_start:nx_end,ny_start)
c      write(*,*)rank,'ny_end:',uindx(nx_start:nx_end,ny_end)

      call syncborder_int(UINDX)
      call syncborder_int(VINDX)
      call syncborder_int(SLHINDX)
c
c      write(*,*)rank,'bnd_x1:',uindx(bnd_x1,ny_start:ny_end)
c      write(*,*)rank,'bnd_x2:',uindx(bnd_x2,ny_start:ny_end)
c      write(*,*)rank,'bnd_y1:',uindx(nx_start:nx_end,bnd_y1)
c      write(*,*)rank,'bnd_y2:',uindx(nx_start:nx_end,bnd_y2)

c      call PetscFinalize(ierr)
c      stop

c      print *, VINDX(max(nx_start,MMM):min(nx_end,MM),
c     &                  max(ny_start,NNN):min(ny_end,NN))
c      print *, rank, max(nx_start,MMM),min(nx_end,MM),
c     &               max(ny_start,NNN),min(ny_end,NN)
c      call PetscFinalize(ierr)
c      stop
      end subroutine

c----------------------------------------------------------------------
c     Form matrice of sea level task
C
C    dU/dt + rb*U - l*V = mg d sl/dx
C
C    dU/dt + rb*V + l*U = ng d sl/dx
C
C             SL          d            d   n
C freesurface*-----  - [  --(HU) +    --(  - * HV ) ] = 0
C              dt         dx           dy  m
C
C    Implicit time scheme:
c
C           U-Uo                             d SLH
C           ----  + Rb * U - l * V = m*g*H * -----
C           tau                              d x
C
C
C           V-Vo                             d SLH
C           ----  + Rb * V + l * U = n*g*H * -----     (1)
C           tau                              d y
C
C             SLH - SLHo       d          d   n
C freesurface*----------  - [  --(U) +    --( - * V ) ] = 0
C              m*tau           dx         dy  m
C
C********************************************************************
C      remark:   freesurface =0 - rigid lid condition
C                freesurface =1 - long gravity waves are avaliable
c----------------------------------------------------------------------
c
C     Form matrice in compress sparse row format (CSR).
C    Use common ordering. the first is u, the second is v , the third is slh.
c    index m -changed first.
c
c  side effect
c     uncyclize operation on LCU and LCV
c
c----------------------------------------------------------------------
      SUBROUTINE FORM_MATRICE_SLH(tau, freesurface)
      IMPLICIT NONE
      integer :: M,N,K, k1, k2, k3
      real    :: tau, SLCU, SLCV
      real*8  :: freesurface
      real*8  :: frictau
      integer :: ileft, iright
      integer :: ierr

c----------------------------------------------------------------------
c
c      Form matrice
c
c---------------------------------------------------------------------

      PRECONDITIONER_FORMED = 0
      IS_FORMED = 0
      MATRICE_PERMUTED = 0

      call MatCreate(PETSC_COMM_WORLD,matrix,ierr)
      call MatSetSizes(matrix,locnozero,locnozero,nvar,nvar,ierr)
      call MatSetFromOptions(matrix,ierr)
      call MatSetUp(matrix,ierr)

c      call MatGetOwnershipRange(matrix,Istart,Iend,ierr)
c      print *, rank, Istart, Iend

c      call PetscFinalize(ierr)
c      stop


      do n = max(ny_start, NNN), min(ny_end, NN)
        do m = max(nx_start, MMM), min(nx_end, MM)
c----------------------------------------------------------------------
c         First equation
c----------------------------------------------------------------------
c         something like: see reports
c
C           U-Uo                             d SLH
C           ----  + Rb * U - l * V = m*g*H * -----
C           tau                              d x
c
c         coriolis on LUH grid
c----------------------------------------------------------------------
          nnz = 0
          ileft  = m-1
          iright = m+1
          if( MMD.ne.0 ) then
              if (m.eq.MMM) then
                  ileft  = MM
                  iright = m+1
              endif
              if (m.eq.MM) then   !cyclic right boundary
                  ileft  = m-1
                  iright = MMM
              end if
          end if !if cyclic condition

          if ( LCU(m,n).gt.0.5 ) then
            !diagonal
            posI(1) = UINDX(m, n)
            nnz = nnz + 1
            frictau = 1.0d0/dble(TAU) + RBOTTOM(M,N)
            A(nnz)  = frictau
	        JA(nnz) = UINDX(m,n)
c            IA(UINDX(m,n)) = nnz

            SLCV=4.0
            !v elements
	        if ( LCV(m,n).gt.0.5 ) then
		        nnz = nnz + 1
	            A(nnz)  = -dble(RLH(m,n)/SLCV)*SWEQ_IMP
	            JA(nnz) = VINDX(m,n)
	        endif
            if ( LCV(m+1,n).gt.0.5.and.LCV(iright,n).gt.0.5 ) then
	            nnz = nnz + 1
	            A(nnz)  = -dble(RLH(m,n)/SLCV)*SWEQ_IMP
	            JA(nnz) = VINDX(iright,n)
	        endif
	        if ( LCV(m,n-1).gt.0.5 ) then
		        nnz = nnz + 1
	            A(nnz)  = -dble(RLH(m,n-1)/SLCV)*SWEQ_IMP
	            JA(nnz) = VINDX(m,n-1)
	        endif
	        if ( LCV(m+1,n-1).gt.0.5.and.LCV(iright,n-1).gt.0.5 ) then
		        nnz = nnz + 1
	            A(nnz)  = -dble(RLH(m,n-1)/SLCV)*SWEQ_IMP
	            JA(nnz) = VINDX(iright,n-1)
	        endif
	        !slh elements
	        if ( LU(m,n).gt.0.5 ) then
		        nnz = nnz + 1
	            A(nnz)  = dble(GRV)/dble(DXT(M,N)*RN)*SWEQ_IMP
	            JA(nnz) = SLHINDX(m,n)
	        endif
	        if ( LU(m+1,n).gt.0.5 .and. LU(iright,n).gt.0.5 ) then
		        nnz = nnz + 1
	            A(nnz)  =-dble(GRV)/dble(DXT(M,N)*RN)*SWEQ_IMP
	            JA(nnz) = SLHINDX(iright,n)
	        endif
            call MatSetValues(matrix, 1, posI, nnz, JA, A,
     &                        INSERT_VALUES, ierr)
          endif
c----------------------------------------------------------------------
c         Second equation
c----------------------------------------------------------------------
c         something like: see reports
c
C           V-Vo                             d SLH
C           ----  + Rb * V + l * U = n*g*H * -----     (1)
C           tau                              d y
c
c         coriolis on LUH grid
c----------------------------------------------------------------------
          nnz = 0
          ileft  = m-1
          iright = m+1
          if( MMD.gt.0) then
              if (m.eq.MMM) then !cyclic left boundary
                  ileft  = MM
                  iright = m+1
              end if
              if (m.eq.MM) then   !cyclic right boundary
                  ileft  = m-1
                  iright = MMM
              end if
          end if !if cyclic condition

          if ( LCV(m,n).gt.0.5 ) then
            !diagonal
            posI(1) = VINDX(m, n)
            nnz = nnz + 1
            frictau = 1.0d0/dble(TAU) + RBOTTOM(M,N)
            A(nnz)  = frictau
	        JA(nnz) = VINDX(m,n)
c            IA(VINDX(m,n)) = nnz

            SLCU=4.0
            !u elements
	        if ( LCU(m,n).gt.0.5 ) then
	              nnz = nnz + 1
                  A(nnz)  =dble(RLH(m,n)/SLCU)*SWEQ_IMP !*dble(DXT(M,N))
!     &                       *dble( HH(M,N))
!     &                      /(dble(DXH(M,N))*dble(HHV(M,N)))
	              JA(nnz) = UINDX(m,n)
            endif
            if ( LCU(m-1,n).gt.0.5 ) then
              if ( LCU(ileft,n).gt.0.5 ) then
		          nnz = nnz + 1
	              A(nnz)  =dble(RLH(m-1,n)/SLCU)*SWEQ_IMP !*dble(DXT(M-1,N))
!     &                         *dble( HH(M-1,N))
!     &                        /(dble(DXH(M,N))*dble(HHV(M,N)))
	              JA(nnz) = UINDX(ileft,n)
              endif
	        endif
	        if ( LCU(m,n+1).gt.0.5 ) then
		        nnz = nnz + 1
                A(nnz)  =dble(RLH(m,n)/SLCU)*SWEQ_IMP !*dble(DXT(M,N+1))
!     &                         *dble( HH(M,N))
!     &                        /(dble(DXH(M,N))*dble(HHV(M,N)))
	            JA(nnz) = UINDX(m,n+1)
	        endif
	        if ( LCU(m-1,n+1).gt.0.5 ) then
              if ( LCU(ileft,n+1).gt.0.5 ) then
                nnz = nnz + 1
  	            A(nnz)  =dble(RLH(m-1,n)/SLCU)*SWEQ_IMP !*dble(DXT(M-1,N+1))
!     &                         *dble( HH(M-1,N))
!     &                        /(dble(DXH(M,N))*dble(HHV(M,N)))
	            JA(nnz) = UINDX(ileft,n+1)
             end if
	        endif

	        !slh elements
	        if ( LU(m,n).gt.0.5 ) then
		        nnz = nnz + 1
	            A(nnz)  = dble(GRV)/dble(DYT(M,N)*RN)*SWEQ_IMP
	            JA(nnz) = SLHINDX(m,n)
	        endif
	        if ( LU(m,n+1).gt.0.5 ) then
		        nnz = nnz + 1
	            A(nnz)  =-dble(GRV)/dble(DYT(M,N)*RN)*SWEQ_IMP
	            JA(nnz) = SLHINDX(m,n+1)
	        endif
            call MatSetValues(matrix, 1, posI, nnz, JA, A,
     &                        INSERT_VALUES, ierr)
          endif
c----------------------------------------------------------------------
c         Third equation
c----------------------------------------------------------------------
c         something like: see reports
c
C
C             SLH - SLHo                d           d    H
C  freesurface*----------  - m(ij)n ([  --(HU/n) +   --(  - * V ) ] = 0
C                  tau                 dx           dy   m
c
c         coriolis on LUH grid
c----------------------------------------------------------------------
          nnz = 0
          ileft  = m-1
          iright = m+1
          if( MMD.gt.0) then
              if (m.eq.MMM) then   !cyclic right boundary
                  ileft  = MM
                  iright = m+1
              end if
              if (m.eq.MM) then   !cyclic right boundary
                  ileft  = m-1
                  iright = MMM
              end if
          end if !if cyclic condition

          if ( LU(m,n).gt.0.5 ) then
            !diagonal
            posI(1) = SLHINDX(m, n)
	        nnz = nnz + 1
            A(nnz)  = 1.0d0/dble(TAU)*freesurface !the same have to be in RHS
	        JA(nnz) = SLHINDX(m,n)
c            IA(SLHINDX(m,n)) = nnz

            !u elements NAN
	        if ( LCU(m,n).gt.0.5 ) then
              nnz = nnz + 1
	          A(nnz)  = -dble(HHU(M,N))*dble(DYH(M,N))*SWEQ_IMP
     &                         /(dble(DX(M,N))*dble(DY(M,N)*RN))
	          JA(nnz) = UINDX(m,n)
	        endif
	        if ( LCU(m-1,n).gt.0.5.and.LCU(ileft,n).gt.0.5 ) then
              nnz = nnz + 1
	          A(nnz)  =  dble(HHU(M-1,N))*dble(DYH(M-1,N))*SWEQ_IMP
     &                         /(dble(DX(M,N))*dble(DY(M,N)*RN))
	          JA(nnz) = UINDX(ileft,n)
	        endif

            !v elements
	        if ( LCV(m,n).gt.0.5 ) then
              nnz = nnz + 1
	          A(nnz)  = -dble(HHV(M,N))*dble(DXH(M,N))*SWEQ_IMP
     &                         /(dble(DX(M,N))*dble(DY(M,N)*RN))
	          JA(nnz) = VINDX(m,n)
	        endif
	        if ( LCV(m,n-1).gt.0.5 ) then
              nnz = nnz + 1
	          A(nnz)  = dble(HHV(M,N-1))*dble(DXH(M,N-1))*SWEQ_IMP
     &                         /(dble(DX(M,N))*dble(DY(M,N)*RN))
	          JA(nnz) = VINDX(m,n-1)
	        endif
            call MatSetValues(matrix, 1, posI, nnz, JA, A,
     &                        INSERT_VALUES, ierr)
          endif
	    enddo
      enddo

      call MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)

c      call MatView(matrix, PETSC_VIEWER_STDOUT_WORLD)

c      IA(nvar+1) = nnz+1
      IS_FORMED = 1

c      call PetscFinalize(ierr)
c      stop

      call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
      call KSPSetOperators(ksp, matrix, matrix, ierr)
      call KSPSetFromOptions(ksp, ierr)
      call KSPSetInitialGuessNonZero(ksp, PETSC_TRUE, ierr)
c      call KSPSetType(ksp, KSPGMRES)

c      print *, 'INIIT KSP'
c      call KSPSetUp(ksp, ierr)
c      CALL KSPSETINITIALGUESSNONZERO(ksp,PETSC_TRUE,IERR)

      END SUBROUTINE

c----------------------------------------------------------------------
c     Form RHS
c     f1 - on LCU grid corresponds to equation for U
c     f2 - on LCV grid corresponds to equation for V
c     U0,V0,SLH0 are changed on exit ( 1/tau*[U0,V0,SLH0] are added)
c     RHS = [f1,f2,0] + 1/tau*[U0,V0, freesurface * SLH0]
c----------------------------------------------------------------------
      SUBROUTINE FORM_RHS(tau, freesurface)
      IMPLICIT NONE

      real    :: tau, SLCU, SLCV
      real*8  :: tau1, freesurface
	  integer :: m,n
      integer :: ierr

      tau1 = 1/tau

      call VecCreateMPI(PETSC_COMM_WORLD, locnozero, nvar,
     &                  SOL, ierr)
      call VecCreateMPI(PETSC_COMM_WORLD, locnozero, nvar,
     &                  RHS, ierr)
c
c      RHS:
c
      do n = max(ny_start, NNN), min(ny_end, NN)
        do m = max(nx_start, MMM), min(nx_end, MM)
         nnz = 0
C----------------------------- 1-ST EQUATION -----------------------------------
         if (LCU(M,N).GT.0.5) then
          nnz = nnz + 1
          pos(nnz) = UINDX(m,n)
          SLCV=4.0
	      valrhs(nnz) = dinx(M,N) + ubrtr(M,N)*tau1
     &       + (
     &           dble(RLH(m,n  )/SLCV)*(vbrtr(M,N  )+vbrtr(M+1,N ))
     &          +dble(RLH(m,n-1)/SLCV)*(vbrtr(M,N-1)+vbrtr(M+1,N-1))
     &          +dble(GRV)/dble(DXT(M,N)*RN)*(slh(M+1,N)-slh(M,N))
     &                                                      )*SWEQ_EXP
          valsol(nnz) = ubrtr(M,N)
         endif

C----------------------------- 2-ND EQUATION -----------------------------------
         if (LCV(M,N).GT.0.5) then
          nnz = nnz + 1
          pos(nnz) = VINDX(m,n)
          SLCU=4.0
	      valrhs(nnz) = diny(M,N) + vbrtr(M,N)*tau1
     &       + (
     &          -dble(RLH(m  ,n)/SLCU)*(ubrtr(M  ,N)+ubrtr(M  ,N+1))
     &          -dble(RLH(m-1,n)/SLCU)*(ubrtr(M-1,N)+ubrtr(M-1,N+1))
     &       +   dble(GRV)/dble(DYT(M,N)*RN)*(slh(M,N+1)-slh(M,N))
     &                                                      )*SWEQ_EXP
          valsol(nnz) = vbrtr(M,N)
         endif
C----------------------------- 3-RD EQUATION -----------------------------------
         if (LU(M,N).GT.0.5) then
          nnz = nnz + 1
          pos(nnz) = SLHINDX(m,n)
          valrhs(nnz) = slh(M,N) * tau1 * freesurface
     &  + (  dble(HHU(M  ,N  ))*dble(DYH(M  ,N  ))*ubrtr(M  ,N  )
     &      -dble(HHU(M-1,N  ))*dble(DYH(M-1,N  ))*ubrtr(M-1,N  )
     &      +dble(HHV(M  ,N  ))*dble(DXH(M  ,N  ))*vbrtr(M  ,N  )
     &      -dble(HHV(M  ,N-1))*dble(DXH(M  ,N-1))*vbrtr(M  ,N-1)
     &                  )  /(dble(DX(M,N))*dble(DY(M,N)*RN))*SWEQ_EXP
          valsol(nnz) = slh(M,N)
         endif

         call VecSetValues(RHS, nnz, pos, valrhs,
     &                     INSERT_VALUES, ierr)
         call VecSetValues(SOL, nnz, pos, valsol,
     &                     INSERT_VALUES, ierr)
        enddo
      enddo

      call VecAssemblyBegin(RHS, ierr)
      call VecAssemblyEnd(RHS, ierr)

      call VecAssemblyBegin(SOL, ierr)
      call VecAssemblyEnd(SOL, ierr)

c      call VecView(RHS, PETSC_VIEWER_STDOUT_WORLD)
c      call VecView(SOL, PETSC_VIEWER_STDOUT_WORLD)

      END SUBROUTINE
c--------------------------------------------------------------

      subroutine solve_system()
      implicit none
      integer :: ierr
      integer :: n, m, k
      real*8 :: w_time

      call start_timer(w_time)
      call KSPSolve(ksp, RHS, SOL, ierr)
      call end_timer(w_time)
c      if (rank.eq.0) print*,'  KSPSolve:',w_time

      call KSPGetIterationNumber(ksp, its, ierr)
      call KSPGetResidualNorm(ksp, norm, ierr)
c      call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)

      if (rank.eq.0) then
      print*,rank,'KSPSolve',w_time,'its',its,'norm',norm,'nvar',nvar
      endif
c      call VecView(SOL, PETSC_VIEWER_STDOUT_WORLD)
c      print *, "its:", its, "norm:", norm

      call start_timer(w_time)
      call VecGetArrayF90(SOL, retarray, ierr)
      k = 0
      do n = max(ny_start, NNN), min(ny_end, NN)
        do m = max(nx_start, MMM), min(nx_end, MM)
            if ( LCU(M,N).GT.0.5 ) then
                ubrtr(M,N) = retarray(k + 1)
                k = k + 1
            endif
            if ( LCV(M,N).GT.0.5 ) then
                vbrtr(M,N) = retarray(k + 1)
                k = k + 1
            endif
            if ( LU(M,N).GT.0.5 ) then
                slh(M,N) = retarray(k + 1)
                k = k + 1
            endif
        enddo
      enddo
      call VecRestoreArrayF90(SOL, retarray, ierr)
      call end_timer(w_time)
      if (rank.eq.0) print*,'  FORM u,v,slh:',w_time

      call start_timer(w_time)
      call syncborder_real8(ubrtr)
      call syncborder_real8(vbrtr)
      call syncborder_real8(slh)
      call end_timer(w_time)
      if (rank.eq.0) print*,'  SYNC:',w_time

      call VecDestroy(SOL, ierr)
      call VecDestroy(RHS, ierr)

      end subroutine

c--------------------------------------------------------------
c
      end module sea_level_no_split

c  Additional function
c-----------------------------------------------------------------------
      function distdot(n,x,ix,y,iy)
      integer n, ix, iy
      real*8 distdot, x(*), y(*), ddot
      external ddot
      distdot = ddot(n,x,ix,y,iy)
      return
c-----end-of-distdot
c-----------------------------------------------------------------------
      end function
C======================================================================
