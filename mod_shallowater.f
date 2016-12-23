! BASE 2D-3D-ARRAYS:
!---------------------------------------------------------------------72
!      COMMON /SEALEVEL/        SLH(NX,NY),  !SEA LEVEL HEIGHT
!     &                        SLH0(NX,NY),
!     &                     RBOTTOM(NX,NY),  !BOTTOM STRESS
!     &        SLHGRX(NX,NY),SLHGRY(NX,NY),  !SEA LEVEL X&Y GRADIENTS
!     &  UBRTR(NX,NY), VBRTR(NX,NY),      !BAROTROPIC VELOCITIES
!     &   DINX(NX,NY),  DINY(NX,NY)
!      REAL*8 SLH, RBOTTOM, SLHGRX, SLHGRY, SLH0
!      REAL*8 UBRTR, VBRTR
!      REAL*8 DINX, DINY
		module mod_shallowater
		implicit none

		real*8, allocatable :: slh(:,:), slh0(:,:)
		real*8, allocatable :: rbottom(:,:)
		real*8, allocatable :: slhgrx(:,:), slhgry(:,:)
		real*8, allocatable :: ubrtr(:,:), vbrtr(:,:)
		real*8, allocatable :: dinx(:,:), diny(:,:)

		contains

		subroutine allocate_shallowater()
		use mpi_parallel_tools
		implicit none

		allocate(slh(bnd_x1:bnd_x2, bnd_y1:bnd_y2))
		allocate(slh0(bnd_x1:bnd_x2, bnd_y1:bnd_y2))
        allocate(ubrtr(bnd_x1:bnd_x2, bnd_y1:bnd_y2))
        allocate(vbrtr(bnd_x1:bnd_x2, bnd_y1:bnd_y2))

        allocate(rbottom(nx_start:nx_end, ny_start:ny_end))
		allocate(slhgrx(nx_start:nx_end, ny_start:ny_end))
		allocate(slhgry(nx_start:nx_end, ny_start:ny_end))
		allocate(dinx(nx_start:nx_end, ny_start:ny_end))
		allocate(diny(nx_start:nx_end, ny_start:ny_end))

		end subroutine

		subroutine deallocate_shallowater()
		implicit none

		deallocate(slh)
		deallocate(slh0)
		deallocate(rbottom)
		deallocate(slhgrx)
		deallocate(slhgry)
		deallocate(ubrtr)
		deallocate(vbrtr)
		deallocate(dinx)
		deallocate(diny)

		end subroutine


		end module
