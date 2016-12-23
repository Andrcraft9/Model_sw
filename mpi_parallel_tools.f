!---------------------module for definition of array dimensions and boundaries-----------------
		module mpi_parallel_tools
#include <petsc/finclude/petscdef.h>

        use petscksp
        implicit none

c        include 'mpif.h'

		integer :: nx_start !first significant point in x-direction
		integer :: nx_end   !last  significant point in x-direction
		integer :: ny_start !first significant point in y-direction
		integer :: ny_end   !last  significant point in y-direction

		integer :: bnd_x1   !left   array boundary in x-direction
		integer :: bnd_x2   !right  array boundary in x-direction
		integer :: bnd_y1   !bottom array boundary in y-direction
		integer :: bnd_y2   !top    array boundary in y-direction

        integer :: rank, procs
        integer :: CART_COMM
        integer, dimension(2) :: p_size, period, p_coord


        contains

        subroutine start_timer(time)
        implicit none

        include 'mpif.h'

        real*8, intent(inout) :: time
        integer :: ierr

        time = MPI_wtime(ierr)
        return
        end subroutine

        subroutine end_timer(time)
        implicit none

        include 'mpif.h'

        real*8, intent(inout) :: time
        real*8 :: outtime
        integer :: ierr

        time = MPI_wtime(ierr) - time
c        call MPI_allreduce(time, outtime, 1, MPI_REAL8,
c     &                     MPI_MAX, CART_COMM, ierr)
c        time = outtime
        return
        end subroutine

        integer function check_p_coord(coord)
            implicit none
            integer, dimension(2), intent(in) :: coord

            check_p_coord = 0

c            write(*,*) coord,all(coord.ge.0),all((p_size-coord).ge.1)
c            print *, coord, p_size - coord, all((p_size-coord).ge.1)

            if (all(coord.ge.0) .and. all((p_size-coord).ge.1)) then
                check_p_coord = 1
            endif
            return
        end function

c-------------------------------------------------------------------------------
        subroutine syncborder_int(fild)
            implicit none
            integer :: fild(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
            integer, dimension(2) :: p_dist, p_src
            integer :: dist_rank, src_rank
            integer :: flag_dist, flag_src
            integer :: ierr
            integer stat(MPI_status_size)

c------------------ send-recv in NY+ -------------------------------------------
            p_dist(1) = p_coord(1)
            p_dist(2) = p_coord(2) + 1
            p_src(1) = p_coord(1)
            p_src(2) = p_coord(2) - 1
            flag_dist = check_p_coord(p_dist)
            flag_src = check_p_coord(p_src)

            if ( (flag_src .eq. 1) .and. (flag_dist .eq. 1) ) then
                call MPI_cart_rank(CART_COMM, p_dist,dist_rank,ierr)
                call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)

                call MPI_sendrecv(fild(nx_start:nx_end, ny_end),
     &                       nx_end - nx_start + 1,
     &                       MPI_integer, dist_rank, 1,
     &                       fild(nx_start:nx_end, bnd_y1),
     &                       nx_end - nx_start + 1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                write(*,*)rank,p_coord,'NY+|s-r',ierr
            else
                if (flag_src .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_src,src_rank,ierr)

                    call MPI_recv(fild(nx_start:nx_end, bnd_y1),
     &                       nx_end - nx_start + 1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,src_rank,'NY+|r',ierr
                endif
                if (flag_dist .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_dist,dist_rank,ierr)

                    call MPI_send(fild(nx_start:nx_end, ny_end),
     &                       nx_end - nx_start + 1,
     &                       MPI_integer, dist_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,dist_rank,'NY+|s',ierr
                endif
            endif
c-------------------------------------------------------------------------------

c------------------ send-recv in NX+ -------------------------------------------
            p_dist(1) = p_coord(1) + 1
            p_dist(2) = p_coord(2)
            p_src(1) = p_coord(1) - 1
            p_src(2) = p_coord(2)

            flag_dist = check_p_coord(p_dist)
            flag_src = check_p_coord(p_src)

            if ( (flag_src .eq. 1) .and. (flag_dist .eq. 1) ) then
                call MPI_cart_rank(CART_COMM, p_dist,dist_rank,ierr)
                call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)

                call MPI_sendrecv(fild(nx_end, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_integer, dist_rank, 1,
     &                       fild(bnd_x1, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                write(*,*)rank,p_coord,'NX+|s-r', ierr
            else
                if (flag_src .eq. 1) then
                    call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)

                    call MPI_recv(fild(bnd_x1, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,src_rank,'NX+|r',ierr
                endif
                if (flag_dist .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_dist,dist_rank,ierr)

                    call MPI_send(fild(nx_end, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_integer, dist_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,dist_rank,'NX+|s',ierr
                endif
            endif
c-------------------------------------------------------------------------------

c------------------ send-recv in NY- -------------------------------------------
            p_dist(1) = p_coord(1)
            p_dist(2) = p_coord(2) - 1
            p_src(1) = p_coord(1)
            p_src(2) = p_coord(2) + 1
            flag_dist = check_p_coord(p_dist)
            flag_src = check_p_coord(p_src)

            if ( (flag_src .eq. 1) .and. (flag_dist .eq. 1) ) then
                call MPI_cart_rank(CART_COMM, p_dist,dist_rank,ierr)
                call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)

                call MPI_sendrecv(fild(nx_start:nx_end, ny_start),
     &                       nx_end - nx_start + 1,
     &                       MPI_integer, dist_rank, 1,
     &                       fild(nx_start:nx_end, bnd_y2),
     &                       nx_end - nx_start + 1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                write(*,*)rank,p_coord,'NY-|s-r',ierr
            else
                if (flag_src .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_src,src_rank,ierr)

                    call MPI_recv(fild(nx_start:nx_end, bnd_y2),
     &                       nx_end - nx_start + 1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,src_rank,'NY-|r',ierr
                endif
                if (flag_dist .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_dist,dist_rank,ierr)

                    call MPI_send(fild(nx_start:nx_end, ny_start),
     &                       nx_end - nx_start + 1,
     &                       MPI_integer, dist_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,dist_rank,'NY-|s',ierr
                endif
            endif
c-------------------------------------------------------------------------------

c------------------ send-recv in NX- -------------------------------------------
            p_dist(1) = p_coord(1) - 1
            p_dist(2) = p_coord(2)
            p_src(1) = p_coord(1) + 1
            p_src(2) = p_coord(2)
            flag_dist = check_p_coord(p_dist)
            flag_src = check_p_coord(p_src)

            if ( (flag_src .eq. 1) .and. (flag_dist .eq. 1) ) then
                call MPI_cart_rank(CART_COMM, p_dist,dist_rank,ierr)
                call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)

                call MPI_sendrecv(fild(nx_start, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_integer, dist_rank, 1,
     &                       fild(bnd_x2, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                write(*,*)rank,p_coord,'NX-|s-r', ierr
            else
                if (flag_src .eq. 1) then
                    call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)

                    call MPI_recv(fild(bnd_x2, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,src_rank,'NX-|r',ierr
                endif
                if (flag_dist .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_dist,dist_rank,ierr)

                    call MPI_send(fild(nx_start, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_integer, dist_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,dist_rank,'NX-|s',ierr
                endif
            endif
c-------------------------------------------------------------------------------

c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
c------------------ Sync edge points (EP) --------------------------------------
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------

c------------------ send-recv EP in NX+,NY+ ------------------------------------
            p_dist(1) = p_coord(1) + 1
            p_dist(2) = p_coord(2) + 1
            p_src(1) = p_coord(1) - 1
            p_src(2) = p_coord(2) - 1
            flag_dist = check_p_coord(p_dist)
            flag_src = check_p_coord(p_src)
            if ( (flag_src .eq. 1) .and. (flag_dist .eq. 1) ) then
                call MPI_cart_rank(CART_COMM, p_dist,dist_rank,ierr)
                call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)
                call MPI_sendrecv(fild(nx_end, ny_end),
     &                       1,
     &                       MPI_integer, dist_rank, 1,
     &                       fild(bnd_x1, bnd_y1),
     &                       1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                write(*,*)rank,p_coord,'NY+|s-r',ierr
            else
                if (flag_src .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_src,src_rank,ierr)
                    call MPI_recv(fild(bnd_x1, bnd_y1),
     &                       1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,src_rank,'NY+|r',ierr
                endif
                if (flag_dist .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_dist,dist_rank,ierr)
                    call MPI_send(fild(nx_end, ny_end),
     &                       1,
     &                       MPI_integer, dist_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,dist_rank,'NY+|s',ierr
                endif
            endif
c-------------------------------------------------------------------------------

c------------------ send-recv EP in NX+,NY- ------------------------------------
            p_dist(1) = p_coord(1) + 1
            p_dist(2) = p_coord(2) - 1
            p_src(1) = p_coord(1) - 1
            p_src(2) = p_coord(2) + 1
            flag_dist = check_p_coord(p_dist)
            flag_src = check_p_coord(p_src)
            if ( (flag_src .eq. 1) .and. (flag_dist .eq. 1) ) then
                call MPI_cart_rank(CART_COMM, p_dist,dist_rank,ierr)
                call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)
                call MPI_sendrecv(fild(nx_end, ny_start),
     &                       1,
     &                       MPI_integer, dist_rank, 1,
     &                       fild(bnd_x1, bnd_y2),
     &                       1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                write(*,*)rank,p_coord,'NX+|s-r', ierr
            else
                if (flag_src .eq. 1) then
                    call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)
                    call MPI_recv(fild(bnd_x1, bnd_y2),
     &                       1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,src_rank,'NX+|r',ierr
                endif
                if (flag_dist .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_dist,dist_rank,ierr)
                    call MPI_send(fild(nx_end, ny_start),
     &                       1,
     &                       MPI_integer, dist_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,dist_rank,'NX+|s',ierr
                endif
            endif
c-------------------------------------------------------------------------------

c------------------ send-recv EP in NX-,NY- ------------------------------------
            p_dist(1) = p_coord(1) - 1
            p_dist(2) = p_coord(2) - 1
            p_src(1) = p_coord(1) + 1
            p_src(2) = p_coord(2) + 1
            flag_dist = check_p_coord(p_dist)
            flag_src = check_p_coord(p_src)

            if ( (flag_src .eq. 1) .and. (flag_dist .eq. 1) ) then
                call MPI_cart_rank(CART_COMM, p_dist,dist_rank,ierr)
                call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)
                call MPI_sendrecv(fild(nx_start, ny_start),
     &                       1,
     &                       MPI_integer, dist_rank, 1,
     &                       fild(bnd_x2, bnd_y2),
     &                       1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                write(*,*)rank,p_coord,'NY-|s-r',ierr
            else
                if (flag_src .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_src,src_rank,ierr)
                    call MPI_recv(fild(bnd_x2, bnd_y2),
     &                       1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,src_rank,'NY-|r',ierr
                endif
                if (flag_dist .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_dist,dist_rank,ierr)
                    call MPI_send(fild(nx_start, ny_start),
     &                       1,
     &                       MPI_integer, dist_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,dist_rank,'NY-|s',ierr
                endif
            endif
c-------------------------------------------------------------------------------

c------------------ send-recv EP in NX-,NY+ ------------------------------------
            p_dist(1) = p_coord(1) - 1
            p_dist(2) = p_coord(2) + 1
            p_src(1) = p_coord(1) + 1
            p_src(2) = p_coord(2) - 1
            flag_dist = check_p_coord(p_dist)
            flag_src = check_p_coord(p_src)
            if ( (flag_src .eq. 1) .and. (flag_dist .eq. 1) ) then
                call MPI_cart_rank(CART_COMM, p_dist,dist_rank,ierr)
                call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)
                call MPI_sendrecv(fild(nx_start, ny_end),
     &                       1,
     &                       MPI_integer, dist_rank, 1,
     &                       fild(bnd_x2, bnd_y1),
     &                       1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                write(*,*)rank,p_coord,'NX-|s-r', ierr
            else
                if (flag_src .eq. 1) then
                    call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)
                    call MPI_recv(fild(bnd_x2, bnd_y1),
     &                       1,
     &                       MPI_integer, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,src_rank,'NX-|r',ierr
                endif
                if (flag_dist .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_dist,dist_rank,ierr)
                    call MPI_send(fild(nx_start, ny_end),
     &                       1,
     &                       MPI_integer, dist_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,dist_rank,'NX-|s',ierr
                endif
            endif
c-------------------------------------------------------------------------------
            return
        end subroutine



c-------------------------------------------------------------------------------
        subroutine syncborder_real8(fild)
            implicit none
            real*8 :: fild(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
            integer, dimension(2) :: p_dist, p_src
            integer :: dist_rank, src_rank
            integer :: flag_dist, flag_src
            integer :: ierr, debg
            integer stat(MPI_status_size)

            debg = 0
c------------------ send-recv in NY+ -------------------------------------------
            p_dist(1) = p_coord(1)
            p_dist(2) = p_coord(2) + 1
            p_src(1) = p_coord(1)
            p_src(2) = p_coord(2) - 1
            flag_dist = check_p_coord(p_dist)
            flag_src = check_p_coord(p_src)

            if ( (flag_src .eq. 1) .and. (flag_dist .eq. 1) ) then
                call MPI_cart_rank(CART_COMM, p_dist,dist_rank,ierr)
                call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)

                call MPI_sendrecv(fild(nx_start:nx_end, ny_end),
     &                       nx_end - nx_start + 1,
     &                       MPI_real8, dist_rank, 1,
     &                       fild(nx_start:nx_end, bnd_y1),
     &                       nx_end - nx_start + 1,
     &                       MPI_real8, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                write(*,*)rank,p_coord,'NY+|s-r',ierr
            else
                if (flag_src .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_src,src_rank,ierr)

                    call MPI_recv(fild(nx_start:nx_end, bnd_y1),
     &                       nx_end - nx_start + 1,
     &                       MPI_real8, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,src_rank,'NY+|r',ierr
                endif

                if (flag_dist .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_dist,dist_rank,ierr)

                    call MPI_send(fild(nx_start:nx_end, ny_end),
     &                       nx_end - nx_start + 1,
     &                       MPI_real8, dist_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,dist_rank,'NY+|s',ierr
                endif
            endif
c-------------------------------------------------------------------------------

c------------------ send-recv in NX+ -------------------------------------------
            p_dist(1) = p_coord(1) + 1
            p_dist(2) = p_coord(2)
            p_src(1) = p_coord(1) - 1
            p_src(2) = p_coord(2)
            flag_dist = check_p_coord(p_dist)
            flag_src = check_p_coord(p_src)

            if ( (flag_src .eq. 1) .and. (flag_dist .eq. 1) ) then
                call MPI_cart_rank(CART_COMM, p_dist,dist_rank,ierr)
                call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)

                call MPI_sendrecv(fild(nx_end, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_real8, dist_rank, 1,
     &                       fild(bnd_x1, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_real8, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                write(*,*)rank,p_coord,'NX+|s-r', ierr
            else
                if (flag_src .eq. 1) then
                    call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)

                    call MPI_recv(fild(bnd_x1, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_real8, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,src_rank,'NX+|r',ierr
                endif

                if (flag_dist .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_dist,dist_rank,ierr)

                    call MPI_send(fild(nx_end, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_real8, dist_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,dist_rank,'NX+|s',ierr
                endif
            endif
c-------------------------------------------------------------------------------

c------------------ send-recv in NY- -------------------------------------------
            p_dist(1) = p_coord(1)
            p_dist(2) = p_coord(2) - 1
            p_src(1) = p_coord(1)
            p_src(2) = p_coord(2) + 1
            flag_dist = check_p_coord(p_dist)
            flag_src = check_p_coord(p_src)

            if ( (flag_src .eq. 1) .and. (flag_dist .eq. 1) ) then
                call MPI_cart_rank(CART_COMM, p_dist,dist_rank,ierr)
                call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)

                call MPI_sendrecv(fild(nx_start:nx_end, ny_start),
     &                       nx_end - nx_start + 1,
     &                       MPI_real8, dist_rank, 1,
     &                       fild(nx_start:nx_end, bnd_y2),
     &                       nx_end - nx_start + 1,
     &                       MPI_real8, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                write(*,*)rank,p_coord,'NY-|s-r',ierr
            else
                if (flag_src .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_src,src_rank,ierr)

                    call MPI_recv(fild(nx_start:nx_end, bnd_y2),
     &                       nx_end - nx_start + 1,
     &                       MPI_real8, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,src_rank,'NY-|r',ierr
                endif

                if (flag_dist .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_dist,dist_rank,ierr)

                    call MPI_send(fild(nx_start:nx_end, ny_start),
     &                       nx_end - nx_start + 1,
     &                       MPI_real8, dist_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,dist_rank,'NY-|s',ierr
                endif
            endif
c-------------------------------------------------------------------------------

c------------------ send-recv in NX- -------------------------------------------
            p_dist(1) = p_coord(1) - 1
            p_dist(2) = p_coord(2)
            p_src(1) = p_coord(1) + 1
            p_src(2) = p_coord(2)
            flag_dist = check_p_coord(p_dist)
            flag_src = check_p_coord(p_src)

            if ( (flag_src .eq. 1) .and. (flag_dist .eq. 1) ) then
                call MPI_cart_rank(CART_COMM, p_dist,dist_rank,ierr)
                call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)

                call MPI_sendrecv(fild(nx_start, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_real8, dist_rank, 1,
     &                       fild(bnd_x2, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_real8, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                write(*,*)rank,p_coord,'NX-|s-r', ierr
            else
                if (flag_src .eq. 1) then
                    call MPI_cart_rank(CART_COMM, p_src, src_rank, ierr)

                    call MPI_recv(fild(bnd_x2, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_real8, src_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,src_rank,'NX-|r',ierr
                endif

                if (flag_dist .eq. 1) then
                    call MPI_cart_rank(CART_COMM,p_dist,dist_rank,ierr)

                    call MPI_send(fild(nx_start, ny_start:ny_end),
     &                       ny_end - ny_start + 1,
     &                       MPI_real8, dist_rank, 1,
     &                       CART_COMM, stat, ierr)
c                    write(*,*)rank,p_coord,dist_rank,'NX-|s',ierr
                endif
            endif
c-------------------------------------------------------------------------------
            return
        end subroutine
		endmodule mpi_parallel_tools
!---------------------end module for definition of array dimensions and boundaries-----------------
