  module mpi_mod

    use, intrinsic :: iso_fortran_env, only : OUTPUT_UNIT

    implicit none

    private

    public :: Calc_tasks_in_x_and_y, Calc_patch_dims, Gather_var2d, Do_halo_exchange, Do_halo_exchange_with_corners, &
        Max_across_mpi_tasks, Sum_across_mpi_tasks, Distribute_var2d, Min_across_mpi_tasks, Convert_mpi_comm_to_f08, &
        Print_cart_info, topology_dim_order

    ! --- Topology Ordering Flag ---
    ! 0 = Default (Fortran order: X is Dim 0, Y is Dim 1)
    ! 1 = External (C/C++ order: Y is Dim 0, X is Dim 1)
    integer :: topology_dim_order = 0
    ! ------------------------------

  contains

    subroutine Calc_patch_dims (nx, ny, px, py, coords, istart, iend, jstart, jend)

      implicit none

      integer, intent (in) :: nx, ny, px, py
      integer, dimension(2) :: coords
      integer, intent (out) :: istart, iend, jstart, jend

      integer :: nx_local, ny_local, x, y, rx, ry


      x = coords(1)
      y = coords(2)

      rx = mod(nx, px)
      ry = mod(ny, py)

        ! size patch
      nx_local = nx / px
      if (x < rx) nx_local = nx_local + 1

      ny_local = ny / py
      if (y < ry) ny_local = ny_local + 1

        ! coords in the grid
      istart = (nx / px) * x + min(x, rx) + 1
      iend = istart + nx_local - 1

      jstart = (ny / py) * y + min(y, ry) + 1
      jend = jstart + ny_local - 1

    end subroutine Calc_patch_dims

    pure subroutine Calc_tasks_in_x_and_y (ntasks, nx, ny, px, py)

      implicit none

      integer, intent (in)  :: ntasks, nx, ny
      integer, intent (out) :: px, py

      integer :: i, j
      real    :: best_ratio, this_ratio, target_ratio,ratio


      best_ratio = huge (best_ratio)
      target_ratio = real (nx) / real (ny)

      px = 1
      py = ntasks

      do i = 1, ntasks
        if (mod (ntasks, i) == 0) then
          j = ntasks / i
          ratio = real (i) / real (j)
          this_ratio = abs (ratio - target_ratio)

          if (this_ratio < best_ratio) then
            best_ratio = this_ratio
            px = i
            py = j
          end if
        end if
      end do

    end subroutine Calc_tasks_in_x_and_y

    subroutine Convert_mpi_comm_to_f08 (comm_old, comm_new)

#ifdef DM_PARALLEL
      use mpi_f08
#endif

      implicit none

      integer, intent (in) :: comm_old
#ifdef DM_PARALLEL
      type(MPI_Comm), intent (out) :: comm_new
#else
      integer, intent (out) :: comm_new
#endif


      comm_new = transfer (comm_old, comm_new)

    end subroutine Convert_mpi_comm_to_f08

    subroutine Distribute_var2d (local_field, ips, ipe, jps, jpe, comm)

      ! Rank 0 has global_field in local_field and it is distributed across tasks. Rank 0 retains only its part

#ifdef DM_PARALLEL
      use mpi
#endif

      implicit none

      real, dimension(:, :), allocatable, intent(in out) :: local_field
      integer, intent(in) :: ips, ipe, jps, jpe
      integer, intent(in) :: comm

      integer :: ierr, rank, nprocs, nx_loc, ny_loc
      integer, allocatable :: all_ips(:), all_ipe(:), all_jps(:), all_jpe(:)
      integer, allocatable :: sendcounts(:), displs(:)
      real, allocatable :: sendbuf(:)
      integer :: r, i, j, k, kbuf
      logical :: am_root

      logical, parameter :: DEBUG_LOCAL = .false.
      real, dimension(:, :), allocatable :: tmp


      if (DEBUG_LOCAL) write (OUTPUT_UNIT, *) 'Entering Distribute_var2d...'

#ifdef DM_PARALLEL
      call MPI_Comm_rank (comm, rank, ierr)
      call MPI_Comm_size (comm, nprocs, ierr)

      am_root = (rank == 0)

      nx_loc = ipe - ips + 1
      ny_loc = jpe - jps + 1

        ! Gather indices on rank 0
      if (rank == 0) then
        allocate (all_ips(nprocs), all_ipe(nprocs))
        allocate (all_jps(nprocs), all_jpe(nprocs))
      end if

      call MPI_Gather(ips, 1, MPI_INTEGER, all_ips, 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_Gather(ipe, 1, MPI_INTEGER, all_ipe, 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_Gather(jps, 1, MPI_INTEGER, all_jps, 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_Gather(jpe, 1, MPI_INTEGER, all_jpe, 1, MPI_INTEGER, 0, comm, ierr)

        ! Rank 0 prepares sendcounts, displacements, and buffer
      if (rank == 0) then
        allocate(sendcounts(nprocs), displs(nprocs))

        do r = 1, nprocs
          sendcounts(r) = (all_ipe(r) - all_ips(r) + 1) * (all_jpe(r) - all_jps(r) + 1)
        end do

        displs(1) = 0
        do r = 2, nprocs
          displs(r) = displs(r - 1) + sendcounts(r - 1)
        end do

        allocate(sendbuf(sum(sendcounts)))
        kbuf = 0
        do r = 1, nprocs
          do j = all_jps(r), all_jpe(r)
            do i = all_ips(r), all_ipe(r)
              kbuf = kbuf + 1
              sendbuf(kbuf) = local_field(i - ips + 1, j - jps + 1)
            end do
          end do
        end do

        allocate (tmp(nx_loc, ny_loc))
        tmp(:, :) = local_field(ips:ipe, jps:jpe)
      end if

      if (rank /= 0) then
        if (allocated (local_field)) deallocate (local_field)
        allocate (local_field(nx_loc, ny_loc))
      end if

        ! Send data to all ranks from rank 0
      call MPI_Scatterv(sendbuf, sendcounts, displs, MPI_REAL, &
          local_field, nx_loc*ny_loc, MPI_REAL, 0, comm, ierr)

      if (rank == 0) then
        if (allocated (sendbuf)) deallocate (sendbuf)
        if (allocated (sendcounts)) deallocate (sendcounts)
        if (allocated (displs)) deallocate (displs)
        if (allocated (all_ips)) deallocate (all_ips, all_ipe, all_jps, all_jpe)

          ! Going from global to local on Rank0
        local_field = tmp
        deallocate (tmp)
      end if
#endif

      if (DEBUG_LOCAL) write (OUTPUT_UNIT, *) 'Leaving Distribute_var2d...'

    end subroutine Distribute_var2d

    subroutine Do_halo_exchange (patch, ims, ime, jms, jme, ips, ipe, jps, jpe, nghost, cart_comm)

    ! It does not update the corners in the halo

#ifdef DM_PARALLEL
      use mpi
#endif

      implicit none

      integer, intent (in) :: cart_comm, ims, ime, jms, jme, ips, ipe, jps, jpe, nghost
      real, dimension(ims:ime, jms:jme), intent (in out) :: patch

      integer :: ierr, nbr_left, nbr_right, nbr_up, nbr_down, tag_base, nx, ny
      integer, dimension(8) :: reqs
      real, dimension(:), allocatable :: sendbuf_right, recvbuf_left, sendbuf_left, recvbuf_right, &
                                         sendbuf_up, recvbuf_down, sendbuf_down, recvbuf_up

      integer :: rank, i, j, k
      integer, dimension(2) :: coords

      ! New variables to handle the topology mapping
      integer :: dim_x, dim_y

#ifdef DM_PARALLEL

      call MPI_Comm_rank(cart_comm, rank, ierr)
      call MPI_Cart_coords(cart_comm, rank, 2, coords, ierr)

      nx = ipe - ips + 1
      ny = jpe - jps + 1

      tag_base = 1000

      ! ---------------------------------------------------------
      ! ADAPT DIMENSIONS BASED ON MODULE FLAG
      ! ---------------------------------------------------------
      if (topology_dim_order == 0) then
          ! Default internal behavior [X, Y]
          dim_x = 0
          dim_y = 1
      else
          ! External framework behavior [Y, X]
          dim_x = 1
          dim_y = 0
      end if

      ! Get neighbor ranks in Cartesian topology using the mapped dimensions
      call MPI_Cart_shift(cart_comm, dim_x, 1, nbr_left, nbr_right, ierr)
      call MPI_Cart_shift(cart_comm, dim_y, 1, nbr_down, nbr_up, ierr)
      ! ---------------------------------------------------------

        ! Allocate buffers
      allocate (sendbuf_right(ny * nghost), recvbuf_left(ny * nghost))
      allocate (sendbuf_left(ny * nghost), recvbuf_right(ny * nghost))
      allocate (sendbuf_up(nx * nghost), recvbuf_down(nx * nghost))
      allocate (sendbuf_down(nx * nghost), recvbuf_up(nx * nghost))

        ! Send RIGHT, receive LEFT
      k = 0
      do j = jps, jpe
        do i = ipe - nghost + 1, ipe
          k = k + 1
          sendbuf_right(k) = patch(i, j)
        end do
      end do
      call MPI_Irecv(recvbuf_left, ny*nghost, MPI_REAL, nbr_left, tag_base + 0, cart_comm, reqs(1), ierr)
      call MPI_Isend(sendbuf_right, ny*nghost, MPI_REAL, nbr_right, tag_base + 0, cart_comm, reqs(2), ierr)

      ! Send LEFT, receive RIGHT
      k = 0
      do j = jps, jpe
        do i = ips, ips + nghost - 1
          k = k + 1
          sendbuf_left(k) = patch(i, j)
        end do
      end do
      call MPI_Irecv (recvbuf_right, ny * nghost, MPI_REAL, nbr_right, tag_base + 1, cart_comm, reqs(3), ierr)
      call MPI_Isend (sendbuf_left, ny * nghost, MPI_REAL, nbr_left, tag_base + 1, cart_comm, reqs(4), ierr)

        ! Send UP, receive DOWN
      k = 0
      do j = jpe - nghost + 1, jpe
        do i = ips, ipe
          k = k + 1
          sendbuf_up(k) = patch(i, j)
        end do
      end do
      call MPI_Irecv (recvbuf_down, nx * nghost, MPI_REAL, nbr_down, tag_base + 2, cart_comm, reqs(5), ierr)
      call MPI_Isend (sendbuf_up, nx * nghost, MPI_REAL, nbr_up, tag_base + 2, cart_comm, reqs(6), ierr)

      ! Send DOWN, receive UP
      k = 0
      do j = jps, jps + nghost - 1
        do i = ips, ipe
          k = k + 1
          sendbuf_down(k) = patch(i, j)
        end do
      end do
      call MPI_Irecv(recvbuf_up, nx * nghost, MPI_REAL, nbr_up, tag_base + 3, cart_comm, reqs(7), ierr)
      call MPI_Isend(sendbuf_down, nx * nghost, MPI_REAL, nbr_down, tag_base + 3, cart_comm, reqs(8), ierr)

      ! Wait for all communications
      call MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE, ierr)

        ! Unpack ghost zones
        ! LEFT ghost
      if (nbr_left /= MPI_PROC_NULL) then
        k = 0
        do j = jps, jpe
          do i = ips - nghost, ips - 1
            k = k + 1
            patch(i, j) = recvbuf_left(k)
          end do
        end do
      end if

        ! RIGHT ghost
      if (nbr_right /= MPI_PROC_NULL) then
        k = 0
        do j = jps, jpe
          do i = ipe + 1, ipe + nghost
            k = k + 1
            patch(i, j) = recvbuf_right(k)
          end do
        end do
      end if

        ! UP ghost (top halo rows)
      if (nbr_up /= MPI_PROC_NULL) then
        k = 0
        do j = jpe + 1, jpe + nghost
          do i = ips, ipe
            k = k + 1
            patch(i, j) = recvbuf_up(k)
          end do
        end do
      end if

        ! DOWN ghost (bottom halo rows)
      if (nbr_down /= MPI_PROC_NULL) then
        k = 0
        do j = jps - nghost, jps - 1
          do i = ips, ipe
            k = k + 1
            patch(i, j) = recvbuf_down(k)
          end do
        end do
      end if

        ! Deallocate buffers
      deallocate (sendbuf_right, recvbuf_left)
      deallocate (sendbuf_left, recvbuf_right)
      deallocate (sendbuf_up, recvbuf_down)
      deallocate (sendbuf_down, recvbuf_up)

#endif

    end subroutine Do_halo_exchange

    subroutine Do_halo_exchange_with_corners (patch, ims, ime, jms, jme, ips, ipe, jps, jpe, nghost, cart_comm)

#ifdef DM_PARALLEL
      use mpi
#endif

      implicit none

      integer, intent (in) :: cart_comm, ims, ime, jms, jme, ips, ipe, jps, jpe, nghost
      real, dimension(ims:ime, jms:jme), intent (inout) :: patch

      integer :: ierr, nbr_left, nbr_right, nbr_up, nbr_down, tag_base
      integer :: nx, ny, rank, i, j, k
      integer, dimension(8) :: reqs
      integer, dimension(2) :: coords
      real, dimension(:), allocatable :: sendbuf_right, recvbuf_left, sendbuf_left, recvbuf_right
      real, dimension(:), allocatable :: sendbuf_up, recvbuf_down, sendbuf_down, recvbuf_up

      ! New variables to handle the topology mapping
      integer :: dim_x, dim_y

#ifdef DM_PARALLEL
      call MPI_Comm_rank(cart_comm, rank, ierr)
      call MPI_Cart_coords(cart_comm, rank, 2, coords, ierr)

      nx = ipe - ips + 1
      ny = jpe - jps + 1
      tag_base = 1000

      ! ---------------------------------------------------------
      ! ADAPT DIMENSIONS BASED ON MODULE FLAG
      ! ---------------------------------------------------------
      if (topology_dim_order == 0) then
          ! Default internal behavior [X, Y]
          dim_x = 0
          dim_y = 1
      else
          ! External framework behavior [Y, X]
          dim_x = 1
          dim_y = 0
      end if

      ! Neighbor ranks using mapped dimensions
      call MPI_Cart_shift (cart_comm, dim_x, 1, nbr_left, nbr_right, ierr)
      call MPI_Cart_shift (cart_comm, dim_y, 1, nbr_down, nbr_up, ierr)
      ! ---------------------------------------------------------

        ! Halo exchange in X direction
      allocate (sendbuf_right(ny * nghost), recvbuf_left(ny * nghost))
      allocate (sendbuf_left(ny * nghost), recvbuf_right(ny * nghost))

        ! Pack right send buffer
      k = 0
      do j = jps, jpe
        do i = ipe - nghost + 1, ipe
          k = k + 1
          sendbuf_right(k) = patch(i, j)
        end do
      end do

        ! Pack left send buffer
      k = 0
      do j = jps, jpe
        do i = ips, ips + nghost - 1
          k = k + 1
          sendbuf_left(k) = patch(i, j)
        end do
      end do

        ! Exchange left/right
      call MPI_Irecv (recvbuf_left, ny * nghost, MPI_REAL, nbr_left, tag_base + 0, cart_comm, reqs(1), ierr)
      call MPI_Isend (sendbuf_right, ny * nghost, MPI_REAL, nbr_right, tag_base + 0, cart_comm, reqs(2), ierr)
      call MPI_Irecv (recvbuf_right, ny * nghost, MPI_REAL, nbr_right, tag_base + 1, cart_comm, reqs(3), ierr)
      call MPI_Isend (sendbuf_left, ny * nghost, MPI_REAL, nbr_left, tag_base + 1, cart_comm, reqs(4), ierr)

      call MPI_Waitall (4, reqs, MPI_STATUSES_IGNORE, ierr)

        ! Unpack left halo
      if (nbr_left /= MPI_PROC_NULL) then
        k = 0
        do j = jps, jpe
          do i = ips - nghost, ips - 1
            k = k + 1
            patch(i, j) = recvbuf_left(k)
          end do
        end do
      end if

        ! Unpack right halo
      if (nbr_right /= MPI_PROC_NULL) then
        k = 0
        do j = jps, jpe
          do i = ipe + 1, ipe + nghost
            k = k + 1
            patch(i, j) = recvbuf_right(k)
          end do
        end do
      end if

      deallocate (sendbuf_right, recvbuf_left, sendbuf_left, recvbuf_right)

        ! Exchange in Y direction (including halos in X direction)
      allocate (sendbuf_up((nx + 2 * nghost) * nghost), recvbuf_down((nx + 2 * nghost) * nghost))
      allocate (sendbuf_down((nx + 2 * nghost) * nghost), recvbuf_up((nx + 2 * nghost) * nghost))

        ! Pack up send buffer
      k = 0
      do j = jpe - nghost + 1, jpe
        do i = ips - nghost, ipe + nghost
          k = k + 1
          sendbuf_up(k) = patch(i, j)
        end do
      end do

        ! Pack down send buffer
      k = 0
      do j = jps, jps + nghost - 1
        do i = ips - nghost, ipe + nghost
          k = k + 1
          sendbuf_down(k) = patch(i, j)
        end do
      end do

        ! Exchange up/down
      call MPI_Irecv (recvbuf_down, (nx + 2 * nghost) * nghost, MPI_REAL, nbr_down, tag_base + 2, cart_comm, reqs(1), ierr)
      call MPI_Isend (sendbuf_up, (nx + 2 * nghost) * nghost, MPI_REAL, nbr_up, tag_base + 2, cart_comm, reqs(2), ierr)
      call MPI_Irecv (recvbuf_up, (nx + 2 * nghost) * nghost, MPI_REAL, nbr_up, tag_base + 3, cart_comm, reqs(3), ierr)
      call MPI_Isend (sendbuf_down, (nx + 2 * nghost) * nghost, MPI_REAL, nbr_down, tag_base + 3, cart_comm, reqs(4), ierr)

      call MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE, ierr)

        ! Unpack down halo
      if (nbr_down /= MPI_PROC_NULL) then
        k = 0
        do j = jps - nghost, jps - 1
          do i = ips - nghost, ipe + nghost
            k = k + 1
            patch(i, j) = recvbuf_down(k)
          end do
        end do
      end if

        ! Unpack up halo
      if (nbr_up /= MPI_PROC_NULL) then
        k = 0
        do j = jpe + 1, jpe + nghost
          do i = ips - nghost, ipe + nghost
            k = k + 1
            patch(i, j) = recvbuf_up(k)
          end do
        end do
      end if

      deallocate (sendbuf_up, recvbuf_down, sendbuf_down, recvbuf_up)
#endif

    end subroutine Do_halo_exchange_with_corners

    subroutine Gather_var2d (cfbm_comm, nx, ny, ifps, ifpe, jfps, jfpe, var2d_local, var2d_global)

#ifdef DM_PARALLEL
      use mpi_f08
#endif

      implicit none

      integer, intent (in) :: cfbm_comm, nx, ny, ifps, ifpe, jfps, jfpe
      real, dimension(ifps:ifpe, jfps:jfpe), intent (in) :: var2d_local
      real, dimension(nx, ny), intent (out) :: var2d_global

      real, dimension(:, :), allocatable :: buf
      integer :: rank, ierr, ntasks, src, i_start, j_start, nx_local, ny_local, nxl, nyl

      integer, dimension(4) :: metadata
#ifdef DM_PARALLEL
      type(MPI_Comm) :: cfbm_comm_f08
#endif


#ifdef DM_PARALLEL
      call Convert_mpi_comm_to_f08 (cfbm_comm, cfbm_comm_f08)

      call Mpi_comm_size (cfbm_comm_f08, ntasks, ierr)
      call Mpi_comm_rank (cfbm_comm_f08, rank, ierr)

      nx_local = ifpe - ifps + 1
      ny_local = jfpe - jfps + 1
      if (rank /= 0) then
        metadata = [ifps, jfps, nx_local, ny_local]
        call Mpi_Send (metadata, 4, MPI_INTEGER, 0, 200, cfbm_comm_f08, ierr)
        call Mpi_Send (var2d_local, nx_local * ny_local, MPI_REAL, 0, 100, cfbm_comm_f08, ierr)
      else
          ! Rank 0 places its own data
        var2d_global(ifps:ifpe, jfps:jfpe) = var2d_local

          ! Receive from others
        do src = 1, ntasks - 1
            ! Receive metadata
          call Mpi_Recv (metadata, 4, MPI_INTEGER, src, 200, cfbm_comm_f08, MPI_STATUS_IGNORE, ierr)
          i_start = metadata(1)
          j_start = metadata(2)
          nxl = metadata(3)
          nyl = metadata(4)
          allocate (buf(nxl, nyl))
          call Mpi_Recv (buf, nxl * nyl, MPI_REAL, src, 100, cfbm_comm_f08, MPI_STATUS_IGNORE, ierr)

            ! Place in global array
          var2d_global(i_start:i_start + nxl - 1, j_start:j_start + nyl - 1) = buf

          deallocate (buf)
        end do
      end if
#endif

    end subroutine Gather_var2d

    subroutine Max_across_mpi_tasks (local_max, cart_comm, global_max)

#ifdef DM_PARALLEL
      use mpi
#endif

      implicit none

      real, intent (in) :: local_max
      integer, intent (in) :: cart_comm
      real, intent (out) :: global_max

      integer :: ierr


#ifdef DM_PARALLEL
      call MPI_Allreduce(local_max, global_max, 1, MPI_REAL, MPI_MAX, cart_comm, ierr)
#else
      global_max = local_max
#endif

    end subroutine Max_across_mpi_tasks

    subroutine Min_across_mpi_tasks (local_min, cart_comm, global_min)

#ifdef DM_PARALLEL
      use mpi
#endif

      implicit none

      real, intent (in) :: local_min
      integer, intent (in) :: cart_comm
      real, intent (out) :: global_min

      integer :: ierr


#ifdef DM_PARALLEL
      call MPI_Allreduce (local_min, global_min, 1, MPI_REAL, MPI_MIN, cart_comm, ierr)
#else
      global_min = local_min
#endif

    end subroutine Min_across_mpi_tasks

    subroutine Print_cart_info (cart_comm)

#ifdef DM_PARALLEL
      use mpi
#endif

      implicit none
      integer, intent(in) :: cart_comm

      integer :: rank, ierr, ndims, topo_type, d
      integer, allocatable, dimension(:) :: dims, coords
      logical, allocatable, dimension(:) :: periods

#ifdef DM_PARALLEL
      ! Get the rank of the calling process
      call MPI_Comm_rank(cart_comm, rank, ierr)

      ! 1. Verify this is actually a Cartesian communicator
      call MPI_Topo_test(cart_comm, topo_type, ierr)

      if (topo_type /= MPI_CART) then
         if (rank == 0) print *, "DEBUG ERROR: The provided communicator is NOT a Cartesian topology."
         return
      end if

      ! 2. Get the number of dimensions
      call MPI_Cartdim_get(cart_comm, ndims, ierr)

      if (ierr == MPI_SUCCESS) then
         allocate(dims(ndims))
         allocate(periods(ndims))
         allocate(coords(ndims))

         ! 3. Extract the full topology information
         call MPI_Cart_get(cart_comm, ndims, dims, periods, coords, ierr)

         ! 4. Print the overall grid configuration (Restricted to Rank 0 to avoid console spam)
         if (rank == 0) then
            print *, "=========================================="
            print *, "       CARTESIAN TOPOLOGY DEBUG           "
            print *, "=========================================="
            print *, "Number of Dimensions: ", ndims
            do d = 1, ndims
               print *, "--- Dimension ", d, " ---"
               print *, "  Grid Size (Ranks):  ", dims(d)
               print *, "  Is Periodic?:       ", periods(d)
            end do
            print *, "------------------------------------------"
         end if

         ! Wait for Rank 0 to finish printing the header
         call MPI_Barrier(cart_comm, ierr)

         ! 5. Print individual rank coordinates
         ! (Restricted to Ranks 0 and 1 for your specific debugging needs)
         if (rank == 0 .or. rank == 1) then
             print *, "DEBUG -> Rank ", rank, " is at Coordinates: ", coords
         end if

         ! Optional: Barrier again to keep standard output clean before the program continues
         call MPI_Barrier(cart_comm, ierr)

         deallocate(dims, periods, coords)
      else
         print *, "DEBUG ERROR [Rank ", rank, "]: Failed to get Cartesian dimensions."
      end if
#endif

    end subroutine Print_cart_info

    subroutine Sum_across_mpi_tasks (local_sum, cart_comm, global_sum)

#ifdef DM_PARALLEL
      use mpi
#endif

      implicit none

      real, intent (in) :: local_sum
      integer, intent (in) :: cart_comm
      real, intent (out) :: global_sum

      integer :: ierr


#ifdef DM_PARALLEL
      call MPI_Allreduce(local_sum, global_sum, 1, MPI_REAL, MPI_SUM, cart_comm, ierr)
#else
      global_sum = local_sum
#endif

    end subroutine Sum_across_mpi_tasks

  end module mpi_mod
