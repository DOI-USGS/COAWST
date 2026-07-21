  module tiles_mod

    use stderrout_mod, only : Stop_simulation

    implicit none

    private

    integer, parameter :: TILE_NONE = 0, TILE_X = 1, TILE_Y = 2, TILE_XY = 3

    public :: Calc_tiles_dims

  contains

    subroutine Calc_tiles_dims (spx, epx, spy, epy, num_tiles, tile_strategy, i_start, i_end, j_start, j_end)

      ! This code is borrowed from module_tiles.F in WRF, specifically SUBROUTINE set_tiles2

      implicit none

      integer, intent (in) :: spx, epx, spy, epy, tile_strategy
      integer, intent (in out) :: num_tiles

      integer, parameter :: MIN_TILES_IN_X = 1, MIN_TILES_IN_Y = 1
      integer, dimension(:), allocatable :: i_start, i_end, j_start, j_end
      integer :: num_tiles_x, num_tiles_y
      integer :: t, ts, te
      integer :: ntiles
      integer :: nt


        ! Consistency check
      if (num_tiles <= 0) then
        num_tiles = 0
        return
      end if

      allocate (i_start(num_tiles))
      allocate (i_end(num_tiles))
      allocate (j_start(num_tiles))
      allocate (j_end(num_tiles))

        ! Calc number of tiles in x and y based on total number of tiles
      select case (tile_strategy)
        case (TILE_NONE, TILE_Y)
          if (num_tiles > (epy - spy + 1)) call Stop_simulation ('Number of tiles is too high for TILE_Y strategy')
          num_tiles_x = 1
          num_tiles_y = num_tiles

        case (TILE_X)
          num_tiles_x = num_tiles
          num_tiles_y = 1

        case (TILE_XY)
          call least_aspect (num_tiles, MIN_TILES_IN_Y, MIN_TILES_IN_X, num_tiles_y, num_tiles_x)

        case default
          call Stop_simulation ('The tile strategy selected is not valid.')

      end select

        ! Calc start and end tile indices
      nt = 1
      do t = 0, num_tiles - 1
        ntiles = t / num_tiles_x
        CALL region_bounds (spy, epy, num_tiles_y, ntiles, ts, te)
        if (ts <= te) then  ! converse happens if number of tiles > number of points in dim
          j_start(nt) = ts
          j_end(nt) = te
          ntiles = mod (t, num_tiles_x)
          call region_bounds (spx, epx, num_tiles_x, ntiles, ts, te)
          if (ts <= te) then  ! converse happens if number of tiles > number of points in dim
            i_start(nt) = ts
            i_end(nt) = te
            nt = nt + 1
          end if
        end if
      end do
      nt = nt - 1

      if (nt /= num_tiles) then
        num_tiles = nt
        if (num_tiles == 0) then
          deallocate (i_start)
          deallocate (j_start)
          deallocate (i_end)
          deallocate (j_end)
        else
          i_start = reshape (i_start, [num_tiles])
          i_end = reshape (i_end, [num_tiles])
          j_start = reshape (j_start, [num_tiles])
          j_end = reshape (j_end, [num_tiles])
        end if
      end if

    end subroutine Calc_tiles_dims

    subroutine least_aspect (nparts, minparts_y, minparts_x, nparts_y, nparts_x)

      implicit none

      integer, intent (in) :: nparts, minparts_y, minparts_x
      integer, intent (out) :: nparts_y, nparts_x

      integer :: x, y, mini


      mini = 2 * nparts
      nparts_y = 1
      nparts_x = nparts

      do y = 1, nparts
        if (mod (nparts, y) == 0) then
          x = nparts / y
          if (abs (y - x) < mini .and. y >= minparts_y .and. x >= minparts_x) then
            mini = abs (y - x)
            nparts_y = y
            nparts_x = x
          end if
        end if
      end do

    end subroutine least_aspect

    pure subroutine region_bounds (region_start, region_end, num_p, p, patch_start, patch_end)

    ! 1-D decomposition routine: Given starting and ending indices of a
    ! vector, the number of patches dividing the vector, and the number of
    ! the patch, give the start and ending indices of the patch within the
    ! vector.  This will work with tiles too.  
    !
    ! Implementation note.  This is
    ! implemented somewhat inefficiently, now, with a loop, so we can use the
    ! locproc function, which returns processor number for a given
    ! index, whereas what we want is index for a given processor number.
    ! With a little thought and a lot of debugging, we can come up with a
    ! direct expression for what we want.  For time being, we loop...
    ! Remember that processor numbering starts with zero.

     implicit none

     integer, intent (in) :: region_start, region_end, num_p, p
     integer, intent (out) :: patch_start, patch_end
     integer :: offset, i


     patch_end = -999999999
     patch_start = 999999999
     offset = region_start
     do i = 0, region_end - offset
       if (locproc (i, region_end - region_start + 1, num_p) == p) then
         patch_end = max (patch_end, i)
         patch_start = min (patch_start, i)
       end if
     end do

     patch_start = patch_start + offset
     patch_end = patch_end + offset

     return
   end subroutine region_bounds

   recursive pure subroutine rlocproc (p, maxi, nproc, ml, mr, ret)

     implicit none

     integer, intent(in)  :: p, maxi, nproc, ml, mr
     integer, intent(out) :: ret

     integer :: width, rem, ret2, bl, br, mid, adjust, p_r, maxi_r, nproc_r, zero


      adjust = 0
      rem = mod (maxi, nproc)
      width = maxi / nproc
      mid = maxi / 2
      if (rem > 0 .and. (((mod (rem, 2) == 0) .or. (rem > 2)) .or. (p <= mid))) then
        width = width + 1
      end if
      if (p <= mid .and. mod (rem, 2) .ne. 0) then
        adjust = adjust + 1
      end if
      bl = max (width, ml)
      br = max (width, mr)
      if (p < bl) then
        ret = 0
      else if (p > maxi - br - 1) then
        ret = nproc - 1
      else
        p_r = p - bl
        maxi_r = maxi - bl - br + adjust
        nproc_r = max (nproc - 2, 1)
        zero = 0
          ! Recursive
        CALL rlocproc (p_r, maxi_r, nproc_r, zero, zero, ret2)
        ret = ret2 + 1
      end if

      return
    end subroutine rlocproc

    integer pure function locproc (i, m, numpart)

      implicit none

      integer, intent(in) :: i, m, numpart

      integer :: retval, ii, im, inumpart, zero


      ii = i
      im = m
      inumpart = numpart
      zero = 0
      CALL rlocproc (ii, im, inumpart, zero, zero, retval)
      locproc = retval

      return
    end function locproc

  end module tiles_mod
