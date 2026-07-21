module NoahmpWriteLandMod

   use mpi
   use netcdf
   use NoahmpIOVarType

   implicit none

   integer, save, private :: ncid, terrain, snowh, shbxy, evbxy, &
                             vegfra, gvfmin, gvfmax, tsk, emiss, &
                             albsfcdirxy, albsfcdifxy, &
                             savxy, sagxy, pahxy, firaxy, hfx, &
                             lh, grdflx, ghbxy, canhsxy, tslb, smois, &
                             tau_ew, tau_ns

contains

   subroutine NoahmpWriteLand(NoahmpIO, filenum, maxblocks)

      implicit none

      type(NoahmpIO_type), intent(inout) :: NoahmpIO
      integer, intent(in) :: filenum, maxblocks

      ! local variables
      integer :: ierr, start(2), count(2), nx, ny, comp2d, nsoil
      character(len=5) :: ts_str
      character(len=1) :: lev_str
      character(len=100) :: dir, filename
      logical :: ex

      if (NoahmpIO%blkid == 0) then
         write (ts_str, '(I5.5)') filenum
         write (lev_str, '(I1.1)') NoahmpIO%LEVEL

         dir = "lnd"//trim(ts_str)
         inquire (file=trim(dir), exist=ex)
         if (.not. ex) then
            call execute_command_line("mkdir -p "//trim(dir), exitstat=ierr)
            if (ierr /= 0) then
               print *, "Failed to create directory: ", trim(dir)
               stop
            end if
         end if

         filename = trim(dir)//"/Level_"//trim(lev_str)//".nc"
         ierr = nf90_create(trim(filename), IOR(NF90_CLOBBER, IOR(NF90_NETCDF4, NF90_MPIIO)), &
                            ncid, comm=NoahmpIO%comm, info=MPI_INFO_NULL)

         if (ierr /= nf90_noerr) then
            print *, "NetCDF create failed: ", trim(nf90_strerror(ierr))
            stop
         end if

         ! Define dimensions and variable
         ierr = nf90_def_dim(ncid, "NX", NoahmpIO%xsglobal, nx)
         ierr = nf90_def_dim(ncid, "NY", NoahmpIO%ysglobal, ny)
         ierr = nf90_def_dim(ncid, "COMP2D", 2, comp2d)
         ierr = nf90_def_dim(ncid, "NSOIL", NoahmpIO%NSOIL, nsoil)
         ierr = nf90_def_var(ncid, "TERRAIN", NF90_FLOAT, (/nx, ny/), terrain)
         ierr = nf90_def_var(ncid, "SNOWH", NF90_FLOAT, (/nx, ny/), snowh)
         ierr = nf90_def_var(ncid, "SHBXY", NF90_FLOAT, (/nx, ny/), shbxy)
         ierr = nf90_def_var(ncid, "EVBXY", NF90_FLOAT, (/nx, ny/), evbxy)
         ierr = nf90_def_var(ncid, "VEGFRA", NF90_FLOAT, (/nx, ny/), vegfra)
         ierr = nf90_def_var(ncid, "GVFMIN", NF90_FLOAT, (/nx, ny/), gvfmin)
         ierr = nf90_def_var(ncid, "GVFMAX", NF90_FLOAT, (/nx, ny/), gvfmax)
         ierr = nf90_def_var(ncid, "TSK", NF90_FLOAT, (/nx, ny/), tsk)
         ierr = nf90_def_var(ncid, "EMISS", NF90_FLOAT, (/nx, ny/), emiss)
         ierr = nf90_def_var(ncid, "ALBSFCDIRXY", NF90_FLOAT, (/nx, comp2d, ny/), albsfcdirxy)
         ierr = nf90_def_var(ncid, "ALBSFCDIFXY", NF90_FLOAT, (/nx, comp2d, ny/), albsfcdifxy)
         ierr = nf90_def_var(ncid, "SAVXY", NF90_FLOAT, (/nx, ny/), savxy)
         ierr = nf90_def_var(ncid, "SAGXY", NF90_FLOAT, (/nx, ny/), sagxy)
         ierr = nf90_def_var(ncid, "PAHXY", NF90_FLOAT, (/nx, ny/), pahxy)
         ierr = nf90_def_var(ncid, "FIRAXY", NF90_FLOAT, (/nx, ny/), firaxy)
         ierr = nf90_def_var(ncid, "HFX", NF90_FLOAT, (/nx, ny/), hfx)
         ierr = nf90_def_var(ncid, "LH", NF90_FLOAT, (/nx, ny/), lh)
         ierr = nf90_def_var(ncid, "GRDFLX", NF90_FLOAT, (/nx, ny/), grdflx)
         ierr = nf90_def_var(ncid, "GHBXY", NF90_FLOAT, (/nx, ny/), ghbxy)
         ierr = nf90_def_var(ncid, "CANHSXY", NF90_FLOAT, (/nx, ny/), canhsxy)
         ierr = nf90_def_var(ncid, "TSLB", NF90_FLOAT, (/nx, nsoil, ny/), tslb)
         ierr = nf90_def_var(ncid, "SMOIS", NF90_FLOAT, (/nx, nsoil, ny/), smois)
         ierr = nf90_def_var(ncid, "TAU_EW", NF90_FLOAT, (/nx, ny/), tau_ew)
         ierr = nf90_def_var(ncid, "TAU_NS", NF90_FLOAT, (/nx, ny/), tau_ns)
         ! End definition mode
         ierr = nf90_enddef(ncid)
      end if

      ! Select portion to write
      start = (/NoahmpIO%xstart-NoahmpIO%xoffset+1, NoahmpIO%ystart-NoahmpIO%yoffset+1/)
      count = (/NoahmpIO%xend-NoahmpIO%xstart+1, NoahmpIO%yend-NoahmpIO%ystart+1/)

      ! Write data
      ierr = nf90_put_var(ncid, terrain, NoahmpIO%TERRAIN, start=start, count=count)
      ierr = nf90_put_var(ncid, snowh, NoahmpIO%SNOWH, start=start, count=count)
      ierr = nf90_put_var(ncid, shbxy, NoahmpIO%SHBXY, start=start, count=count)
      ierr = nf90_put_var(ncid, evbxy, NoahmpIO%EVBXY, start=start, count=count)
      ierr = nf90_put_var(ncid, vegfra, NoahmpIO%VEGFRA, start=start, count=count)
      ierr = nf90_put_var(ncid, gvfmin, NoahmpIO%GVFMIN, start=start, count=count)
      ierr = nf90_put_var(ncid, gvfmax, NoahmpIO%GVFMAX, start=start, count=count)
      ierr = nf90_put_var(ncid, tsk, NoahmpIO%TSK, start=start, count=count)
      ierr = nf90_put_var(ncid, emiss, NoahmpIO%EMISS, start=start, count=count)
      ierr = nf90_put_var(ncid, albsfcdirxy, NoahmpIO%ALBSFCDIRXY, start=(/start(1), 1, start(2)/), &
                                                                   count=(/count(1), 2, count(2)/))
      ierr = nf90_put_var(ncid, albsfcdifxy, NoahmpIO%ALBSFCDIFXY, start=(/start(1), 1, start(2)/), &
                                                                   count=(/count(1), 2, count(2)/))
      ierr = nf90_put_var(ncid, savxy, NoahmpIO%SAVXY, start=start, count=count)
      ierr = nf90_put_var(ncid, sagxy, NoahmpIO%SAGXY, start=start, count=count)
      ierr = nf90_put_var(ncid, pahxy, NoahmpIO%PAHXY, start=start, count=count)
      ierr = nf90_put_var(ncid, firaxy, NoahmpIO%FIRAXY, start=start, count=count)
      ierr = nf90_put_var(ncid, hfx, NoahmpIO%HFX, start=start, count=count)
      ierr = nf90_put_var(ncid, lh, NoahmpIO%LH, start=start, count=count)
      ierr = nf90_put_var(ncid, grdflx, NoahmpIO%GRDFLX, start=start, count=count)
      ierr = nf90_put_var(ncid, ghbxy, NoahmpIO%GHBXY, start=start, count=count)
      ierr = nf90_put_var(ncid, canhsxy, NoahmpIO%CANHSXY, start=start, count=count)
      ierr = nf90_put_var(ncid, tslb, NoahmpIO%TSLB, start=(/start(1), 1, start(2)/), &
                                                     count=(/count(1), NoahmpIO%NSOIL, count(2)/))
      ierr = nf90_put_var(ncid, smois, NoahmpIO%SMOIS, start=(/start(1), 1, start(2)/), &
                                                       count=(/count(1), NoahmpIO%NSOIL, count(2)/))
      ierr = nf90_put_var(ncid, tau_ew, NoahmpIO%TAU_EW, start=start, count=count)
      ierr = nf90_put_var(ncid, tau_ns, NoahmpIO%TAU_NS, start=start, count=count)

      if (NoahmpIO%blkid == (maxblocks-1)) then
         ! Close file
         ierr = nf90_close(ncid)
      end if

   end subroutine NoahmpWriteLand

end module NoahmpWriteLandMod
