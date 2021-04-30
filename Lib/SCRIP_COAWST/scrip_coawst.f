!======================================================================
!     This is a wrapper developed at USGS, Woods Hole Coastal and Marine
!     Science Center. This routine calls the SCRIP interpolation 
!     package routines that are the driver for computing
!     interpolation addresses and weights between two grids.
!
!---- Written by John C. Warner-------------------------------------- 
!-----         Tarandeep S. Kalra -------------------------------------
!--------------Date: 08/04/2015----------------------------------------
!======================================================================

      program scrip_coawst

      use kinds_mod
      use constants
      use iounits 
      use scripwrap_mod 
      use read_swan 
      use read_roms 
      use read_wrf
      use read_hydro
      use read_ww3
      use create_masks

      implicit none
#ifdef MPI
      include 'mpif.h'
      integer (kind=int_kind) :: Nnodes
#endif
      integer (kind=int_kind) :: MyComm
!     local variables
      character(len=char_len) :: inputfile
      integer (kind=int_kind) :: OutThread
!
#ifdef MPI
!  Initialize MPI execution environment.
!
      CALL mpi_init (MyError)
!
!  Get rank of the local process in the group associated with the
!  comminicator.
!
      CALL mpi_comm_size (MPI_COMM_WORLD, Nnodes, MyError)
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)
      MyComm=MPI_COMM_WORLD
      OutThread=MyRank
#else
      MyComm=0
      OutThread=0
#endif
!
!     Reading input file 
      call getarg(1,inputfile)
!     inputfile = arg1
!     write(*,*) 'infiel si ', inputfile, arg1
!     write(stdout,*) "---------------------------------------------"
!     write(stdout,*) "Enter the input file for SCRIP_COAWST package"
!     write(stdout,*) "---------------------------------------------"
!     read(stdin,*) inputfile
      inputfile = adjustl(inputfile)
      open(unit=iunit,file=inputfile,status='old',form='formatted')
        call read_inputs( MyComm )
      close(iunit)

!     Read the information from different grids 
      if ((Ngrids_swan>0).and.(Ngrids_ww3>0)) then 
        write(stdout,*) 'Only use WW3 or SWAN, not both '
        write(stdout,*) 'Num WW3 grids = ', Ngrids_ww3
        write(stdout,*) 'Num SWAN grids = ', Ngrids_swan
      end if

      if (Ngrids_roms>0) then
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Reading ROMS grids'
        END IF
        call load_roms_grid( MyComm )
      end if
      if (Ngrids_swan>0) then 
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Reading SWAN grids'
        END IF
        call load_swan_grid( MyComm )
      end if
      if (Ngrids_ww3>0) then 
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Reading WW3 grids'
        END IF
        call load_ww3_grid()
      end if
      if (Ngrids_wrf>0) then
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Reading WRF grids'
        END IF
        call load_wrf_grid( MyComm )
      end if
      if (Ngrids_hyd>0) then
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Reading HYD grids'
        END IF
        call load_hydro_grid( MyComm )
      end if

!     Calculate interpolation weights between combinations of grids
      if ((Ngrids_roms>0).and.(Ngrids_swan>0)) then 
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Calling ocn2wav_mask'
        END IF
        call ocn2wav_mask(MyComm)
      end if
      if ((Ngrids_roms>0).and.(Ngrids_ww3>0)) then 
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Calling ocn2ww3_mask'
        END IF
        call ocn2ww3_mask(MyComm)
      end if
      if ((Ngrids_roms>0).and.(Ngrids_wrf>0))  then 
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Calling ocn2atm_mask'
        END IF
        call ocn2atm_mask(MyComm)
      end if

      if ((Ngrids_swan>0).and.(Ngrids_roms>0)) then 
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Calling wav2ocn_mask'
        END IF
        call wav2ocn_mask(MyComm)
      end if 
      if ((Ngrids_swan>0).and.(Ngrids_wrf>0)) then
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Calling wav2atm_mask'
        END IF
        call wav2atm_mask(MyComm)
      end if
      if ((Ngrids_ww3>0).and.(Ngrids_roms>0)) then 
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Calling ww32ocn_mask'
        END IF
        call ww32ocn_mask(MyComm)
      end if 
      if ((Ngrids_ww3>0).and.(Ngrids_wrf>0)) then
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Calling ww32atm_mask'
        END IF
        call ww32atm_mask(MyComm)
      end if

      if ((Ngrids_wrf>0).and.(Ngrids_roms>0)) then
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Calling atm2ocn_mask'
        END IF
        call atm2ocn_mask(MyComm)
      end if
      if ((Ngrids_wrf>0).and.(Ngrids_swan>0)) then
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Calling atm2wav_mask'
        END IF
        call atm2wav_mask(MyComm)
      end if
      if ((Ngrids_wrf>0).and.(Ngrids_ww3>0)) then
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Calling atm2ww3_mask'
        END IF
        call atm2ww3_mask(MyComm)
      end if

      if ((Ngrids_hyd>0).and.(Ngrids_roms>0)) then
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Calling hyd2ocn_mask'
        END IF
        call hyd2ocn_mask(MyComm)
      end if
      if ((Ngrids_roms>0).and.(Ngrids_hyd>0)) then
        IF (OutThread.eq.0) THEN
          write(stdout,*) 'Calling ocn2hyd_mask'
        END IF
        call ocn2hyd_mask(MyComm)
      end if

#ifdef MPI
      CALL mpi_finalize (MyError)
#endif
      STOP
      end program scrip_coawst

!======================================================================

      subroutine read_inputs( MyComm )
      use kinds_mod
      use iounits
      use scripwrap_mod
 
      implicit none
 
      integer (kind=int_kind), intent(in) :: MyComm

#ifdef MPI
      include 'mpif.h'
      integer (kind=int_kind) :: MyRank, Nnodes, MyError
#endif

!     local variables
      integer(int_kind) :: i

!     Allocate input variables
      namelist /inputs/ Ngrids_roms, Ngrids_swan, Ngrids_ww3,           &
     &                  Ngrids_wrf,  Ngrids_hyd,                        &
     &                  roms_grids, swan_coord, swan_bath,              &
     &                  swan_numx, swan_numy, cartesian,                &
     &                  ww3_xcoord, ww3_ycoord, ww3_bath,               &
     &                  ww3_numx, ww3_numy,                             &
     &                  wrf_grids, parent_grid_ratio, parent_id,        &
     &                  output_ncfile, hydro_grids

#ifdef MPI
      CALL mpi_comm_rank (MyComm, MyRank, MyError)
      IF (MyRank.eq.0) THEN
#endif
      write(stdout,*)"================================================"
      write(stdout,*) ' Read input_file for SCRIP_COAWST Wrapper '
      write(stdout,*)"================================================"
#ifdef MPI
      END IF
      CALL mpi_barrier (MyComm, MyError)
#endif

      read (iunit, inputs)

!     Total number of grid interpolation combinations 
      Ngrids_comb_total=(Ngrids_roms*Ngrids_swan +                      &
     &                   Ngrids_roms*Ngrids_wrf  +                      &
     &                   Ngrids_swan*Ngrids_wrf  +                      & 
     &                   Ngrids_roms*Ngrids_ww3  +                      &
     &                   Ngrids_roms*Ngrids_hyd  +                      &
     &                   Ngrids_ww3*Ngrids_wrf)*2

#ifdef MPI
      CALL mpi_comm_rank (MyComm, MyRank, MyError)
      IF (MyRank.eq.0) THEN
#endif
      write(stdout,*) "Ngrid_roms=",Ngrids_roms
      write(stdout,*) "Ngrid_swan=",Ngrids_swan
      write(stdout,*) "Ngrid_ww3= ",Ngrids_ww3
      write(stdout,*) "Ngrid_wrf =",Ngrids_wrf
      write(stdout,*) "Ngrid_hyd =",Ngrids_hyd

      write(stdout,*) "Common netcdf file is: ",output_ncfile
      write(stdout,*) "Number of netcdf subgroups ", Ngrids_comb_total

      do i = 1,Ngrids_roms 
        write(*,10)"Input ROMS grid",i,"=",TRIM(ADJUSTL(roms_grids(i)))
      end do
      do i = 1,Ngrids_swan
        write(*,10)"Input SWAN grid",i,"=", TRIM(ADJUSTL(swan_coord(i)))
        write(*,11)"Input SWAN bath",i,"=", TRIM(ADJUSTL(swan_bath(i)))
        write(*,12)"Cartesian input In meters=", cartesian(i)
      end do 
      do i = 1,Ngrids_ww3
        write(*,10)"Input WW3 xgrid",i,"=", TRIM(ADJUSTL(ww3_xcoord(i)))
        write(*,10)"Input WW3 ygrid",i,"=", TRIM(ADJUSTL(ww3_ycoord(i)))
        write(*,11)"Input WW3 bath",i,"=", TRIM(ADJUSTL(ww3_bath(i)))
      end do 
      do i = 1,Ngrids_wrf
        write(*,10)"Input WRF grid ", i,"=", TRIM(ADJUSTL(wrf_grids(i)))
      end do 
      do i = 1,Ngrids_hyd
        write(*,10)"Input Hydro grid ", i,"=",                          &
     &              TRIM(ADJUSTL(hydro_grids(i)))
      end do 
      write(stdout,*)"================================================"

#ifdef MPI
      END IF
      CALL mpi_barrier (MyComm, MyError)
#endif

 10   FORMAT(A16, 1X, I1, A3, 1X, A)
 11   FORMAT(A16, 1X, I1, A3, 1X, A)
 12   FORMAT(A27, 1X, I20)

      end subroutine read_inputs 
