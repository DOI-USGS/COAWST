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
      use create_masks

      implicit none

!     local variables
      character(len=char_len)   :: inputfile

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
        call read_inputs()
      close(iunit)

!     Read the information from different grids 
      if (Ngrids_swan>0) then 
        call load_swan_grid()
      end if
      if (Ngrids_roms>0) then
        call load_roms_grid()
      end if
      if (Ngrids_wrf>0) then
        call load_wrf_grid()
      end if

!     Calculate interpolation weights between combinations of grids
      if ((Ngrids_roms>0).and.(Ngrids_swan>0)) then 
        call ocn2wav_mask()
      end if
      if ((Ngrids_roms>0).and.(Ngrids_wrf>0))  then 
        call ocn2atm_mask()
      end if

      if ((Ngrids_swan>0).and.(Ngrids_roms>0)) then 
        call wav2ocn_mask()
      end if 
      if ((Ngrids_swan>0).and.(Ngrids_wrf>0)) then
        call wav2atm_mask()
      end if

      if ((Ngrids_wrf>0).and.(Ngrids_roms>0)) then
        call atm2ocn_mask()
      end if
      if ((Ngrids_wrf>0).and.(Ngrids_swan>0)) then
        call atm2wav_mask()
      end if

      end program scrip_coawst

!======================================================================

      subroutine read_inputs()
      use kinds_mod
      use iounits
      use scripwrap_mod
 
      implicit none

!     local variables
      integer(int_kind) :: i

!     Allocate input variables
      namelist /inputs/ Ngrids_roms, Ngrids_swan, Ngrids_wrf,           &
     &                  roms_grids, swan_coord, swan_bath,              &
     &                  swan_numx, swan_numy, cartesian,                &
     &                  wrf_grids, parent_grid_ratio, parent_id,
     &                  output_ncfile

      write(stdout,*)"================================================"
      write(stdout,*) ' Read input_file for SCRIP_COAWST Wrapper '
      write(stdout,*)"================================================"
      read (iunit, inputs)

      write(stdout,*) "Ngrid_roms=",Ngrids_roms
      write(stdout,*) "Ngrid_swan=",Ngrids_swan
      write(stdout,*) "Ngrid_wrf =",Ngrids_wrf

!     Total number of grid interpolation combinations 
      Ngrids_comb_total=(Ngrids_roms*Ngrids_swan + 
     &                   Ngrids_roms*Ngrids_wrf  + 
     &                   Ngrids_swan*Ngrids_wrf)*2 

      write(stdout,*) "Common netcdf file",output_ncfile
      write(stdout,*) "Number of netcdf subgroups", Ngrids_comb_total

!      allocate(roms_grids(Ngrids_roms))
      do i = 1,Ngrids_roms 
        write(*,10)"Input ROMS grid",i,"=",TRIM(ADJUSTL(roms_grids(i)))
      end do
      do i = 1,Ngrids_swan
        write(*,10)"Input SWAN grid",i,"=", TRIM(ADJUSTL(swan_coord(i)))
        write(*,11)"Input SWAN bath",i,"=", TRIM(ADJUSTL(swan_bath(i)))
        write(*,12)"Cartesian input In meters=", cartesian(i)
      end do 
      do i = 1,Ngrids_wrf
        write(*,10)"Input WRF grid ", i,"=", TRIM(ADJUSTL(wrf_grids(i)))
      end do 
      write(stdout,*)"================================================"
 10   FORMAT(A16, 1X, I1, A3, 1X, A)
 11   FORMAT(A16, 1X, I1, A3, 1X, A)
 12   FORMAT(A27, 1X, I20)

      end subroutine read_inputs 
