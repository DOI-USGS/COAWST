module SnowInputSnicarMod

!!! read in required SNICAR snow albedo parameter datasets
!!! This module should be called in host model driver but archived here in NoahMP src cold

  use netcdf
  use Machine
  use NoahmpIOVarType, only : NoahmpIO_type

#ifdef _PARALLEL_
  use mpi
#endif

  implicit none

contains

  subroutine SnowInputSnicar(NoahmpIO)

! ------------------------ Code history -----------------------------------
! Code: T.-S. Lin, C. He, et al. (2025, JHM)
! -------------------------------------------------------------------------

    implicit none

    type(NoahmpIO_type), intent(inout) :: NoahmpIO

! local variables
    character(len=30) :: name
    integer           :: rank
    integer           :: ierr,iret
    integer           :: ncid, varid

! --------------------------------------------------------------------
    associate(                                                              &
              idx_Mie_snw_mx          => NoahmpIO%idx_Mie_snw_mx           ,& ! in,  number of effective radius indices used in Mie lookup table [idx]
              idx_T_max               => NoahmpIO%idx_T_max                ,& ! in,  maxiumum temperature index used in aging lookup table [idx]
              idx_Tgrd_max            => NoahmpIO%idx_Tgrd_max             ,& ! in,  maxiumum temperature gradient index used in aging lookup table [idx]
              idx_rhos_max            => NoahmpIO%idx_rhos_max             ,& ! in,  maxiumum snow density index used in aging lookup table [idx]
              snicar_numrad_snw       => NoahmpIO%snicar_numrad_snw        ,& ! in,  wavelength bands used in SNICAR snow albedo calculation
              snicar_snw_optics       => NoahmpIO%SNICAR_SNOWOPTICS_OPT    ,& ! in,  snow optics type using different refractive index databases in SNICAR
              snicar_dust_optics      => NoahmpIO%SNICAR_DUSTOPTICS_OPT    ,& ! in,  dust optics type for SNICAR snow albedo calculation
              snicar_solarspec        => NoahmpIO%SNICAR_SOLARSPEC_OPT     ,& ! in,  type of downward solar radiation spectrum for SNICAR snow albedo calculation
              snicar_optic_flnm       => NoahmpIO%snicar_optic_flnm        ,& ! in,  filename for SNICAR optics parameters
              snicar_age_flnm         => NoahmpIO%snicar_age_flnm          ,& ! in,  filename for snow aging parameters
              ss_alb_snw_drc          => NoahmpIO%ss_alb_snw_drc           ,& ! out, Mie single scatter albedos for direct-beam ice  
              asm_prm_snw_drc         => NoahmpIO%asm_prm_snw_drc          ,& ! out, asymmetry parameter of direct-beam ice
              ext_cff_mss_snw_drc     => NoahmpIO%ext_cff_mss_snw_drc      ,& ! out, mass extinction coefficient for direct-beam ice [m2/kg]
              ss_alb_snw_dfs          => NoahmpIO%ss_alb_snw_dfs           ,& ! out, Mie single scatter albedos for diffuse ice
              asm_prm_snw_dfs         => NoahmpIO%asm_prm_snw_dfs          ,& ! out, asymmetry parameter of diffuse ice
              ext_cff_mss_snw_dfs     => NoahmpIO%ext_cff_mss_snw_dfs      ,& ! out, mass extinction coefficient for diffuse ice [m2/kg]
              ss_alb_bc1              => NoahmpIO%ss_alb_bc1               ,& ! out, Mie single scatter albedos for hydrophillic BC
              asm_prm_bc1             => NoahmpIO%asm_prm_bc1              ,& ! out, asymmetry parameter for hydrophillic BC
              ext_cff_mss_bc1         => NoahmpIO%ext_cff_mss_bc1          ,& ! out, mass extinction coefficient for hydrophillic BC [m2/kg]
              ss_alb_bc2              => NoahmpIO%ss_alb_bc2               ,& ! out, Mie single scatter albedos for hydrophobic BC
              asm_prm_bc2             => NoahmpIO%asm_prm_bc2              ,& ! out, asymmetry parameter for hydrophobic BC
              ext_cff_mss_bc2         => NoahmpIO%ext_cff_mss_bc2          ,& ! out, mass extinction coefficient for hydrophobic BC [m2/kg]
              ss_alb_oc1              => NoahmpIO%ss_alb_oc1               ,& ! out, Mie single scatter albedos for hydrophillic OC
              asm_prm_oc1             => NoahmpIO%asm_prm_oc1              ,& ! out, asymmetry parameter for hydrophillic OC
              ext_cff_mss_oc1         => NoahmpIO%ext_cff_mss_oc1          ,& ! out, mass extinction coefficient for hydrophillic OC [m2/kg]
              ss_alb_oc2              => NoahmpIO%ss_alb_oc2               ,& ! out, Mie single scatter albedos for hydrophobic OC
              asm_prm_oc2             => NoahmpIO%asm_prm_oc2              ,& ! out, asymmetry parameter for hydrophobic OC
              ext_cff_mss_oc2         => NoahmpIO%ext_cff_mss_oc2          ,& ! out, mass extinction coefficient for hydrophobic OC [m2/kg]
              ss_alb_dst1             => NoahmpIO%ss_alb_dst1              ,& ! out, Mie single scatter albedos for dust species 1 
              asm_prm_dst1            => NoahmpIO%asm_prm_dst1             ,& ! out, asymmetry parameter for dust species 1
              ext_cff_mss_dst1        => NoahmpIO%ext_cff_mss_dst1         ,& ! out, mass extinction coefficient for dust species 1 [m2/kg]
              ss_alb_dst2             => NoahmpIO%ss_alb_dst2              ,& ! out, Mie single scatter albedos for dust species 2
              asm_prm_dst2            => NoahmpIO%asm_prm_dst2             ,& ! out, asymmetry parameter for dust species 2
              ext_cff_mss_dst2        => NoahmpIO%ext_cff_mss_dst2         ,& ! out, mass extinction coefficient for dust species 2 [m2/kg]
              ss_alb_dst3             => NoahmpIO%ss_alb_dst3              ,& ! out, Mie single scatter albedos for dust species 3
              asm_prm_dst3            => NoahmpIO%asm_prm_dst3             ,& ! out, asymmetry parameter for dust species 3
              ext_cff_mss_dst3        => NoahmpIO%ext_cff_mss_dst3         ,& ! out, mass extinction coefficient for dust species 3 [m2/kg]
              ss_alb_dst4             => NoahmpIO%ss_alb_dst4              ,& ! out, Mie single scatter albedos for dust species 4
              asm_prm_dst4            => NoahmpIO%asm_prm_dst4             ,& ! out, asymmetry parameter for dust species 4
              ext_cff_mss_dst4        => NoahmpIO%ext_cff_mss_dst4         ,& ! out, mass extinction coefficient for dust species 4 [m2/kg]
              ss_alb_dst5             => NoahmpIO%ss_alb_dst5              ,& ! out, Mie single scatter albedos for dust species 5
              asm_prm_dst5            => NoahmpIO%asm_prm_dst5             ,& ! out, asymmetry parameter for dust species 5
              ext_cff_mss_dst5        => NoahmpIO%ext_cff_mss_dst5         ,& ! out, mass extinction coefficient for dust species 5 [m2/kg]
              flx_wgt_dir             => NoahmpIO%flx_wgt_dir              ,& ! out, downward direct solar radiation spectral weights for wavelength band
              flx_wgt_dif             => NoahmpIO%flx_wgt_dif              ,& ! out, downward diffuse solar radiation spectral weights for wavelength band
              snowage_tau             => NoahmpIO%snowage_tau              ,& ! out, Snow aging parameters retrieved from lookup table [hour]
              snowage_kappa           => NoahmpIO%snowage_kappa            ,& ! out, Snow aging parameters retrieved from lookup table [unitless]
              snowage_drdt0           => NoahmpIO%snowage_drdt0             & ! out, Snow aging parameters retrieved from lookup table [m2 kg-1 hr-1]
             )
! ----------------------------------------------------------------------

#ifdef _PARALLEL_  
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    if (ierr /= MPI_SUCCESS) stop "MPI_COMM_RANK"
#else
    rank = 0
#endif

    !=================== read in SNICAR snow and aerosol optics parameters =========================

    ! Open the NetCDF file.
    if (rank == 0) write(*,'("Snicar SnowOptics init: ''", A, "''")') trim(snicar_optic_flnm)
#ifdef _PARALLEL_
    ierr = nf90_open_par(snicar_optic_flnm, NF90_NOWRITE, MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
#else
    ierr = nf90_open(snicar_optic_flnm, NF90_NOWRITE, ncid)
#endif
    if (ierr /= 0) then
       write(*,'("read_snicar_data:  Problem opening file: ''", A, "''")') trim(snicar_optic_flnm)
#ifdef _PARALLEL_
       call mpi_finalize(ierr)
       if (ierr /= 0) write(*, '("Problem with MPI_finalize.")')
#endif
       stop
    endif

    if (snicar_numrad_snw==5) then

       ! mid-latitude winter spectrum
       if (snicar_solarspec == 1) then

          ! flux weights/spectrum
          name = "flx_wgt_dir5_mlw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, flx_wgt_dir, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "flx_wgt_dif5_mlw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, flx_wgt_dif, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! BC species 1 Mie parameters, uncoated BC, same as bc2 before BC-snow internal mixing
          name = "ss_alb_bcphob_dif_mlw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_bcphob_dif_mlw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_bcphob_dif_mlw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! BC species 2 Mie parameters, uncoated BC
          name = "ss_alb_bcphob_dif_mlw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_bcphob_dif_mlw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_bcphob_dif_mlw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! OC species 1 Mie parameters, uncoated OC, same as oc2 before OC-snow internal mixing
          name = "ss_alb_ocphob_dif_mlw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_ocphob_dif_mlw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_ocphob_dif_mlw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! OC species 2 Mie parameters, uncoated OC
          name = "ss_alb_ocphob_dif_mlw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_ocphob_dif_mlw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_ocphob_dif_mlw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

           ! ice refractive index options
          if (snicar_snw_optics == 1) then  ! Warren (1984)
            name = "ss_alb_ice_wrn84_dir_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn84_dir_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn84_dir_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_wrn84_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn84_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn84_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_snw_optics == 2) then ! Warren and Brandt (2008)
            name = "ss_alb_ice_wrn08_dir_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn08_dir_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn08_dir_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_wrn08_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn08_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn08_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_snw_optics == 3) then ! Picard et al (2016)
            name = "ss_alb_ice_pic16_dir_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_pic16_dir_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_pic16_dir_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_pic16_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_pic16_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_pic16_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          endif

          ! dust optical properties
          if(snicar_dust_optics == 1) then ! Saharan dust (Balkanski et al., 2007, central hematite)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_sah_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_sah_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_sah_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_sah_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust02_sah_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_sah_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_sah_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust03_sah_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_sah_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_sah_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust04_sah_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_sah_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_sah_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_sah_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_sah_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_dust_optics == 2) then  ! San Juan Mountains, CO (Skiles et al, 2017)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_col_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_col_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_col_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_col_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
       
            name = "asm_prm_dust02_col_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_col_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_col_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust03_col_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_col_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_col_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust04_col_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_col_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_col_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_col_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_col_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_dust_optics == 3) then  ! Greenland (Polashenski et al., 2015, central absorptivity)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_gre_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_gre_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_gre_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_gre_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust02_gre_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_gre_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_gre_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust03_gre_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_gre_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_gre_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust04_gre_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_gre_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_gre_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_gre_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_gre_dif_mlw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          endif

       ! mid-latitude summer spectrum
       elseif (snicar_solarspec == 2) then
          ! flux weights/spectrum
          name = "flx_wgt_dir5_mls"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, flx_wgt_dir, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "flx_wgt_dif5_mls"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, flx_wgt_dif, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! BC species 1 Mie parameters, uncoated BC, same as bc2 before BC-snow internal mixing
          name = "ss_alb_bcphob_dif_mls"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_bcphob_dif_mls"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_bcphob_dif_mls"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! BC species 2 Mie parameters, uncoated BC
          name = "ss_alb_bcphob_dif_mls"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_bcphob_dif_mls"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_bcphob_dif_mls"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! OC species 1 Mie parameters, uncoated OC, same as oc2 before OC-snow internal mixing
          name = "ss_alb_ocphob_dif_mls"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_ocphob_dif_mls"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_ocphob_dif_mls"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! OC species 2 Mie parameters, uncoated OC
          name = "ss_alb_ocphob_dif_mls"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_ocphob_dif_mls"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_ocphob_dif_mls"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          if (snicar_snw_optics == 1) then  ! Warren (1984)
            name = "ss_alb_ice_wrn84_dir_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn84_dir_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn84_dir_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_wrn84_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn84_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn84_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_snw_optics == 2) then  ! Warren and Brandt (2008)
            name = "ss_alb_ice_wrn08_dir_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn08_dir_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn08_dir_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_wrn08_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn08_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn08_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_snw_optics == 3) then  ! Picard et al (2016)
            name = "ss_alb_ice_pic16_dir_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_pic16_dir_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_pic16_dir_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_pic16_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_pic16_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_pic16_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          endif

          if (snicar_dust_optics == 1) then ! Saharan dust (Balkanski et al., 2007, central hematite)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_sah_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_sah_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_sah_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_sah_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
    
            name = "asm_prm_dust02_sah_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_sah_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_sah_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
    
            name = "asm_prm_dust03_sah_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_sah_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_sah_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
    
            name = "asm_prm_dust04_sah_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_sah_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_sah_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_sah_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_sah_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_dust_optics == 2) then  ! San Juan Mountains, CO (Skiles et al, 2017)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_col_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_col_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_col_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_col_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust02_col_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_col_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_col_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust03_col_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_col_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_col_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust04_col_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_col_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_col_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_col_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_col_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_dust_optics == 3) then  ! Greenland (Polashenski et al., 2015, central absorptivity)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_gre_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_gre_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_gre_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_gre_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust02_gre_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_gre_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_gre_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust03_gre_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_gre_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_gre_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust04_gre_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_gre_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_gre_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_gre_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_gre_dif_mls"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          endif

       ! sub-Arctic winter spectrum
       elseif (snicar_solarspec == 3) then
          ! flux weights/spectrum
          name = "flx_wgt_dir5_saw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, flx_wgt_dir, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "flx_wgt_dif5_saw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, flx_wgt_dif, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! BC species 1 Mie parameters, uncoated BC, same as bc2 before BC-snow internal mixing
          name = "ss_alb_bcphob_dif_saw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_bcphob_dif_saw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_bcphob_dif_saw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! BC species 2 Mie parameters, uncoated BC
          name = "ss_alb_bcphob_dif_saw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_bcphob_dif_saw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_bcphob_dif_saw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! OC species 1 Mie parameters, uncoated OC, same as oc2 before OC-snow internal mixing
          name = "ss_alb_ocphob_dif_saw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_ocphob_dif_saw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_ocphob_dif_saw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! OC species 2 Mie parameters, uncoated OC
          name = "ss_alb_ocphob_dif_saw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_ocphob_dif_saw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_ocphob_dif_saw"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          if (snicar_snw_optics == 1) then  ! Warren (1984)
            name = "ss_alb_ice_wrn84_dir_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn84_dir_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn84_dir_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_wrn84_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn84_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn84_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_snw_optics == 2) then  ! Warren and Brandt (2008)
            name = "ss_alb_ice_wrn08_dir_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn08_dir_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn08_dir_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_wrn08_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn08_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn08_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_snw_optics == 3) then  ! Picard et al (2016)
            name = "ss_alb_ice_pic16_dir_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_pic16_dir_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_pic16_dir_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_pic16_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_pic16_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_pic16_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          endif

          if (snicar_dust_optics == 1) then ! Saharan dust (Balkanski et al., 2007, central hematite)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_sah_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_sah_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_sah_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_sah_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
    
            name = "asm_prm_dust02_sah_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_sah_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_sah_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
    
            name = "asm_prm_dust03_sah_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_sah_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_sah_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
    
            name = "asm_prm_dust04_sah_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_sah_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_sah_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_sah_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_sah_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_dust_optics == 2) then  ! San Juan Mountains, CO (Skiles et al, 2017)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_col_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_col_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_col_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_col_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust02_col_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_col_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_col_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust03_col_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_col_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_col_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust04_col_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_col_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_col_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_col_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_col_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_dust_optics == 3) then  ! Greenland (Polashenski et al., 2015, central absorptivity)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_gre_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_gre_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_gre_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_gre_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust02_gre_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_gre_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_gre_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust03_gre_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_gre_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_gre_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust04_gre_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_gre_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_gre_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_gre_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_gre_dif_saw"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          endif

       ! sub-Arctic summer spectrum
       elseif (snicar_solarspec == 4) then
          ! flux weights/spectrum
          name = "flx_wgt_dir5_sas"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, flx_wgt_dir, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "flx_wgt_dif5_sas"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, flx_wgt_dif, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! BC species 1 Mie parameters, uncoated BC, same as bc2 before BC-snow internal mixing
          name = "ss_alb_bcphob_dif_sas"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_bcphob_dif_sas"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_bcphob_dif_sas"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! BC species 2 Mie parameters, uncoated BC
          name = "ss_alb_bcphob_dif_sas"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_bcphob_dif_sas"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_bcphob_dif_sas"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! OC species 1 Mie parameters, uncoated OC, same as oc2 before OC-snow internal mixing
          name = "ss_alb_ocphob_dif_sas"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_ocphob_dif_sas"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_ocphob_dif_sas"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! OC species 2 Mie parameters, uncoated OC
          name = "ss_alb_ocphob_dif_sas"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_ocphob_dif_sas"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_ocphob_dif_sas"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          if (snicar_snw_optics == 1) then  ! Warren (1984)
            name = "ss_alb_ice_wrn84_dir_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn84_dir_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn84_dir_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_wrn84_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn84_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn84_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_snw_optics == 2) then  ! Warren and Brandt (2008)
            name = "ss_alb_ice_wrn08_dir_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn08_dir_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn08_dir_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_wrn08_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn08_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn08_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_snw_optics == 3) then  ! Picard et al (2016)
            name = "ss_alb_ice_pic16_dir_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_pic16_dir_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_pic16_dir_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_pic16_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_pic16_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_pic16_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          endif

          if (snicar_dust_optics == 1) then ! Saharan dust (Balkanski et al., 2007, central hematite)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_sah_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_sah_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_sah_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_sah_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
    
            name = "asm_prm_dust02_sah_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_sah_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_sah_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
    
            name = "asm_prm_dust03_sah_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_sah_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_sah_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
    
            name = "asm_prm_dust04_sah_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_sah_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_sah_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_sah_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_sah_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_dust_optics == 2) then  ! San Juan Mountains, CO (Skiles et al, 2017)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_col_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_col_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_col_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_col_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust02_col_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_col_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_col_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust03_col_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_col_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_col_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust04_col_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_col_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_col_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_col_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_col_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_dust_optics == 3) then  ! Greenland (Polashenski et al., 2015, central absorptivity)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_gre_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_gre_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_gre_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_gre_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust02_gre_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_gre_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_gre_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust03_gre_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_gre_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_gre_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust04_gre_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_gre_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_gre_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_gre_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_gre_dif_sas"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          endif

       ! Summit,Greenland,summer spectrum
       elseif (snicar_solarspec == 5) then
          ! flux weights/spectrum
          name = "flx_wgt_dir5_smm"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, flx_wgt_dir, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "flx_wgt_dif5_smm"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, flx_wgt_dif, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! BC species 1 Mie parameters, uncoated BC, same as bc2 before BC-snow internal mixing
          name = "ss_alb_bcphob_dif_smm"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_bcphob_dif_smm"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_bcphob_dif_smm"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! BC species 2 Mie parameters, uncoated BC
          name = "ss_alb_bcphob_dif_smm"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_bcphob_dif_smm"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_bcphob_dif_smm"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! OC species 1 Mie parameters, uncoated OC, same as oc2 before OC-snow internal mixing
          name = "ss_alb_ocphob_dif_smm"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_ocphob_dif_smm"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_ocphob_dif_smm"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! OC species 2 Mie parameters, uncoated OC
          name = "ss_alb_ocphob_dif_smm"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_ocphob_dif_smm"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_ocphob_dif_smm"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          if (snicar_snw_optics == 1) then  ! Warren (1984)
            name = "ss_alb_ice_wrn84_dir_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn84_dir_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn84_dir_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_wrn84_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn84_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn84_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_snw_optics == 2) then  ! Warren and Brandt (2008)
            name = "ss_alb_ice_wrn08_dir_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn08_dir_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn08_dir_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_wrn08_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn08_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn08_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_snw_optics == 3) then  ! Picard et al (2016)
            name = "ss_alb_ice_pic16_dir_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_pic16_dir_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_pic16_dir_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_pic16_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_pic16_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_pic16_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          endif

          if (snicar_dust_optics == 1) then ! Saharan dust (Balkanski et al., 2007, central hematite)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_sah_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_sah_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_sah_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_sah_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
    
            name = "asm_prm_dust02_sah_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_sah_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_sah_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
    
            name = "asm_prm_dust03_sah_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_sah_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_sah_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
    
            name = "asm_prm_dust04_sah_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_sah_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_sah_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_sah_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_sah_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_dust_optics == 2) then  ! San Juan Mountains, CO (Skiles et al, 2017)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_col_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_col_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_col_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_col_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust02_col_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_col_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_col_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust03_col_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_col_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_col_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust04_col_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_col_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_col_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_col_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_col_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_dust_optics == 3) then  ! Greenland (Polashenski et al., 2015, central absorptivity)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_gre_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_gre_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_gre_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_gre_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust02_gre_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_gre_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_gre_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust03_gre_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_gre_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_gre_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust04_gre_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_gre_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_gre_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_gre_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_gre_dif_smm"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          endif
       
       ! High Mountain summer spectrum
       elseif (snicar_solarspec == 6) then
          ! flux weights/spectrum
          name = "flx_wgt_dir5_hmn"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, flx_wgt_dir, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "flx_wgt_dif5_hmn"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, flx_wgt_dif, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! BC species 1 Mie parameters, uncoated BC, same as bc2 before BC-snow internal mixing
          name = "ss_alb_bcphob_dif_hmn"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_bcphob_dif_hmn"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_bcphob_dif_hmn"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_bc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! BC species 2 Mie parameters, uncoated BC
          name = "ss_alb_bcphob_dif_hmn"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_bcphob_dif_hmn"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_bcphob_dif_hmn"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_bc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! OC species 1 Mie parameters, uncoated OC, same as oc2 before OC-snow internal mixing
          name = "ss_alb_ocphob_dif_hmn"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_ocphob_dif_hmn"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_ocphob_dif_hmn"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_oc1, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          ! OC species 2 Mie parameters, uncoated OC
          name = "ss_alb_ocphob_dif_hmn"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ss_alb_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "asm_prm_ocphob_dif_hmn"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, asm_prm_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          name = "ext_cff_mss_ocphob_dif_hmn"
          iret = nf90_inq_varid(ncid,  name,  varid)
          if (iret == 0) then
            ierr = nf90_get_var(ncid, varid, ext_cff_mss_oc2, start=(/1/), count=(/snicar_numrad_snw/))
          else
            write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
          endif

          if (snicar_snw_optics == 1) then  ! Warren (1984)
            name = "ss_alb_ice_wrn84_dir_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn84_dir_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn84_dir_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_wrn84_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn84_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn84_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_snw_optics == 2) then  ! Warren and Brandt (2008)
            name = "ss_alb_ice_wrn08_dir_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn08_dir_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn08_dir_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_wrn08_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_wrn08_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_wrn08_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_snw_optics == 3) then  ! Picard et al (2016)
            name = "ss_alb_ice_pic16_dir_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_pic16_dir_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_pic16_dir_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ss_alb_ice_pic16_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_ice_pic16_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_ice_pic16_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          endif

          if (snicar_dust_optics == 1) then ! Saharan dust (Balkanski et al., 2007, central hematite)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_sah_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_sah_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_sah_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_sah_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
    
            name = "asm_prm_dust02_sah_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_sah_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_sah_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
    
            name = "asm_prm_dust03_sah_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_sah_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_sah_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif
    
            name = "asm_prm_dust04_sah_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_sah_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_sah_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_sah_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_sah_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_dust_optics == 2) then  ! San Juan Mountains, CO (Skiles et al, 2017)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_col_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_col_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_col_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_col_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust02_col_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_col_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_col_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust03_col_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_col_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_col_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust04_col_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_col_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_col_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_col_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_col_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          elseif (snicar_dust_optics == 3) then  ! Greenland (Polashenski et al., 2015, central absorptivity)
            ! dust species 1 Mie parameters
            name = "ss_alb_dust01_gre_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust01_gre_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust01_gre_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 2 Mie parameters
            name = "ss_alb_dust02_gre_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust02_gre_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust02_gre_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 3 Mie parameters
            name = "ss_alb_dust03_gre_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust03_gre_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust03_gre_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 4 Mie parameters
            name = "ss_alb_dust04_gre_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust04_gre_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust04_gre_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            ! dust species 5 Mie parameters
            name = "ss_alb_dust05_gre_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "asm_prm_dust05_gre_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

            name = "ext_cff_mss_dust05_gre_dif_hmn"
            iret = nf90_inq_varid(ncid,  name,  varid)
            if (iret == 0) then
              ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
            else
              write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
            endif

          endif
      
       endif ! end of snicar_solarspec
 
    endif ! end snicar five bands

    if (snicar_numrad_snw==480) then
   
       ! BC species 1 Mie parameters, uncoated BC, same as bc2 before BC-snow internal mixing
       name = "ss_alb_bcphob"
       iret = nf90_inq_varid(ncid,  name,  varid)
       if (iret == 0) then
         ierr = nf90_get_var(ncid, varid, ss_alb_bc1, start=(/1/), count=(/snicar_numrad_snw/))
       else
         write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
       endif

       name = "asm_prm_bcphob"
       iret = nf90_inq_varid(ncid,  name,  varid)
       if (iret == 0) then
         ierr = nf90_get_var(ncid, varid, asm_prm_bc1, start=(/1/), count=(/snicar_numrad_snw/))
       else
         write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
       endif

       name = "ext_cff_mss_bcphob"
       iret = nf90_inq_varid(ncid,  name,  varid)
       if (iret == 0) then
         ierr = nf90_get_var(ncid, varid, ext_cff_mss_bc1, start=(/1/), count=(/snicar_numrad_snw/))
       else
         write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
       endif

       ! BC species 2 Mie parameters, uncoated BC
       name = "ss_alb_bcphob"
       iret = nf90_inq_varid(ncid,  name,  varid)
       if (iret == 0) then
         ierr = nf90_get_var(ncid, varid, ss_alb_bc2, start=(/1/), count=(/snicar_numrad_snw/))
       else
         write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
       endif

       name = "asm_prm_bcphob"
       iret = nf90_inq_varid(ncid,  name,  varid)
       if (iret == 0) then
         ierr = nf90_get_var(ncid, varid, asm_prm_bc2, start=(/1/), count=(/snicar_numrad_snw/))
       else
         write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
       endif

       name = "ext_cff_mss_bcphob"
       iret = nf90_inq_varid(ncid,  name,  varid)
       if (iret == 0) then
         ierr = nf90_get_var(ncid, varid, ext_cff_mss_bc2, start=(/1/), count=(/snicar_numrad_snw/))
       else
         write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
       endif

       ! OC species 1 Mie parameters, uncoated OC, same as oc2 before OC-snow internal mixing
       name = "ss_alb_ocphob"
       iret = nf90_inq_varid(ncid,  name,  varid)
       if (iret == 0) then
         ierr = nf90_get_var(ncid, varid, ss_alb_oc1, start=(/1/), count=(/snicar_numrad_snw/))
       else
         write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
       endif

       name = "asm_prm_ocphob"
       iret = nf90_inq_varid(ncid,  name,  varid)
       if (iret == 0) then
         ierr = nf90_get_var(ncid, varid, asm_prm_oc1, start=(/1/), count=(/snicar_numrad_snw/))
       else
         write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
       endif

       name = "ext_cff_mss_ocphob"
       iret = nf90_inq_varid(ncid,  name,  varid)
       if (iret == 0) then
         ierr = nf90_get_var(ncid, varid, ext_cff_mss_oc1, start=(/1/), count=(/snicar_numrad_snw/))
       else
         write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
       endif

       ! OC species 2 Mie parameters, uncoated OC
       name = "ss_alb_ocphob"
       iret = nf90_inq_varid(ncid,  name,  varid)
       if (iret == 0) then
         ierr = nf90_get_var(ncid, varid, ss_alb_oc2, start=(/1/), count=(/snicar_numrad_snw/))
       else
         write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
       endif

       name = "asm_prm_ocphob"
       iret = nf90_inq_varid(ncid,  name,  varid)
       if (iret == 0) then
         ierr = nf90_get_var(ncid, varid, asm_prm_oc2, start=(/1/), count=(/snicar_numrad_snw/))
       else
         write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
       endif

       name = "ext_cff_mss_ocphob"
       iret = nf90_inq_varid(ncid,  name,  varid)
       if (iret == 0) then
         ierr = nf90_get_var(ncid, varid, ext_cff_mss_oc2, start=(/1/), count=(/snicar_numrad_snw/))
       else
         write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
       endif

       ! snow optical properties derived from different ice refractive index dataset
       ! same value for direct and diffuse due to high spectral res without spectra averaging in database
       if (snicar_snw_optics == 1) then  ! Warren (1984)
         name = "ss_alb_ice_wrn84"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_ice_wrn84"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_ice_wrn84"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ss_alb_ice_wrn84"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_ice_wrn84"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_ice_wrn84"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

       elseif (snicar_snw_optics == 2) then  ! Warren and Brandt (2008)
         name = "ss_alb_ice_wrn08"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_ice_wrn08"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_ice_wrn08"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ss_alb_ice_wrn08"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_ice_wrn08"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_ice_wrn08"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

       elseif (snicar_snw_optics == 3) then  ! Picard et al (2016)
         name = "ss_alb_ice_pic16"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_ice_pic16"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_ice_pic16"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_drc, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ss_alb_ice_pic16"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_ice_pic16"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_ice_pic16"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_snw_dfs, start=(/1,1/), count=(/idx_Mie_snw_mx,snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

       endif

      ! dust optical properties
       if (snicar_dust_optics == 1) then  ! Saharan dust (Balkanski et al., 2007, central hematite)
         ! dust species 1 Mie parameters
         name = "ss_alb_dust01_sah"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_dust01_sah"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_dust01_sah"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         ! dust species 2 Mie parameters
         name = "ss_alb_dust02_sah"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_dust02_sah"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_dust02_sah"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         ! dust species 3 Mie parameters
         name = "ss_alb_dust03_sah"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_dust03_sah"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_dust03_sah"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         ! dust species 4 Mie parameters
         name = "ss_alb_dust04_sah"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_dust04_sah"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_dust04_sah"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         ! dust species 5 Mie parameters
         name = "ss_alb_dust05_sah"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_dust05_sah"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_dust05_sah"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

       elseif (snicar_dust_optics == 2) then  ! San Juan Mountains, CO (Skiles et al, 2017)
         ! dust species 1 Mie parameters
         name = "ss_alb_dust01_col"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_dust01_col"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_dust01_col"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         ! dust species 2 Mie parameters
         name = "ss_alb_dust02_col"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_dust02_col"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_dust02_col"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         ! dust species 3 Mie parameters
         name = "ss_alb_dust03_col"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_dust03_col"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_dust03_col"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         ! dust species 4 Mie parameters
         name = "ss_alb_dust04_col"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_dust04_col"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_dust04_col"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         ! dust species 5 Mie parameters
         name = "ss_alb_dust05_col"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_dust05_col"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_dust05_col"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

       elseif (snicar_dust_optics == 3) then  ! Greenland (Polashenski et al., 2015, central absorptivity)
         ! dust species 1 Mie parameters
         name = "ss_alb_dust01_gre"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_dst1, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_dust01_gre"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_dst1, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_dust01_gre"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst1, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         ! dust species 2 Mie parameters
         name = "ss_alb_dust02_gre"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_dst2, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_dust02_gre"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_dst2, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_dust02_gre"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst2, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         ! dust species 3 Mie parameters
         name = "ss_alb_dust03_gre"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_dst3, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_dust03_gre"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_dst3, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_dust03_gre"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst3, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         ! dust species 4 Mie parameters
         name = "ss_alb_dust04_gre"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_dst4, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_dust04_gre"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_dst4, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_dust04_gre"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst4, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         ! dust species 5 Mie parameters
         name = "ss_alb_dust05_gre"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ss_alb_dst5, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "asm_prm_dust05_gre"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, asm_prm_dst5, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "ext_cff_mss_dust05_gre"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, ext_cff_mss_dst5, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

       endif

       ! downward solar radiation spectral weights for 480-band
       if (snicar_solarspec == 1) then     ! mid-latitude winter
         name = "flx_wgt_dir480_mlw"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, flx_wgt_dir, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

         name = "flx_wgt_dif480_mlw"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, flx_wgt_dif, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

       elseif (snicar_solarspec == 2) then ! mid-latitude summer
         name = "flx_wgt_dir480_mls"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, flx_wgt_dir, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif
    
         name = "flx_wgt_dif480_mls"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, flx_wgt_dif, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

       elseif (snicar_solarspec == 3) then ! sub-Arctic winter
         name = "flx_wgt_dir480_saw"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, flx_wgt_dir, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif
    
         name = "flx_wgt_dif480_saw"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, flx_wgt_dif, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

       elseif (snicar_solarspec == 4) then ! sub-Arctic summer
         name = "flx_wgt_dir480_sas"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, flx_wgt_dir, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif
    
         name = "flx_wgt_dif480_sas"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, flx_wgt_dif, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

       elseif (snicar_solarspec == 5) then ! Summit,Greenland,summer
         name = "flx_wgt_dir480_smm"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, flx_wgt_dir, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif
    
         name = "flx_wgt_dif480_smm"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, flx_wgt_dif, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

 
       elseif (snicar_solarspec == 6) then ! High Mountain summer
         name = "flx_wgt_dir480_hmn"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, flx_wgt_dir, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif
    
         name = "flx_wgt_dif480_hmn"
         iret = nf90_inq_varid(ncid,  name,  varid)
         if (iret == 0) then
           ierr = nf90_get_var(ncid, varid, flx_wgt_dif, start=(/1/), count=(/snicar_numrad_snw/))
         else
           write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
         endif

       endif

   endif ! end of 480-band read in

    ! Close the NetCDF file
    ierr = nf90_close(ncid)
    if (ierr /= 0) stop "MODULE_NOAHLSM_HRLDAS_INPUT:  read_snicar_data:  NF90_CLOSE"

    write(*,*) 'Successfully read snow optical properties'
    write(*,*) 'band numbers', snicar_numrad_snw
    write (*,*) 'SNICAR: Mie single scatter albedos for direct-beam ice, rds=100um: ', &
             ss_alb_snw_drc(71,1), ss_alb_snw_drc(71,2), ss_alb_snw_drc(71,3),     &
             ss_alb_snw_drc(71,4), ss_alb_snw_drc(71,5)
    write (*,*) 'SNICAR: Mie single scatter albedos for diffuse ice, rds=100um: ',     &
             ss_alb_snw_dfs(71,1), ss_alb_snw_dfs(71,2), ss_alb_snw_dfs(71,3),     &
             ss_alb_snw_dfs(71,4), ss_alb_snw_dfs(71,5)
    write (*,*) 'SNICAR: Mie single scatter albedos for hydrophillic BC (before): ', &
             ss_alb_bc1(1), ss_alb_bc1(2), ss_alb_bc1(3), ss_alb_bc1(4), ss_alb_bc1(5)
    write (*,*) 'SNICAR: Mie single scatter albedos for hydrophobic BC: ', &
             ss_alb_bc2(1), ss_alb_bc2(2), ss_alb_bc2(3), ss_alb_bc2(4), ss_alb_bc2(5)
    write (*,*) 'SNICAR: Mie single scatter albedos for dust species 1: ', &
             ss_alb_dst1(1), ss_alb_dst1(2), ss_alb_dst1(3), ss_alb_dst1(4), ss_alb_dst1(5)
    write (*,*) 'SNICAR: Mie single scatter albedos for dust species 2: ', &
             ss_alb_dst2(1), ss_alb_dst2(2), ss_alb_dst2(3), ss_alb_dst2(4), ss_alb_dst2(5)
    write (*,*) 'SNICAR: Mie single scatter albedos for dust species 3: ', &
             ss_alb_dst3(1), ss_alb_dst3(2), ss_alb_dst3(3), ss_alb_dst3(4), ss_alb_dst3(5)
    write (*,*) 'SNICAR: Mie single scatter albedos for dust species 4: ', &
             ss_alb_dst4(1), ss_alb_dst4(2), ss_alb_dst4(3), ss_alb_dst4(4), ss_alb_dst4(5)


    !=================== read in SNICAR aging parameters =========================

    ! Open the NetCDF file.
    if (rank == 0) write(*,'("Snicar SnowAge init: ''", A, "''")') trim(snicar_age_flnm)
#ifdef _PARALLEL_
    ierr = nf90_open_par(snicar_age_flnm, NF90_NOWRITE, MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
#else
    ierr = nf90_open(snicar_age_flnm, NF90_NOWRITE, ncid)
#endif
    if (ierr /= 0) then
       write(*,'("read_snicar_data:  Problem opening file: ''", A, "''")') trim(snicar_age_flnm)
#ifdef _PARALLEL_
       call mpi_finalize(ierr)
       if (ierr /= 0) write(*, '("Problem with MPI_finalize.")')
#endif
       stop
    endif

    name = "tau"
    iret = nf90_inq_varid(ncid,  name,  varid)
    if (iret == 0) then
      ierr = nf90_get_var(ncid, varid, snowage_tau, start=(/1,1,1/), count=(/idx_rhos_max,idx_Tgrd_max,idx_T_max/))
    else
      write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
    endif

    name = "kappa"
    iret = nf90_inq_varid(ncid,  name,  varid)
    if (iret == 0) then
      ierr = nf90_get_var(ncid, varid, snowage_kappa, start=(/1,1,1/), count=(/idx_rhos_max,idx_Tgrd_max,idx_T_max/))
    else
      write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
    endif

    name = "drdsdt0"
    iret = nf90_inq_varid(ncid,  name,  varid)
    if (iret == 0) then
      ierr = nf90_get_var(ncid, varid, snowage_drdt0, start=(/1,1,1/), count=(/idx_rhos_max,idx_Tgrd_max,idx_T_max/))
    else
      write(*,*) "MODULE_HRLDAS_NETCDF_IO:  Problem finding variable '"//trim(name)//"' in NetCDF file."
    endif

    ! Close the NetCDF file
    ierr = nf90_close(ncid)
    if (ierr /= 0) stop "MODULE_NOAHLSM_HRLDAS_INPUT:  read_snicar_data:  NF90_CLOSE"

    write(*,*) 'Successfully read snow aging properties'

    ! print some diagnostics:
    write (*,*) 'SNICAR: snowage tau for T=263.15K, dTdz = 100 K/m, rhos = 150 kg/m3: ', snowage_tau(3,11,9)
    write (*,*) 'SNICAR: snowage kappa for T=263.15K, dTdz = 100 K/m, rhos = 150 kg/m3: ', snowage_kappa(3,11,9)
    write (*,*) 'SNICAR: snowage dr/dt_0 for T=263.15K, dTdz = 100 K/m, rhos = 150 kg/m3: ', snowage_drdt0(3,11,9)

    end associate

  end subroutine SnowInputSnicar

end module SnowInputSnicarMod

