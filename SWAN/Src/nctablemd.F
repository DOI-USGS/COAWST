    module nctablemd
!/ ------------------------------------------------------------------- /
!/
!/
!/    08-dec-2009 : initial revision (sander.hulst@bmtargoss.com)
!/
!/    copyright 2009 bmtargoss
!/       no unauthorized use without permission.
!/
!  1. purpose :
!
!     nctable is a support module for agioncmd. it defines the supported
!     netcdf variables
!
!  2.attributes
!     --name-------------type-------description-----------------------
!       nctable          ncttype    struct array with fields
!         %standard_name char       see remarks a
!         %long_name     char       "
!         %units         char       "
!         %nctype        integer    nc_short, nc_byte etc.
!         %datamin       real       minimum expected data value
!         %datamax       real       maximum expected data value
!         %varid         int        nc90_var_id
        use netcdf
        implicit none
        private
        public                                       :: nctable, nctable_record
        public                                       :: get_nctable_record, set_nctable_convention_nautical
        public                                       :: pi
        real,    parameter                           :: pi  = 3.1415927
        real,    parameter                           :: r2d = 180. / pi
        real,    parameter                           :: d2r = pi / 180.
        integer, parameter                           :: variddum = -1
        logical                                      :: nautical = .true.
        type nctable_record
            character(128)                           :: name, standard_name
            character(128)                           :: long_name
            character(24)                            :: units
            integer                                  :: nctype
            real                                     :: datamin, datamax
            integer                                  :: varid
        end type nctable_record

        type (nctable_record)                        :: nctable(42)

        data nctable(1) / nctable_record("depth",&
                            "sea_floor_depth_below_sea_level",&
                            "depth below mean sea level", &
                            "m", NF90_FLOAT, 0, 0, variddum) /, &
             nctable(2) / nctable_record("xcur",&
                            "eastward_sea_water_velocity",&
                            "x component of current", &
                            "m s-1", NF90_SHORT, -20., 20., variddum) /, &
             nctable(3) / nctable_record("ycur",&
                            "northward_sea_water_velocity",&
                            "y component of current", &
                            "m s-1", NF90_SHORT, -20., 20., variddum) /, &
             nctable(4) / nctable_record("xwnd",&
                            "eastward_wind",&
                            "U-Component of Wind", &
                            "m s-1", NF90_SHORT, -100., 100., variddum) /, &
             nctable(5) / nctable_record("ywnd",&
                            "northward_wind",&
                            "V-Component of Wind", &
                            "m s-1", NF90_SHORT, -100., 100., variddum) /, &
             nctable(6) / nctable_record("astd",&
                            "none",&
                            "air sea temperature difference (k)", &
                            "k", NF90_BYTE, 0, 0, variddum) /, &
             nctable(7) / nctable_record("ustar",&
                            "none",&
                            "friction velocity", &
                            "m s-1", NF90_SHORT, 0, 20, variddum) /,&
             nctable(8) / nctable_record("hs",&
                            "sea_surface_wave_significant_height",&
                            "hs", &
                            "m", NF90_SHORT, 0., 50., variddum) /, &
             nctable(9) / nctable_record("L",&
                            "none",&
                            "L", &
                            "m", NF90_SHORT, 0, 0, variddum) /, &
             nctable(10) / nctable_record("theta0",&
                            "sea_surface_wave_from_direction",&
                            "theta0", &
                            "degrees", NF90_SHORT, 0., 360., variddum) /, &
             nctable(11) / nctable_record("tmm10",&
                            "sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment",&
                            "tm-10", &
                            "s", NF90_SHORT, 0., 50., variddum) /, &
             nctable(12) / nctable_record("thetam",&
                            "none",&
                            "thetam", &
                            "degrees", NF90_BYTE, 0, 0, variddum) /, &
             nctable(13) / nctable_record("tp",&
                            "none",&
                            "tp", &
                            "s", NF90_SHORT, 0., 50., variddum) /, &
             nctable(14) / nctable_record("thetap",&
                            "none",&
                            "thetap", &
                            "degrees", NF90_SHORT, 0., 360., variddum) /,&
             nctable(15) / nctable_record("fpl",&
                            "none",&
                            "peak frequency of wind sea part of spectrum", &
                            "s-1", NF90_BYTE, 0., 10., variddum) /,&
             nctable(16) / nctable_record("dpl",&
                            "none",&
                            "direction of peak frequency of wind sea part of spectrum", &
                            "degrees", NF90_BYTE, 0., 360., variddum) /,&
             nctable(17) / nctable_record("sign_wave_height_partitions",&
                            "none",&
                            "hs of partitions. does this mean that i have an extra dimension???", &
                            "m", NF90_SHORT, 0., 30., variddum) /, &
             nctable(18) / nctable_record("tp_partitions",&
                            "none",&
                            "relative peak periods of partitions of the spectrum", &
                            "s", NF90_BYTE, 0., 24., variddum) /,&
             nctable(19) / nctable_record("tp_length_partitions",&
                            "none",&
                            "peak wave lengths of partitions of the spectrum", &
                            "s", NF90_SHORT, 0, 0, variddum) /,&
             nctable(20) / nctable_record("theta0_partitions",&
                            "none",&
                            "mean wave direction of partitions of spectrum", &
                            "degrees", NF90_BYTE, 0., 360., variddum) /,&
             nctable(21) / nctable_record("spread_partitions",&
                            "none",&
                            "directional spread of partition of spectrum cf", &
                            "degrees", NF90_BYTE, 0, 81, variddum) /, &
             nctable(22) / nctable_record("wind_sea_fraction_partitions",&
                            "none",&
                            "wind sea fraction of partitions of spectrum", &
                            "1", NF90_SHORT, 0., 1., variddum) /,&
             nctable(23) / nctable_record("wind_sea_fraction",&
                            "none",&
                            "wind sea fraction of entire spectrum", &
                            "1", NF90_SHORT, 0., 1., variddum) /,&
             nctable(24) / nctable_record("npartitions",&
                            "none",&
                            "number of partitions found in spectrum", &
                            "1", NF90_BYTE, 0, 0, variddum) /,&
             nctable(25) / nctable_record("source_term_timestep",&
                            "none",&
                            "average timestep in the source term integration", &
                            "s", NF90_SHORT, 0, 0, variddum) /,&
             nctable(26) / nctable_record("cut_off_frequency",&
                            "none",&
                            "cut-off frequency", &
                            "s-1", NF90_SHORT, 0., 40., variddum) /,&
             nctable(27) / nctable_record("icec",&
                            "sea_ice_area_fraction",&
                            "ice cover", &
                            "1", NF90_BYTE, 0., 1., variddum) /,&
             nctable(28) / nctable_record("ssh",&
                            "sea_surface_height",&
                            "SSH", &
                            "m", NF90_SHORT, -15., 15., variddum) /, &
             nctable(29) / nctable_record("landmask",&
                            "land_binary_mask",&
                            "land cover (1=land, 0=sea)", &
                            "m", nf90_byte, 0, 0, variddum) /, &
             nctable(30) / nctable_record("tm02",&
                            "sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment", &
                            "tm02", &
                            "s", NF90_SHORT, 0., 50., variddum) /, &
             nctable(31) / nctable_record("tm01",&
                            "sea_surface_wave_mean_period_from_variance_spectral_density_first_frequency_moment", &
                            "tm01", &
                            "s", NF90_SHORT, 0., 50., variddum) /, &
             nctable(32) / nctable_record("tps",&
                            "none",&
                            "tps", &
                            "s", NF90_SHORT, 0., 50., variddum) /, &
             nctable(33) / nctable_record("spread",&
                            "none",&
                            "directional spreading", &
                            "degrees", NF90_SHORT, 0., 360., variddum) /, &
             nctable(34) / nctable_record("rtm01",&
                            "none", &
                            "rtm01", &
                            "s", NF90_SHORT, 0., 50., variddum) /, &
             nctable(35) / nctable_record("hswe",&
                            "sea_surface_swell_wave_significant_height",&
                            "wave height of swell part", &
                            "m", NF90_SHORT, 0., 50., variddum) /, &
             nctable(36) / nctable_record("rtmm10",&
                            "none",&
                            "rtm-10", &
                            "s", NF90_SHORT, 0., 50., variddum) /, &
             nctable(37) / nctable_record("botl",&
                            "none",&
                            "depth below still water level", &
                            "m", NF90_FLOAT, 0, 0, variddum) /, &
             nctable(38) / nctable_record("ubot",&
                            "none",&
                            "orbital velocity near bottom", &
                            "m s-1", NF90_SHORT, 0, 15, variddum) /, &
             nctable(39) / nctable_record("urms",&
                            "none",&
                            "rms of orbital velocity near bottom", &
                            "m s-1", NF90_SHORT, 0, 5, variddum) /, &
             nctable(40) / nctable_record("dhs",&
                            "none",&
                            "dHs", &
                            "m", NF90_SHORT, 0, 25, variddum) /, &
             nctable(41) / nctable_record("dtm",&
                            "none",&
                            "dTm", &
                            "s", NF90_SHORT, 0, 25, variddum) /, &
             nctable(42) / nctable_record("cdrag",&
                            "none",&
                            "Cdrag", &
                            "1", NF90_SHORT, 0, 1, variddum) /

    contains
        subroutine get_nctable_record(varname, trecord, found)
            character(len=*),           intent( in)     :: varname
            type(nctable_record),       intent(out)     :: trecord
            logical, optional,          intent(out)     :: found
            character(128)                              :: tmpname
            integer                                     :: i
            if ( present(found) ) found = .false.
            do i=1,size(nctable)
                if ( nctable(i)%name == varname ) then
                    trecord = nctable(i)
                    if ( .not.nautical ) then
                        tmpname = trecord%standard_name
                        call nautical_to_cartesian_name(tmpname, trecord%standard_name)
                    end if
                    if ( present(found) ) found = .true.
                    return
                end if
            end do
            if ( i > size(nctable) ) then
                print "(A)", "Could not find ", trim(varname), " in record table"
            end if
        end subroutine get_nctable_record

        subroutine set_nctable_convention_nautical(isnautical)
            logical,                    intent( in)    :: isnautical
            nautical = isnautical
        end subroutine set_nctable_convention_nautical

        subroutine nautical_to_cartesian_name(nname, cname)
            character(128),             intent( in)    :: nname
            character(128),             intent(out)    :: cname
            integer                                    :: i
            character(20)                              :: search = 'from_direction'
            i = index(nname, trim(search))
            if ( i < 1 ) then
                cname = nname
            else
                if ( i+len(trim(search)) < len(trim(nname)) ) then
                    cname = nname(1:i-1) // 'to_direction' // nname(i+len(trim(search)):len(trim(nname)))
                else
                    cname = nname(1:i-1) // 'to_direction'
                end if
            end if
        end subroutine nautical_to_cartesian_name
    end module nctablemd
