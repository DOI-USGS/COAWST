!/ ------------------------------------------------------------------- /
!/
!/
!/    03-dec-2009 : initial revision (sander.hulst@bmtargoss.com)
!/
!/    copyright 2009 bmtargoss
!/       no unauthorized use without permission.
!/
!  1. purpose :
!
!     agioncmd is a module that provides methods to write cf1.5 compliant
!     netcdf files containing integrated parameters and (non-cf compliant)
!     spectra.
!     furthermore, some methods are available to read data from cf-1.5 and
!     coads compliant files. The supported time units are limited to
!     seconds or days since 1970-01-01 (so no foam, but that's not a geographic
!     grid anyway, it's a mercator). If this is a problem, then the UDUNITS
!     utScan() and utIsTime() can be used to make this determination.
!
!     For more information:
!       http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.5
!
!  2. method :
!
!  3.a.types
!     --name-------------type-------description-----------------------
!       mapgrid          intgtype   struct with fields
!         %longitude     r.a        x-coordinates (dim)
!         %latitude      r.a        y-coordinates (dim)
!         %nx            int        number of columns
!         %ny            int        number of rows
!         %dx            real       horizontal cellsize
!         %dy            real       vertical cellsize
!         %lunit         char       likely unit (meter, degrees)
!
!       mmapgrid         mmapgtype  two dimensional lon lat
!         %longitude     r.a        x-coordinates (dim)
!         %latitude      r.a        y-coordinates (dim)
!         %nx            int        number of columns
!         %ny            int        number of rows
!         %dx            real       horizontal cellsize
!         %dy            real       vertical cellsize
!         %alpc          real       angle of the grid
!         %lunit         char       likely unit (meter, degrees)
!
!       recordaxe         ttype      struct with fields
!         %content       int.a      seconds since 1-jan-1970 / run
!         %delta         int        timestep (s) / 1
!         %ncontent      int        number of times / runs
!         %nstatm        bool       non-stationary mode (time /run switch)
!         %dimid         int        netcdf dimension id
!         %varid         int        netcdf variable id
!
!       spcgrid          spcgtype   struct with fields
!         %points        int        number of spectra (points dim)
!         %frequency     r.a        frequency axe (dim)
!         %direction     r.a.       directional axe (dim)
!         %longitude     r.a        x coordinates points
!         %latitude      r.a        y coordinates points
!         %lunit         char       likely unit (meter, degrees)
!
!       unstrgrid        unstrgtype struct with fields
!         %npoints        int        number of spectra (points dim)
!         %longitude     r.a        x coordinates points
!         %latitude      r.a        y coordinates points
!         %lunit         char       likely unit (meter, degrees)
!
!       nctable          ncttype    struct array from nctablemd.f90
!         %standard_name char       see remarks a
!         %long_name     char       "
!         %units         char       "
!         %nctype        integer    agnc_short, agnc_byte etc.
!         %datamin       real       minimum expected data value
!         %datamax       real       maximum expected data value
!         %varid         int        nc90_var_id (remark e)
!
!       fill_values      fvtype     struct with fields
!         %byte          int*1
!         %short         int*2
!         %float         real
!         %double        real*8
!
!  3.a.variables
!       ncid             int        identifier of opened / created file
!       pi               real       3.1415927
!       r2d              real       180. / pi
!       d2r              real       pi / 180.
!       e_code           int        error code (0: no error, 1: warning >1: exits)
!       e_msg            char       error string
!
!  4.procedures
!       ---- create new netcdf files ----
!       create_ncfile(ncfile, ncid)
!       agnc_set_mapgrid(ncid, mapgrid)
!       agnc_set_spcgrid(ncid, spcgrid)
!       agnc_set_time(ncid, recordaxe)
!       set_monthly(reftime, recordaxe)
!       agnc_define_mapvariable(ncid, name)
!       agnc_define_pntvariables
!       ------------- or ---------------
!       open_ncfile( ncfile, 'read/write' )
!       --------- get data -------------
!       agnc_get_mapgrid(ncid, mapgrid )
!       agnc_get_spcgrid(ncid, spcgrid )
!       agnc_get_recordaxe(ncid, ttype )
!       agnc_get_var(ncid, standard_name, result, <noautoscale>)
!       <agnc_get_mapvar(ncid, standard_name, result, lon=lon, lat=lat, time=time)>
!       <agnc_get_spcvar(time, points, result)>
!       ------ insert/append data ------
!       agnc_add_mapdata(standard_name, time, data)
!       agnc_add_spcvariables(time, e1, th1, sp1 )
!       ------------- utils -------------
!       agnc_get_varid_by_name(ncid, searchfor, varid)
!       close_ncfile(ncid)
!       nccheck(status)
!       axe_is_regular(x, nx, dx)
!       derive_mapgrid(longitude, latitude, mapgrid)
!       derive_mapgrid(longitude, latitude, points, frequency, direction, spcgrid)
!       find_time_index(ttype, time)
!       compute_1d_spectra(energy, e1, th1, sp1, undef)
!
!  5. used by :
!
!     ww3_outnc.ftn (wavewatch)
!     swoutp (swan)
!
!  6. exit codes
!     30) Netcdf F90 library error
!     31) Uknown nc format in nctable
!     32) Error while reading grid from file
!
!  7. remarks :
!     a. http://cf-pcmdi.llnl.gov/documents/cf-standard-names/
!     b. int16 = round((data - offset) / scale)
!        scale = (datamax - datamin) / (2^n -2)
!        offset = (datamax - datamin)/2
!     c. officially, netcdf apps should un-/pack data on request
!        however, i didn't feel like writing the interface on
!        http://www.lahey.com/lookat90.htm, interface to deal with
!        the different data storage types.
!     d. this module is fortran 95 due to the usage of allocatable
!        arrays in derived data type (eg, structs)
!     e. set when opening an existing file or creation of a variable in
!        a new file
!
      module agioncmd
        use netcdf
        use nctablemd
        implicit none
        real, parameter                              :: AGIONCMD_VERSION = 1.2
        real, parameter                              :: AGNC_DUMMY = NF90_FILL_FLOAT
        integer(kind=1), parameter                   :: AGNC_FILL_BYTE = -2**7
        integer(kind=2), parameter                   :: AGNC_FILL_SHORT= -2**15
        integer(kind=1)                              :: e_code
        character(1024)                              :: e_msg
        type mapgrid_type
            real, dimension(:), allocatable          :: longitude, latitude
            real, dimension(:,:), allocatable        :: xmesh, ymesh
            integer                                  :: nx, ny, lon_dimid, lon_varid, &
                                                        lat_dimid, lat_varid
            real                                     :: dx, dy, angle = 0.
            character(len=10)                        :: lunit = 'degrees'
        end type mapgrid_type

        type mmapgrid_type
            real, dimension(:,:), allocatable        :: longitude, latitude
            integer                                  :: nx, ny, lon_dimid, lon_varid, &
                                                        lat_dimid, lat_varid
            real                                     :: alpc = 0.
            character(len=10)                        :: lunit = 'degrees'
        end type mmapgrid_type

        type unstrgrid_type
            real, dimension(:), allocatable          :: longitude, latitude
            integer                                  :: npoints, pnt_dimid, lon_varid, &
                                                        lat_varid
            character(len=10)                        :: lunit = 'degrees'
        end type unstrgrid_type

        type spcgrid_type
            real, dimension(:), allocatable          :: frequency, direction, &
                                                        longitude, &
                                                        latitude
            integer                                  :: frq_dimid, frq_varid, &
                                                        dir_dimid, dir_varid, pnt_dimid, &
                                                        lon_varid, lat_varid, &
                                                        nfreq, ndir, npoints
            logical                                  :: relative
            character(len=10)                        :: lunit = 'degrees'
            character(len=10)                        :: convention = 'nautical'
        end type spcgrid_type

        type recordaxe_type
            integer(kind=4), dimension(:), allocatable  :: content
            integer(kind=4)                             :: delta
            integer                                     :: ncontent, dimid, varid
            logical                                     :: nstatm = .true.
        end type recordaxe_type

        ! Given the dimension size and/or type, select proper procedure/function
        interface agnc_get_var
            module procedure agnc_get_var1d_integer, agnc_get_var2d_integer, &
                             agnc_get_var3d_integer, agnc_get_var4d_integer, &
                             agnc_get_var1d_float,   agnc_get_var1d_double, &
                             agnc_get_var2d_float,   agnc_get_var2d_double, &
                             agnc_get_var3d_float,   agnc_get_var3d_double, &
                             agnc_get_var4d_float,   agnc_get_var4d_double
        end interface agnc_get_var
        interface agnc_get_spcvars
            module procedure agnc_get_spcvars_vector, agnc_get_spcvars_scalar
        end interface agnc_get_spcvars

        interface axe_is_regular
            module procedure axe_is_regular_integer, &
                             axe_is_regular_float, &
                             axe_is_regular_float_2d, &
                             axe_is_regular_double
        end interface axe_is_regular
        interface agnc_add_mapdata
            module procedure agnc_add_mapdata2d_float, agnc_add_mapdata3d_float, &
                             agnc_add_mapdata2d_double, agnc_add_mapdata3d_double
        end interface agnc_add_mapdata
        interface agnc_add_spcdata
            module procedure agnc_add_spcdata_3d, agnc_add_spcdata_2d
        end interface agnc_add_spcdata
        interface agnc_add_pntdata
            module procedure agnc_add_pntdata1d_float, agnc_add_pntdata2d_float, &
                             agnc_add_pntdata3d_float, &
                             agnc_add_pntdata1d_double, agnc_add_pntdata2d_double, &
                             agnc_add_pntdata3d_double
        end interface agnc_add_pntdata
        interface seconds_since_epoch
            module procedure seconds_since_epoch_datevec, seconds_since_epoch_twoint
        end interface seconds_since_epoch
        interface datevec
            module procedure datevec_from_twoint, datevec_from_epoch
        end interface datevec
      contains
!
! ---- create new netcdf files ----
!       name                                            status
!       create_ncfile(ncfile, ncid)                     done
!       agnc_define_mapgrid(ncid, mapgrid)              done
!       agnc_define_spcgrid(ncid, spcgrid, Escale)      done
!       agnc_define_unstrgrid(ncid, unstrgrid)            done
!       agnc_define_mapvariable(ncid, name)             done
!       agnc_define_pntvariables                        done(moved to define_spcgrid)
!
        subroutine create_ncfile( ncfile, ncid, recordaxe, mapgrid, mmapgrid, spcgrid, unstrgrid, Escale, meteorological )
            character(len=*),     intent(in)                    :: ncfile
            integer,              intent(out)                   :: ncid
            type(recordaxe_type), intent(inout)                 :: recordaxe
            type(mapgrid_type),   intent(inout), optional       :: mapgrid
            type(mmapgrid_type),  intent(inout), optional       :: mmapgrid
            type(spcgrid_type),   intent(inout), optional       :: spcgrid
            type(unstrgrid_type) ,intent(inout), optional       :: unstrgrid
            logical              ,intent(   in), optional       :: Escale
            logical              ,intent(   in), optional       :: meteorological
            character(100)                                      :: history, dconv = 'cartesian'
            logical                                             :: do_scale

            if ( present(meteorological) ) then
                if ( meteorological ) dconv = 'meteorological'
            end if

            do_scale = .true.
            if ( present(Escale) ) then
                do_scale = Escale
            end if

            write(history,'(A,F3.1)'), 'Created with agioncmd version ', AGIONCMD_VERSION
            call nccheck ( nf90_create( ncfile, NF90_NOCLOBBER + NF90_64BIT_OFFSET, ncid) )
            if ( recordaxe%nstatm ) then
                call nccheck ( nf90_def_dim( ncid, 'time', NF90_UNLIMITED, recordaxe%dimid ) );
                call nccheck ( nf90_def_var( ncid, 'time', NF90_INT, recordaxe%dimid, recordaxe%varid ) );
                call nccheck ( nf90_put_att( ncid, recordaxe%varid, 'units', 'seconds since 1970-01-01') )
                call nccheck ( nf90_put_att( ncid, recordaxe%varid, 'calendar', 'julian') )
                call nccheck ( nf90_put_att( ncid, recordaxe%varid, 'long_name', 'time') )
            else
                call nccheck ( nf90_def_dim( ncid, 'run', NF90_UNLIMITED, recordaxe%dimid ) );
                call nccheck ( nf90_def_var( ncid, 'run', NF90_INT, recordaxe%dimid, recordaxe%varid ) );
                call nccheck ( nf90_put_att( ncid, recordaxe%varid, 'long_name', 'run number') )
                call nccheck ( nf90_put_att( ncid, recordaxe%varid, 'units', '1' ) )
            end if

            ! global attributes
            call nccheck ( nf90_put_att( ncid, NF90_GLOBAL, 'Conventions', 'CF-1.5') )
            call nccheck ( nf90_put_att( ncid, NF90_GLOBAL, 'History', history) )
            call nccheck ( nf90_put_att( ncid, NF90_GLOBAL, 'Directional convention', dconv) )

            if ( present(mapgrid) ) then
                call agnc_define_mapgrid(ncid, mapgrid)
                ! Restore convention for maps
                call nccheck ( nf90_put_att( ncid, NF90_GLOBAL, 'Conventions', 'CF-1.5') )
            end if

            if ( present(mmapgrid) )  call agnc_define_mmapgrid(ncid, mmapgrid)
            if ( present(spcgrid) )   call agnc_define_spcgrid(ncid, spcgrid, do_scale)
            if ( present(unstrgrid) ) call agnc_define_unstrgrid(ncid, unstrgrid)

            ! reserve extra space in the header for adding attributes
            ! later without reorganizing the file structure. 4kb should be
            ! enough for most attributes and negligible compared to the
            ! total file size.
            call nccheck( nf90_enddef(ncid, h_minfree=4096) )
            ! Set the file back to definition mode so variables can be defined later on
            call nccheck( nf90_redef(ncid) )

        end subroutine create_ncfile

        subroutine coordinate_names_units( lunit, xname, xunit, yname, yunit )
            character(len=10),         intent(   in):: lunit
            character(len=20),         intent(  out):: xname, xunit
            character(len=20),         intent(  out):: yname, yunit
             if ( trim(lunit) == 'degrees') then
                xname = 'longitude'
                xunit = 'degrees_east'
                yname = 'latitude'
                yunit = 'degrees_north'
            else
                xname = 'x'
                xunit = lunit
                yname = 'y'
                yunit = lunit
            end if
        end subroutine coordinate_names_units

        subroutine agnc_define_mapgrid( ncid, mapgrid )
            integer,            intent( in)         :: ncid
            type (mapgrid_type),intent( inout)      :: mapgrid
            character(len=20)                       :: xname, xunit
            character(len=20)                       :: yname, yunit
            call coordinate_names_units( mapgrid%lunit, xname, xunit, yname, yunit )
            call nccheck ( nf90_def_dim( ncid, xname, mapgrid%nx, mapgrid%lon_dimid ) );
            call nccheck ( nf90_def_var( ncid, xname, NF90_FLOAT, mapgrid%lon_dimid, mapgrid%lon_varid ) );
            call nccheck ( nf90_put_att( ncid, mapgrid%lon_varid, 'units', xunit) )
            call nccheck ( nf90_put_att( ncid, mapgrid%lon_varid, 'long_name', xname) )

            call nccheck ( nf90_def_dim( ncid, yname, mapgrid%ny, mapgrid%lat_dimid ) );
            call nccheck ( nf90_def_var( ncid, yname, NF90_FLOAT, mapgrid%lat_dimid, mapgrid%lat_varid ) );
            call nccheck ( nf90_put_att( ncid, mapgrid%lat_varid, 'units', yunit) )
            call nccheck ( nf90_put_att( ncid, mapgrid%lat_varid, 'long_name', yname) )
        end subroutine agnc_define_mapgrid

        subroutine agnc_define_mmapgrid( ncid, mapgrid )
            integer,             intent( in)        :: ncid
            type (mmapgrid_type),intent( inout)     :: mapgrid
            character(len=20)                       :: xname, xunit
            character(len=20)                       :: yname, yunit
            call coordinate_names_units( mapgrid%lunit, xname, xunit, yname, yunit )

            call nccheck ( nf90_def_dim( ncid, 'xc', mapgrid%nx, mapgrid%lon_dimid ) );
            call nccheck ( nf90_def_dim( ncid, 'yc', mapgrid%ny, mapgrid%lat_dimid ) );

            call nccheck ( nf90_def_var( ncid, xname, NF90_FLOAT, &
                                         (/ mapgrid%lon_dimid, mapgrid%lat_dimid /), &
                                         mapgrid%lon_varid ) )
            call nccheck ( nf90_def_var( ncid, yname, NF90_FLOAT, &
                                         (/ mapgrid%lon_dimid, mapgrid%lat_dimid /), &
                                         mapgrid%lat_varid ) )

            call nccheck ( nf90_put_att( ncid, mapgrid%lon_varid, 'units', xunit) )
            call nccheck ( nf90_put_att( ncid, mapgrid%lon_varid, 'long_name', xname) )
            call nccheck ( nf90_put_att( ncid, mapgrid%lat_varid, 'units', yunit) )
            call nccheck ( nf90_put_att( ncid, mapgrid%lat_varid, 'long_name', yname) )
        end subroutine agnc_define_mmapgrid

        subroutine agnc_define_unstrgrid( ncid, unstrgrid )
            integer,            intent( in)         :: ncid
            type (unstrgrid_type),intent( inout)    :: unstrgrid
            character(len=20)                       :: xname, xunit
            character(len=20)                       :: yname, yunit
            call coordinate_names_units( unstrgrid%lunit, xname, xunit, yname, yunit )

            call nccheck ( nf90_def_dim( ncid, 'points', unstrgrid%npoints, unstrgrid%pnt_dimid ) )

            call nccheck ( nf90_def_var( ncid, xname, NF90_FLOAT, unstrgrid%pnt_dimid, unstrgrid%lon_varid ) );
            call nccheck ( nf90_put_att( ncid, unstrgrid%lon_varid, 'units', xunit) )
            call nccheck ( nf90_put_att( ncid, unstrgrid%lon_varid, 'long_name', xname) )

            call nccheck ( nf90_def_var( ncid, yname, NF90_FLOAT, unstrgrid%pnt_dimid, unstrgrid%lat_varid ) );
            call nccheck ( nf90_put_att( ncid, unstrgrid%lat_varid, 'units', yunit) )
            call nccheck ( nf90_put_att( ncid, unstrgrid%lat_varid, 'long_name', yname) )

        end subroutine agnc_define_unstrgrid

        subroutine agnc_define_spcgrid( ncid, spcgrid, Escale )
            integer,            intent( in)         :: ncid
            logical, optional,  intent( in)         :: Escale
            type (spcgrid_type),intent( inout)      :: spcgrid
            integer                                 :: evarid, ra_dimid, dvarid
            type(nctable_record)                    :: trecord
            logical                                 :: do_scale
            character(len=20)                       :: xname, xunit
            character(len=20)                       :: yname, yunit

            call coordinate_names_units( spcgrid%lunit, xname, xunit, yname, yunit )

            do_scale = .true.
            if ( present(Escale) ) do_scale = Escale

            !
            ! frequency
            !
            call nccheck ( nf90_def_dim( ncid, 'frequency', spcgrid%nfreq, spcgrid%frq_dimid ) );
            call nccheck ( nf90_def_var( ncid, 'frequency', NF90_FLOAT, spcgrid%frq_dimid, spcgrid%frq_varid ) );
            call nccheck ( nf90_put_att( ncid, spcgrid%frq_varid, 'units', 's-1') )
            call nccheck ( nf90_put_att( ncid, spcgrid%frq_varid, 'standard_name', 'wave_frequency') )

            !
            ! points dimension
            !
            call nccheck ( nf90_def_dim( ncid, 'points', spcgrid%npoints, spcgrid%pnt_dimid ) )

            !
            ! longitude
            !
            call nccheck ( nf90_def_var( ncid, xname, NF90_FLOAT, spcgrid%pnt_dimid, spcgrid%lon_varid ) );
            call nccheck ( nf90_put_att( ncid, spcgrid%lon_varid, 'units', xunit) )
            call nccheck ( nf90_put_att( ncid, spcgrid%lon_varid, 'long_name', xname//' of spectral point') )

            !
            ! latitude
            !
            call nccheck ( nf90_def_var( ncid, yname, NF90_FLOAT, spcgrid%pnt_dimid, spcgrid%lat_varid ) );
            call nccheck ( nf90_put_att( ncid, spcgrid%lat_varid, 'units', yunit) )
            call nccheck ( nf90_put_att( ncid, spcgrid%lat_varid, 'long_name', yname//' of spectral point') )

            call agnc_get_recordaxe_ids( ncid, ra_dimid )

            if ( spcgrid%ndir > 0 ) then
                !
                ! 2D direction axe
                !
                call nccheck ( nf90_def_dim( ncid, 'direction', spcgrid%ndir, spcgrid%dir_dimid ) )
                call nccheck ( nf90_def_var( ncid, 'direction', NF90_FLOAT, spcgrid%dir_dimid, spcgrid%dir_varid ) );
                call nccheck ( nf90_put_att( ncid, spcgrid%dir_varid, 'units', 'radians') )
                call nccheck ( nf90_put_att( ncid, spcgrid%dir_varid, 'long_name', 'direction') )
                call nccheck ( nf90_put_att( ncid, spcgrid%dir_varid, 'standard_name', 'sea_surface_wave_to_direction') )

                !
                ! 2D density spectrum
                !
                if ( do_scale) then
                    call nccheck ( nf90_def_var( ncid, 'density', NF90_SHORT, &
                                            (/ spcgrid%frq_dimid, spcgrid%dir_dimid, spcgrid%pnt_dimid, ra_dimid /), dvarid ) )
                    call nccheck ( nf90_put_att( ncid, dvarid, 'long_name', &
                            'Density. Note that you have to multiply with scale_density in order to get final values') )
                    call nccheck ( nf90_put_att( ncid, dvarid, '_FillValue', NF90_FILL_SHORT) )

                else
                    call nccheck ( nf90_def_var( ncid, 'density', NF90_FLOAT, &
                                            (/ spcgrid%frq_dimid, spcgrid%dir_dimid, spcgrid%pnt_dimid, ra_dimid /), dvarid ) )
                    call nccheck ( nf90_put_att( ncid, dvarid, 'long_name', &
                            'Density') )
                    call nccheck ( nf90_put_att( ncid, dvarid, '_FillValue', NF90_FILL_FLOAT) )
                end if
                call nccheck ( nf90_put_att( ncid, dvarid, 'units', 'm2 s rad-1') )
                call nccheck ( nf90_put_att( ncid, dvarid, 'standard_name', &
                            'sea_surface_wave_directional_variance_spectral_density') )

                if ( spcgrid%relative ) &
                     call nccheck ( nf90_put_att( ncid, dvarid, 'relative_to_current', 'true') )

                !
                ! scale_density
                !
                if ( do_scale ) then
                    call nccheck ( nf90_def_var( ncid, 'scale_density', NF90_FLOAT, &
                                            (/ spcgrid%pnt_dimid, ra_dimid /), dvarid ) )
                    call nccheck ( nf90_put_att( ncid, dvarid, 'long_name', &
                                        'Multiply density values with this scale factor') )
                    call nccheck ( nf90_put_att( ncid, dvarid, 'units', '1') )
                end if
            else
                if ( do_scale) then
                    call nccheck ( nf90_def_var( ncid, 'energy_1d', NF90_SHORT, &
                                            (/ spcgrid%frq_dimid, spcgrid%pnt_dimid,ra_dimid /), evarid ) )
                    call nccheck ( nf90_put_att( ncid, evarid, '_FillValue', NF90_FILL_SHORT) )
                else
                    call nccheck ( nf90_def_var( ncid, 'energy_1d', NF90_FLOAT, &
                                            (/ spcgrid%frq_dimid, spcgrid%pnt_dimid,ra_dimid /), evarid ) )
                    call nccheck ( nf90_put_att( ncid, evarid, '_FillValue', NF90_FILL_FLOAT) )
                end if
                call nccheck ( nf90_put_att( ncid, evarid, 'units', 'm2 s') )
                call nccheck ( nf90_put_att( ncid, evarid, 'standard_name', 'sea_surface_wave_variance_spectral_density') )

                if ( spcgrid%relative ) &
                     call nccheck ( nf90_put_att( ncid, dvarid, 'relative_to_current', 'true') )

                !
                ! scale_density
                !
                if ( do_scale ) then
                    call nccheck ( nf90_put_att( ncid, evarid, 'long_name', &
                            '1D energy. Note that you have to multiply with scale_energy_1d in order to get final values') )

                    call nccheck ( nf90_def_var( ncid, 'scale_energy_1d', NF90_FLOAT, &
                                            (/ spcgrid%pnt_dimid, ra_dimid /), evarid ) )
                    call nccheck ( nf90_put_att( ncid, evarid, 'long_name', &
                            'Multiply energy_1d values with this scale factor') )
                    call nccheck ( nf90_put_att( ncid, evarid, 'units', &
                            '1') )
                else
                    call nccheck ( nf90_put_att( ncid, evarid, 'long_name', &
                            '1D energy') )
                end if

                !
                ! theta_1d
                !
                call nccheck ( nf90_def_var( ncid, 'theta_1d', NF90_BYTE, &
                                            (/ spcgrid%frq_dimid, spcgrid%pnt_dimid,ra_dimid /), evarid ) )
                call nccheck ( nf90_put_att( ncid, evarid, 'units', 'degree') )
                call nccheck ( nf90_put_att( ncid, evarid, 'long_name', 'principal wave direction') )
                call nccheck ( nf90_put_att( ncid, evarid, '_FillValue', AGNC_FILL_BYTE) )
                call nccheck ( nf90_put_att( ncid, evarid, 'scale_factor', 360. / (2**8 -2)) )
                call nccheck ( nf90_put_att( ncid, evarid, 'add_offset', 360. / 2.) )

                !
                ! spread_1d
                !
                call nccheck ( nf90_def_var( ncid, 'spread_1d', NF90_BYTE, &
                                            (/ spcgrid%frq_dimid, spcgrid%pnt_dimid,ra_dimid /), evarid ) )
                call nccheck ( nf90_put_att( ncid, evarid, 'units', 'degree') )
                call nccheck ( nf90_put_att( ncid, evarid, 'long_name', &
                    'Longuet-Higgins short-crestedness parameter (s in cos(theta/2)^2s)') )
                call nccheck ( nf90_put_att( ncid, evarid, '_FillValue', AGNC_FILL_BYTE) )

            end if

        end subroutine agnc_define_spcgrid

        subroutine agnc_define_mapvariable (ncid, varname)
            integer,                intent( in)         :: ncid
            character(len=*),       intent( in)         :: varname
            type(nctable_record)                        :: trecord
            integer                                     :: ra_dimid, lon_dimid, lat_dimid,&
                                                           varid

            call get_nctable_record(varname, trecord)

            call agnc_get_recordaxe_ids(ncid, ra_dimid)
            call agnc_get_mapgrid_ids(ncid, lon_dimid, lat_dimid)

            call nccheck ( nf90_def_var( ncid, trecord%name , trecord%nctype, &
                                            (/ lon_dimid, lat_dimid, ra_dimid /), varid ) )

            call agcn_add_variable_attributes( ncid, varid, trecord )

        end subroutine agnc_define_mapvariable

        subroutine agnc_define_pntvariable (ncid, varname)
            integer,                intent( in)         :: ncid
            character(len=*),       intent( in)         :: varname
            type(nctable_record)                        :: trecord
            integer                                     :: pnt_dimid, ra_dimid, varid

            call get_nctable_record(varname, trecord)

            call agnc_get_recordaxe_ids(ncid, ra_dimid)

            ! both unstructured as spectral grid have the "point" dimension
            call agnc_get_spcgrid_ids(ncid, pnt_dimid=pnt_dimid)

            call nccheck ( nf90_def_var( ncid, trecord%name , trecord%nctype, &
                                            (/ pnt_dimid, ra_dimid /), varid ) )

            call agcn_add_variable_attributes( ncid, varid, trecord )

        end subroutine agnc_define_pntvariable

        subroutine agcn_add_variable_attributes( ncid, varid, trecord )
            integer,                intent(in) :: ncid, varid
            type(nctable_record),   intent(in) :: trecord
            real(kind=4)                       :: add_offset, scale_factor, rng
            integer                            :: tmpid
            character(len=nf90_max_name)       :: lonname

            call nccheck ( nf90_put_att( ncid, varid, 'units', trecord%units) )
            if ( trim(trecord%standard_name) /= 'none' ) &
                call nccheck ( nf90_put_att( ncid, varid, 'standard_name', trecord%standard_name) )

            call nccheck ( nf90_put_att( ncid, varid, 'long_name', trecord%long_name) )

            ! CF-1.2+ prescibes that for a multi-dimension coordinate field a coordinates atribute

            if ( nf90_inq_dimid(ncid, 'xc', tmpid) == nf90_noerr ) then
                ! multi-dimension field because of the hard coded dimension name. The coordinates
                ! are either in the longitude/latitude or x,y variables
                if ( nf90_inq_varid(ncid, 'x', tmpid) == nf90_noerr ) then
                    call nccheck ( nf90_put_att( ncid, varid, 'coordinates', 'x y') )
                else
                    call nccheck ( nf90_put_att( ncid, varid, 'coordinates', 'longitude latitude') )
                end if
            end if

            ! http://www.unidata.ucar.edu/software/netcdf/docs/BestPractices.html
            !
            ! In either the signed or unsigned case, an alternate formula may be used for the add_offset
            ! and scale_factor packing parameters that reserves a packed value for a special value, such
            ! as an indicator of missing data. For example, to reserve the minimum packed value (-2n - 1)
            ! for use as a special value in the case of signed packed values:
            !
            !    scale_factor =(dataMax - dataMin) / (2n - 2)
            !    add_offset = (dataMax + dataMin) / 2
            !
            ! Note that the NF90_FILL_ values are one off in comparison with the formulea used here
            !
            ! if datamin equals datamax in an integer kind of way, do not apply scaling
            rng = (trecord%datamax - trecord%datamin)
            add_offset = (trecord%datamin + trecord%datamax) * .5
            if ( rng /= 0 ) then
                select case ( trecord%nctype )
                case (NF90_BYTE)
                    call nccheck ( nf90_put_att( ncid, varid, '_FillValue', AGNC_FILL_BYTE) )
                    scale_factor = rng/(2**8-2)
                case (NF90_SHORT)
                    call nccheck ( nf90_put_att( ncid, varid, '_FillValue', AGNC_FILL_SHORT) )
                    scale_factor = rng/(2**16-2)
                case default
                    call agnc_set_error( 'The nctable contains an invalid agnc_TYPE')
                end select
                call nccheck ( nf90_put_att( ncid, varid, 'scale_factor', scale_factor) )
                call nccheck ( nf90_put_att( ncid, varid, 'add_offset', add_offset) )
            end if
            select case ( trecord%nctype )
            case (NF90_BYTE)
                call nccheck ( nf90_put_att( ncid, varid, '_FillValue', AGNC_FILL_BYTE) )
            case (NF90_SHORT)
                call nccheck ( nf90_put_att( ncid, varid, '_FillValue', AGNC_FILL_SHORT) )
            end select

        end subroutine agcn_add_variable_attributes

!       ------------- or ---------------
!       name                                            status
!       open_ncfile( ncfile, 'read/write', &
!                    recordaxe, monthly )       done
!
        subroutine open_ncfile( ncfile, permission, ncid, recordaxe, monthly )
        !
        ! recordaxe, is optional. When permission is "read", the information is read from the file.
        ! When permission is write, a check is made whether the recordaxe and grid is set in
        ! the dimension variables in the file.
            character(len=*),                intent(in)      :: ncfile, permission
            integer         ,                intent(out)     :: ncid
            type(recordaxe_type) , optional, intent(inout)   :: recordaxe
            logical              , optional, intent(   in)   :: monthly
            type(recordaxe_type)                             :: recordaxe_tmp
            if ( permission == "read") then
                call nccheck( nf90_open(ncfile, nf90_nowrite, ncid)  )
            else
                call nccheck( nf90_open(ncfile, nf90_write, ncid)  )
            end if

            if ( present(recordaxe)) then
                call agnc_get_recordaxe(ncid, recordaxe_tmp)
                ! Check if read time axe is regular and does not only contain _FillValues
                if ( axe_is_regular(recordaxe_tmp%content, recordaxe_tmp%ncontent, recordaxe_tmp%delta, .true.) ) then
                    recordaxe = recordaxe_tmp
                else
                    if ( axe_is_regular(recordaxe%content, recordaxe%ncontent, recordaxe%delta, .true.)) then
                        if ( present(monthly) ) then
                            call agnc_set_recordaxe(ncid, recordaxe, monthly)
                        else
                            call agnc_set_recordaxe(ncid, recordaxe)
                        end if
                    else
                        ! If you're writing a program, you end up here if you didn't provide a
                        ! regular grid. Check that dt and nt are set
                        call agnc_set_error( 'Could not read valid record axe from file')
                    end if
                end if
            end if
        end subroutine open_ncfile

!
!       --------- set (quasi-) dimension variables -------------
!       name                                            status
!       agnc_set_recordaxe (ncid, recordaxe, month)       done
!       agnc_set_mapgrid (ncid, mapgrid)                  done
!       agnc_set_mmapgrid (ncid, mmapgrid)                done
!       agnc_set_unstrgrid (ncid, unstrgrid)              done
!       agnc_set_spcgrid (ncid, spcgrid)                  done

        subroutine agnc_set_recordaxe (ncid, recordaxe, monthly)
        ! When the optional monthly is set to true, the time axis in the netcdf file
        ! is predeclared for the month of trecordaxe%content(1)
            integer,              intent (   in)                  :: ncid
            type (recordaxe_type), intent (inout)                 :: recordaxe
            logical, optional   , intent (   in)                  :: monthly
            character(len=nf90_max_name)                          :: name
            integer                                               :: t1(6), t2(6), i
            if ( recordaxe%varid < 1 ) call agnc_get_recordaxe_ids(ncid, recordaxe%dimid, recordaxe%varid)

            ! Some sanity checks
            if (recordaxe%ncontent == 0) then
                call agnc_set_error( 'Cannot set empty record axis')
            end if
            if ( present(monthly) ) then
                if (monthly) then
                    t1 = datevec(recordaxe%content(1))
                    t1(3) = 1
                    t1(4:6) = (/0,0,0/)
                    t2=t1
                    t2(2) = t2(2)+1
                    recordaxe%ncontent = (seconds_since_epoch(t2) - seconds_since_epoch(t1)) / recordaxe%delta
                    deallocate(recordaxe%content)
                    allocate(recordaxe%content(recordaxe%ncontent))
                    recordaxe%content = (/ (i, i=seconds_since_epoch(t1), seconds_since_epoch(t2) - recordaxe%delta, recordaxe%delta) /)
                end if
            end if
            if (.not. axe_is_regular(recordaxe%content, recordaxe%ncontent, recordaxe%delta) ) then
                call close_ncfile(ncid)
                call agnc_set_error( 'Failed to append data to this file. It would result in a irregular record axis');
            end if

            call nccheck ( nf90_put_var(ncid, recordaxe%varid, recordaxe%content) )

        end subroutine agnc_set_recordaxe

        subroutine agnc_set_unstrgrid (ncid, unstrgrid)
            integer,                intent (   in)                  :: ncid
            type (unstrgrid_type),  intent (inout)                  :: unstrgrid
            character(len=nf90_max_name)                            :: name
            if ( unstrgrid%lon_varid < 1 ) then
                call agnc_get_unstrgrid_ids(ncid, unstrgrid%pnt_dimid,        &
                                    unstrgrid%lon_varid, unstrgrid%lat_varid);
            end if
            call nccheck ( nf90_put_var(ncid, unstrgrid%lon_varid, unstrgrid%longitude) )
            call nccheck ( nf90_put_var(ncid, unstrgrid%lat_varid, unstrgrid%latitude) )
        end subroutine agnc_set_unstrgrid

        subroutine agnc_set_mapgrid (ncid, mapgrid)
            integer,              intent (   in)                  :: ncid
            type (mapgrid_type),  intent (inout)                  :: mapgrid
            character(len=nf90_max_name)                          :: name
            if ( mapgrid%lon_varid < 1 ) then
                call agnc_get_mapgrid_ids(ncid, mapgrid%lon_dimid, mapgrid%lat_dimid, &
                                              mapgrid%lon_varid, mapgrid%lat_varid);
            end if
            call nccheck ( nf90_put_var(ncid, mapgrid%lon_varid, mapgrid%longitude) )
            call nccheck ( nf90_put_var(ncid, mapgrid%lat_varid, mapgrid%latitude) )
        end subroutine agnc_set_mapgrid

        subroutine agnc_set_mmapgrid (ncid, mapgrid)
            integer,              intent (   in)                  :: ncid
            type (mmapgrid_type), intent (inout)                  :: mapgrid
            character(len=nf90_max_name)                          :: name
            if ( mapgrid%lon_varid < 1 ) then
                call agnc_get_mapgrid_ids(ncid, mapgrid%lon_dimid, mapgrid%lat_dimid, &
                                              mapgrid%lon_varid, mapgrid%lat_varid);
            end if
            call nccheck ( nf90_put_var(ncid, mapgrid%lon_varid, mapgrid%longitude) )
            call nccheck ( nf90_put_var(ncid, mapgrid%lat_varid, mapgrid%latitude) )
        end subroutine agnc_set_mmapgrid


        subroutine agnc_set_spcgrid (ncid, spcgrid)
            integer,              intent (   in)                  :: ncid
            type (spcgrid_type),  intent (inout)                  :: spcgrid
            character(len=nf90_max_name)                          :: name
            if ( spcgrid%frq_varid < 1 ) then
                call agnc_get_spcgrid_ids(ncid, spcgrid%frq_dimid, spcgrid%frq_varid, &
                                              spcgrid%pnt_dimid, &
                                              spcgrid%dir_dimid, spcgrid%dir_varid, &
                                              spcgrid%lon_varid, spcgrid%lat_varid)
            end if

            call nccheck ( nf90_put_var(ncid, spcgrid%frq_varid, spcgrid%frequency) )
            call nccheck ( nf90_put_var(ncid, spcgrid%lon_varid, spcgrid%longitude) )
            call nccheck ( nf90_put_var(ncid, spcgrid%lat_varid, spcgrid%latitude) )
            if (spcgrid%ndir > 0) &
                call nccheck ( nf90_put_var(ncid, spcgrid%dir_varid, spcgrid%direction) )

        end subroutine agnc_set_spcgrid

!       --------- get data -------------
!       name                                                        status
!       agnc_get_mapgrid(ncid, mapgrid )                            done
!       agnc_get_spcgrid(ncid, spcgrid )                            done
!       agnc_get_recordaxe(ncid, recordaxe )                          done
!       agnc_get_mapvar(ncid, name, values, noautoscale)            done
!       agnc_get_spcvars(ncid, times, points, energy, theta, spr)   done
!

        subroutine agnc_get_mapgrid(ncid, mapgrid)
            ! Structure:
            !   find dimension id of dimenstion named longitude, lon or x
            !   get meta information on this dimension
            !   assume the name of the dimension corresponds with a variable
            !   read variable
            !   set results in mapgrid struct
            !   check that geographic axes are regular
            integer, intent ( in )              :: ncid
            type (mapgrid_type), intent (out)   :: mapgrid

            call agnc_get_mapgrid_ids(ncid, mapgrid%lon_dimid, mapgrid%lat_dimid, &
                                          mapgrid%lon_varid, mapgrid%lat_varid, &
                                          mapgrid%nx, mapgrid%ny)

!             if ( size(shape(mapgrid%longitude)) == 1) then
                allocate(mapgrid%longitude(mapgrid%nx))
                call nccheck ( nf90_get_var(ncid, mapgrid%lon_varid, mapgrid%longitude) );
                mapgrid%dx = mapgrid%longitude(2) - mapgrid%longitude(1)
                if (.not. axe_is_regular(mapgrid%longitude, mapgrid%nx, mapgrid%dx) ) then
                    call agnc_set_error( 'Longitude axe does not appear to be regular')
                end if
!             else
!                 allocate(mapgrid%longitude(mapgrid%ny, mapgrid%nx))
!                 call nccheck ( nf90_get_var(ncid, mapgrid%lon_varid, mapgrid%longitude) );
!             end if

!             if ( size((shape(mapgrid%latitude)) == 1) then
                allocate(mapgrid%latitude(mapgrid%ny))
                call nccheck ( nf90_get_var(ncid, mapgrid%lat_varid, mapgrid%latitude) )
                mapgrid%dy = mapgrid%latitude(2) - mapgrid%latitude(1)
                if (.not. axe_is_regular(mapgrid%latitude, mapgrid%ny, mapgrid%dy) ) then
                    call agnc_set_error( 'Latitude axe does not appear to be regular')
                end if
!             else
!                 allocate(mapgrid%latitude(mapgrid%ny, mapgrid%nx))
!                 call nccheck ( nf90_get_var(ncid, mapgrid%lat_varid, mapgrid%latitude) )
!             end if
        end subroutine agnc_get_mapgrid

        subroutine agnc_get_unstrgrid(ncid, unstrgrid)
            ! Structure:
            !   find dimension id of dimenstion named longitude, lon or x
            !   get meta information on this dimension
            !   assume the name of the dimension corresponds with a variable
            !   read variable
            !   set results in mapgrid struct
            !   check that geographic axes are regular
            integer, intent ( in )              :: ncid
            type (unstrgrid_type), intent (out) :: unstrgrid

            call agnc_get_unstrgrid_ids(ncid, unstrgrid%pnt_dimid, &
                                          unstrgrid%lon_varid, unstrgrid%lat_varid, &
                                          unstrgrid%npoints)

            allocate(unstrgrid%longitude(unstrgrid%npoints))
            call nccheck ( nf90_get_var(ncid, unstrgrid%lon_varid, unstrgrid%longitude) );

            allocate(unstrgrid%latitude(unstrgrid%npoints))
            call nccheck ( nf90_get_var(ncid, unstrgrid%lat_varid, unstrgrid%latitude) )

        end subroutine agnc_get_unstrgrid

        subroutine agnc_get_spcgrid(ncid, spcgrid)
            integer, intent ( in )              :: ncid
            type (spcgrid_type), intent (out)   :: spcgrid

            call agnc_get_spcgrid_ids(ncid, spcgrid%frq_dimid, spcgrid%frq_varid, &
                                          spcgrid%pnt_dimid, &
                                          spcgrid%dir_dimid, spcgrid%dir_varid, &
                                          spcgrid%lon_varid, spcgrid%lat_varid, &
                                          spcgrid%nfreq, spcgrid%npoints, &
                                          spcgrid%ndir)

            if ( spcgrid%frq_varid > 0 ) then
                allocate(spcgrid%frequency(spcgrid%nfreq))
                call nccheck ( nf90_get_var(ncid, spcgrid%frq_varid, spcgrid%frequency) )
                allocate(spcgrid%longitude(spcgrid%npoints))
                call nccheck ( nf90_get_var(ncid, spcgrid%lon_varid, spcgrid%longitude) )
                allocate(spcgrid%latitude(spcgrid%npoints))
                call nccheck ( nf90_get_var(ncid, spcgrid%lat_varid, spcgrid%latitude) )
                if ( spcgrid%dir_varid > 0 ) then
                    allocate(spcgrid%direction(spcgrid%ndir))
                    call nccheck ( nf90_get_var(ncid, spcgrid%dir_varid, spcgrid%direction) )
                end if
            end if
        end subroutine agnc_get_spcgrid

        subroutine agnc_get_recordaxe(ncid, recordaxe)
            integer, intent ( in )                  :: ncid
            type (recordaxe_type), intent(out)      :: recordaxe
            character(len=nf90_max_name)            :: units, raname
            real(kind=8), allocatable, dimension(:) :: dvalues
            ! find time dimid
            call agnc_get_recordaxe_ids(ncid, recordaxe%dimid, recordaxe%varid, recordaxe%ncontent)

            if ( recordaxe%varid < 0 ) then
                call agnc_set_error( 'Could not find record axe variable')
            end if
            call nccheck ( nf90_get_att(ncid, recordaxe%varid, 'units', units) )
            call nccheck ( nf90_inquire_variable(ncid, recordaxe%varid, raname) )
            allocate(recordaxe%content(recordaxe%ncontent))

            if (raname /= 'run') then
                ! ouchouch, ughlyughly. A better but more complicated solution would
                ! be to link in the UDUNITS package. However, they just switched from
                ! version 1.3 to version 2.0 and removed the Perl and Fortan interface.
                !
                ! I leave it as it is now, but this needs to be changed sooner or later
                !
                ! TODO
                select case (units(1:21))
                case ("seconds since 1970-01")
                    call agnc_get_var(ncid, recordaxe%varid, recordaxe%content)
                case ("days since 1970-01-01")
                    allocate(dvalues(recordaxe%ncontent))
                    call agnc_get_var(ncid, recordaxe%varid, dvalues)
                    recordaxe%content = nint(dvalues * 86400)
                    deallocate(dvalues)
                case default
                    call agnc_set_error( 'Unknown time unit'//units(1:21))
                end select
                recordaxe%nstatm = .true.
            else
                call agnc_get_var(ncid, recordaxe%varid, recordaxe%content)
                recordaxe%nstatm = .false.
            end if

            if (recordaxe%ncontent > 1) then
                recordaxe%delta = recordaxe%content(2) - recordaxe%content(1)
            else
                ! You'll have set dt yourself when nt is 1
                recordaxe%delta = 0
            end if
            if (.not. axe_is_regular(recordaxe%content, recordaxe%ncontent, recordaxe%delta)) then
                call agnc_set_error( 'record axe does not appear to be regular')
            end if
        end subroutine agnc_get_recordaxe

        subroutine agnc_get_recordaxe_ids(ncid, dimid, varid, ncontent, raname)
            integer,                             intent( in)           :: ncid
            integer,                             intent(out)           :: dimid
            integer,                   optional, intent(out)           :: varid, ncontent
            character(len=nf90_max_name),   optional, intent(out)      :: raname
            character*6                                                :: dnames(2), dname
            integer                                                    :: i
            logical                                                    :: dim_found

            dim_found = .false.
            dnames = (/ 'time  ', 'run   ' /)

            ! Note that the assumption is that the record axe dimension has a correspomding
            ! variable

            do i = 1, 2
                if ( nf90_inq_dimid(ncid, dnames(i), dimid) == nf90_noerr ) then
                    dname = dnames(i)
                    dim_found = .true.
                    exit
                end if
            end do

            if ( .not. dim_found ) then
                call agnc_set_error( 'record axe could not be found')
            end if

            if (present(ncontent)) then
                call nccheck ( nf90_inquire_dimension(ncid, dimid, dname, ncontent) )
            end if

            if ( present(varid) ) call agnc_get_varid_by_name(ncid, dname, varid)
            if ( present(raname)) raname = dname

        end subroutine agnc_get_recordaxe_ids

        subroutine agnc_get_mapgrid_ids(ncid, londid, latdid, lonvid, latvid, nx, ny, lonname, latname)
            integer             , intent( in)                 :: ncid
            integer             , intent(out)                 :: londid, latdid
            integer, optional   , intent(out)                 :: lonvid, latvid, nx, ny
            character(len=nf90_max_name), optional, intent(out)    :: lonname, latname
            ! find longitude dimid
            if ( nf90_inq_dimid(ncid, "longitude", londid) /= nf90_noerr ) then
                if ( nf90_inq_dimid(ncid, "x", londid) /= nf90_noerr ) then
                    if ( nf90_inq_dimid(ncid, "xc", londid) /= nf90_noerr ) then
                        call nccheck ( nf90_inq_dimid(ncid, "lon", londid) )
                    end if
                end if
            end if
            if ( present(nx) ) then
                call nccheck ( nf90_inquire_dimension(ncid, londid, lonname, nx) )
            else
                call nccheck ( nf90_inquire_dimension(ncid, londid, lonname) )
            end if
            if ( present(lonvid)) call nccheck ( nf90_inq_varid(ncid, lonname, lonvid) )

            ! find latitude dimid
            if ( nf90_inq_dimid(ncid, "latitude", latdid) /= nf90_noerr ) then
                if ( nf90_inq_dimid(ncid, "y", latdid) /= nf90_noerr ) then
                    if ( nf90_inq_dimid(ncid, "yc", latdid) /= nf90_noerr ) then
                        call nccheck ( nf90_inq_dimid(ncid, "lat", latdid) )
                    end if
                end if
            end if
            if ( present(ny) ) then
                call nccheck ( nf90_inquire_dimension(ncid, latdid, latname, ny) )
            else
                call nccheck ( nf90_inquire_dimension(ncid, latdid, latname) )
            end if
            if ( present(latvid)) call nccheck ( nf90_inq_varid(ncid, latname, latvid) )
        end subroutine agnc_get_mapgrid_ids

        subroutine agnc_get_unstrgrid_ids(ncid, pnt_dimid, lonvid, latvid, npoints)
            integer             , intent( in)           :: ncid
            integer, optional   , intent(out)           :: pnt_dimid
            integer, optional   , intent(out)           :: lonvid, latvid, npoints
            integer                                     :: londid, latdid
            character(len=nf90_max_name)                :: lonname, latname

            if ( present(pnt_dimid) ) then
                if ( nf90_inq_dimid(ncid, "points", pnt_dimid) /= nf90_noerr ) then
                    call agnc_set_error( "No 'points' dimension found")
                end if

                if ( present(npoints) ) then
                    call nccheck ( nf90_inquire_dimension(ncid, pnt_dimid, len=npoints) )
                end if
            end if

            if ( present(lonvid) ) then
                ! find longitude varid
                if ( nf90_inq_varid(ncid, "longitude", lonvid) /= nf90_noerr ) then
                    if ( nf90_inq_varid(ncid, "x", lonvid) /= nf90_noerr ) then
                        call nccheck ( nf90_inq_varid(ncid, "lon", lonvid) )
                    end if
                end if
            end if

            if ( present( latvid) ) then
                ! find latitude varid
                if ( nf90_inq_varid(ncid, "latitude", latvid) /= nf90_noerr ) then
                    if ( nf90_inq_varid(ncid, "y", latvid) /= nf90_noerr ) then
                        call nccheck ( nf90_inq_varid(ncid, "lat", latvid) )
                    end if
                end if
            end if

        end subroutine agnc_get_unstrgrid_ids

        subroutine agnc_get_spcgrid_ids(ncid, frq_dimid, frq_varid, pnt_dimid, dir_dimid, dir_varid, &
                                      lon_varid, lat_varid, nfreq, npoints, ndir)
            integer             , intent( in)           :: ncid
            integer, optional   , intent(out)           :: frq_dimid, frq_varid, &
                                                           pnt_dimid, &
                                                           lon_varid, lat_varid, &
                                                           dir_dimid, dir_varid
            integer, optional   , intent(out)           :: nfreq, npoints, ndir


            if ( present(frq_dimid) ) then
                if ( nf90_inq_dimid(ncid, "frequency", frq_dimid) /= nf90_noerr ) then
                    print *, "No 'frequency' dimension found"
                    frq_dimid = -1
                    if ( present( frq_varid ) ) frq_varid = -1
                    return
                end if

                if ( present(nfreq) ) then
                    call nccheck ( nf90_inquire_dimension(ncid, frq_dimid, len=nfreq) )
                end if
            end if

            if ( present(frq_varid) ) then
              if ( nf90_inq_varid(ncid, "frequency", frq_varid) /= nf90_noerr ) then
                call agnc_set_error( "No 'frequency' variable found")
                frq_varid = -1
                return
              end if
            end if

            if ( present(pnt_dimid) ) then
                if ( nf90_inq_dimid(ncid, "points", pnt_dimid) /= nf90_noerr ) then
                    call agnc_set_error( "No 'points' dimension found")
                    return
                end if

                if ( present(npoints) ) then
                    call nccheck ( nf90_inquire_dimension(ncid, pnt_dimid, len=npoints) )
                end if
            end if

            if ( present(lon_varid) ) then
                ! find longitude dimid
                if ( nf90_inq_varid(ncid, "longitude", lon_varid) /= nf90_noerr ) then
                    if ( nf90_inq_varid(ncid, "x", lon_varid) /= nf90_noerr ) then
                        call nccheck ( nf90_inq_varid(ncid, "lon", lon_varid) )
                    end if
                end if
            end if
            ! find latitude dimid
            if ( present(lat_varid) ) then
                ! find latitude dimid
                if ( nf90_inq_varid(ncid, "latitude", lat_varid) /= nf90_noerr ) then
                    if ( nf90_inq_varid(ncid, "y", lat_varid) /= nf90_noerr ) then
                        call nccheck ( nf90_inq_varid(ncid, "lat", lat_varid) )
                    end if
                end if
            end if

            if ( present(dir_dimid) ) then
                if ( nf90_inq_dimid(ncid, "direction", dir_dimid) /= nf90_noerr ) then
                    dir_dimid = -1;
                    dir_varid = -1;
                    if (present(ndir)) ndir=0;
                    return
                end if

                if ( present(ndir) ) then
                    call nccheck ( nf90_inquire_dimension(ncid, dir_dimid, len=ndir) )
                end if
            end if

            if ( present(dir_varid) ) then
                call nccheck ( nf90_inq_varid(ncid, "direction", dir_varid) )
            end if

        end subroutine agnc_get_spcgrid_ids

        subroutine agnc_get_var1d_integer(ncid, varid, values)
            integer,                            intent( in)            :: ncid, varid
            integer(kind=4),dimension(:),       intent(inout)          :: values
            if (size(values) < 1) return
            call nccheck ( nf90_get_var(ncid, varid, values ) )
        end subroutine agnc_get_var1d_integer

        subroutine agnc_get_var2d_integer(ncid, varid, values)
            integer,                            intent( in)            :: ncid, varid
            integer(kind=4),dimension(:,:),     intent(inout)          :: values
            if (size(values) < 1) return
            call nccheck ( nf90_get_var(ncid, varid, values ) )
        end subroutine agnc_get_var2d_integer

        subroutine agnc_get_var3d_integer(ncid, varid, values)
            integer,                            intent( in)            :: ncid, varid
            integer(kind=4),dimension(:,:,:),   intent(inout)          :: values
            if (size(values) < 1) return
            call nccheck ( nf90_get_var(ncid, varid, values ) )
        end subroutine agnc_get_var3d_integer

        subroutine agnc_get_var4d_integer(ncid, varid, values)
            integer,                            intent( in)            :: ncid, varid
            integer(kind=4),dimension(:,:,:,:), intent(inout)          :: values
            if (size(values) < 1) return
            call nccheck ( nf90_get_var(ncid, varid, values ) )
        end subroutine agnc_get_var4d_integer

        subroutine agnc_get_var1d_float(ncid, varid, values, noautoscale)
            integer,                            intent( in)            :: ncid, varid
            real(kind=4),dimension(:),          intent(inout)          :: values
            logical, optional,                  intent( in)            :: noautoscale
            real(kind=4)                                               :: add_offset, scale_factor, fill_value
            integer                                                    :: status
            real                                                       :: eps = 1e-8
            if (size(values) < 1) return
            call nccheck ( nf90_get_var(ncid, varid, values ) )
            if ( present(noautoscale) ) then
                if ( noautoscale ) return
            end if
            call get_scalies_float(ncid, varid, add_offset, scale_factor, fill_value)
            where ( abs( values - fill_value ) > eps ) values = values * scale_factor + add_offset
        end subroutine agnc_get_var1d_float

        subroutine agnc_get_var2d_float(ncid, varid, values, noautoscale)
            integer,                            intent( in)            :: ncid, varid
            real(kind=4),dimension(:,:),        intent(inout)          :: values
            logical, optional,                  intent( in)            :: noautoscale
            real(kind=4)                                               :: add_offset, scale_factor, fill_value
            integer                                                    :: status
            real                                                       :: eps = 1e-8
            if (size(values) < 1) return
            call nccheck ( nf90_get_var(ncid, varid, values ) )
            if ( present(noautoscale) ) then
                if ( noautoscale ) return
            end if
            call get_scalies_float(ncid, varid, add_offset, scale_factor, fill_value)
            where ( abs( values - fill_value ) > eps ) values = values * scale_factor + add_offset
        end subroutine agnc_get_var2d_float

        subroutine agnc_get_var3d_float(ncid, varid, values, noautoscale)
            integer,                            intent( in)            :: ncid, varid
            real(kind=4),dimension(:,:,:),      intent(inout)          :: values
            logical, optional,                  intent( in)            :: noautoscale
            real(kind=4)                                               :: add_offset, scale_factor, fill_value
            real                                                       :: eps = 1e-8
            if (size(values) < 1) return
            call nccheck ( nf90_get_var(ncid, varid, values ) )
            if ( present(noautoscale) ) then
                if ( noautoscale ) return
            end if
            call get_scalies_float(ncid, varid, add_offset, scale_factor, fill_value)
            where ( abs( values - fill_value ) > eps ) values = values * scale_factor + add_offset

        end subroutine agnc_get_var3d_float

        subroutine agnc_get_var4d_float(ncid, varid, values, noautoscale)
            integer,                            intent( in)            :: ncid, varid
            real(kind=4),dimension(:,:,:,:),    intent(inout)          :: values
            logical, optional,                  intent( in)            :: noautoscale
            real(kind=4)                                               :: add_offset, scale_factor, fill_value
            real                                                       :: eps = 1e-8
            if (size(values) < 1) return
            call nccheck ( nf90_get_var(ncid, varid, values ) )
            if ( present(noautoscale) ) then
                if ( noautoscale ) return
            end if
            call get_scalies_float(ncid, varid, add_offset, scale_factor, fill_value)
            where ( abs( values - fill_value ) > eps ) values = values * scale_factor + add_offset

        end subroutine agnc_get_var4d_float
!
! double get_var functions
!
        subroutine agnc_get_var1d_double(ncid, varid, values, noautoscale)
            integer,                            intent( in)            :: ncid, varid
            real(kind=8),dimension(:),          intent(inout)          :: values
            logical, optional,                  intent( in)            :: noautoscale
            real(kind=8)                                               :: add_offset, scale_factor, fill_value
            integer                                                    :: status
            real                                                       :: eps = 1e-8
            if (size(values) < 1) return
            call nccheck ( nf90_get_var(ncid, varid, values ) )
            if ( present(noautoscale) ) then
                if ( noautoscale ) return
            end if
            call get_scalies_double(ncid, varid, add_offset, scale_factor, fill_value)
            where ( abs( values - fill_value ) > eps ) values = values * scale_factor + add_offset
        end subroutine agnc_get_var1d_double

        subroutine agnc_get_var2d_double(ncid, varid, values, noautoscale)
            integer,                            intent( in)            :: ncid, varid
            real(kind=8),dimension(:,:),        intent(inout)          :: values
            logical, optional,                  intent( in)            :: noautoscale
            real(kind=8)                                               :: add_offset, scale_factor, fill_value
            integer                                                    :: status
            real                                                       :: eps = 1e-8
            if (size(values) < 1) return
            call nccheck ( nf90_get_var(ncid, varid, values ) )
            if ( present(noautoscale) ) then
                if ( noautoscale ) return
            end if
            call get_scalies_double(ncid, varid, add_offset, scale_factor, fill_value)
            where ( abs( values - fill_value ) > eps ) values = values * scale_factor + add_offset
        end subroutine agnc_get_var2d_double

        subroutine agnc_get_var3d_double(ncid, varid, values, noautoscale)
            integer,                            intent( in)            :: ncid, varid
            real(kind=8),dimension(:,:,:),      intent(inout)          :: values
            logical, optional,                  intent( in)            :: noautoscale
            real(kind=8)                                               :: add_offset, scale_factor, fill_value
            real                                                       :: eps = 1e-8
            if (size(values) < 1) return
            call nccheck ( nf90_get_var(ncid, varid, values ) )
            if ( present(noautoscale) ) then
                if ( noautoscale ) return
            end if
            call get_scalies_double(ncid, varid, add_offset, scale_factor, fill_value)
            where ( abs( values - fill_value ) > eps ) values = values * scale_factor + add_offset

        end subroutine agnc_get_var3d_double

        subroutine agnc_get_var4d_double(ncid, varid, values, noautoscale)
            integer,                            intent( in)            :: ncid, varid
            real(kind=8),dimension(:,:,:,:),    intent(inout)          :: values
            logical, optional,                  intent( in)            :: noautoscale
            real(kind=8)                                               :: add_offset, scale_factor, fill_value
            real                                                       :: eps = 1e-8
            if (size(values) < 1) return
            call nccheck ( nf90_get_var(ncid, varid, values ) )
            if ( present(noautoscale) ) then
                if ( noautoscale ) return
            end if
            call get_scalies_double(ncid, varid, add_offset, scale_factor, fill_value)
            where ( abs( values - fill_value ) > eps ) values = values * scale_factor + add_offset

        end subroutine agnc_get_var4d_double

        subroutine get_scalies_float(ncid, varid, add_offset, scale_factor, fill_value)
            integer,                            intent( in)            :: ncid, varid
            real(kind=4),                       intent(out)            :: add_offset, scale_factor, fill_value
            integer                                                    :: xtype
            call nccheck ( nf90_inquire_variable(ncid, varid, xtype=xtype))
            if ( nf90_get_att(ncid, varid, "_FillValue", fill_value) /= NF90_NOERR ) then
                select case (xtype)
                case (NF90_BYTE)
                    fill_value = NF90_FILL_BYTE
                case (NF90_SHORT)
                    fill_value = NF90_FILL_SHORT
                case (NF90_INT)
                    fill_value = NF90_FILL_INT
                case default
                    fill_value = NF90_FILL_FLOAT
                end select
            end if
            if ( nf90_get_att(ncid, varid, "add_offset", add_offset) /= NF90_NOERR ) add_offset = 0.
            if ( nf90_get_att(ncid, varid, "scale_factor", scale_factor) /= NF90_NOERR ) scale_factor = 1.

        end subroutine get_scalies_float

        subroutine get_scalies_double(ncid, varid, add_offset, scale_factor, fill_value)
            integer,                            intent( in)            :: ncid, varid
            real(kind=8),                       intent(out)            :: add_offset, scale_factor, fill_value
            integer                                                    :: xtype
            call nccheck ( nf90_inquire_variable(ncid, varid, xtype=xtype))
            if ( nf90_get_att(ncid, varid, "_FillValue", fill_value) /= NF90_NOERR ) then
                select case (xtype)
                case (NF90_BYTE)
                    fill_value = NF90_FILL_BYTE
                case (NF90_SHORT)
                    fill_value = NF90_FILL_SHORT
                case (NF90_INT)
                    fill_value = NF90_FILL_INT
                case default
                    fill_value = NF90_FILL_DOUBLE
                end select
            end if
            if ( nf90_get_att(ncid, varid, "add_offset", add_offset) /= NF90_NOERR ) add_offset = 0.
            if ( nf90_get_att(ncid, varid, "scale_factor", scale_factor) /= NF90_NOERR ) scale_factor = 1.


        end subroutine get_scalies_double

        subroutine agnc_get_spcvars_by_point(ncid, points, energy, theta, spr)
            integer,                           intent( in)      :: ncid
            integer,         dimension(:),     intent( in)      :: points
            real(kind=4),    dimension(:,:,:), intent(out)      :: energy, theta, spr
            real(kind=4),    dimension(:,:),     allocatable    :: scale_energy
            integer                                             :: i, varid_scale_energy, &
                                                                   varid_energy, varid_theta, &
                                                                   varid_spr
            real(kind=4)                                        :: sf_theta, sf_spr, &
                                                                   ao_theta, ao_spr, &
                                                                   fv_theta, fv_spr
            type(spcgrid_type)                                  :: spcgrid
            type (recordaxe_type)                                :: recordaxe
            ! energy(nfreq, npoint, ntime)
            call agnc_collect_spcmeta(ncid, spcgrid, recordaxe, varid_scale_energy, varid_energy,&
                                      varid_theta, varid_spr, ao_theta, sf_theta, fv_theta,&
                                      ao_spr, sf_spr, fv_spr)
            allocate(scale_energy(size(points), recordaxe%ncontent))
            do i=1, size(points)
                call nccheck ( nf90_get_var(ncid, varid_scale_energy, scale_energy(i,:), (/points(i), 1/), &
                                    (/1, recordaxe%ncontent/) ) )
                call nccheck ( nf90_get_var(ncid, varid_energy, energy(:,i,:), (/1, points(i), 1/), &
                                    (/spcgrid%nfreq, 1, recordaxe%ncontent/) ) )
                call nccheck ( nf90_get_var(ncid, varid_theta, theta(:,i,:), (/1, points(i), 1/), &
                                    (/spcgrid%nfreq, 1, recordaxe%ncontent/) ) )
                call nccheck ( nf90_get_var(ncid, varid_spr, spr(:,i,:), (/1, points(i), 1/), &
                                    (/spcgrid%nfreq, 1, recordaxe%ncontent/) ) )
            end do
            where ( abs(energy - AGNC_FILL_SHORT) > epsilon(1.) .and. abs(energy - NF90_FILL_FLOAT) > epsilon(1.) )
                energy = nint(energy) * spread(scale_energy,1,spcgrid%nfreq)
            elsewhere
                energy = AGNC_DUMMY
            end where
            where ( abs(theta - fv_theta) > epsilon(1.) .and. abs(theta - NF90_FILL_FLOAT) > epsilon(1.) )
                theta = theta * sf_theta + ao_theta
            elsewhere
                theta = AGNC_DUMMY
            end where
            where ( abs(spr - fv_spr) > epsilon(1.) .and. abs(spr - NF90_FILL_FLOAT) > epsilon(1.) )
                spr = spr * sf_spr + ao_spr
            elsewhere
                spr = AGNC_DUMMY
            end where
            deallocate(scale_energy)
        end subroutine agnc_get_spcvars_by_point

        subroutine agnc_get_spcvars_vector(ncid, time_indices, energy, theta, spr)
            integer,                           intent( in)      :: ncid
            integer,         dimension(:),     intent( in)      :: time_indices
            real(kind=4),    dimension(:,:,:), intent(out)      :: energy, theta, spr
            real(kind=4),    dimension(:,:),   allocatable      :: scale_energy
            integer                                             :: i, varid_scale_energy, &
                                                                   varid_energy, varid_theta, &
                                                                   varid_spr
            real(kind=4)                                        :: sf_theta, sf_spr, &
                                                                   ao_theta, ao_spr, &
                                                                   fv_theta, fv_spr
            type(spcgrid_type)                                  :: spcgrid
            type (recordaxe_type)                                :: recordaxe
            ! energy(nfreq, npoint, ntime)
            call agnc_collect_spcmeta(ncid, spcgrid, recordaxe, varid_scale_energy, varid_energy,&
                                      varid_theta, varid_spr, ao_theta, sf_theta, fv_theta,&
                                      ao_spr, sf_spr, fv_spr)
            allocate(scale_energy(spcgrid%npoints, size(time_indices)))

            do i=1, size(time_indices)
                call nccheck ( nf90_get_var(ncid, varid_scale_energy, scale_energy(:,i), (/1, time_indices(i)/), &
                                    (/spcgrid%npoints, 1/) ) )
                call nccheck ( nf90_get_var(ncid, varid_energy, energy(:,:,i), (/1, 1, time_indices(i)/), &
                                    (/spcgrid%nfreq, spcgrid%npoints, 1/) ) )
                call nccheck ( nf90_get_var(ncid, varid_theta, theta(:,:,i), (/1, 1, time_indices(i)/), &
                                    (/spcgrid%nfreq, spcgrid%npoints, 1/) ) )
                call nccheck ( nf90_get_var(ncid, varid_spr, spr(:,:,i), (/1, 1, time_indices(i)/), &
                                    (/spcgrid%nfreq, spcgrid%npoints, 1/) ) )
            end do

            where ( abs(energy - AGNC_FILL_SHORT) > epsilon(1.) .and. abs(energy - NF90_FILL_FLOAT) > epsilon(1.) )
                energy = nint(energy) * spread(scale_energy,1,spcgrid%nfreq)
            elsewhere
                energy = AGNC_DUMMY
            end where

            where ( abs(theta - fv_theta) > epsilon(1.) .and. abs(theta - NF90_FILL_FLOAT) > epsilon(1.) )
                theta = theta * sf_theta + ao_theta
            elsewhere
                theta = AGNC_DUMMY
            end where

            where ( abs(spr - fv_spr) > epsilon(1.) .and. abs(spr - NF90_FILL_FLOAT) > epsilon(1.) )
                spr = spr * sf_spr + ao_spr
            elsewhere
                spr = AGNC_DUMMY
            end where
            deallocate(scale_energy)
        end subroutine agnc_get_spcvars_vector

        subroutine agnc_get_spcvars_scalar(ncid, time_index, energy, theta, spr)
            integer,                           intent( in)      :: ncid
            integer,                           intent( in)      :: time_index
            real(kind=4),    dimension(:,:),   intent(out)      :: energy, theta, spr
            real(kind=4),    dimension(:),     allocatable      :: scale_energy
            integer                                             :: i, varid_scale_energy, &
                                                                   varid_energy, varid_theta, &
                                                                   varid_spr
            real(kind=4)                                        :: sf_theta, sf_spr, &
                                                                   ao_theta, ao_spr, &
                                                                   fv_theta, fv_spr
            type(spcgrid_type)                                  :: spcgrid
            type (recordaxe_type)                                :: recordaxe
            ! energy(nfreq, npoint, ntime)
            call agnc_collect_spcmeta(ncid, spcgrid, recordaxe, varid_scale_energy, varid_energy,&
                                      varid_theta, varid_spr, ao_theta, sf_theta, fv_theta,&
                                      ao_spr, sf_spr, fv_spr)
            allocate(scale_energy(spcgrid%npoints))

            call nccheck ( nf90_get_var(ncid, varid_scale_energy, scale_energy, (/1, time_index/), &
                                (/spcgrid%npoints, 1/) ) )
            call nccheck ( nf90_get_var(ncid, varid_energy, energy, (/1, 1, time_index/), &
                                (/spcgrid%nfreq, spcgrid%npoints, 1/) ) )
            call nccheck ( nf90_get_var(ncid, varid_theta, theta, (/1, 1, time_index/), &
                                (/spcgrid%nfreq, spcgrid%npoints, 1/) ) )
            call nccheck ( nf90_get_var(ncid, varid_spr, spr, (/1, 1, time_index/), &
                                (/spcgrid%nfreq, spcgrid%npoints, 1/) ) )

            where ( abs(energy - AGNC_FILL_SHORT) > epsilon(1.) .and. abs(energy - NF90_FILL_FLOAT) > epsilon(1.) )
                energy = nint(energy) * spread(scale_energy,1,spcgrid%nfreq)
            elsewhere
                energy = AGNC_DUMMY
            end where
            where ( abs(theta - fv_theta) > epsilon(1.) .and. abs(theta - NF90_FILL_FLOAT) > epsilon(1.) )
                theta = theta * sf_theta + ao_theta
            elsewhere
                theta = AGNC_DUMMY
            end where
            where ( abs(spr - fv_spr) > epsilon(1.) .and. abs(spr - NF90_FILL_FLOAT) > epsilon(1.) )
                spr = spr * sf_spr + ao_spr
            elsewhere
                spr = AGNC_DUMMY
            end where
            deallocate(scale_energy)
        end subroutine agnc_get_spcvars_scalar

!
!       ------ insert/append data ------
!
!       agnc_add_mapdata(ncid, varname, time, values)     done
!       agnc_add_spcdata(ncid, e1, th1, spr1)             todo
        subroutine agnc_add_mapdata2d_float(ncid, varname, ti, values, dummyvalue, skip_range_error)
            ! if the optional dummyvalue is given, the values are searched for this value
            ! and replaced with the _FillValue
            ! ti is the index where the data write should be started
            integer,                        intent( in)         :: ncid, ti
            character(len=*),               intent( in)         :: varname
            real(kind=4), optional,         intent( in)         :: dummyvalue
            logical,      optional                              :: skip_range_error
            real(kind=4), dimension(:,:)                        :: values
            real(kind=4)                                        :: add_offset, scale_factor, &
                                                                   fill_value
            integer                                             :: varid, nf90_stat

            call agnc_get_varid_by_name(ncid, varname, varid)

            call get_scalies_float(ncid, varid, add_offset, scale_factor, fill_value)
            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            where ( abs(values - fill_value) > epsilon(1.) ) &
                values = nint( (values - add_offset) / scale_factor)
            nf90_stat = nf90_put_var(ncid, varid, values, (/1, 1, ti/))
            if ( nf90_stat /= NF90_NOERR ) then
                if ( values_in_range_float(ncid, varid, reshape(values, (/size(values)/)), &
                        fill_value, scale_factor, add_offset, varname) ) then
                    call nccheck (nf90_stat)
                else
                    if (.not. present(skip_range_error)) call agnc_set_error( nf90_strerror(nf90_stat))
                end if
            end if
        end subroutine agnc_add_mapdata2d_float

        subroutine agnc_add_mapdata3d_float(ncid, varname, ti, values, dummyvalue, skip_range_error)
            ! if the optional dummyvalue is given, the values are searched for this value
            ! and replaced with the _FillValue
            ! ti is the index where the data write should be started
            integer,                        intent( in)         :: ncid, ti
            character(len=*),               intent( in)         :: varname
            real(kind=4), optional,         intent( in)         :: dummyvalue
            logical,      optional                              :: skip_range_error
            real(kind=4), dimension(:,:,:)                      :: values
            real(kind=4)                                        :: add_offset, scale_factor, &
                                                                   fill_value
            integer                                             :: varid, nf90_stat
            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_float(ncid, varid, add_offset, scale_factor, fill_value)

            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            where ( abs(values - fill_value) > epsilon(1.) ) &
                values = nint( (values - add_offset) / scale_factor)
            nf90_stat = nf90_put_var(ncid, varid, values, (/1, 1, ti/))
            if ( nf90_stat /= NF90_NOERR ) then
                if ( values_in_range_float(ncid, varid, reshape(values, (/size(values)/)), &
                        fill_value, scale_factor, add_offset, varname) ) then
                    call agnc_set_error( nf90_strerror(nf90_stat))
                else
                    if (.not. present(skip_range_error)) call agnc_set_error( nf90_strerror(nf90_stat))
                end if
            end if
        end subroutine agnc_add_mapdata3d_float

        subroutine agnc_add_mapdata2d_double(ncid, varname, ti, values, dummyvalue, skip_range_error)
            ! if the optional dummyvalue is given, the values are searched for this value
            ! and replaced with the _FillValue
            ! ti is the index where the data write should be started
            integer,                        intent( in)         :: ncid, ti
            character(len=*),               intent( in)         :: varname
            real(kind=8), optional,         intent( in)         :: dummyvalue
            logical,      optional                              :: skip_range_error
            real(kind=8), dimension(:,:)                        :: values
            real(kind=8)                                        :: add_offset, scale_factor, &
                                                                   fill_value
            integer                                             :: varid, nf90_stat
            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_double(ncid, varid, add_offset, scale_factor, fill_value)

            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            where ( abs(values - fill_value) > epsilon(1.) ) &
                values = nint( (values - add_offset) / scale_factor)

            nf90_stat = nf90_put_var(ncid, varid, values, (/1, 1, ti/))
            if ( nf90_stat /= NF90_NOERR ) then
                if ( values_in_range_double(ncid, varid, reshape(values, (/size(values)/)), &
                        fill_value, scale_factor, add_offset, varname) ) then
                    call agnc_set_error( nf90_strerror(nf90_stat))
                else
                    if (.not. present(skip_range_error)) call agnc_set_error( nf90_strerror(nf90_stat))
                end if
            end if
        end subroutine agnc_add_mapdata2d_double

        subroutine agnc_add_mapdata3d_double(ncid, varname, ti, values, dummyvalue, skip_range_error)
            ! if the optional dummyvalue is given, the values are searched for this value
            ! and replaced with the _FillValue
            ! ti is the index where the data write should be started
            integer,                        intent( in)         :: ncid, ti
            character(len=*),               intent( in)         :: varname
            real(kind=8), optional,         intent( in)         :: dummyvalue
            logical,      optional                              :: skip_range_error
            real(kind=8), dimension(:,:,:)                      :: values
            real(kind=8)                                        :: add_offset, scale_factor, &
                                                                   fill_value
            integer                                             :: varid, nf90_stat
            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_double(ncid, varid, add_offset, scale_factor, fill_value)

            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            where ( abs(values - fill_value) > epsilon(1.) ) &
                values = nint( (values - add_offset) / scale_factor)

            nf90_stat = nf90_put_var(ncid, varid, values, (/1, 1, ti/))
            if ( nf90_stat /= NF90_NOERR ) then
                if ( values_in_range_double(ncid, varid, reshape(values, (/size(values)/)), &
                        fill_value, scale_factor, add_offset, varname) ) then
                    call agnc_set_error( nf90_strerror(nf90_stat))
                else
                    if (.not. present(skip_range_error)) call agnc_set_error( nf90_strerror(nf90_stat))
                end if
            end if
        end subroutine agnc_add_mapdata3d_double

        subroutine agnc_add_spcdata_density(ncid, ti, e2d)
            integer,                            intent( in) :: ncid, ti
            real(kind=4),   dimension(:,:)                  :: e2d
            real(kind=4),   dimension(:,:,:), allocatable   :: density
            real(kind=4),   dimension(:), allocatable       :: scale_density
            integer                                         :: varid_scale_density, &
                                                               varid_density
            real(kind=4)                                    :: fv_density
            integer                                         :: ip, ith, ik, isp
            type(spcgrid_type)                              :: spcgrid
            type(recordaxe_type)                            :: recordaxe
            logical                                         :: do_scale

            call agnc_collect_spcmeta_2d(ncid, spcgrid, recordaxe, varid_density, &
                                        varid_scale_density, fv_density)

            do_scale = .false.
            if ( varid_scale_density > 0 ) do_scale = .true.

            allocate(density(spcgrid%nfreq, spcgrid%ndir, spcgrid%npoints))

            ! Note the assumption that AGNC_FILL_SHORT is negative
            if ( do_scale ) then
                where ( abs(e2d - AGNC_DUMMY) < epsilon(1.) ) e2d = AGNC_FILL_SHORT
                allocate(scale_density(spcgrid%npoints))
                scale_density = maxval(e2d, 1) / (2**15-1)
                where (scale_density < epsilon(1.) ) scale_density=epsilon(1.)
            else
                where ( abs(e2d - AGNC_DUMMY) < epsilon(1.) ) e2d = NF90_FILL_FLOAT
            end if

            do ip=1, spcgrid%npoints
                do ith=1, spcgrid%ndir
                    do ik=1, spcgrid%nfreq
                        isp = ith + (ik - 1) * spcgrid%ndir
                        if ( do_scale .and. abs(e2d(isp,ip) - AGNC_FILL_SHORT) > epsilon(1.) ) then
                            density(ik,ith,ip) = nint(e2d(isp,ip) / scale_density(ip))
                        else
                            density(ik,ith,ip) = e2d(isp,ip)
                        end if
                    end do
                end do
            end do

            call nccheck ( nf90_put_var(ncid, varid_density, density, (/1,1,1,ti/)) )
            deallocate(density)

            if (do_scale) then
                call nccheck ( nf90_put_var(ncid, varid_scale_density, scale_density, (/1, ti/)) )
                deallocate(scale_density)
            end if

        end subroutine agnc_add_spcdata_density

        subroutine agnc_add_spcdata_3d(ncid, ti, energy, theta, spr)
            integer,                           intent( in)      :: ncid, ti
            real(kind=4),    dimension(:,:,:)                   :: energy, theta, spr
            real(kind=4),    dimension(:,:),   allocatable      :: scale_energy
            integer                                             :: varid_scale_energy, &
                                                                   varid_energy, varid_theta, &
                                                                   varid_spr
            real(kind=4)                                        :: sf_theta, sf_spr, &
                                                                   ao_theta, ao_spr, &
                                                                   fv_theta, fv_spr
            type(spcgrid_type)                                  :: spcgrid
            type (recordaxe_type)                               :: recordaxe
            logical                                             :: do_scale

            ! energy(nfreq, npoint, ntime)
            call agnc_collect_spcmeta(ncid, spcgrid, recordaxe, varid_scale_energy, varid_energy,&
                                      varid_theta, varid_spr, ao_theta, sf_theta, fv_theta,&
                                      ao_spr, sf_spr, fv_spr)

            do_scale = .false.
            if ( varid_scale_energy > 0 ) then
                do_scale = .true.

                allocate(scale_energy(spcgrid%npoints, size(energy,3)))

                where ( abs(energy - AGNC_DUMMY) > epsilon(1.) ) energy = AGNC_FILL_SHORT

                ! Note the assumption that AGNC_FILL_SHORT is negative
                scale_energy = maxval(energy, 1) / (2**15-1)

                where (scale_energy < epsilon(1.) ) scale_energy=epsilon(1.)

                where ( abs(energy - AGNC_FILL_SHORT) > epsilon(1.) ) &
                    energy = nint(energy / spread(scale_energy,1,spcgrid%nfreq))

                call nccheck ( nf90_put_var(ncid, varid_scale_energy, scale_energy, (/1, ti/)) )
            else
                where ( abs(energy - AGNC_DUMMY) < epsilon(1.) ) energy = NF90_FILL_FLOAT
            end if

            where ( abs(theta - AGNC_DUMMY) > epsilon(1.) )
                theta = nint((theta - ao_theta) / sf_theta)
            elsewhere
                theta = AGNC_FILL_BYTE
            end where

            where ( abs(spr - AGNC_DUMMY) > epsilon(1.) )
                spr = nint((spr - ao_spr) / sf_spr)
            elsewhere
                spr = AGNC_FILL_BYTE
            end where

            call nccheck ( nf90_put_var(ncid, varid_energy, energy, (/1,1,ti/)) )
            call nccheck ( nf90_put_var(ncid, varid_theta, theta, (/1,1,ti/)) )
            call nccheck ( nf90_put_var(ncid, varid_spr, spr, (/1,1,ti/)) )

            if ( do_scale ) deallocate(scale_energy)

        end subroutine agnc_add_spcdata_3d

        subroutine agnc_add_spcdata_2d(ncid, ti, energy, theta, spr)
            integer,                        intent(   in)       :: ncid, ti
            real(kind=4),    dimension(:,:),intent(inout)       :: energy, theta, spr
            real(kind=4),    dimension(:),  allocatable         :: scale_energy
            integer                                             :: varid_scale_energy, &
                                                                   varid_energy, varid_theta, &
                                                                   varid_spr, ip, ik
            real(kind=4)                                        :: sf_theta, sf_spr, &
                                                                   ao_theta, ao_spr, &
                                                                   fv_theta, fv_spr
            type(spcgrid_type)                                  :: spcgrid
            type (recordaxe_type)                               :: recordaxe
            logical                                             :: do_scale

            ! energy(nfreq, npoint)
            call agnc_collect_spcmeta(ncid, spcgrid, recordaxe, varid_scale_energy, varid_energy,&
                                      varid_theta, varid_spr, ao_theta, sf_theta, fv_theta,&
                                      ao_spr, sf_spr, fv_spr)

            do_scale = .false.
            if ( varid_scale_energy > 0 ) then
                do_scale = .true.

                allocate(scale_energy(spcgrid%npoints))

                where ( abs(energy - AGNC_DUMMY) <= epsilon(1.) ) energy = AGNC_FILL_SHORT

                ! Note the assumption that AGNC_FILL_SHORT is negative
                scale_energy = maxval(energy, 1) / (2**15-1)
                where (scale_energy < epsilon(1.) ) scale_energy=epsilon(1.)

                do ip=1, spcgrid%npoints
                    do ik=1, spcgrid%nfreq
                        if ( abs(energy(ik,ip) - AGNC_FILL_SHORT) > epsilon(1.)) energy(ik,ip) = nint(energy(ik,ip) / scale_energy(ip))
                    end do
                end do

                call nccheck ( nf90_put_var(ncid, varid_scale_energy, scale_energy, (/1, ti/)) )
            else
                where ( abs(energy - AGNC_DUMMY) <= epsilon(1.) ) energy = NF90_FILL_FLOAT
            end if

            where ( abs(theta - AGNC_DUMMY) > epsilon(1.) )
                theta = nint((theta - ao_theta) / sf_theta)
            elsewhere
                theta = AGNC_FILL_BYTE
            end where

            where ( abs(spr - AGNC_DUMMY) > epsilon(1.) )
                spr = nint((spr - ao_spr) / sf_spr)
            elsewhere
                spr = AGNC_FILL_BYTE
            end where

            call nccheck ( nf90_put_var(ncid, varid_energy, energy, (/1,1,ti/)) )
            call nccheck ( nf90_put_var(ncid, varid_theta, theta, (/1,1,ti/)) )
            call nccheck ( nf90_put_var(ncid, varid_spr, spr, (/1,1,ti/)) )

            if ( do_scale ) deallocate(scale_energy)

        end subroutine agnc_add_spcdata_2d

        subroutine agnc_add_pntdata1d_float(ncid, varname, ti, values, dummyvalue, skip_range_error)
            ! if the optional dummyvalue is given, the values are searched for this value
            ! and replaced with the _FillValue
            ! ti is the index where the data write should be started
            integer,                        intent( in)         :: ncid, ti
            character(len=*),               intent( in)         :: varname
            real(kind=4), optional,         intent( in)         :: dummyvalue
            logical,      optional                              :: skip_range_error
            real(kind=4), dimension(:)                          :: values
            real(kind=4)                                        :: add_offset, scale_factor, &
                                                                   fill_value
            integer                                             :: varid, nf90_stat

            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_float(ncid, varid, add_offset, scale_factor, fill_value)

            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            where ( abs(values - fill_value) > epsilon(1.) ) &
                values = nint( (values - add_offset) / scale_factor)

            nf90_stat = nf90_put_var(ncid, varid, values, (/1, ti/))
            if (  nf90_stat /= NF90_NOERR ) then
                if ( values_in_range_float(ncid, varid, reshape(values, (/size(values)/)), &
                        fill_value, scale_factor, add_offset, varname) ) then
                    call agnc_set_error( nf90_strerror(nf90_stat))
                else
                    if (.not. present(skip_range_error)) call agnc_set_error( nf90_strerror(nf90_stat))
                end if
            end if
        end subroutine agnc_add_pntdata1d_float

        subroutine agnc_add_pntdata2d_float(ncid, varname, ti, values, dummyvalue, skip_range_error)
            ! if the optional dummyvalue is given, the values are searched for this value
            ! and replaced with the _FillValue
            ! ti is the index where the data write should be started
            integer,                        intent( in)         :: ncid, ti
            character(len=*),               intent( in)         :: varname
            real(kind=4), optional,         intent( in)         :: dummyvalue
            logical,      optional                              :: skip_range_error
            real(kind=4), dimension(:,:)                        :: values
            real(kind=4)                                        :: add_offset, scale_factor, &
                                                                   fill_value
            integer                                             :: varid, nf90_stat

            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_float(ncid, varid, add_offset, scale_factor, fill_value)

            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            where ( abs(values - fill_value) > epsilon(1.) ) &
                values = nint( (values - add_offset) / scale_factor)
                        nf90_stat = nf90_put_var(ncid, varid, values, (/1, ti/))

            if (  nf90_stat /= NF90_NOERR ) then
                if ( values_in_range_float(ncid, varid, reshape(values, (/size(values)/)), &
                        fill_value, scale_factor, add_offset, varname) ) then
                    call agnc_set_error( nf90_strerror(nf90_stat))
                else
                    if (.not. present(skip_range_error)) call agnc_set_error( nf90_strerror(nf90_stat))
                end if
            end if

        end subroutine agnc_add_pntdata2d_float

        subroutine agnc_add_pntdata3d_float(ncid, varname, ti, values, dummyvalue, skip_range_error)
            ! if the optional dummyvalue is given, the values are searched for this value
            ! and replaced with the _FillValue
            ! ti is the index where the data write should be started
            integer,                        intent( in)         :: ncid, ti
            character(len=*),               intent( in)         :: varname
            real(kind=4), optional,         intent( in)         :: dummyvalue
            logical,      optional                              :: skip_range_error
            real(kind=4), dimension(:,:,:)                      :: values
            real(kind=4)                                        :: add_offset, scale_factor, &
                                                                   fill_value
            integer                                             :: varid, nf90_stat

            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_float(ncid, varid, add_offset, scale_factor, fill_value)

            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            where ( abs(values - fill_value) > epsilon(1.) ) &
                values = nint( (values - add_offset) / scale_factor)

            nf90_stat = nf90_put_var(ncid, varid, values, (/1, ti/))
            if (  nf90_stat /= NF90_NOERR ) then
                if ( values_in_range_float(ncid, varid, reshape(values, (/size(values)/)), &
                        fill_value, scale_factor, add_offset, varname) ) then
                    call agnc_set_error( nf90_strerror(nf90_stat))
                else
                    if (.not. present(skip_range_error)) call agnc_set_error( nf90_strerror(nf90_stat))
                end if
            end if

        end subroutine agnc_add_pntdata3d_float

        subroutine agnc_add_pntdata1d_double(ncid, varname, ti, values, dummyvalue, skip_range_error)
            ! if the optional dummyvalue is given, the values are searched for this value
            ! and replaced with the _FillValue
            ! ti is the index where the data write should be started
            integer,                        intent( in)         :: ncid, ti
            character(len=*),               intent( in)         :: varname
            real(kind=8), optional,         intent( in)         :: dummyvalue
            logical,      optional                              :: skip_range_error
            real(kind=8), dimension(:)                          :: values
            real(kind=8)                                        :: add_offset, scale_factor, &
                                                                   fill_value
            integer                                             :: varid, nf90_stat
            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_double(ncid, varid, add_offset, scale_factor, fill_value)

            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            where ( abs(values - fill_value) > epsilon(1.) ) &
                values = nint( (values - add_offset) / scale_factor)

            nf90_stat = nf90_put_var(ncid, varid, values, (/1, ti/))
            if (  nf90_stat /= NF90_NOERR ) then
                if ( values_in_range_double(ncid, varid, reshape(values, (/size(values)/)), &
                        fill_value, scale_factor, add_offset, varname) ) then
                    call agnc_set_error( nf90_strerror(nf90_stat))
                else
                    if (.not. present(skip_range_error)) call agnc_set_error( nf90_strerror(nf90_stat))
                end if
            end if

        end subroutine agnc_add_pntdata1d_double

        subroutine agnc_add_pntdata2d_double(ncid, varname, ti, values, dummyvalue, skip_range_error)
            ! if the optional dummyvalue is given, the values are searched for this value
            ! and replaced with the _FillValue
            ! ti is the index where the data write should be started
            integer,                        intent( in)         :: ncid, ti
            character(len=*),               intent( in)         :: varname
            real(kind=8), optional,         intent( in)         :: dummyvalue
            logical,      optional                              :: skip_range_error
            real(kind=8), dimension(:,:)                        :: values
            real(kind=8)                                        :: add_offset, scale_factor, &
                                                                   fill_value
            integer                                             :: varid, nf90_stat
            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_double(ncid, varid, add_offset, scale_factor, fill_value)

            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            where ( abs(values - fill_value) > epsilon(1.) ) &
                values = nint( (values - add_offset) / scale_factor)

            nf90_stat = nf90_put_var(ncid, varid, values, (/1, ti/))
            if (  nf90_stat /= NF90_NOERR ) then
                if ( values_in_range_double(ncid, varid, reshape(values, (/size(values)/)), &
                        fill_value, scale_factor, add_offset, varname) ) then
                    call agnc_set_error( nf90_strerror(nf90_stat))
                else
                    if (.not. present(skip_range_error)) call agnc_set_error( nf90_strerror(nf90_stat))
                end if
            end if

        end subroutine agnc_add_pntdata2d_double

        subroutine agnc_add_pntdata3d_double(ncid, varname, ti, values, dummyvalue, skip_range_error)
            ! if the optional dummyvalue is given, the values are searched for this value
            ! and replaced with the _FillValue
            ! ti is the index where the data write should be started
            integer,                        intent( in)         :: ncid, ti
            character(len=*),               intent( in)         :: varname
            real(kind=8), optional,         intent( in)         :: dummyvalue
            logical,      optional                              :: skip_range_error
            real(kind=8), dimension(:,:,:)                      :: values
            real(kind=8)                                        :: add_offset, scale_factor, &
                                                                   fill_value
            integer                                             :: varid, nf90_stat
            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_double(ncid, varid, add_offset, scale_factor, fill_value)

            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            where ( abs(values - fill_value) > epsilon(1.) ) &
                values = nint( (values - add_offset) / scale_factor)

            nf90_stat = nf90_put_var(ncid, varid, values, (/1, ti/))
            if (  nf90_stat /= NF90_NOERR ) then
                if ( values_in_range_double(ncid, varid, reshape(values, (/size(values)/)), &
                        fill_value, scale_factor, add_offset, varname) ) then
                    call agnc_set_error( nf90_strerror(nf90_stat))
                else
                    if (.not. present(skip_range_error)) call agnc_set_error( nf90_strerror(nf90_stat))
                end if
            end if

        end subroutine agnc_add_pntdata3d_double
!
!       ------------- utils -------------
!       name                                            status
!       agnc_get_varid_by_name(ncid, searchfor, varid)  done
!       close_ncfile(ncid)                              done
!       nccheck                                         done
!       axe_is_regular(x, nx, dx)                       done
!       derive_mapgrid(longitude, latitude, mapgrid)    todo
!       derive_mapgrid(longitude, latitude, &   todo
!                points, frequency, direction, spcgrid)
!       compute_1d_spectra(energy, e1, th1, sp1, undef) todo

        subroutine agnc_get_varid_by_name(ncid, searchfor, varid)
            ! search for name in:
            !  1) name of the variable
            !  2) standard_name of the variable
            !  3) short_name of the variable
            integer,            intent(in)          :: ncid
            character(len=*),   intent(in)          :: searchfor
            integer,            intent(out)         :: varid
            integer                                 :: ndims,nvars
            character(len=NF90_MAX_NAME)            :: name

            call nccheck ( nf90_inquire(ncid, ndims, nvars) );
            do varid=1,nvars
                call nccheck (nf90_inquire_variable(ncid, varid, name) )
                if ( trim(name) == trim(searchfor) ) exit

                if ( nf90_get_att(ncid, varid, "standard_name", name) == NF90_NOERR ) then
                    if ( trim(name) == trim(searchfor) ) exit
                else

                    if (nf90_get_att(ncid, varid, "short_name", name) == NF90_NOERR ) then
                        if ( trim(name) == trim(searchfor) ) exit
                    end if
                end if
            end do
            if ( varid > nvars ) then
                varid = -1;
            end if
        end subroutine agnc_get_varid_by_name

        subroutine close_ncfile(ncid)
            integer         , intent(in)       :: ncid
            call nccheck( nf90_close(ncid) )
        end subroutine close_ncfile

        function axe_is_regular_integer(x, nx, dx, fill_check) result (isregular)
            integer,   intent(in)           :: x(:),dx,nx
            logical,   intent(in), optional :: fill_check
            logical                         :: isregular
            integer                         :: dx2(nx-1)

            isregular = .false.
            dx2 = x(2:nx) - x(1:nx-1)
            if ( any( (dx2 - dx) /= 0 )) return

            if ( present(fill_check) ) then
                if  ( fill_check .and. all( (x - NF90_FILL_INT) == 0) ) return
            end if
            isregular = .true.

        end function axe_is_regular_integer

        function axe_is_regular_float(x, nx, dx, fill_check) result (isregular)
            real(kind=4),intent(in)             :: x(:),dx
            integer,     intent(in)             :: nx
            logical,     intent(in), optional   :: fill_check
            logical                             :: isregular
            real(kind=4)                        :: dx2(nx-1)

            isregular = .false.
            dx2 = x(2:nx) - x(1:nx-1)
            ! Should be, epsilon(1.0). However, our own aden grid does not comply
            if ( any( abs(dx2 - dx) > 1e-4 )) return
            if ( present(fill_check) ) then
                if ( fill_check .and. all( abs(x - NF90_FILL_FLOAT) < epsilon(1.0)) ) return
            end if
            isregular = .true.

        end function axe_is_regular_float

        function axe_is_regular_float_2d(x, nx, dx, fill_check) result (isregular)
            real(kind=4),intent(in)             :: x(:,:),dx
            integer,     intent(in)             :: nx
            logical,     intent(in), optional   :: fill_check
            logical                             :: isregular
            real(kind=4)                        :: dx2(nx-1)

            ! The reason for 2D meshes is that we can use the SWAN rotated grid feature
            !
            ! I have never used this before so I would really know how to deal with it
            ! as I'm not sure how x0, x1, and alpha are implemented.
            isregular = .true.

        end function axe_is_regular_float_2d


        function axe_is_regular_double(x, nx, dx, fill_check) result(isregular)
            real(kind=8),intent(in)             :: x(:),dx
            integer,     intent(in)             :: nx
            logical,     intent(in), optional   :: fill_check
            logical                             :: isregular
            real(kind=8)                        :: dx2(nx-1)

            isregular = .false.
            dx2 = x(2:nx) - x(1:nx-1)
            ! Should be, epsilon(1.0d0). However, our own aden grid does not comply
            if ( any( abs(dx2 - dx) > 1e-4 )) return
            if ( present(fill_check) ) then
                if ( fill_check .and. all( abs(x - NF90_FILL_DOUBLE) < epsilon(1.0d0)) ) return
            end if
            isregular = .true.

        end function axe_is_regular_double

        subroutine nccheck(status)
            integer, intent ( in ) :: status

            if (status /= nf90_noerr) then
                write(6,*) trim(nf90_strerror(status))
                STOP 1
            end if
        end subroutine nccheck

        subroutine agnc_set_error(msg)
            character(len=*), intent (in) :: msg

            write(6,*) trim(msg)
            STOP 2

        end subroutine agnc_set_error

        function values_in_range_float(ncid, varid, ovalues, fill_value, scale_factor, add_offset, varname) result (inrange)
            integer,                    intent(   in)         :: ncid, varid
            real(kind=4), dimension(:), intent(   in)         :: ovalues
            real(kind=4),               intent(   in)         :: fill_value, scale_factor, add_offset
            character(len=*),           intent(   in)         :: varname
            ! local
            real(kind=4), dimension(:), allocatable           :: values
            integer                                           :: xtype
            logical                                           :: inrange
            real(kind=4)                                      :: bmax, bmin, thmin, thmax, acmin, acmax
            inrange = .true.
            write (6, *), "Checking the values of ", varname
            allocate(values(size(ovalues)))
            values = ovalues
            ! Get variable meta information
            call nccheck ( nf90_inquire_variable(ncid, varid, xtype = xtype))
            ! Compute theoretic min and max values given the type and offset / scale
            select case (xtype)
            case (NF90_BYTE)
                bmin = -127.
                bmax = 127.
            case (NF90_SHORT)
                bmin = -32767.
                bmax =  32767.
            case (NF90_INT)
                bmin = -2147483648.
                bmax = 2147483647.
            case default
                write(6,*) "Cannot check the values of this variable as the netcdf type is not byte, short or int"
                deallocate(values)
                inrange = .false.
                return
            end select
            thmin = bmin * scale_factor + add_offset
            thmax = bmax * scale_factor + add_offset
            ! compute original values
            where (abs(values - fill_value) > epsilon(1.)) values = values * scale_factor + add_offset
            acmin = minval(values, 1, abs(values - fill_value) > epsilon(1.))
            acmax = maxval(values, 1, abs(values - fill_value) > epsilon(1.))
            deallocate(values)
            if (acmin < thmin .or. acmax > thmax) then
                write (6, *) "The values are outside the boundaries provided in nctablemd.f90"
                print *, "Theoretical min / max: ", thmin, "/", thmax
                print *, "Actual      min / max: ", acmin, "/", acmax
                inrange = .false.
            end if
        end function values_in_range_float

        function values_in_range_double(ncid, varid, ovalues, fill_value, scale_factor, add_offset, varname) result (inrange)
            integer,                    intent(   in)         :: ncid, varid
            real(kind=8), dimension(:), intent(   in)         :: ovalues
            real(kind=8),               intent(   in)         :: fill_value, scale_factor, add_offset
            character(len=*),           intent(   in)         :: varname
            ! local
            real(kind=8), dimension(:), allocatable           :: values
            integer                                           :: xtype
            logical                                           :: inrange
            real(kind=8)                                      :: bmax, bmin, thmin, thmax, acmin, acmax
            inrange = .true.
            write (6, *), "Checking the values of ", varname
            allocate(values(size(ovalues)))
            values = ovalues
            ! Get variable meta information
            call nccheck ( nf90_inquire_variable(ncid, varid, xtype = xtype))
            ! Compute theoretic min and max values given the type and offset / scale
            select case (xtype)
            case (NF90_BYTE)
                bmin = -127.
                bmax = 127.
            case (NF90_SHORT)
                bmin = -32767.
                bmax =  32767.
            case (NF90_INT)
                bmin = -2147483648.
                bmax = 2147483647.
            case default
                write(6,*) "Cannot check the values of this variable as the netcdf type is not byte, short or int"
                deallocate(values)
                inrange = .false.
                return
            end select
            thmin = bmin * scale_factor + add_offset
            thmax = bmax * scale_factor + add_offset
            ! compute original values
            where (abs(values - fill_value) > epsilon(1.)) values = values * scale_factor + add_offset
            acmin = minval(values, 1, abs(values - fill_value) > epsilon(1.))
            acmax = maxval(values, 1, abs(values - fill_value) > epsilon(1.))
            deallocate(values)
            if (acmin < thmin .or. acmax > thmax) then
                write (6, *) "The values are outside the boundaries provided in nctablemd.f90"
                print *, "Theoretical min / max: ", thmin, "/", thmax
                print *, "Actual      min / max: ", acmin, "/", acmax
                inrange = .false.
            end if
        end function values_in_range_double

        function timeindex(tarr, t) result (ti)
            integer, dimension(:), intent( in)          :: tarr
            integer,               intent( in)          :: t
            integer                                     :: ti
            if ( t >= minval(tarr) .and. t <= maxval(tarr) ) then
                ! somewhere in the existing time axe  .and. (t <= max(tarr,1))
                ti = minloc( abs(tarr-t), 1)
            else
                if ( t < minval(tarr) ) then
                    ! smaller then the existing time axe. This should not be allowed
                    ti = -1
                else
                    ! larger. Extend time axe by one
                    ti = 0
                end if
            end if
        end function timeindex

        subroutine agnc_collect_spcmeta(ncid, spcgrid, recordaxe, varid_scale_energy, varid_energy,&
                                      varid_theta, varid_spr, ao_theta, sf_theta, fv_theta,&
                                      ao_spr, sf_spr, fv_spr)
            integer,                           intent( in)      :: ncid
            integer                                             :: varid_scale_energy, &
                                                                   varid_energy, varid_theta, &
                                                                   varid_spr
            real(kind=4)                                        :: sf_theta, sf_spr, &
                                                                   ao_theta, ao_spr, &
                                                                   fv_theta, fv_spr
            type(spcgrid_type)                                  :: spcgrid
            type (recordaxe_type)                                :: recordaxe
            call agnc_get_spcgrid(ncid, spcgrid)
            call agnc_get_recordaxe(ncid, recordaxe)
            call agnc_get_varid_by_name(ncid, "scale_energy_1d", varid_scale_energy)
            call agnc_get_varid_by_name(ncid, "energy_1d", varid_energy)
            call agnc_get_varid_by_name(ncid, "theta_1d", varid_theta)
            call agnc_get_varid_by_name(ncid, "spread_1d", varid_spr)
            call get_scalies_float(ncid, varid_theta, ao_theta, sf_theta, fv_theta)
            call get_scalies_float(ncid, varid_spr, ao_spr, sf_spr, fv_spr)
        end subroutine agnc_collect_spcmeta

        subroutine agnc_collect_spcmeta_2d(ncid, spcgrid, recordaxe, varid_density, varid_scale_density, &
                                         fv_density)
            integer,                           intent( in)      :: ncid
            integer,                           intent(out)      :: varid_scale_density, &
                                                                   varid_density
            real(kind=4)                                        :: ao_density, sf_density
            real(kind=4),                      intent(out)      :: fv_density
            type(spcgrid_type),                intent(out)      :: spcgrid
            type (recordaxe_type),              intent(out)      :: recordaxe

            call agnc_get_spcgrid(ncid, spcgrid)
            call agnc_get_recordaxe(ncid, recordaxe)
            call agnc_get_varid_by_name(ncid, "density", varid_density)
            call agnc_get_varid_by_name(ncid, "scale_density", varid_scale_density)
            call get_scalies_float(ncid, varid_density, ao_density, sf_density, fv_density)
        end subroutine agnc_collect_spcmeta_2d

        function datevec_from_epoch( t_in ) result (datevec)
            integer, intent( in)    :: t_in
            integer                 :: t, year, month ,day, s
            integer                 :: datevec(6)
            ! number of days per month
            integer, save           :: ndpm(12)
            data    ndpm / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
            t = t_in
            datevec = (/0, 0, 0, 0, 0, 0/)
            datevec(6) = mod(t, 60)
            t=t-datevec(6)
            datevec(5) = mod(t,3600)/60
            t = t - datevec(5)*60
            datevec(4) = mod(t,86400)/3600
            t = t-datevec(4)*3600
            do year = 1970, 2038
                if ( (mod(year,400) == 0) .or. ( (mod(year,4) == 0) .and. (mod(year,100) /= 0 )) ) then
                    s = 366*86400
                else
                    s = 365*86400
                end if
                if (t < s) then
                    exit
                else
                    t = t-s
                end if
            end do
            datevec(1) = year
            do month = 1, 12
                if ( ((mod(year,400) == 0) .or. ( (mod(year,4) == 0) &
                      .and. (mod(year,100) /= 0) )) .and. month == 2 ) then
                    s = 29 * 86400
                else
                    s = ndpm(month) * 86400
                end if
                if (t < s) then
                    exit
                else
                    t = t-s
                end if
            end do
            datevec(2) = month
            datevec(3) = t/86400+1
        end function datevec_from_epoch

        function datevec_from_twoint(date, time) result (datevec)
            integer, intent( in)    :: date, time
            integer                 :: datevec(6)
            datevec(1)    = date / 10000
            datevec(2)    = mod(date,10000) / 100
            datevec(3)    = mod(date,100)
            datevec(4)    = time / 10000
            datevec(5)    = mod(time,10000) / 100
            datevec(6)    = mod(time,100)
        end function datevec_from_twoint

        integer function seconds_since_epoch_twoint ( date, time ) result (seconds_since_epoch)
            integer, intent( in)    :: date, time
            seconds_since_epoch = seconds_since_epoch_datevec( datevec_from_twoint(date, time) )
        end function seconds_since_epoch_twoint

        integer function seconds_since_epoch_datevec ( datevec ) result (seconds_since_epoch)
        !     convert date time in (int yyyymmdd, hhmmss) format to seconds since 1-jan-1970
            integer                 :: datevec(6)
            integer                 :: year, month
            ! number of days per month
            integer, save           :: ndpm(12)
            data    ndpm / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
            if (datevec(2) > 12) then
                datevec(1) = datevec(1) + datevec(2)/12
                datevec(2) = mod(datevec(2),12)
            end if
            seconds_since_epoch = 0
            do year = 1970, datevec(1)-1
                if ( (mod(year,400) == 0) .or. ( (mod(year,4) == 0) .and. (mod(year,100) /= 0 )) ) then
                    seconds_since_epoch = seconds_since_epoch + 366*86400
                else
                    seconds_since_epoch = seconds_since_epoch + 365*86400
                end if
            end do
            ! year = nyears
            do month = 1, datevec(2) - 1
                if ( ((mod(year,400) == 0) .or. ( (mod(year,4) == 0) &
                      .and. (mod(year,100) /= 0) )) .and. month == 2 ) then
                    seconds_since_epoch = seconds_since_epoch + (ndpm(month) + 1) * 86400
                else
                    seconds_since_epoch = seconds_since_epoch + ndpm(month) * 86400
                end if
            end do
            seconds_since_epoch = seconds_since_epoch + (datevec(3)-1)*86400 + datevec(4)*3600 &
                    + datevec(5) * 60 + datevec(6)
        end function seconds_since_epoch_datevec

    end module agioncmd
