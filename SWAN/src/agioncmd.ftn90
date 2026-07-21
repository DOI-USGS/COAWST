!/ ------------------------------------------------------------------- /
!/
!/
!/    version | date        | description
!/    1.0     | 01-nov-2006 | netcdf support in wavewatch III 2.22
!/    1.1     | 03-dec-2009 | netcdf support in wavewatch III 3.04
!/    1.2     | 01-sep-2012 | netcdf support in SWAN 40.91
!/    1.3     | 01-may-2013 | overhauled to simplify maintenance
!/    1.4     | 11-apr-2014 | bug fixes
!/                          | reduced memory consumption
!/
!/
!  1. purpose :
!
!     agioncmd is a module that provides methods to write cf1.5 compliant
!     netcdf files containing integrated parameters and spectra.
!     Furthermore, some methods are available to read data from cf-1.5 and
!     coads compliant files.
!
!     For more information:
!       http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.5
!
!  2. method :
!
!  3.a.types
!     --name-------------type-------description-----------------------
!       mapgrid                     data on a lon / lat grid
!         %longitude     r.a        x-coordinates
!         %latitude      r.a        y-coordinates
!         %nx            int        number of columns
!         %ny            int        number of rows
!         %dx            real       horizontal cellsize
!         %dy            real       vertical cellsize
!         %lunit         char       likely unit (meter, degrees)
!         %alpc          real       angle of the grid
!         %mdc           bool       multi dimensional coordinates
!
!       recordaxe                   derived type with fields
!         %content       int.a      seconds since 1-jan-1970 / run
!         %delta         int        timestep (s) / 1
!         %ncontent      int        number of times / runs
!         %nstatm        bool       non-stationary mode (time /run switch)
!         %dimid         int        netcdf dimension id
!         %varid         int        netcdf variable id
!
!       pointgrid                   data along a point dimension
!         %points        int        1:npoints
!         %longitude     r.a        x coordinates points
!         %latitude      r.a        y coordinates points
!         %lunit         char       likely unit (meter, degrees)
!         %npoints
!
!       spcgrid                     spectral grid
!         %frequency     r.a        frequency axe
!         %direction     r.a.       directional axe
!         %nfreq         int        number of frequencies
!         %ndir          int        number of directions
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
!
!  4.procedures
!       ---- create new netcdf files ----
!       create_ncfile(ncfile, ncid)
!       agnc_define_mapgrid
!       agnc_define_pntgrid
!       agnc_define_spcgrid
!       agnc_define_mapvariable
!       agnc_define_pntvariables
!       agnc_set_mapgrid(ncid, mapgrid)
!       agnc_set_pntgrid(ncid, pntgrid)
!       agnc_set_spcgrid(ncid, spcgrid)
!       agnc_set_recordaxe   (ncid, recordaxe)
!       agnc_add_variable_attributes
!       agnc_add_coordinates_attribute
!       ------------- or ---------------
!       open_ncfile
!       --------- get grids -------------
!       agnc_get_mapgrid
!       agnc_get_pntgrid
!       agnc_get_spcgrid
!       agnc_get_recordaxe
!       ------ insert/append data ------
!       agnc_add_mapdata
!       agnc_add_pntdata
!       agnc_add_spcdata
!       agnc_add_density
!       ------------- utils -------------
!       coordinate_names_units
!       agnc_get_varid_by_name
!       close_ncfile
!       nccheck
!       axe_is_regular
!       compute_1d_spectra
!
!  5. used by :
!
!     ww3_outnc.ftn (wavewatch)
!     swoutp (swan)
!
!  6. exit codes
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
        real, parameter                              :: AGIONCMD_VERSION = 1.5
        real, parameter                              :: AGNC_DUMMY = NF90_FILL_FLOAT
        integer(kind=1), parameter                   :: AGNC_FILL_BYTE = -2**7
        integer(kind=2), parameter                   :: AGNC_FILL_SHORT= -2**15

        type mapgrid_type
            real, dimension(:,:), allocatable        :: longitude, latitude
            character(len=10)                        :: lunit = 'degrees'
            integer                                  :: nx = 0, ny = 0, &
                                                        lon_dimid, lat_dimid, &
                                                        lon_varid, lat_varid
            real                                     :: dx, dy, alpc = 0.
            logical                                  :: mdc = .false.

        end type mapgrid_type

        type pntgrid_type
            real, dimension(:), allocatable          :: longitude, latitude
            integer, dimension(:), allocatable       :: ips
            character(len=10)                        :: lunit = 'degrees'
            integer                                  :: pnt_dimid, &
                                                        lon_varid, lat_varid, &
                                                        npoints = 0, ips_varid = 0, &
                                                        xdimlen = 0, ydimlen = 0
        end type pntgrid_type

        type spcgrid_type
            real, dimension(:), allocatable          :: frequency, direction
            logical                                  :: relative
            integer                                  :: frq_dimid, frq_varid, &
                                                        dir_dimid, dir_varid, &
                                                        nfreq = 0, ndir = 0
        end type spcgrid_type

        type recordaxe_type
            integer*8, dimension(:), allocatable     :: content
            integer                                  :: delta
            integer                                  :: ncontent = 0, dimid, varid
            logical                                  :: nstatm = .true.
        end type recordaxe_type

        ! Given the dimension size and/or type, select proper procedure/function
        interface agnc_get_griddef
            module procedure agnc_get_griddef_spc, agnc_get_griddef_map, &
                             agnc_get_griddef_pnt
        end interface agnc_get_griddef
        interface axe_is_regular
            module procedure axe_is_regular_integer, &
                             axe_is_regular_integer64, &
                             axe_is_regular_float, &
                             axe_is_regular_float_2d, &
                             axe_is_regular_double
        end interface axe_is_regular
        interface agnc_add_mapdata
            module procedure agnc_add_mapdata_float, agnc_add_mapdata_double
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
        interface timeindex
            module procedure timeindex32, timeindex64
        end interface timeindex
      contains
!
! ---- create new netcdf files ----
!       name
!       agnc_define_mapgrid
!       agnc_define_pntgrid
!       agnc_define_spcgrid
!       agnc_define_mapvariable
!       agnc_define_pntvariables
!       agnc_set_mapgrid(ncid, mapgrid)
!       agnc_set_pntgrid(ncid, pntgrid)
!       agnc_set_spcgrid(ncid, spcgrid)
!       agnc_set_recordaxe   (ncid, recordaxe)

        subroutine create_ncfile( ncfile, ncid, recordaxe, mapgrid, pntgrid, spcgrid, Escale, nautical )
            character(len=*),     intent(in)                    :: ncfile
            integer,              intent(out)                   :: ncid
            type(recordaxe_type), intent(inout)                 :: recordaxe
            type(mapgrid_type),   intent(inout), optional       :: mapgrid
            type(spcgrid_type),   intent(inout), optional       :: spcgrid
            type(pntgrid_type),   intent(inout), optional       :: pntgrid
            logical              ,intent(   in), optional       :: Escale
            logical              ,intent(   in), optional       :: nautical
            character(len=100)                                  :: history, dconv
            logical                                             :: do_scale

            dconv = 'nautical'
            if ( present(nautical) ) then
                if ( .not.nautical ) then
                    dconv = 'cartesian'
                    ! nctablemd::nautical defaults to true
                    call set_nctable_convention_nautical(.false.)
                end if
            end if

            do_scale = .true.
            if ( present(Escale) ) then
                do_scale = Escale
            end if

            write(history,'(a,f3.1)') 'Created with agioncmd version ', AGIONCMD_VERSION
            call nccheck ( nf90_create( ncfile, NF90_NOCLOBBER + NF90_64BIT_OFFSET, ncid) )
            if ( recordaxe%nstatm ) then
                call nccheck ( nf90_def_dim( ncid, 'time', NF90_UNLIMITED, recordaxe%dimid ) );
                call nccheck ( nf90_def_var( ncid, 'time', NF90_INT, recordaxe%dimid, recordaxe%varid ) )
                call nccheck ( nf90_put_att( ncid, recordaxe%varid, 'units', 'seconds since 1970-01-01') )
                call nccheck ( nf90_put_att( ncid, recordaxe%varid, 'calendar', 'gregorian') )
                call nccheck ( nf90_put_att( ncid, recordaxe%varid, 'standard_name', 'time') )
                call nccheck ( nf90_put_att( ncid, recordaxe%varid, 'long_name', 'time') )
            else
                call nccheck ( nf90_def_dim( ncid, 'run', NF90_UNLIMITED, recordaxe%dimid ) );
                call nccheck ( nf90_def_var( ncid, 'run', NF90_INT, recordaxe%dimid, recordaxe%varid ) )
                call nccheck ( nf90_put_att( ncid, recordaxe%varid, 'long_name', 'run number') )
                call nccheck ( nf90_put_att( ncid, recordaxe%varid, 'units', '1' ) )
            end if

            ! global attributes
            call nccheck ( nf90_put_att( ncid, NF90_GLOBAL, 'Conventions', 'CF-1.5') )
            call nccheck ( nf90_put_att( ncid, NF90_GLOBAL, 'History', history) )
            call nccheck ( nf90_put_att( ncid, NF90_GLOBAL, 'Directional_convention', dconv) )

            if ( present(mapgrid) ) then
                call agnc_define_mapgrid(ncid, mapgrid)
                if ( present(spcgrid) ) call agnc_define_spcgrid(ncid, spcgrid, do_scale, dconv, mapgrid=mapgrid)
            else
                if ( present(pntgrid) ) then
                    call agnc_define_pntgrid(ncid, pntgrid)
                    if ( present(spcgrid) ) call agnc_define_spcgrid(ncid, spcgrid, do_scale, dconv, pntgrid=pntgrid)
                else
                    call agnc_set_error('Provide a point grid or a map grid to create_ncfile')
                end if
            end if

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

            if ( mapgrid%mdc ) then
                call nccheck ( nf90_def_dim( ncid, 'xc', mapgrid%nx, mapgrid%lon_dimid ) );
                call nccheck ( nf90_def_dim( ncid, 'yc', mapgrid%ny, mapgrid%lat_dimid ) );

                call nccheck ( nf90_def_var( ncid, xname, NF90_FLOAT, &
                                            (/ mapgrid%lon_dimid, mapgrid%lat_dimid /), &
                                            mapgrid%lon_varid ) )
                call nccheck ( nf90_def_var( ncid, yname, NF90_FLOAT, &
                                            (/ mapgrid%lon_dimid, mapgrid%lat_dimid /), &
                                            mapgrid%lat_varid ) )
                call nccheck ( nf90_put_att( ncid, mapgrid%lon_varid, '_FillValue', NF90_FILL_FLOAT) )
                call nccheck ( nf90_put_att( ncid, mapgrid%lat_varid, '_FillValue', NF90_FILL_FLOAT) )
            else
                call nccheck ( nf90_def_dim( ncid, xname, mapgrid%nx, mapgrid%lon_dimid ) );
                call nccheck ( nf90_def_var( ncid, xname, NF90_FLOAT, mapgrid%lon_dimid, mapgrid%lon_varid ) );
                call nccheck ( nf90_def_dim( ncid, yname, mapgrid%ny, mapgrid%lat_dimid ) );
                call nccheck ( nf90_def_var( ncid, yname, NF90_FLOAT, mapgrid%lat_dimid, mapgrid%lat_varid ) );
            end if

            call nccheck ( nf90_put_att( ncid, mapgrid%lon_varid, 'units', xunit) )
            call nccheck ( nf90_put_att( ncid, mapgrid%lon_varid, 'long_name', xname) )
            if ( trim(xname) == 'longitude' ) then
                call nccheck ( nf90_put_att( ncid, mapgrid%lon_varid, 'standard_name', 'longitude') )
            end if

            call nccheck ( nf90_put_att( ncid, mapgrid%lat_varid, 'units', yunit) )
            call nccheck ( nf90_put_att( ncid, mapgrid%lat_varid, 'long_name', yname) )
            if ( trim(yname) == 'latitude' ) then
                call nccheck ( nf90_put_att( ncid, mapgrid%lat_varid, 'standard_name', 'latitude') )
            end if

        end subroutine agnc_define_mapgrid


        subroutine agnc_define_spcgrid( ncid, spcgrid, do_scale, dconv, mapgrid, pntgrid )
            integer,              intent(in   )                 :: ncid
            type(spcgrid_type),   intent(inout)                 :: spcgrid
            logical,              intent( in)                   :: do_scale
            character(len=100),   intent( in)                   :: dconv
            type(mapgrid_type),   intent(inout), optional       :: mapgrid
            type(pntgrid_type),   intent(inout), optional       :: pntgrid

            integer                                             :: evarid, ra_dimid, &
                                                                   nf90_var_type
            logical                                             :: nautical

            nautical = .false.
            if ( dconv == 'nautical') nautical = .true.

            if ( do_scale ) then
                nf90_var_type  = NF90_SHORT
            else
                nf90_var_type  = NF90_FLOAT
            end if

            !
            ! frequency
            !
            call nccheck ( nf90_def_dim( ncid, 'frequency', spcgrid%nfreq, spcgrid%frq_dimid ) );
            call nccheck ( nf90_def_var( ncid, 'frequency', NF90_FLOAT, spcgrid%frq_dimid, spcgrid%frq_varid ) );
            call nccheck ( nf90_put_att( ncid, spcgrid%frq_varid, 'units', 's-1') )
            call nccheck ( nf90_put_att( ncid, spcgrid%frq_varid, 'standard_name', 'wave_frequency') )

            call nccheck ( nf90_put_att( ncid, spcgrid%frq_varid, 'flow',  spcgrid%frequency(1) ) )
            call nccheck ( nf90_put_att( ncid, spcgrid%frq_varid, 'fhigh', spcgrid%frequency(spcgrid%nfreq) ) )
            call nccheck ( nf90_put_att( ncid, spcgrid%frq_varid, 'msc',   spcgrid%nfreq-1 ) )

            call agnc_get_recordaxe_ids( ncid, ra_dimid )

            if ( spcgrid%ndir > 0 ) then
                !
                ! 2D direction axe
                !
                call nccheck ( nf90_def_dim( ncid, 'direction', spcgrid%ndir, spcgrid%dir_dimid ) )
                call nccheck ( nf90_def_var( ncid, 'direction', NF90_FLOAT, spcgrid%dir_dimid, spcgrid%dir_varid ) );
                call nccheck ( nf90_put_att( ncid, spcgrid%dir_varid, 'units', 'radians') )
                call nccheck ( nf90_put_att( ncid, spcgrid%dir_varid, 'long_name', 'direction') )
                if ( nautical ) then
                    call nccheck ( nf90_put_att( ncid, spcgrid%dir_varid, 'standard_name', 'sea_surface_wave_from_direction') )
                else
                    call nccheck ( nf90_put_att( ncid, spcgrid%dir_varid, 'standard_name', 'sea_surface_wave_to_direction') )
                end if
                call nccheck ( nf90_put_att( ncid, spcgrid%dir_varid, 'mdc',   spcgrid%ndir ) )

                !
                ! 2D density spectrum
                !
                if ( present(mapgrid) ) then
                    call nccheck ( nf90_def_var( ncid, 'density', nf90_var_type, &
                                            (/ spcgrid%dir_dimid, spcgrid%frq_dimid, mapgrid%lon_dimid, mapgrid%lat_dimid, ra_dimid /), evarid ) )
                else
                    call nccheck ( nf90_def_var( ncid, 'density', nf90_var_type, &
                                        (/ spcgrid%dir_dimid, spcgrid%frq_dimid, pntgrid%pnt_dimid, ra_dimid /), evarid ) )
                end if

                if ( do_scale ) then
                    call nccheck ( nf90_put_att( ncid, evarid, '_FillValue',NF90_FILL_SHORT) )
                else
                    call nccheck ( nf90_put_att( ncid, evarid, '_FillValue',NF90_FILL_FLOAT) )
                end if
                call nccheck ( nf90_put_att( ncid, evarid, 'long_name', 'density') )
                call nccheck ( nf90_put_att( ncid, evarid, 'units', 'm2 s rad-1') )
                call nccheck ( nf90_put_att( ncid, evarid, 'standard_name', &
                            'sea_surface_wave_directional_variance_spectral_density') )

                if ( do_scale) then
                    call nccheck ( nf90_put_att( ncid, evarid, 'note', &
                            'multiply with scale_density') )
                end if

                if ( spcgrid%relative ) then
                    call nccheck ( nf90_put_att( ncid, evarid, 'relative_to_current', 'true') )
                else
                    call nccheck ( nf90_put_att( ncid, evarid, 'relative_to_current', 'false') )
                end if

                !
                ! scale_density
                !
                if ( do_scale ) then
                    if ( present(mapgrid) ) then
                        call nccheck ( nf90_def_var( ncid, 'scale_density', NF90_FLOAT, &
                                                (/ mapgrid%lon_dimid, mapgrid%lat_dimid, ra_dimid /), evarid ) )
                    else
                        call nccheck ( nf90_def_var( ncid, 'scale_density', NF90_FLOAT, &
                                                (/ pntgrid%pnt_dimid, ra_dimid /), evarid ) )
                    end if
                    call nccheck ( nf90_put_att( ncid, evarid, 'note', &
                                        'multiply density values with this scale factor') )
                    call nccheck ( nf90_put_att( ncid, evarid, 'units', '1') )
                end if
            else
                if ( present(mapgrid) ) then
                    call nccheck ( nf90_def_var( ncid, 'energy_1d', nf90_var_type, &
                                            (/ spcgrid%frq_dimid, mapgrid%lon_dimid, mapgrid%lat_dimid, ra_dimid /), evarid ) )
                else
                    call nccheck ( nf90_def_var( ncid, 'energy_1d', nf90_var_type, &
                                            (/ spcgrid%frq_dimid, pntgrid%pnt_dimid, ra_dimid /), evarid ) )
                end if
                if ( do_scale ) then
                    call nccheck ( nf90_put_att( ncid, evarid, '_FillValue',NF90_FILL_SHORT) )
                else
                    call nccheck ( nf90_put_att( ncid, evarid, '_FillValue',NF90_FILL_FLOAT) )
                end if
                call nccheck ( nf90_put_att( ncid, evarid, 'long_name', 'energy') )
                call nccheck ( nf90_put_att( ncid, evarid, 'units', 'm2 s') )
                call nccheck ( nf90_put_att( ncid, evarid, 'standard_name', 'sea_surface_wave_variance_spectral_density') )

                if ( spcgrid%relative ) then
                    call nccheck ( nf90_put_att( ncid, evarid, 'relative_to_current', 'true') )
                else
                    call nccheck ( nf90_put_att( ncid, evarid, 'relative_to_current', 'false') )
                end if
                call agnc_add_coordinates_attribute(ncid, evarid)

                !
                ! scale_density
                !
                if ( do_scale ) then
                    call nccheck ( nf90_put_att( ncid, evarid, 'note', &
                            'multiply with scale_energy_1d') )

                    if ( present(mapgrid) ) then
                        call nccheck ( nf90_def_var( ncid, 'scale_energy_1d', NF90_FLOAT, &
                                                (/ mapgrid%lon_dimid, mapgrid%lat_dimid, ra_dimid /), evarid ) )
                    else
                        call nccheck ( nf90_def_var( ncid, 'scale_energy_1d', NF90_FLOAT, &
                                                (/ pntgrid%pnt_dimid, ra_dimid /), evarid ) )
                    end if
                    call nccheck ( nf90_put_att( ncid, evarid, 'long_name', &
                            'Multiply energy_1d values with this scale factor') )
                    call nccheck ( nf90_put_att( ncid, evarid, 'units', &
                            '1') )
                end if

                !
                ! theta_1d
                !
                if ( present(mapgrid) ) then
                    call nccheck ( nf90_def_var( ncid, 'theta_1d', NF90_BYTE, &
                                                (/ spcgrid%frq_dimid, mapgrid%lon_dimid, mapgrid%lat_dimid, ra_dimid /), evarid ) )
                else
                    call nccheck ( nf90_def_var( ncid, 'theta_1d', NF90_BYTE, &
                                                (/ spcgrid%frq_dimid, pntgrid%pnt_dimid, ra_dimid /), evarid ) )

                end if
                call nccheck ( nf90_put_att( ncid, evarid, 'units', 'degree') )
                call nccheck ( nf90_put_att( ncid, evarid, 'long_name', 'principal wave direction') )
                call nccheck ( nf90_put_att( ncid, evarid, '_FillValue', AGNC_FILL_BYTE) )
                call nccheck ( nf90_put_att( ncid, evarid, 'scale_factor', 360. / (2**8 -2)) )
                call nccheck ( nf90_put_att( ncid, evarid, 'add_offset', 360. / 2.) )
                call agnc_add_coordinates_attribute(ncid, evarid)

                !
                ! spread_1d
                !
                if ( present(mapgrid) ) then
                    call nccheck ( nf90_def_var( ncid, 'spread_1d', NF90_BYTE, &
                                                (/ spcgrid%frq_dimid, mapgrid%lon_dimid, mapgrid%lat_dimid, ra_dimid /), evarid ) )
                else
                    call nccheck ( nf90_def_var( ncid, 'spread_1d', NF90_BYTE, &
                                                (/ spcgrid%frq_dimid, pntgrid%pnt_dimid, ra_dimid /), evarid ) )
                end if
                call nccheck ( nf90_put_att( ncid, evarid, 'units', 'degree') )
                call nccheck ( nf90_put_att( ncid, evarid, 'long_name', &
                    'Longuet-Higgins short-crestedness parameter (s in cos(theta/2)^2s)') )
                call nccheck ( nf90_put_att( ncid, evarid, '_FillValue', AGNC_FILL_BYTE) )
                call agnc_add_coordinates_attribute(ncid, evarid)

            end if

        end subroutine agnc_define_spcgrid

        subroutine agnc_define_pntgrid( ncid, pntgrid )
            integer,            intent(in   )       :: ncid
            type (pntgrid_type),intent(inout)       :: pntgrid
            character(len=20)                       :: xname, xunit
            character(len=20)                       :: yname, yunit
            call coordinate_names_units( pntgrid%lunit, xname, xunit, yname, yunit )
            call nccheck ( nf90_def_dim( ncid, 'points', pntgrid%npoints, pntgrid%pnt_dimid ) )

            call nccheck ( nf90_def_var( ncid, xname, NF90_FLOAT, pntgrid%pnt_dimid, pntgrid%lon_varid ) );
            call nccheck ( nf90_put_att( ncid, pntgrid%lon_varid, 'units', xunit) )
            call nccheck ( nf90_put_att( ncid, pntgrid%lon_varid, 'long_name', xname) )

            if ( trim(xname) == 'longitude' ) then
                call nccheck ( nf90_put_att( ncid, pntgrid%lon_varid, 'standard_name', 'longitude') )
            end if

            call nccheck ( nf90_def_var( ncid, yname, NF90_FLOAT, pntgrid%pnt_dimid, pntgrid%lat_varid ) );
            call nccheck ( nf90_put_att( ncid, pntgrid%lat_varid, 'units', yunit) )
            call nccheck ( nf90_put_att( ncid, pntgrid%lat_varid, 'long_name', yname) )
            if ( trim(yname) == 'latitude' ) then
                call nccheck ( nf90_put_att( ncid, pntgrid%lat_varid, 'standard_name', 'latitude') )
            end if

            ! http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.5/ch08s02.html
            ! states that compression by gathering consists of having a pointlist with one dimensional
            ! indices in the two dimensional grid (lon, lat). In the example lon and lat are arrays sized
            ! nx and ny. This code delibertly breaks that rule because it is also used to write unstructured
            ! data and along a user defined pointlist. It makes more sense to have a longitude / latitude with the
            ! size of the pointlist.
            if ( allocated(pntgrid%ips) ) then
                call nccheck ( nf90_put_att( ncid, pntgrid%lon_varid, 'dimlen', pntgrid%xdimlen) )
                call nccheck ( nf90_put_att( ncid, pntgrid%lat_varid, 'dimlen', pntgrid%ydimlen) )
                call nccheck ( nf90_def_var( ncid, 'ips', NF90_INT, pntgrid%pnt_dimid, pntgrid%ips_varid ) );
                call nccheck ( nf90_put_att( ncid, pntgrid%ips_varid, 'units', '1') )
                call nccheck ( nf90_put_att( ncid, pntgrid%ips_varid, 'computed_as', '(i - 1) * dimlen1 + j') )
                call nccheck ( nf90_put_att( ncid, pntgrid%ips_varid, 'description', 'i is index in first spatial dimension, dimlen1 its length') )
                call nccheck ( nf90_put_att( ncid, pntgrid%ips_varid, 'long_name', 'one-dimensional index') )
            end if

        end subroutine agnc_define_pntgrid

        subroutine agnc_define_mapvariable (ncid, varname)
            integer,                intent( in)         :: ncid
            character(len=*),       intent( in)         :: varname
            type(nctable_record)                        :: trecord
            integer                                     :: ra_dimid, lon_dimid, lat_dimid,&
                                                           varid

            call get_nctable_record(varname, trecord)
            call agnc_get_recordaxe_ids(ncid, ra_dimid)
            call agnc_get_ll_dimids(ncid, lon_dimid, lat_dimid)

            call nccheck ( nf90_def_var( ncid, trecord%name , trecord%nctype, &
                                            (/ lon_dimid, lat_dimid, ra_dimid /), varid ) )

            call agnc_add_variable_attributes( ncid, varid, trecord )

        end subroutine agnc_define_mapvariable

        subroutine agnc_define_pntvariable (ncid, varname)
            integer,                intent( in)         :: ncid
            character(len=*),       intent( in)         :: varname
            type(nctable_record)                        :: trecord
            integer                                     :: pnt_dimid, ra_dimid, varid

            call get_nctable_record(varname, trecord)

            call agnc_get_recordaxe_ids(ncid, ra_dimid)

            ! both unstructured as spectral grid have the "point" dimension
            call nccheck ( nf90_inq_dimid(ncid, 'points', pnt_dimid) )

            call nccheck ( nf90_def_var( ncid, trecord%name , trecord%nctype, &
                                            (/ pnt_dimid, ra_dimid /), varid ) )

            call agnc_add_variable_attributes( ncid, varid, trecord )

        end subroutine agnc_define_pntvariable

        subroutine agnc_add_variable_attributes( ncid, varid, trecord )
            integer,                intent(in) :: ncid, varid
            type(nctable_record),   intent(in) :: trecord
            real(kind=4)                       :: add_offset, scale_factor, rng
            integer                            :: tmpid
            character(len=nf90_max_name)       :: lonname

            call nccheck ( nf90_put_att( ncid, varid, 'units', trecord%units) )
            if ( trim(trecord%standard_name) /= 'none' ) &
                call nccheck ( nf90_put_att( ncid, varid, 'standard_name', trecord%standard_name) )

            call nccheck ( nf90_put_att( ncid, varid, 'long_name', trecord%long_name) )

            call agnc_add_coordinates_attribute(ncid, varid)

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
                    scale_factor = rng/(2**8-2)
                case (NF90_SHORT)
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
            case (NF90_FLOAT)
                call nccheck ( nf90_put_att( ncid, varid, '_FillValue', NF90_FILL_FLOAT) )
            end select

        end subroutine agnc_add_variable_attributes

        subroutine agnc_add_coordinates_attribute(ncid, varid)
            integer,                intent(in) :: ncid, varid
            integer                            :: tmpid
            ! CF-1.2+ prescibes a "coordinates" attribute on multi-dimension fields

            if ( nf90_inq_dimid(ncid, 'xc', tmpid) == nf90_noerr ) then
                ! multi-dimension field because of the hard coded dimension name. The coordinates
                ! are either in the longitude/latitude or x,y variables
                if ( nf90_inq_varid(ncid, 'x', tmpid) == nf90_noerr ) then
                    call nccheck ( nf90_put_att( ncid, varid, 'coordinates', 'x y') )
                else
                    call nccheck ( nf90_put_att( ncid, varid, 'coordinates', 'longitude latitude') )
                end if
            end if

        end subroutine agnc_add_coordinates_attribute

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
            integer                                          :: oldmode
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
!       name
!       agnc_set_recordaxe
!       agnc_set_mapgrid
!       agnc_set_pntgrid
!       agnc_set_spcgrid (ncid, spcgrid)                  done

        subroutine agnc_set_recordaxe (ncid, recordaxe, monthly)
        ! When the optional monthly is set to true, the time axis in the netcdf file
        ! is predeclared for the month of trecordaxe%content(1)
            integer,               intent (   in)                 :: ncid
            type (recordaxe_type), intent (inout)                 :: recordaxe
            logical, optional   ,  intent (   in)                 :: monthly
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

        subroutine agnc_set_pntgrid (ncid, pntgrid)
            integer,                intent (   in)                  :: ncid
            type (pntgrid_type),    intent (inout)                  :: pntgrid
            character(len=nf90_max_name)                            :: name

            if ( pntgrid%lon_varid < 1 ) then
                call agnc_get_griddef(ncid, pntgrid)
            end if
            call nccheck ( nf90_put_var(ncid, pntgrid%lon_varid, pntgrid%longitude) )
            call nccheck ( nf90_put_var(ncid, pntgrid%lat_varid, pntgrid%latitude) )
             if ( pntgrid%ips_varid > 0 ) then
                call nccheck ( nf90_put_var(ncid, pntgrid%ips_varid, pntgrid%ips) )
            end if
        end subroutine agnc_set_pntgrid

        subroutine agnc_set_mapgrid (ncid, mapgrid)
            integer,              intent (   in)                  :: ncid
            type (mapgrid_type),  intent (inout)                  :: mapgrid
            character(len=nf90_max_name)                          :: name

            if ( mapgrid%lon_varid < 1 ) then
                call agnc_get_griddef(ncid, mapgrid)
            end if
            call nccheck ( nf90_put_var(ncid, mapgrid%lon_varid, mapgrid%longitude) )
            if ( mapgrid%mdc ) then
                call nccheck ( nf90_put_var(ncid, mapgrid%lat_varid, mapgrid%latitude) )
            else
                call nccheck ( nf90_put_var(ncid, mapgrid%lat_varid, &
                    reshape(mapgrid%latitude, (/mapgrid%ny/)) ))
            end if
        end subroutine agnc_set_mapgrid


        subroutine agnc_set_spcgrid (ncid, spcgrid)
            integer,              intent (   in)                  :: ncid
            type (spcgrid_type),  intent (inout)                  :: spcgrid
            character(len=nf90_max_name)                          :: name

            if ( spcgrid%frq_varid < 1 ) then
                call agnc_get_griddef(ncid, spcgrid)
            end if
            call nccheck ( nf90_put_var(ncid, spcgrid%frq_varid, spcgrid%frequency) )
            if (spcgrid%ndir > 0) &
                call nccheck ( nf90_put_var(ncid, spcgrid%dir_varid, spcgrid%direction) )

        end subroutine agnc_set_spcgrid

!       --------- get data -------------
!       agnc_get_mapgrid
!       agnc_get_pntgrid
!       agnc_get_spcgrid
!       agnc_get_recordaxe
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
            integer                             :: ndims
            real                                :: xlen, ylen
            real, dimension(:), allocatable     :: lattmp

            call agnc_get_griddef(ncid, mapgrid)

            call nccheck ( nf90_inquire_variable(ncid, mapgrid%lon_varid, ndims=ndims) )
            if ( ndims == 1 ) then
                allocate(mapgrid%longitude(mapgrid%nx, 1))
                allocate(mapgrid%latitude (1,mapgrid%ny))
                allocate(lattmp(mapgrid%ny))

                call nccheck ( nf90_get_var(ncid, mapgrid%lon_varid, mapgrid%longitude) )
                call nccheck ( nf90_get_var(ncid, mapgrid%lat_varid, lattmp) )
                mapgrid%latitude(1,:) = lattmp
                deallocate(lattmp)

                xlen = mapgrid%longitude(mapgrid%nx,1) - mapgrid%longitude(1,1)
                ylen = mapgrid%latitude(1,mapgrid%ny)  - mapgrid%latitude(1,1)

                mapgrid%dx   = xlen / ( mapgrid%nx - 1 )
                mapgrid%dy   = ylen / ( mapgrid%ny - 1 )

                mapgrid%alpc = 0.
            else
                allocate(mapgrid%longitude(mapgrid%nx, mapgrid%ny))
                allocate(mapgrid%latitude (mapgrid%nx, mapgrid%ny))
                call nccheck ( nf90_get_var(ncid, mapgrid%lat_varid, mapgrid%latitude) )
                call nccheck ( nf90_get_var(ncid, mapgrid%lon_varid, mapgrid%longitude) )

                xlen = mapgrid%longitude(mapgrid%nx,1) - mapgrid%longitude(1,1)
                ylen = mapgrid%latitude(mapgrid%nx,1)  - mapgrid%latitude(1,1)

                mapgrid%alpc = atan2(xlen, ylen)
                ! greatcircle??
                mapgrid%dx   = ( xlen**2 + ylen**2  )**0.5 / ( mapgrid%nx - 1 )

                xlen = mapgrid%longitude(1,mapgrid%ny) - mapgrid%longitude(1,1)
                ylen = mapgrid%latitude(1,mapgrid%ny)  - mapgrid%latitude(1,1)
                mapgrid%dy = ( xlen**2 + ylen**2)**0.5 / ( mapgrid%ny - 1 )
            end if

        end subroutine agnc_get_mapgrid

        subroutine agnc_get_pntgrid(ncid, pntgrid)
            integer, intent ( in )            :: ncid
            type (pntgrid_type), intent (out) :: pntgrid

            call agnc_get_griddef(ncid, pntgrid)

            allocate(pntgrid%longitude(pntgrid%npoints))
            call nccheck ( nf90_get_var(ncid, pntgrid%lon_varid, pntgrid%longitude) );

            allocate(pntgrid%latitude(pntgrid%npoints))
            call nccheck ( nf90_get_var(ncid, pntgrid%lat_varid, pntgrid%latitude) )

            if ( pntgrid%ips_varid /= -1 ) then
                allocate(pntgrid%ips(pntgrid%npoints))
                call nccheck ( nf90_get_var(ncid, pntgrid%ips_varid, pntgrid%ips) )
            end if

        end subroutine agnc_get_pntgrid

        subroutine agnc_get_spcgrid(ncid, spcgrid)
            integer, intent ( in )              :: ncid
            type (spcgrid_type), intent (out)   :: spcgrid

            call agnc_get_griddef(ncid, spcgrid)

            if ( spcgrid%frq_varid > 0 ) then
                allocate(spcgrid%frequency(spcgrid%nfreq))
                call nccheck ( nf90_get_var(ncid, spcgrid%frq_varid, spcgrid%frequency) )
                if ( spcgrid%dir_varid > 0 ) then
                    allocate(spcgrid%direction(spcgrid%ndir))
                    call nccheck ( nf90_get_var(ncid, spcgrid%dir_varid, spcgrid%direction) )
                end if
            end if
        end subroutine agnc_get_spcgrid

        subroutine agnc_get_recordaxe(ncid, recordaxe)
            integer, intent ( in )                  :: ncid
            type (recordaxe_type), intent(out)      :: recordaxe
            character(len=nf90_max_name)            :: units, raname, emsg
            integer*8, allocatable, dimension(:)    :: dvalues
            integer                                 :: factor, dvec(6)

            ! find time dimid
            call agnc_get_recordaxe_ids(ncid, recordaxe%dimid, recordaxe%varid, recordaxe%ncontent)

            if ( recordaxe%varid < 0 ) then
                call agnc_set_error( 'Could not find record axe variable')
            end if
            call nccheck ( nf90_get_att(ncid, recordaxe%varid, 'units', units) )
            call nccheck ( nf90_inquire_variable(ncid, recordaxe%varid, raname) )
            allocate(recordaxe%content(recordaxe%ncontent))

            if (raname /= 'run') then

                call agnc_scan_time(units, factor, dvec, emsg)
                if ( factor < 0 ) call agnc_set_error(emsg)

                allocate(dvalues(recordaxe%ncontent))
                call nccheck ( nf90_get_var(ncid, recordaxe%varid, dvalues) )
                recordaxe%content = nint(dvalues * dble(factor) + dble(seconds_since_epoch(dvec)),8)
                deallocate(dvalues)
                recordaxe%nstatm = .true.
            else
                call nccheck ( nf90_get_var(ncid, recordaxe%varid, recordaxe%content) )
                recordaxe%nstatm = .false.
            end if

            if (recordaxe%ncontent > 1) then
                recordaxe%delta = recordaxe%content(2) - recordaxe%content(1)
                if (.not. axe_is_regular(recordaxe%content, recordaxe%ncontent, recordaxe%delta)) then
                    call agnc_set_error( 'record axe does not appear to be regular')
                end if
            else
                ! You'll have set dt yourself when nt is 1
                recordaxe%delta = 0
            end if
        end subroutine agnc_get_recordaxe

        subroutine agnc_get_recordaxe_ids(ncid, dimid, varid, ncontent, raname)
            integer,                                intent( in)        :: ncid
            integer,                                intent(out)        :: dimid
            integer,                      optional, intent(out)        :: varid, ncontent
            character(len=nf90_max_name), optional, intent(out)        :: raname
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

        subroutine agnc_get_ll_dimids(ncid, lon_dimid, lat_dimid)
            integer             , intent( in)                   :: ncid
            integer             , intent(out)                   :: lon_dimid, lat_dimid

            if ( nf90_inq_dimid(ncid, "longitude", lon_dimid) /= nf90_noerr ) then
                if ( nf90_inq_dimid(ncid, "x", lon_dimid) /= nf90_noerr ) then
                    if ( nf90_inq_dimid(ncid, "xc", lon_dimid) /= nf90_noerr ) then
                        call nccheck ( nf90_inq_dimid(ncid, "lon", lon_dimid) )
                    end if
                end if
            end if
            if ( nf90_inq_dimid(ncid, "latitude", lat_dimid) /= nf90_noerr ) then
                if ( nf90_inq_dimid(ncid, "y", lat_dimid) /= nf90_noerr ) then
                    if ( nf90_inq_dimid(ncid, "yc", lat_dimid) /= nf90_noerr ) then
                        call nccheck ( nf90_inq_dimid(ncid, "lat", lat_dimid) )
                    end if
                end if
            end if

        end subroutine agnc_get_ll_dimids

        subroutine agnc_get_ll_varids(ncid, lon_varid, lat_varid)
            integer            , intent(in   )                 :: ncid
            integer            , intent(  out)                 :: lon_varid, lat_varid

            character(len=nf90_max_name)                       :: lonnames(3), latnames(3)
            integer                                            :: n
            lonnames = (/'longitude', &
                         'lon      ', &
                         'x        '/)
            latnames = (/'latitude ', &
                         'lat      ', &
                         'y        '/)

            do n=1,3
                call agnc_get_varid_by_name(ncid, trim(lonnames(n)), lon_varid)
                if ( lon_varid /= -1 ) exit
            end do
            if ( lon_varid == -1 ) call agnc_set_error('Could not find horizontal coordinate axe')

            ! find latitude varid
            do n=1,3
                call agnc_get_varid_by_name(ncid, trim(latnames(n)), lat_varid)
                if ( lat_varid /= -1 ) exit
            end do
            if ( lat_varid == -1 ) call agnc_set_error('Could not find vertical coordinate axe')

        end subroutine agnc_get_ll_varids

        subroutine agnc_get_griddef_map(ncid, mapgrid)
            integer            , intent(in   )                 :: ncid
            type(mapgrid_type) , intent(  out)                 :: mapgrid

            call agnc_get_ll_dimids(ncid, mapgrid%lon_dimid, mapgrid%lat_dimid)
            call nccheck ( nf90_inquire_dimension(ncid, mapgrid%lon_dimid, len=mapgrid%nx) )
            call nccheck ( nf90_inquire_dimension(ncid, mapgrid%lat_dimid, len=mapgrid%ny) )
            call agnc_get_ll_varids(ncid, mapgrid%lon_varid, mapgrid%lat_varid)

        end subroutine agnc_get_griddef_map

        subroutine agnc_get_griddef_pnt(ncid, pntgrid)
            integer            , intent(in   )                 :: ncid
            type(pntgrid_type) , intent(  out)                 :: pntgrid

            call nccheck ( nf90_inq_dimid(ncid, 'points', pntgrid%pnt_dimid) )
            call nccheck ( nf90_inquire_dimension(ncid, pntgrid%pnt_dimid, len=pntgrid%npoints) )

            call agnc_get_ll_varids(ncid, pntgrid%lon_varid, pntgrid%lat_varid)

            if ( nf90_inq_varid(ncid, 'ips', pntgrid%ips_varid) /= nf90_noerr ) then
                pntgrid%ips_varid = -1
            else
                call nccheck ( nf90_get_att(ncid, pntgrid%lon_varid, 'dimlen', pntgrid%xdimlen) )
                call nccheck ( nf90_get_att(ncid, pntgrid%lat_varid, 'dimlen', pntgrid%ydimlen) )
            end if

        end subroutine agnc_get_griddef_pnt

        subroutine agnc_get_griddef_spc(ncid, spcgrid)
            integer            , intent(in   )                 :: ncid
            type(spcgrid_type) , intent(  out)                 :: spcgrid

            ! spectral grid
            call nccheck ( nf90_inq_dimid(ncid, "frequency", spcgrid%frq_dimid) )
            call nccheck ( nf90_inquire_dimension(ncid, spcgrid%frq_dimid, len=spcgrid%nfreq) )
            call nccheck ( nf90_inq_varid(ncid, "frequency", spcgrid%frq_varid) )

            ! direction only available when frequency is set
            if ( nf90_inq_dimid(ncid, "direction", spcgrid%dir_dimid) == nf90_noerr ) then
                call nccheck ( nf90_inquire_dimension(ncid, spcgrid%dir_dimid, len=spcgrid%ndir) )
                call nccheck ( nf90_inq_varid(ncid, "direction", spcgrid%dir_varid) )
            else
                spcgrid%dir_dimid = -1
                spcgrid%dir_varid = -1
                spcgrid%ndir      =  0
            end if

        end subroutine agnc_get_griddef_spc

        subroutine get_scalies_float(ncid, varid, add_offset, scale_factor, fill_value, xtype)
            integer,                            intent( in)            :: ncid, varid
            real(kind=4),                       intent(out)            :: add_offset, scale_factor, fill_value
            integer, optional,                  intent(out)            :: xtype
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

        subroutine get_scalies_double(ncid, varid, add_offset, scale_factor, fill_value, xtype)
            integer,                            intent( in)            :: ncid, varid
            real(kind=8),                       intent(out)            :: add_offset, scale_factor, fill_value
            integer, optional,                  intent(out)            :: xtype
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

!
!       ------ insert/append data ------
!
!       agnc_add_mapdata(ncid, varname, time, values)     done
!       agnc_add_spcdata(ncid, e1, th1, spr1)             done

        subroutine agnc_add_mapdata_float(ncid, varname, ti, values, dummyvalue, skip_range_error)
            ! if the optional dummyvalue is given, the values are searched for this value
            ! and replaced with the _FillValue
            ! ti is the index where the data write should be started
            integer,                        intent( in)         :: ncid, ti
            character(len=*),               intent( in)         :: varname
            real(kind=4), optional,         intent( in)         :: dummyvalue
            logical,      optional                              :: skip_range_error
            real(kind=4), dimension(:)                          :: values
            real(kind=4), dimension(:,:), allocatable           :: vloc
            real(kind=4)                                        :: add_offset, scale_factor, &
                                                                   fill_value
            type(mapgrid_type)                                  :: mapgrid
            integer                                             :: varid, xtype, chunksize, &
                                                                   sx, cx, sy, cy, &
                                                                   lxi, lyi, nf90_stat
            ! integrated parameters, chunks can be fairly large
            chunksize = 256

            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_float(ncid, varid, add_offset, scale_factor, fill_value, xtype)
            call agnc_get_griddef(ncid, mapgrid)

            !
            ! _FillValue
            !
            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            !
            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            !
            if ( xtype == NF90_SHORT .or. xtype == NF90_BYTE .or. xtype == NF90_INT ) then
                where ( abs(values - fill_value) > epsilon(1.) ) &
                    values = nint( (values - add_offset) / scale_factor)
            end if

            !
            ! Chunked write
            !
            allocate(vloc(chunksize, chunksize))

            do sx=1,mapgrid%nx,chunksize
                cx = minval( (/chunksize, mapgrid%nx-sx+1/) )
                do sy=1,mapgrid%ny,chunksize
                    cy = minval( (/chunksize, mapgrid%ny-sy+1/) )

                    vloc = 0.
                    do lxi=1,cx
                        do lyi=1,cy
                            vloc(lxi,lyi) = values((sx + lxi - 1) + (sy + lyi - 2) * mapgrid%nx )
                        end do
                    end do

                    nf90_stat = nf90_put_var(ncid, varid, vloc(1:cx,1:cy), (/sx,sy,ti/), (/cx, cy,1/))

                    !
                    ! If the write failed, check the values against the limits configured in nctablemd.ftn90
                    !
                    if ( nf90_stat /= NF90_NOERR ) then
                        if ( values_in_range_float(ncid, varid, values, &
                                fill_value, scale_factor, add_offset, varname) ) then
                            call agnc_set_error( nf90_strerror(nf90_stat))
                        else
                            if (.not. present(skip_range_error)) call agnc_set_error( nf90_strerror(nf90_stat))
                        end if
                    end if
                end do
            end do
            deallocate(vloc)

        end subroutine agnc_add_mapdata_float

        subroutine agnc_add_mapdata_double(ncid, varname, ti, values, dummyvalue, skip_range_error)
            integer,                        intent( in)         :: ncid, ti
            character(len=*),               intent( in)         :: varname
            real(kind=8), optional,         intent( in)         :: dummyvalue
            logical,      optional                              :: skip_range_error
            real(kind=8), dimension(:)                          :: values
            real(kind=8), dimension(:,:), allocatable           :: vloc
            real(kind=8)                                        :: add_offset, scale_factor, &
                                                                   fill_value
            type(mapgrid_type)                                  :: mapgrid
            integer                                             :: varid, xtype, chunksize, &
                                                                   sx, cx, sy, cy, &
                                                                   lxi, lyi, nf90_stat
            ! integrated parameters, chunks can be fairly large
            chunksize = 128

            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_double(ncid, varid, add_offset, scale_factor, fill_value, xtype)
            call agnc_get_griddef(ncid, mapgrid)

            !
            ! _FillValue
            !
            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            !
            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            !
            if ( xtype == NF90_SHORT .or. xtype == NF90_BYTE .or. xtype == NF90_INT ) then
                where ( abs(values - fill_value) > epsilon(1.) ) &
                    values = nint( (values - add_offset) / scale_factor)
            end if

            !
            ! Chunked write
            !
            allocate(vloc(chunksize, chunksize))

            do sx=1,mapgrid%nx,chunksize
                cx = minval( (/chunksize, mapgrid%nx-sx+1/) )
                do sy=1,mapgrid%ny,chunksize
                    cy = minval( (/chunksize, mapgrid%ny-sy+1/) )

                    vloc = 0.
                    do lxi=1,cx
                        do lyi=1,cy
                            vloc(lxi,lyi) = values((sx + lxi - 1) + (sy + lyi - 2) * mapgrid%nx )
                        end do
                    end do

                    nf90_stat = nf90_put_var(ncid, varid, vloc(1:cx,1:cy), (/sx,sy,ti/), (/cx, cy,1/))

                    !
                    ! If the write failed, check the values against the limits configured in nctablemd.ftn90
                    !
                    if ( nf90_stat /= NF90_NOERR ) then
                        if ( values_in_range_double(ncid, varid, values, &
                                fill_value, scale_factor, add_offset, varname) ) then
                            call agnc_set_error( nf90_strerror(nf90_stat))
                        else
                            if (.not. present(skip_range_error)) call agnc_set_error( nf90_strerror(nf90_stat))
                        end if
                    end if
                end do
            end do
            deallocate(vloc)

       end subroutine agnc_add_mapdata_double

       subroutine agnc_add_spcdata_density(ncid, ti, e2d, spc_as_map)
            integer,                          intent(in   ) :: ncid, ti
            real(kind=4),   dimension(:,:),   intent(inout) :: e2d
            logical,                          intent(in   ) :: spc_as_map

            real(kind=4),   dimension(:,:,:), allocatable   :: density
            real(kind=4),   dimension(:,:,:,:), allocatable :: edloc
            real(kind=4),   dimension(:),     allocatable   :: scale_density
            integer                                         :: varid_scale_density, &
                                                               varid_density
            real(kind=4)                                    :: fv_density
            integer                                         :: ip, ith, ik, isp, msc, &
                                                               chunksize, sx, cx, sy, cy, &
                                                               lxi, lyi
            type(spcgrid_type)                              :: spcgrid
            type(pntgrid_type)                              :: pntgrid
            type(mapgrid_type)                              :: mapgrid
            type(recordaxe_type)                            :: recordaxe
            logical                                         :: do_scale
            chunksize = 64
            do_scale = .false.

            call agnc_get_griddef(ncid, spcgrid)
            call agnc_get_recordaxe(ncid, recordaxe)
            call agnc_collect_spcmeta_2d(ncid, varid_density, varid_scale_density, fv_density)

            if ( spc_as_map ) then
                call agnc_get_griddef(ncid, mapgrid)
                msc = mapgrid%nx * mapgrid%ny
            else
                call agnc_get_griddef(ncid, pntgrid)
                msc = pntgrid%npoints
            end if

            ! should I scale and write in a big chuncked loop over msc?
            ! Allocating the density array creates a new copy for all 2D spectral points in the heap.
            allocate(density(spcgrid%ndir, spcgrid%nfreq, msc))

            !
            ! - scale and/or replace energy dummies with _FillValue
            ! - if scaling requested, write scale factors to netcdf
            !
            if ( varid_scale_density > 0 ) then
                do_scale = .true.
                allocate(scale_density(msc))

                ! Note the assumption that AGNC_FILL_SHORT is negative
                scale_density = maxval(e2d, 1) / (2**15-1)

                do ip=1, msc
                    if (scale_density(ip) < epsilon(1.) ) scale_density(ip) = epsilon(1.)
                    do ik=1, spcgrid%nfreq
                        do ith=1, spcgrid%ndir
                            isp = ith + (ik - 1) * spcgrid%ndir
                            if ( abs(e2d(isp,ip) - AGNC_DUMMY) > epsilon(1.) ) then
                                density(ith,ik,ip) = nint(e2d(isp,ip) / scale_density(ip))
                            else
                                density(ith,ik,ip) = AGNC_FILL_SHORT
                            end if
                        end do
                    end do
                end do

                if ( spc_as_map ) then
                    call nccheck ( nf90_put_var(ncid, varid_scale_density, &
                                    reshape(scale_density, (/mapgrid%nx, mapgrid%ny, 1/)), (/1, 1, ti/)) )
                else
                    call nccheck ( nf90_put_var(ncid, varid_scale_density, scale_density, (/1, ti/)) )
                end if
            else
                ! replace unscaled energy dummies with _FillValue
                do ip=1, msc
                    do ik=1, spcgrid%nfreq
                        do ith=1, spcgrid%ndir
                            isp = ith + (ik - 1) * spcgrid%ndir
                            if ( abs(e2d(isp,ip) - AGNC_DUMMY) > epsilon(1.) ) then
                                density(ith,ik,ip) = e2d(isp,ip)
                            else
                                density(ith,ik,ip) = NF90_FILL_FLOAT
                            end if
                        end do
                    end do
                end do
            end if

            !
            ! chunked write to netcdf
            !
            if ( spc_as_map ) then
                allocate(edloc(spcgrid%ndir, spcgrid%nfreq, chunksize, chunksize))

                do sx=1,mapgrid%nx,chunksize
                    cx = minval( (/chunksize, mapgrid%nx-sx+1/) )

                    do sy=1,mapgrid%ny,chunksize
                        cy = minval( (/chunksize, mapgrid%ny-sy+1/) )
                        edloc = 0.
                        do lxi=1,cx
                            do lyi=1,cy
                                edloc(:,:,lxi,lyi) = density(:,:, (sx + lxi - 1) + (sy + lyi - 2) * mapgrid%nx )
                            end do
                        end do
                        call nccheck ( nf90_put_var(ncid, varid_density, edloc(:,:,1:cx,1:cy), &
                                                                         (/1,1,sx,sy,ti/),     &
                                                                         (/spcgrid%ndir, spcgrid%nfreq, cx, cy,1/)) )
                    end do
                end do
                deallocate(edloc)
            else
                do ip=1,msc, chunksize**2
                    ! abuse sx and cx for size and end index in pointlist
                    sx = minval( (/chunksize**2, msc-ip+1/) )
                    cx = ip + sx - 1

                    call nccheck ( nf90_put_var(ncid, varid_density, density(:,:,ip:cx), &
                                                                      (/1,1,ip,ti/),      &
                                                                      (/spcgrid%ndir, spcgrid%nfreq, sx, 1/) ) )
                end do
            end if

            deallocate(density)

            if (do_scale) deallocate(scale_density)

        end subroutine agnc_add_spcdata_density

        subroutine agnc_add_spcdata_3d(ncid, ti_start, energy, theta, spr, spc_as_map)
            integer,                           intent(in   )    :: ncid, ti_start
            real(kind=4),    dimension(:,:,:), intent(inout)    :: energy, theta, spr
            logical,                           intent(in   )    :: spc_as_map

            ! local variables
            integer                                             :: s(3), ti

            s  = shape(energy)
            do ti = ti_start, ti_start + s(3) - 1
                call agnc_add_spcdata_2d(ncid, ti, reshape(energy(:,:, ti), (/s(1), s(2)/)), &
                                                   reshape( theta(:,:, ti), (/s(1), s(2)/)), &
                                                   reshape(   spr(:,:, ti), (/s(1), s(2)/)), spc_as_map)
            end do

        end subroutine agnc_add_spcdata_3d

        subroutine agnc_add_spcdata_2d(ncid, ti, energy, theta, spr, spc_as_map)
            integer,                           intent(in   )    :: ncid, ti
            real(kind=4),    dimension(:,:)                     :: energy, theta, spr
            logical,                           intent(in   )    :: spc_as_map

            ! local variables
            real(kind=4),    dimension(:),       allocatable    :: scale_energy
            real(kind=4),    dimension(:,:,:),   allocatable    :: eloc, tloc, sloc
            integer                                             :: varid_scale_energy, &
                                                                   varid_energy, varid_theta, &
                                                                   varid_spr, msc, nt, &
                                                                   ip, ik, ir, &
                                                                   chunksize, sx, cx, sy, cy, &
                                                                   lxi, lyi
            real(kind=4)                                        :: sf_theta, sf_spr, &
                                                                   ao_theta, ao_spr, &
                                                                   fv_theta, fv_spr
            type(spcgrid_type)                                  :: spcgrid
            type(pntgrid_type)                                  :: pntgrid
            type(mapgrid_type)                                  :: mapgrid
            type(recordaxe_type)                                :: recordaxe
            logical                                             :: do_scale
            chunksize = 64
            do_scale = .false.

            call agnc_get_griddef(ncid, spcgrid)
            call agnc_get_recordaxe(ncid, recordaxe)
            call agnc_collect_spcmeta(ncid, &
                                      varid_scale_energy, varid_energy,         &
                                      varid_theta, ao_theta, sf_theta, fv_theta,&
                                      varid_spr,   ao_spr,   sf_spr,   fv_spr)

            if ( spc_as_map ) then
                call agnc_get_griddef(ncid, mapgrid)
                msc = mapgrid%nx * mapgrid%ny
            else
                call agnc_get_griddef(ncid, pntgrid)
                msc = pntgrid%npoints
            end if
            if ( msc /= size(energy,2) ) then
                write(6,*) '#points of this run and the netcdf file is different'
                STOP 2
            end if

            if ( varid_scale_energy > 0 ) then
                do_scale = .true.

                allocate(scale_energy(msc))

                ! Note the assumption that AGNC_DUMMY is negative (or zero)
                scale_energy = maxval(energy, 1) / (2**15-1)

                do ip=1, msc
                    ! scale or replace energy dummies with _FillValue
                    if (scale_energy(ip) < epsilon(1.) ) scale_energy(ip) = epsilon(1.)
                    do ik=1, spcgrid%nfreq
                        if ( abs(energy(ik,ip) - AGNC_DUMMY) > epsilon(1.)) then
                            energy(ik,ip) = nint(energy(ik,ip) / scale_energy(ip))
                        else
                            energy(ik,ip) = AGNC_FILL_SHORT
                        end if
                    end do
                end do

                if ( spc_as_map ) then
                    ! make a clone of scale_energy. Reshape is acceptable.....
                    call nccheck ( nf90_put_var(ncid, varid_scale_energy, &
                                    reshape(scale_energy, (/mapgrid%nx, mapgrid%ny, 1/)), (/1, 1, ti/)) )
                else
                    call nccheck ( nf90_put_var(ncid, varid_scale_energy, scale_energy, (/1, ti/)) )
                end if
            else
                ! replace unscaled energy dummies with _FillValue
                do ip=1, msc
                    do ik=1, spcgrid%nfreq
                        if ( abs(energy(ik,ip) - AGNC_DUMMY) < 2*epsilon(1.)) &
                          energy(ik,ip) = NF90_FILL_FLOAT
                    end do
                end do
            end if

            ! replace direction and spread dummies with _FillValue
            do ip=1, msc
                do ik=1, spcgrid%nfreq
                    if ( abs(theta(ik,ip) - AGNC_DUMMY) > epsilon(1.) ) then
                        theta(ik,ip) = nint((theta(ik,ip) - ao_theta) / sf_theta)
                    else
                        theta(ik,ip) = AGNC_FILL_BYTE
                    end if
                    if ( abs(spr(ik,ip) - AGNC_DUMMY) > epsilon(1.) ) then
                        spr(ik,ip) = nint((spr(ik,ip) - ao_spr) / sf_spr)
                    else
                        spr(ik,ip) = AGNC_FILL_BYTE
                    end if
                end do
            end do

            !
            ! chunked write to netcdf
            !
            if ( spc_as_map ) then
                allocate(eloc(spcgrid%nfreq, chunksize, chunksize))
                allocate(tloc(spcgrid%nfreq, chunksize, chunksize))
                allocate(sloc(spcgrid%nfreq, chunksize, chunksize))

                do sx=1,mapgrid%nx,chunksize
                    cx = minval( (/chunksize, mapgrid%nx-sx+1/) )

                    do sy=1,mapgrid%ny,chunksize
                        cy = minval( (/chunksize, mapgrid%ny-sy+1/) )
                        eloc = 0.
                        tloc = 0.
                        sloc = 0.
                        do lxi=1,cx
                            do lyi=1,cy
                                eloc(:,lxi,lyi) = energy(:, (sx + lxi - 1) + (sy + lyi - 2) * mapgrid%nx)
                                tloc(:,lxi,lyi) = theta (:, (sx + lxi - 1) + (sy + lyi - 2) * mapgrid%nx)
                                sloc(:,lxi,lyi) = spr   (:, (sx + lxi - 1) + (sy + lyi - 2) * mapgrid%nx)
                            end do
                        end do
                        call nccheck ( nf90_put_var(ncid, varid_energy, eloc(:,1:cx,1:cy), &
                                                                         (/1,sx,sy,ti/),   &
                                                                         (/spcgrid%nfreq, cx, cy,1/)) )
                        call nccheck ( nf90_put_var(ncid, varid_theta , tloc(:,1:cx,1:cy), &
                                                                         (/1,sx,sy,ti/),   &
                                                                         (/spcgrid%nfreq, cx, cy,1/)) )
                        call nccheck ( nf90_put_var(ncid, varid_spr   , sloc(:,1:cx,1:cy), &
                                                                         (/1,sx,sy,ti/),   &
                                                                         (/spcgrid%nfreq, cx, cy,1/)) )
                    end do
                end do
                deallocate(eloc)
                deallocate(tloc)
                deallocate(sloc)
            else
                do ip=1,msc, chunksize**2
                    ! abuse sx and cx for size and end index in pointlist
                    sx = minval( (/chunksize**2, msc-ip+1/) )
                    cx = ip + sx - 1

                    call nccheck ( nf90_put_var(ncid, varid_energy, energy(:,ip:cx), &
                                                                     (/1,ip,ti/),    &
                                                                     (/spcgrid%nfreq, sx,  1/) ))
                    call nccheck ( nf90_put_var(ncid, varid_theta,  theta(:,ip:cx),  &
                                                                     (/1,ip,ti/),    &
                                                                     (/spcgrid%nfreq, sx,  1/) ))
                    call nccheck ( nf90_put_var(ncid, varid_spr,    spr(:,ip:cx),    &
                                                                     (/1,ip,ti/),    &
                                                                     (/spcgrid%nfreq, sx,  1/) ))
                end do
            end if
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
            integer                                             :: varid, nf90_stat, xtype

            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_float(ncid, varid, add_offset, scale_factor, fill_value, xtype)

            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            if ( xtype == NF90_SHORT .or. xtype == NF90_BYTE .or. xtype == NF90_INT ) then
                where ( abs(values - fill_value) > epsilon(1.) ) &
                    values = nint( (values - add_offset) / scale_factor)
            end if

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
            integer                                             :: varid, nf90_stat, xtype

            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_float(ncid, varid, add_offset, scale_factor, fill_value, xtype)

            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            if ( xtype == NF90_SHORT .or. xtype == NF90_BYTE .or. xtype == NF90_INT ) then
                where ( abs(values - fill_value) > epsilon(1.) ) &
                    values = nint( (values - add_offset) / scale_factor)
                            nf90_stat = nf90_put_var(ncid, varid, values, (/1, ti/))
            end if

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
            integer                                             :: varid, nf90_stat, xtype

            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_float(ncid, varid, add_offset, scale_factor, fill_value, xtype)

            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            if ( xtype == NF90_SHORT .or. xtype == NF90_BYTE .or. xtype == NF90_INT ) then
                where ( abs(values - fill_value) > epsilon(1.) ) &
                    values = nint( (values - add_offset) / scale_factor)
            end if

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
            integer                                             :: varid, nf90_stat, xtype
            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_double(ncid, varid, add_offset, scale_factor, fill_value, xtype)

            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            if ( xtype == NF90_SHORT .or. xtype == NF90_BYTE .or. xtype == NF90_INT ) then
                where ( abs(values - fill_value) > epsilon(1.) ) &
                    values = nint( (values - add_offset) / scale_factor)
            end if

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
            integer                                             :: varid, nf90_stat, xtype
            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_double(ncid, varid, add_offset, scale_factor, fill_value, xtype)

            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            if ( xtype == NF90_SHORT .or. xtype == NF90_BYTE .or. xtype == NF90_INT ) then
                where ( abs(values - fill_value) > epsilon(1.) ) &
                    values = nint( (values - add_offset) / scale_factor)
            end if

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
            integer                                             :: varid, nf90_stat, xtype
            call agnc_get_varid_by_name(ncid, varname, varid)
            call get_scalies_double(ncid, varid, add_offset, scale_factor, fill_value, xtype)

            if ( present(dummyvalue) ) then
                where (abs(values - dummyvalue) < epsilon(1.)) values = fill_value
            end if

            ! packed_data_value = nint((unpacked_data_value - add_offset) / scale_factor)
            if ( xtype == NF90_SHORT .or. xtype == NF90_BYTE .or. xtype == NF90_INT ) then
                where ( abs(values - fill_value) > epsilon(1.) ) &
                    values = nint( (values - add_offset) / scale_factor)
            end if

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
!       compute_1d_spectra(energy, e1, th1, sp1, undef) done

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
                if ( trim(name) == trim(searchfor) ) return

                if ( nf90_get_att(ncid, varid, "standard_name", name) == NF90_NOERR ) then
                    if ( trim(name) == trim(searchfor) ) return
                end if

                if (nf90_get_att(ncid, varid, "short_name", name) == NF90_NOERR ) then
                    if ( trim(name) == trim(searchfor) ) return
                end if
                if (nf90_get_att(ncid, varid, "long_name", name) == NF90_NOERR ) then
                    if ( trim(name) == trim(searchfor) ) return
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
            integer, intent(in)                :: x(:),dx,nx
            logical, intent(in), optional      :: fill_check
            logical                            :: isregular
            integer, dimension(:), allocatable :: dx2

            isregular = .false.
            if ( nx > 0 ) then
                if ( present(fill_check) ) then
                    if  ( fill_check .and. all( (x - NF90_FILL_INT) == 0) ) return
                end if
                if ( nx > 1 ) then
                    allocate(dx2(nx-1))
                    dx2 = x(2:nx) - x(1:nx-1)
                    if ( any( (dx2 - dx) /= 0 )) then
                        deallocate(dx2)
                        return
                    end if
                    deallocate(dx2)
                end if
                isregular = .true.
            end if

        end function axe_is_regular_integer

        function axe_is_regular_integer64(x, nx, dx, fill_check) result (isregular)
            integer*8, intent(in)                :: x(:)
            integer,   intent(in)                :: dx,nx
            logical,   intent(in), optional      :: fill_check
            logical                              :: isregular
            integer,   dimension(:), allocatable :: dx2

            isregular = .false.
            if ( nx > 0 ) then
                if ( present(fill_check) ) then
                    if  ( fill_check .and. all( (x - NF90_FILL_INT) == 0) ) return
                end if
                if ( nx > 1 ) then
                    allocate(dx2(nx-1))
                    dx2 = x(2:nx) - x(1:nx-1)
                    if ( any( (dx2 - dx) /= 0 )) then
                        deallocate(dx2)
                        return
                    end if
                    deallocate(dx2)
                end if
                isregular = .true.
            end if

        end function axe_is_regular_integer64

        function axe_is_regular_float(x, nx, dx, fill_check) result (isregular)
            real(kind=4),intent(in)                :: x(:),dx
            integer,     intent(in)                :: nx
            logical,     intent(in), optional      :: fill_check
            logical                                :: isregular
            real(kind=4),dimension(:), allocatable :: dx2

            isregular = .false.
            if ( nx > 0 ) then
                if ( present(fill_check) ) then
                    if  ( fill_check .and. all( (x - NF90_FILL_INT) == 0) ) return
                end if
                if ( nx > 1 ) then
                    allocate(dx2(nx-1))
                    dx2 = x(2:nx) - x(1:nx-1)
                    if ( any( abs(dx2 - dx) > 1e-4 )) then
                        deallocate(dx2)
                        return
                    end if
                    deallocate(dx2)
                end if
                isregular = .true.
            end if

        end function axe_is_regular_float

        function axe_is_regular_float_2d(x, nx, dx, fill_check) result (isregular)
            real(kind=4),intent(in)             :: x(:,:),dx
            integer,     intent(in)             :: nx
            logical,     intent(in), optional   :: fill_check
            logical                             :: isregular
            real(kind=4)                        :: dx2

            ! The reason for 2D meshes is that we can use the SWAN rotated grid feature
            !
            ! I have never used this before so I would really know how to deal with it
            ! as I'm not sure how x0, x1, and alpha are implemented.
            isregular = .true.

        end function axe_is_regular_float_2d


        function axe_is_regular_double(x, nx, dx, fill_check) result(isregular)
            real(kind=8),intent(in)                :: x(:),dx
            integer,     intent(in)                :: nx
            logical,     intent(in), optional      :: fill_check
            logical                                :: isregular
            real(kind=8),dimension(:), allocatable :: dx2

            isregular = .false.
            if ( nx > 0 ) then
                if ( present(fill_check) ) then
                    if  ( fill_check .and. all( (x - NF90_FILL_INT) == 0) ) return
                end if
                if ( nx > 1 ) then
                    allocate(dx2(nx-1))
                    dx2 = x(2:nx) - x(1:nx-1)
                    if ( any( abs(dx2 - dx) > 1e-4 )) then
                        deallocate(dx2)
                        return
                    end if
                    deallocate(dx2)
                end if
                isregular = .true.
            end if

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

        function values_in_range_float(ncid, varid, values, fill_value, scale_factor, add_offset, varname) result (inrange)
            integer,                    intent(   in)         :: ncid, varid
            real(kind=4), dimension(:), intent(   in)         :: values
            real(kind=4),               intent(   in)         :: fill_value, scale_factor, add_offset
            character(len=*),           intent(   in)         :: varname

            ! local
            integer                                           :: xtype, ind1d, fillv_int, value
            logical                                           :: inrange
            real(kind=4)                                      :: bmax, bmin
            inrange = .true.

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
                bmin = -2147483647.
                bmax = 2147483647.
            case default
                return
            end select

            if ( (fill_value - bmin - 1) > epsilon(1.) ) then
                ! This is not the usual NF90_FILL_XXXX value, so compute the integer value of the fill value
                fillv_int = nint( (fill_value - add_offset) / scale_factor )
            else
                fillv_int = nint(bmin - 1)
            end if

            do ind1d=1,size(values)
                value = nint(values(ind1d))
                if ( value /= fillv_int .and. value .lt. bmin .or. value .gt. bmax ) then
                    write(6, '("The value in 1d index ", I12, " lies outside the boundaries provided in nctablemd.f90")') ind1d
                    write(6, *) "fill value   int / float: ", fill_value, fill_value * scale_factor + add_offset
                    write(6, *) "Theoretical    min / max: ", &
                                               bmin * scale_factor + add_offset, &
                                               bmax * scale_factor + add_offset
                    write(6, *) "Actual value int / float: ", value, " / ", &
                                                   values(ind1d) * scale_factor + add_offset
                    inrange = .false.
                    return
                end if
            end do

        end function values_in_range_float

        function values_in_range_double(ncid, varid, values, fill_value, scale_factor, add_offset, varname) result (inrange)
            integer,                    intent(   in)         :: ncid, varid
            real(kind=8), dimension(:), intent(   in)         :: values
            real(kind=8),               intent(   in)         :: fill_value, scale_factor, add_offset
            character(len=*),           intent(   in)         :: varname

            ! local
            integer                                           :: xtype, ind1d, fillv_int, value
            logical                                           :: inrange
            real(kind=8)                                      :: bmax, bmin
            inrange = .true.

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
                bmin = -2147483647.
                bmax = 2147483647.
            case default
                return
            end select

            if ( (fill_value - bmin - 1) > epsilon(1.) ) then
                ! This is not the usual NF90_FILL_XXXX value, so compute the integer value of the fill value
                fillv_int = nint( (fill_value - add_offset) / scale_factor )
            else
                fillv_int = nint(bmin - 1)
            end if

            do ind1d=1,size(values)
                value = nint(values(ind1d))
                if ( value /= fillv_int .and. value .lt. bmin .or. value .gt. bmax ) then
                    write(6, '("The value in 1d index ", I12, " lies outside the boundaries provided in nctablemd.f90")') ind1d
                    write(6, *) "fill value               : ", fill_value
                    write(6, *) "Theoretical     min / max: ", &
                                               bmin * scale_factor + add_offset, &
                                               bmax * scale_factor + add_offset
                    write(6, *) "Actual values int / float: ", value, " / ", &
                                                   values(ind1d) * scale_factor + add_offset
                    inrange = .false.
                    return
                end if
            end do

        end function values_in_range_double

        function timeindex32(tarr, t) result (ti)
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
        end function timeindex32

        function timeindex64(tarr, t) result (ti)
            integer*8, dimension(:), intent( in)        :: tarr
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
        end function timeindex64

        subroutine agnc_collect_spcmeta(ncid, varid_scale_energy, varid_energy, &
                                        varid_theta, ao_theta, sf_theta, fv_theta, &
                                        varid_spr,   ao_spr,   sf_spr,   fv_spr)
            integer,                           intent( in)      :: ncid
            integer,                           intent(out)      :: varid_scale_energy, &
                                                                   varid_energy, varid_theta, &
                                                                   varid_spr
            real(kind=4),                      intent(out)      :: sf_theta, sf_spr, &
                                                                   ao_theta, ao_spr, &
                                                                   fv_theta, fv_spr
            call agnc_get_varid_by_name(ncid, "scale_energy_1d", varid_scale_energy)
            call agnc_get_varid_by_name(ncid, "energy_1d", varid_energy)
            call agnc_get_varid_by_name(ncid, "theta_1d", varid_theta)
            call agnc_get_varid_by_name(ncid, "spread_1d", varid_spr)
            call get_scalies_float(ncid, varid_theta, ao_theta, sf_theta, fv_theta)
            call get_scalies_float(ncid, varid_spr, ao_spr, sf_spr, fv_spr)
        end subroutine agnc_collect_spcmeta

        subroutine agnc_collect_spcmeta_2d(ncid, varid_density, varid_scale_density, fv_density)
            integer,                           intent( in)      :: ncid
            integer,                           intent(out)      :: varid_scale_density, &
                                                                   varid_density
            real(kind=4)                                        :: ao_density, sf_density
            real(kind=4),                      intent(out)      :: fv_density

            call agnc_get_varid_by_name(ncid, "density", varid_density)
            call agnc_get_varid_by_name(ncid, "scale_density", varid_scale_density)
            call get_scalies_float(ncid, varid_density, ao_density, sf_density, fv_density)
        end subroutine agnc_collect_spcmeta_2d

        function datevec_from_epoch( t_in ) result (datevec)
            integer*8, intent( in)  :: t_in
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

        function seconds_since_epoch_twoint ( date, time ) result (seconds_since_epoch)
            integer, intent( in)    :: date, time
            integer*8               :: seconds_since_epoch
            seconds_since_epoch = seconds_since_epoch_datevec( datevec_from_twoint(date, time) )
        end function seconds_since_epoch_twoint

        function seconds_since_epoch_datevec ( datevec ) result (seconds_since_epoch)
        !     convert date time in (int yyyymmdd, hhmmss) format to seconds since 1-jan-1970
            integer                 :: datevec(6)
            integer*8               :: seconds_since_epoch
            integer                 :: year, month, y1, y2, f
            ! number of days per month
            integer, save           :: ndpm(12)
            data    ndpm / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
            seconds_since_epoch = 0
            y1 = 1970
            y2 = 1970
            f  = 1

            if (datevec(2) > 12) then
                datevec(1) = datevec(1) + datevec(2)/12
                datevec(2) = mod(datevec(2),12)
            end if

            if (datevec(1) < 1970) then
                y1 = datevec(1)+1
                f  = -1
            else
                y2 = datevec(1)-1
            end if

            do year = y1, y2
                if ( (mod(year,400) == 0) .or. ( (mod(year,4) == 0) .and. (mod(year,100) /= 0 )) ) then
                    seconds_since_epoch = seconds_since_epoch + 366*86400 * f
                else
                    seconds_since_epoch = seconds_since_epoch + 365*86400 * f
                end if
            end do

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

        subroutine agnc_scan_time(ustr, factor, dvec, emsg)
            character(len=*),      intent( in) :: ustr
            integer,               intent(out) :: factor, dvec(6)
            character(len=*),      intent(out) :: emsg
            integer                            :: j, slen, nspace
            character                          :: fstr*20, dlim1, dlim2

            factor = -1
            emsg   = ' '
            slen   = len(trim(ustr))
            j    = index(ustr, 'since')
            if ( j < 1 ) then
                emsg = 'Keyword "since" not found in unit string '//trim(ustr)
                return
            end if

            fstr = trim(ustr(1:j-1))
            if      ( fstr .eq. 'second'  .or. &
                    fstr .eq. 'seconds' .or. &
                    fstr .eq. 'sec'     .or.     &
                    fstr .eq. 's' ) then
            factor = 1
            elseif ( fstr .eq. 'minute'  .or. &
                    fstr .eq. 'minutes' .or. &
                    fstr .eq. 'min'     .or.     &
                    fstr .eq. 'm' ) then
            factor = 60
            elseif ( fstr .eq. 'day'  .or. &
                    fstr .eq. 'days' .or. &
                    fstr .eq. 'd' ) then
            factor = 86400
            elseif ( fstr .eq. 'year'  .or. &
                    fstr .eq. 'years' ) then
            ! int32(365.242198781 * 86400)
            factor = 31556926
            elseif ( fstr .eq. 'week'  .or. &
                    fstr .eq. 'weeks' ) then
            ! int32(365.242198781 * 86400 / 52)
            factor = 606864
            elseif ( fstr .eq. 'month'  .or. &
                    fstr .eq. 'months' ) then
            ! int32(365.242198781 * 86400 / 12)
            factor = 2629744
            else
                write(emsg,'("Time unit ", A, " not supported")') trim(fstr)
                return
            end if
            ! step to the character after the 'since ', start of yyyy-mm-dd, mm and dd
            ! do not have to be zero padded...
            j   = j + 6
            read(ustr(j:slen),'(I4,A,I2,A,I2)') dvec(1), dlim1, dvec(2), dlim2, dvec(3)

            ! find the starting charter of the time
            nspace = index(ustr(j:slen), ' ')

            if (nspace > 0 .and. nspace < slen) then
                j  = j+nspace
                read(ustr(j:slen),'(I2,A,I2,A,I2)') dvec(4), dlim1, dvec(5), dlim2, dvec(6)
            else
                ! no time provided
                dvec(4:6) = 0
            end if

        end subroutine agnc_scan_time

    end module agioncmd
