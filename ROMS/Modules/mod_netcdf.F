#include "cppdefs.h"
      MODULE mod_netcdf
!
!svn $Id: mod_netcdf.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This MODULE contains all NetCDF variables definitions. It also      !
!  contains several variables and routines to facilitate  generic      !
!  manipulations of NetCDF data.  Usually, several NetCDF library      !
!  calls are required to inquire and read a dimension or variable.     !
!  These routines provide a single interface for such operations.      !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE netcdf
!
      implicit none
!
      INTERFACE netcdf_get_fvar
        MODULE PROCEDURE netcdf_get_fvar_0d
        MODULE PROCEDURE netcdf_get_fvar_1d
        MODULE PROCEDURE netcdf_get_fvar_2d
        MODULE PROCEDURE netcdf_get_fvar_3d
        MODULE PROCEDURE netcdf_get_fvar_4d
      END INTERFACE netcdf_get_fvar
!
      INTERFACE netcdf_get_lvar
        MODULE PROCEDURE netcdf_get_lvar_0d
        MODULE PROCEDURE netcdf_get_lvar_1d
      END INTERFACE netcdf_get_lvar
!
      INTERFACE netcdf_get_ivar
        MODULE PROCEDURE netcdf_get_ivar_0d
        MODULE PROCEDURE netcdf_get_ivar_1d
        MODULE PROCEDURE netcdf_get_ivar_2d
      END INTERFACE netcdf_get_ivar
!
      INTERFACE netcdf_get_svar
        MODULE PROCEDURE netcdf_get_svar_0d
        MODULE PROCEDURE netcdf_get_svar_1d
        MODULE PROCEDURE netcdf_get_svar_2d
        MODULE PROCEDURE netcdf_get_svar_3d
      END INTERFACE netcdf_get_svar
!
      INTERFACE netcdf_get_time
        MODULE PROCEDURE netcdf_get_time_0d
        MODULE PROCEDURE netcdf_get_time_1d
      END INTERFACE netcdf_get_time
!
      INTERFACE netcdf_put_fvar
#ifdef SINGLE_PRECISION
        MODULE PROCEDURE netcdf_put_fvar_0dp
        MODULE PROCEDURE netcdf_put_fvar_1dp
        MODULE PROCEDURE netcdf_put_fvar_2dp
#endif
        MODULE PROCEDURE netcdf_put_fvar_0d
        MODULE PROCEDURE netcdf_put_fvar_1d
        MODULE PROCEDURE netcdf_put_fvar_2d
        MODULE PROCEDURE netcdf_put_fvar_3d
        MODULE PROCEDURE netcdf_put_fvar_4d
      END INTERFACE netcdf_put_fvar
!
      INTERFACE netcdf_put_ivar
        MODULE PROCEDURE netcdf_put_ivar_0d
        MODULE PROCEDURE netcdf_put_ivar_1d
        MODULE PROCEDURE netcdf_put_ivar_2d
      END INTERFACE netcdf_put_ivar
!
      INTERFACE netcdf_put_lvar
        MODULE PROCEDURE netcdf_put_lvar_0d
        MODULE PROCEDURE netcdf_put_lvar_1d
        MODULE PROCEDURE netcdf_put_lvar_2d
      END INTERFACE netcdf_put_lvar
!
      INTERFACE netcdf_put_svar
        MODULE PROCEDURE netcdf_put_svar_0d
        MODULE PROCEDURE netcdf_put_svar_1d
        MODULE PROCEDURE netcdf_put_svar_2d
        MODULE PROCEDURE netcdf_put_svar_3d
      END INTERFACE netcdf_put_svar
!
!  Switch to debug creating, opening, and closing of NetCDF file IDs
!  to monitor 'Too many open files' error. In the Unix environment
!  there is a limit to the number of open files. Use 'ulimit -a' or
!  'ulimit -S -n' to check.
!
#ifdef CHECK_OPEN_FILES
      logical :: Ldebug_ncid = .TRUE.
#else
      logical :: Ldebug_ncid = .FALSE.
#endif
      integer :: DBout = 1000           ! debug file 'fort.1000'
!
!  Local dimension parameters.
!
      integer, parameter :: Matts = 50  ! maximum number of attributes
      integer, parameter :: Mdims = 50  ! maximum number of dimensions
      integer, parameter :: Mvars = 900 ! maximum number of variables
      integer, parameter :: NvarD = 5   ! number of variable dimensions
      integer, parameter :: NvarA = 50  ! number of variable attributes
!
!  Generic information about current NetCDF for all dimensions and
!  all variables.
!
      integer :: n_dim                  ! number of dimensions
      integer :: n_var                  ! number of variables
      integer :: n_gatt                 ! number of global attributes
      integer :: ncformat               ! file format number (version)
      integer :: rec_id                 ! unlimited dimension ID
      integer :: rec_size               ! unlimited dimension value
      integer :: att_kind(Matts)        ! attribute data type
      integer :: dim_id(Mdims)          ! dimensions ID
      integer :: dim_size(Mdims)        ! dimensions value
      integer :: var_id(Mvars)          ! variables ID
      integer :: var_natt(Mvars)        ! variables number of attributes
      integer :: var_flag(Mvars)        ! Variables water points flag
      integer :: var_type(Mvars)        ! variables external data type
      integer :: var_ndim(Mvars)        ! variables number of dimensions
      integer :: var_dim(NvarD,Mvars)   ! variables dimensions ID
!
      character (len=100) :: att_name(Matts)     ! attribute names
      character (len=100) :: dim_name(Mdims)     ! dimension names
      character (len=100) :: var_name(Mvars)     ! variable names
!
!  Generic information about requested current variable.
!
      integer :: n_vdim                 ! number of variable dimensions
      integer :: n_vatt                 ! number of variable attributes
      integer :: var_kind               ! external data type
      integer :: var_Dids(NvarD)        ! dimensions ID
      integer :: var_Dsize(NvarD)       ! dimensions values
      integer :: var_Aint(NvarA)        ! attribute integer values
      real(r8) :: var_Afloat(NvarA)     ! attribute float values
!
      character (len=100)  :: var_Aname(NvarA)   ! Attribute names
      character (len=100)  :: var_Dname(NvarD)   ! dimension names
      character (len=1024) :: var_Achar(NvarA)   ! Attribute char values
!
!  External data representation for floating-point variables.
!
#ifdef OUT_DOUBLE
      integer, parameter :: NF_FOUT = nf90_double
#else
      integer, parameter :: NF_FOUT = nf90_real
#endif
#ifdef RST_SINGLE
      integer, parameter :: NF_FRST = nf90_real
#else
      integer, parameter :: NF_FRST = nf90_double
#endif
#ifdef DOUBLE_PRECISION
      integer, parameter :: NF_TYPE = nf90_double
#else
      integer, parameter :: NF_TYPE = nf90_real
#endif
!
!  External data representation for floating-point time and depth
!  variables. It is set to double precision for accuaracy in both
!  single and douple precision numerical kernel.
!
      integer, parameter :: NF_TOUT = nf90_double
!
!  Netcdf file basic creation mode flag. Its value is further modified
!  in "netcdf_create" using the IOR function to include additional
!  flags, like "nf90_clobber" to overwrite existing file.
!
#if defined HDF5 || defined PARALLEL_IO
      integer :: CMODE = nf90_netcdf4      ! NetCDF-4/HDF5 format file
#else
      integer :: CMODE = nf90_64bit_offset ! NetCDF 64-bit offset format
!                                            (large file support)
#endif
!
!  Set NetCDF-4/HDF5 deflate (file compression) parameters.
!
      integer :: shuffle = 1
      integer :: deflate = 1
      integer :: deflate_level = 1

#ifdef PARALLEL_IO
!
!  Set parallel I/O access for processing non-tiled and tiled data.
!  Two modes of parallel I/O access are possible: Independent and
!  Collective. Independent I/O access means that processing do not
!  depend on or be affected by other parallel processes (nodes).
!  This is the case in non-tiled I/O variables. Contrarily,
!  Collective I/O access implies that all parallel processess
!  participate during processing. This is the case for tiled
!  variables: each node in the group reads/writes their own tile
!  data when parallel I/O is activated.
!
!  New versions of the NetCDF4/HDF5 libraries require that all
!  variables with unlimited dimension collective I/O access including
!  scalar variables like "ocean_time".
!
!  It is more efficient to use collective access for both non-tiled
!  and tiled variables.
!
!  Note: the parallel access flags nf90_independent and nf90_collective
!        were missing in module netcdf.mod in early versions of the
!        NetCDF 4.x library. Currently,
!
!            nf_independent = 0,    nf90_independent = 0
!            nf_collective  = 1,    nf90_collective  = 1
!
      integer, parameter :: IO_nontiled_access = nf90_independent
      integer, parameter :: IO_tiled_access    = nf90_collective
#endif
!
      CONTAINS
!
      SUBROUTINE netcdf_get_dim (ng, model, ncname, ncid, DimName,      &
     &                           DimSize, DimID)
!
!=======================================================================
!                                                                      !
!  This routine inquires a NetCDF file dimensions names and values.    !
!  All the dimension information is stored in the module variables     !
!  declared above.  In addition, if a particular dimension name is     !
!  provided, this routine returns the requested information in the     !
!  optional arguments.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     DimName      Requested dimension name (string, OPTIONAL)         !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     DimSize      Size of requested dimension (integer, OPTIONAL)     !
!     DimID        ID of requested dimension (integer, OPTIONAL)       !
!                                                                      !
!  Other information stored in this module:                            !
!                                                                      !
!     n_dim        Number of dimensions                                !
!     n_var        Number of variables                                 !
!     n_gatt       Number of global attributes                         !
!     rec_id       Unlimited dimension ID                              !
!     rec_size     Size of unlimited dimension                         !
!     dim_name     Dimensions name (1:n_dim)                           !
!     dim_id       Dimensions ID (1:n_dim)                             !
!     dim_size     Dimensions value (1:n_dim)                          !
!                                                                      !
!     WARNING: This is information is rewritten during each CALL.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcasts
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in), optional :: DimName
!
      integer, intent(out), optional :: DimSize
      integer, intent(out), optional :: DimID
!
!  Local variable declarations.
!
      logical :: foundit
!
      integer :: my_ncid, i, j, status
      integer :: myID, myValue

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(5) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_dim"
!
!-----------------------------------------------------------------------
!  Inquire about the NetCDF dimensions (names and values).
!-----------------------------------------------------------------------
!
!  Initialize.
!
      n_dim=0
      n_var=0
      n_gatt=0
      ncformat=-1
      rec_id=-1
      rec_size=0
      dim_id=0
      dim_size=0
      DO i=1,Mdims
        DO j=1,LEN(dim_name(1))
          dim_name(i)(j:j)=' '
        END DO
      END DO
!
!  Open file for reading.
!
      IF (.not.PRESENT(ncid)) THEN
         CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
         IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Inquire file.
!
      IF (InpThread) THEN
#ifdef HDF5
        status=nf90_inquire(my_ncid, n_dim, n_var, n_gatt, rec_id,      &
     &                      ncformat)
#else
        status=nf90_inquire(my_ncid, n_dim, n_var, n_gatt, rec_id)
#endif
        IF ((status.eq.nf90_noerr).and.(n_dim.le.Mdims)) THEN
#if !defined PARALLEL_IO  && defined DISTRIBUTE
          ibuffer(1)=n_dim
          ibuffer(2)=n_var
          ibuffer(3)=n_gatt
          ibuffer(4)=rec_id
#endif
!
!  Inquire about dimensions: names, ID, and size.
!
          rec_size=-1
          DO i=1,n_dim
            dim_id(i)=i
            status=nf90_inquire_dimension(my_ncid, dim_id(i),           &
     &                                    dim_name(i), dim_size(i))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              WRITE (stdout,10) dim_id(i), TRIM(ncname),                &
     &                          TRIM(SourceFile), nf90_strerror(status)
              exit_flag=2
              ioerror=status
              EXIT
            END IF
            IF (dim_id(i).eq.rec_id) THEN
              rec_size=dim_size(i)
#if !defined PARALLEL_IO  && defined DISTRIBUTE
              ibuffer(5)=rec_size
#endif
            END IF
          END DO
        ELSE
          IF (n_dim.gt.Mdims) THEN
            WRITE (stdout,20) ' Mdims = ', Mdims, n_dim
            exit_flag=2
            ioerror=0
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,30) TRIM(ncname), TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast dimension to all processors in the group.
!
      CALL mp_bcasti (ng, model, exit_flag)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, ibuffer)
        n_dim=ibuffer(1)
        n_var=ibuffer(2)
        n_gatt=ibuffer(3)
        rec_id=ibuffer(4)
        rec_size=ibuffer(5)
        CALL mp_bcasti (ng, model, dim_id)
        CALL mp_bcasti (ng, model, dim_size)
        CALL mp_bcasts (ng, model, dim_name)
      END IF
#endif
!
!  Load requested information.
!
      IF (exit_flag.eq.NoError) THEN
        foundit=.FALSE.
        IF (PRESENT(DimName)) THEN
          DO i=1,n_dim
            IF (TRIM(dim_name(i)).eq.TRIM(DimName)) THEN
              foundit=.TRUE.
              myID=dim_id(i)
              myValue=dim_size(i)
            END IF
          END DO
          IF (foundit) THEN
            IF (PRESENT(DimSize)) THEN
              DimSize=myValue
            END IF
            IF (PRESENT(DimID)) THEN
              DimID=myID
            END IF
          ELSE
            WRITE (stdout,40) TRIM(DimName), TRIM(ncname)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_DIM - error while reading dimension ID:',  &
     &        2x,i3,/,18x,'in input file:',2x,a,/,18x,'call from:',     &
     &        2x,a,/,18x,a)
  20  FORMAT (/,' NETCDF_GET_DIM - too small dimension parameter,',a,   &
     &        2i5,/,18x,'change file  mod_netcdf.F  and recompile')
  30  FORMAT (/,' NETCDF_GET_DIM - unable to inquire about contents',   &
     &          ' of input NetCDF file:',2x,a,/,18x,'call from:',2x,a)
  40  FORMAT (/,' NETCDF_GET_DIM - requested dimension: ',a,/18x,       &
     &        'not found in input file:',2x,a,/,18x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_dim
!
      SUBROUTINE netcdf_check_dim (ng, model, ncname, ncid)
!
!=======================================================================
!                                                                      !
!  This routine inquires a NetCDF file dimensions names and values.    !
!  It checks the file dimensions against model dimension parameters    !
!  for consitency.  All the dimensions information is stored in the    !
!  module variables declared above.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!                                                                      !
!  On output the following information is stored in this module:       !
!                                                                      !
!     n_dim        Number of dimensions                                !
!     n_var        Number of variables                                 !
!     n_gatt       Number of global attributes                         !
!     rec_id       Unlimited dimension ID                              !
!     rec_size     Size of unlimited dimension                         !
!     dim_name     Dimensions name (1:n_dim)                           !
!     dim_id       Dimensions ID (1:n_dim)                             !
!     dim_size     Dimensions value (1:n_dim)                          !
!                                                                      !
!     WARNING: This is information is rewritten during each CALL.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
!
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: i, status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_check_dim"
!
!-----------------------------------------------------------------------
!  Inquire about the NetCDF dimensions (names and values).
!-----------------------------------------------------------------------
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_get_dim (ng, model, ncname)
      ELSE
        CALL netcdf_get_dim (ng, model, ncname,                         &
     &                       ncid = ncid)
      END IF
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Check dimensions for consistency.
!-----------------------------------------------------------------------
!
      DO i=1,n_dim
        SELECT CASE (TRIM(ADJUSTL(dim_name(i))))
          CASE ('xi_psi')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_psi) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_psi,        &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_psi')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_psi) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_psi,       &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('xi_rho')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_rho) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_rho,        &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_rho')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_rho) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_rho,       &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('xi_u')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_u) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_u,          &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_u')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_u) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_u,         &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('xi_v')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_v) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_v,          &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_v')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_v) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_v,         &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('IorJ')
            IF (dim_size(i).ne.IOBOUNDS(ng)%IorJ) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%IorJ,          &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
#ifdef SOLVE3D
          CASE ('s_rho')
            IF (dim_size(i).ne.N(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i), N(ng),         &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
#endif
#ifndef RBL4DVAR_FCT_SENSITIVITY
# ifdef ADJUST_BOUNDARY
          CASE ('obc_adjust')
            IF (dim_size(i).ne.Nbrec(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i), Nbrec(ng),     &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
          CASE ('frc_adjust')
            IF (dim_size(i).ne.Nfrec(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i), Nfrec(ng),     &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
# endif
#endif
        END SELECT
      END DO

 10   FORMAT (/,' NETCDF_CHECK_DIM - inconsistent size of dimension: ', &
     &        a,2x,2i5,/,20x,'in file: ',a)
!
      RETURN
      END SUBROUTINE netcdf_check_dim
!
      SUBROUTINE netcdf_check_var (ng, model, ncname, ncid)
!
!=======================================================================
!                                                                      !
!  This routine inquires the NetCDF file variables and check if the    !
!  values of few of them are consitent with the parameters provided    !
!  in input scripts. All the variables information is stored in the    !
!  module variables declared above.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!                                                                      !
!  On output the following information is stored in this module:       !
!                                                                      !
!     n_dim        Number of dimensions                                !
!     n_var        Number of variables                                 !
!     n_gatt       Number of global attributes                         !
!     rec_id       Unlimited dimension ID                              !
!     var_name     Variables name (1:n_var)                            !
!     var_id       Variables ID (1:n_var)                              !
!     var_natt     Variables number of attributes (1:n_var)            !
!     var_flag     Variables flag [1=full field, -1=water points only] !
!     var_type     Variables external data type (1:n_var)              !
!     var_ndim     Variables number of dimensions (1:n_var)            !
!     var_dim      Variables dimensions ID (:,1:n_var)                 !
!                                                                      !
!     WARNING: This is information is rewritten during each CALL.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
#ifdef FOUR_DVAR
      USE strings_mod, ONLY : find_string
#endif
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
!
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: IDmod, Npts, i, ib, ic, j, j1, j2, status
      integer :: IvarS
!
      real(r8), parameter :: RoundOff = 1.0e-7_r8

      real(r8) :: FvarS, FvarV(50), VarVal
!
      character (len=80) :: text

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_check_var"
!
!-----------------------------------------------------------------------
!  Inquire about the NetCDF variables.
!-----------------------------------------------------------------------
!
!  Limit model identifier. The profiling is limited to iNLM, iTLM, iRPM,
!  and iADM.
!
      IF ((model.lt.1).or.(model.gt.4)) THEN
        IDmod=iNLM
      ELSE
        IDmod=model
      END IF
!
! Inquire about all variables.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_inq_var (ng, IDmod, ncname)
      ELSE
        CALL netcdf_inq_var (ng, IDmod, ncname,                         &
     &                       ncid = ncid)
      END IF
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Check several important variables for consistency.
!-----------------------------------------------------------------------
!
      DO i=1,n_var
        SELECT CASE (TRIM(ADJUSTL(var_name(i))))

#ifdef SOLVE3D
          CASE ('Vtransform')
            IF (.not.PRESENT(ncid)) THEN
              CALL netcdf_get_ivar (ng, IDmod, ncname, var_name(i),     &
     &                              IvarS)
            ELSE
              CALL netcdf_get_ivar (ng, IDmod, ncname, var_name(i),     &
     &                              IvarS,                              &
     &                              ncid = ncid)
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            IF (IvarS.ne.Vtransform(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(var_name(i)),          &
     &                                      IvarS, Vtransform(ng),      &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
          CASE ('Vstretching')
            IF (.not.PRESENT(ncid)) THEN
              CALL netcdf_get_ivar (ng, IDmod, ncname, var_name(i),     &
     &                              IvarS)
            ELSE
              CALL netcdf_get_ivar (ng, IDmod, ncname, var_name(i),     &
     &                              IvarS,                              &
     &                              ncid = ncid)
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            IF (IvarS.ne.Vstretching(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(var_name(i)),          &
     &                                      IvarS, Vstretching(ng),     &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
          CASE ('hc')
            IF (.not.PRESENT(ncid)) THEN
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS)
            ELSE
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS,                              &
     &                              ncid = ncid)
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            IF (ABS(hc(ng)-FvarS).gt.RoundOff) THEN
              IF (Master) WRITE (stdout,20) TRIM(var_name(i)),          &
     &                                      FvarS, hc(ng),              &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
          CASE ('theta_s')
            IF (.not.PRESENT(ncid)) THEN
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS)
            ELSE
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS,                              &
     &                              ncid = ncid)
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            IF (ABS(theta_s(ng)-FvarS).gt.RoundOff) THEN
              IF (Master) WRITE (stdout,20) TRIM(var_name(i)),          &
     &                                      FvarS, theta_s(ng),         &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
          CASE ('theta_b')
            IF (.not.PRESENT(ncid)) THEN
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS)
            ELSE
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS,                              &
     &                              ncid = ncid)
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            IF (ABS(theta_b(ng)-FvarS).gt.RoundOff) THEN
              IF (Master) WRITE (stdout,20) TRIM(var_name(i)),          &
     &                                      FvarS, theta_b(ng),         &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
          CASE ('Tcline')
            IF (.not.PRESENT(ncid)) THEN
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS)
            ELSE
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS,                              &
     &                              ncid = ncid)
            END IF
            IF (FoundError(exit_flag, NoError,__LINE__, MyFile)) RETURN

            IF (ABS(Tcline(ng)-FvarS).gt.RoundOff) THEN
              IF (Master) WRITE (stdout,20) TRIM(var_name(i)),          &
     &                                      FvarS, Tcline(ng),          &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
#endif
#ifdef FOUR_DVAR
          CASE ('Hgamma')
            IF ((model.eq.5).or.(model.eq.10).or.(model.eq.11)) THEN
              IF (.not.find_string(var_name,n_var,'HgammaM',ic)         &
     &            .and.(model.eq.5).and.(NSA.eq.2)) THEN
                VarVal=Hgamma(2)
# ifdef ADJUST_BOUNDARY
              ELSE IF (.not.find_string(var_name,n_var,'HgammaB',ic)    &
     &                 .and.(model.eq.10)) THEN
                VarVal=Hgamma(3)
# endif
# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
              ELSE IF (.not.find_string(var_name,n_var,'HgammaF',ic)    &
     &                 .and.(model.eq.11)) THEN
                VarVal=Hgamma(4)
# endif
              ELSE                   ! Backward compatability logic
                VarVal=Hgamma(1)     ! for a single Hgamma value
              END IF
              IF (.not.PRESENT(ncid)) THEN
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarS)
              ELSE
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarS,                            &
     &                                ncid = ncid)
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (ABS(VarVal-FvarS).gt.RoundOff) THEN
                IF (Master) WRITE (stdout,20) TRIM(var_name(i)),        &
     &                                        FvarS, VarVal,            &
     &                                        TRIM(ncname)
                exit_flag=5
                EXIT
              END IF
            END IF
# ifdef WEAK_CONSTRAINT
          CASE ('HgammaM')
            IF (model.eq.5) THEN
              IF (.not.PRESENT(ncid)) THEN
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarS)
              ELSE
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarS,                            &
     &                                ncid = ncid)
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (ABS(Hgamma(2)-FvarS).gt.RoundOff) THEN
                IF (Master) WRITE (stdout,20) TRIM(var_name(i)),        &
     &                                        FvarS, Hgamma(2),         &
     &                                        TRIM(ncname)
                exit_flag=5
                EXIT
              END IF
            END IF
# endif
# ifdef ADJUST_BOUNDARY
          CASE ('HgammaB')
            IF (model.eq.10) THEN
              IF (.not.PRESENT(ncid)) THEN
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarS)
              ELSE
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarS,                            &
     &                                ncid = ncid)
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (ABS(Hgamma(3)-FvarS).gt.RoundOff) THEN
                IF (Master) WRITE (stdout,20) TRIM(var_name(i)),        &
     &                                        FvarS, Hgamma(3),         &
     &                                        TRIM(ncname)
                exit_flag=5
                EXIT
              END IF
            END IF
# endif
# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          CASE ('HgammaF')
            IF (model.eq.11) THEN
              IF (.not.PRESENT(ncid)) THEN
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarS)
              ELSE
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarS,                            &
     &                                ncid = ncid)
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (ABS(Hgamma(4)-FvarS).gt.RoundOff) THEN
                IF (Master) WRITE (stdout,20) TRIM(var_name(i)),        &
     &                                        FvarS, Hgamma(4),         &
     &                                        TRIM(ncname)
                exit_flag=5
                EXIT
              END IF
            END IF
# endif
# ifdef SOLVE3D
          CASE ('Vgamma')
            IF ((model.eq.5).or.(model.eq.10)) THEN
              IF (.not.find_string(var_name,n_var,'VgammaM',ic)         &
     &            .and.(model.eq.5).and.(NSA.eq.2)) THEN
                VarVal=Vgamma(2)
# ifdef ADJUST_BOUNDARY
              ELSE IF (.not.find_string(var_name,n_var,'VgammaB',ic)    &
     &                 .and.(model.eq.10)) THEN
                VarVal=Vgamma(3)
# endif
              ELSE                   ! Backward compatability logic
                VarVal=Vgamma(1)     ! for a single Vgamma value
              END IF
              IF (.not.PRESENT(ncid)) THEN
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarS)
              ELSE
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarS,                            &
     &                                ncid = ncid)
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (ABS(FvarS-VarVal).gt.RoundOff) THEN
                IF (Master) WRITE (stdout,20) TRIM(var_name(i)),        &
     &                                        FvarS, VarVal,            &
     &                                        TRIM(ncname)
                exit_flag=5
                EXIT
              END IF
            END IF
#  ifdef WEAK_CONSTRAINT
          CASE ('VgammaM')
            IF (model.eq.5) THEN
              IF (.not.PRESENT(ncid)) THEN
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarS)
              ELSE
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarS,                            &
     &                                ncid = ncid)
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (ABS(FvarS-Vgamma(2)).gt.RoundOff) THEN
                IF (Master) WRITE (stdout,20) TRIM(var_name(i)),        &
     &                                        FvarS, Vgamma(2),         &
     &                                        TRIM(ncname)
                exit_flag=5
                EXIT
              END IF
            END IF
#  endif
#  ifdef ADJUST_BOUNDARY
          CASE ('VgammaB')
            IF (model.eq.5) THEN
              IF (.not.PRESENT(ncid)) THEN
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarS)
              ELSE
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarS,                            &
     &                                ncid = ncid)
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (ABS(FvarS-Vgamma(3)).gt.RoundOff) THEN
                IF (Master) WRITE (stdout,20) TRIM(var_name(i)),        &
     &                                        FvarS, Vgamma(3),         &
     &                                        TRIM(ncname)
                exit_flag=5
                EXIT
              END IF
            END IF
#  endif
# endif
          CASE ('Hdecay')
            IF ((model.eq.5).or.(model.eq.11)) THEN
              Npts=UBOUND(Hdecay,DIM=2)
              IF (.not.PRESENT(ncid)) THEN
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarV,                            &
     &                                start = (/1/), total = (/Npts/))
              ELSE
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarV,                            &
     &                                ncid = ncid,                      &
     &                                start = (/1/), total = (/Npts/))
              END IF

              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (model.eq.5) THEN
                j1=1                   ! first state variable to check
# ifdef SOLVE3D
                j2=isTvar(NT(ng))      ! last  state variable to check
# else
                j2=3                   ! last  state variable to check
# endif
              ELSE IF (model.eq.11) THEN
# ifdef SOLVE3D
                j1=isTvar(NT(ng))+1    ! first state variable to check
# else
                j1=4                   ! first state variable to check
# endif
                j2=Npts                ! last  state variable to check
              END IF

              DO j=j1,j2
                IF (ABS(Hdecay(1,j,ng)-FvarV(j)).gt.RoundOff) THEN
                  text=TRIM(var_name(i))//                              &
     &                 '(1,'//TRIM(Vname(1,idSvar(j)))//')'
                  IF (Master) WRITE (stdout,20) TRIM(text),             &
     &                                          FvarV(j),               &
     &                                          Hdecay(1,j,ng),         &
     &                                          TRIM(ncname)
                  exit_flag=5
                  EXIT
                END IF
              END DO
            END IF
# ifdef SOLVE3D
          CASE ('Vdecay')
            IF ((model.eq.5).or.(model.eq.11)) THEN
              Npts=UBOUND(Vdecay,DIM=2)
              IF (.not.PRESENT(ncid)) THEN
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarV,                            &
     &                                start = (/1/), total = (/Npts/))
              ELSE
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarV,                            &
     &                                ncid = ncid,                      &
     &                                start = (/1/), total = (/Npts/))
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (model.eq.5) THEN
                j1=1                   ! first state variable to check
                j2=isTvar(NT(ng))      ! last  state variable to check
              ELSE IF (model.eq.11) THEN
                j1=isTvar(NT(ng))+1    ! first state variable to check
                j2=Npts                ! last  state variable to check
              END IF

              DO j=j1,j2
                IF (ABS(Vdecay(1,j,ng)-FvarV(j)).gt.RoundOff) THEN
                  text=TRIM(var_name(i))//                              &
     &                 '(1,'//TRIM(Vname(1,idSvar(j)))//')'
                  IF (Master) WRITE (stdout,20) TRIM(text),             &
     &                                          FvarV(j),               &
     &                                          Vdecay(1,j,ng),         &
     &                                          TRIM(ncname)
                  exit_flag=5
                  EXIT
                END IF
              END DO
            END IF
# endif
# ifdef WEAK_CONSTRAINT
          CASE ('HdecayM')
            IF ((model.eq.5).and.(NSA.eq.2)) THEN
              Npts=UBOUND(Hdecay,DIM=2)
              IF (.not.PRESENT(ncid)) THEN
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarV,                            &
     &                                start = (/1/), total = (/Npts/))
              ELSE
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarV,                            &
     &                                ncid = ncid,                      &
     &                                start = (/1/), total = (/Npts/))
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              j1=1                     ! first state variable to check
#  ifdef SOLVE3D
              j2=isTvar(NT(ng))        ! last  state variable to check
#  else
              j2=3                     ! last  state variable to check
#  endif
              DO j=j1,j2
                IF (ABS(Hdecay(NSA,j,ng)-FvarV(j)).gt.RoundOff) THEN
                  text=TRIM(var_name(i))//                              &
     &                 '(2,'//TRIM(Vname(1,idSvar(j)))//')'
                  IF (Master) WRITE (stdout,20) TRIM(text),             &
     &                                          FvarV(j),               &
     &                                          Hdecay(NSA,j,ng),       &
     &                                          TRIM(ncname)
                  exit_flag=5
                  EXIT
                END IF
              END DO
            END IF
#  ifdef SOLVE3D
          CASE ('VdecayM')
            IF ((model.eq.5).and.(NSA.eq.2)) THEN
              Npts=UBOUND(Vdecay,DIM=2)
              IF (.not.PRESENT(ncid)) THEN
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarV,                            &
     &                                start = (/1/), total = (/Npts/))
              ELSE
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarV,                            &
     &                                ncid = ncid,                      &
     &                                start = (/1/), total = (/Npts/))
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              j1=1                     ! first state variable to check
              j2=isTvar(NT(ng))        ! last  state variable to check

              DO j=j1,j2
                IF (ABS(Vdecay(NSA,j,ng)-FvarV(j)).gt.RoundOff) THEN
                  text=TRIM(var_name(i))//                              &
     &                 '(2,'//TRIM(Vname(1,idSvar(j)))//')'
                  IF (Master) WRITE (stdout,20) TRIM(text),             &
     &                                          FvarV(j),               &
     &                                          Vdecay(NSA,j,ng),       &
     &                                          TRIM(ncname)
                  exit_flag=5
                  EXIT
                END IF
              END DO
            END IF
#  endif
# endif
# ifdef ADJUST_BOUNDARY
          CASE ('HdecayB')
            IF (model.eq.10) THEN
#  ifdef SOLVE3D
              Npts=isTvar(NT(ng))
#  else
              Npts=3
#  endif
              IF (.not.PRESENT(ncid)) THEN
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarV,                            &
     &                                start = (/1,1/),                  &
     &                                total = (/Npts,4/))
              ELSE
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarV,                            &
     &                                ncid = ncid,                      &
     &                                start = (/1,1/),                  &
     &                                total = (/Npts,4/))
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              ic=0
              j1=1                     ! first state variable to check
              j2=Npts                  ! last  state variable to check

              DO ib=1,4
                DO j=j1,j2
                  ic=ic+1
                  IF (Lobc(ib,j,ng)) THEN
                    WRITE (text,"(a,'(',i1,',',a,')')")                 &
     &                          TRIM(var_name(i)), ib,                  &
     &                          TRIM(Vname(1,idSvar(j)))
                    IF (ABS(HdecayB(j,ib,ng)-FvarV(ic)).gt.             &
     &                  RoundOff) THEN
                      IF (Master) WRITE (stdout,20) TRIM(text),         &
     &                                              FvarV(ic),          &
     &                                              HdecayB(j,ib,ng),   &
     &                                              TRIM(ncname)
                      exit_flag=5
                      EXIT
                    END IF
                  END IF
                END DO
              END DO
            END IF
#  ifdef SOLVE3D
          CASE ('VdecayB')
            IF (model.eq.10) THEN
              Npts=isTvar(NT(ng))
              IF (.not.PRESENT(ncid)) THEN
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarV,                            &
     &                                start = (/1,1/),                  &
     &                                total = (/Npts,4/))
              ELSE
                CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),   &
     &                                FvarV,                            &
     &                                ncid = ncid,                      &
     &                                start = (/1,1/),                  &
     &                                total = (/Npts,4/))
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              ic=0
              j1=1                     ! first state variable to check
              j2=Npts                  ! last  state variable to check

              DO ib=1,4
                DO j=j1,j2
                  ic=ic+1
                  IF (Lobc(ib,j,ng)) THEN
                    WRITE (text,"(a,'(',i1,',',a,')')")                 &
     &                          TRIM(var_name(i)), ib,                  &
     &                          TRIM(Vname(1,idSvar(j)))
                    IF (ABS(VdecayB(j,ib,ng)-FvarV(ic)).gt.             &
     &                  RoundOff) THEN
                      IF (Master) WRITE (stdout,20) TRIM(text),         &
     &                                              FvarV(ic),          &
     &                                              VdecayB(j,ib,ng),   &
     &                                              TRIM(ncname)
                      exit_flag=5
                      EXIT
                    END IF
                  END IF
                END DO
              END DO
            END IF
#  endif
# endif
#endif
        END SELECT
      END DO

 10   FORMAT (/,' NETCDF_CHECK_VAR - inconsistent value of variable: ', &
     &        a,2x,2i5,/,20x,'in file: ',a)
 20   FORMAT (/,' NETCDF_CHECK_VAR - inconsistent value of variable: ', &
     &        a,2x,2(1pe14.6),/,20x,'in file: ',a)
!
      RETURN
      END SUBROUTINE netcdf_check_var
!
      SUBROUTINE netcdf_inq_var (ng, model, ncname, ncid, myVarName,    &
     &                           SearchVar, VarID, nVarDim, nVarAtt)
!
!=======================================================================
!                                                                      !
!  This routine inquires a NetCDF file dimensions names and values.    !
!  All the dimension information is stored in the module variables     !
!  declared above.  In addition,  if a particular variable name is     !
!  provided, this routine returns the requested information in the     !
!  optional arguments.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     myVarName    Requested variable name (string, OPTIONAL)          !
!     SearchVar    Switch used when searching a variable over          !
!                    multiple NetCDF files (logical, OPTIONAL)         !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     VarID        ID of requested variable (integer, OPTIONAL)        !
!     nVarDim      Number of variable dimensions (integer, OPTIONAL)   !
!     nVarAtt      Number of variable attributes (integer, OPTIONAL)   !
!                                                                      !
!  Other information stored in this module:                            !
!                                                                      !
!     n_dim        Number of dimensions                                !
!     n_var        Number of variables                                 !
!     n_gatt       Number of global attributes                         !
!     rec_id       Unlimited dimension ID                              !
!     var_name     Variables name (1:n_var)                            !
!     var_id       Variables ID (1:n_var)                              !
!     var_natt     Variables number of attributes (1:n_var)            !
!     var_flag     Variables flag [1=full field, -1=water points only] !
!     var_type     Variables external data type (1:n_var)              !
!     var_ndim     Variables number of dimensions (1:n_var)            !
!     var_dim      Variables dimensions ID (:,1:n_var)                 !
!                                                                      !
!  If the OPTIONAL argument myVarName is provided, the following       !
!  information for requested variable is also stored:                  !
!                                                                      !
!     n_vatt       Number of variable attributes                       !
!     n_vdim       Number of variable dimensions                       !
!     var_kind     Variable external data type                         !
!     var_Aname    Variable attribute names (1:n_vatt)                 !
!     var_Achar    Variable string attribute values (1:n_vatt)         !
!     var_Afloat   Variable float attribute values (1:n_vatt)          !
!     var_Aint     Variable integer attribute values (1:n_vatt)        !
!     var_Dids     Variable dimensions ID (1:n_vdim)                   !
!     var_Dname    Variable dimensions name (1:n_vdim)                 !
!     var_Dsize    Variable dimensions size (1:n_vdim)                 !
!                                                                      !
!     WARNING: This is information is rewritten during each CALL.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcastf, mp_bcasts
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in), optional :: myVarName
!
      logical, intent(out), optional :: SearchVar
!
      integer, intent(in), optional :: ncid

      integer, intent(out), optional :: VarID
      integer, intent(out), optional :: nVarDim
      integer, intent(out), optional :: nVarAtt
!
!  Local variable declarations.
!
      logical :: foundit, WriteError
!
      integer :: i, j, status
      integer :: att_id, my_Alen, my_Atype, my_id, my_ncid

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(4) :: ibuffer
#endif
!
      character (len=1024) :: text

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_inq_var"
!
!-----------------------------------------------------------------------
!  Inquire about the NetCDF dimensions (names and values).
!-----------------------------------------------------------------------
!
!  Initialize.
!
      n_dim=0
      n_var=0
      n_gatt=0
      rec_id=-1
      att_kind=-1
      var_id=0
      var_natt=0
      var_flag=0
      var_type=0
      var_ndim=0
      var_dim=0
      status=nf90_noerr
      DO i=1,Matts
        DO j=1,LEN(att_name(1))
          att_name(i)(j:j)=' '
        END DO
      END DO
      DO i=1,Mvars
        DO j=1,LEN(var_name(1))
          var_name(i)(j:j)=' '
        END DO
      END DO
!
!  Open file for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Inquire NetCDF file.
!
      IF (InpThread) THEN
        status=nf90_inquire(my_ncid, n_dim, n_var, n_gatt, rec_id)
        IF ((status.eq.nf90_noerr).and.(n_var.le.Mvars)) THEN
#if !defined PARALLEL_IO  && defined DISTRIBUTE
          ibuffer(1)=n_dim
          ibuffer(2)=n_var
          ibuffer(3)=n_gatt
          ibuffer(4)=rec_id
#endif
!
!  Inquire global attribute names and their external data type.
!
          DO i=1,MIN(Matts,n_gatt)
            att_id=i
            status=nf90_inq_attname(my_ncid, nf90_global, att_id,       &
     &                              att_name(i))
            IF (status.eq.nf90_noerr) THEN
              status=nf90_inquire_attribute(my_ncid, nf90_global,       &
     &                                      TRIM(att_name(i)),          &
     &                                      xtype = att_kind(i),        &
     &                                      attnum = att_id)
              IF (status.ne.nf90_noerr) THEN
                WRITE (stdout,5) i, TRIM(ncname), TRIM(SourceFile),     &
     &                           nf90_strerror(status)
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            ELSE
              WRITE (stdout,5) i, TRIM(ncname), TRIM(SourceFile),       &
     &                         nf90_strerror(status)
              exit_flag=2
              ioerror=status
              EXIT
            END IF
          END DO
!
!  Inquire about variables: name, ID, dimensions, data type, and number
!  of attributes.
!
          IF (status.eq.nf90_noerr) THEN
            DO i=1,n_var
              var_id(i)=i
              var_flag(i)=1
              status=nf90_inquire_variable(my_ncid, var_id(i),          &
     &                                     var_name(i), var_type(i),    &
     &                                     var_ndim(i), var_dim(:,i),   &
     &                                     var_natt(i))
              IF (status.eq.nf90_noerr) THEN
                DO j=1,MIN(NvarA,var_natt(i))
                  status=nf90_inq_attname(my_ncid, var_id(i), j,        &
     &                                    var_Aname(j))
                  IF (status.eq.nf90_noerr) THEN
                    IF (TRIM(var_Aname(j)).eq.'water_points'.and.       &
     &                  (var_ndim(i).gt.0)) THEN
                      var_flag(i)=-1
                    END IF
                  ELSE
                    WRITE (stdout,10) j, TRIM(var_name(i)),             &
     &                                TRIM(ncname), TRIM(SourceFile),   &
     &                                nf90_strerror(status)
                    exit_flag=2
                    ioerror=status
                    EXIT
                  END IF
                END DO
              ELSE
                WRITE (stdout,20) var_id(i), TRIM(ncname),              &
     &                            TRIM(SourceFile),                     &
     &                            nf90_strerror(status)
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            END DO
          END IF
        ELSE
          IF (n_var.gt.Mvars) THEN
            WRITE (stdout,30) 'Mvars = ', Mvars, n_var
            exit_flag=2
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast dimension to all processors in the group.
!
      CALL mp_bcasti (ng, model, exit_flag)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, ibuffer)
        n_dim=ibuffer(1)
        n_var=ibuffer(2)
        n_gatt=ibuffer(3)
        rec_id=ibuffer(4)
        CALL mp_bcasti (ng, model, att_kind)
        CALL mp_bcasti (ng, model, var_id)
        CALL mp_bcasti (ng, model, var_flag)
        CALL mp_bcasti (ng, model, var_type)
        CALL mp_bcasti (ng, model, var_ndim)
        CALL mp_bcasti (ng, model, var_natt)
        CALL mp_bcasti (ng, model, var_dim)
        CALL mp_bcasts (ng, model, att_name)
        CALL mp_bcasts (ng, model, var_name)
      END IF
#endif
!
!  Load requested requested variable information.
!
      IF (exit_flag.eq.NoError) THEN
        foundit=.FALSE.
        IF (PRESENT(myVarName)) THEN
          var_Dids=-1
          var_Dsize=0
          var_Aint=0
          var_Afloat=0.0_r8
          DO i=1,NvarA
            DO j=1,LEN(var_Aname(1))
              var_Aname(i)(j:j)=' '
            END DO
            DO j=1,LEN(var_Achar(1))
              var_Achar(i)(j:j)=' '
            END DO
          END DO
          DO i=1,NvarD
            DO j=1,LEN(var_Dname(1))
              var_Dname(i)(j:j)=' '
            END DO
          END DO
!
          DO i=1,n_var
            IF (TRIM(var_name(i)).eq.TRIM(myVarName)) THEN
              foundit=.TRUE.
              my_id=var_id(i)
              n_vatt=var_natt(i)
              n_vdim=var_ndim(i)
              DO j=1,n_vdim
                var_Dids(j)=var_dim(j,i)
              END DO
              var_kind=var_type(i)
            END IF
          END DO
          IF (foundit) THEN
            IF (PRESENT(VarID)) THEN
              VarID=my_id
            END IF
            IF (PRESENT(nVarDim)) THEN
              nVarDim=n_vdim
            END IF
            IF (PRESENT(nVarAtt)) THEN
              nVarAtt=n_vatt
            END IF
          END IF
!
!  If founded requested variable, inquire about is dimensions and
!  attributes.
!
          IF (foundit.and.InpThread) THEN
            DO i=1,n_vdim
              status=nf90_inquire_dimension(my_ncid, var_Dids(i),       &
     &                                      var_Dname(i), var_Dsize(i))
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                WRITE (stdout,50) i, TRIM(myVarName), TRIM(ncname),     &
     &                            TRIM(SourceFile),                     &
     &                            nf90_strerror(status)
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            END DO
            IF (status.eq.nf90_noerr) THEN
              DO i=1,MIN(NvarA, n_vatt)
                status=nf90_inq_attname(my_ncid, my_id, i, var_Aname(i))
                IF (status.eq.nf90_noerr) THEN
                  status=nf90_inquire_attribute(my_ncid, my_id,         &
     &                                          TRIM(var_Aname(i)),     &
     &                                          xtype = my_Atype,       &
     &                                          len = my_Alen)
                  IF (status.eq.nf90_noerr) THEN
                    IF ((my_Alen.eq.1).and.                             &
     &                  (my_Atype.eq.NF90_INT)) THEN
                      status=nf90_get_att(my_ncid, my_id,               &
     &                                    TRIM(var_Aname(i)),           &
     &                                    var_Aint(i))
                      IF (FoundError(status, nf90_noerr,                &
     &                               __LINE__, MyFile)) THEN
                        WRITE (stdout,60) 'integer',                    &
     &                                    TRIM(var_Aname(i)),           &
     &                                    TRIM(myVarName),              &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile),             &
     &                                    nf90_strerror(status)
                        exit_flag=2
                        ioerror=status
                        EXIT
                      END IF
                    ELSE IF ((my_Alen.eq.1).and.                        &
     &                       ((my_Atype.eq.NF90_FLOAT ).or.             &
     &                        (my_Atype.eq.NF90_DOUBLE))) THEN
                      status=nf90_get_att(my_ncid, my_id,               &
     &                                    TRIM(var_Aname(i)),           &
     &                                    var_Afloat(i))
                      IF (FoundError(status, nf90_noerr,                &
     &                               __LINE__, MyFile)) THEN
                        WRITE (stdout,60) 'float',                      &
     &                                    TRIM(var_Aname(i)),           &
     &                                    TRIM(myVarName),              &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile),             &
     &                                    nf90_strerror(status)
                        exit_flag=2
                        ioerror=status
                        EXIT
                      END IF
                    ELSE IF (my_Atype.eq.NF90_CHAR) THEN
                      status=nf90_get_att(my_ncid, my_id,               &
     &                                    TRIM(var_Aname(i)),           &
     &                                    text(1:my_Alen))
                      IF (FoundError(status, nf90_noerr,                &
     &                               __LINE__, MyFile)) THEN
                        WRITE (stdout,60) 'string',                     &
     &                                    TRIM(var_Aname(i)),           &
     &                                    TRIM(myVarName),              &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile),             &
     &                                    nf90_strerror(status)
                        exit_flag=2
                        ioerror=status
                        EXIT
                      END IF
                      var_Achar(i)=text(1:my_Alen)
                    END IF
                  ELSE
                    WRITE (stdout,10) i, TRIM(myVarName), TRIM(ncname), &
     &                                TRIM(SourceFile),                 &
     &                                nf90_strerror(status)
                    exit_flag=2
                    ioerror=status
                    EXIT
                  END IF
                ELSE
                  WRITE (stdout,70) i, TRIM(myVarName), TRIM(ncname),   &
     &                              TRIM(SourceFile)
                  exit_flag=4
                  ioerror=status
                  EXIT
                END IF
              END DO
            END IF
          END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast requested variable information to all processors in the
!  group.
!
          IF (foundit) THEN
            CALL mp_bcasti (ng, model, exit_flag)
            IF (exit_flag.eq.NoError) THEN
              CALL mp_bcasts (ng, model, var_Dname)
              IF (n_vdim.gt.0) THEN
                CALL mp_bcasti (ng, model, var_Dsize)
              END IF
              CALL mp_bcasts (ng, model, var_Aname)
              CALL mp_bcasts (ng, model, var_Achar)
              IF (n_vatt.gt.0) THEN
                CALL mp_bcasti (ng, model, var_Aint)
                CALL mp_bcastf (ng, model, var_Afloat)
              END IF
            END IF
          END IF
#endif
!
!  Ignore error message if requested variable not found when searching
!  over multiple input NetCDF files.
!
          IF (PRESENT(SearchVar)) THEN
            SearchVar=foundit
            WriteError=.FALSE.
          ELSE
            WriteError=.TRUE.
          END IF
          IF (.not.foundit.and.WriteError) THEN
            IF (Master) WRITE (stdout,80) TRIM(myVarName),              &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
   5  FORMAT (/,' NETCDF_INQ_VAR - error while inquiring global',       &
     &        ' attribute: ',i2.2,/,18x,'in input file:',2x,a,          &
     &        /,18x,'call from:',2x,a,/,18x,a)
  10  FORMAT (/,' NETCDF_INQ_VAR - error while inquiring attribute ',   &
     &        i1,' for variable: ',a,/,18x,'in input file:',2x,a,       &
     &        /,18x,'call from:',2x,a,/18x,a)
  20  FORMAT (/,' NETCDF_INQ_VAR - error while inquiring variable ID:', &
     &        2x,i3,/,18x,'in input file:',2x,a,/,18x,'call from:',     &
     &        2x,a,/,18x,a)
  30  FORMAT (/,' NETCDF_INQ_VAR - too small dimension parameter, ',a,  &
     &        2i5,/,18x,'change file  mod_netcdf.F  and recompile')
  40  FORMAT (/,' NETCDF_INQ_VAR - unable to inquire about contents',   &
     &        ' of input NetCDF file:',2x,a,/,18x,'call from:',2x,a,    &
     &        /,18x,a)
  50  FORMAT (/,' NETCDF_INQ_VAR - error while inquiring dimension ',   &
     &        i1,' for variable:',2x,a,/,18x,'in input file:',2x,a,     &
     &        /,18x,'call from:',2x,a,/,18x,a)
  60  FORMAT (/,' NETCDF_INQ_VAR - error while reading ',a,             &
     &        'attribute:',1x,a,' for variable ',a,/,18x,               &
     &        'in input file:',2x,a,/,18x,'call from:',2x,a,/,18x,a)
  70  FORMAT (/,' NETCDF_INQ_VAR - unable to inquire name of ',         &
     &        'attribute ',i1,' for variable ',a,/,18x,                 &
     &        'in input file:',2x,a,/,18x,'call from:',2x,a,/,18x,a)
  80  FORMAT (/,' NETCDF_INQ_VAR - requested variable:',2x,a,/18x,      &
     &        'not found in input file:',2x,a,/,18x,'call from:',2x,a,  &
     &        /,18x,a)
!
      RETURN
      END SUBROUTINE netcdf_inq_var
!
      SUBROUTINE netcdf_inq_varid (ng, model, ncname, myVarName,        &
     &                             ncid, VarID)
!
!=======================================================================
!                                                                      !
!  This routine inquires ID of requested NetCDF variable.              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Requested variable name (string)                    !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     VarID        ID of requested variable (integer)                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      integer, intent(out) :: VarID
!
!  Local variable declarations.
!
      integer :: status
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(3) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_inq_varid"
!
!-----------------------------------------------------------------------
!  Inquire ID of requested variable.
!-----------------------------------------------------------------------
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncid, TRIM(myVarName), VarID)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=3
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      ibuffer(3)=VarID
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      VarID=ibuffer(3)
#endif
!
  10  FORMAT (/,' NETCDF_INQ_VARID - error while inquiring ID for ',    &
     &        'variable:',2x,a,/,20x,'in input file:',2x,a,/,           &
     &        20x,'call from:',2x,a,/,20x,a)
!
      RETURN
      END SUBROUTINE netcdf_inq_varid
!
      SUBROUTINE netcdf_get_fatt (ng, model, ncname, varid, AttName,    &
     &                            AttValue, foundit, ncid)
!
!=======================================================================
!                                                                      !
!  This routine gets requested variable floating-point attribute(s).   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     varid        Variable ID for variable attribute or               !
!                    NF90_GLOBAL for a global attribute (integer)      !
!     AttName      Attribute name to read (string array)               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     AttValue     Attribute value (real array)                        !
!     foundit      Switch (T/F) activated when the requested           !
!                    attribute is found (logical array)                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, varid

      integer, intent(in), optional :: ncid
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: AttName(:)
!
      logical, intent(out) :: foundit(:)
!
      real(r8), intent(out) :: AttValue(:)
!
!  Local variable declarations.
!
      integer :: i, j, my_natts, my_ncid, natts, status

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
      real(r8), allocatable :: buffer(:)
#endif
!
      character (len=40) :: my_Aname
      character (len=40) :: my_Vname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_fatt"
!
!-----------------------------------------------------------------------
!  Inquire ID of requested variable.
!-----------------------------------------------------------------------
!
!  Get number of variable attributes to process.
!
      natts=UBOUND(AttName, DIM=1)
      DO i=1,natts
        foundit(i)=.FALSE.
        AttValue(i)=0.0_r8
      END DO
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      IF (.not.allocated(buffer)) THEN
        allocate ( buffer(2*natts+1) )
      END IF
#endif
!
!  If appropriate, open file for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Inquire about requested attribute value.
!
      IF (InpThread) THEN
        IF (varid.eq.nf90_global) THEN
          status=nf90_inquire(my_ncid,                                  &
     &                        nAttributes = my_natts)
        ELSE
          status=nf90_inquire_variable(my_ncid, varid,                  &
     &                                 name = my_Vname,                 &
     &                                 nAtts = my_natts)
        END IF
        IF (status.eq.nf90_noerr) THEN
          DO j=1,my_natts
            status=nf90_inq_attname(my_ncid, varid, j, my_Aname)
            IF (status.eq.nf90_noerr) THEN
              DO i=1,natts
                IF (TRIM(my_Aname).eq.TRIM(AttName(i))) THEN
                  status=nf90_get_att(my_ncid, varid, TRIM(AttName(i)), &
     &                                AttValue(i))
                  IF (FoundError(status, nf90_noerr,                    &
     &                           __LINE__, MyFile)) THEN
                    IF (Master) WRITE (stdout,10) TRIM(AttName(i)),     &
     &                                            TRIM(my_Vname),       &
     &                                            TRIM(ncname),         &
     &                                            TRIM(SourceFile),     &
     &                                            nf90_strerror(status)
                    exit_flag=2
                    ioerror=status
                  END IF
                  foundit(i)=.TRUE.
                  EXIT
                END IF
              END DO
            ELSE
              IF (Master) WRITE (stdout,20) j,                          &
     &                                      TRIM(my_Vname),             &
     &                                      TRIM(ncname),               &
     &                                      TRIM(SourceFile),           &
     &                                      nf90_strerror(status)
              exit_flag=2
              ioerror=status
              EXIT
            END IF
          END DO
#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
          buffer=0.0_r8
          DO i=1,natts
            buffer(i)=AttValue(i)
            IF (foundit(i)) THEN
               buffer(i+natts)=1.0_r8
            END IF
          END DO
          buffer(2*natts+1)=REAL(ioerror, r8)
#endif
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(my_Vname),                 &
     &                                  TRIM(ncname),                   &
     &                                  TRIM(SourceFile),               &
     &                                  nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast information to all processors in the group. Pack all
!  data into a real array for efficient communications.
!
      CALL mp_bcasti (ng, model, exit_flag)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcastf (ng, model, buffer)
        DO i=1,natts
          AttValue(i)=buffer(i)
          IF (buffer(i+natts).gt.0.0_r8) THEN
            foundit(i)=.TRUE.
          END IF
        END DO
        ioerror=INT(buffer(2*natts+1))
        IF (allocated(buffer)) THEN
          deallocate (buffer)
        END IF
      END IF
#endif
!
!  If applicable, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FATT - error while reading attribute:',    &
     &        1x,a,'for variable',1x,a,/,19x,'in input file:',2x,a,/,   &
     &        19x,'call from:',2x,a,/,19x,a)
  20  FORMAT (/,' NETCDF_GET_FATT - error while inquiring attribute:',  &
     &        1x,i2.2,'for variable',1x,a,/,19x,'in input file:',2x,a,  &
     &        /,19x,'call from:',2x,a,/,19x,a)
  30  FORMAT (/,' NETCDF_GET_FATT - error while inquiring number of ',  &
     &        'attributes for variable:',1x,a,/,19x,'in input file:',   &
     &        2x,a,/,19x,'call from:',2x,a,/,19x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_fatt
!
      SUBROUTINE netcdf_get_satt (ng, model, ncname, varid, AttName,    &
     &                            AttValue, foundit, ncid)
!
!=======================================================================
!                                                                      !
!  This routine gets requested variable string attribute(s).           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     varid        Variable ID for variable attribute or               !
!                    NF90_GLOBAL for a global attribute (integer)      !
!     AttName      Attribute name to read (string array)               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     AttValue     Attribute value (string array)                      !
!     foundit      Switch (T/F) activated when the requested           !
!                    attribute is found (logical array)                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcasts
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, varid

      integer, intent(in), optional :: ncid
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: AttName(:)
!
      logical, intent(out) :: foundit(:)
!
      character (len=*), intent(out) :: AttValue(:)
!
!  Local variable declarations.
!
      integer :: i, j, my_natts, my_ncid, natts, status

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, allocatable :: ibuffer(:)
#endif
!
      character (len=40) :: my_Aname
      character (len=40) :: my_Vname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_satt"
!
!-----------------------------------------------------------------------
!  Inquire ID of requested variable.
!-----------------------------------------------------------------------
!
!  Get number of variable attributes to process and initialize.
!
      natts=UBOUND(AttName, DIM=1)
      DO i=1,natts
        foundit(i)=.FALSE.
        AttValue(i)=' '
      END DO
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      IF (.not.allocated(ibuffer)) THEN
        allocate ( ibuffer(natts+1) )
      END IF
#endif
!
!  If appropriate, open file for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Inquire about requested attribute value.
!
      IF (InpThread) THEN
        IF (varid.eq.nf90_global) THEN
          status=nf90_inquire(my_ncid,                                  &
     &                        nAttributes = my_natts)
        ELSE
          status=nf90_inquire_variable(my_ncid, varid,                  &
     &                                 name = my_Vname,                 &
     &                                 nAtts = my_natts)
        END IF
        IF (status.eq.nf90_noerr) THEN
          DO j=1,my_natts
            status=nf90_inq_attname(my_ncid, varid, j, my_Aname)
            IF (status.eq.nf90_noerr) THEN
              DO i=1,natts
                IF (TRIM(my_Aname).eq.TRIM(AttName(i))) THEN
                  status=nf90_get_att(my_ncid, varid, TRIM(AttName(i)), &
     &                                AttValue(i))
                  IF (FoundError(status, nf90_noerr,                    &
     &                           __LINE__, MyFile)) THEN
                    IF (Master) WRITE (stdout,10) TRIM(AttName(i)),     &
     &                                            TRIM(my_Vname),       &
     &                                            TRIM(ncname),         &
     &                                            TRIM(SourceFile),     &
     &                                            nf90_strerror(status)
                    exit_flag=2
                    ioerror=status
                  END IF
                  foundit(i)=.TRUE.
                  EXIT
                END IF
              END DO
            ELSE
              IF (Master) WRITE (stdout,20) j,                          &
     &                                      TRIM(my_Vname),             &
     &                                      TRIM(ncname),               &
     &                                      TRIM(SourceFile),           &
     &                                      nf90_strerror(status)
              exit_flag=2
              ioerror=status
              EXIT
            END IF
          END DO
#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
          ibuffer=0
          DO i=1,natts
            IF (foundit(i)) THEN
              ibuffer(i)=1
            END IF
          END DO
          ibuffer(natts+1)=ioerror
#endif
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(my_Vname),                 &
     &                                  TRIM(ncname),                   &
     &                                  TRIM(SourceFile),               &
     &                                  nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast information to all processors in the group. Pack all
!  data into a real array for efficient communications.
!
      CALL mp_bcasti (ng, model, exit_flag)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, ibuffer)
        DO i=1,natts
          IF (ibuffer(i).gt.0) THEN
            foundit(i)=.TRUE.
          END IF
        END DO
        ioerror=ibuffer(natts+1)
        CALL mp_bcasts (ng, model, AttValue)
      END IF
      IF (allocated(ibuffer)) THEN
        deallocate (ibuffer)
      END IF
#endif
!
!  If applicable, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_SATT - error while reading attribute:',    &
     &        1x,a,'for variable',1x,a,/,19x,'in input file:',2x,a,/,   &
     &        19x,'call from:',2x,a,/,19x,a)
  20  FORMAT (/,' NETCDF_GET_SATT - error while inquiring attribute:',  &
     &        1x,i2.2,'for variable',1x,a,/,19x,'in input file:',2x,a,  &
     &        /,19x,'call from:',2x,a,/,19x,a)
  30  FORMAT (/,' NETCDF_GET_SATT - error while inquiring number of ',  &
     &        'attributes for variable:',1x,a,/,19x,'in input file:',   &
     &        2x,a,/,19x,'call from:',2x,a,/19x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_satt
!
      SUBROUTINE netcdf_get_fvar_0d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total, broadcast,     &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point scalar variable from    !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!     broadcast    Switch to broadcast read values from root to all    !
!                    members of the communicator in distributed-       !
!                    memory applications (logical, OPTIONAL,           !
!                    default=TRUE)                                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (real)                         !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(1))    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: broadcast
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val

      real(r8), intent(out) :: A
!
!  Local variable declarations.
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      logical :: DoBroadcast
!
#endif

      integer :: my_ncid, status, varid

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      real(r8), dimension(1) :: my_A
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_fvar_0d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Determine if the read data is only needed and allocated by the master
!  node ins distribute-memory applications.
!
      IF (.not.PRESENT(broadcast)) THEN
        DoBroadcast=.TRUE.
      ELSE
        DoBroadcast=broadcast
      END IF
#endif
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, my_A, start, total)
            A=my_A(1)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (DoBroadcast.and.(exit_flag.eq.NoError)) THEN
        CALL mp_bcastf (ng, model, A)
      END IF
#endif
!
!  Compute minimum and maximum values of read variable. Notice that
!  the same read value is assigned since a scalar variable was
!  processed.
!
      IF (PRESENT(min_val)) THEN
        min_val=A
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=A
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FVAR_0D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_FVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_fvar_0d
!
      SUBROUTINE netcdf_get_fvar_1d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total, broadcast,     &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 1D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!     broadcast    Switch to broadcast read values from root to all    !
!                    members of the communicator in distributed-       !
!                    memory applications (logical, OPTIONAL,           !
!                    default=TRUE)                                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(0:))   !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(:,1))  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: broadcast
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val

      real(r8), intent(out) :: A(:)
!
!  Local variable declarations.
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      logical :: DoBroadcast
#endif
      logical, dimension(3) :: foundit
!
      integer :: i, my_ncid, status, varid

      integer, dimension(1) :: Asize
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      real(r8) :: Afactor, Aoffset, Aspval

      real(r8), parameter :: Aepsilon = 1.0E-8_r8

      real(r8), dimension(3) :: AttValue
!
      character (len=12), dimension(3) :: AttName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_fvar_1d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point 1D-array variable.
!-----------------------------------------------------------------------
!
      IF (PRESENT(start).and.PRESENT(total)) THEN
        Asize(1)=1
        DO i=1,SIZE(total)              ! this logic is for the case
          Asize(1)=Asize(1)*total(i)    ! of reading multidimensional
        END DO                          ! data into a compact 1D array
      ELSE
        Asize(1)=UBOUND(A, DIM=1)
      END IF
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Determine if the read data is only needed and allocated by the master
!  node ins distribute-memory applications.
!
      IF (.not.PRESENT(broadcast)) THEN
        DoBroadcast=.TRUE.
      ELSE
        Dobroadcast=broadcast
      END IF
#endif
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (DoBroadcast.and.(exit_flag.eq.NoError)) THEN
        CALL mp_bcastf (ng, model, A)
      END IF
#endif
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '

      CALL netcdf_get_fatt (ng, model, ncname, varid, AttName,          &
     &                      AttValue, foundit,                          &
     &                      ncid = my_ncid)

      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=AttValue(1)
        END IF

        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=AttValue(2)
        END IF

        IF (.not.foundit(3)) THEN
          Aspval=spval_check
        ELSE
          Aspval=AttValue(3)
        END IF

        DO i=1,Asize(1)                       ! zero out missing values
          IF ((foundit(3).and.(ABS(A(i)-Aspval).lt.Aepsilon)).or.       &
     &        (.not.foundit(3).and.(ABS(A(i)).ge.ABS(Aspval)))) THEN
            A(i)=0.0_r8
          END IF
        END DO

        IF (foundit(1)) THEN                  ! scale data
          DO i=1,Asize(1)
            A(i)=Afactor*A(i)
          END DO
        END IF

        IF (foundit(2)) THEN                  ! add data offset
          DO i=1,Asize(1)
            A(i)=A(i)+Aoffset
          END DO
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FVAR_1D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_FVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_fvar_1d
!
      SUBROUTINE netcdf_get_fvar_2d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total, broadcast,     &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 2D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!     broadcast    Switch to broadcast read values from root to all    !
!                    members of the communicator in distributed-       !
!                    memory applications (logical, OPTIONAL,           !
!                    default=TRUE)                                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 2D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(0:,:)) !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(0:,0:))!
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(:,:,1))!
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: broadcast
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val

      real(r8), intent(out) :: A(:,:)
!
!  Local variable declarations.
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      logical :: Dobroadcast
#endif
      logical, dimension(3) :: foundit
!
      integer :: i, j, my_ncid, status, varid

      integer, dimension(2) :: Asize

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      real(r8) :: Afactor, Aoffset, Aspval

      real(r8), parameter :: Aepsilon = 1.0E-8_r8

      real(r8), dimension(3) :: AttValue
!
      character (len=12), dimension(3) :: AttName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_fvar_2d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point 2D-array variable.
!-----------------------------------------------------------------------
!
      IF (PRESENT(start).and.PRESENT(total)) THEN
        Asize(1)=total(1)
        Asize(2)=total(2)
      ELSE
        Asize(1)=UBOUND(A, DIM=1)
        Asize(2)=UBOUND(A, DIM=2)
      END IF
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Determine if the read data is only needed and allocated by the master
!  node ins distribute-memory applications.
!
      IF (.not.PRESENT(broadcast)) THEN
        DoBroadcast=.TRUE.
      ELSE
        DoBroadcast=broadcast
      END IF
#endif
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (DoBroadcast.and.(exit_flag.eq.NoError)) THEN
        CALL mp_bcastf (ng, model, A)
      END IF
#endif
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '

      CALL netcdf_get_fatt (ng, model, ncname, varid, AttName,          &
     &                      AttValue, foundit,                          &
     &                      ncid = my_ncid)

      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=AttValue(1)
        END IF

        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=AttValue(2)
        END IF

        IF (.not.foundit(3)) THEN
          Aspval=spval_check
        ELSE
          Aspval=AttValue(3)
        END IF

        DO j=1,Asize(2)                       ! zero out missing values
          DO i=1,Asize(1)
            IF ((foundit(3).and.(ABS(A(i,j)-Aspval).lt.Aepsilon)).or.   &
     &          (.not.foundit(3).and.(ABS(A(i,j)).ge.ABS(Aspval)))) THEN
              A(i,j)=0.0_r8
            END IF
          END DO
        END DO

        IF (foundit(1)) THEN                  ! scale data
          DO j=1,Asize(2)
            DO i=1,Asize(1)
              A(i,j)=Afactor*A(i,j)
            END DO
          END DO
        END IF

        IF (foundit(2)) THEN                  ! add data offset
          DO j=1,Asize(2)
            DO i=1,Asize(1)
              A(i,j)=A(i,j)+Aoffset
            END DO
          END DO
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FVAR_2D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_FVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_fvar_2d
!
      SUBROUTINE netcdf_get_fvar_3d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total, broadcast,     &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 3D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!     broadcast    Switch to broadcast read values from root to all    !
!                    members of the communicator in distributed-       !
!                    memory applications (logical, OPTIONAL,           !
!                    default=TRUE)                                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 3D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: broadcast
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val

      real(r8), intent(out) :: A(:,:,:)
!
!  Local variable declarations.
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      logical :: DoBroadcast
#endif
      logical, dimension(3) :: foundit
!
      integer :: i, j, k, my_ncid, status, varid

      integer, dimension(3) :: Asize

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      real(r8) :: Afactor, Aoffset, Aspval

      real(r8), parameter :: Aepsilon = 1.0E-8_r8

      real(r8), dimension(3) :: AttValue
!
      character (len=12), dimension(3) :: AttName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_fvar_3d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point 2D-array variable.
!-----------------------------------------------------------------------
!
      IF (PRESENT(start).and.PRESENT(total)) THEN
        Asize(1)=total(1)
        Asize(2)=total(2)
        Asize(3)=total(3)
      ELSE
        Asize(1)=UBOUND(A, DIM=1)
        Asize(2)=UBOUND(A, DIM=2)
        Asize(3)=UBOUND(A, DIM=3)
      END IF
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Determine if the read data is only needed and allocated by the master
!  node ins distribute-memory applications.
!
      IF (.not.PRESENT(broadcast)) THEN
        DoBroadcast=.TRUE.
      ELSE
        DoBroadcast=broadcast
      END IF
#endif
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (DoBroadcast.and.(exit_flag.eq.NoError)) THEN
        CALL mp_bcastf (ng, model, A)
      END IF
#endif
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '

      CALL netcdf_get_fatt (ng, model, ncname, varid, AttName,          &
     &                      AttValue, foundit,                          &
     &                      ncid = my_ncid)

      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=AttValue(1)
        END IF

        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=AttValue(2)
        END IF

        IF (.not.foundit(3)) THEN
          Aspval=spval_check
        ELSE
          Aspval=AttValue(3)
        END IF

        DO k=1,Asize(3)                       ! zero out missing values
          DO j=1,Asize(2)
            DO i=1,Asize(1)
              IF ((foundit(3).and.                                      &
     &             (ABS(A(i,j,k)-Aspval).lt.Aepsilon)).or.              &
     &            (.not.foundit(3).and.                                 &
     &             (ABS(A(i,j,k)).ge.ABS(Aspval)))) THEN
                A(i,j,k)=0.0_r8
              END IF
            END DO
          END DO
        END DO

        IF (foundit(1)) THEN                  ! scale data
          DO k=1,Asize(3)
            DO j=1,Asize(2)
              DO i=1,Asize(1)
                A(i,j,k)=Afactor*A(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (foundit(2)) THEN                  ! add data offset
          DO k=1,Asize(3)
            DO j=1,Asize(2)
              DO i=1,Asize(1)
                A(i,j,k)=A(i,j,k)+Aoffset
              END DO
            END DO
          END DO
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FVAR_3D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_FVAR_3D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_fvar_3d
!
      SUBROUTINE netcdf_get_fvar_4d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total, broadcast,     &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 4D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!     broadcast    Switch to broadcast read values from root to all    !
!                    members of the communicator in distributed-       !
!                    memory applications (logical, OPTIONAL,           !
!                    default=TRUE)                                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 4D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: broadcast
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val

      real(r8), intent(out) :: A(:,:,:,:)
!
!  Local variable declarations.
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      logical :: DoBroadcast
#endif
      logical, dimension(3) :: foundit
!
      integer :: i, j, k, l, my_ncid, status, varid

      integer, dimension(4) :: Asize

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      real(r8) :: Afactor, Aoffset, Aspval

      real(r8), parameter :: Aepsilon = 1.0E-8_r8

      real(r8), dimension(3) :: AttValue
!
      character (len=12), dimension(3) :: AttName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_fvar_4d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point 2D-array variable.
!-----------------------------------------------------------------------
!
      IF (PRESENT(start).and.PRESENT(total)) THEN
        Asize(1)=total(1)
        Asize(2)=total(2)
        Asize(3)=total(3)
        Asize(4)=total(4)
      ELSE
        Asize(1)=UBOUND(A, DIM=1)
        Asize(2)=UBOUND(A, DIM=2)
        Asize(3)=UBOUND(A, DIM=3)
        Asize(4)=UBOUND(A, DIM=4)
      END IF
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Determine if the read data is only needed and allocated by the master
!  node ins distribute-memory applications.
!
      IF (.not.PRESENT(broadcast)) THEN
        DoBroadcast=.TRUE.
      ELSE
        DoBroadcast=broadcast
      END IF
#endif
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (DoBroadcast.and.(exit_flag.eq.NoError)) THEN
        CALL mp_bcastf (ng, model, A)
      END IF
#endif
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '

      CALL netcdf_get_fatt (ng, model, ncname, varid, AttName,          &
     &                      AttValue, foundit,                          &
     &                      ncid = my_ncid)

      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=AttValue(1)
        END IF

        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=AttValue(2)
        END IF

        IF (.not.foundit(3)) THEN
          Aspval=spval_check
        ELSE
          Aspval=AttValue(3)
        END IF

        DO l=1,Asize(4)                       ! zero out missing values
          DO k=1,Asize(3)
            DO j=1,Asize(2)
              DO i=1,Asize(1)
                IF ((foundit(3).and.                                    &
     &               (ABS(A(i,j,k,l)-Aspval).lt.Aepsilon)).or.          &
     &              (.not.foundit(3).and.                               &
     &               (ABS(A(i,j,k,l)).ge.ABS(Aspval)))) THEN
                  A(i,j,k,l)=0.0_r8
                END IF
              END DO
            END DO
          END DO
        END DO

        IF (foundit(1)) THEN                  ! scale data
          DO l=1,Asize(4)
            DO k=1,Asize(3)
              DO j=1,Asize(2)
                DO i=1,Asize(1)
                  A(i,j,k,l)=Afactor*A(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF

        IF (foundit(2)) THEN                  ! add data offset
          DO l=1,Asize(4)
            DO k=1,Asize(3)
              DO j=1,Asize(2)
                DO i=1,Asize(1)
                  A(i,j,k,l)=A(i,j,k,l)+Aoffset
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FVAR_4D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_FVAR_4D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_fvar_4d
!
      SUBROUTINE netcdf_get_lvar_0d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested logical scalar variable from specified !
!  NetCDF file.  The variable can be stored as an interger (0 or 1) or !
!  as a character ('T' or 'F').   Reading a character variable is very !
!  inefficient in parallel I/O.                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (logical)                      !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_lvar (ng, iNLM, 'file.nc', 'VarName', lvar)       !
!    CALL netcdf_get_lvar (ng, iNLM, 'file.nc', 'VarName', lvar(1))    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcastl
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      logical, intent(out) :: A
!
!  Local variable declarations.
!
      integer :: my_ncid, my_type, status, varid

      integer :: AI
      integer, dimension(1) :: my_AI

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=1) :: Achar

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_lvar_0d"
!
!-----------------------------------------------------------------------
!  Read in an integer scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          status=nf90_inquire_variable(my_ncid, varid,                  &
     &                                 xtype = my_type)
          IF (status.eq.nf90_noerr) THEN
            IF (my_type.eq.nf90_int) THEN
              IF (PRESENT(start).and.PRESENT(total)) THEN
                status=nf90_get_var(my_ncid, varid, my_AI, start, total)
                AI=my_AI(1)
              ELSE
                status=nf90_get_var(my_ncid, varid, AI)
              END IF
              IF (status.eq.nf90_noerr) THEN
                IF (AI.eq.0) THEN
                  A=.FALSE.
                ELSE
                  A=.TRUE.
                END IF
              END IF
            ELSE IF (my_type.eq.nf90_char) THEN
              IF (PRESENT(start).and.PRESENT(total)) THEN
                status=nf90_get_var(my_ncid, varid, Achar, start, total)
              ELSE
                status=nf90_get_var(my_ncid, varid, Achar)
              END IF
              IF (status.eq.nf90_noerr) THEN
                A=.FALSE.
                IF ((Achar.eq.'t').or.(Achar.eq.'T')) THEN
                  A=.TRUE.
                END IF
              END IF
            END IF
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile), nf90_strerror(status)
              exit_flag=2
              ioerror=status
            END IF
          ELSE
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcastl (ng, model, A)
      END IF
#endif
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_LVAR_0D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_LVAR_0D - error while inquiring type for ',&
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  30  FORMAT (/,' NETCDF_GET_LVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_lvar_0d
!
      SUBROUTINE netcdf_get_lvar_1d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested logical 1D-array variable from         !
!  specified  NetCDF file.  The  variable  can be stored as an         !
!  interger (0 or 1) or  as a character ('T' or 'F').  Reading         !
!  a character variable is very inefficient in parallel I/O.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable (logical)                    !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_lvar (ng, iNLM, 'file.nc', 'VarName', lvar)       !
!    CALL netcdf_get_lvar (ng, iNLM, 'file.nc', 'VarName', lvar(0:))   !
!    CALL netcdf_get_lvar (ng, iNLM, 'file.nc', 'VarName', lvar(:,1))  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcastl
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      logical, intent(out) :: A(:)
!
!  Local variable declarations.
!
      integer :: i, my_ncid, my_type, status, varid

      integer, dimension(SIZE(A,1)) :: AI

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=1), dimension(SIZE(A,1)) :: Achar

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_lvar_1d"
!
!-----------------------------------------------------------------------
!  Read in an integer scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          status=nf90_inquire_variable(my_ncid, varid,                  &
     &                                 xtype = my_type)
          IF (status.eq.nf90_noerr) THEN
            IF (my_type.eq.nf90_int) THEN
              IF (PRESENT(start).and.PRESENT(total)) THEN
                status=nf90_get_var(my_ncid, varid, AI, start, total)
              ELSE
                status=nf90_get_var(my_ncid, varid, AI)
              END IF
              IF (status.eq.nf90_noerr) THEN
                DO i=1,SIZE(A,1)
                  IF (AI(i).eq.0) THEN
                    A(i)=.FALSE.
                  ELSE
                    A(i)=.TRUE.
                  END IF
                END DO
              END IF
            ELSE IF (my_type.eq.nf90_char) THEN
              IF (PRESENT(start).and.PRESENT(total)) THEN
                status=nf90_get_var(my_ncid, varid, Achar, start, total)
              ELSE
                status=nf90_get_var(my_ncid, varid, Achar)
              END IF
              IF (status.eq.nf90_noerr) THEN
                DO i=1,SIZE(A,1)
                  A(i)=.FALSE.
                  IF ((Achar(i).eq.'t').or.(Achar(i).eq.'T')) THEN
                    A(i)=.TRUE.
                  END IF
                END DO
              END IF
            END IF
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile), nf90_strerror(status)
              exit_flag=2
              ioerror=status
            END IF
          ELSE
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcastl (ng, model, A)
      END IF
#endif
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_LVAR_1D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_LVAR_1D - error while inquiring type for ',&
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  30  FORMAT (/,' NETCDF_GET_LVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_lvar_1d
!
      SUBROUTINE netcdf_get_ivar_0d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested integer scalar variable from specified !
!  NetCDF file.                                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (integer)                      !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar)       !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(1))    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      integer, intent(out) :: A
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

      integer, dimension(1) :: my_A

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_ivar_0d"
!
!-----------------------------------------------------------------------
!  Read in an integer scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, my_A, start, total)
            A=my_A(1)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, A)
      END IF
#endif
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_IVAR_0D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_IVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_ivar_0d
!
      SUBROUTINE netcdf_get_ivar_1d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested integer 1D-array variable from         !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable (integer)                    !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar)       !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(0:))   !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(:,1))  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      integer, intent(out) :: A(:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_ivar_1d"
!
!-----------------------------------------------------------------------
!  Read in an integer 1D-array variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, A)
      END IF
#endif
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_IVAR_1D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_IVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_ivar_1d
!
      SUBROUTINE netcdf_get_ivar_2d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested integer 2D-array variable from         !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 2D-array variable (integer)                    !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar)       !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(0:,:)) !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(0:,0:))!
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(:,:,1))!
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      integer, intent(out) :: A(:,:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_ivar_2d"
!
!-----------------------------------------------------------------------
!  Read in an integer 2D-array variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, A)
      END IF
#endif
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_IVAR_2D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_IVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_ivar_2d
!
      SUBROUTINE netcdf_get_svar_0d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested string scalar variable from specified  !
!  NetCDF file. The CDL of the scalar variable has one-dimension in    !
!  the NetCDF file for the number of characters:                       !
!                                                                      !
!     char string(Nchars)                                         CDL  !
!                                                                      !
!     character (len=Nchars) :: string                            F90  !
!                                                                      !
!  to read a scalar string use:                                        !
!                                                                      !
!     start = (/1/)                                                    !
!     total = (/Nchars/)                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (string)                       !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_svar (ng, iNLM, 'file.nc', 'VarName', ivar)       !
!    CALL netcdf_get_svar (ng, iNLM, 'file.nc', 'VarName', ivar(1))    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcasts
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      character (len=*), intent(out) :: A
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_svar_0d"
!
!-----------------------------------------------------------------------
!  Read in a string scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasts (ng, model, A)
      END IF
#endif
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_SVAR_0D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_SVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_svar_0d
!
      SUBROUTINE netcdf_get_svar_1d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested string 1D-array variable or array      !
!  element from specified NetCDF file. The CDL of the 1D-array         !
!  variable has two-dimensions in the NetCDF file, and the first       !
!  dimension is the number of characters:                              !
!                                                                      !
!     char string(dim1, Nchars)                                   CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1)                      F90  !
!                                                                      !
!  to read a single array element at location (i) use:                 !
!                                                                      !
!     start = (/1, i/)                                                 !
!     total = (/Nchars, 1/)                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable or array element (string)    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcasts
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      character (len=*), intent(out) :: A(:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_svar_1d"
!
!-----------------------------------------------------------------------
!  Read in a string 1D-array or array element.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasts (ng, model, A)
      END IF
#endif
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_SVAR_1D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_SVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_svar_1d
!
      SUBROUTINE netcdf_get_svar_2d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested string 2D-array variable or array      !
!  element from specified NetCDF file. The CDL of the 2D-array         !
!  variable has three-dimensions in the NetCDF file, and the first     !
!  dimension is the number of characters:                              !
!                                                                      !
!     char string(dim2, dim1, Nchars)                             CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1, dim2)                F90  !
!                                                                      !
!  to read a single array element at location (i,j) use:               !
!                                                                      !
!     start = (/1, i, j/)                                              !
!     total = (/Nchars, 1, 1/)                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (3D vector integer, OPTIONAL)        !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (3D vector integer, OPTIONAL)           !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 2D-array variable or array element (string)    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcasts
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      character (len=*), intent(out) :: A(:,:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_svar_2d"
!
!-----------------------------------------------------------------------
!  Read in a string 2D-array or array element.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasts (ng, model, A)
      END IF
#endif
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_SVAR_2D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_SVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_svar_2d
!
      SUBROUTINE netcdf_get_svar_3d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested string 3D-array variable or array      !
!  element from specified NetCDF file. The CDL of the 3D-array         !
!  variable has four-dimensions in the NetCDF file, and the first      !
!  dimension is the number of characters:                              !
!                                                                      !
!     char string(dim3, dim2, dim1, Nchars)                       CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1, dim2, dim3)          F90  !
!                                                                      !
!  to write a single array element at location (i,j,k) use:            !
!                                                                      !
!     start = (/1, i, j, k/)                                           !
!     total = (/Nchars, 1, 1, 1/)                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (4D vector      !
!                    integer, OPTIONAL)                                !
!     total        Number of data values to be read along each         !
!                    dimension (3D vector integer, OPTIONAL)           !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 3D-array variable or element (string)          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcasts
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      character (len=*), intent(out) :: A(:,:,:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_svar_3d"
!
!-----------------------------------------------------------------------
!  Read in a string 3D-array or array element.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasts (ng, model, A)
      END IF
#endif
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_SVAR_3D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_SVAR_3D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_svar_3d
!
      SUBROUTINE netcdf_get_time_0d (ng, model, ncname, myVarName,      &
     &                               Rdate, A,                          &
     &                               ncid, start, total,                &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested time scalar variable from specified    !
!  NetCDF file. If the "units" attribute of the form:                  !
!                                                                      !
!    'time-units since YYYY-MM-DD hh:mm:ss'                            !
!                                                                      !
!  is different than provided reference date "Rdate", it converts to   !
!  elapsed time since "Rdate".                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     Rdate        Reference date (real; [1] seconds, [2] days)        !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (real)                         !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(1))    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE dateclock_mod,  ONLY : datenum, datestr, time_units
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError, lowercase
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(dp), intent(in) :: Rdate(2)

      real(dp), intent(out), optional :: min_val
      real(dp), intent(out), optional :: max_val

      real(dp), intent(out) :: A
!
!  Local variable declarations.
!
      logical :: JulianOffset = .FALSE.
      logical :: Ldebug = .FALSE.

      logical, dimension(1) :: got_units
      logical, dimension(2) :: foundit
!
      integer :: ind, lstr, my_ncid, status, varid
      integer :: year, month, day, hour, minutes

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
      real(dp) :: Afactor, Aoffset, my_Rdate(2), seconds
      real(dp) :: dnum_old, dnum_new

      real(dp), dimension(1) :: my_A
      real(r8), dimension(2) :: AttValue
!
      character (len=12) :: AttName(2)
      character (len=22) :: dstr_old, dstr_new
      character (len=40) :: UnitsAtt(1), UnitsValue(1)
      character (len=40) :: Units, Ustring

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_time_0d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, my_A, start, total)
            A=my_A(1)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcastf (ng, model, A)
      END IF
#endif
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '

      CALL netcdf_get_fatt (ng, model, ncname, varid, AttName,          &
     &                      AttValue, foundit,                          &
     &                      ncid = my_ncid)

      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=REAL(AttValue(1),dp)
        END IF

        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=REAL(AttValue(2),dp)
        END IF

        IF (foundit(1)) THEN                  ! scale data
          A=Afactor*A
        END IF

        IF (foundit(2)) THEN                  ! add data offset
          A=A+Aoffset
          IF (time_ref.eq.-2) JulianOffset=.TRUE.
        END IF
      END IF
!
!  Get time variable "units" attribute and convert to elapsed time
!  since reference date. If Julian Day Number (days or seconds) and
!  'add_offset' attribute,
!
      UnitsAtt(1)='units'

      CALL netcdf_get_satt (ng, model, ncname, varid, UnitsAtt,         &
     &                      UnitsValue, got_units,                      &
     &                      ncid = my_ncid)

      IF (exit_flag.eq.NoError) THEN
        IF (got_units(1)) THEN
          Units=TRIM(lowercase(UnitsValue(1)))
          lstr=LEN_TRIM(Units)
          ind=INDEX(Units,'since')
          IF (ind.gt.0) THEN
            CALL time_units (TRIM(Units), year, month, day, hour,       &
     &                       minutes, seconds)
            CALL datenum (my_Rdate, year, month, day, hour, minutes,    &
     &                    seconds)
            IF (Rdate(1).ne.my_Rdate(1)) THEN
              Ustring=Units(1:ind-2)
              SELECT CASE (TRIM(Ustring))
                CASE ('second', 'seconds')
                  IF (Ldebug) THEN
                    IF (JulianOffset) THEN
                      dnum_old=A
                    ELSE
                      dnum_old=my_Rdate(2)+A
                    END IF
                    CALL datestr (dnum_old, .FALSE., dstr_old)
                  END IF
                  IF (JulianOffset) THEN
                    A=A-Rdate(2)              ! 'add_offset' added above
                  ELSE
                    A=(my_Rdate(2)+A)-Rdate(2)
                  END IF
                  IF (Ldebug) THEN
                    dnum_new=Rdate(2)+A
                    CALL datestr (dnum_new, .FALSE., dstr_new)
                  END IF
                CASE ('day', 'days')
                  IF (Ldebug) THEN
                    IF (JulianOffset) THEN
                      dnum_old=A
                    ELSE
                      dnum_old=my_Rdate(1)+A
                    END IF
                    CALL datestr (dnum_old, .TRUE., dstr_old)
                  END IF
                  IF (JulianOffset) THEN
                    A=A-Rdate(1)              ! 'add_offset' added above
                  ELSE
                    A=(my_Rdate(1)+A)-Rdate(1)
                  END IF
                  IF (Ldebug) THEN
                    dnum_new=Rdate(1)+A
                    CALL datestr (dnum_new, .TRUE., dstr_new)
                  END IF
              END SELECT
            END IF
          END IF
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable. Notice that
!  the same read value is assigned since a scalar variable was
!  processed.
!
      IF (PRESENT(min_val)) THEN
        min_val=A
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=A
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_TIME_0D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_TIME_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_time_0d
!
      SUBROUTINE netcdf_get_time_1d (ng, model, ncname, myVarName,      &
     &                               Rdate, A,                          &
     &                               ncid, start, total,                &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested time 1D-array variable from specified  !
!  NetCDF file.  If the "units" attribute of the form:                 !
!                                                                      !
!    'time-units since YYYY-MM-DD hh:mm:ss'                            !
!                                                                      !
!  is different than provided reference date "Rdate", it converts to   !
!  elapsed time since "Rdate".                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    time variable name (string)                         !
!     Rdate        Reference date (real; [1] seconds, [2] days)        !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array time variable (real)                  !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(0:))   !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(:,1))  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE dateclock_mod,  ONLY : datenum, datestr, time_units
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError, lowercase
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(dp), intent(in) :: Rdate(2)

      real(dp), intent(out), optional :: min_val
      real(dp), intent(out), optional :: max_val

      real(dp), intent(out) :: A(:)
!
!  Local variable declarations.
!
      logical :: JulianOffset = .FALSE.
      logical :: Ldebug = .FALSE.

      logical, dimension(1) :: got_units
      logical, dimension(2) :: foundit
!
      integer :: i, ind, lstr, my_ncid, status, varid
      integer :: year, month, day, hour, minutes

      integer, dimension(1) :: Asize
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      real(dp) :: Afactor, Aoffset, my_Rdate(2), seconds
      real(dp) :: dnum_old, dnum_new

      real(r8), dimension(2) :: AttValue
!
      character (len=12) :: AttName(2)
      character (len=22) :: dstr_old, dstr_new
      character (len=40) :: UnitsAtt(1), UnitsValue(1)
      character (len=40) :: Units, Ustring

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_get_time_1d"
!
!-----------------------------------------------------------------------
!  Read in a time 1D-array variable.
!-----------------------------------------------------------------------
!
      IF (PRESENT(start).and.PRESENT(total)) THEN
        Asize(1)=1
        DO i=1,SIZE(total)              ! this logic is for the case
          Asize(1)=Asize(1)*total(i)    ! of reading multidimensional
        END DO                          ! data into a compact 1D array
      ELSE
        Asize(1)=UBOUND(A, DIM=1)
      END IF
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in time variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile), nf90_strerror(status)
          exit_flag=2
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcastf (ng, model, A)
      END IF
#endif
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '

      CALL netcdf_get_fatt (ng, model, ncname, varid, AttName,          &
     &                      AttValue, foundit,                          &
     &                      ncid = my_ncid)

      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=REAL(AttValue(1),dp)
        END IF

        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=REAL(AttValue(2),dp)
        END IF

        IF (foundit(1)) THEN                  ! scale data
          DO i=1,Asize(1)
            A(i)=Afactor*A(i)
          END DO
        END IF

        IF (foundit(2)) THEN                  ! add data offset
          DO i=1,Asize(1)
            A(i)=A(i)+Aoffset
          END DO
          IF (time_ref.eq.-2) JulianOffset=.TRUE.
        END IF
      END IF
!
!  Get time variable "units" attribute and convert to elapsed time
!  since reference date.
!
      UnitsAtt(1)='units'

      CALL netcdf_get_satt (ng, model, ncname, varid, UnitsAtt,         &
     &                      UnitsValue, got_units,                      &
     &                      ncid = my_ncid)

      IF (exit_flag.eq.NoError) THEN
        IF (got_units(1)) THEN
          Units=TRIM(lowercase(UnitsValue(1)))
          lstr=LEN_TRIM(Units)
          ind=INDEX(Units,'since')
          IF (ind.gt.0) THEN
            CALL time_units (TRIM(Units), year, month, day, hour,       &
     &                       minutes, seconds)
            CALL datenum (my_Rdate, year, month, day, hour, minutes,    &
     &                    seconds)
            IF (Rdate(1).ne.my_Rdate(1)) THEN
              Ustring=Units(1:ind-2)
              SELECT CASE (TRIM(Ustring))
                CASE ('second', 'seconds')
                  IF (Ldebug) THEN
                    IF (JulianOffset) THEN
                      dnum_old=A(1)
                    ELSE
                      dnum_old=my_Rdate(2)+A(1)
                    END IF
                    CALL datestr (dnum_old, .FALSE., dstr_old)
                  END IF
                  IF (JulianOffset) THEN
                    DO i=1,Asize(1)
                      A(i)=A(i)-Rdate(2)      ! 'add_offset' added above
                    END DO
                  ELSE
                    DO i=1,Asize(1)
                      A(i)=(my_Rdate(2)+A(i))-Rdate(2)
                    END DO
                  END IF
                  IF (Ldebug) THEN
                    dnum_new=Rdate(2)+A(1)
                    CALL datestr (dnum_new, .FALSE., dstr_new)
                  END IF
                CASE ('day', 'days')
                  IF (Ldebug) THEN
                    IF (JulianOffset) THEN
                      dnum_old=A(1)
                    ELSE
                      dnum_old=my_Rdate(1)+A(1)
                    END IF
                    CALL datestr (dnum_old, .TRUE., dstr_old)
                  END IF
                  IF (JulianOffset) THEN
                    DO i=1,Asize(1)
                      A(i)=A(i)-Rdate(1)      ! 'add_offset' added above
                    END DO
                  ELSE
                    DO i=1,Asize(1)
                      A(i)=(my_Rdate(1)+A(i))-Rdate(1)
                    END DO
                  END IF
                  IF (Ldebug) THEN
                    dnum_new=Rdate(1)+A(1)
                    CALL datestr (dnum_new, .TRUE., dstr_new)
                  END IF
              END SELECT
            END IF
          END IF
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_TIME_1D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
  20  FORMAT (/,' NETCDF_GET_TIME_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_get_time_1d

#ifdef SINGLE_PRECISION
!
      SUBROUTINE netcdf_put_fvar_0dp (ng, model, ncname, myVarName, A,  &
     &                                start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a double-precision floating-point scalar        !
!  variable into a NetCDF file. If the NetCDF ID is not provided,      !
!  it opens the file, writes data, and then closes the file.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (double precision)      !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
# if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid
!
      real(dp), intent(in) :: A
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

# if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
# endif
!
      real(r8), dimension(1) :: my_A
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_fvar_0dp"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((start(1).eq.0).and.(total(1).eq.0)) THEN
            status=nf90_put_var(my_ncid, my_varid, A)
          ELSE
            my_A(1)=A
            status=nf90_put_var(my_ncid, my_varid, my_A, start, total)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

# if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
# endif
!
!  Close input file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_0DP - error while inquiring ID for ', &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_0DP - error while writing variable:', &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_fvar_0dp
!
      SUBROUTINE netcdf_put_fvar_1dp (ng, model, ncname, myVarName, A,  &
     &                                start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a double-precision floating-point 1D-array      !
!  variable into a file. If the NetCDF ID is not provided, it opens    !
!  the file, writes data, and then closes the file.                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (double precision)      !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
# if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid
!
      real(dp), intent(in) :: A(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

# if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_fvar_1dp"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

# if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
# endif
!
!  Close input file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_1DP - error while inquiring ID for ', &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_1DP - error while writing variable:', &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_fvar_1dp
!
      SUBROUTINE netcdf_put_fvar_2dp (ng, model, ncname, myVarName, A,  &
     &                                start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a double-precision floating-point 2D-array      !
!  variable into a file. If the NetCDF ID is not provided, it opens    !
!  the file, writes data, and then closes the file.                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (double precision)      !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
# if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid

      real(dp), intent(in) :: A(:,:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

# if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_fvar_2dp"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

# if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
# endif
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_2DP - error while inquiring ID for ', &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/22x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_2DP - error while writing variable:', &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_fvar_2dp
#endif
!
      SUBROUTINE netcdf_put_fvar_0d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point scalar variable into a NetCDF  !
!  file.  If the NetCDF ID is not provided, it opens the file, writes  !
!  data, and then closes the file.                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid
!
      real(r8), intent(in) :: A
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      real(r8), dimension(1) :: my_A
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_fvar_0d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((start(1).eq.0).and.(total(1).eq.0)) THEN
            status=nf90_put_var(my_ncid, my_varid, A)
          ELSE
            my_A(1)=A
            status=nf90_put_var(my_ncid, my_varid, my_A, start, total)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
!  Close input file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_0D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_fvar_0d
!
      SUBROUTINE netcdf_put_fvar_1d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point 1D-array variable into a file. !
!  If the NetCDF ID is not provided,  it opens the file,  writes data, !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid
!
      real(r8), intent(in) :: A(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_fvar_1d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
!  Close input file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_1D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_fvar_1d
!
      SUBROUTINE netcdf_put_fvar_2d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point 2D-array variable into a file. !
!  If the NetCDF ID is not provided,  it opens the file,  writes data, !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid
!
      real(r8), intent(in) :: A(:,:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_fvar_2d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_2D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_fvar_2d
!
      SUBROUTINE netcdf_put_fvar_3d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point 3D-array variable into a file. !
!  If the NetCDF ID is not provided,  it opens the file,  writes data, !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid
!
      real(r8), intent(in) :: A(:,:,:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_fvar_3d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_3D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_3D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_fvar_3d
!
      SUBROUTINE netcdf_put_fvar_4d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point 4D-array variable into a file. !
!  If the NetCDF ID is not provided,  it opens the file,  writes data, !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid
!
      real(r8), intent(in) :: A(:,:,:,:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_fvar_4d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_4D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_4D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_fvar_4d
!
      SUBROUTINE netcdf_put_ivar_0d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a integer scalar variable into a NetCDF file.   !
!  If the NetCDF ID is not provided, it opens the file, writes data,   !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (integer)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid

      integer, intent(in) :: A
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

      integer, dimension(1) :: my_A

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_ivar_0d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((start(1).eq.0).and.(total(1).eq.0)) THEN
            status=nf90_put_var(my_ncid, my_varid, A)
          ELSE
            my_A(1)=A
            status=nf90_put_var(my_ncid, my_varid, my_A, start, total)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
!  Close file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_IVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_IVAR_0D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_ivar_0d
!
      SUBROUTINE netcdf_put_ivar_1d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a integer 1D-array variable into a file. If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (integer)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid

      integer, intent(in) :: A(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_ivar_1d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
!  Close file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_IVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_IVAR_1D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_ivar_1d
!
      SUBROUTINE netcdf_put_ivar_2d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a integer 2D-array variable into a file. If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (integer)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid

      integer, intent(in) :: A(:,:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_ivar_2d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
!  Close file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_IVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_IVAR_2D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_ivar_2d
!
      SUBROUTINE netcdf_put_lvar_0d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a logical scalar variable into a NetCDF file.   !
!  If the NetCDF ID is not provided, it opens the file, writes data,   !
!  and then closes the file.                                           !
!                                                                      !
!  The input logical data is converted to integer such that .FALSE.    !
!  is interpreted as zero, and .TRUE. is interpreted as one.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (logical)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid
!
      logical, intent(in) :: A
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

      integer :: AI
      integer, dimension(1) :: my_AI

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_lvar_0d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Convert logical data to integer: .FALSE. is interpreted as zero, and
!  .TRUE. is interpreted as one.
!
        IF (A) THEN
          AI=1
        ELSE
          AI=0
        END IF
!
!  Write out logical data as integers.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((start(1).eq.0).and.(total(1).eq.0)) THEN
            status=nf90_put_var(my_ncid, my_varid, AI)
          ELSE
            my_AI(1)=AI
            status=nf90_put_var(my_ncid, my_varid, my_AI, start, total)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
!  Close file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_LVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_LVAR_0D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_lvar_0d
!
      SUBROUTINE netcdf_put_lvar_1d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a logical 1D-array variable into a file. If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  The input logical data is converted to integer such that .FALSE.    !
!  is interpreted as zero, and .TRUE. is interpreted as one.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (logical)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid
!
      logical, intent(in) :: A(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: i, my_ncid, my_varid, status

      integer, dimension(SIZE(A,1)) :: AI

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_lvar_1d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Convert logical data to integer: .FALSE. is interpreted as zero, and
!  .TRUE. is interpreted as one.
!
        DO i=1,SIZE(A,1)
          IF (A(i)) THEN
            AI(i)=1
          ELSE
            AI(i)=0
          END IF
        END DO
!
!  Write out logical data as integers.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, AI, start, total)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
!  Close file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_LVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_LVAR_1D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_lvar_1d
!
      SUBROUTINE netcdf_put_lvar_2d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a logical 2D-array variable into a file. If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  The input logical data is converted to integer such that .FALSE.    !
!  is interpreted as zero, and .TRUE. is interpreted as one.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (logical)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid
!
      logical, intent(in) :: A(:,:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: i, j, my_ncid, my_varid, status

      integer, dimension(SIZE(A,1),SIZE(A,2)) :: AI

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_lvar_2d"
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Convert logical data to integer: .FALSE. is interpreted as zero, and
!  .TRUE. is interpreted as one.
!
        DO j=1,SIZE(A,2)
          DO i=1,SIZE(A,1)
            IF (A(i,j)) THEN
              AI(i,j)=1
            ELSE
              AI(i,j)=0
            END IF
          END DO
        END DO
!
!  Write out logical data as integers.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, AI, start, total)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
!  Close file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_LVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_LVAR_2D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_lvar_2d
!
      SUBROUTINE netcdf_put_svar_0d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a string scalar variable into a file.  If       !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file. The CDL of the scalar variable has        !
!  one-dimension in the NetCDF file for the number of characters:      !
!                                                                      !
!     char string(Nchars)                                         CDL  !
!                                                                      !
!     character (len=Nchars) :: string                            F90  !
!                                                                      !
!  to write a scalar string use:                                       !
!                                                                      !
!     start = (/1/)                                                    !
!     total = (/Nchars/)                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (string)                !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (1D vector   !
!                    integer)                                          !
!     total        Number of data values to be written along each      !
!                    dimension (1D vector integer)                     !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid
!
      character (len=*), intent(in) :: A
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_svar_0d"
!
!-----------------------------------------------------------------------
!  Write out a scalar string.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_SVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_SVAR_0D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_svar_0d
!
      SUBROUTINE netcdf_put_svar_1d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a string 1D-array variable into a file.  If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file. The CDL of the 1D-array variable has      !
!  two-dimensions in the NetCDF file, and the first dimension is       !
!  the number of characters:                                           !
!                                                                      !
!     char string(dim1, Nchars)                                   CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1)                      F90  !
!                                                                      !
!  to write a single array element at location (i) use:                !
!                                                                      !
!     start = (/1, i/)                                                 !
!     total = (/Nchars, 1/)                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (1D string array)       !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (2D vector   !
!                    integer)                                          !
!     total        Number of data values to be written along each      !
!                    dimension (2D vector integer)                     !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid
!
      character (len=*), intent(in) :: A(:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_svar_1d"
!
!-----------------------------------------------------------------------
!  Write out a string 1D array or array element.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_SVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_SVAR_1D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_svar_1d
!
      SUBROUTINE netcdf_put_svar_2d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a string 2D-array variable into a file.  If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file. The CDL of the 2D-array variable has      !
!  three-dimensions in the NetCDF file, and the first dimension is     !
!  the number of characters:                                           !
!                                                                      !
!     char string(dim2, dim1, Nchars)                             CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1, dim2)                F90  !
!                                                                      !
!  to write a single array element at location (i,j) use:              !
!                                                                      !
!     start = (/1, i, j/)                                              !
!     total = (/Nchars, 1, 1/)                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (2D string array)       !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (3D vector   !
!                    integer)                                          !
!     total        Number of data values to be written along each      !
!                    dimension (3D vector integer)                     !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid
!
      character (len=*), intent(in) :: A(:,:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_svar_2d"
!
!-----------------------------------------------------------------------
!  Write out a string array element or full 2D array.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_SVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_SVAR_2D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_svar_2d
!
      SUBROUTINE netcdf_put_svar_3d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a string 3D-array variable into a file.  If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file. The CDL of the 3D-array variable has      !
!  four-dimensions in the NetCDF file, and the first dimension is      !
!  the number of characters:                                           !
!                                                                      !
!     char string(dim3, dim2, dim1, Nchars)                       CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1, dim2, dim3)          F90  !
!                                                                      !
!  to write a single array element at location (i,j,k) use:            !
!                                                                      !
!     start = (/1, i, j, k/)                                           !
!     total = (/Nchars, 1, 1, 1/)                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (3D string array)       !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (4D vector   !
!                    integer)                                          !
!     total        Number of data values to be written along each      !
!                    dimension (4D vector integer)                     !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid
!
      character (len=*), intent(in) :: A(:,:,:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_put_svar_3d"
!
!-----------------------------------------------------------------------
!  Write out a string array element or full 3D array.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile), nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_SVAR_3D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a,/,22x,a)
  20  FORMAT (/,' NETCDF_PUT_SVAR_3D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a, &
     &        /,22x,a)
!
      RETURN
      END SUBROUTINE netcdf_put_svar_3d
!
      SUBROUTINE netcdf_close (ng, model, ncid, ncname, Lupdate)
!
!=======================================================================
!                                                                      !
!  This routine closes requested NetCDF file. If appropriate, it       !
!  also performs additional processing, like updating the global       !
!  attributes, before closing the file.                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncid         NetCDF file ID (integer)                            !
!     ncname       NetCDF file name (string, OPTIONAL)                 !
!     Lupdate      Update global attribute (logical, OPTIONAl)         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(inout) :: ncid
!
      logical, intent(in), optional :: Lupdate
!
      character (len=*), intent(in), optional :: ncname
!
!  Local variable declarations.
!
#ifdef BIOLOGY
      logical :: my_Lupdate
!
#endif
      integer :: i, status
#ifdef BIOLOGY
      integer :: AttLen, is, ie, lstr
#endif
#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(3) :: ibuffer
#endif
!
      character (len=200) :: my_ncname
#ifdef BIOLOGY
      character (len=512) :: bio_file
#endif

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_close"
!
!-----------------------------------------------------------------------
!  If open, close requested NetCDF file.
!-----------------------------------------------------------------------
!
      IF (OutThread.and.(ncid.ne.-1)) THEN
        DO i=1,LEN(my_ncname)
          my_ncname(i:i)=' '
        END DO

        IF (.not.PRESENT(ncname)) THEN
!
!  Get filename, if any. It will be nice if there is a function in
!  the NetCDF library to do this. Fortunately, the filename is
!  written as a global attribute.
!
          status=nf90_get_att(ncid, nf90_global, 'file', my_ncname)
        ELSE
          my_ncname=TRIM(ncname)
        END IF

#ifdef BIOLOGY
!
!  Determine updating value of biology header files global attribute.
!  This is only possible in output files. An error occurs in input
!  files open for reading only. This allows to use ROMS input files
!  with the "bio_file" attribute.
!
        IF (.not.PRESENT(Lupdate)) THEN
          my_Lupdate=.TRUE.
        ELSE
          my_Lupdate=Lupdate
        END IF
!
!  Update global attribute with the biology header files used.
!
        IF (my_Lupdate) THEN
          is=1
          status=nf90_inquire_attribute(ncid, nf90_global, 'bio_file',  &
     &                                  len = AttLen)
          IF (status.eq.nf90_noerr) THEN
            DO i=1,512
              bio_file(i:i)=' '
            END DO
            DO i=1,4
              lstr=LEN_TRIM(BIONAME(i))
              IF (lstr.gt.0) THEN
                ie=is+lstr-1
                bio_file(is:ie)=TRIM(BIONAME(i))
                is=ie+1
                bio_file(is:is)=','
                is=is+2
              END IF
            END DO
            lstr=LEN_TRIM(bio_file)-1
            IF (lstr.gt.0) THEN
              status=nf90_put_att(ncid, nf90_global, 'bio_file',        &
     &                            bio_file(1:lstr))
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                WRITE (stdout,10) 'bio_file', TRIM(my_ncname),          &
     &                            TRIM(SourceFile),                     &
     &                            nf90_strerror(status)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF
        END IF
#endif
!
!  Close requested NetCDF file.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_close(ncid)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,20) ncid, TRIM(my_ncname), TRIM(SourceFile),  &
     &                        nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
          IF (Ldebug_ncid) THEN
            WRITE (DBout,'(a,1x," <= ",i8,2(2x,a))') kernel(model)//    &
     &                                               ' CLOSE', ncid,    &
     &                                               TRIM(my_ncname),   &
     &                                               TRIM(SourceFile)
            CALL my_flush (DBout)
          END IF
          ncid=-1
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      ibuffer(3)=ncid
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      ncid=ibuffer(3)
#endif
!
  10  FORMAT (/,' NETCDF_CLOSE - error while writing global ',          &
     &        'attribute:',2x,a,/,16x,'file:',2x,a,/,16x,               &
     &        'call from:',2x,a,/,16x,a)
  20  FORMAT (/,' NETCDF_CLOSE - error during closing of file, ',       &
     &        'ncid = ',i3,/,16x,'file:',2x,a,/,16x,'call from:',2x,a,  &
     &        /,16x,a)
!
      RETURN
      END SUBROUTINE netcdf_close
!
      SUBROUTINE netcdf_create (ng, model, ncname, ncid)
!
!=======================================================================
!                                                                      !
!  This routine creates a new NetCDF file. Currently, it only creates  !
!  file for serial or parallel I/O access.                             !
!                                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       Name of the new NetCDF file (string)                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(out) :: ncid
!
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: my_cmode, status
      integer :: OldFillMode
#if !defined PARALLEL_IO && defined DISTRIBUTE
      integer :: ibuffer(3)
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_create"
!
!-----------------------------------------------------------------------
!  Create requested NetCDF file.
!-----------------------------------------------------------------------
!
#if defined PARALLEL_IO && defined DISTRIBUTE

!  Create a netCDF-4/HDF5  format file. Since nf90_clobber=0, then
!  OR(nf90_clobber, nf90_netcdf4)=nf90_netcdf.  Therefore, the default
!  is to overwrite an existing file on disk with the same file name.
!  This is what we usually do anyway.
!
      my_cmode=IOR(CMODE, nf90_mpiio)
      my_cmode=IOR(my_cmode, nf90_clobber)
      status=nf90_create(path = TRIM(ncname),                           &
     &                   cmode = my_cmode,                              &
     &                   ncid = ncid,                                   &
     &                   comm = OCN_COMM_WORLD,                         &
     &                   info = MP_INFO)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile),   &
     &                                nf90_strerror(status)
        exit_flag=3
        ioerror=status
      END IF
!
!  Set NOFILL mode to enhance performance.
!
      IF (exit_flag.eq.NoError) THEN
        status=nf90_set_fill(ncid, nf90_nofill, OldFillMode)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile), &
     &                                  nf90_strerror(status)
          exit_flag=3
          ioerror=status
        END IF
      END IF
#else
      my_cmode=IOR(nf90_clobber, CMODE)
# if defined DISTRIBUTE        && !defined HDF5       && \
    (defined CONCURRENT_KERNEL || defined DISJOINTED)
      my_cmode=IOR(my_cmode, nf90_share)
# endif
      IF (OutThread) THEN
        status=nf90_create(TRIM(ncname), my_cmode, ncid)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile),             &
     &                      nf90_strerror(status)
          exit_flag=3
          ioerror=status
        END IF
        IF (Ldebug_ncid) THEN
            WRITE (DBout,'(a," <= ",i8,2(2x,a))') kernel(model)//       &
     &                                            ' CREATE', ncid,      &
     &                                            TRIM(ncname),         &
     &                                            TRIM(SourceFile)
          CALL my_flush (DBout)
        END IF
!
!  Set NOFILL mode to enhance performance.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_set_fill(ncid, nf90_nofill, OldFillMode)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) WRITE (stdout,10) TRIM(ncname),                 &
     &                                    TRIM(SourceFile),             &
     &                                    nf90_strerror(status)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
# ifdef DISTRIBUTE
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      ibuffer(3)=ncid
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      ncid=ibuffer(3)
# endif
#endif
!
  10  FORMAT (/,' NETCDF_CREATE - unable to create output NetCDF ',     &
#if defined PARALLEL_IO && defined DISTRIBUTE
     &        'parallel I/O file:',/,17x,a,/,17x,'call from:',2x,a,     &
     &        /,17x,a)
#else
     &        'file:',/,17x,a,/,17x,'call from:',2x,a,/,17x,a)
#endif
!
      RETURN
      END SUBROUTINE netcdf_create
!
      SUBROUTINE netcdf_enddef (ng, model, ncname, ncid)
!
!=======================================================================
!                                                                      !
!  This routine ends definition in an  open NetCDF  dataset.  The      !
!  changes made in define mode  are checked and committed to disk      !
!  if no errors occurred. The dataset is then placed in data mode,     !
!  so variable data can be read or written.                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       Name of the new NetCDF file (string)                !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!  Note:   There is not need to call this function when processing     !
!          NetCDF-4 format type file.                                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
!
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: status

#if !defined PARALLEL_IO && defined DISTRIBUTE
      integer :: ibuffer(2)
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_enddef"
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
#if defined PARALLEL_IO && defined DISTRIBUTE

! There is not need to call this function when processing NetCDF-4
! format type file.
!
#else
      IF (OutThread) THEN
        status=nf90_enddef(ncid)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile), &
     &                                  nf90_strerror(status)
          exit_flag=3
          ioerror=status
        END IF
      END IF

# ifdef DISTRIBUTE
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
# endif
!
  10  FORMAT (/,' NETCDF_ENDDEF - unable to end definition mode for',   &
     &        ' file:',/,17x,a,/,17x,'call from:',2x,a,/,17x,a)
#endif
!
      RETURN
      END SUBROUTINE netcdf_enddef
!
      SUBROUTINE netcdf_open (ng, model, ncname, omode, ncid)
!
!=======================================================================
!                                                                      !
!  This routine opens an existing NetCDF file for access. Currently,   !
!  it only opens file for serial or parallel I/O access.               !
!                                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       Name of the new NetCDF file (string)                !
!     omode        Open mode flag:                                     !
!                    omode = 0,  read-only access (NF90_NOWRITE)       !
!                    omode = 1,  read and write access (NF90_WRITE)    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, omode

      integer, intent(out) :: ncid
!
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: my_omode, status

#if !defined PARALLEL_IO && defined DISTRIBUTE
      integer :: ibuffer(3)
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_open"
!
!-----------------------------------------------------------------------
!  Create requested NetCDF file.
!-----------------------------------------------------------------------
!
#if defined PARALLEL_IO && defined DISTRIBUTE
      SELECT CASE (omode)
        CASE (0)
# ifdef NF_OPEN_SERIAL
          status=nf90_open(TRIM(ncname), nf90_nowrite, ncid)
# else
          status=nf90_open(path = TRIM(ncname),                         &
     &                     mode = IOR(nf90_nowrite, nf90_mpiio),        &
     &                     ncid = ncid,                                 &
     &                     comm = OCN_COMM_WORLD,                       &
     &                     info = MP_INFO)
# endif
        CASE (1)
# ifdef NF_OPEN_SERIAL
          status=nf90_open(TRIM(ncname), nf90_write, ncid)
# else
          status=nf90_open(path = TRIM(ncname),                         &
     &                     mode = IOR(nf90_write, nf90_mpiio),          &
     &                     ncid = ncid,                                 &
     &                     comm = OCN_COMM_WORLD,                       &
     &                     info = MP_INFO)
# endif
        CASE DEFAULT
# ifdef NF_OPEN_SERIAL
          status=nf90_open(TRIM(ncname), nf90_nowrite, ncid)
# else
          status=nf90_open(path = TRIM(ncname),                         &
     &                     mode = IOR(nf90_nowrite, nf90_mpiio),        &
     &                     ncid = ncid,                                 &
     &                     comm = OCN_COMM_WORLD,                       &
     &                     info = MP_INFO)
# endif
      END SELECT
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile),   &
     &                                nf90_strerror(status)
        exit_flag=3
        ioerror=status
      END IF
#else
      IF (InpThread) THEN
        SELECT CASE (omode)
          CASE (0)
            my_omode=nf90_nowrite
# if defined DISTRIBUTE        && !defined HDF5       && \
    (defined CONCURRENT_KERNEL || defined DISJOINTED)
            my_omode=IOR(my_omode, nf90_share)
# endif
            status=nf90_open(TRIM(ncname), my_omode, ncid)
          CASE (1)
            my_omode=nf90_write
# if defined DISTRIBUTE        && !defined HDF5       && \
    (defined CONCURRENT_KERNEL || defined DISJOINTED)
            my_omode=IOR(my_omode, nf90_share)
# endif
            status=nf90_open(TRIM(ncname), my_omode, ncid)
          CASE DEFAULT
            my_omode=nf90_nowrite
# if defined DISTRIBUTE        && !defined HDF5       && \
    (defined CONCURRENT_KERNEL || defined DISJOINTED)
            my_omode=IOR(my_omode, nf90_share)
# endif
            status=nf90_open(TRIM(ncname), my_omode, ncid)
        END SELECT
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile),             &
     &                      nf90_strerror(status)
          exit_flag=3
          ioerror=status
        END IF
        IF (Ldebug_ncid) THEN
           WRITE (DBout,'(a,2x," => ",i8,2(2x,a))') kernel(model)//     &
     &                                              ' OPEN', ncid,      &
     &                                              TRIM(ncname),       &
     &                                              TRIM(SourceFile)
          CALL my_flush (DBout)
        END IF
# ifdef DISTRIBUTE
        ibuffer(1)=exit_flag
        ibuffer(2)=ioerror
        ibuffer(3)=ncid
# endif
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      ncid=ibuffer(3)
# endif
#endif
!
  10  FORMAT (/,' NETCDF_OPEN - unable to open existing NetCDF ',       &
#if defined PARALLEL_IO && defined DISTRIBUTE
     &        'file for parallel access:',/,15x,a,/,15x,                &
     &        'call from:',2x,a,/,15x,a)
#else
     &        'file:',/,15x,a,/,15x,'call from:',2x,a,/,15x,a)
#endif
!
      RETURN
      END SUBROUTINE netcdf_open
!
      SUBROUTINE netcdf_redef (ng, model, ncname, ncid)
!
!=======================================================================
!                                                                      !
!  This routine puts an open NetCDF dataset into define mode, so       !
!  dimensions, variables, and attributes can be added or renamed       !
!  an attributes can be deleted.                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       Name of the new NetCDF file (string)                !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!  Note:   There is not need to call this function when processing     !
!          NetCDF-4 format type files.                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
!
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: status

#if !defined PARALLEL_IO && defined DISTRIBUTE
      integer :: ibuffer(2)
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_redef"
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_redef(ncid)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile), &
     &                                  nf90_strerror(status)
          exit_flag=3
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO && defined DISTRIBUTE
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_REDEF - unable to put in definition mode',     &
     &        ' file:',/,16x,a,/,16x,'call from:',2x,a,/,16x,a)
!
      RETURN
      END SUBROUTINE netcdf_redef
!
      SUBROUTINE netcdf_sync (ng, model, ncname, ncid)
!
!=======================================================================
!                                                                      !
!  This routine synchronize to disk requested NetCDF file with         !
!  in-memory buffer to make data available to other  processes         !
!  immediately after it is written.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       Name of the new NetCDF file (string)                !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!  NOTE:                                                               !
!                                                                      !
!  Nowadays, it is recommended to have the writer and readers open     !
!  the NetCDF file with the NF90_SHARE flag to improve performance.    !
!  See NetCDF user manual for more details.                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#if !defined PARALLEL_IO && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
!
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: status

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer :: ibuffer(2)
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", netcdf_sync"
!
!-----------------------------------------------------------------------
!  Synchronize requested NetCDF file.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_sync(ncid)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile),             &
     &                      nf90_strerror(status)
          exit_flag=3
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO && defined DISTRIBUTE
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_SYNC - unable to synchronize to disk file:',   &
     &        /,15x,a,/,15x,'call from:',2x,a,/15x,a)
!
      RETURN
      END SUBROUTINE netcdf_sync
!
      END MODULE mod_netcdf
