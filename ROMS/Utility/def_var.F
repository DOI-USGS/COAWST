#include "cppdefs.h"
      MODULE def_var_mod
!
!svn $Id: def_var.F 806 2016-06-10 01:01:57Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine defines the requested NetCDF variable.                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer).                       !
!     model        Calling model identifier (integer).                 !
!     ncid         NetCDF file ID (integer).                           !
!     Vtype        NetCDF variable type (integer).                     !
!     nVdim        Number of variable dimensions (integer; 0=scalar).  !
!     Vdim         Dimensions IDs for this variable (integer vector).  !
!     Aval         Attribute values (real vector):                     !
!                    Aval(1)   =>  Add offset value                    !
!                    Aval(2)   =>  Valid minimum value                 !
!                    Aval(3)   =>  Valid maximum value                 !
!                    Aval(4)   =>  Missing value                       !
!                    Aval(5)   =>  C-grid variable type                !
!                    Aval(6)   =>  Fill value                          !
!     Vinfo        Variable information (character array):             !
!                    Vinfo( 1) =>  Variable name                       !
!                    Vinfo( 2) =>  Variable "longname" attribute       !
!                    Vinfo( 3) =>  Variable "units" attribute          !
!                    Vinfo( 4) =>  Variable "calendar" attribute       !
!                    Vinfo( 5) =>  Variable "valid_min" attribute      !
!                    Vinfo( 6) =>  Variable "valid_max" attribute      !
!                    Vinfo( 7) =>  Variable "option_T" attribute       !
!                    Vinfo( 8) =>  Variable "option_F" attribute       !
!                    Vinfo( 9) =>  Variable "option_0" attribute       !
!                    Vinfo(10) =>  Variable "option_1" attribute       !
!                    Vinfo(11) =>  Variable "negative_value" attribute !
!                    Vinfo(12) =>  Variable "positive_value" attribute !
!                    Vinfo(13) =>  Variable "cycle" attribute          !
!                    Vinfo(14) =>  Variable "field" attribute          !
!                    Vinfo(15) =>  Variable "positions" attribute      !
!                    Vinfo(16) =>  Variable "time" attribute           !
!                    Vinfo(17) =>  Variable "missing_value" attribute  !
!                    Vinfo(18) =>  Variable "add_offset" attribute     !
!                    Vinfo(19) =>  Variable "size_class" attribute     !
!                    Vinfo(20) =>  Variable "water_points" attribute   !
!                    Vinfo(21) =>  Variable "standard_name" attribute  !
!                    Vinfo(22) =>  Variable "coordinates" attribute    !
!                    Vinfo(23) =>  Variable "formula_terms" attribute  !
!                    Vinfo(24) =>  Variable "_FillValue" attribute     !
!                    Vinfo(25) =>  Variable "positive" attribute       !
!     ncname       NetCDF file name.                                   !
!     SetFillVal   Logical switch to set fill value in land areas      !
!                    (optional).                                       !
!     SetParAccess Logical switch to set parallel I/O access flag to   !
!                    either collective or independent (optional).      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     def_var      Error flag (integer).                               !
!     Vid          NetCDF variable ID (integer).                       !
!                                                                      !
!  Notice that arrays "Aval" and "Vinfo" is destroyed on output to     !
!  facilitate the definition of the next variable.                     !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      FUNCTION def_var (ng, model, ncid, Vid, Vtype, nVdim, Vdim,       &
     &                  Aval, Vinfo, ncname, SetFillVal, SetParAccess)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

#if !defined PARALLEL_OUT  && defined DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal
      logical, intent(in), optional :: SetParAccess

      integer, intent (in) :: ng, model, ncid, Vtype, nVdim

      integer, dimension(:), intent (in) :: Vdim

      integer, intent (out) :: Vid

      real(r8), dimension(:), intent(inout) :: Aval

      character (len=*), intent(in) :: ncname
      character (len=*), intent(inout) :: Vinfo(25)
!
!  Local variable declarations.
!
#ifdef MASKING
      logical :: LandFill
#endif
#if defined PARALLEL_OUT && defined DISTRIBUTE
      logical :: Ltiled
#endif
      integer :: i, j, latt, status

      integer :: def_var

      character (len=  5) location
      character (len=160) text
!
!-----------------------------------------------------------------------
!  Define requested variable and its attributes.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  Define variable.
!
        IF (exit_flag.eq.NoError) THEN
          IF (LEN_TRIM(Vinfo(1)).gt.0) THEN
            IF ((nVdim.eq.1).and.(Vdim(1).eq.0)) THEN
              status=nf90_def_var(ncid, TRIM(Vinfo(1)), Vtype,          &
     &                            varid = Vid)
            ELSE
              status=nf90_def_var(ncid, TRIM(Vinfo(1)), Vtype,          &
     &                            Vdim(1:nVdim), Vid)
            END IF
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF

#if !defined PARALLEL_OUT && (defined HDF5 && defined DEFLATE)
!
!  Define deflate (file compresion) parameters. Notice that deflation
!  cannot be used in parallel I/O for writing data. This is because
!  the compression makes it impossible for the HDF5 library to exactly
!  map the data to the disk location. However, deflated data can be
!  read with parallel I/O
!
        IF (exit_flag.eq.NoError) THEN
          IF (LEN_TRIM(Vinfo(1)).gt.0) THEN
            IF ((nVdim.gt.1).and.(Vdim(1).ne.0)) THEN
              status=nf90_def_var_deflate (ncid, Vid, shuffle, deflate, &
     &                                     deflate_level)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,20) TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF
        END IF
#endif
!
!  Define special attributes for SGRID conventions variable "grid".
!
        IF (exit_flag.eq.NoError) THEN
          IF (TRIM(Vinfo(1)).eq.'grid') THEN
            status=nf90_put_att(ncid, Vid, 'cf_role',                   &
     &                          'grid_topology')
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'cf_role',                  &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            status=nf90_put_att(ncid, Vid, 'topology_dimension',        &
     &                          (/2/))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'topology_dimension',       &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            status=nf90_put_att(ncid, Vid, 'node_dimensions',           &
     &                          'xi_psi eta_psi')
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'node_dimensions',          &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            text='xi_rho: xi_psi (padding: both) '//                    &
     &           'eta_rho: eta_psi (padding: both)'
            status=nf90_put_att(ncid, Vid, 'face_dimensions',           &
     &                          TRIM(text))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'face_dimensions',          &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            text='xi_u: xi_psi eta_u: eta_psi (padding: both)'
            status=nf90_put_att(ncid, Vid, 'edge1_dimensions',          &
     &                          TRIM(text))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'edge1_dimensions',         &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            text='xi_v: xi_psi (padding: both) eta_v: eta_psi'
            status=nf90_put_att(ncid, Vid, 'edge2_dimensions',          &
     &                          TRIM(text))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'edge2_dimensions',         &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            IF (spherical) THEN
              status=nf90_put_att(ncid, Vid, 'node_coordinates',        &
     &                            'lon_psi lat_psi')
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,30) 'node_coordinates',       &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF

              status=nf90_put_att(ncid, Vid, 'face_coordinates',        &
     &                            'lon_rho lat_rho')
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,30) 'face_coordinates',       &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF

              status=nf90_put_att(ncid, Vid, 'edge1_coordinates',       &
     &                            'lon_u lat_u')
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,30) 'edge1_coordinates',      &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF

              status=nf90_put_att(ncid, Vid, 'edge2_coordinates',       &
     &                            'lon_v lat_v')
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,30) 'edge2_coordinates',      &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              status=nf90_put_att(ncid, Vid, 'node_coordinates',        &
     &                            'x_psi y_psi')
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,30) 'node_coordinates',       &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF

              status=nf90_put_att(ncid, Vid, 'face_coordinates',        &
     &                            'x_rho y_rho')
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,30) 'face_coordinates',       &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF

              status=nf90_put_att(ncid, Vid, 'edge1_coordinates',       &
     &                            'x_u y_u')
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,30) 'edge1_coordinates',      &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF

              status=nf90_put_att(ncid, Vid, 'edge2_coordinates',       &
     &                            'x_v y_v')
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,30) 'edge2_coordinates',      &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
#ifdef SOLVE3D
!
            status=nf90_put_att(ncid, Vid, 'vertical_dimensions',       &
     &                          's_rho: s_w (padding: none)')
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'vertical_dimensions',      &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
#endif
          END IF
        END IF

!
!  Define "longname" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(2))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'long_name',                 &
     &                          Vinfo(2)(1:latt))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'long_name',                &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  Define "size_class" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(19))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'size_class',                &
     &                          Vinfo(19)(1:latt))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'size_class',               &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "units" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(3))
          IF (latt.gt.0) THEN
            IF (TRIM(Vinfo(3)).ne.'nondimensional') THEN
              status=nf90_put_att(ncid, Vid, 'units',                   &
     &                            Vinfo(3)(1:latt))
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,30) 'units',                  &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF
        END IF
!
!  If applicable, define "calendar" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(4))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'calendar',                  &
     &                          Vinfo(4)(1:latt))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'calendar',                 &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "valid_min" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(5))
          IF (latt.gt.0) THEN
            IF (Vtype.eq.nf90_int) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(5)),            &
     &                            INT(Aval(2)))
#ifndef NO_4BYTE_REALS
            ELSE IF (Vtype.eq.nf90_float) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(5)),            &
     &                            REAL(Aval(2),r4))
#endif
            ELSE
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(5)),            &
     &                            Aval(2))
            END IF
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) TRIM(Vinfo(5)),             &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
            Aval(2)=0.0_r8
          END IF
        END IF
!
!  If applicable, define "valid_max" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(6))
          IF (latt.gt.0) THEN
            IF (Vtype.eq.nf90_int) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(6)),            &
     &                            INT(Aval(3)))
#ifndef NO_4BYTE_REALS
            ELSE IF (Vtype.eq.nf90_float) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(6)),            &
     &                            REAL(Aval(3),r4))
#endif
            ELSE
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(6)),            &
     &                            Aval(3))
            END IF
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) TRIM(Vinfo(6)),             &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
            Aval(3)=0.0_r8
          END IF
        END IF
!
!  If applicable, define "flag_values" and "flag_meanings" attributes
!  for logical variables.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((LEN_TRIM(Vinfo(7)).gt.0).and.                            &
     &        (LEN_TRIM(Vinfo(8)).gt.0)) THEN
            text='T, F'
            latt=LEN_TRIM(text)
            status=nf90_put_att(ncid, Vid, 'flag_values',               &
     &                          text(1:latt))
            IF (status.eq.nf90_noerr) THEN
              text=TRIM(Vinfo(7))//' '//TRIM(Vinfo(8))
              latt=LEN_TRIM(text)
              status=nf90_put_att(ncid, Vid, 'flag_meanings',           &
     &                            text(1:latt))
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,30) 'flag_meanings (T/F)',    &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,30) 'flag_values (T/F)',        &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "flag_values" and "flag_meanings" attributes
!  for integer and floating point variables.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((LEN_TRIM(Vinfo( 9)).gt.0).and.                           &
     &        (LEN_TRIM(Vinfo(10)).gt.0)) THEN
            IF (Vtype.eq.nf90_int) THEN
              status=nf90_put_att(ncid, Vid, 'flag_values',             &
     &                            (/0,1/))
#ifndef NO_4BYTE_REALS
            ELSE IF (Vtype.eq.nf90_float) THEN
              status=nf90_put_att(ncid, Vid, 'flag_values',             &
     &                            (/0.0_r4, 1.0_r4/))
#endif
            ELSE
              status=nf90_put_att(ncid, Vid, 'flag_values',             &
     &                            (/0.0_r8, 1.0_r8/))
            END IF
            IF (status.eq.nf90_noerr) THEN
              text=TRIM(Vinfo(9))//' '//TRIM(Vinfo(10))
              latt=LEN_TRIM(text)
              status=nf90_put_att(ncid, Vid, 'flag_meanings',           &
     &                            text(1:latt))
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,30) 'flag_meanings',          &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,30) 'flag_values',              &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "negative_value" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(11))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'negative_value',            &
     &                          Vinfo(11)(1:latt))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'negative_value',           &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "positive_value" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(12))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'positive_value',            &
     &                          Vinfo(12)(1:latt))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'positive_value',           &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define CF "positive" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(25))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'positive',                  &
     &                          Vinfo(25)(1:latt))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'positive',                 &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "cycle" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(13))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'cycle',                     &
     &                          Vinfo(13)(1:latt))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'cycle',                    &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "positions" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(15))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'positions',                 &
     &                          Vinfo(15)(1:latt))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'positions',                &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "time" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(16))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'time',                      &
     &                          Vinfo(16)(1:latt))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'time',                     &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "missing_value" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(17))
          IF (latt.gt.0) THEN
            IF (Vtype.eq.nf90_int) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(17)),           &
     &                            INT(Aval(4)))
#ifndef NO_4BYTE_REALS
            ELSE IF (Vtype.eq.nf90_float) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(17)),           &
     &                            REAL(Aval(4),r4))
#endif
            ELSE
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(17)),           &
     &                            Aval(4))
            END IF
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) TRIM(Vinfo(17)),            &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
            Aval(4)=0.0_r8
          END IF
        END IF
!
!  If applicable, define "add_offset" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(18))
          IF (latt.gt.0) THEN
            IF (Vtype.eq.nf90_int) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(18)),           &
     &                            INT(Aval(1)))
#ifndef NO_4BYTE_REALS
            ELSE IF (Vtype.eq.nf90_float) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(18)),           &
     &                            REAL(Aval(1),r4))
#endif
            ELSE
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(18)),           &
     &                            Aval(1))
            END IF
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) TRIM(Vinfo(18)),            &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
            Aval(1)=0.0_r8
          END IF
        END IF
!
!  If applicable, define "_FillValue" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(24))
          IF (latt.gt.0) THEN
            IF (Vtype.eq.nf90_int) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(24)),           &
     &                            INT(Aval(6)))
#ifndef NO_4BYTE_REALS
            ELSE IF (Vtype.eq.nf90_float) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(24)),           &
     &                            REAL(Aval(6),r4))
#endif
            ELSE
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(24)),           &
     &                            Aval(6))
            END IF
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) TRIM(Vinfo(24)),            &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
            Aval(6)=0.0_r8
          END IF
        END IF
!
!  If applicable, define "water_points" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(20))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'water_points',              &
     &                          Vinfo(20)(1:latt))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'water_points',             &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "standard_name" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(21))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'standard_name',             &
     &                          Vinfo(21)(1:latt))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'standard_name',            &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "coordinates" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(22))
          IF (latt.gt.0) THEN
            IF (spherical) THEN
              IF (INT(Aval(5)).eq.r2dvar) THEN
                text='lon_rho lat_rho'
                location='face'
              ELSE IF (INT(Aval(5)).eq.r3dvar) THEN
                text='lon_rho lat_rho s_rho'
                location='face'
              ELSE IF (INT(Aval(5)).eq.w3dvar) THEN
                text='lon_rho lat_rho s_w'
                location='face'
              ELSE IF (INT(Aval(5)).eq.b3dvar) THEN
                text='lon_rho lat_rho'
                location='face'
#ifdef SPECTRAL_LIGHT
              ELSE IF (INT(Aval(5)).eq.s3dvar) THEN
                text='lon_rho lat_rho'
#endif
              ELSE IF (INT(Aval(5)).eq.u2dvar) THEN
                text='lon_u lat_u'
                location='edge1'
              ELSE IF (INT(Aval(5)).eq.u3dvar) THEN
                text='lon_u lat_u s_rho'
                location='edge1'
              ELSE IF (INT(Aval(5)).eq.-u3dvar) THEN
                text='lon_u lat_u s_w'
                location='edge1'
              ELSE IF (INT(Aval(5)).eq.v2dvar) THEN
                text='lon_v lat_v'
                location='edge2'
              ELSE IF (INT(Aval(5)).eq.v3dvar) THEN
                text='lon_v lat_v s_rho'
                location='edge2'
              ELSE IF (INT(Aval(5)).eq.-v3dvar) THEN
                text='lon_v lat_v s_w'
                location='edge2'
              ELSE IF (INT(Aval(5)).eq.p2dvar) THEN
                text='lon_psi lat_psi'
                location='node'
              ELSE IF (INT(Aval(5)).eq.p3dvar) THEN
                text='lon_psi lat_psi s_rho'
                location='node'
              END IF
            ELSE
              IF (INT(Aval(5)).eq.r2dvar) THEN
                text='x_rho y_rho'
                location='face'
              ELSE IF (INT(Aval(5)).eq.r3dvar) THEN
                text='x_rho y_rho s_rho'
                location='face'
              ELSE IF (INT(Aval(5)).eq.w3dvar) THEN
                text='x_rho y_rho s_w'
                location='face'
              ELSE IF (INT(Aval(5)).eq.b3dvar) THEN
                text='x_rho y_rho'
                location='face'
#ifdef SPECTRAL_LIGHT
              ELSE IF (INT(Aval(5)).eq.s3dvar) THEN
                text='x_rho y_rho'
#endif
              ELSE IF (INT(Aval(5)).eq.u2dvar) THEN
                text='x_u y_u'
                location='edge1'
              ELSE IF (INT(Aval(5)).eq.u3dvar) THEN
                text='x_u y_u s_rho'
                location='edge1'
              ELSE IF (INT(Aval(5)).eq.-u3dvar) THEN
                text='x_u y_u s_w'
                location='edge1'
              ELSE IF (INT(Aval(5)).eq.v2dvar) THEN
                text='x_v y_v'
                location='edge2'
              ELSE IF (INT(Aval(5)).eq.v3dvar) THEN
                text='x_v y_v s_rho'
                location='edge2'
              ELSE IF (INT(Aval(5)).eq.-v3dvar) THEN
                text='x_v y_v s_w'
                location='edge2'
              ELSE IF (INT(Aval(5)).eq.p2dvar) THEN
                text='x_psi y_psi'
                location='node'
              ELSE IF (INT(Aval(5)).eq.p3dvar) THEN
                text='x_psi y_psi s_rho'
                location='node'
              END IF
            END IF

            status=nf90_put_att(ncid, Vid, 'grid',                      &
     &                          'grid')
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'grid',                     &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            status=nf90_put_att(ncid, Vid, 'location',                  &
     &                          TRIM(location))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'location',                 &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            latt=LEN_TRIM(text)
            IF (((ABS(INT(Aval(5))).le.4).and.(nVdim.gt.2)).or.         &
     &          ((ABS(INT(Aval(5))).gt.4).and.(nVdim.gt.3))) THEN
              text=text(1:latt)//' ocean_time'
              latt=LEN_TRIM(text)
            END IF
            status=nf90_put_att(ncid, Vid, TRIM(Vinfo(22)),             &
     &                          text(1:latt))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'coordinates',              &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
            Aval(5)=0.0_r8
          END IF
        END IF
!
!  If applicable, define "formula_terms" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(23))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'formula_terms',             &
     &                          Vinfo(23)(1:latt))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'formula_terms',            &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "field" attribute (always last).
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(14))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'field',                     &
     &                          Vinfo(14)(1:latt))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) 'field',                    &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
      END IF

#if defined MASKING && !defined WRITE_WATER
!
!  If land/sea masking, define "_FillValue" attribute since masking
!  areas are overwrite with special value (spval) during output.
!  Notice that the coordinate attribute is used to check which
!  variables need the "_FillValue" attribute.
!

        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(22))
          IF (PRESENT(SetFillVal)) THEN
            LandFill=SetFillVal
          ELSE
            LandFill=(latt.gt.0).and.(nVdim.gt.2)
          END IF
          IF (LandFill) THEN
            IF (Vtype.eq.nf90_double) THEN
              status=nf90_put_att(ncid, Vid, '_FillValue',              &
     &                            spval)
# ifndef NO_4BYTE_REALS
            ELSE IF (Vtype.eq.nf90_float) THEN
              status=nf90_put_att(ncid, Vid, '_FillValue',              &
     &                            REAL(spval,r4))
# endif
            END IF
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,30) '_FillValue',               &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
#endif
#if defined PARALLEL_OUT && defined DISTRIBUTE
!
!  Set parallel access for tiled and non-tiled variables. The optional
!  argument "SetParAccess" is only available for non-tiled variables.
!  Otherwise, it is assumed that defined variables are tiled. That is,
!  they are variables (2D/3D) with parallel domain decomposition.
!
        IF (exit_flag.eq.NoError) THEN
          IF (PRESENT(SetParAccess)) THEN
            Ltiled=SetParAccess
          ELSE
            Ltiled=.TRUE.
          END IF
          IF (Ltiled) THEN
            status=nf90_var_par_access(ncid, Vid, IO_tiled_access)
          ELSE
            status=nf90_var_par_access(ncid, Vid, IO_nontiled_access)
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vinfo(1)),               &
     &                                    TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif
!
!  Clean information variables.
!
      DO i=1,SIZE(Vinfo)
        DO j=1,LEN(Vinfo(1))
          Vinfo(i)(j:j)=' '
        END DO
      END DO
!
!  Set error flag
!
#if !defined PARALLEL_OUT  && defined DISTRIBUTE
      CALL mp_bcasti (ng, model, exit_flag)
#endif
      def_var=exit_flag
!
 10   FORMAT (/,' DEF_VAR - unable to define variable: ',a,/,           &
     &        11x,'in NetCDF file: ',a)
 20   FORMAT (/,' DEF_VAR - error while setting deflate parameters',    &
     &        ' for variable: ',a,/,11x,'in NetCDF file: ',a)
 30   FORMAT (/,'DEF_VAR - error while defining attribute: ',a,         &
     &        ' for variable: ',a,/,11x,'in NetCDF file: ',a)
#if defined PARALLEL_OUT && defined DISTRIBUTE
 40   FORMAT (/,'DEF_VAR - error while setting parallel access flag',   &
     &        ' for variable: ',a,/,11x,'in NetCDF file: ',a)
#endif

      RETURN
      END FUNCTION def_var

      END MODULE def_var_mod
