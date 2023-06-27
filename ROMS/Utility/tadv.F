#include "cppdefs.h"
      MODULE tadv_mod
#ifdef SOLVE3D
!
!git $Id$
!svn $Id: tadv.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines are used to process tracer advection switches        !
!  structure:                                                          !
!                                                                      !
!     tadv_putatt   Writes activated keyword strings into specified    !
!                     output NetCDF file global attribute.             !
!
!     tadv_report   Reports to standard output activated keyword       !
!                     strings.                                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      INTERFACE tadv_putatt
        MODULE PROCEDURE tadv_putatt_nf90
#if defined PIO_LIB && defined DISTRIBUTE
        MODULE PROCEDURE tadv_putatt_pio
#endif
      END INTERFACE tadv_putatt
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tadv_putatt_nf90 (ng, ncid, ncname, aname,             &
     &                             Hadv, Vadv, status)
!***********************************************************************
!                                                                      !
!  This routine writes tracer advection scheme keywords strings into   !
!  specified NetCDF file global attribute when using the standard      !
!  NetCDF-3 or NetCDF-4 library.                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng            Nested grid number (integer)                       !
!     model         Calling model identifier (integer)                 !
!     ncid          NetCDF file ID (integer)                           !
!     ncname        NetCDF filename (character)                        !
!     aname         NetCDF global attribute name (character)           !
!     Hadv          Horizontal advection type structure, TYPE(T_ADV)   !
!     Vadv          Vertical   advection type structure, TYPE(T_ADV)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag     Error flag (integer) stored in MOD_SCALARS         !
!     ioerror       NetCDF return code (integer) stored in MOD_IOUNITS !
!     status        NetCDF return code (integer)                       !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE strings_mod, ONLY : FoundError
!
      implicit none
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, ncid
      integer, intent(out) :: status
!
      character (*), intent(in) :: ncname
      character (*), intent(in) :: aname
!
      TYPE(T_ADV), intent(in) :: Hadv(MAXVAL(NT),Ngrids)
      TYPE(T_ADV), intent(in) :: Vadv(MAXVAL(NT),Ngrids)
!
! Local variable declarations
!
      integer :: i, ie, is, itrc, lvar, lstr, nTvar
!
      character (len=   1) :: newline
      character (len=  13) :: Hstring, Vstring
      character (len=  17) :: frmt
      character (len=  70) :: line
      character (len=2816) :: tadv_att

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", tadv_putatt_nf90"
!
!-----------------------------------------------------------------------
!  Write lateral boundary conditions global attribute.
!-----------------------------------------------------------------------
!
!  Determine maximum length of state variable length.
!
      nTvar=MAXVAL(NT)
      lvar=0
      DO itrc=1,nTvar
        lvar=MAX(lvar, LEN_TRIM(Vname(1,idTvar(itrc))))
      END DO
      WRITE (frmt,10) MAX(10,lvar)+4
  10  FORMAT ("(a,':',t",i2.2,',a,a,a)')
!
!  Initialize attribute.
!
      newline=CHAR(10)            ! Line Feed (LF) character for
      lstr=LEN_TRIM(newline)      ! attribute clarity with "ncdump"
      DO i=1,LEN(tadv_att)
        tadv_att(i:i)=' '
      END DO
      tadv_att(1:lstr)=newline(1:lstr)
      is=lstr+1
      WRITE (line,frmt) 'ADVECTION',                                    &
     &                  'HORIZONTAL   ',                                &
     &                  'VERTICAL     ',                                &
     &                   newline(1:lstr)
      lstr=LEN_TRIM(line)
      ie=is+lstr
      tadv_att(is:ie)=line(1:lstr)
      is=ie
!
!  Check if the local string "tadv_att" is big enough to store the
!  tracer advection scheme global attribute.
!
      lstr=(nTvar+1)*(26+lvar+4)+1
      IF (LEN(tadv_att).lt.lstr) THEN
        IF (Master) THEN
          WRITE (stdout,20) LEN(tadv_att), lstr
  20      FORMAT (/,' TADV_PUTATT_NF90 - Length of local string ',      &
     &            ' tadv_att too small',/,20x,'Current = ',i5,          &
     &            '  Needed = ',i5)
        END IF
        exit_flag=5
        RETURN
      END IF
!
!  Build attribute string.
!
      DO itrc=1,nTvar
        IF (Hadv(itrc,ng)%AKIMA4) THEN
          Hstring='Akima4'
        ELSE IF (Hadv(itrc,ng)%CENTERED2) THEN
          Hstring='Centered2'
        ELSE IF (Hadv(itrc,ng)%CENTERED4) THEN
          Hstring='Centered4'
        ELSE IF (Hadv(itrc,ng)%HSIMT) THEN
          Hstring='HSIMT'
        ELSE IF (Hadv(itrc,ng)%MPDATA) THEN
          Hstring='MPDATA'
        ELSE IF (Hadv(itrc,ng)%SPLINES) THEN
          Hstring='Splines'
        ELSE IF (Hadv(itrc,ng)%SPLIT_U3) THEN
          Hstring='Split_U3'
        ELSE IF (Hadv(itrc,ng)%UPSTREAM3) THEN
          Hstring='Upstream3'
        END IF
!
        IF (Vadv(itrc,ng)%AKIMA4) THEN
          Vstring='Akima4'
        ELSE IF (Vadv(itrc,ng)%CENTERED2) THEN
          Vstring='Centered2'
        ELSE IF (Vadv(itrc,ng)%CENTERED4) THEN
          Vstring='Centered4'
        ELSE IF (Vadv(itrc,ng)%HSIMT) THEN
          Vstring='HSIMT'
        ELSE IF (Vadv(itrc,ng)%MPDATA) THEN
          Vstring='MPDATA'
        ELSE IF (Vadv(itrc,ng)%SPLINES) THEN
          Vstring='Splines'
        ELSE IF (Vadv(itrc,ng)%SPLIT_U3) THEN
          Vstring='Split_U3'
        ELSE IF (Vadv(itrc,ng)%UPSTREAM3) THEN
          Vstring='Upstream3'
        END IF
        IF (itrc.eq.nTvar) newline=' '
        WRITE (line,frmt) TRIM(Vname(1,idTvar(itrc))),                  &
     &                    Hstring, Vstring,                             &
     &                    newline
        lstr=LEN_TRIM(line)
        ie=is+lstr
        tadv_att(is:ie)=line(1:lstr)
        is=ie
      END DO
!
!  Write attribute to NetCDF file.
!
      status=nf90_put_att(ncid, nf90_global, TRIM(aname),               &
     &                    TRIM(tadv_att))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      RETURN
      END SUBROUTINE tadv_putatt_nf90

#if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE tadv_putatt_pio (ng, pioFile, ncname, aname,           &
     &                            Hadv, Vadv, status)
!***********************************************************************
!                                                                      !
!  This routine writes tracer advection scheme keywords strings into   !
!  specified NetCDF file global attribute when using the NCAR          !
!  Parallel-IO library.                                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng            Nested grid number (integer)                       !
!     model         Calling model identifier (integer)                 !
!     pioFile       PIO file descriptor structure, TYPE(File_desc_t)   !
!                     pioFile%fh         file handler                  !
!                     pioFile%iosystem   IO system descriptor (struct) !
!     ncname        PIO filename (character)                           !
!     aname         PIO global attribute name (character)              !
!     Hadv          Horizontal advection type structure, TYPE(T_ADV)   !
!     Vadv          Vertical   advection type structure, TYPE(T_ADV)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag     Error flag (integer) stored in MOD_SCALARS         !
!     ioerror       NetCDF return code (integer) stored in MOD_IOUNITS !
!     status        NetCDF return code (integer)                       !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE pio
!
      USE strings_mod, ONLY : FoundError
!
      implicit none
!
! Imported variable declarations.
!
      integer, intent(in) :: ng
      integer, intent(out) :: status
!
      character (*), intent(in) :: ncname
      character (*), intent(in) :: aname
!
      TYPE (File_desc_t), intent(in) :: pioFile
      TYPE(T_ADV), intent(in) :: Hadv(MAXVAL(NT),Ngrids)
      TYPE(T_ADV), intent(in) :: Vadv(MAXVAL(NT),Ngrids)
!
! Local variable declarations
!
      integer :: i, ie, is, itrc, lvar, lstr, nTvar
!
      character (len=   1) :: newline
      character (len=  13) :: Hstring, Vstring
      character (len=  17) :: frmt
      character (len=  70) :: line
      character (len=2816) :: tadv_att

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", tadv_putatt_pio"
!
!-----------------------------------------------------------------------
!  Write lateral boundary conditions global attribute.
!-----------------------------------------------------------------------
!
!  Determine maximum length of state variable length.
!
      nTvar=MAXVAL(NT)
      lvar=0
      DO itrc=1,nTvar
        lvar=MAX(lvar, LEN_TRIM(Vname(1,idTvar(itrc))))
      END DO
      WRITE (frmt,10) MAX(10,lvar)+4
  10  FORMAT ("(a,':',t",i2.2,',a,a,a)')
!
!  Initialize attribute.
!
      newline=CHAR(10)            ! Line Feed (LF) character for
      lstr=LEN_TRIM(newline)      ! attribute clarity with "ncdump"
      DO i=1,LEN(tadv_att)
        tadv_att(i:i)=' '
      END DO
      tadv_att(1:lstr)=newline(1:lstr)
      is=lstr+1
      WRITE (line,frmt) 'ADVECTION',                                    &
     &                  'HORIZONTAL   ',                                &
     &                  'VERTICAL     ',                                &
     &                   newline(1:lstr)
      lstr=LEN_TRIM(line)
      ie=is+lstr
      tadv_att(is:ie)=line(1:lstr)
      is=ie
!
!  Check if the local string "tadv_att" is big enough to store the
!  tracer advection scheme global attribute.
!
      lstr=(nTvar+1)*(26+lvar+4)+1
      IF (LEN(tadv_att).lt.lstr) THEN
        IF (Master) THEN
          WRITE (stdout,20) LEN(tadv_att), lstr
  20      FORMAT (/,' TADV_PUTATT_PIO - Length of local string ',       &
     &            'tadv_att too small',/,19x,'Current = ',i5,           &
     &            '  Needed = ',i5)
        END IF
        exit_flag=5
        RETURN
      END IF
!
!  Build attribute string.
!
      DO itrc=1,nTvar
        IF (Hadv(itrc,ng)%AKIMA4) THEN
          Hstring='Akima4'
        ELSE IF (Hadv(itrc,ng)%CENTERED2) THEN
          Hstring='Centered2'
        ELSE IF (Hadv(itrc,ng)%CENTERED4) THEN
          Hstring='Centered4'
        ELSE IF (Hadv(itrc,ng)%HSIMT) THEN
          Hstring='HSIMT'
        ELSE IF (Hadv(itrc,ng)%MPDATA) THEN
          Hstring='MPDATA'
        ELSE IF (Hadv(itrc,ng)%SPLINES) THEN
          Hstring='Splines'
        ELSE IF (Hadv(itrc,ng)%SPLIT_U3) THEN
          Hstring='Split_U3'
        ELSE IF (Hadv(itrc,ng)%UPSTREAM3) THEN
          Hstring='Upstream3'
        END IF
!
        IF (Vadv(itrc,ng)%AKIMA4) THEN
          Vstring='Akima4'
        ELSE IF (Vadv(itrc,ng)%CENTERED2) THEN
          Vstring='Centered2'
        ELSE IF (Vadv(itrc,ng)%CENTERED4) THEN
          Vstring='Centered4'
        ELSE IF (Vadv(itrc,ng)%HSIMT) THEN
          Vstring='HSIMT'
        ELSE IF (Vadv(itrc,ng)%MPDATA) THEN
          Vstring='MPDATA'
        ELSE IF (Vadv(itrc,ng)%SPLINES) THEN
          Vstring='Splines'
        ELSE IF (Vadv(itrc,ng)%SPLIT_U3) THEN
          Vstring='Split_U3'
        ELSE IF (Vadv(itrc,ng)%UPSTREAM3) THEN
          Vstring='Upstream3'
        END IF
        IF (itrc.eq.nTvar) newline=' '
        WRITE (line,frmt) TRIM(Vname(1,idTvar(itrc))),                  &
     &                    Hstring, Vstring,                             &
     &                    newline
        lstr=LEN_TRIM(line)
        ie=is+lstr
        tadv_att(is:ie)=line(1:lstr)
        is=ie
      END DO
!
!  Write attribute to NetCDF file.
!
      status=PIO_put_att(pioFile, PIO_global, TRIM(aname),              &
     &                   TRIM(tadv_att))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      RETURN
      END SUBROUTINE tadv_putatt_pio
#endif

!
!***********************************************************************
      SUBROUTINE tadv_report (iunit, model, Hadv, Vadv, Lwrite)
!***********************************************************************
!                                                                      !
!  This routine reports horizontal and vertical advection scheme for   !
!  each tracer variable.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     iunit         Output logical unit (integer)                      !
!     model         Calling model identifier (integer)                 !
!     Hadv          Horizontal advection type structure, TYPE(T_ADV)   !
!     Vadv          Vertical   advection type structure, TYPE(T_ADV)   !
!     Lwrite        Switch to report information to standard output    !
!                     unit or file (logical)                           !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
! Imported variable declarations.
!
      logical, intent(in) :: Lwrite
!
      integer, intent(in) :: iunit, model
!
      TYPE(T_ADV), intent(in) :: Hadv(MAXVAL(NT),Ngrids)
      TYPE(T_ADV), intent(in) :: Vadv(MAXVAL(NT),Ngrids)
!
! Local variable declarations
!
      integer :: i, itrc, ng
!
      character (len=11) :: Hstring(MAXVAL(NT),Ngrids)
      character (len=11) :: Vstring(MAXVAL(NT),Ngrids)
!
!-----------------------------------------------------------------------
!  Tracer horizontal and vertical switches.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO itrc=1,NT(ng)
          IF (Hadv(itrc,ng)%AKIMA4) THEN
            Hstring(itrc,ng)='Akima4'
          ELSE IF (Hadv(itrc,ng)%CENTERED2) THEN
            Hstring(itrc,ng)='Centered2'
          ELSE IF (Hadv(itrc,ng)%CENTERED4) THEN
            Hstring(itrc,ng)='Centered4'
          ELSE IF (Hadv(itrc,ng)%HSIMT) THEN
            Hstring(itrc,ng)='HSIMT'
          ELSE IF (Hadv(itrc,ng)%MPDATA) THEN
            Hstring(itrc,ng)='MPDATA'
          ELSE IF (Hadv(itrc,ng)%SPLINES) THEN
            Hstring(itrc,ng)='Splines'
          ELSE IF (Hadv(itrc,ng)%SPLIT_U3) THEN
            Hstring(itrc,ng)='Split_U3'
          ELSE IF (Hadv(itrc,ng)%UPSTREAM3) THEN
            Hstring(itrc,ng)='Upstream3'
          END IF
!
          IF (Vadv(itrc,ng)%AKIMA4) THEN
            Vstring(itrc,ng)='Akima4'
          ELSE IF (Vadv(itrc,ng)%CENTERED2) THEN
            Vstring(itrc,ng)='Centered2'
          ELSE IF (Vadv(itrc,ng)%CENTERED4) THEN
            Vstring(itrc,ng)='Centered4'
          ELSE IF (Vadv(itrc,ng)%HSIMT) THEN
            Vstring(itrc,ng)='HSIMT'
          ELSE IF (Vadv(itrc,ng)%MPDATA) THEN
            Vstring(itrc,ng)='MPDATA'
          ELSE IF (Vadv(itrc,ng)%SPLINES) THEN
            Vstring(itrc,ng)='Splines'
          ELSE IF (Vadv(itrc,ng)%SPLIT_U3) THEN
            Vstring(itrc,ng)='Split_U3'
          ELSE IF (Vadv(itrc,ng)%UPSTREAM3) THEN
            Vstring(itrc,ng)='Upstream3'
          END IF
        END DO
      END DO
!
!  Report.
!
      IF (Master.and.Lwrite) THEN
        DO itrc=1,MAXVAL(NT)
          DO ng=1,Ngrids
            IF (ng.eq.1) THEN
              WRITE (iunit,10) TRIM(Vname(1,idTvar(itrc))), ng,         &
     &                         TRIM(Hstring(itrc,ng)),                  &
     &                         TRIM(Vstring(itrc,ng))
            ELSE
              WRITE (iunit,20) ng,                                      &
     &                         TRIM(Hstring(itrc,ng)),                  &
     &                         TRIM(Vstring(itrc,ng))
            END IF
          END DO
        END DO
!
        IF (model.eq.iNLM) THEN
          WRITE (iunit,'(1x)')
          WRITE (iunit,30) 'Akima4',                                    &
     &          'Fourth-order Akima advection'
          WRITE (iunit,30) 'Centered2',                                 &
     &          'Second-order centered differences advection'
          WRITE (iunit,30) 'Centered4',                                 &
     &          'Fourth-order centered differences advection'
          WRITE (iunit,30) 'HSIMT',                                     &
     &          'Third High-order Spatial Inteporlation at Middle '//   &
     &          'Time Advection with TVD limiter'
          WRITE (iunit,30) 'MPDATA',                                    &
     &          'Multidimensional Positive Definite Advection '//       &
     &          'Algorithm, recursive method'
          WRITE (iunit,30) 'Splines',                                   &
     &          'Conservative Parabolic Splines Reconstruction '//      &
     &          'Advection (only vertical; not recommended)'
          WRITE (iunit,30) 'Split_U3',                                  &
     &          'Split third-order Upstream Advection'
          WRITE (iunit,30) 'Upstream3',                                 &
     &          'Third-order Upstream-biased Advection '//              &
     &          '(only horizontal)'
          WRITE (iunit,'(1x)')
        END IF
      END IF
!
!  Check switches for consistency.
!
      IF (model.eq.iNLM) THEN
        DO ng=1,Ngrids
          DO i=1,NT(ng)
            IF (.not.Vadv(i,ng)%MPDATA.and.Hadv(i,ng)%MPDATA) THEN
              IF (Master) THEN
                WRITE (iunit,40) TRIM(Vname(1,idTvar(i))), ng,          &
     &                           TRIM(Hstring(i,ng)),                   &
     &                'must be specified for both advective terms'
              END IF
              exit_flag=5
              RETURN
            ELSE IF (.not.Hadv(i,ng)%MPDATA.and.Vadv(i,ng)%MPDATA) THEN
              IF (Master) THEN
                WRITE (iunit,40) TRIM(Vname(1,idTvar(i))), ng,          &
     &                           TRIM(Vstring(i,ng)),                   &
     &                'must be specified for both advective terms'
              END IF
              exit_flag=5
              RETURN
            ELSE IF (.not.Vadv(i,ng)%HSIMT.and.Hadv(i,ng)%HSIMT) THEN
              IF (Master) THEN
                WRITE (iunit,40) TRIM(Vname(1,idTvar(i))), ng,          &
     &                           TRIM(Hstring(i,ng)),                   &
     &                'must be specified for both advective terms'
              END IF
              exit_flag=5
              RETURN
            ELSE IF (.not.Hadv(i,ng)%HSIMT.and.Vadv(i,ng)%HSIMT) THEN
              IF (Master) THEN
                WRITE (iunit,40) TRIM(Vname(1,idTvar(i))), ng,          &
     &                           TRIM(Vstring(i,ng)),                   &
     &                'must be specified for both advective terms'
              END IF
              exit_flag=5
              RETURN
            ELSE IF (Hadv(i,ng)%SPLINES) THEN
              IF (Master) THEN
                WRITE (iunit,40) TRIM(Vname(1,idTvar(i))), ng,          &
     &                           TRIM(Hstring(i,ng)),                   &
     &                'is only available for the vertical term'

              END IF
              exit_flag=5
              RETURN
            ELSE IF (Vadv(i,ng)%UPSTREAM3) THEN
              IF (Master) THEN
                WRITE (iunit,40) TRIM(Vname(1,idTvar(i))), ng,          &
     &                           TRIM(Vstring(i,ng)),                   &
     &                'is only available for the horizontal term'

              END IF
              exit_flag=5
              RETURN
            END IF
          END DO
        END DO
      ELSE IF (model.eq.iADM) THEN
        DO ng=1,Ngrids
          DO i=1,NT(ng)
            IF (Hadv(i,ng)%MPDATA.or.Hadv(i,ng)%HSIMT) THEN
              IF (Master) THEN
                WRITE (iunit,40) TRIM(Vname(1,idTvar(i))), ng,          &
     &                           TRIM(Hstring(i,ng)),                   &
     &                'is not supported in adjoint-based algorithms'
              END IF
              exit_flag=5
              RETURN
            ELSE IF (Vadv(i,ng)%MPDATA.or.Vadv(i,ng)%HSIMT) THEN
              IF (Master) THEN
                WRITE (iunit,40) TRIM(Vname(1,idTvar(i))), ng,          &
     &                           TRIM(Vstring(i,ng)),                   &
     &                'is not supported in adjoint-based algorithms'
              END IF
              exit_flag=5
              RETURN
            ELSE IF (Hadv(i,ng)%SPLINES) THEN
              IF (Master) THEN
                WRITE (iunit,40) TRIM(Vname(1,idTvar(i))), ng,          &
     &                           TRIM(Hstring(i,ng)),                   &
     &                'is only available for the vertical term'

              END IF
              exit_flag=5
              RETURN
            ELSE IF (Vadv(i,ng)%UPSTREAM3) THEN
              IF (Master) THEN
                WRITE (iunit,40) TRIM(Vname(1,idTvar(i))), ng,          &
     &                           TRIM(Vstring(i,ng)),                   &
     &                'is only available for the horizontal term'

              END IF
              exit_flag=5
              RETURN
            END IF
          END DO
        END DO
      END IF
!
 10   FORMAT (/,1x,a,t26,i2,t31,a,t50,a)
 20   FORMAT (t26,i2,t31,a,t50,a)
 30   FORMAT (1x,a,t13,a)
 40   FORMAT (/,'TADV_REPORT - Illegal tracer advection scheme for ''', &
     &        a,''' in grid: ',i0,/,14x,'''',a,'''',1x,a,'.',/)
!
      RETURN
      END SUBROUTINE tadv_report
#endif
      END MODULE tadv_mod
