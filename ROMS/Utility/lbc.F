#include "cppdefs.h"
!
!svn $Id: lbc.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines are used to process lateral boundary condition       !
!  structure:                                                          !
!                                                                      !
!     lbc_getatt    Reads activated keyword strings from specified     !
!                     input NetCDF file and checks input script        !
!                     values during restart for consistency.           !
!                                                                      !
!     lbc_putatt    Writes activated keyword strings into specified    !
!                     output NetCDF file global attribute.             !
!                                                                      !
!     lbc_report    Reports to standard output activated keyword       !
!                     strings.                                         !
!                                                                      !
!=======================================================================
!
!
!***********************************************************************
      SUBROUTINE lbc_getatt (ng, model, ncid, ncname, aname, S)
!***********************************************************************
!                                                                      !
!  If restart,  this routine reads lateral boundary conditions         !
!  keywords strings from restart NetCDF file global attribute.         !
!  Then, it checks their values with those provided from input         !
!  script for consistency.                                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng            Nested grid number (integer).                      !
!     model         Calling model identifier (integer).                !
!     ncid          NetCDF file ID (integer).                          !
!     ncname        NetCDF file name (character).                      !
!     aname         NetCDF global attribute name (character).          !
!     S             Derived type structure, TYPE(T_LBC)                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag     Error flag (integer) stored in MOD_SCALARS         !
!     ioerror       NetCDF return code (integer) stored in MOD_IOUNITS !
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
#if !defined PARALLEL_IN  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcasts
#endif
!
      implicit none
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
      character (*), intent(in) :: ncname
      character (*), intent(in) :: aname

      TYPE(T_LBC), intent(in) :: S(4,nLBCvar,Ngrids)
!
! Local variable declarations
!
      integer :: i, ibry, ie, ifield, is, ne, lstr, lvar, status
#if !defined PARALLEL_IN  && defined DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif

      character (len=   7) :: string(4)
      character (len=   8) :: B(4)
      character (len=  40) :: BryVar1, BryVar2
      character (len=  70) :: Bstring, line
      character (len=2816) :: lbc_att
!
!-----------------------------------------------------------------------
!  If restart, read and check lateral boundary conditions global
!  attribute.
!-----------------------------------------------------------------------
!
!  Read in global attribute.
!
      IF (InpThread) THEN
        status=nf90_get_att(ncid, nf90_global, TRIM(aname), lbc_att)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(aname), TRIM(ncname),                  &
     &                      TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IN  && defined DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasts (ng, model, lbc_att)
      END IF
#endif
!
!  Check keyword values from global attribute and compare against
!  logical switches in lateral boundary condition structure.
!
      B(iwest )='western '
      B(isouth)='southern'
      B(ieast )='eastern '
      B(inorth)='northern'
      DO i=1,LEN(Bstring)
        Bstring(i:i)=' '
      END DO
      DO i=1,LEN(line)
        line(i:i)=' '
      END DO
!
      DO ifield=1,nLBCvar
        BryVar1=TRIM(Vname(1,idBvar(ifield))) // ':'
        is=INDEX(lbc_att, TRIM(BryVar1))
        IF (ifield.lt.nLBCvar) THEN
          BryVar2=TRIM(Vname(1,idBvar(ifield+1))) // ':'
          ie=INDEX(lbc_att, TRIM(BryVar2))-1
        ELSE
          ie=LEN_TRIM(lbc_att)
        END IF
        IF ((is.gt.0).and.(ie.gt.0).and.(ie.gt.is)) THEN
          line=lbc_att(is:ie)
          is=INDEX(line, ':')+1               ! find colon location
          ie=INDEX(line, CHAR(10))-1          ! find Line Feed location
          IF (ie.le.0) THEN
            ie=LEN_TRIM(line)                 ! last attribute entry
          END IF
          Bstring=TRIM(ADJUSTL(line(is:ie)))  ! extract boundary string
          ne=MIN(LEN_TRIM(Bstring), 28)       ! north keyword < 7 chars
          string(1)=Bstring( 1: 7)            ! west  keyword
          string(2)=Bstring( 8:14)            ! south keyword
          string(3)=Bstring(15:21)            ! east  keyword
          string(4)=Bstring(22:ne)            ! north keyword
          DO ibry=1,4
            SELECT CASE (TRIM(string(ibry)))
              CASE ('Cha')
                IF (.not.S(ibry,ifield,ng)%Chapman_implicit) THEN
                  IF (Master) THEN
                    WRITE (stdout,20) B(ibry),                          &
     &                    TRIM(Vname(1,idBvar(ifield))),                &
     &                    TRIM(string(ibry)),                           &
     &                    'S(',ibry,ifield,ng,')%Chapman_implicit',     &
     &                    S(ibry,ifield,ng)%Chapman_implicit,           &
     &                    TRIM(ncname)
                  END IF
                  exit_flag=5
                END IF
              CASE ('Che')
                IF (.not.S(ibry,ifield,ng)%Chapman_explicit) THEN
                  IF (Master) THEN
                    WRITE (stdout,20) B(ibry),                          &
     &                    TRIM(Vname(1,idBvar(ifield))),                &
     &                    TRIM(string(ibry)),                           &
     &                    'S(',ibry,ifield,ng,')%Chapman_explicit',     &
     &                    S(ibry,ifield,ng)%Chapman_explicit,           &
     &                    TRIM(ncname)
                  END IF
                  exit_flag=5
                END IF

              CASE ('Cla')
                IF (.not.S(ibry,ifield,ng)%clamped) THEN
                  IF (Master) THEN
                    WRITE (stdout,20) B(ibry),                          &
     &                    TRIM(Vname(1,idBvar(ifield))),                &
     &                    TRIM(string(ibry)),                           &
     &                    'S(',ibry,ifield,ng,')%clamped',              &
     &                    S(ibry,ifield,ng)%clamped,                    &
     &                    TRIM(ncname)
                  END IF
                  exit_flag=5
                END IF
              CASE ('Clo')
                IF (.not.S(ibry,ifield,ng)%closed) THEN
                  IF (Master) THEN
                    WRITE (stdout,20) B(ibry),                          &
     &                    TRIM(Vname(1,idBvar(ifield))),                &
     &                    TRIM(string(ibry)),                           &
     &                    'S(',ibry,ifield,ng,')%closed',               &
     &                    S(ibry,ifield,ng)%closed,                     &
     &                    TRIM(ncname)
                  END IF
                  exit_flag=5
                END IF
              CASE ('Fla')
                IF (.not.S(ibry,ifield,ng)%Flather) THEN
                  IF (Master) THEN
                    WRITE (stdout,20) B(ibry),                          &
     &                    TRIM(Vname(1,idBvar(ifield))),                &
     &                    TRIM(string(ibry)),                           &
     &                    'S(',ibry,ifield,ng,')%Flather',              &
     &                    S(ibry,ifield,ng)%Flather,                    &
     &                    TRIM(ncname)
                  END IF
                  exit_flag=5
                END IF
              CASE ('Gra')
                IF (.not.S(ibry,ifield,ng)%gradient) THEN
                  IF (Master) THEN
                    WRITE (stdout,20) B(ibry),                          &
     &                    TRIM(Vname(1,idBvar(ifield))),                &
     &                    TRIM(string(ibry)),                           &
     &                    'S(',ibry,ifield,ng,')%gradient',             &
     &                    S(ibry,ifield,ng)%gradient,                   &
     &                    TRIM(ncname)
                  END IF
                  exit_flag=5
                END IF
              CASE ('Mix')
                IF (.not.S(ibry,ifield,ng)%mixed) THEN
                  IF (Master) THEN
                    WRITE (stdout,20) B(ibry),                          &
     &                                TRIM(Vname(1,idBvar(ifield))),    &
     &                                TRIM(string(ibry)),               &
     &                                'S(',ibry,ifield,ng,')%mixed',    &
     &                                S(ibry,ifield,ng)%mixed,          &
     &                                TRIM(ncname)
                  END IF
                  exit_flag=5
                END IF
              CASE ('Nes')
                IF (.not.S(ibry,ifield,ng)%nested) THEN
                  IF (Master) THEN
                    WRITE (stdout,20) B(ibry),                          &
     &                    TRIM(Vname(1,idBvar(ifield))),                &
     &                    TRIM(string(ibry)),                           &
     &                    'S(',ibry,ifield,ng,')%nested',               &
     &                    S(ibry,ifield,ng)%nested,                     &
     &                    TRIM(ncname)
                  END IF
                  exit_flag=5
                END IF
              CASE ('Per')
                IF (.not.S(ibry,ifield,ng)%periodic) THEN
                  IF (Master) THEN
                    WRITE (stdout,20) B(ibry),                          &
     &                    TRIM(Vname(1,idBvar(ifield))),                &
     &                    TRIM(string(ibry)),                           &
     &                    'S(',ibry,ifield,ng,')%periodic',             &
     &                    S(ibry,ifield,ng)%periodic,                   &
     &                    TRIM(ncname)
                  END IF
                  exit_flag=5
                END IF
              CASE ('Rad')
                IF (.not.S(ibry,ifield,ng)%radiation) THEN
                  IF (Master) THEN
                    WRITE (stdout,20) B(ibry),                          &
     &                    TRIM(Vname(1,idBvar(ifield))),                &
     &                    TRIM(string(ibry)),                           &
     &                    'S(',ibry,ifield,ng,')%radiation',            &
     &                    S(ibry,ifield,ng)%radiation,                  &
     &                    TRIM(ncname)
                  END IF
                  exit_flag=5
                END IF
              CASE ('RadNud')
                IF (.not.(S(ibry,ifield,ng)%radiation.and.              &
     &                    S(ibry,ifield,ng)%nudging)) THEN
                  IF (Master) THEN
                    WRITE (stdout,20) B(ibry),                          &
     &                    TRIM(Vname(1,idBvar(ifield))),                &
     &                    TRIM(string(ibry)),                           &
     &                    'S(',ibry,ifield,ng,')%radiation',            &
     &                    S(ibry,ifield,ng)%radiation,                  &
     &                    TRIM(ncname)
                  END IF
                  exit_flag=5
                END IF
              CASE ('Red')
                IF (.not.S(ibry,ifield,ng)%reduced) THEN
                  IF (Master) THEN
                    WRITE (stdout,20) B(ibry),                          &
     &                    TRIM(Vname(1,idBvar(ifield))),                &
     &                    TRIM(string(ibry)),                           &
     &                    'S(',ibry,ifield,ng,')%reduced',              &
     &                    S(ibry,ifield,ng)%reduced,                    &
     &                    TRIM(ncname)
                  END IF
                  exit_flag=5
                END IF
              CASE ('Shc')
                IF (.not.S(ibry,ifield,ng)%Shchepetkin) THEN
                  IF (Master) THEN
                    WRITE (stdout,20) B(ibry),                          &
     &                    TRIM(Vname(1,idBvar(ifield))),                &
     &                    TRIM(string(ibry)),                           &
     &                    'S(',ibry,ifield,ng,')%Shchepetkin',          &
     &                    S(ibry,ifield,ng)%Shchepetkin,                &
     &                    TRIM(ncname)
                  END IF
                  exit_flag=5
                END IF
              CASE DEFAULT
                IF (Master) THEN
                  WRITE (stdout,30) B(ibry),                            &
     &                              TRIM(Vname(1,idBvar(ifield))),      &
     &                              TRIM(string(ibry)), TRIM(ncname)
                END IF
                exit_flag=5
            END SELECT
          END DO
        END IF
      END DO

  10  FORMAT (/,' LBC_GETATT - error while reading global attribute:',  &
     &        2x,a,/,14x,'in restart file:',2x,a,/,14x,'call from:',    &
     &        2x,a,                                                     &
     &        /,14x,'Probably global attribute was not found ...',      &
     &        /,14x,'restart file needs to be generated by ROMS ',      &
     &        'version 3.6 or higher',                                  &
     &        /,14x,'Alternatively, you may use  NO_LBC_ATT  at your '  &
     &        'own risk!')
  20  FORMAT (/,' LBC_GETATT - inconsistent ',a,' lateral boundary ',   &
     &        'condition for variable: ',2x,a,                          &
     &        /,14x,'restart file LBC keyword = ',1x,a,                 &
     &        /,14x,'but assigned structure switch: ',                  &
     &         1x,a,i1,',',i2,',',i1,a,' = ',l1,                        &
     &        /,14x,'check input script LBC keyword for consitency ...',&
     &        /,14x,'restart file:',2x,a)
  30  FORMAT (/,' LBC_GETATT - inconsistent ',a,' boundary for ',       &
     &        'variable: ',a,2x,'Keyword = ',a,/,14x,'in input file:',  &
     &        2x,a)

      RETURN
      END SUBROUTINE lbc_getatt
!
!***********************************************************************
      SUBROUTINE lbc_putatt (ng, ncid, ncname, aname, S, status)
!***********************************************************************
!                                                                      !
!  This routine writes lateral boundary conditions keywords strings    !
!  into specified output NetCDF file global attribute.                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng            Nested grid number (integer).                      !
!     ncid          NetCDF file ID (integer).                          !
!     ncname        NetCDF file name (character).                      !
!     aname         NetCDF global attribute name (character).          !
!     S             Derived type structure, TYPE(T_LBC)                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag     Error flag (integer) stored in MOD_SCALARS         !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
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
      implicit none
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, ncid
      character (*), intent(in) :: ncname
      character (*), intent(in) :: aname

      TYPE(T_LBC), intent(in) :: S(4,nLBCvar,Ngrids)

      integer, intent(out) :: status
!
! Local variable declarations
!
      integer :: i, ibry, ie, ifield, is, lstr, lvar

      character (len=   1) :: newline
      character (len=   7) :: string(4)
      character (len=  21) :: frmt
      character (len=  70) :: line
      character (len=2816) :: lbc_att
!
!-----------------------------------------------------------------------
!  Write lateral boundary conditions global attribute.
!-----------------------------------------------------------------------
!
!  Determine maximum length of state variable length.
!
      lvar=0
      DO ifield=1,nLBCvar
        lvar=MAX(lvar, LEN_TRIM(Vname(1,idBvar(ifield))))
      END DO
      WRITE (frmt,10) lvar+4
  10  FORMAT ("(a,':',t",i2.2,',a,a,a,a,a)')
!
!  Initialize attribute.
!
      newline=CHAR(10)            ! Line Feed (LF) character for
      lstr=LEN_TRIM(newline)      ! attribute clarity with "ncdump"
      DO i=1,LEN(lbc_att)
        lbc_att(i:i)=' '
      END DO
      lbc_att(1:lstr)=newline(1:lstr)
      is=lstr+1
      WRITE (line,frmt) 'EDGE',                                         &
     &                  'WEST   ',                                      &
     &                  'SOUTH  ',                                      &
     &                  'EAST   ',                                      &
     &                  'NORTH  ',                                      &
     &                   newline(1:lstr)
      lstr=LEN_TRIM(line)
      ie=is+lstr
      lbc_att(is:ie)=line(1:lstr)
      is=ie
!
!  Check if the local string "lbc_att" is big enough to store the
!  lateral boundary conditions global attribute.
!
      lstr=(nLBCvar+1)*(29+lvar+4)+1
      IF (LEN(lbc_att).lt.lstr) THEN
        IF (Master) THEN
          WRITE (stdout,20) LEN(lbc_att), lstr
  20      FORMAT (/,' LBC_PUTATT - Length of local string lbc_att ',    &
     &            'too small',/,14x,'Current = ',i5,'  Needed = ',i5)
        END IF
        exit_flag=5
        RETURN
      END IF
!
!  Build attribute string.
!
      DO ifield=1,nLBCvar
        DO ibry=1,4
          IF (S(ibry,ifield,ng)%Chapman_explicit) THEN
            string(ibry)='Che    '
          ELSE IF (S(ibry,ifield,ng)%Chapman_implicit) THEN
            string(ibry)='Cha    '
          ELSE IF (S(ibry,ifield,ng)%clamped) THEN
            string(ibry)='Cla    '
          ELSE IF (S(ibry,ifield,ng)%closed) THEN
            string(ibry)='Clo    '
          ELSE IF (S(ibry,ifield,ng)%Flather) THEN
            string(ibry)='Fla    '
          ELSE IF (S(ibry,ifield,ng)%gradient) THEN
            string(ibry)='Gra    '
          ELSE IF (S(ibry,ifield,ng)%mixed) THEN
            string(ibry)='Mix    '
          ELSE IF (S(ibry,ifield,ng)%nested) THEN
            string(ibry)='Nes    '
          ELSE IF (S(ibry,ifield,ng)%periodic) THEN
            string(ibry)='Per    '
          ELSE IF (S(ibry,ifield,ng)%radiation) THEN
            IF (S(ibry,ifield,ng)%nudging) THEN
              string(ibry)='RadNud '
            ELSE
              string(ibry)='Rad    '
            END IF
          ELSE IF (S(ibry,ifield,ng)%reduced) THEN
            string(ibry)='Red    '
          ELSE IF (S(ibry,ifield,ng)%Shchepetkin) THEN
            string(ibry)='Shc    '
          END IF
        END DO
        IF (ifield.eq.nLBCvar) newline=' '
        WRITE (line,frmt) TRIM(Vname(1,idBvar(ifield))),                &
     &                    string(iwest),                                &
     &                    string(isouth),                               &
     &                    string(ieast),                                &
     &                    string(inorth),                               &
     &                    newline
        lstr=LEN_TRIM(line)
        ie=is+lstr
        lbc_att(is:ie)=line(1:lstr)
        is=ie
      END DO
!
!  Write attribute to NetCDF file.
!
      status=nf90_put_att(ncid, nf90_global, TRIM(aname),               &
     &                    TRIM(lbc_att))

      RETURN
      END SUBROUTINE lbc_putatt
!
!***********************************************************************
      SUBROUTINE lbc_report (iunit, ifield, S)
!***********************************************************************
!                                                                      !
!  This routine reports lateral boundary conditions to requested       !
!  state variable.                                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     iunit         Output logical unit (integer)                      !
!     ifield        Lateral boundary variable index (integer)          !
!     S             Derived type structure, TYPE(T_LBC)                !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
! Imported variable declarations.
!
      integer, intent(in) :: ifield, iunit

      TYPE(T_LBC), intent(in) :: S(4,nLBCvar,Ngrids)
!
! Local variable declarations
!
      integer :: ibry, ng

      character (len=11) :: string(4,Ngrids)
!
!-----------------------------------------------------------------------
!  Report lateral boundary conditions.
!-----------------------------------------------------------------------
!
!  Build information strings.
!
      DO ng=1,Ngrids
        DO ibry=1,4
          IF (S(ibry,ifield,ng)%Chapman_explicit) THEN
            string(ibry,ng)='Chapman Exp'
          ELSE IF (S(ibry,ifield,ng)%Chapman_implicit) THEN
            string(ibry,ng)='Chapman Imp'
          ELSE IF (S(ibry,ifield,ng)%clamped) THEN
            string(ibry,ng)='Clamped    '
          ELSE IF (S(ibry,ifield,ng)%closed) THEN
            string(ibry,ng)='Closed     '
          ELSE IF (S(ibry,ifield,ng)%Flather) THEN
            string(ibry,ng)='Flather    '
          ELSE IF (S(ibry,ifield,ng)%gradient) THEN
            string(ibry,ng)='Gradient   '
          ELSE IF (S(ibry,ifield,ng)%mixed) THEN
            string(ibry,ng)='Mixed      '
          ELSE IF (S(ibry,ifield,ng)%nested) THEN
            string(ibry,ng)='Nested     '
          ELSE IF (S(ibry,ifield,ng)%periodic) THEN
            string(ibry,ng)='Periodic   '
          ELSE IF (S(ibry,ifield,ng)%radiation) THEN
            IF (S(ibry,ifield,ng)%nudging) THEN
              string(ibry,ng)='Rad + Nud  '
            ELSE
              string(ibry,ng)='Radiation  '
            END IF
          ELSE IF (S(ibry,ifield,ng)%reduced) THEN
            string(ibry,ng)='Reduced    '
          ELSE IF (S(ibry,ifield,ng)%Shchepetkin) THEN
            string(ibry,ng)='Shchepetkin'
          END IF
        END DO
      END DO
!
!  Report.
!
      DO ng=1,Ngrids
        IF (ng.eq.1) THEN
          WRITE (iunit,10) TRIM(Vname(1,idBvar(ifield))), ng,           &
     &                     TRIM(string(1,ng)),                          &
     &                     TRIM(string(2,ng)),                          &
     &                     TRIM(string(3,ng)),                          &
     &                     TRIM(string(4,ng))
        ELSE
          WRITE (iunit,20) ng,                                          &
     &                     TRIM(string(1,ng)),                          &
     &                     TRIM(string(2,ng)),                          &
     &                     TRIM(string(3,ng)),                          &
     &                     TRIM(string(4,ng))
        END IF
      END DO
!
!  Check periodic bounadry conditions: both east and west or north and
!  need to be specified.
!
      DO ng=1,Ngrids
        IF (.not.S(iwest,ifield,ng)%periodic.and.                       &
     &           S(ieast,ifield,ng)%periodic) THEN
          WRITE (iunit,30) 'Western Edge boundary',                     &
     &                     TRIM(Vname(1,idBvar(ifield))), ng
          exit_flag=5
        ELSE IF (.not.S(ieast,ifield,ng)%periodic.and.                  &
     &                S(iwest,ifield,ng)%periodic) THEN
          WRITE (iunit,30) 'Eastern Edge boundary',                     &
     &                     TRIM(Vname(1,idBvar(ifield))), ng
          exit_flag=5
        ELSE IF (.not.S(inorth,ifield,ng)%periodic.and.                 &
     &                S(isouth,ifield,ng)%periodic) THEN
          WRITE (iunit,30) 'Northern Edge boundary',                    &
     &                     TRIM(Vname(1,idBvar(ifield))), ng
          exit_flag=5
        ELSE IF (.not.S(isouth,ifield,ng)%periodic.and.                 &
     &                S(inorth,ifield,ng)%periodic) THEN
          WRITE (iunit,30) 'Southern Edge boundary',                    &
     &                     TRIM(Vname(1,idBvar(ifield))), ng
          exit_flag=5
        END IF
      END DO

 10   FORMAT (/,1x,a,t26,i2,t31,a,t44,a,t57,a,t70,a)
 20   FORMAT (t26,i2,t31,a,t44,a,t57,a,t70,a)
 30   FORMAT (/,' LBC_REPORT - illegal configuration: The ',a,          &
     &        ' needs to be periodic too!',/,14x,'Variable: ',a,3x,     &
     &        'Grid = ',i2.2)

      RETURN
      END SUBROUTINE lbc_report
