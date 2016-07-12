#include "cppdefs.h"
#ifdef FLOATS
      SUBROUTINE read_FltPar (model, inp, out, Lwrite)
!
!svn $Id: read_fltpar.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads and reports floats input parameters.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_floats
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      logical :: find_file

      integer :: Npts, Nval
      integer :: i, j, igrid, mc, nc, ng, status

      integer, dimension(Ngrids) :: ncount, nentry

      integer, allocatable :: Fcoor(:,:), Fcount(:,:), Ftype(:,:)

      integer :: decode_line, load_i, load_l, load_r

      real(r8) :: xfloat, yfloat, zfloat

      real(r8), dimension(200) :: Rval

      real(r8), allocatable :: Ft0(:,:), Fx0(:,:), Fy0(:,:), Fz0(:,:)
      real(r8), allocatable :: Fdt(:,:), Fdx(:,:), Fdy(:,:), Fdz(:,:)

      character (len=1  ), parameter :: blank = ' '

      character (len=35 ) :: frmt
      character (len=40 ) :: KeyWord
      character (len=256) :: fname, line
      character (len=256), dimension(200) :: Cval
!
!-----------------------------------------------------------------------
!  Read in initial float locations.
!-----------------------------------------------------------------------
!
!  Allocate floats parameters that do not depend on the number of
!  floats, Nfloats(ng).
!
      CALL allocate_floats (.FALSE.)
!
!  Notice I added one when allocating local scratch arrays to avoid
!  out of bounds in some compilers when reading the last blank line
!  which signal termination of input data.
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=20,END=30) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Lfloats')
              Npts=load_l(Nval, Cval, Ngrids, Lfloats)
            CASE ('Fprint')
              Npts=load_l(Nval, Cval, Ngrids, Fprint)
            CASE ('FRREC')
              Npts=load_i(Nval, Rval, Ngrids, frrec)
            CASE ('FBIONAM')
              DO i=1,LEN(fbionam)
                fbionam(i:i)=blank
              END DO
              fbionam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('NFLOATS')
              Npts=load_i(Nval, Rval, Ngrids, Nfloats)
            CASE ('POS')
              Npts=Nfloats(1)+1
              IF (Ngrids.gt.1) Npts=MAXVAL(Nfloats)+1
              IF (.not.allocated(Fcoor)) THEN
                allocate ( Fcoor (Npts,Ngrids) )
              END IF
              IF (.not.allocated(Fcount)) THEN
                allocate ( Fcount(Npts,Ngrids) )
              END IF
              IF (.not.allocated(Ftype)) THEN
                allocate ( Ftype (Npts,Ngrids) )
              END IF
              IF (.not.allocated(Ft0)) THEN
                allocate ( Ft0(Npts,Ngrids) )
              END IF
              IF (.not.allocated(Fx0)) THEN
                allocate ( Fx0(Npts,Ngrids) )
              END IF
              IF (.not.allocated(Fy0)) THEN
                allocate ( Fy0(Npts,Ngrids) )
              END IF
              IF (.not.allocated(Fz0)) THEN
                allocate ( Fz0(Npts,Ngrids) )
              END IF
              IF (.not.allocated(Fdt)) THEN
                allocate ( Fdt(Npts,Ngrids) )
              END IF
              IF (.not.allocated(Fdx)) THEN
                allocate ( Fdx(Npts,Ngrids) )
              END IF
              IF (.not.allocated(Fdy)) THEN
                allocate ( Fdy(Npts,Ngrids) )
              END IF
              IF (.not.allocated(Fdz)) THEN
                allocate ( Fdz(Npts,Ngrids) )
              END IF
              CALL allocate_floats (.TRUE.) ! allocate DRIFTER structure
              ncount(1:Ngrids)=0
              nentry(1:Ngrids)=0
              DO WHILE (.TRUE.)
                READ (inp,*,ERR=30,END=30) igrid,                       &
     &                                   Fcoor (nentry(igrid)+1,igrid), &
     &                                   Ftype (nentry(igrid)+1,igrid), &
     &                                   Fcount(nentry(igrid)+1,igrid), &
     &                                   Ft0(nentry(igrid)+1,igrid),    &
     &                                   Fx0(nentry(igrid)+1,igrid),    &
     &                                   Fy0(nentry(igrid)+1,igrid),    &
     &                                   Fz0(nentry(igrid)+1,igrid),    &
     &                                   Fdt(nentry(igrid)+1,igrid),    &
     &                                   Fdx(nentry(igrid)+1,igrid),    &
     &                                   Fdy(nentry(igrid)+1,igrid),    &
     &                                   Fdz(nentry(igrid)+1,igrid)
                IF (igrid.gt.Ngrids) THEN
                  IF (Lwrite) WRITE (out,60) fposnam
                  exit_flag=4
                  RETURN
                END IF
                ncount(igrid)=ncount(igrid)+                            &
     &                        Fcount(nentry(igrid)+1,igrid)
                nentry(igrid)=nentry(igrid)+1
              END DO
          END SELECT
        END IF
      END DO
  20  IF (Master) WRITE (out,70) line
      exit_flag=4
      RETURN
  30  CONTINUE
!
!  Turn off the processing of floats if not running long enough to
!  create a floats file (LdefFLT=.FALSE. because nFLT < ntimes or
!  nFLT = 0 when nrrec = 0).
!
      DO ng=1,Ngrids
        IF (.not.LdefFLT(ng).and.Lfloats(ng)) THEN
          Lfloats(ng)=.FALSE.
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lfloats(ng).and.Fprint(ng)) THEN
            IF (ncount(ng).ne.Nfloats(ng)) THEN
              WRITE (stdout,80) ncount(ng), Nfloats(ng)
              exit_flag=4
              RETURN
            END IF
            WRITE (out,90) ng
            DO i=1,nentry(ng)
              IF (.not.spherical.and.(Fcoor(i,ng).eq.0)) THEN
                frmt='(i1,i2,i5,f10.4,2f8.2,f8.2,4f9.3)'
              ELSE
                frmt='(i1,i2,i5,f10.4,3f8.2,4f9.3)'
              END IF
              WRITE (out,frmt) Fcoor(i,ng), Ftype(i,ng), Fcount(i,ng),  &
     &                         Ft0(i,ng), Fx0(i,ng), Fy0(i,ng),         &
     &                         Fz0(i,ng), Fdt(i,ng), Fdx(i,ng),         &
     &                         Fdy(i,ng), Fdz(i,ng)
            END DO
            WRITE (out,100) Nfloats(ng),                                &
     &            'Nfloats',                                            &
     &            'Number of float trajectories to compute.'
          END IF
        END DO
      END IF

# ifdef FLOAT_BIOLOGY
      fname=fbionam
      IF (.not.find_file(1, fname, 'FBIONAM')) GO TO 40
      GO TO 50
 40   IF (Master) WRITE (out,110) TRIM(fname)
      exit_flag=4
      RETURN
 50   CONTINUE
      IF (Lwrite) WRITE (out,120) 'biological floats behavior File:  ', &
     &                            TRIM(fname)
# endif
!
!-----------------------------------------------------------------------
!  Process initial float locations.
!-----------------------------------------------------------------------
!
!  Set time of float release (seconds after model initialization) and
!  initial float horizontal positions (grid units).  Fill the initial
!  vertical level or depth position.
!
      DO ng=1,Ngrids
        mc=0
        nc=0
        IF (Lfloats(ng)) THEN
          DO i=1,nentry(ng)
            IF (Fcount(i,ng).eq.1) THEN
              nc=nc+1
              DRIFTER(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng))*day2sec
              DRIFTER(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
              DRIFTER(ng)%Ftype(nc)=Ftype(i,ng)
              IF (Fcoor(i,ng).eq.0) THEN
                DRIFTER(ng)%Tinfo(ixgrd,nc)=Fx0(i,ng)
                DRIFTER(ng)%Tinfo(iygrd,nc)=Fy0(i,ng)
              ELSE
                mc=mc+1
                DRIFTER(ng)%Flon(mc)=Fx0(i,ng)
                DRIFTER(ng)%Flat(mc)=Fy0(i,ng)
                DRIFTER(ng)%Findex(mc)=nc
              END IF
            ELSE IF (Fcount(i,ng).gt.1) THEN
              DO j=1,Fcount(i,ng)
                nc=nc+1
                IF (Fdt(i,ng).gt.0.0_r8) THEN
                  DRIFTER(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng)+        &
     &                                         REAL(j-1,r8)*Fdt(i,ng))* &
     &                                        day2sec
                  DRIFTER(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
                  DRIFTER(ng)%Ftype(nc)=Ftype(i,ng)
                  IF (Fcoor(i,ng).eq.0) THEN
                    DRIFTER(ng)%Tinfo(ixgrd,nc)=Fx0(i,ng)
                    DRIFTER(ng)%Tinfo(iygrd,nc)=Fy0(i,ng)
                  ELSE
                    mc=mc+1
                    DRIFTER(ng)%Flon(mc)=Fx0(i,ng)
                    DRIFTER(ng)%Flat(mc)=Fy0(i,ng)
                    DRIFTER(ng)%Findex(mc)=nc
                  END IF
                ELSE
                  DRIFTER(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng))*day2sec
                  IF (Fdz(i,ng).eq.0.0_r8) THEN
                    DRIFTER(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
                  ELSE
                    IF (Fz0(i,ng).gt.0.0_r8) THEN
                      zfloat=Fz0(i,ng)+REAL(j-1,r8)*Fdz(i,ng)
                      DRIFTER(ng)%Tinfo(izgrd,nc)=MIN(MAX(0.0_r8,       &
     &                                                    zfloat),      &
     &                                                REAL(N(ng),r8))
                    ELSE
                      DRIFTER(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)+            &
     &                                            REAL(j-1,r8)*Fdz(i,ng)
                    END IF
                  END IF
                  DRIFTER(ng)%Ftype(nc)=Ftype(i,ng)
                  IF (Fcoor(i,ng).eq.0) THEN
                    xfloat=Fx0(i,ng)+REAL(j-1,r8)*Fdx(i,ng)
                    yfloat=Fy0(i,ng)+REAL(j-1,r8)*Fdy(i,ng)
                    DRIFTER(ng)%Tinfo(ixgrd,nc)=xfloat
                    DRIFTER(ng)%Tinfo(iygrd,nc)=yfloat
                  ELSE
                    mc=mc+1
                    DRIFTER(ng)%Flon(mc)=Fx0(i,ng)+                     &
     &                                   REAL(j-1,r8)*Fdx(i,ng)
                    DRIFTER(ng)%Flat(mc)=Fy0(i,ng)+                     &
     &                                   REAL(j-1,r8)*Fdy(i,ng)
                    DRIFTER(ng)%Findex(mc)=nc
                  END IF
                END IF
              END DO
            END IF
          END DO
          DRIFTER(ng)%Findex(0)=mc
        END IF
      END DO
!
!  Deallocate local arrays.
!
      IF (allocated(Fcoor))  deallocate ( Fcoor )
      IF (allocated(Fcount)) deallocate ( Fcount )
      IF (allocated(Ftype))  deallocate ( Ftype )
      IF (allocated(Ft0))    deallocate ( Ft0 )
      IF (allocated(Fx0))    deallocate ( Fx0 )
      IF (allocated(Fy0))    deallocate ( Fy0 )
      IF (allocated(Fz0))    deallocate ( Fz0 )
      IF (allocated(Fdt))    deallocate ( Fdt )
      IF (allocated(Fdx))    deallocate ( Fdx )
      IF (allocated(Fdy))    deallocate ( Fdy )
      IF (allocated(Fdz))    deallocate ( Fdz )
!
  60  FORMAT (/,' READ_FltPar - Error while reading floats',            &
     &          ' locations in input script: ',a)
  70  FORMAT (/,' READ_FltPar - Error while processing line: ',/,a)
  80  FORMAT (/,' READ_FltPar - Inconsistent number of floats to',      &
     &          ' process: ', 2i6,/,18x,'change input script.')
  90  FORMAT (/,/,' Floats Initial Locations, Grid: ',i2.2,             &
     &        /,  ' ==================================',/,/,            &
     &        15x,'Ft0',5x,'Fx0',5x,'Fy0',5x,'Fz0',                     &
     &        6x,'Fdt',6x,'Fdx',6x,'Fdy',6x,'Fdz',/)
 100  FORMAT (/,1x,i10,2x,a,t30,a)
 110  FORMAT (/,' READ_FltPar - could not find input file: ', a)
 120  FORMAT (/,2x,a,a)

      RETURN
      END SUBROUTINE read_FltPar
#else
      SUBROUTINE read_FltPar
      END SUBROUTINE read_FltPar
#endif
