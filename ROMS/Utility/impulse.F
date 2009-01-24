#include "cppdefs.h"
      MODULE impulse_mod

#if defined ADJOINT && defined IMPULSE
!
!svn $Id: impulse.F 584 2008-03-18 20:01:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutines read in requested adjoint solution from input      !
!  NetCDF (saved at nADJ time-step intervals), convert to impulse      !
!  forcing, and then writes to output  NetCDF  file in  ascending      !
!  time order since it is processed by the TL and RP models.           !
!                                                                      !
!  The impulse forcing, F,  is computed from snapshots of adjoint      !
!  solution, ad_F, as follows:                                         !
!                                                                      !
!  n     0        1        2        3              N-1        N        !
!        |________|________|________|___......._____|_________|        !
!  iic   0      1*nADJ   2*nADJ   3*nADJ      (N-1)*nADJ    N*nADJ     !
!                                                                      !
!  They are only computed at the time-steps where the  solution is     !
!  saved:                                                              !
!                                                                      !
!  At the first time-step, iic=ntdtstr=0, n=0:                         !
!                                                                      !
!        F(0) = C1 * ad_F(0) + C2 * ad_F(1)                            !
!                                                                      !
!  at the last time-step, iic=ntend=N*nADJ, n=N:                       !
!                                                                      !
!        F(N) = C2 * ad_F(N-1) + C1 * ad_F(N)                          !
!                                                                      !
!  otherwise,                                                          !
!                                                                      !
!        F(n) = C2 * ad_F(n-1) + C3 * Ad_F(n) + C2 * ad_F(n+1)         !
!                                                                      !
!  where,                                                              !
!                                                                      !
!           N = number of NetCDF records to process.                   !
!                                                                      !
!        nADJ = number of time-steps between adjoint snapshots         !
!                                                                      !
!          C1 = (nADJ + 1) * (2 * nADJ + 1) / (6 * nADJ)               !
!                                                                      !
!          C2 = (nADJ * nADJ - 1) / (6 * nADJ)                         !
!                                                                      !
!          C3 = (2 * nADJ * nADJ + 1) / (3 * nADJ)                     !
!                                                                      !
!======================================================================!
!
      implicit none

      PRIVATE
      PUBLIC :: impulse

      CONTAINS

      SUBROUTINE impulse (ng, tile, model, outer_impulse, INPncname)
!
!=======================================================================
!                                                                      !
!  This subroutine reads in requested adjoint solution from input      !
!  NetCDF, converts to impulse forcing, and then writes to output      !
!  NetCDF file.                                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng              Nested grid number.                              !
!     tile            Domain partition.                                !
!     model           Calling model identifier.                        !
!     outer_impulse   Switch to process outer loop impulses.           !
!     INPncname       Input adjoint solution NetCDF file name.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: outer_impulse

      integer, intent(in) :: ng, tile, model

      character (len=*), intent(in) :: INPncname
!
!  Local variable declarations.
!
      logical, dimension(NV) :: get_var
      logical, dimension(NV) :: have_var

      integer :: LBi, UBi, LBj, UBj
      integer :: INPncid, INPvid
      integer :: Nrec, Nsave, gtype, status, varid
      integer :: i
      integer :: save(1), start(4), total(4)

      real(r8) :: inp_time(1)

      character (len=6 ) :: string
      character (len=15) :: tvarnam

# include "set_bounds.h"
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Determine variables to read and their availability.
!-----------------------------------------------------------------------
!
!  Set model identification string.
!
      string=' FRC: '
!
!  Determine variables to read.
!
      CALL checkvars (ng, model, INPncname, string, Nrec, NV, tvarnam,  &
     &                get_var, have_var)
!
!-----------------------------------------------------------------------
!  Read adjoint solution and convert to impulse forcing.  Then, write
!  impulse forcing into output NetCDF file.
!-----------------------------------------------------------------------
!
!  Open input NetCDF file.
!
      IF (InpThread) THEN
        status=nf90_open(TRIM(INPncname), nf90_nowrite, INPncid)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) string, TRIM(INPncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Read in input file time-step size and number of time-steps between
!  snapshots.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(INPncid, 'nADJ', varid)
        start(1)=ng
        total(1)=1
        status=nf90_get_var(INPncid, varid, save, start, total)       
        Nsave=save(1)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'nADJ', TRIM(INPncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, Nsave, 1)
# endif
!
!  Process free-surface impulse.
!
      IF (get_var(idFsur)) THEN
        IF (InpThread) THEN
          status=nf90_inq_varid(INPncid, TRIM(Vname(1,idFsur)), INPvid)
          gtype=vflag(varid)*r2dvar
        END IF
        CALL impulse2d (ng, model, gtype, outer_impulse, Nrec, Nsave,   &
     &                  INPncname, INPncid, INPvid,                     &
     &                  TLFname(ng), ncTLFid(ng), tlfVid(idFsur,ng),    &
     &                  IstrR, IendR, JstrR, JendR,                     &
     &                  LBi, UBi, LBj, UBj,                             &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
     &                  OCEAN(ng) % ad_zeta)
        IF (exit_flag.ne.NoError) RETURN
      END IF

# ifndef SOLVE3D
!
!  Process 2D U-momentum component impulse.
!
      IF (get_var(idUbar)) THEN
        IF (InpThread) THEN
          status=nf90_inq_varid(INPncid, TRIM(Vname(1,idUbar)), INPvid)
          gtype=vflag(INPvid)*u2dvar
        END IF
        CALL impulse2d (ng, model, gtype, outer_impulse, Nrec, Nsave,   &
     &                  INPncname, INPncid, INPvid,                     &
     &                  TLFname(ng), ncTLFid(ng), tlfVid(idUbar,ng),    &
     &                  Istr, IendR, JstrR, JendR,                      &
     &                  LBi, UBi, LBj, UBj,                             &
#  ifdef MASKING
     &                  GRID(ng) % umask,                               &
#  endif
     &                  OCEAN(ng) % ad_ubar)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Process 2D V-momentum component impulse.
!
      IF (get_var(idVbar)) THEN
        IF (InpThread) THEN
          status=nf90_inq_varid(INPncid, TRIM(Vname(1,idVbar)), INPvid)
          gtype=vflag(varid)*v2dvar
        END IF
        CALL impulse2d (ng, model, gtype, outer_impulse, Nrec, Nsave,   &
     &                  INPncname, INPncid, INPvid,                     &
     &                  TLFname(ng), ncTLFid(ng), tlfVid(idVbar,ng),    &
     &                  IstrR, IendR, Jstr, JendR,                      &
     &                  LBi, UBi, LBj, UBj,                             &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  OCEAN(ng) % ad_vbar)
        IF (exit_flag.ne.NoError) RETURN
      END IF

# else
!
!  Process 3D U-momentum component impulse. Notice that "ad_t" is used
!  here as a temporary array since it has three-time levels.
!
      IF (get_var(idUvel)) THEN
        IF (InpThread) THEN
          status=nf90_inq_varid(INPncid, TRIM(Vname(1,idUvel)), INPvid)
          gtype=vflag(varid)*u3dvar
        END IF
        CALL impulse3d (ng, model, gtype, outer_impulse, Nrec, Nsave,   &
     &                  INPncname, INPncid, INPvid,                     &
     &                  TLFname(ng), ncTLFid(ng), tlfVid(idUvel,ng),    &
     &                  Istr, IendR, JstrR, JendR,                      &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
#  ifdef MASKING
     &                  GRID(ng) % umask,                               &
#  endif
     &                  OCEAN(ng) % ad_t(:,:,:,:,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Process 3D V-momentum component impulse. Notice that "ad_t" is used
!  here as a temporary array since it has three-time levels.
!
      IF (get_var(idVvel)) THEN
        IF (InpThread) THEN
          status=nf90_inq_varid(INPncid, TRIM(Vname(1,idVvel)), INPvid)
          gtype=vflag(varid)*v3dvar
        END IF
        CALL impulse3d (ng, model, gtype, outer_impulse, Nrec, Nsave,   &
     &                  INPncname, INPncid, INPvid,                     &
     &                  TLFname(ng), ncTLFid(ng), tlfVid(idVvel,ng),    &
     &                  IstrR, IendR, Jstr, JendR,                      &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  OCEAN(ng) % ad_t(:,:,:,:,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Process tracer type variables impulses.
!
      DO i=1,NT(ng)
        IF (get_var(idTvar(i))) THEN
          IF (InpThread) THEN
            status=nf90_inq_varid(INPncid, TRIM(Vname(1,idTvar(i))),    &
     &                            INPvid)
            gtype=vflag(varid)*r3dvar
          END IF
          CALL impulse3d (ng, model, gtype, outer_impulse, Nrec, Nsave, &
     &                    INPncname, INPncid, INPvid,                   &
     &                    TLFname(ng), ncTLFid(ng), tlfTid(i,ng),       &
     &                    IstrR, IendR, JstrR, JendR,                   &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % ad_t(:,:,:,:,i))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
# endif
!
!  Write out time in backward order.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(INPncid, TRIM(Vname(1,idtime)), INPvid)
        DO i=1,Nrec        
          start(1)=i
          total(1)=1
          status=nf90_get_var(INPncid,INPvid, inp_time, start, total)
          start(1)=Nrec-i+1
          total(1)=1
          status=nf90_put_var(ncTLFid(ng), tlfVid(idtime,ng),           &
     &                        inp_time, start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) TRIM(Vname(1,idtime)), Nrec-i+1,          &
     &                        TRIM(TLFname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize impulse NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      IF (InpThread) THEN
        status=nf90_sync(ncTLFid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,40)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        WRITE (stdout,50) Nrec, TRIM(TLFname(ng))
      END IF

  10  FORMAT (/,' IMPULSE - unable to open input NetCDF file: ',a)
  20  FORMAT (/,' IMPULSE - error while reading variable: ',a,          &
     &        /,12x,'in input NetCDF file: ',a)
  30  FORMAT (/,' IMPULSE - error while writing variable: ',a,2x,       &
     &        'at time record = ',i3,/,12x,'into NetCDF file: ',a)
  40  FORMAT (/,' IMPULSE - unable to synchronize impulse NetCDF to ',  &
     &        'disk.')
  50  FORMAT (6x,'IMPULSE   - wrote convolved adjoint impulses, ',      &
     &           'records: 001 to ',i3.3,/,18x,'file: ',a)
      
      RETURN
      END SUBROUTINE impulse

      SUBROUTINE impulse2d (ng, model, gtype, outer_impulse,            &
     &                      Nrec, Nsave,                                &
     &                      INPncname, INPncid, INPvid,                 &
     &                      OUTncname, OUTncid, OUTvid,                 &
     &                      Imin, Imax, Jmin, Jmax,                     &
     &                      LBi, UBi, LBj, UBj,                         &
# ifdef MASKING
     &                      mask,                                       &
# endif
     &                      ad_A)
!
!=======================================================================
!                                                                      !
!  This routine reads in requested 2D adjoint variable from input      !
!  NetCDF file, converts to impulse,  and then writes into output      !
!  NetCDF file.  The impulses are written in time ascending order      !
!  since they will be processed by forward TL and RP models.           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_netcdf
      USE mod_scalars
!
!  Imported variable declarations.
!
      logical, intent(in) :: outer_impulse

      integer, intent(in) :: ng, model, gtype, Nrec, Nsave
      integer, intent(in) :: INPncid, INPvid, OUTncid, OUTvid
      integer, intent(in) :: Imin, Imax, Jmin, Jmax
      integer, intent(in) :: LBi, UBi, LBj, UBj

      character(len=*), intent(in) :: INPncname
      character(len=*), intent(in) :: OUTncname
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: mask(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: mask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, status
      integer :: Irec, Iold, Inow, Inxt, Iout, Isav

      integer :: nf_fread2d
      integer :: nf_fwrite2d

      integer, dimension(4) :: Vsize(4)

      real(r8) :: C1, C2, C3, Fmin, Fmax, cff, scale

      real(r8), dimension(LBi:UBi,LBj:UBj) :: Awrk

      character(len=40) :: ncvname
!
!-----------------------------------------------------------------------
!  Initialize time level indices and impulse focing coefficients.
!-----------------------------------------------------------------------
!
      Inow=1
      Inxt=2
      Iold=3
      cff=REAL(Nsave,r8)
      IF (outer_impulse) THEN
        C1=(cff+1.0_r8)*(2.0_r8*cff+1.0_r8)/(6.0_r8*cff)
        C2=(cff*cff-1.0_r8)/(6.0_r8*cff)
        C3=(2.0_r8*cff*cff+1.0_r8)/(3.0_r8*cff)
      ELSE
        C1=1.0_r8
        C2=0.0_r8
        C3=1.0_r8
      END IF
      scale=1.0_r8
!
!  Set Vsize to zero to deactivate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Initialize work array.
!
      Awrk=0.0_r8
!
!-----------------------------------------------------------------------
!  Process first adjoint time record.
!-----------------------------------------------------------------------
!
      Irec=1
!
!  Read initial time record.
!
      status=nf_fread2d(ng, model, INPncid, INPvid, Irec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  mask(LBi,LBj),                                  &
# endif
     &                  ad_A(LBi,LBj,Inow))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          status=nf90_inquire_variable(INPncid, INPvid, ncvname)
          WRITE (stdout,10) TRIM(ncvname), Irec, TRIM(INPncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF
!
!  Read next time record.
!
      IF (Nrec.gt.1) THEN
        status=nf_fread2d(ng, model, INPncid, INPvid, Irec+1, gtype,    &
     &                    Vsize, LBi, UBi, LBj, UBj,                    &
     &                    scale, Fmin, Fmax,                            &
# ifdef MASKING
     &                    mask(LBi,LBj),                                &
# endif
     &                    ad_A(LBi,LBj,Inxt))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            status=nf90_inquire_variable(INPncid, INPvid, ncvname)
            WRITE (stdout,10) TRIM(ncvname), Irec, TRIM(INPncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      ELSE
        C2=0.0_r8
      END IF
!
!  Compute impulse forcing.
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          Awrk(i,j)=C1*ad_A(i,j,Inow)+                                  &
     &              C2*ad_A(i,j,Inxt)
        END DO
      END DO
!
!  Write out initial adjoint impulse forcing to last record.
!
      Iout=Nrec
      status=nf_fwrite2d(ng, model, OUTncid, OUTvid, Iout, gtype,       &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   mask(LBi,LBj),                                 &
# endif
     &                   Awrk(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          status=nf90_inquire_variable(OUTncid, OUTvid, ncvname)
          WRITE (stdout,20) TRIM(ncvname), Irec, TRIM(OUTncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Process intermediate adjoint time records.
!-----------------------------------------------------------------------
!
      DO Irec=2,Nrec-1
        Isav=Iold
        Iold=Inow
        Inow=Inxt
        Inxt=Isav
!
!  Read in next time record.
!
        status=nf_fread2d(ng, model, INPncid, INPvid, Irec+1, gtype,    &
     &                    Vsize, LBi, UBi, LBj, UBj,                    &
     &                    scale, Fmin, Fmax,                            &
# ifdef MASKING
     &                    mask(LBi,LBj),                                &
# endif
     &                    ad_A(LBi,LBj,Inxt))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            status=nf90_inquire_variable(INPncid, INPvid, ncvname)
            WRITE (stdout,10) TRIM(ncvname), Irec, TRIM(INPncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        END IF
!
!  Compute impulse forcing.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            Awrk(i,j)=C2*ad_A(i,j,Iold)+                                &
     &                C3*ad_A(i,j,Inow)+                                &
     &                C2*ad_A(i,j,Inxt)
          END DO
        END DO
!
!  Write out impulse forcing record in backward order so it is in
!  ascending time order.
!
        Iout=Iout-1
        status=nf_fwrite2d(ng, model, OUTncid, OUTvid, Iout, gtype,     &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     mask(LBi,LBj),                               &
# endif
     &                     Awrk(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            status=nf90_inquire_variable(OUTncid, OUTvid, ncvname)
            WRITE (stdout,20) TRIM(ncvname), Irec, TRIM(OUTncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Process last adjoint time record.
!-----------------------------------------------------------------------
!
      IF (Nrec.gt.1) THEN
        Irec=Nrec
        Iold=Inow
        Inow=Inxt
!
!  Compute impulse forcing.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            Awrk(i,j)=C2*ad_A(i,j,Iold)+                                &
     &                C1*ad_A(i,j,Inow)

          END DO
        END DO
!
!  Write out last impulse forcing into first record.
!
        Iout=1
        status=nf_fwrite2d(ng, model, OUTncid, OUTvid, Iout, gtype,     &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     mask(LBi,LBj),                               &
# endif
     &                     Awrk(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            status=nf90_inquire_variable(OUTncid, OUTvid, ncvname)
              WRITE (stdout,20) TRIM(ncvname), Irec, TRIM(OUTncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

  10  FORMAT (/,' IMPULSE2D - error while reading variable: ',a,2x,     &
     &        'at time record = ',i3,/,17x,'in input NetCDF file: ',a)
  20  FORMAT (/,' IMPULSE2D - error while writing variable: ',a,2x,     &
     &        'at time record = ',i3,/,17x,'into NetCDF file: ',a)

      RETURN
      END SUBROUTINE impulse2d

# ifdef SOLVE3D
      SUBROUTINE impulse3d (ng, model, gtype, outer_impulse,            &
     &                      Nrec, Nsave,                                &
     &                      INPncname, INPncid, INPvid,                 &
     &                      OUTncname, OUTncid, OUTvid,                 &
     &                      Imin, Imax, Jmin, Jmax,                     &
     &                      LBi, UBi, LBj, UBj, LBk, UBk,               &
#  ifdef MASKING
     &                      mask,                                       &
#  endif
     &                      ad_A)
!
!=======================================================================
!                                                                      !
!  This routine reads in requested 3D adjoint variable from input      !
!  NetCDF file, converts to impulse,  and then writes into output      !
!  NetCDF file.  The impulses are written in time ascending order      !
!  since they will be processed by forward TL and RP models.           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_netcdf
      USE mod_scalars
!
!  Imported variable declarations.
!
      logical, intent(in) :: outer_impulse

      integer, intent(in) :: ng, model, gtype, Nrec, Nsave
      integer, intent(in) :: INPncid, INPvid, OUTncid, OUTvid
      integer, intent(in) :: Imin, Imax, Jmin, Jmax
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk

      character(len=*), intent(in) :: INPncname
      character(len=*), intent(in) :: OUTncname
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: mask(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: mask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk,3)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, k, status
      integer :: Irec, Iold, Inow, Inxt, Iout, Isav

      integer :: nf_fread3d
      integer :: nf_fwrite3d

      integer, dimension(4) :: Vsize(4)

      real(r8) :: C1, C2, C3, Fmin, Fmax, cff, scale

      real(r8), dimension(LBi:UBi,LBj:UBj,LBk:UBk) :: Awrk

      character(len=40) :: ncvname
!
!-----------------------------------------------------------------------
!  Initialize time level indices and impulse focing coefficients.
!-----------------------------------------------------------------------
!
      Inow=1
      Inxt=2
      Iold=3
      cff=REAL(Nsave,r8)
      IF (outer_impulse) THEN
        C1=(cff+1.0_r8)*(2.0_r8*cff+1.0_r8)/(6.0_r8*cff)
        C2=(cff*cff-1.0_r8)/(6.0_r8*cff)
        C3=(2.0_r8*cff*cff+1.0_r8)/(3.0_r8*cff)
      ELSE
        C1=1.0_r8
        C2=0.0_r8
        C3=1.0_r8
      END IF
      scale=1.0_r8
!
!  Set Vsize to zero to deactivate interpolation of input data to model
!  grid in "nf_fread3d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Initialize work array.
!
      Awrk=0.0_r8
!
!-----------------------------------------------------------------------
!  Process first adjoint time record.
!-----------------------------------------------------------------------
!
      Irec=1
!
!  Read initial time record.
!
      status=nf_fread3d(ng, model, INPncid, INPvid, Irec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj, LBk, UBk,            &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  mask(LBi,LBj),                                  &
#  endif
     &                  ad_A(LBi,LBj,LBk,Inow))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          status=nf90_inquire_variable(INPncid, INPvid, ncvname)
          WRITE (stdout,10) TRIM(ncvname), Irec, TRIM(INPncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF
!
!  Read next time record.
!
      IF (Nrec.gt.1) THEN
        status=nf_fread3d(ng, model, INPncid, INPvid, Irec+1, gtype,    &
     &                    Vsize, LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                    scale, Fmin, Fmax,                            &
#  ifdef MASKING
     &                    mask(LBi,LBj),                                &
#  endif
     &                    ad_A(LBi,LBj,LBk,Inxt))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            status=nf90_inquire_variable(INPncid, INPvid, ncvname)
            WRITE (stdout,10) TRIM(ncvname), Irec+1, TRIM(INPncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      ELSE
        C2=0.0_r8
      END IF
!
!  Compute impulse forcing.
!
      DO k=LBk,UBk
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            Awrk(i,j,k)=C1*ad_A(i,j,k,Inow)+                            &
     &                  C2*ad_A(i,j,k,Inxt)
          END DO
        END DO
      END DO
!
!  Write out initial impulse forcing record to last record.
!
      Iout=Nrec
      status=nf_fwrite3d(ng, model, OUTncid, OUTvid, Iout, gtype,       &
     &                   LBi, UBi, LBj, UBj, LBk, UBk, scale,           &
#  ifdef MASKING
     &                   mask(LBi,LBj),                                 &
#  endif
     &                   Awrk(LBi,LBj,LBk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          status=nf90_inquire_variable(OUTncid, OUTvid, ncvname)
          WRITE (stdout,20) TRIM(ncvname), Irec, TRIM(OUTncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Process intermediate adjoint time records.
!-----------------------------------------------------------------------
!
      DO Irec=2,Nrec-1
        Isav=Iold
        Iold=Inow
        Inow=Inxt
        Inxt=Isav
!
!  Read in next time record.
!
        status=nf_fread3d(ng, model, INPncid, INPvid, Irec+1, gtype,    &
     &                    Vsize, LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                    scale, Fmin, Fmax,                            &
#  ifdef MASKING
     &                    mask(LBi,LBj),                                &
#  endif
     &                    ad_A(LBi,LBj,LBk,Inxt))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            status=nf90_inquire_variable(INPncid, INPvid, ncvname)
            WRITE (stdout,10) TRIM(ncvname), Irec+1, TRIM(INPncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        END IF
!
!  Compute impulse forcing.
!
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              Awrk(i,j,k)=C2*ad_A(i,j,k,Iold)+                          &
     &                    C3*ad_A(i,j,k,Inow)+                          &
     &                    C2*ad_A(i,j,k,Inxt)
            END DO
          END DO
        END DO
!
!  Write out impulse forcing record in backward order so it is in
!  ascending time order.
!
        Iout=Iout-1
        status=nf_fwrite3d(ng, model, OUTncid, OUTvid, Iout, gtype,     &
     &                     LBi, UBi, LBj, UBj, LBk, UBk, scale,         &
#  ifdef MASKING
     &                     mask(LBi,LBj),                               &
#  endif
     &                     Awrk(LBi,LBj,LBk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            status=nf90_inquire_variable(OUTncid, OUTvid, ncvname)
            WRITE (stdout,20) TRIM(ncvname), Irec, TRIM(OUTncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Process last adjoint time record.
!-----------------------------------------------------------------------
!
      IF (Nrec.gt.1) THEN
        Irec=Nrec
        Iold=Inow
        Inow=Inxt
!
!  Compute impulse forcing.
!
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              Awrk(i,j,k)=C2*ad_A(i,j,k,Iold)+                          &
     &                    C1*ad_A(i,j,k,Inow)
            END DO
          END DO
        END DO
!
!  Write out last impulse forcing into first record
!
        Iout=1
        status=nf_fwrite3d(ng, model, OUTncid, OUTvid, 1, gtype,        &
     &                     LBi, UBi, LBj, UBj, LBk, UBk, scale,         &
#  ifdef MASKING
     &                     mask(LBi,LBj),                               &
#  endif
     &                     Awrk(LBi,LBj,LBk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            status=nf90_inquire_variable(OUTncid, OUTvid, ncvname)
            WRITE (stdout,20) TRIM(ncvname), Irec, TRIM(OUTncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

  10  FORMAT (/,' IMPULSE3D - error while reading variable: ',a,2x,     &
     &        'at time record = ',i3,/,17x,'in input NetCDF file: ',a)
  20  FORMAT (/,' IMPULSE3D - error while writing variable: ',a,2x,     &
     &        'at time record = ',i3,/,17x,'into NetCDF file: ',a)

      RETURN
      END SUBROUTINE impulse3d
# endif
#endif
      END MODULE impulse_mod
