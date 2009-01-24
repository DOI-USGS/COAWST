#include "cppdefs.h"
#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
      SUBROUTINE stats_modobs (ng)
!
!svn $Id: stats_modobs.F 645 2008-07-16 19:26:48Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine computes several statistical quantities between     !
!  model and observations:                                             !
!                                                                      !
!     CC         Cross-Correlation                                     !
!     MB         Model Bias                                            !
!     MSE        Mean Squared Error                                    !
!     SDE        Standard Deviation Error                              !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Oke, P.R., J.S. Allen, R.N. Miller, G.D. Egbert, J.A. Austin,       !
!    J.A. Barth, T.J. Boyd, P.M. Kosro, and M.D. Levine, 2002: A       !
!    Modeling Study of the Three-Dimensional Continental Shelf         !
!    Circulation off Oregon. Part I: Model-Data Comparison, J.         !
!    Phys. Oceanogr., 32, 1360-1382.                                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: i, ic, iobs, status, varid
      integer, dimension(2) :: start, total
      integer, dimension(0:NstateVar(ng)) :: Ncount, is

      integer, allocatable :: Obs_type(:)

      real(r8) :: cff1, cff2
      real(r8), parameter :: LARGE = 1.0e+10_r8

      real(r8), dimension(0:NstateVar(ng)) :: CC, MB, MSE, SDE
      real(r8), dimension(0:NstateVar(ng)) :: mod_min, mod_max
      real(r8), dimension(0:NstateVar(ng)) :: mod_mean, mod_std
      real(r8), dimension(0:NstateVar(ng)) :: obs_min, obs_max
      real(r8), dimension(0:NstateVar(ng)) :: obs_mean, obs_std

      real(r8), allocatable :: mod_value(:)
      real(r8), allocatable :: Obs_scale(:)
      real(r8), allocatable :: Obs_value(:)

      character (len=11), dimension(NstateVar(ng)) :: text, svar_name
!
!-----------------------------------------------------------------------
!  Read in model and observations data.
!-----------------------------------------------------------------------
!
      IF (InpThread) THEN
        start(1)=1
        total(1)=Ndatum(ng)
!
!  Allocate working arrays.
!
        IF (.not.allocated(mod_value)) allocate (mod_value(Ndatum(ng)))
        IF (.not.allocated(obs_scale)) allocate (obs_scale(Ndatum(ng)))
        IF (.not.allocated(obs_type))  allocate (obs_type (Ndatum(ng)))
        IF (.not.allocated(obs_value)) allocate (obs_value(Ndatum(ng)))
!
!  Read in observation type identifier.
!
        status=nf90_get_var(ncOBSid(ng), obsVid(idOtyp,ng),             &
     &                      obs_type, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idOtyp)), TRIM(OBSname(ng))
          exit_flag=2
          ioerror=status
          RETURN
        END IF
!
!  Read in observation values.
!
        status=nf90_get_var(ncOBSid(ng), obsVid(idOval,ng),             &
     &                      obs_value, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idOval)), TRIM(OBSname(ng))
          exit_flag=2
          ioerror=status
          RETURN
        END IF
!
!  Read in observation screening flag.
!
        status=nf90_get_var(ncMODid(ng), modVid(idObsS,ng),             &
     &                      obs_scale, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idObsS)), TRIM(MODname(ng))
          exit_flag=2
          ioerror=status
          RETURN
        END IF
!
!  Read in model values at observation locations.
!
# ifdef S4DVAR
        start(2)=2
        total(2)=1
# endif
# if defined W4DVAR || defined IOMS
        status=nf90_get_var(ncMODid(ng), modVid(idTLmo,ng),             &
     &                      mod_value, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTLmo)), TRIM(MODname(ng))
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# else
        status=nf90_get_var(ncMODid(ng), modVid(idNLmo,ng),             &
     &                      mod_value, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idNLmo)), TRIM(MODname(ng))
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Compute model and observations comparison statistics.
!-----------------------------------------------------------------------
!
!  Initialize.
!
        DO i=0,NstateVar(ng)
          CC(i)=0.0_r8
          MB(i)=0.0_r8
          MSE(i)=0.0_r8
          SDE(i)=0.0_r8
          mod_min(i)=LARGE
          mod_max(i)=-LARGE
          mod_mean(i)=0.0_r8
          obs_min(i)=LARGE
          obs_max(i)=-LARGE
          obs_mean(i)=0.0_r8
          mod_std(i)=0.0_r8
          obs_std(i)=0.0_r8
          Ncount(i)=0
        END DO
!
!  Compute model and observations mean per each state variable.
!
        DO iobs=1,Ndatum(ng)
          IF (obs_scale(iobs).gt.0.0_r8) THEN
            i=obs_type(iobs)
            Ncount(i)=Ncount(i)+1
            mod_min(i)=MIN(mod_min(i),mod_value(iobs))
            obs_min(i)=MIN(obs_min(i),obs_value(iobs))
            mod_max(i)=MAX(mod_max(i),mod_value(iobs))
            obs_max(i)=MAX(obs_max(i),obs_value(iobs))
            mod_mean(i)=mod_mean(i)+mod_value(iobs)
            obs_mean(i)=obs_mean(i)+obs_value(iobs)
          END IF
        END DO
        DO i=1,NstateVar(ng)
          IF (Ncount(i).gt.0) THEN
            mod_mean(i)=mod_mean(i)/REAL(Ncount(i),r8)
            obs_mean(i)=obs_mean(i)/REAL(Ncount(i),r8)
          END IF
        END DO
!
!  Compute standard deviation and cross-correlation between model and
!  observations (CC).
!
        DO iobs=1,Ndatum(ng)
          IF (obs_scale(iobs).gt.0.0_r8) THEN
            i=obs_type(iobs)
            cff1=mod_value(iobs)-mod_mean(i)
            cff2=obs_value(iobs)-obs_mean(i)
            mod_std(i)=mod_std(i)+cff1*cff1
            obs_std(i)=obs_std(i)+cff2*cff2
            CC(i)=CC(i)+cff1*cff2
          END IF
        END DO
        DO i=1,NstateVar(ng)
          IF (Ncount(i).gt.1) THEN
            mod_std(i)=SQRT(mod_std(i)/REAL(Ncount(i)-1,r8))
            obs_std(i)=SQRT(obs_std(i)/REAL(Ncount(i)-1,r8))
            CC(i)=(CC(i)/REAL(Ncount(i),r8))/(mod_std(i)*obs_std(i))
          END IF
        END DO
!
!  Compute model bias (MB), standard deviation error (SDE), and mean
!  squared error (MSE).
!
        DO i=1,NstateVar(ng)
          IF (Ncount(i).gt.0) THEN
            MB(i)=mod_mean(i)-obs_mean(i)
            SDE(i)=mod_std(i)-obs_std(i)
            MSE(i)=MB(i)*MB(i)+                                         &
     &             SDE(i)*SDE(i)+                                       &
     &             2.0_r8*mod_std(i)*obs_std(i)*(1.0_r8-CC(i))            
          END IF
        END DO
!
!  Write comparison statistics to NetCDF file.
!
        start(1)=1
        total(1)=NstateVar(ng)+1

        status=nf90_inq_varid(ncMODid(ng), 'Nobs', varid)
        status=nf90_put_var(ncMODid(ng), varid, Ncount(0:),             &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'Nobs', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf90_inq_varid(ncMODid(ng), 'obs_mean', varid)
        status=nf90_put_var(ncMODid(ng), varid, obs_mean(0:),           &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'obs_mean', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf90_inq_varid(ncMODid(ng), 'obs_std', varid)
        status=nf90_put_var(ncMODid(ng), varid, obs_std(0:),            &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'obs_std', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf90_inq_varid(ncMODid(ng), 'model_mean', varid)
        status=nf90_put_var(ncMODid(ng), varid, mod_mean(0:),           &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'model_mean', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf90_inq_varid(ncMODid(ng), 'model_std', varid)
        status=nf90_put_var(ncMODid(ng), varid, mod_std(0:),            &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'model_std', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf90_inq_varid(ncMODid(ng), 'model_bias', varid)
        status=nf90_put_var(ncMODid(ng), varid, MB(0:),                 &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'model_bias', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf90_inq_varid(ncMODid(ng), 'SDE', varid)
        status=nf90_put_var(ncMODid(ng), varid, SDE(0:),                &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'SDE', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf90_inq_varid(ncMODid(ng), 'CC', varid)
        status=nf90_put_var(ncMODid(ng), varid, CC(0:),                 &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'CC', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf90_inq_varid(ncMODid(ng), 'MSE', varid)
        status=nf90_put_var(ncMODid(ng), varid, MSE(0:),                &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'MSE', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Synchronize NetCDF file to disk.
!
        status=nf90_sync(ncMODid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Report model and observations comparison statistics.
!
        ic=0
        DO i=1,NstateVar(ng)
          svar_name(i)='           '
          text(i)='           '
          IF (Ncount(i).gt.0) THEN
            ic=ic+1
            is(ic)=i
            svar_name(ic)=TRIM(Vname(1,idSvar(i)))
            text(ic)='-----------'
          END IF
        END DO
        WRITE (stdout,40)
        WRITE (stdout,50) (svar_name(i),i=1,ic)
        WRITE (stdout,60) (text(i),i=1,ic)
        WRITE (stdout,70) 'Observation Min   ', (obs_min (is(i)),i=1,ic)
        WRITE (stdout,70) 'Observation Max   ', (obs_max (is(i)),i=1,ic)
        WRITE (stdout,70) 'Observation Mean  ', (obs_mean(is(i)),i=1,ic)
        WRITE (stdout,70) 'Observation STD   ', (obs_std (is(i)),i=1,ic)
        WRITE (stdout,70) 'Model Min         ', (mod_min (is(i)),i=1,ic)
        WRITE (stdout,70) 'Model Max         ', (mod_max (is(i)),i=1,ic)
        WRITE (stdout,70) 'Model Mean        ', (mod_mean(is(i)),i=1,ic)
        WRITE (stdout,70) 'Model STD         ', (mod_std (is(i)),i=1,ic)
        WRITE (stdout,70) 'Model Bias        ', (MB(is(i)),i=1,ic)
        WRITE (stdout,70) 'STD Error         ', (SDE(is(i)),i=1,ic)
        WRITE (stdout,70) 'Cross-Correlation ', (CC(is(i)),i=1,ic)
        WRITE (stdout,70) 'Mean Squared Error', (MSE(is(i)),i=1,ic)
        WRITE (stdout,80) 'Observation Count ', (Ncount(is(i)),i=1,ic)
!
!  Allocate working arrays.
!
        IF (allocated(mod_value)) deallocate (mod_value)
        IF (allocated(obs_scale)) deallocate (obs_scale)
        IF (allocated(obs_type))  deallocate (obs_type)
        IF (allocated(obs_value)) deallocate (obs_value)
      END IF
!
 10   FORMAT (/,' STATS_4DVAR - error while reading variable: ',a,/     &
     &        15x,'in NetCDF file: ',a)
 20   FORMAT (/,' STATS_4DVAR - error while writing variable: ',a,/,    &
     &        12x,'into NetCDF file: ',a)
 30   FORMAT (/,' STATS_4DVAR - unable to synchronize 4DVAR',           &
     &        1x,'NetCDF file to disk.')
 40   FORMAT (/,' 4DVAR Model-Observations Comparison Statistics:',/)
 50   FORMAT (t22,5(a11,1x))
 60   FORMAT (t22,5(a11,1x),/)
 70   FORMAT (a,3x,5(1p,e11.4,0p,1x))
 80   FORMAT (a,3x,5(i11,1x))
#else
      SUBROUTINE stats_modobs
#endif
      RETURN
      END SUBROUTINE stats_modobs
