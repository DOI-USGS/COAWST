#include "cppdefs.h"

!  1. Select double or single. The double sum was taking too long,
!     stick with single sum.
#undef  DOUBLESUM
#define SINGLESUM
!
#ifdef SINGLESUM
#  define INTCDF
#  undef  FCDF
#endif
!  2. Select direct or integrated hilbert.
#define DHILBERT
#undef  IHILBERT
!  3. this hscale may be needed for ihilbert.
#undef HSCALE

! SINGLESUM: SINGLE SUMMATION METHOD
! DOUBLESUM: DOUBLE SUMMATION METHOD
! INTCDF:    INTEGRATED DIRECTIONAL DISTRIBUTION
! FCDF:      FREQUENCY-DEPENDENT DIRECTIONAL DISTRIBUTIONS
! DHILBERT:  HILBERT TAKEN OVER DIRECTIONAL AMPZETA
! IHILBERT:  HILBERT TAKEN OVER INTEGRATED AMPZETA
! HSCALE:    SCALE ENERGY BY ENERGY IN SWAN SPECTRUM

      MODULE mod_inwave_swan
!
!svn $Id: swan_reader.F 1336 2008-01-24 02:45:56Z jcwarner $
!=======================================================================
!                                                                      !
!  NDIR        number of directional dimensions in Fourier components  !
!  NFREQ       number of frequency dimensions in Fourier components    !
!  NDIRSW      number of directional dimensions in SWAN spectra        !
!  NFREQSW     number of frequency dimensions in SWAN spectra          !
!  Insteps     number of time steps                                    !
!  DUR         duration of the time series to generate                 !
!  DT          time step in the time series of the free surface        !
!              elevation signal                                        !
!  FMAX        Maximum frequency to consider when recontructing the    !
!              free surface signal - this is equivalent to the Nyquist !
!              frequency                                               !
!  FMIN        Minimum frequency to consider when recontructing the    !
!              free surface signal                                     !
!  DF          Frequency interval for the frequency array              !
!  TREP        Representative period                                   !
!  DDIR        Directional interval for the directions array           !
!  FP          Peak frequency                                          !
!  MAINANG     Short waves main propagation angle                      !
!  MAINANGI    Short waves main propagation angle index                !
!  HM0GEW      Significant wave height computed from the               !
!              interpolated spectra                                    !
!  DANG        
!  INT         
!  DF_FFT      Frequency increment for the positive Fourier Components !
!  FSW         Frequency array                                         !
!  F           Frequency array                                         !
!  DIR         Directional array read from swan                        !
!  THETA       Directional array                                       !
!  SF          Spectral density function read from swan                !
!  SDSW        Directional spectral density function read from swan    !
!  SD          Directional spectral density function                   !
!  SDD         Spectral density function integrated over frequency     !
!  SDF         Spectral density function integrated over direction     !
!  PHASE       Random phase for each frequency- direction component    !
!  AMP         Amplitude of the Fourier components                     !
!  AMPZETA     Amplitude of the free surface elevation for each        !
!              directional component                                   !
!  AMPZETA_TOT Amplitude of the free surface elevation                 !
!  POSITION    positive frequencies of the Fourier Components          !
!  DIR_3_INT   index of the directional array for each freq            !
!=======================================================================
!
!
!=======================================================================
!                                                                      !
!  CONTAINS THE FOLLOWING SUBROUTINES:                                 !
!                                                                      !
!     inwave_swan  : main driver of boundary condition computations    !
!                    from swan output spectra                          !
!     swan_reader  : read swan output spectra                          !
!     array_gen    : generate the arrays necessary to compute the free !
!                    surface elevation time series from the            !
!                    directional wave spectra derived from swan        !
!     random_phase : assign random phases to each frequency-direction  !
!                    component                                         !
!     tpdcalc    : This subroutine computes the respresentative period !
!     amplitudes   : compute the amplitude for each frequency-direction!
!                    components                                        !
!     FourierComp  : compute the Fourier components                    !
!=======================================================================
!

#ifdef INWAVE_SWAN_COUPLING

        USE mod_kinds
        USE mod_inwave_params
        USE math_tools
        USE mod_iounits

        implicit none

        TYPE SHORT_WAVE

         integer             :: ndir, nfreq, mainangi
         integer             :: ndirsw, nfreqsw
         integer             :: Insteps, Swancircle
         integer             :: SstrR, SendR

         real(r8)            :: dur, idt
         real(r8)            :: fmax, fmin, df
         real(r8)            :: ddir, fp, exc
         real(r8)            :: mainang, hm0gew, dang, integ
         real(r8)            :: df_fft, ampmax
         real(r8)            :: SpecTimeIso(2)
         real(r8)            :: SpecTime(2)
         real(r8)            :: SpecTimedt
         real(r8)            :: h0
         real(r8)            :: dtifac

         integer, pointer  :: posang(:)
         integer, pointer  :: dir_3_int(:)

         real(r8), pointer :: fsw(:)
         real(r8), pointer :: f(:)
         real(r8), pointer :: dir(:)
         real(r8), pointer :: theta(:)
         real(r8), pointer :: SF(:)

         real(r8), pointer :: SDSW(:,:)
         real(r8), pointer :: SD(:,:)
         real(r8), pointer :: SDD(:)
         real(r8), pointer :: SDF(:)
         real(r8), pointer :: phase(:,:)
# if defined SINGLESUM
         real(r8), pointer :: amp(:)
# elif defined DOUBLESUM
         real(r8), pointer :: amp(:,:)
# endif
         real(r8), pointer :: zeta(:,:)
         real(r8), pointer :: Ampzeta(:,:,:)
         complex(fftkind), dimension(:), allocatable :: Ampzeta_tot
         real(r8), pointer :: bndwave_z(:,:)
         real(r8), pointer :: bndwave_u(:,:)
         real(r8), pointer :: bndwave_v(:,:)
         real(r8), pointer :: position(:)
         real(r8), pointer :: dist_x(:)
         real(r8), pointer :: dist_y(:)
         real(r8), pointer :: dist_xy(:)
         real(r8), pointer :: angle_edge(:)
         real(r8), pointer :: waven(:)

        END TYPE SHORT_WAVE

        TYPE (SHORT_WAVE), allocatable :: WAVES(:)

      CONTAINS
!
!***********************************************************************
      SUBROUTINE allocate_inwave_swan (ng)
!***********************************************************************
!
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_iounits
!
!  Local variable declarations.
!
      integer, intent(in) :: ng

!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( WAVES(Ngrids) )

      RETURN
      END SUBROUTINE allocate_inwave_swan

!
!***********************************************************************
      SUBROUTINE inwave_swan_run (ng, first, tile)
!***********************************************************************
!
!=======================================================================
!                                                                      !
!  Computes the free surface elevation time series from the swan       !
!  directional spectra                                                 !
!                                                                      !
!  On Input:                                                           !
!     ng         Number of grids                                       !
!                                                                      !
!  On Output:                                                          !
!     Ampzeta    Free surface elevation time series for each           !
!                directional bin                                       !
!                                                                      !
!=======================================================================
!

      USE mod_iounits
      USE mod_scalars
      USE mod_inwave_params

!  Imported variable declarations.
!
      integer, intent(in) :: ng, first, tile

!-----------------------------------------------------------------------
! Read swan output spectra
!-----------------------------------------------------------------------
!
      call swan_reader (ng, first, tile)
!
!-----------------------------------------------------------------------
! Compute the boundary water depth and grid distances
!-----------------------------------------------------------------------
!
      call grid_metrics (ng, tile, first)
!
!-----------------------------------------------------------------------
! Generate the arrays 
!-----------------------------------------------------------------------
!
      call array_gen (ng, first)
!
!-----------------------------------------------------------------------
! Generate random phases for each frequency-direction component
!-----------------------------------------------------------------------
!
      call random_phase (ng)
!
!-----------------------------------------------------------------------
! Compute the amplitudes for each frequency-direction component
!-----------------------------------------------------------------------
!
      call amplitudes (ng)
!
!-----------------------------------------------------------------------
! Compute the wave numbers
!-----------------------------------------------------------------------
!
      call compute_k (ng, tile, first)
!
!-----------------------------------------------------------------------
! Compute the Fourier Components
!-----------------------------------------------------------------------
!
      call FourierComp (ng, tile, first)
!
!-----------------------------------------------------------------------
! Compute the boundwave
!-----------------------------------------------------------------------
!
      call boundwave (ng, tile, first)
!
!-----------------------------------------------------------------------
!
      RETURN

      END SUBROUTINE inwave_swan_run

!
!***********************************************************************
      SUBROUTINE array_gen (ng, first)
!***********************************************************************
!
!=======================================================================
!                                                                      !
!  Generates the arrays necesary to compute the free surface           !
!  elevation time series from the swan directional spectra             !
!                                                                      !
!=======================================================================
!
      USE mod_inwave_bound
      USE mod_inwave_params
      USE mod_inwave_vars
      USE mod_parallel
      USE mod_scalars
      USE interpolate_mod

      implicit none

!  Imported variable declarations.
!
      integer, intent(in) :: ng, first

!  Local variable declarations.

      logical             :: rectangular
      integer             :: dd, ff, k, offset, i, tidx
      real(r8), parameter :: IJspv = 0.0_r8
      real(r8)            :: cff
      real(r8)            :: my_min, my_max
# ifdef INWAVE_IGBCSOURCE
      real(r8)            :: cff1, cff2, cff3, cff4
      real(r8)            :: twopi, otwopi
      real(r8)            :: kh, cwc, G1, cgr
      real(r8), allocatable :: k(:)
      real(r8), allocatable :: kwc(:)
# endif
      real(r8), allocatable :: angle(:,:), Iout(:,:), Jout(:,:)
      real(r8), allocatable :: fsw_2d(:,:), theta_2d(:,:)
      real(r8), allocatable :: f_2d(:,:), wd_2d(:,:), SDSW_circle(:,:)

# ifdef INWAVE_IGBCSOURCE
      real(r8), parameter :: khmax=100.0_r8
      real(r8), parameter :: eps = 1.0E-10_r8
      real(r8), parameter :: maxErr = 0.1_r8
      real(r8), parameter :: alpha1 = 0.0006_r8
      twopi=2.0_r8*pi
      otwopi=1.0_r8/(2.0_r8*pi)
# endif

!
!-----------------------------------------------------------------------
!  Define the Nyquist frequency and minimum frequency and 
!  the number of directional bins
!  Determine the time dimensions of duration and num steps.
!  The Ampzeta time series will repeat every 1/df time steps.
!  For old way:
!  If you keep fmax=0.5 then the number of freqs will be 1/2 the time length.
!  We want this for the bound wave.
!  The Ampzeta time series will repeat every 1/df time steps.
!-----------------------------------------------------------------------
!
      IF (first.eq.1) THEN
# ifdef SINGLESUM
! new way, use a 3600s time array at dt=1 sec
        WAVES(ng)%ndir=ND  !use this to force AC/bnd computation grid be user size
!       WAVES(ng)%ndir=36  !use this to force AC/bnd computation grid be 10deg.
        WAVES(ng)%idt=1.0_r8
        WAVES(ng)%dur=3600.0_r8
        WAVES(ng)%Insteps=nint(WAVES(ng)%dur/WAVES(ng)%idt)
        WAVES(ng)%df=1.0_r8/(WAVES(ng)%Insteps*WAVES(ng)%idt)
        WAVES(ng)%fmin=WAVES(ng)%df
        WAVES(ng)%nfreq=WAVES(ng)%Insteps/2                     !-1
        WAVES(ng)%fmax=REAL(WAVES(ng)%nfreq,r8)*WAVES(ng)%df
        WAVES(ng)%dtifac=1.0_r8/5.0_r8
# else
! old way
        WAVES(ng)%nfreq=200
        WAVES(ng)%fmax=0.5_r8  !! Tmin = 2 s, this is ~ the Nyquist f
        WAVES(ng)%df=WAVES(ng)%fmax/WAVES(ng)%nfreq
        WAVES(ng)%ndir=ND
        WAVES(ng)%fmin=WAVES(ng)%df
        WAVES(ng)%idt=1.0_r8
        WAVES(ng)%dur=1./WAVES(ng)%df
        WAVES(ng)%Insteps=nint(WAVES(ng)%dur/WAVES(ng)%idt)
# endif
!
        IF (Master) THEN
          WRITE(stdout,*) 'Computing InWave boundary forcing'
          WRITE(stdout,*) 'Freqs min max are : ',WAVES(ng)%fmin,        &
     &                                           WAVES(ng)%fmax
        END IF
!
!-----------------------------------------------------------------------
! Create the frequency and directional arrays for the fft
!-----------------------------------------------------------------------
!
        allocate (WAVES(ng)%f(WAVES(ng)%nfreq))
        allocate (WAVES(ng)%dir(WAVES(ng)%ndir))
        allocate (WAVES(ng)%SD(WAVES(ng)%nfreq,WAVES(ng)%ndir))
        allocate (WAVES(ng)%SDD(WAVES(ng)%ndir))
        allocate (WAVES(ng)%SDF(WAVES(ng)%nfreq))
        DO dd=1,WAVES(ng)%ndir
          DO ff=1,WAVES(ng)%nfreq
            WAVES(ng)%SD(ff,dd)=0.0_r8
          ENDDO
          WAVES(ng)%SDD(dd)=0.0_r8
          WAVES(ng)%dir(dd)=WAVEG(ng)%wd(dd)  !use this to set dirs to be ND
!                                             !use these next lines for
!                                             !internal computation grid.
!         WAVES(ng)%dir(dd)=REAL(dd-1,r8)*                              &
!    &                       360.0_r8/REAL(WAVES(ng)ndir,r8)*pi/180.0_r8
        ENDDO
        DO ff=1,WAVES(ng)%nfreq
          WAVES(ng)%SDF(ff)=0.0_r8
        ENDDO
        DO ff=1,WAVES(ng)%nfreq
          WAVES(ng)%f(ff)=REAL(ff-1,r8)*WAVES(ng)%df+WAVES(ng)%fmin
        END DO
      END IF
!
!-----------------------------------------------------------------------
! Interpolate from the SWAN 2D spectral grid to the 2D spectral grid that
! we predefined in the ini file.
!
!  Set up 2d gridded freq and dir arrays using the SWAN data.
!  Also here if SWAN was computed on a full circle, then 
!  we mirror the SWAN data from -360 to +720. this allows
!  user to define a smaller computational grid if needed.
!
      rectangular=.TRUE.
      IF (WAVES(ng)%Swancircle.eq.1) THEN
        offset=3
      ELSE
        offset=1
      ENDIF
      allocate (angle(1:WAVES(ng)%nfreqsw,1:WAVES(ng)%ndirsw*offset))
      allocate (fsw_2d(1:WAVES(ng)%nfreqsw,1:WAVES(ng)%ndirsw*offset))
      allocate (theta_2d(1:WAVES(ng)%nfreqsw,1:WAVES(ng)%ndirsw*offset))
      allocate (SDSW_circle(1:WAVES(ng)%nfreqsw,                        &
     &                      1:WAVES(ng)%ndirsw*offset))
      DO ff=1,WAVES(ng)%nfreqsw
        DO dd=1,WAVES(ng)%ndirsw*offset
          angle(ff,dd)=0.0_r8
          fsw_2d(ff,dd)=0.0_r8
          theta_2d(ff,dd)=0.0_r8
          SDSW_circle(ff,dd)=0.0_r8
        END DO
      END DO
      IF (WAVES(ng)%Swancircle.eq.1) THEN
        DO ff=1,WAVES(ng)%nfreqsw
          DO dd=1,WAVES(ng)%ndirsw*offset
            fsw_2d(ff,dd)=WAVES(ng)%fsw(ff)
            IF (dd.le.WAVES(ng)%ndirsw) THEN
              k=dd
              theta_2d(ff,dd)=WAVES(ng)%theta(k)-360.0_r8*pi/180.0_r8
            ELSEIF (dd.le.WAVES(ng)%ndirsw*2) THEN
              k=dd-WAVES(ng)%ndirsw
              theta_2d(ff,dd)=WAVES(ng)%theta(k)
            ELSE
              k=dd-WAVES(ng)%ndirsw*2
              theta_2d(ff,dd)=WAVES(ng)%theta(k)+360.0_r8*pi/180.0_r8
            END IF
            SDSW_circle(ff,dd)=WAVES(ng)%SDSW(ff,k)
          END DO
        END DO
      ELSE
        DO ff=1,WAVES(ng)%nfreqsw
          DO dd=1,WAVES(ng)%ndirsw*offset
            fsw_2d(ff,dd)=WAVES(ng)%fsw(ff)
            IF (dd.le.WAVES(ng)%ndirsw) THEN
              k=dd
              theta_2d(ff,dd)=WAVES(ng)%theta(k)
            END IF
            SDSW_circle(ff,dd)=WAVES(ng)%SDSW(ff,k)
          END DO
        END DO
      ENDIF
!
!  Set up 2d gridded freq and dir arrays to be used in the computations
!  of Ampzeta and Boundwave. SWAN 2dspec is interpolated to these arrays.
!
      allocate (Iout(1:WAVES(ng)%nfreq,1:WAVES(ng)%ndir))
      allocate (Jout(1:WAVES(ng)%nfreq,1:WAVES(ng)%ndir))
      allocate (f_2d(1:WAVES(ng)%nfreq,1:WAVES(ng)%ndir))
      allocate (wd_2d(1:WAVES(ng)%nfreq,1:WAVES(ng)%ndir))
      DO ff=1,WAVES(ng)%nfreq
        DO dd=1,WAVES(ng)%ndir
          Iout(ff,dd)=0.0_r8
          Jout(ff,dd)=0.0_r8
        END DO
      END DO
      DO ff=1,WAVES(ng)%nfreq
        DO dd=1,WAVES(ng)%ndir
          f_2d(ff,dd)=WAVES(ng)%f(ff)
          wd_2d(ff,dd)=WAVES(ng)%dir(dd)
        END DO
      END DO
!
      CALL hindices (ng, 1, WAVES(ng)%nfreqsw,                          &
     &                   1, WAVES(ng)%ndirsw*offset,                    &
     &                   1, WAVES(ng)%nfreqsw,                          &
     &                   1, WAVES(ng)%ndirsw*offset,                    &
     &                   angle, fsw_2d, theta_2d,                       &
     &                   1, WAVES(ng)%nfreq, 1, WAVES(ng)%ndir,         &
     &                   1, WAVES(ng)%nfreq, 1, WAVES(ng)%ndir,         &
     &                   f_2d, wd_2d,                                   &
     &                   Iout, Jout,                                    &
     &                   IJspv, rectangular)
      CALL linterp2d (ng, 1, WAVES(ng)%nfreqsw,                         &
     &                    1, WAVES(ng)%ndirsw*offset,                   &
     &                    fsw_2d, theta_2d, SDSW_circle,                &
     &                    1, WAVES(ng)%nfreq, 1, WAVES(ng)%ndir,        &
     &                    1, WAVES(ng)%nfreq, 1, WAVES(ng)%ndir,        &
     &                    Iout, Jout,                                   &
     &                    f_2d, wd_2d,                                  &
     &                    WAVES(ng)%SD,                                 &
                          my_min, my_max)
      deallocate(angle, Iout, Jout)
      deallocate(fsw_2d, theta_2d, SDSW_circle, f_2d, wd_2d)
!
!  Here the array WAVES(ng)%SD now has the spectrum from SWAN
!  interpolated to the computational grid.
!
# ifdef INWAVE_IGBCSOURCE
!  We can now add the IG boundary data to this energy array.
!  Ardhuin et al., 2014. Numerical model for free IG, OM 77, 20-32.
!    
      DO dd=1,WAVES(ng)%ndir
        DO ff=1,WAVES(ng)%nfreq
          IF (ff.gt.1.0_r8/400.0_r8) THEN
            kh=MIN(WAVES(ng)%waven(ff)*WAVES(ng)%h0,khmax)
            cwc=SQRT(g/(WAVES(ng)%waven(ff)+eps)*TANH(kh))
            G1=2.0_r8*kh/(SINH(2.0_r8*kh)+eps)
            cgr=0.5_r8*cwc*(1.0_r8+G1)
!
            cff1=1.0_r8/(cgr*twopi*WAVES(ng)%f(ff)*WAVES(ng)%df)
            cff2=(0.25*WAVES(ng)%hm0gew*WAVEG(ng)%Trep**2)**2
            cff3=MIN(1.0_r8,0.015_r8/WAVES(ng)%f(ff))
            cff4=1.2_r8*alpha1**2*WAVES(ng)%waven(ff)*g*g*              &
     &           cff1*cff2*cff3**1.5
            WAVES(ng)%SD(ff,dd)=WAVES(ng)%SD(ff,dd)+                    &
     &                          cff4*otwopi*WAVEG(ng)%pd
          ENDIF
        ENDDO
      ENDDO
# endif
!
!  Integrate energy over frequency
!
      DO dd=1,WAVES(ng)%ndir
        cff=0.0_r8
        DO ff=1,WAVES(ng)%nfreq
          cff=cff+WAVES(ng)%SD(ff,dd)
        ENDDO
        WAVES(ng)%SDD(dd)=cff
      ENDDO
!
!  Integrate frequency spectrum over direction 
!
      WAVES(ng)%integ=0.0_r8
      DO dd=1,WAVES(ng)%ndir
        WAVES(ng)%integ=WAVES(ng)%integ+WAVES(ng)%SDD(dd)
      ENDDO
      DO ff=1,WAVES(ng)%nfreq
        cff=0.0_r8
        DO dd=1,WAVES(ng)%ndir
          cff=cff+WAVES(ng)%SD(ff,dd)
        ENDDO
        WAVES(ng)%SDF(ff)=cff
      ENDDO
!
      IF (first.eq.1) THEN
!
!  Set directional arrays for model boundary to equal user input dirs.
!
        WAVEB(ng)%ND_bnd=ND
        allocate(WAVEB(ng)%WD_bnd(WAVEB(ng)%ND_bnd))
        DO dd=1,WAVEB(ng)%ND_bnd
          WAVEB(ng)%WD_bnd(dd)=WAVEG(ng)%wd(dd)
        END DO
!
!  Allocate and init the computational arrays.
!
        allocate(WAVES(ng)%Ampzeta(WAVES(ng)%SstrR:WAVES(ng)%SendR,     &
     &                             WAVES(ng)%Insteps,WAVES(ng)%ndir))
        allocate(WAVES(ng)%Ampzeta_tot(WAVES(ng)%Insteps))
!
!-----------------------------------------------------------------------
! Create frequency and directional arrays for the spectra
!-----------------------------------------------------------------------
!
        allocate (WAVES(ng)%position(WAVES(ng)%nfreq))
        allocate (WAVES(ng)%phase(WAVES(ng)%nfreq,WAVES(ng)%ndir))
# if defined SINGLESUM
        allocate (WAVES(ng)%amp(WAVES(ng)%nfreq))
        allocate (WAVES(ng)%dir_3_int(WAVES(ng)%nfreq))
# elif defined DOUBLESUM
        allocate (WAVES(ng)%amp(WAVES(ng)%nfreq,WAVES(ng)%ndir))
# endif
      END IF
!
      DO i=WAVES(ng)%SstrR,WAVES(ng)%SendR
        DO tidx=1,WAVES(ng)%Insteps
          DO dd=1,WAVES(ng)%ndir
            WAVES(ng)%Ampzeta(i,tidx,dd)=0.0_r8
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE array_gen

!
!***********************************************************************
      SUBROUTINE random_phase (ng)
!***********************************************************************
!
      USE mod_parallel
      USE mod_scalars

      implicit none

!  Imported variable declarations.
!
      integer, intent(in) :: ng
      
!  Local variable declarations.

      integer :: ff, dd, k, Npts, MyError
      integer :: i, clock
      real(r8) :: twopi
# ifdef DISTRIBUTE
      real(r8), allocatable :: wrk(:)
# endif
      call random_seed
      call random_number(WAVES(ng)%phase)
      twopi=2.0_r8*pi
# ifdef DISTRIBUTE
      IF (Master) THEN
# endif
        DO ff=1,WAVES(ng)%nfreq
          DO dd=1,WAVES(ng)%ndir
            WAVES(ng)%phase(ff,dd)=WAVES(ng)%phase(ff,dd)*twopi
          END DO
        END DO
# ifdef DISTRIBUTE
      END IF
# endif
# ifdef DISTRIBUTE
!
!  Scatter phase to all the nodes.
!
      Npts=WAVES(ng)%nfreq*WAVES(ng)%ndir
      allocate(wrk(Npts))
      IF (Master) THEN
        k=0
        DO ff=1,WAVES(ng)%nfreq
          DO dd=1,WAVES(ng)%ndir
            k=k+1
            wrk(k)=WAVES(ng)%phase(ff,dd)
          END DO
        END DO
      END IF
      CALL MPI_BCAST(wrk, Npts, MP_FLOAT, 0,                            &
     &               OCN_COMM_WORLD, MyError)
      k=0
      DO ff=1,WAVES(ng)%nfreq
        DO dd=1,WAVES(ng)%ndir
          k=k+1
          WAVES(ng)%phase(ff,dd)=wrk(k)
        END DO
      END DO
!
      deallocate(wrk)
# endif

      RETURN
      END SUBROUTINE random_phase

!
!***********************************************************************
      SUBROUTINE tpdcalc(Sf,f,frep)
!***********************************************************************
!
!  Here is only a calculation of the peak period.
!
      USE mod_inwave_bound

      implicit none

      real(r8), intent(in)  :: Sf(:), f(:)
      real(r8), pointer     :: temp(:)
      real(r8)              :: frep

      allocate(temp(size(Sf)))
      temp=0.0_r8
      where (Sf>0.8_r8*maxval(Sf))
       temp=1.0_r8
      end where

      frep=sum(temp*Sf*f)/sum(temp*Sf)
      deallocate(temp)

      RETURN
      END SUBROUTINE tpdcalc

!
!***********************************************************************
      SUBROUTINE amplitudes (ng)
!***********************************************************************
!
      USE mod_scalars
      USE mod_inwave_params
      USE mod_inwave_vars

      implicit none

!  Imported variable declarations.
!
      integer, intent(in) :: ng

!  Local variable declarations.
!
      integer                         :: ff, dd 
      real(r8)                        :: cff, cff1, cff2
!
!-----------------------------------------------------------------------
! Compute the amplitude for the Fourier components
! For each frequency in the spectra there is one amplitude
!-----------------------------------------------------------------------
!
      cff1=WAVES(ng)%df
      cff2=WAVEG(ng)%pd

# if defined DOUBLESUM
      WAVES(ng)%ampmax=0.0_r8
      DO dd=1,WAVES(ng)%ndir
        DO ff=1,WAVES(ng)%nfreq
          WAVES(ng)%amp(ff,dd)=sqrt(2.0_r8*WAVES(ng)%SD(ff,dd)*         &
     &                              cff1*cff2)
          IF (WAVES(ng)%amp(ff,dd).gt.WAVES(ng)%ampmax) THEN
            WAVES(ng)%ampmax=WAVES(ng)%amp(ff,dd)
          END IF
        END DO
      END DO
# elif defined SINGLESUM
      DO ff=1,WAVES(ng)%nfreq
        cff=0.0_r8
        DO dd=1,WAVES(ng)%ndir
          cff=cff+(WAVES(ng)%SD(ff,dd)*cff2)
        END DO
        WAVES(ng)%amp(ff)= sqrt(2.0_r8*cff*cff1)
      END DO
# endif

! Assigning a position in the spectral frequency array
! to each Fourier component

      DO ff=1,WAVES(ng)%nfreq
        WAVES(ng)%position(ff)=ff
      END DO

      RETURN
      END SUBROUTINE amplitudes
!
!***********************************************************************
      SUBROUTINE grid_metrics (ng, tile, first)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_inwave_params
      USE mod_grid

      implicit none

!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, first

!  Local variable declarations.
!
# include "tile.h"

      CALL grid_metrics_tile (ng, tile, first,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS)
      RETURN
      END SUBROUTINE grid_metrics
!
!***********************************************************************
      SUBROUTINE grid_metrics_tile (ng, tile, first,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_ncparam
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_inwave_params
      USE mod_inwave_vars
      USE mod_grid
      USE mod_ocean
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_gather2d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, first
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer                         :: i, j, m
      integer                         :: ff, dd, my_rank
      integer                         :: is, jcenter
      integer                         :: Npts, p1
      real(r8)                        :: x_mid, y_mid
      real(r8), allocatable           :: dist_local(:)

      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: Awrk
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: h_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: zeta_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: xr_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: yr_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: pm_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: pn_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: angle_local

# include "set_bounds.h"
!
!  Gather the depth and zeta to determine total depth at boundary point.
!
      DO j=0,Mm(ng)+1
        DO i=0,Lm(ng)+1
          h_local(i,j)=0.0_r8
          zeta_local(i,j)=0.0_r8
          xr_local(i,j)=0.0_r8
          yr_local(i,j)=0.0_r8
          pm_local(i,j)=0.0_r8
          pn_local(i,j)=0.0_r8
          angle_local(i,j)=0.0_r8
        END DO
      END DO

# ifdef DISTRIBUTE
      CALL mp_gather2d (ng, iNLM, LBi, UBi, LBj, UBj,                   &
     &                  0, r2dvar, 1.0_r8,                              &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  GRID(ng)%h, Npts, Awrk)
      CALL mp_bcastf (ng, iNLM, Awrk)
      p1=0
      DO j=0,Mm(ng)+1
        DO i=0,Lm(ng)+1
          p1=p1+1
          h_local(i,j)=Awrk(p1)
        END DO
      END DO
!
      CALL mp_gather2d (ng, iNLM, LBi, UBi, LBj, UBj,                   &
     &                  0, r2dvar, 1.0_r8,                              &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  OCEAN(ng)%zeta(:,:,1), Npts, Awrk)
      CALL mp_bcastf (ng, iNLM, Awrk)
      p1=0
      DO j=0,Mm(ng)+1
        DO i=0,Lm(ng)+1
          p1=p1+1
          zeta_local(i,j)=Awrk(p1)
        END DO
      END DO
!
      CALL mp_gather2d (ng, iNLM, LBi, UBi, LBj, UBj,                   &
     &                  0, r2dvar, 1.0_r8,                              &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  GRID(ng) % xr, Npts, Awrk)
      CALL mp_bcastf (ng, iNLM, Awrk)
      p1=0
      DO j=0,Mm(ng)+1
        DO i=0,Lm(ng)+1
          p1=p1+1
          xr_local(i,j)=Awrk(p1)
        END DO
      END DO
!
      CALL mp_gather2d (ng, iNLM, LBi, UBi, LBj, UBj,                   &
     &                  0, r2dvar, 1.0_r8,                              &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  GRID(ng) % yr, Npts, Awrk)
      CALL mp_bcastf (ng, iNLM, Awrk)
      p1=0
      DO j=0,Mm(ng)+1
        DO i=0,Lm(ng)+1
          p1=p1+1
          yr_local(i,j)=Awrk(p1)
        END DO
      END DO
!
      CALL mp_gather2d (ng, iNLM, LBi, UBi, LBj, UBj,                   &
     &                  0, r2dvar, 1.0_r8,                              &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  GRID(ng) % pm, Npts, Awrk)
      CALL mp_bcastf (ng, iNLM, Awrk)
      p1=0
      DO j=0,Mm(ng)+1
        DO i=0,Lm(ng)+1
          p1=p1+1
          pm_local(i,j)=Awrk(p1)
        END DO
      END DO
!
      CALL mp_gather2d (ng, iNLM, LBi, UBi, LBj, UBj,                   &
     &                  0, r2dvar, 1.0_r8,                              &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  GRID(ng) % pn, Npts, Awrk)
      CALL mp_bcastf (ng, iNLM, Awrk)
      p1=0
      DO j=0,Mm(ng)+1
        DO i=0,Lm(ng)+1
          p1=p1+1
          pn_local(i,j)=Awrk(p1)
        END DO
      END DO
!
      CALL mp_gather2d (ng, iNLM, LBi, UBi, LBj, UBj,                   &
     &                  0, r2dvar, 1.0_r8,                              &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  GRID(ng) % angler, Npts, Awrk)
      CALL mp_bcastf (ng, iNLM, Awrk)
      p1=0
      DO j=0,Mm(ng)+1
        DO i=0,Lm(ng)+1
          p1=p1+1
          angle_local(i,j)=Awrk(p1)
        END DO
      END DO
# else
      DO j=0,Mm(ng)+1
        DO i=0,Lm(ng)+1
          h_local(i,j)=GRID(ng)%h(i,j)
          zeta_local(i,j)=OCEAN(ng)%zeta(i,j,1)
          xr_local(i,j)=GRID(ng)%xr(i,j)
          yr_local(i,j)=GRID(ng)%yr(i,j)
          pm_local(i,j)=GRID(ng)%pm(i,j)
          pn_local(i,j)=GRID(ng)%pn(i,j)
          angle_local(i,j)=GRID(ng)%angle(i,j)
        END DO
      END DO
# endif
!
      IF (((LBC(iwest,isAC3d,ng)%acquire).and.                         &
     &    (DOMAIN(ng)%Western_Edge(tile))).or.                         &
     &    ((LBC(ieast,isAC3d,ng)%acquire).and.                         &
     &    (DOMAIN(ng)%Eastern_Edge(tile)))) THEN
!
!  Find the boundary cell closest to the center of the boundary
!
        IF (LBC(iwest,isAC3d,ng)%acquire) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            write(*,*) 'Computing western AC forcing'
          END IF
          is=0
        ELSE
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            write(*,*) 'Computing eastern AC forcing'
          END IF
          is=Lm(ng)+1
        END IF
        WAVES(ng)%SstrR=JstrR
        WAVES(ng)%SendR=JendR

        IF (first.eq.1) THEN
          allocate(WAVES(ng)%angle_edge(0:Mm(ng)+1))
        END IF
        DO j=0,Mm(ng)+1
          WAVES(ng)%angle_edge(j)=angle_local(is,j)
        END DO

        x_mid=0.5_r8*(xr_local(is,Mm(ng)+1)+xr_local(is,0))
        y_mid=0.5_r8*(yr_local(is,Mm(ng)+1)+yr_local(is,0))

        allocate(dist_local(0:Mm(ng)+1))
        DO j=0,Mm(ng)+1
          dist_local(j)=(yr_local(is,j)-y_mid)
        END DO
        jcenter=minloc(abs(dist_local), DIM = 1)
        WAVES(ng)%h0=h_local(is,jcenter)+zeta_local(is,jcenter)
        IF (DOMAIN(ng)%SouthWest_Corner(tile).or.                       &
     &      DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            write(*,*) 'Mean offshore water depth for AC is ',          &
     &                  WAVES(ng)%h0
        END IF
        deallocate(dist_local)

        IF (first.eq.1) THEN
          allocate(WAVES(ng)%dist_x(0:Mm(ng)+1))
          allocate(WAVES(ng)%dist_y(0:Mm(ng)+1))
          allocate(WAVES(ng)%dist_xy(0:Mm(ng)+1))
        END IF

        WAVES(ng)%dist_x(0)=0.0_r8
        WAVES(ng)%dist_y(0)=0.0_r8
        DO j=1,Mm(ng)+1
!         WAVES(ng)%dist_x(j)=xr_local(is,j)-xr_local(is,jcenter)
!         WAVES(ng)%dist_y(j)=yr_local(is,j)-yr_local(is,jcenter)
          WAVES(ng)%dist_x(j)=0.0_r8
          WAVES(ng)%dist_y(j)=WAVES(ng)%dist_y(j-1)+                    &
     &                        1.0_r8/pn_local(is,j)
        END DO
!try a new method here
        WAVES(ng)%dist_xy(0)=0.0_r8
        DO j=0,Mm(ng)+1
          WAVES(ng)%dist_xy(j)=sqrt( (xr_local(is,j)-xr_local(is,0))**2+&
     &                               (yr_local(is,j)-yr_local(is,0))**2)
        END DO
      END IF
!
      IF (((LBC(isouth,isAC3d,ng)%acquire).and.                         &
     &    (DOMAIN(ng)%Southern_Edge(tile))).or.                         &
     &    ((LBC(inorth,isAC3d,ng)%acquire).and.                         &
     &    (DOMAIN(ng)%Northern_Edge(tile)))) THEN
!
!  Find the boundary cell closest to the center of the boundary
!
        IF (LBC(isouth,isAC3d,ng)%acquire) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            write(*,*) 'Computing Southern AC forcing'
          END IF
          is=0
        ELSE
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            write(*,*) 'Computing Northern AC forcing'
          END IF
          is=Mm(ng)+1
        END IF
        WAVES(ng)%SstrR=IstrR
        WAVES(ng)%SendR=IendR
        
        IF (first.eq.1) THEN
          allocate(WAVES(ng)%angle_edge(0:Lm(ng)+1))
        END IF
        DO i=0,Lm(ng)+1
          WAVES(ng)%angle_edge(i)=angle_local(i,is)
        END DO

        x_mid=0.5_r8*(xr_local(Lm(ng)+1,is)+xr_local(0,is))
        y_mid=0.5_r8*(yr_local(Lm(ng)+1,is)+yr_local(0,is))

        allocate(dist_local(0:Lm(ng)+1))
        DO i=0,Lm(ng)+1
          dist_local(i)=(xr_local(i,is)-x_mid)
        END DO

        jcenter=minloc(abs(dist_local), DIM = 1)
        WAVES(ng)%h0=h_local(jcenter,is)+zeta_local(jcenter,is)

        IF (DOMAIN(ng)%SouthWest_Corner(tile).or.                       &
     &      DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            write(*,*) 'Mean offshore water depth for AC is ',          &
     &                  WAVES(ng)%h0
        END IF

        deallocate(dist_local)
        IF (first.eq.1) THEN
          allocate(WAVES(ng)%dist_x(0:Lm(ng)+1))
          allocate(WAVES(ng)%dist_y(0:Lm(ng)+1))
          allocate(WAVES(ng)%dist_xy(0:Lm(ng)+1))
        END IF

        WAVES(ng)%dist_x(0)=0.0_r8
        WAVES(ng)%dist_y(0)=0.0_r8
        DO i=1,Lm(ng)+1
!         WAVES(ng)%dist_x(i)=xr_local(i,is)-xr_local(jcenter,is)
!         WAVES(ng)%dist_y(i)=yr_local(i,is)-yr_local(jcenter,is)
          WAVES(ng)%dist_x(i)=WAVES(ng)%dist_x(i-1)+                    &
     &                        1.0_r8/pm_local(i,is)
          WAVES(ng)%dist_y(i)=0.0_r8
        END DO
!try a new method here
        WAVES(ng)%dist_xy(0)=0.0_r8
        DO i=1,Lm(ng)+1
          WAVES(ng)%dist_xy(i)=sqrt( (xr_local(i,is)-xr_local(0,is))**2+&
     &                               (yr_local(i,is)-yr_local(0,is))**2)
        END DO
      END IF

      RETURN
      END SUBROUTINE grid_metrics_tile
!
!
!***********************************************************************
      SUBROUTINE compute_k (ng, tile, first)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_inwave_params
      USE mod_grid

      implicit none

!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, first

!  Local variable declarations.
!
# include "tile.h"

      CALL compute_k_tile (ng, tile, first,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     WAVES(ng) % nfreq)
      RETURN
      END SUBROUTINE compute_k
!
!***********************************************************************
      SUBROUTINE compute_k_tile (ng, tile, first,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nfreq)
!***********************************************************************
!
      USE mod_param
      USE mod_inwave_params
      USE mod_inwave_vars
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, first
      integer, intent(in) :: nfreq
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer                         :: ff,error
      real(r8)                        :: cff1, cff2
      real(r8)                        :: twopi, otwopi
      real(r8)                        :: F, FDER, tanhkh
      real(r8)                        :: L0, k0, k1, kh, wr
      real(r8), parameter :: maxErr = 0.1_r8

# include "set_bounds.h"

      twopi=2.0_r8*pi
      otwopi=1.0_r8/(2.0_r8*pi)
!
!  Compute the wave numbers.
!
!  Only compute Ampzeta for the boundary tiles.
!
      IF (first.eq.1) THEN
        allocate (WAVES(ng)%waven(nfreq))
      END IF
      WAVES(ng)%waven=0.0_r8
!
      IF (  ((LBC(iwest,isAC3d,ng)%acquire).and.                        &
     &       (DOMAIN(ng)%Western_Edge(tile)) ).or.                      &
     &      ((LBC(ieast,isAC3d,ng)%acquire).and.                        &
     &       (DOMAIN(ng)%Eastern_Edge(tile))).or.                       &
     &      ((LBC(isouth,isAC3d,ng)%acquire).and.                       &
     &       (DOMAIN(ng)%Southern_Edge(tile))).or.                      &
     &      ((LBC(inorth,isAC3d,ng)%acquire).and.                       &
     &       (DOMAIN(ng)%Northern_Edge(tile)))) THEN
        DO ff=1,nfreq
          L0=g*otwopi*(1.0_r8/WAVES(ng)%f(ff))**2
          k0=twopi/L0
          error=100.0_r8
          wr=twopi*WAVES(ng)%f(ff)
          DO WHILE(error.gt.maxErr)
            kh=k0*WAVES(ng)%h0
            tanhkh=TANH(kh)
            cff1=wr**2.0_r8
            cff2=-g*k0*tanhkh
            F=cff1+cff2
            cff1=-g*tanhkh
            cff2=-g*kh/COSH(kh)**2.0_r8
            FDER=cff1+cff2
            k1=k0-F/FDER
            error=100.0_r8*ABS((k1-k0)/k0)
            k0=k1
          END DO
          WAVES(ng)%waven(ff)=k0
        END DO
      END IF
      RETURN
      END SUBROUTINE compute_k_tile
!
!
!***********************************************************************
      SUBROUTINE FourierComp (ng, tile, first)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_inwave_params
      USE mod_inwave_vars
      USE mod_grid

      implicit none

!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, first

!  Local variable declarations.
!
#  include "tile.h"

      CALL FourierComp_tile (ng, tile, first,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       ND, WAVES(ng) % SstrR, WAVES(ng) % SendR,  &
     &                       WAVEG(ng) % wd,                            &
     &                       WAVES(ng) % Insteps,                       &
     &                       WAVES(ng) % nfreq,                         &
     &                       WAVES(ng) % ndir,                          &
     &                       WAVES(ng) % amp,                           &
#  ifdef SINGLESUM
#   ifdef INTCDF
     &                       WAVES(ng) % SDD,                           &
#   endif
#   ifdef FCDF
     &                       WAVES(ng) % SD,                            &
#   endif
#  endif
     &                       WAVES(ng) % dir_3_int,                     &
     &                       WAVES(ng) % waven,                         &
     &                       WAVES(ng) % phase)
      RETURN
      END SUBROUTINE FourierComp
!
!***********************************************************************
      SUBROUTINE FourierComp_tile (ng, tile, first,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             ND, SstrR, SendR, wd,                &
     &                             Insteps, nfreq, ndir,                &
     &                             amp,                                 &
#  ifdef SINGLESUM
#   ifdef INTCDF
     &                             SDD,                                 &
#   endif
#   ifdef FCDF
     &                             SD,                                  &
#   endif
#  endif
     &                             dir_3_int, waven,                    &
     &                             phase)
!***********************************************************************
!
!     USE mod_ncparam
      USE mod_param
      USE mod_parallel
      USE mod_scalars
!     USE mod_inwave_params
      USE mod_inwave_vars
      USE mod_grid
      USE mod_ocean
#  ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_gather2d
#  endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, first
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: ND, SstrR, SendR
      integer, intent(in) :: Insteps, nfreq, ndir
      integer, intent(inout) :: dir_3_int(nfreq)

#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: wd(ND)
#   ifdef SINGLESUM
      real(r8), intent(in) :: amp(nfreq)
#   elif defined DOUBLESUM
      real(r8), intent(in) :: amp(nfreq,ndir)
#   endif
#   ifdef SINGLESUM
#    ifdef INTCDF
      real(r8), intent(in) :: SDD(ndir)
#    endif
#    ifdef FCDF
      real(r8), intent(in) :: SD(nfreq,ndir)
#    endif
#   endif
      real(r8), intent(in) :: waven(nfreq)
      real(r8), intent(in) :: phase(nfreq,ndir)
#  else
      real(r8), intent(in) :: wd(ND)
#   ifdef SINGLESUM
      real(r8), intent(in) :: amp(nfreq)
#   elif defined DOUBLESUM
      real(r8), intent(in) :: amp(nfreq,ndir)
#   endif
#   ifdef SINGLESUM
#    ifdef INTCDF
      real(r8), intent(in) :: SDD(ndir)
#    endif
#    ifdef FCDF
      real(r8), intent(in) :: SD(nfreq,ndir)
#    endif
#   endif
      real(r8), intent(in) :: waven(nfreq)
      real(r8), intent(in) :: phase(nfreq,ndir)
#  endif
!
!  Local variable declarations.
!
      integer                         :: i, ii, j, k, m, tidx, nn
      integer                         :: ff, dd, MyError
      integer                         :: index1, index2, numdirs
      real(r8)                        :: cff, cff1, cff2, DIR3, kwx, kwy
      real(r8)                        :: twopi, otwopi, x, cff3, kedge
      real(r8)                        :: theta_pos, fac1, fac2
#  ifdef GAUSSIAN_DIST
      real(r8)                        :: rsq, meanx, sigma, xx
      real(r8)                        :: r1, r2, fac, tmp
#  endif
#  ifdef INTCDF
      real(r8), allocatable           :: SDD_cumul(:)
#  elif defined FCDF
      real(r8), allocatable           :: SDD_cumul(:,:)
#  endif
#  ifdef HSCALE
      real(r8), allocatable           :: Hsig_check(:)
#  endif
      real(r8), allocatable           :: atemp(:)

      real(r8), parameter :: maxErr = 0.1_r8
      real(r8), parameter :: alphad = 0.2_r8

#  include "set_bounds.h"

      twopi=2.0_r8*pi
      otwopi=1.0_r8/(2.0_r8*pi)
!
!  If singlesum, then we need to compute dir for each freq.
!
#  ifdef SINGLESUM
!
!  Create a cumulative distribution function (CDF) of energy.
!  then we select directional bins based on uniform distribution.
!
#   ifdef INTCDF
!     Create a CDF integrated over frequencies. 
!     SDD is directional energy spectra (but not x pd).
      allocate(SDD_cumul(ndir))
      DO dd=1,ndir
        SDD_cumul(dd)=0.0_r8
      END DO
      DO dd=2,ndir
          SDD_cumul(dd)=SDD_cumul(dd-1)+SDD(dd)
      END DO
!
!     Scale SDD_cumul so it goes from ~0.0 to 1.0
!
      DO dd=1,ndir
        IF (SDD_cumul(ndir).gt.0.0_r8) THEN
          SDD_cumul(dd)=SDD_cumul(dd)/SDD_cumul(ndir)
        END IF
      END DO
#   elif defined FCDF
      allocate(SDD_cumul(nfreq,ndir))
      DO ff=1,nfreq
        DO dd=1,ndir
          SDD_cumul(ff,dd)=0.0_r8
        END DO
        DO dd=2,ndir
          SDD_cumul(ff,dd)=SDD_cumul(ff,dd-1)+SD(ff,dd)
        END DO
      END DO
!
!  Scale SDD_cumul so it goes from ~0.0 to 1.0
!
      DO ff=1,nfreq
        DO dd=1,ndir
          IF (SDD_cumul(ff,ndir).gt.0.0_r8) THEN
            SDD_cumul(ff,dd)=SDD_cumul(ff,dd)/SDD_cumul(ff,ndir)
          END IF
        END DO
      END DO
#   endif
!
      call random_seed
      DO ff=1,nfreq
!
!  Determine which direction to use for this freq.
!  Set default to 1.
        m=1
        dir_3_int(ff)=m
#   ifdef DISTRIBUTE
        IF (Master) THEN
#   endif
        call random_number(cff1)
#   ifdef GAUSSIAN_DIST
!
!  Here we can try a normal distribution instead of uniform.
!
        rsq=2.0_r8
        meanx=0.5_r8
        sigma=0.30_r8
        DO WHILE ((rsq.ge.1.0_r8).or.(rsq.eq.0.0_r8))
          call random_number(xx)
          r1=2.0_r8*xx-1.0_r8
          call random_number(xx)
          r2=2.0_r8*xx-1.0_r8
          rsq=r1*r1+r2*r2
        END DO
        fac=sqrt(-2.0_r8*log(rsq)/rsq);
        tmp=r2*fac;
        cff1=tmp*sigma+meanx;
        if (cff1.le.0.0_r8) cff1=cff1+1.0_r8
        if (cff1.ge.1.0_r8) cff1=cff1-1.0_r8
#   endif
!
#   ifdef INTCDF
        DO dd=1,ndir
          IF (cff1.lt.SDD_cumul(dd)) THEN
            m=dd
            EXIT
          END IF
        END DO
        dir_3_int(ff)=m
#   elif defined FCDF
        DO dd=1,ndir
          IF (cff1.lt.SDD_cumul(ff,dd)) THEN
            m=dd
            EXIT
          END IF
        END DO
        dir_3_int(ff)=m
#   endif
#   ifdef DISTRIBUTE
        END IF
#   endif
      END DO
#   ifdef DISTRIBUTE
!
!  Scatter rand num from master to all nodes.
!
      CALL MPI_BCAST(dir_3_int, nfreq,                                  &
     &               MPI_INTEGER, 0, OCN_COMM_WORLD, MyError)
#   endif
#  endif
!
!  Only compute Ampzeta for the boundary tiles.
!
      IF (  ((LBC(iwest,isAC3d,ng)%acquire).and.                        &
     &       (DOMAIN(ng)%Western_Edge(tile)) ).or.                      &
     &      ((LBC(ieast,isAC3d,ng)%acquire).and.                        &
     &       (DOMAIN(ng)%Eastern_Edge(tile))).or.                       &
     &      ((LBC(isouth,isAC3d,ng)%acquire).and.                       &
     &       (DOMAIN(ng)%Southern_Edge(tile))).or.                      &
     &      ((LBC(inorth,isAC3d,ng)%acquire).and.                       &
     &       (DOMAIN(ng)%Northern_Edge(tile)))) THEN
!
!  Compute the amplitude time series for the boundary AC forcing.
!  Loop over all frequencies
!
        DO ff=1,nfreq
!jcw      DO ff=18,900
          DO i=SstrR,SendR
!
!  If DOUBLESUM, double summation loop over all directions.
!  If SINGLESUM, check if the random direction assigned above.
!
            DO dd=1,ndir
#  ifdef SINGLESUM
              IF (dd.eq.dir_3_int(ff)) THEN
#  endif
#  ifdef VARY_ACBC
                ii=i
#  else
                ii=1                     ! use same location for
                                         ! distx disty and angle_edge
#  endif

                DIR3=1.5_r8*pi-WAVES(ng)%dir(dd)-                       &
     &                         WAVES(ng)%angle_edge(ii)
!
                kwx=waven(ff)*cos(DIR3)
                kwy=waven(ff)*sin(DIR3)
                kedge=kwx*WAVES(ng)%dist_x(ii)+kwy*WAVES(ng)%dist_y(ii)
                DO tidx=1,Insteps
                  cff=REAL((tidx-1),r8)*WAVES(ng)%idt
#  ifdef SINGLESUM
                  cff2=amp(ff)*                                         &
     &                 COS(twopi*WAVES(ng)%f(ff)*cff-                   &
     &                 kedge+phase(ff,dd))
#  elif defined DOUBLESUM
                  IF (amp(ff,dd).gt.0.10_r8*WAVES(ng)%ampmax) THEN
                    cff2=amp(ff,dd)*                                    &
     &                      COS(twopi*WAVES(ng)%f(ff)*cff-              &
     &                      kedge+phase(ff,dd))
                   ELSE
                    cff2=0.0_r8
                  END IF
#  endif
                  WAVES(ng)%Ampzeta(i,tidx,dd)=                         &
     &                          WAVES(ng)%Ampzeta(i,tidx,dd)+cff2
                END DO
#  ifdef SINGLESUM
              END IF
#  endif
            END DO
          END DO
        END DO
!
!  Take the Hilbert transform of the amplitude time series 
#  ifdef DHILBERT
!  Maintain the directional distribution
!
        DO i=SstrR,SendR
          DO dd=1,ndir
            DO tidx=1,Insteps
              WAVES(ng)%Ampzeta_tot(tidx)=WAVES(ng)%Ampzeta(i,tidx,dd)
            END DO
            call hilbert(WAVES(ng)%Ampzeta_tot,                         &
     &                   size(WAVES(ng)%Ampzeta_tot))
            DO tidx=1,Insteps
              cff1=(ABS(WAVES(ng)%Ampzeta_tot(tidx)))**2
              WAVES(ng)%Ampzeta(i,tidx,dd)=0.5_r8*cff1*g*rho0
            END DO
          END DO
        END DO
#  elif defined IHILBERT
!  Integrated 
!
#   ifdef HSCALE
        allocate(Hsig_check(SstrR:SendR))
#   endif
        DO i=SstrR,SendR
          DO tidx=1,Insteps
            WAVES(ng)%Ampzeta_tot(tidx)=0.0_r8
            DO dd=1,ndir
              WAVES(ng)%Ampzeta_tot(tidx)=WAVES(ng)%Ampzeta_tot(tidx)+  &
     &                                    WAVES(ng)%Ampzeta(i,tidx,dd)
            END DO
          END DO
          call hilbert(WAVES(ng)%Ampzeta_tot,                           &
     &                 size(WAVES(ng)%Ampzeta_tot))
!
!  Partition back out to all directions
!
#   ifdef HSCALE
          cff=0.0_r8
#   endif
          DO dd=1,ndir
            cff2=SDD(dd)/WAVES(ng)%integ
            DO tidx=1,Insteps
              cff1=((ABS(WAVES(ng)%Ampzeta_tot(tidx)))**2)*cff2     
              WAVES(ng)%Ampzeta(i,tidx,dd)=0.5_r8*cff1*g*rho0
#   ifdef HSCALE
              cff=cff+cff1
#   endif
            END DO
          END DO
#   ifdef HSCALE
          Hsig_check(i)=4.0_r8*sqrt(cff/Insteps)
#   endif
        END DO
!
#   ifdef HSCALE
!
!  Scale energy so the time-averaged Hm0 is equivalent to SWAN spectra 
!
        DO i=SstrR,SendR
          cff=WAVES(ng)%hm0gew/Hsig_check(i)
          DO dd=1,ndir
            DO tidx=1,Insteps
              WAVES(ng)%Ampzeta(i,tidx,dd)=                             &
     &               WAVES(ng)%Ampzeta(i,tidx,dd)*cff
            END DO
          END DO
        END DO
        deallocate(Hsig_check)
#   endif
#  endif
#  ifdef SINGLESUM
        deallocate(SDD_cumul)
#  endif
      END IF
      RETURN
      END SUBROUTINE FourierComp_tile

!***********************************************************************
      SUBROUTINE boundwave (ng, tile, first)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, first, tile
!
!  Local variable declarations.
!
!     integer :: Insteps

# include "tile.h"

      CALL boundwave_tile (ng, tile, first,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     WAVES(ng) % SstrR, WAVES(ng) % SendR,        &
     &                     WAVES(ng) % Insteps,                         &
     &                     WAVES(ng) % nfreq,                           &
     &                     WAVES(ng) % ndir,                            &
     &                     WAVES(ng) % dir_3_int,                       &
     &                     WAVES(ng) % waven,                           &
     &                     WAVES(ng) % h0)
      RETURN
      END SUBROUTINE boundwave
!
!***********************************************************************
      SUBROUTINE boundwave_tile (ng, tile, first,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          SstrR, SendR,                           &
     &                          Insteps, nfreq, ndir,                   &
     &                          dir_3_int,                              &
     &                          waven,                                  &
     &                          h0)
!                                                                       !
!                compute the wave envelope and the associated bound     !
!                wave using a double summation technique.               !
!                (Hasselman, 1962; Herbers et al., 1994;                !
!                 Van Dongeren et al., 2003).                           !
!***********************************************************************

      USE mod_boundary
      USE mod_grid
      USE mod_ocean
      USE mod_parallel
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_inwave_vars
      USE math_tools
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcastf, mp_gather2d
# endif
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, first
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: SstrR, SendR
      integer, intent(in) :: Insteps, nfreq, ndir
      integer, intent(in) :: dir_3_int(nfreq)
      real(r8), intent(in) :: h0

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: waven(nfreq)
# else
      real(r8), intent(in) :: waven(nfreq)
# endif
!
!  Local variable declarations.
!
      integer                         :: i, j, jj, t, is, p1, p2, f1, f2
      integer                         :: n1, n2, nt1
      integer                         :: mdir, nf, nn
      integer                         :: mdir1
      integer                         :: m, ier, lensav, lenwrk
      real(r8)                        :: cff, cff1, cff2, cff3, cff4
      real(r8)                        :: twopi, otwopi, og2, p
      real(r8)                        :: error, fmin, fmax, deltaf
      real(r8)                        :: fac1, fac2, fac3, fac4, ofac4
      real(r8)                        :: ont1, oh0
# if !defined SINGLESUM
      real(r8)                        :: DDtheta, k3, Z_bw
      real(r8)                        :: k3_x, k3_y
      real(r8)                        :: D1, D2, D3, D4, DTOT
      real(r8)                        :: E3, A3
# endif
      real(r8)                        :: D1a, D3a
      real(r8)                        :: DIR3, fac5
      real(r8), dimension(WAVES(ng)%Insteps) :: temp4

      real(r8), allocatable           :: owaven(:)
      real(r8), allocatable           :: ocoshkh(:)
      real(r8), allocatable           :: tempz(:,:)
      real(r8), allocatable           :: tempu(:,:)
      real(r8), allocatable           :: tempv(:,:)
# if defined SINGLESUM
      real(r8), allocatable           :: DDtheta(:)
      real(r8), allocatable           :: k3(:)
      real(r8), allocatable           :: Z_bw(:,:)
      real(r8), allocatable           :: bw_ang(:,:)
      real(r8), allocatable           :: k3_x(:,:)
      real(r8), allocatable           :: k3_y(:,:)
      real(r8), allocatable           :: theta_3(:,:)
      real(r8), allocatable           :: cg_3(:,:)
      real(r8), allocatable           :: ff(:)
      real(r8), allocatable           :: DDf(:)
      real(r8), allocatable           :: D1(:)
      real(r8), allocatable           :: D2(:)
      real(r8), allocatable           :: D3(:)
      real(r8), allocatable           :: D4(:)
      real(r8), allocatable           :: DTOT(:)
      real(r8), allocatable           :: E3(:)
      real(r8), allocatable           :: A3(:,:)
!     real(r8), allocatable           :: Comptemp(:,:)
!     real(r8), allocatable           :: Comptempd(:,:)

      complex(fftkind), dimension(:,:), allocatable :: Comptemp
      complex(fftkind), dimension(:,:), allocatable :: Comptempd
      complex(fftkind), dimension(:), allocatable :: A, AA
      complex(fftkind), dimension(:), allocatable :: Comp_int

      real(r8), allocatable           :: zz(:)
      real(r8), allocatable           :: dpi(:)
      real(r8), allocatable           :: comp1(:)
      real(kind=8),dimension(:),allocatable    :: work, wsave
# endif

      real(r8), parameter :: maxErr = 0.1_r8
      real(r8), parameter :: eps = 1.0E-10_r8
      real(r8), parameter :: fcutoff = 1.0_r8/30.0_r8

# include "set_bounds.h"
      twopi=2.0_r8*pi
      otwopi=1.0_r8/(2.0_r8*pi)
      fmin=1.0_r8/400.0_r8
      fmax=1.0_r8/30.0_r8
      og2=1.0_r8/(g**2)
      oh0=1.0_r8/h0
!
      IF (  ((LBC(iwest,isAC3d,ng)%acquire).and.                        &
     &       (DOMAIN(ng)%Western_Edge(tile)) ).or.                      &
     &      ((LBC(ieast,isAC3d,ng)%acquire).and.                        &
     &       (DOMAIN(ng)%Eastern_Edge(tile))).or.                       &
     &      ((LBC(isouth,isAC3d,ng)%acquire).and.                       &
     &       (DOMAIN(ng)%Southern_Edge(tile))).or.                      &
     &      ((LBC(inorth,isAC3d,ng)%acquire).and.                       &
     &       (DOMAIN(ng)%Northern_Edge(tile)))) THEN
!
!  Allocate boundwave and local wavenumber arrays.
!
        IF (first.eq.1) THEN
          allocate(WAVES(ng)%bndwave_z(SstrR:SendR,Insteps))
          allocate(WAVES(ng)%bndwave_u(SstrR:SendR,Insteps))
          allocate(WAVES(ng)%bndwave_v(SstrR:SendR,Insteps))
        END IF
        DO i=SstrR,SendR
          DO t=1,Insteps
            WAVES(ng)%bndwave_z(i,t)=0.0_r8
            WAVES(ng)%bndwave_u(i,t)=0.0_r8
            WAVES(ng)%bndwave_v(i,t)=0.0_r8
          END DO
        END DO
!
        allocate (owaven(nfreq))
!
! Compute the wave number of each freq bin
!
        DO i=1,nfreq
          owaven(i)=1.0_r8/waven(i)
        END DO
!
!  make some temp arrays
!
        allocate (ocoshkh(nfreq))
        DO i=1,nfreq
          fac2=cosh(waven(i)*h0)
          ocoshkh(i)=1.0_r8/fac2
        END DO
!
! Compute the energy transfer for each pair of frequency bins.
!
        allocate (tempz(SstrR:SendR,Insteps))
        allocate (tempu(SstrR:SendR,Insteps))
        allocate (tempv(SstrR:SendR,Insteps))
        DO t=1,Insteps
          temp4(t)=REAL(t-1,r8)*WAVES(ng)%idt
          DO j=SstrR,SendR
            tempz(j,t)=0.0_r8
            tempu(j,t)=0.0_r8
            tempv(j,t)=0.0_r8
          END DO
        END DO
        fac4=WAVES(ng)%df
        ofac4=1.0_r8/WAVES(ng)%df
        fac5=WAVEG(ng)%pd
        og2=1.0_r8/(g**2)
!
!  Do the big loops to compute bound.
!
# ifdef SINGLESUM
!
!  New way to compute bound wave with less loops.  This is a single 
!  summation approach.
!
!  Loop over all frequencies
!
! m to f1, k to nf
        nt1=WAVES(ng)%Insteps
        ont1=1.0_r8/REAL(Insteps,r8)
        nf=nfreq

        allocate (DDtheta(nf-1))
        allocate (k3(nf-1))
        allocate (Z_bw(nf-1,nf-1))
        allocate (bw_ang(nf-1,nf-1))
        allocate (k3_x(nf-1,nf-1))
        allocate (k3_y(nf-1,nf-1))
        allocate (theta_3(nf-1,nf-1))
        allocate (cg_3(nf-1,nf-1))
        allocate (ff(nf-1))
        allocate (DDf(nf-1))
        allocate (D1(nf-1))
        allocate (D2(nf-1))
        allocate (D3(nf-1))
        allocate (D4(nf-1))
        allocate (DTOT(nf-1))
        allocate (E3(nf-1))
        allocate (A3(nf-1,nf-1))
        allocate (Comptemp(nf-1,nf-1))
        allocate (Comptempd(nf-1,nf-1))
        allocate (Comp_int(nf))
        allocate (A(nt1))
        allocate (zz(nt1))

        DDtheta=0.0_r8
        k3=0.0_r8
        Z_bw=0.0_r8
        bw_ang=0.0_r8
        k3_x=0.0_r8
        k3_y=0.0_r8
        theta_3=0.0_r8
        cg_3=0.0_r8
        ff=0.0_r8
        DDf=0.0_r8
        D1=0.0_r8
        D2=0.0_r8
        D3=0.0_r8
        D4=0.0_r8
        DTOT=0.0_r8
        E3=0.0_r8
        A3=0.0_r8
        Comptemp=0.0_r8
        Comptempd=0.0_r8
        Comp_int=0.0_r8
        A=0.0_r8
        zz=0.0_r8
!
        DO f1=1,nf-1
          allocate (dpi(1:nf-f1))
          DTOT=0.0_r8
          E3=0.0_r8
          DDtheta=0.0_r8
          k3=0.0_r8
          ff=0.0_r8
          DDf=0.0_r8
          D1=0.0_r8
          D2=0.0_r8
          D3=0.0_r8
          D4=0.0_r8
!
!  Report percent done
!
          IF ((DOMAIN(ng)%SouthWest_Corner(tile).and.                   &
              (LBC(isouth,isAC3d,ng)%acquire.or.                        &
               LBC(iwest,isAC3d,ng)%acquire)) .or.                      &
             (DOMAIN(ng)%NorthEast_Corner(tile).and.                    &
              (LBC(inorth,isAC3d,ng)%acquire.or.                        &
               LBC(ieast,isAC3d,ng)%acquire))) THEN
            IF (MOD(f1*1000/nfreq,100).eq.0) THEN
              write(*,*) 'percent done for computing bndwave is ',      &
     &                    f1*100/nfreq
            END IF
          END IF
!
! Compute the freq dif
!
          deltaf=REAL(f1,r8)*WAVES(ng)%df
!
!         direction difference

!         DDtheta(1:nf-f1)=(WAVEG(ng)%wd(dir_3_int(f1+1:nf))-           &
!    &                      WAVEG(ng)%wd(dir_3_int(1:nf-f1)))
          DDtheta(1:nf-f1)=(WAVES(ng)%dir(dir_3_int(f1+1:nf))-          &
     &                      WAVES(ng)%dir(dir_3_int(1:nf-f1)))
!
          dpi=1.0_r8*pi
!
          DDtheta(1:nf-f1)=ABS(DDtheta(1:nf-f1))+dpi
!
!         wave number, -dpi is needed becasue we added it above.
          k3(1:nf-f1)=sqrt(waven(f1+1:nf)**2+waven(1:nf-f1)**2-          &
     &                   2.0_r8*(waven(f1+1:nf)*waven(1:nf-f1))*         &
     &                   (cos(DDtheta(1:nf-f1)-dpi)))
!
!  Phase of the bound wave. Use the m value from INTCDF or FCDF above.
!
          allocate (comp1(1:nf))
          DO i=1,nf
            mdir1=dir_3_int(i)
            comp1(i)=WAVES(ng)%phase(i,mdir1)
          END DO
          Z_bw(f1,1:nf-f1)=comp1(f1+1:nf)-comp1(1:nf-f1)+pi
          deallocate (comp1)
!
!  bw_ang is in radians and nautical convention
!
!          bw_ang(f1,1:nf-f1)=atan2(                                     &
!     &             waven(f1+1:nf)*sin(WAVEG(ng)%wd(                     &
!     &                               dir_3_int(f1+1:nf)))-              &
!     &             waven(1:nf-f1)*sin(WAVEG(ng)%wd(                     &
!     &                               dir_3_int(1:nf-f1))),              &
!     &             waven(f1+1:nf)*cos(WAVEG(ng)%wd(                     &
!     &                               dir_3_int(f1+1:nf)))-              &
!     &             waven(1:nf-f1)*cos(WAVEG(ng)%wd(                     &
!     &                               dir_3_int(1:nf-f1))))
          bw_ang(f1,1:nf-f1)=atan2(                                      &
     &             waven(f1+1:nf)*sin(WAVES(ng)%dir(                     &
     &                               dir_3_int(f1+1:nf)))-               &
     &             waven(1:nf-f1)*sin(WAVES(ng)%dir(                     &
     &                               dir_3_int(1:nf-f1))),               &
     &             waven(f1+1:nf)*cos(WAVES(ng)%dir(                     &
     &                               dir_3_int(f1+1:nf)))-               &
     &             waven(1:nf-f1)*cos(WAVES(ng)%dir(                     &
     &                               dir_3_int(1:nf-f1))))
!
!  Compute wave numbers in x and y
!  The 1.5pi is beasue dir angles are nautical conv, from the direction.
!
          k3_x(f1,1:nf-f1)=k3(1:nf-f1)*cos(1.5_r8*pi-bw_ang(f1,1:nf-f1))
          k3_y(f1,1:nf-f1)=k3(1:nf-f1)*sin(1.5_r8*pi-bw_ang(f1,1:nf-f1))
          cg_3(f1,1:nf-f1)=twopi*deltaf/k3(1:nf-f1)
!
!  Compute df*df.  watch minus signs, it is -f1 + f2
          ff(1:nf-f1)=WAVES(ng)%f(1:nf-f1)*WAVES(ng)%f(f1+1:nf)
          DDf(1:nf-f1)=-WAVES(ng)%f(1:nf-f1)+WAVES(ng)%f(f1+1:nf)
!
!  Compute the D terms
!
! ??      ocoshkh=1.0_r8/cosh(k*h0) - several places below D4 D1
          D1(1:nf-f1)=g*waven(f1+1:nf)*waven(1:nf-f1)*                  &
     &                 cos(DDtheta(1:nf-f1))/                           &
     &                (8.0_r8*pi**2*ff(1:nf-f1))
!  correct for bottom pressure.
          D1(1:nf-f1)=D1(1:nf-f1)*                                      &
     &                 MIN(cosh(k3(1:nf-f1)*h0),100.0_r8)  /            &
     &                 (cosh(waven(f1+1:nf)*h0)*cosh(waven(1:nf-f1)*h0))
          D2(1:nf-f1)=-g*DDf(1:nf-f1)/((g*k3(1:nf-f1)*                  &
     &                tanh(k3(1:nf-f1)*h0)-                             &
     &                (twopi)**2*DDf(1:nf-f1)**2)*ff(1:nf-f1)+eps)
          D3(1:nf-f1)=DDf(1:nf-f1)*((twopi**4*(ff(1:nf-f1))**2)*og2-    &
     &                waven(f1+1:nf)*waven(1:nf-f1)*                    &
     &                cos(DDtheta(1:nf-f1))) 
          D4(1:nf-f1)=-0.5*((-WAVES(ng)%f(1:nf-f1)*waven(f1+1:nf)**2)/  &
     &                      (cosh(waven(f1+1:nf)*h0))**2+               &
     &                      ( WAVES(ng)%f(f1+1:nf)*waven(1:nf-f1)**2)/  &
     &                      (cosh(waven(1:nf-f1)*h0))**2)
!
!  Compute the double summation.
!
          DTOT(1:nf-f1)=D1(1:nf-f1)+D2(1:nf-f1)*                        &
     &                 (D3(1:nf-f1)+D4(1:nf-f1))
!
!  do not consider frequencies higher than fcutoff.
!
          WHERE (WAVES(ng)%f(1:nf-1).le.deltaf)
            DTOT(:)=0.0_r8
          END WHERE
          IF (deltaf.le.fmin) THEN
            DTOT(:)=0.0_r8
          END IF
          IF (deltaf.ge.fmax) THEN
            DTOT(:)=0.0_r8
          END IF
!
!  compute Energy and Amplitude of bound wave.
!
          E3(1:nf-f1)=2.0_r8*(DTOT(1:nf-f1)**2)*                        &
     &              (WAVES(ng)%SDF(f1+1:nf)*WAVES(ng)%SDF(1:nf-f1)*     &
     &               fac5*fac5*fac4)
          A3(f1,1:nf-f1)=sqrt(2.0_r8*E3(1:nf-f1)*fac4)

          deallocate (dpi)
        END DO
!
!  Compute 
!
        Comptemp=0.5_r8*A3*exp(-1.0_r8*CMPLX(0.0_r8,Z_bw))
!
!  Compute this contribution to the bndwave (water level) and associated 
!  currents; currents depend on direction
!  Here we need to include the phase difference along the boundary;
!
!  These are for the fourier computation below.
!
        p=log(real(nt1))/log(2.)
        nn=ceiling(p)
        nn=2**nn
        lenwrk=2*nn
        lensav=2*nn+int(log(real(nn,kind=8))/log(2.0D+00))+4
        allocate ( work(1:lenwrk) )
        allocate ( wsave(1:lensav) )
        allocate ( AA(nn) )
        work=0.
        wsave=0.
        call cfft1i ( nt1, wsave, lensav, ier )
        DO j=SstrR,SendR
#  ifdef VARY_ACBC
          jj=j
#  else
          jj=1                     ! use same location for
                                   ! distx disty and angle_edge
#  endif
!
!  first get the bound_z term.
          Comptempd=Comptemp*exp(-1.0_r8*CMPLX(0.0_r8,                  &
     &                           k3_y*WAVES(ng)%dist_y(jj)+             &
     &                           k3_x*WAVES(ng)%dist_x(jj)))
!
!  Determine Fourier coefficients
!
!  we add all the components of same frequency
!  need to start at 2 for the fft  conj symmetry
           Comp_int(2:nf)=sum(Comptempd,2)
!
!  Compute the Fourier Coefficients
!
          A=0.0_r8
          AA=0.0_r8
          A(1:nf)=Comp_int
          A(nt1-nf+2:nt1)=conjg(A(nf:2:-1))
          AA(1:nt1)=A
!
          work=0.
          call cfft1b (nt1, 1, AA, nn, wsave, lensav, work, lenwrk, ier)
!
          zz=0.0_r8
          zz=real(AA(1:nt1),r8)
          tempz(j,:)=zz
!  now do u
          Comptempd=Comptemp*exp(-1.0_r8*CMPLX(0.0_r8,                  &
     &                           k3_y*WAVES(ng)%dist_y(jj)+             &
     &                           k3_x*WAVES(ng)%dist_x(jj)))*           &
     &                           cg_3*cos(1.5_r8*pi-bw_ang)
!
!  Determine Fourier coefficients
!
!  we add all the components of same frequency
!  need to start at 2 for the fft  conj symmetry
           Comp_int(2:nf)=sum(Comptempd,2)
!
!  Compute the Fourier Coefficients
!
          A=0.0_r8
          AA=0.0_r8
          A(1:nf)=Comp_int
          A(nt1-nf+2:nt1)=conjg(A(nf:2:-1))
          AA(1:nt1)=A
!
          work=0.
          call cfft1b (nt1, 1, AA, nn, wsave, lensav, work, lenwrk, ier)
!
          zz=0.0_r8
          zz=real(AA(1:nt1),r8)

          tempu(j,:)=zz
!  now do v
          Comptempd=Comptemp*exp(-1.0_r8*CMPLX(0.0_r8,                  &
     &                           k3_y*WAVES(ng)%dist_y(jj)+             &
     &                           k3_x*WAVES(ng)%dist_x(jj)))*           &
     &                           cg_3*sin(1.5_r8*pi-bw_ang)
!
!  Determine Fourier coefficients
!
!  we add all the components of same frequency
!  need to start at 2 for the fft  conj symmetry
           Comp_int(2:nf)=sum(Comptempd,2)
!
!  Compute the Fourier Coefficients
!
          A=0.0_r8
          AA=0.0_r8
          A(1:nf)=Comp_int
          A(nt1-nf+2:nt1)=conjg(A(nf:2:-1))
          AA(1:nt1)=A
!
          work=0.
          call cfft1b (nt1, 1, AA, nn, wsave, lensav, work, lenwrk, ier)
!
          zz=0.0_r8
          zz=real(AA(1:nt1),r8)

          tempv(j,:)=zz
        END DO
        deallocate (work, wsave)
!
!  The currents need to be rotated to the roms grid plus they are
!  defined in the rho points. Keep the u v as fluxes and convert in
!  the bndwave bc routine.
!
        DO j=SstrR,SendR
#  ifdef VARY_ACBC
          jj=j
#  else
          jj=1                     ! use same location for
                                   ! distx disty and angle_edge
#  endif
          WAVES(ng)%bndwave_z(j,:)=tempz(j,:)
          WAVES(ng)%bndwave_u(j,:)=                                     &
     &                   (tempu(j,:)*cos(WAVES(ng)%angle_edge(jj))+     &
     &                    tempv(j,:)*sin(WAVES(ng)%angle_edge(jj)))
          WAVES(ng)%bndwave_v(j,:)=                                     &
     &                   (-tempu(j,:)*sin(WAVES(ng)%angle_edge(jj))+    &
     &                     tempv(j,:)*cos(WAVES(ng)%angle_edge(jj)))
        END DO
!
        deallocate (DDtheta, k3, Z_bw, bw_ang)
        deallocate (k3_x, k3_y, cg_3, ff, DDf, D1, D2, D3, D4)
        deallocate (DTOT, E3, A3)
        deallocate (Comptemp, Comptempd, Comp_int)
        deallocate (A, AA, zz)
!
# else
!
! ORIGINAL_BOUND  with 4 loops, use for double summation.
!
        DO f1=1,nfreq-1
!
!  Report percent done
!
          IF ((DOMAIN(ng)%SouthWest_Corner(tile).and.                   &
              (LBC(isouth,isAC3d,ng)%acquire.or.                        &
               LBC(iwest,isAC3d,ng)%acquire)) .or.                      &
              (DOMAIN(ng)%NorthEast_Corner(tile).and.                   &
              (LBC(inorth,isAC3d,ng)%acquire.or.                        &
               LBC(ieast,isAC3d,ng)%acquire))) THEN
            IF (MOD(f1*1000/nfreq,100).eq.0) THEN
              write(*,*) 'percent done for computing bndwave is ',      &
     &                    f1*100/nfreq
            END IF
          END IF

          DO f2=f1+1,nfreq
!
! Compute the freq dif
!
            DDf=WAVES(ng)%f(f2)-WAVES(ng)%f(f1)
!
!  Limit the freq diff range from 1/400 to 1/30, that is
!  df=1/400=.0025 to df=1/30=0.033
!
            IF ((DDf.ge.fmin).and.(DDf.le.fmax)) THEN
              fac1=WAVES(ng)%f(f1)*WAVES(ng)%f(f2)
              D1a=g*k(f1)*k(f2)/(8.0_r8*pi**2*fac1)*                    &
     &            ocoshkh(f1)*ocoshkh(f2)
              D3a=((2.0_r8*pi)**4*(fac1)**2)*og2
              D4=-0.5*(-WAVES(ng)%f(f1)*k(f2)**2*ocoshkh(f2)**2+        &
     &                  WAVES(ng)%f(f2)*k(f1)**2*ocoshkh(f1)**2)

              DO n1=1,WAVES(ng)%ndir
                DO n2=1,WAVES(ng)%ndir
!
! Compute the dir diff, wave number, and phase of bound wave.
!
                  DDtheta=abs(WAVES(ng)%dir(n1)-WAVES(ng)%dir(n2))
                  k3=sqrt(k(f1)**2+k(f2)**2-                            &
     &               2.0_r8*k(f1)*k(f2)*COS(DDtheta))
                  Z_bw=WAVES(ng)%phase(f2,n2)-WAVES(ng)%phase(f1,n1)+pi

!
!  Compute the double summation.
!
                  D1=D1a*COS(DDtheta+pi)*MIN(cosh(k3*h0),1.0E10_r8)
                  D2=-g*(DDf)/((g*k3*tanh(k3*h0)-(2.0_r8*pi)**2*DDf**2)*&
     &                   fac1+eps)
                  D3=DDf*(D3a-k(f1)*k(f2)*COS(DDtheta+pi))
                  DTOT=D1+D2*(D3+D4)
!
!  Compute energy of bound wave
!
!   jcw is this **2 for both
                  E3=2.0_r8*DTOT**2*                                    &
     &              WAVES(ng)%SD(f1,n1)*WAVES(ng)%SD(f2,n2)*fac4*fac5**2
                  A3=SQRT(2.0_r8*E3*fac4)
!
!  Compute the direction of the bound wave (this is in radians)
!
                  cff1=k(f2)*sin(WAVES(ng)%dir(n2))-                    &
     &                 k(f1)*sin(WAVES(ng)%dir(n1))
                  cff2=k(f2)*cos(WAVES(ng)%dir(n2))-                    &
     &                 k(f1)*cos(WAVES(ng)%dir(n1))
                  DIR3=atan(cff1/cff2)
                  k3_x=k3*cos(1.5_r8*pi-DIR3)  !DIR3 is nautical coords
                  k3_y=k3*sin(1.5_r8*pi-DIR3)
!
!  Compute this contribution to the bndwave (water level) and associated 
!  currents; currents depend on direction
!  Here we need to include the phase difference along the boundary; this
!  needs to be done; bndwave needs to depend on the boudary (and so does
!  the free surface elevation)
!
#  ifdef VARY_ACBC
!  The implicit time do loop here allows compiler to optimize it.
                  DO j=SstrR,SendR
                    cff2=k3_x*WAVES(ng)%dist_x(j)+                      &
     &                   k3_y*WAVES(ng)%dist_y(j)
                    DO t=1,WAVES(ng)%Insteps
                      tempz(j,t)=tempz(j,t)+A3*COS(-2.0_r8*pi*DDf*      &
     &                                             temp4(t)+Z_bw+cff2)
                    END DO
                  END DO
#  else
                  tempz=tempz+A3*COS(-2.0_r8*pi*DDf*temp4+Z_bw)
#  endif
                END DO
              END DO
            END IF
          END DO
        END DO
!
!  The currents need to be rotated to the roms grid plus they are
!  defined in the rho points
!  Here we assume these waves are shallow progressive water waves
!
        DO j=SstrR,SendR
#  ifdef VARY_ACBC
          WAVES(ng)%bndwave_z(j,:)=WAVES(ng)%bndwave_z(j,:)+            &
     &                             tempz(j,:)
          WAVES(ng)%bndwave_u(j,:)=WAVES(ng)%bndwave_u(j,:)+            &
     &                             tempz(j,:)*oh0*sqrt(g*h0)
          WAVES(ng)%bndwave_v(j,:)=WAVES(ng)%bndwave_v(j,:)+            &
     &                             tempz(j,:)*oh0*sqrt(g*h0)
#  else
          WAVES(ng)%bndwave_z(:)=WAVES(ng)%bndwave_z(:)+tempz(:)
          WAVES(ng)%bndwave_u(:)=WAVES(ng)%bndwave_u(:)+                &
     &                           tempz(:)*oh0*sqrt(g*h0)
          WAVES(ng)%bndwave_v(:)=WAVES(ng)%bndwave_v(:)+                &
     &                           tempz(:)*oh0*sqrt(g*h0)
#  endif
        END DO
# endif
        deallocate (tempz, tempu, tempv)
        deallocate (owaven, ocoshkh)
      END IF

      RETURN
      END SUBROUTINE boundwave_tile

#endif

      END MODULE mod_inwave_swan
