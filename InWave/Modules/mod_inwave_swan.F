#include "cppdefs.h"

#ifdef VARY_ACBC
#  define DOUBLESUM
#  undef  SINGLESUM
#  ifdef SINGLESUM
#    undef  INTCDF
#    define FCDF
#  endif
#  define DHILBERT
#  undef  IHILBERT
#  undef  HSCALE
#endif

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
!  PHASE       Random phase for each frequency- direction component    !
!  AMP         Amplitude of the Fourier components                     !
!  AMPZETA     Amplitude of the free surface elevation for each        !
!              directional component                                   !
!  AMPZETA_TOT Amplitude of the free surface elevation                 !
!  POSITION    positive frequencies of the Fourier Components          !
!  CompFn      Fourier components                                      !
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
         real(r8)            :: dur, idt
         real(r8)            :: fmax, fmin, df
         real(r8)            :: ddir,fp
         real(r8)            :: mainang, hm0gew, dang,int
         real(r8)            :: df_fft
         real(r8)            :: SpecTimeIso(2)
         real(r8)            :: SpecTime(2)
         real(r8)            :: SpecTimedt

         real(r8), pointer :: fsw(:)
         real(r8), pointer :: f(:)
         real(r8), pointer :: dir(:)
         real(r8), pointer :: theta(:)
         real(r8), pointer :: SF(:)

         real(r8), pointer :: SDSW(:,:)
         real(r8), pointer :: SD(:,:)
         real(r8), pointer :: SDD(:)
         real(r8), pointer :: phase(:,:)
# ifdef SINGLESUM
         real(r8), pointer :: amp(:)
# elif defined DOUBLESUM
         real(r8), pointer :: amp(:,:)
# endif
         real(r8), pointer :: zeta(:,:)
# ifdef VARY_ACBC
         real(r8), pointer :: Ampzeta(:,:,:)
         complex(fftkind), dimension(:), allocatable :: Ampzeta_tot
         real(r8), pointer :: bndwave_z(:,:)
         real(r8), pointer :: bndwave_u(:,:)
         real(r8), pointer :: bndwave_v(:,:)
# else
         real(r8), pointer :: Ampzeta(:,:)
         real(r8), pointer :: Ampzeta_tot(:)
         real(r8), pointer :: bndwave_z(:)
         real(r8), pointer :: bndwave_u(:)
         real(r8), pointer :: bndwave_v(:)
# endif
         real(r8), pointer :: position(:)

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
      integer             :: dd, ff, k, offset, i
      real(r8), parameter :: IJspv = 0.0_r8
      real(r8)            :: cff
      real(r8)            :: my_min, my_max
      real(r8), allocatable :: angle(:,:), Iout(:,:), Jout(:,:)
      real(r8), allocatable :: fsw_2d(:,:), theta_2d(:,:)
      real(r8), allocatable :: f_2d(:,:), wd_2d(:,:), SDSW_circle(:,:)
!
!-----------------------------------------------------------------------
! Define the Nyquist frequency and minimum frequency and 
! the number of directional bins
!-----------------------------------------------------------------------
!
      IF (first.eq.1) THEN
!        WAVES(ng)%nfreq=2000
        WAVES(ng)%nfreq=200
        WAVES(ng)%fmax=0.5_r8  !! Tmin = 2 s, this is ~ the Nyquist f
        WAVES(ng)%df=WAVES(ng)%fmax/WAVES(ng)%nfreq
        WAVES(ng)%ndir=ND
        WAVES(ng)%fmin=WAVES(ng)%df
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
        allocate (WAVES(ng)%SD(WAVES(ng)%nfreq,WAVES(ng)%ndir))
        allocate (WAVES(ng)%SDD(WAVES(ng)%ndir))
        DO dd=1,WAVES(ng)%ndir
          DO ff=1,WAVES(ng)%nfreq
            WAVES(ng)%SD(ff,dd)=0.0_r8
          ENDDO
          WAVES(ng)%SDD(dd)=0.0_r8
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
        END DO
      END DO
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
!
!  Set up 2d gridded freq and dir arrays for user defined computation grid.
!
      allocate (Iout(1:WAVES(ng)%nfreq,1:ND))
      allocate (Jout(1:WAVES(ng)%nfreq,1:ND))
      allocate (f_2d(1:WAVES(ng)%nfreq,1:ND))
      allocate (wd_2d(1:WAVES(ng)%nfreq,1:ND))
      DO ff=1,WAVES(ng)%nfreq
        DO dd=1,ND
          Iout(ff,dd)=0.0_r8
          Jout(ff,dd)=0.0_r8
        END DO
      END DO
      DO ff=1,WAVES(ng)%nfreq
        DO dd=1,ND
          f_2d(ff,dd)=WAVES(ng)%f(ff)
          wd_2d(ff,dd)=WAVEG(ng)%wd(dd)
        END DO
      END DO
!
      CALL hindices (ng, 1, WAVES(ng)%nfreqsw,                          &
     &                   1, WAVES(ng)%ndirsw*offset,                    &
     &                   1, WAVES(ng)%nfreqsw,                          &
     &                   1, WAVES(ng)%ndirsw*offset,                    &
     &                   angle, fsw_2d, theta_2d,                       &
     &                   1, WAVES(ng)%nfreq, 1, ND,                     &
     &                   1, WAVES(ng)%nfreq, 1, ND,                     &
     &                   f_2d, wd_2d,                                   &
     &                   Iout, Jout,                                    &
     &                   IJspv, rectangular)
      CALL linterp2d (ng, 1, WAVES(ng)%nfreqsw,                         &
     &                    1, WAVES(ng)%ndirsw*offset,                   &
     &                    fsw_2d, theta_2d, SDSW_circle,                &
     &                    1, WAVES(ng)%nfreq, 1, ND,                    &
     &                    1, WAVES(ng)%nfreq, 1, ND,                    &
     &                    Iout, Jout,                                   &
     &                    f_2d, wd_2d,                                  &
     &                    WAVES(ng)%SD,                                 &
                          my_min, my_max)
      deallocate(angle, Iout, Jout)
      deallocate(fsw_2d, theta_2d, SDSW_circle, f_2d, wd_2d)
# ifdef SINGLESUM
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
      WAVES(ng)%int=0.0_r8
      DO dd=1,WAVES(ng)%ndir
        WAVES(ng)%int=WAVES(ng)%int+WAVES(ng)%SDD(dd)
      ENDDO
# endif
!
      IF (first.eq.1) THEN
!
!  Set bound directional arrays to be equal to the computational dirs.
!
        WAVEB(ng)%ND_bnd=ND
        allocate(WAVEB(ng)%WD_bnd(WAVEB(ng)%ND_bnd))
        DO dd=1,WAVEB(ng)%ND_bnd
          WAVEB(ng)%WD_BND(dd)=WAVEG(ng)%wd(dd)
        END DO
!
!  Determine the time dimensions of duration and num steps.
!  The Ampzeta time series will repeat every 1/df time steps.
!
        WAVES(ng)%dur=1./WAVES(ng)%df
!       WAVES(ng)%idt=4.0_r8
        WAVES(ng)%idt=1.0_r8
!       WAVES(ng)%Insteps=nint(WAVES(ng)%dur/dt(ng))
!       WAVES(ng)%Insteps=nint(WAVES(ng)%dur/WAVES(ng)%idt)+1
        WAVES(ng)%Insteps=nint(WAVES(ng)%dur/WAVES(ng)%idt)
!
!  Allocate and init the computational arrays.
!
# if !defined VARY_ACBC
        allocate(WAVES(ng)%Ampzeta(WAVES(ng)%Insteps,WAVES(ng)%ndir))
# endif
        allocate(WAVES(ng)%Ampzeta_tot(WAVES(ng)%Insteps))
!
!-----------------------------------------------------------------------
! Create frequency and directional arrays for the spectra
!-----------------------------------------------------------------------
!
        allocate (WAVES(ng)%position(WAVES(ng)%nfreq))
        allocate (WAVES(ng)%phase(WAVES(ng)%nfreq,WAVES(ng)%ndir))
# ifdef SINGLESUM
        allocate (WAVES(ng)%amp(WAVES(ng)%nfreq))
# elif defined DOUBLESUM
        allocate (WAVES(ng)%amp(WAVES(ng)%nfreq,WAVES(ng)%ndir))
# endif
      END IF
!
# if !defined VARY_ACBC
      DO dd=1,WAVES(ng)%ndir
        DO i=1,WAVES(ng)%Insteps
          WAVES(ng)%Ampzeta(i,dd)=0.0_r8
        END DO
      END DO

# endif
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
!      integer, dimension(:), allocatable :: seed
      real(r8) :: twopi
# ifdef DISTRIBUTE
      real(r8), allocatable :: wrk(:)
# endif

!      call random_seed(size=1)
!      allocate(seed(1))
!      call system_clock(count=clock)
!      seed = clock + 37 * (/(i-1,i=1,n)/)
!      call random_seed(put=seed)
!      deallocate(seed)

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
      USE mod_inwave_bound

      implicit none

      real(r8), intent(in)  :: Sf(:), f(:)
      real(r8), pointer     :: temp(:)
      real(r8)              :: frep

      allocate(temp(size(Sf)))
      temp=0.0_r8
! CH note: do we need this limitation of the energy?
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

# ifdef DOUBLESUM
      DO dd=1,WAVES(ng)%ndir
        DO ff=1,WAVES(ng)%nfreq
          WAVES(ng)%amp(ff,dd)= sqrt(2.0_r8*WAVES(ng)%SD(ff,dd)*cff1*cff2)
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

#ifdef VARY_ACBC
!
!***********************************************************************
      SUBROUTINE FourierComp (ng, tile, first)
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

      CALL FourierComp_tile (ng, tile, first,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS)
      RETURN
      END SUBROUTINE FourierComp
!
!***********************************************************************
      SUBROUTINE FourierComp_tile (ng, tile, first,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS)
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
      integer                         :: i, j, m, tidx
      integer                         :: ff, dd
      integer                         :: SstrR, SendR, is, jcenter
      integer                         :: MyError, Npts, p1
      real(r8)                        :: cff, cff1, cff2, DIR3, kwx, kwy
      real(r8)                        :: twopi, otwopi, x, cff3
      real(r8)                        :: error, kedge
      real(r8)                        :: F, FDER, tanhkh
      real(r8)                        :: L0, k0, k1, kh, wr
      real(r8)                        :: h0, oh0 
      real(r8)                        :: x_mid, y_mid, dx_min
      real(r8), allocatable           :: k(:), angle_edge(:)
      real(r8), allocatable           :: dist_local(:), dist_x(:), dist_y(:)
      real(r8), allocatable           :: dx_local(:)
# ifdef INTCDF
      real(r8), allocatable           :: SDD_cumul(:)
# elif defined FCDF
      real(r8), allocatable           :: SDD_cumul(:,:)
# endif
# ifdef HSCALE
      real(r8), allocatable           :: Hsig_check(:)
# endif
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: Awrk
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: h_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: zeta_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: xr_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: yr_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: pm_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: pn_local

      real(r8), parameter :: maxErr = 0.1_r8
      real(r8), parameter :: alphad = 0.2_r8

# include "set_bounds.h"

      twopi=2.0_r8*pi
      otwopi=1.0_r8/(2.0_r8*pi)
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
        END DO
      END DO

# ifdef DISTRIBUTE
      CALL mp_gather2d (ng, iNLM, LBi, UBi, LBj, UBj,                   &
     &                  0, r2dvar, 1.0_r8,                              &
     &                  GRID(ng) % rmask,                               &
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
     &                  GRID(ng) % rmask,                               &
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
     &                  GRID(ng) % rmask,                               &
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
     &                  GRID(ng) % rmask,                               &
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
     &                  GRID(ng) % rmask,                               &
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
     &                  GRID(ng) % rmask,                               &
     &                  GRID(ng) % pn, Npts, Awrk)
      CALL mp_bcastf (ng, iNLM, Awrk)
      p1=0
      DO j=0,Mm(ng)+1
        DO i=0,Lm(ng)+1
          p1=p1+1
          pn_local(i,j)=Awrk(p1)
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
        END DO
      END DO
# endif
!
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
        SstrR=JstrR
        SendR=JendR

        allocate(angle_edge(SstrR:SendR))
        DO j=SstrR,SendR
          angle_edge(j)=GRID(ng)%angler(is,j)
        END DO

        x_mid=0.5_r8*(xr_local(is,Mm(ng)+1)+xr_local(is,0))
        y_mid=0.5_r8*(yr_local(is,Mm(ng)+1)+yr_local(is,0))

        allocate(dist_local(0:Mm(ng)+1))
        allocate(dx_local(0:Mm(ng)+1))

        DO j=0,Mm(ng)+1
          dist_local(j)=(yr_local(is,j)-y_mid)
          dx_local(j)=1.0_r8/pn_local(is,j)
        END DO

        jcenter=minloc(abs(dx_local), DIM = 1)
        dx_min=dx_local(jcenter)
        jcenter=minloc(abs(dist_local), DIM = 1)
        h0=h_local(is,jcenter)+zeta_local(is,jcenter)
        oh0=1.0_r8/h0

        IF (DOMAIN(ng)%SouthWest_Corner(tile).or.                       &
     &      DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            write(*,*) 'Mean offshore water depth for AC is ', h0
        END IF
        deallocate(dist_local, dx_local)
        allocate(dist_x(0:Mm(ng)+1),dist_y(0:Mm(ng)+1))        

        DO j=0,Mm(ng)+1
          dist_x(j)=xr_local(is,j)-xr_local(is,jcenter)
          dist_y(j)=yr_local(is,j)-yr_local(is,jcenter)
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
        SstrR=IstrR
        SendR=IendR
        allocate(angle_edge(SstrR:SendR))
        DO i=SstrR,SendR
          angle_edge(i)=GRID(ng)%angler(i,is)
        END DO

        x_mid=0.5_r8*(xr_local(Lm(ng)+1,is)+xr_local(0,is))
        y_mid=0.5_r8*(yr_local(Lm(ng)+1,is)+yr_local(0,is))

        allocate(dist_local(0:Lm(ng)+1))
        allocate(dx_local(0:Lm(ng)+1))

        DO i=0,Lm(ng)+1
          dist_local(i)=(xr_local(i,is)-x_mid)
          dx_local(i)=1.0_r8/pm_local(i,is)
        END DO

        jcenter=minloc(abs(dx_local), DIM = 1)
        dx_min=dx_local(jcenter)
        jcenter=minloc(abs(dist_local), DIM = 1)
        h0=h_local(jcenter,is)+zeta_local(jcenter,is)
        oh0=1.0_r8/h0
        IF (DOMAIN(ng)%SouthWest_Corner(tile).or.                       &
     &      DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            write(*,*) 'Mean offshore water depth for AC is ', h0
        END IF

        deallocate(dist_local, dx_local)
        allocate(dist_x(0:Lm(ng)+1),dist_y(0:Lm(ng)+1))        

        DO i=0,Lm(ng)+1
          dist_x(i)=xr_local(i,is)-xr_local(jcenter,is)
          dist_y(i)=yr_local(i,is)-yr_local(jcenter,is)
        END DO
      END IF

      IF (  ((LBC(iwest,isAC3d,ng)%acquire).and.                        &
     &       (DOMAIN(ng)%Western_Edge(tile)) ).or.                      &
     &      ((LBC(ieast,isAC3d,ng)%acquire).and.                        &
     &       (DOMAIN(ng)%Eastern_Edge(tile))).or.                       &
     &     ((LBC(isouth,isAC3d,ng)%acquire).and.                        &
     &    (DOMAIN(ng)%Southern_Edge(tile))).or.                         &
     &    ((LBC(inorth,isAC3d,ng)%acquire).and.                         &
     &    (DOMAIN(ng)%Northern_Edge(tile)))) THEN
!
!  Compute the wave numbers.
!
        allocate (k(WAVES(ng)%nfreq))
        DO ff=1,WAVES(ng)%nfreq
          L0=g*otwopi*(1.0_r8/WAVES(ng)%f(ff))**2
          k0=twopi/L0
          error=100.0_r8
          wr=twopi*WAVES(ng)%f(ff)
          DO WHILE(error.gt.maxErr)
            kh=k0*h0
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
          k(ff)=k0
        END DO

        IF (first.eq.1) THEN
          allocate(WAVES(ng)%Ampzeta(SstrR:SendR,                       &
     &             WAVES(ng)%Insteps,WAVES(ng)%ndir))
        END IF
!
! Initialize the final array.
!
        DO tidx=1,WAVES(ng)%Insteps
          DO i=SstrR,SendR
            DO dd=1,WAVES(ng)%ndir
              WAVES(ng)%Ampzeta(i,tidx,dd)=0.0_r8
            END DO
          END DO
        END DO
!
!  Choose single or double summation method.
!
# ifdef SINGLESUM
!
!  Create a cumulative distribution function (CDF) of energy.
!
#   ifdef INTCDF
!  Create a CDF integrated over frequencies. 
!  SDD is directional energy spectra (but not x pd).
        allocate(SDD_cumul(WAVES(ng)%ndir))
        SDD_cumul(1)=0.0_r8
        DO dd=2,WAVES(ng)%ndir
! CH note: do we need this energy check?
!          cff2=WAVES(ng)%SDD(dd)/WAVES(ng)%int
!          IF (cff2.gt.0.10_r8) THEN
            SDD_cumul(dd)=SDD_cumul(dd-1)+WAVES(ng)%SDD(dd)
!          ELSE
!            SDD_cumul(dd)=SDD_cumul(dd-1)
!          END IF
        END DO
!  Scale SDD_cumul so it goes from ~0.0 to 1.0
        DO dd=1,WAVES(ng)%ndir
          SDD_cumul(dd)=SDD_cumul(dd)/SDD_cumul(WAVES(ng)%ndir)
        END DO
#   elif defined FCDF
        allocate(SDD_cumul(WAVES(ng)%nfreq,WAVES(ng)%ndir))
        DO ff=1,WAVES(ng)%nfreq
          SDD_cumul(ff,1)=0.0_r8
          DO dd=2,WAVES(ng)%ndir
!            cff2=WAVES(ng)%SD(ff,dd)/WAVES(ng)%SDD(ff)
!            IF (cff2.gt.0.10_r8) THEN
              SDD_cumul(ff,dd)=SDD_cumul(ff,dd-1)+WAVES(ng)%SD(ff,dd)
!            ELSE
!              SDD_cumul(ff,dd)=SDD_cumul(ff,dd-1)
!            END IF
          END DO
        END DO
!  Scale SDD_cumul so it goes from ~0.0 to 1.0
        DO ff=1,WAVES(ng)%nfreq
          DO dd=1,WAVES(ng)%ndir
            SDD_cumul(ff,dd)=SDD_cumul(ff,dd)/SDD_cumul(ff,WAVES(ng)%ndir)
          END DO
        END DO
#   endif
# endif
!
!  Compute the amplitude time series for the boundary AC forcing.
!  Loop over all frequencies
!
        DO ff=1,WAVES(ng)%nfreq
# ifdef SINGLESUM
!
!  Determine which direction to use for this freq.
!
          call random_number(cff1)
#   ifdef INTCDF
          DO dd=1,WAVES(ng)%ndir
            IF (cff1.lt.SDD_cumul(dd)) THEN
              m=dd
              EXIT
            END IF
          END DO
#   elif defined FCDF
          DO dd=1,WAVES(ng)%ndir
            IF (cff1.lt.SDD_cumul(ff,dd)) THEN
              m=dd
              EXIT
            END IF
          END DO
#   endif
# endif
          DO i=SstrR,SendR
!
!  If DOUBLESUM, double summation loop over all directions.
!  If SINGLESUM, check if the random direction assigned above.
!
            DO dd=1,WAVES(ng)%ndir
! CH note: removed the energy check here
# ifdef SINGLESUM
              IF (m.eq.dd) THEN
# endif
                DIR3=1.5_r8*pi-WAVEG(ng)%wd(dd)-angle_edge(i)
                kwx=k(ff)*cos(DIR3)
                kwy=k(ff)*sin(DIR3)
                kedge=kwx*dist_x(i)+kwy*dist_y(i)
                DO tidx=1,WAVES(ng)%Insteps
                  cff=REAL((tidx-1),r8)*WAVES(ng)%idt
# ifdef SINGLESUM
                  cff2=WAVES(ng)%amp(ff)*                               &
     &                      COS(twopi*WAVES(ng)%f(ff)*cff-              &
     &                      kedge+WAVES(ng)%phase(ff,dd))               
# elif defined DOUBLESUM
                  cff2=WAVES(ng)%amp(ff,dd)*                            &
     &                      COS(twopi*WAVES(ng)%f(ff)*cff-              &
     &                      kedge+WAVES(ng)%phase(ff,dd))               
# endif
                  WAVES(ng)%Ampzeta(i,tidx,dd)=                         &
     &                          WAVES(ng)%Ampzeta(i,tidx,dd)+cff2

                END DO
# ifdef SINGLESUM
              END IF
# endif
            END DO
          END DO
        END DO
!
!  Take the Hilbert transform of the amplitude time series 
# ifdef DHILBERT
!  Maintain the directional distribution
!
        DO i=SstrR,SendR
          DO dd=1,WAVES(ng)%ndir
            DO tidx=1,WAVES(ng)%Insteps
              WAVES(ng)%Ampzeta_tot(tidx)=WAVES(ng)%Ampzeta(i,tidx,dd)
            END DO
            call hilbert(WAVES(ng)%Ampzeta_tot,                         &
     &                   size(WAVES(ng)%Ampzeta_tot))
            DO tidx=1,WAVES(ng)%Insteps
              cff1=(ABS(WAVES(ng)%Ampzeta_tot(tidx)))**2
              WAVES(ng)%Ampzeta(i,tidx,dd)=0.5_r8*cff1*g*rho0
            END DO
          END DO
        END DO
# elif defined IHILBERT
!  Integrated 
!
# ifdef HSCALE
        allocate(Hsig_check(SstrR:SendR))
# endif
        DO i=SstrR,SendR
          DO tidx=1,WAVES(ng)%Insteps
            WAVES(ng)%Ampzeta_tot(tidx)=0.0_r8
            DO dd=1,WAVES(ng)%ndir
              WAVES(ng)%Ampzeta_tot(tidx)=WAVES(ng)%Ampzeta_tot(tidx)+  &
     &                                    WAVES(ng)%Ampzeta(i,tidx,dd)
            END DO
          END DO
          call hilbert(WAVES(ng)%Ampzeta_tot,                           &
     &                 size(WAVES(ng)%Ampzeta_tot))
!
!  Partition back out to all directions
!
# ifdef HSCALE
        cff=0.0_r8
# endif
          DO dd=1,WAVES(ng)%ndir
              cff2=WAVES(ng)%SDD(dd)/WAVES(ng)%int
            DO tidx=1,WAVES(ng)%Insteps
              cff1=((ABS(WAVES(ng)%Ampzeta_tot(tidx)))**2)*cff2     
              WAVES(ng)%Ampzeta(i,tidx,dd)=0.5_r8*cff1*g*rho0
# ifdef HSCALE
              cff=cff+cff1
# endif
            END DO
          END DO
# ifdef HSCALE
          Hsig_check(i)=4.0_r8*sqrt(cff/WAVES(ng)%Insteps)
# endif
        END DO
!
# ifdef HSCALE
!
!  Scale energy so the time-averaged Hm0 is equivalent to SWAN spectra 
!
        DO i=SstrR,SendR
          cff=WAVES(ng)%hm0gew/Hsig_check(i)
          DO dd=1,WAVES(ng)%ndir
            DO tidx=1,WAVES(ng)%Insteps
              WAVES(ng)%Ampzeta(i,tidx,dd)=                             &
     &               WAVES(ng)%Ampzeta(i,tidx,dd)*cff
            END DO
          END DO
        END DO
        deallocate(Hsig_check)
# endif
# endif
# ifdef SINGLESUM
        deallocate(SDD_cumul)
# endif
        deallocate(k, angle_edge, dist_x, dist_y)
      END IF
      RETURN
      END SUBROUTINE FourierComp_tile
#else
!
!***********************************************************************
      SUBROUTINE FourierComp (ng, tile, first)
!***********************************************************************
!
      USE mod_scalars

      implicit none

!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, first

!  Local variable declarations.
!
      integer                         :: i, ff, dd
      real(r8)                        :: cff
      real(r8)                        :: twopi 

      twopi=2.0_r8*pi

      IF (  ((LBC(iwest,isAC3d,ng)%acquire).and.                        &
     &       (DOMAIN(ng)%Western_Edge(tile)) ).or.                      &
     &      ((LBC(ieast,isAC3d,ng)%acquire).and.                        &
     &       (DOMAIN(ng)%Eastern_Edge(tile))).or.                       &
     &     ((LBC(isouth,isAC3d,ng)%acquire).and.                        &
     &    (DOMAIN(ng)%Southern_Edge(tile))).or.                         &
     &    ((LBC(inorth,isAC3d,ng)%acquire).and.                         &
     &    (DOMAIN(ng)%Northern_Edge(tile)))) THEN
! COMPUTES THE FOURIER COMPONENTS
      DO i=1,WAVES(ng)%Insteps
        WAVES(ng)%CompFn(i)=0.0_r8
      END DO
      DO dd=1,ND
        DO ff=1,WAVES(ng)%nfreq
          DO i=1,WAVES(ng)%Insteps
            cff=REAL((i-1),r8)*WAVES(ng)%idt
            WAVES(ng)%CompFn(i)=WAVES(ng)%CompFn(i)+                    &
     &                          WAVES(ng)%amp(ff,dd)*                   &
     &                          COS(twopi*WAVES(ng)%f(ff)*cff+          &
     &                          WAVES(ng)%phase(ff,dd))
          END DO
        END DO
      END DO
!
      call hilbert(WAVES(ng)%CompFn,size(WAVES(ng)%CompFn))
!
      WAVES(ng)%Ampzeta_tot(:)=abs(WAVES(ng)%CompFn)

!
!-----------------------------------------------------------------------
! Compute the wave energy, scale for each direction
!-----------------------------------------------------------------------
!
      DO dd=1,ND
        DO i=1,WAVES(ng)%Insteps
          WAVES(ng)%Ampzeta(i,dd)=0.5_r8*g*rho0*                         &
     &                                (WAVES(ng)%Ampzeta_tot(i)**2)*     &
     &                                 WAVES(ng)%SDD(dd)/                &
     &                                 WAVES(ng)%int
        ENDDO
      ENDDO
          END IF

      RETURN
      END SUBROUTINE FourierComp
#endif
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
     &                     WAVES(ng)%Insteps)
      RETURN
      END SUBROUTINE boundwave
!
!***********************************************************************
      SUBROUTINE boundwave_tile (ng, tile, first,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Insteps)
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
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcastf, mp_gather2d
# endif
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, first, Insteps
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

!
!  Local variable declarations.
!
      integer                         :: i, j, t, is, p1, p2, f1, f2
      integer                         :: n1, n2, SstrR, SendR
      integer                         :: Npts, jcenter
      real(r8)                        :: cff, cff1, cff2, cff3, cff4
      real(r8)                        :: twopi, otwopi, og2
      real(r8)                        :: error, fmin, fmax
      real(r8)                        :: fac1, fac2, fac3, fac4
      real(r8)                        :: L0, k0, k1, kh, wr
      real(r8)                        :: F, FDER, tanhkh
      real(r8)                        :: DDf, A3, Z_bw
      real(r8)                        :: DDtheta, k3, k3_x, k3_y
      real(r8)                        :: D1, D2, D3, D4, D1a, D3a, DTOT
      real(r8)                        :: h0, oh0, E3, DIR3, fac5
      real(r8)                        :: x_mid, y_mid
      real(r8), dimension(WAVES(ng)%Insteps) :: junk4
!      real(r8), dimension(0:Lm(ng)+1,WAVES(ng)%Insteps) :: junk1
!      real(r8), dimension(WAVES(ng)%Insteps,0:Lm(ng)+1) :: junk1
!
!      real(r8), dimension(0:Lm(ng)+1) :: dist_x, dist_y
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: Awrk
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: h_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: zeta_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: xr_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: yr_local

      real(r8), allocatable           :: k(:), ocoshkh(:)
      real(r8), allocatable           :: dist_local(:)
      real(r8), allocatable           :: junk1(:,:)
      real(r8), allocatable           :: dist_x(:), dist_y(:)

!# ifdef DISTRIBUTE
!      real(r8), allocatable :: wrk(:)
!# endif
      real(r8), parameter :: maxErr = 0.1_r8
      real(r8), parameter :: eps = 1.0E-10_r8

# include "set_bounds.h"
        twopi=2.0_r8*pi
        otwopi=1.0_r8/(2.0_r8*pi)
        fmin=1.0_r8/400.0_r8
        fmax=1.0_r8/30.0_r8
!
!  Gather the depth and zeta to determine total depth at boundary point.
!
        DO j=0,Mm(ng)+1
          DO i=0,Lm(ng)+1
            h_local(i,j)=0.0_r8
            zeta_local(i,j)=0.0_r8
            xr_local(i,j)=0.0_r8
            yr_local(i,j)=0.0_r8
          END DO
        END DO

# ifdef DISTRIBUTE
        CALL mp_gather2d (ng, iNLM, LBi, UBi, LBj, UBj,                 &
     &                    0, r2dvar, 1.0_r8,                            &
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng)%h, Npts, Awrk)
        CALL mp_bcastf (ng, iNLM, Awrk)
!
        p1=0
        DO j=0,Mm(ng)+1
          DO i=0,Lm(ng)+1
            p1=p1+1
            h_local(i,j)=Awrk(p1)
          END DO
        END DO
!
        CALL mp_gather2d (ng, iNLM, LBi, UBi, LBj, UBj,                 &
     &                    0, r2dvar, 1.0_r8,                            &
     &                    GRID(ng) % rmask,                             &
     &                    OCEAN(ng)%zeta(:,:,1), Npts, Awrk)
        CALL mp_bcastf (ng, iNLM, Awrk)
!
        p1=0
        DO j=0,Mm(ng)+1
          DO i=0,Lm(ng)+1
            p1=p1+1
            zeta_local(i,j)=Awrk(p1)
          END DO
        END DO
!
        CALL mp_gather2d (ng, iNLM, LBi, UBi, LBj, UBj,                 &
     &                    0, r2dvar, 1.0_r8,                            &
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % xr, Npts, Awrk)
        CALL mp_bcastf (ng, iNLM, Awrk)
!
        p1=0
        DO j=0,Mm(ng)+1
          DO i=0,Lm(ng)+1
            p1=p1+1
            xr_local(i,j)=Awrk(p1)
          END DO
        END DO
!
        CALL mp_gather2d (ng, iNLM, LBi, UBi, LBj, UBj,                 &
     &                    0, r2dvar, 1.0_r8,                            &
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % yr, Npts, Awrk)
        CALL mp_bcastf (ng, iNLM, Awrk)
!
        p1=0
        DO j=0,Mm(ng)+1
          DO i=0,Lm(ng)+1
            p1=p1+1
            yr_local(i,j)=Awrk(p1)
          END DO
        END DO
# else
        DO j=0,Mm(ng)+1
          DO i=0,Lm(ng)+1
            h_local(i,j)=GRID(ng)%h(i,j)
            zeta_local(i,j)=OCEAN(ng)%zeta(i,j,1)
            xr_local(i,j)=GRID(ng)%xr(i,j)
            yr_local(i,j)=GRID(ng)%yr(i,j)
          END DO
        END DO
# endif
!
      IF (  ((LBC(iwest,isAC3d,ng)%acquire).and.                        &
     &       (DOMAIN(ng)%Western_Edge(tile)) ).or.                      &
     &      ((LBC(ieast,isAC3d,ng)%acquire).and.                        &
     &       (DOMAIN(ng)%Eastern_Edge(tile))).or.                       &
     &     ((LBC(isouth,isAC3d,ng)%acquire).and.                        &
     &    (DOMAIN(ng)%Southern_Edge(tile))).or.                         &
     &    ((LBC(inorth,isAC3d,ng)%acquire).and.                         &
     &    (DOMAIN(ng)%Northern_Edge(tile)))) THEN

        IF (((LBC(iwest,isAC3d,ng)%acquire).and.                       &
     &      (DOMAIN(ng)%Western_Edge(tile))).or.                       &
     &      ((LBC(ieast,isAC3d,ng)%acquire).and.                       &
     &      (DOMAIN(ng)%Eastern_Edge(tile)))) THEN
!
!  Find the boundary cell closest to the center of the boundary
!
          IF (LBC(iwest,isAC3d,ng)%acquire) THEN
            IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
              write(*,*) 'Computing western bound wave'
            END IF
            is=0
          ELSE
            IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
              write(*,*) 'Computing eastern bound wave'
            END IF
            is=Lm(ng)+1
          END IF
          SstrR=JstrR
          SendR=JendR

          x_mid=0.5_r8*(xr_local(is,Mm(ng)+1)+xr_local(is,0))
          y_mid=0.5_r8*(yr_local(is,Mm(ng)+1)+yr_local(is,0))

          allocate(dist_local(0:Mm(ng)+1))

          DO j=0,Mm(ng)+1
            dist_local(j)=(yr_local(is,j)-y_mid)
          END DO

          jcenter=minloc(abs(dist_local), DIM = 1)
          h0=h_local(is,jcenter)+zeta_local(is,jcenter)
          oh0=1.0_r8/h0

          IF (DOMAIN(ng)%SouthWest_Corner(tile).or.                     &
     &        DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            write(*,*) 'Mean offshore water depth for bndwave is ', h0
          END IF
          allocate(dist_x(0:Mm(ng)+1),dist_y(0:Mm(ng)+1))        

          DO j=0,Mm(ng)+1
            dist_x(j)=xr_local(is,j)-xr_local(is,jcenter)
            dist_y(j)=yr_local(is,j)-yr_local(is,jcenter)
          END DO
        END IF
!
        IF (((LBC(isouth,isAC3d,ng)%acquire).and.                       &
     &      (DOMAIN(ng)%Southern_Edge(tile))).or.                       &
     &      ((LBC(inorth,isAC3d,ng)%acquire).and.                       &
     &      (DOMAIN(ng)%Northern_Edge(tile)))) THEN
!
!  Find the boundary cell closest to the center of the boundary
!
          IF (LBC(isouth,isAC3d,ng)%acquire) THEN
            IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
               write(*,*) 'Computing Southern bound wave'
            END IF
            is=0
          ELSE
            IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
              write(*,*) 'Computing Northern bound wave'
            END IF
            is=Mm(ng)+1
          END IF
          SstrR=IstrR
          SendR=IendR

          x_mid=0.5_r8*(xr_local(Lm(ng)+1,is)+xr_local(0,is))
          y_mid=0.5_r8*(yr_local(Lm(ng)+1,is)+yr_local(0,is))

          allocate(dist_local(0:Lm(ng)+1))

          DO i=0,Lm(ng)+1
            dist_local(i)=(xr_local(i,is)-x_mid)
          END DO

          jcenter=minloc(abs(dist_local), DIM = 1)
          h0=h_local(jcenter,is)+zeta_local(jcenter,is)
          oh0=1.0_r8/h0

          IF (DOMAIN(ng)%SouthWest_Corner(tile).or.                     &
     &        DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            write(*,*) 'Mean offshore water depth for bndwave is ', h0
          END IF        
          allocate(dist_x(0:Lm(ng)+1),dist_y(0:Lm(ng)+1))        

          DO i=0,Lm(ng)+1
            dist_x(i)=xr_local(i,is)-xr_local(jcenter,is)
            dist_y(i)=yr_local(i,is)-yr_local(jcenter,is)
          END DO
        END IF
!
!  Allocate boundwave and local wavenumber arrays.
!
# ifdef VARY_ACBC
        IF (first.eq.1) THEN
          allocate(WAVES(ng)%bndwave_z(SstrR:SendR,WAVES(ng)%Insteps))
          allocate(WAVES(ng)%bndwave_u(SstrR:SendR,WAVES(ng)%Insteps))
          allocate(WAVES(ng)%bndwave_v(SstrR:SendR,WAVES(ng)%Insteps))
        END IF
        DO i=SstrR,SendR
          DO t=1,WAVES(ng)%Insteps
            WAVES(ng)%bndwave_z(i,t)=0.0_r8
            WAVES(ng)%bndwave_u(i,t)=0.0_r8
            WAVES(ng)%bndwave_v(i,t)=0.0_r8
          END DO
        END DO
# else
        IF (first.eq.1) THEN
          allocate(WAVES(ng)%bndwave_z(WAVES(ng)%Insteps))
          allocate(WAVES(ng)%bndwave_u(WAVES(ng)%Insteps))
          allocate(WAVES(ng)%bndwave_v(WAVES(ng)%Insteps))
        END IF
        DO t=1,WAVES(ng)%Insteps
          WAVES(ng)%bndwave_z(t)=0.0_r8
          WAVES(ng)%bndwave_u(t)=0.0_r8
          WAVES(ng)%bndwave_v(t)=0.0_r8
        END DO
# endif
!
        allocate (k(WAVES(ng)%nfreq))
!
! Compute the wave number of each freq bin
!
        DO i=1,WAVES(ng)%nfreq
          L0=g*otwopi*(1.0_r8/WAVES(ng)%f(i))**2.0_r8
          k0=twopi/L0
          error=100.0_r8
          wr=twopi*WAVES(ng)%f(i)
          DO WHILE(error.gt.maxErr)
            kh=k0*h0
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
          k(i)=k0
        END DO
!
!make some temp arrays
!
        allocate (ocoshkh(WAVES(ng)%nfreq))
        DO i=1,WAVES(ng)%nfreq
          fac2=cosh(k(i)*h0)
          ocoshkh(i)=1.0_r8/fac2
        END DO
!
! Compute the energy transfer for each pair of frequency bins.
!
        allocate (junk1(SstrR:SendR,WAVES(ng)%Insteps))
        DO t=1,WAVES(ng)%Insteps
          junk4(t)=REAL(t-1,r8)*WAVES(ng)%idt
          DO j=SstrR,SendR
            junk1(j,t)=0.0_r8
          END DO
        END DO
        fac4=WAVES(ng)%df
        fac5=WAVEG(ng)%pd
        og2=1.0_r8/(g**2)
!  Do the double summation loop.
!
        DO f1=1,WAVES(ng)%nfreq-1
!
!  Report percent done
!
          IF (DOMAIN(ng)%SouthWest_Corner(tile).or.                     &
     &        DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            cff2=REAL(f1,r8)/REAL(WAVES(ng)%nfreq,r8)
            IF (MOD(INT(cff2*100.0_r8),10).eq.0) THEN
              write(*,*) 'percent done for computing bndwave is ',      &
     &                    INT(cff2*100.0_r8)
            END IF
          END IF

          DO f2=f1+1,WAVES(ng)%nfreq
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
              D4=-0.5*(-WAVES(ng)%f(f1)*k(f2)**2*ocoshkh(f2)**2+               &
     &                  WAVES(ng)%f(f2)*k(f1)**2*ocoshkh(f1)**2)

              DO n1=1,WAVES(ng)%ndir
                DO n2=1,WAVES(ng)%ndir
!
! Compute the dir diff, wave number, and phase of bound wave.
!
                  DDtheta=abs(WAVEG(ng)%wd(n1)-WAVEG(ng)%wd(n2))
                  k3=sqrt(k(f1)**2+k(f2)**2-                            &
     &               2.0_r8*k(f1)*k(f2)*COS(DDtheta))
                  Z_bw=Waves(ng)%phase(f2,n2)-Waves(ng)%phase(f1,n1)+pi

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
                  cff1=k(f2)*sin(WAVEG(ng)%wd(n2))-                     &
     &                 k(f1)*sin(WAVEG(ng)%wd(n1))
                  cff2=k(f2)*cos(WAVEG(ng)%wd(n2))-                     &
     &                 k(f1)*cos(WAVEG(ng)%wd(n1))
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
# ifdef VARY_ACBC
!  The implicit time do loop here allows compiler to optimize it.
                  DO j=SstrR,SendR
                    cff2=k3_x*dist_x(j)+k3_y*dist_y(j)
!                     junk1(j,:)=junk1(j,:)+A3*COS(-2.0_r8*pi*DDf*junk4(:)+Z_bw+cff2)
                    DO t=1,WAVES(ng)%Insteps
                      junk1(j,t)=junk1(j,t)+A3*COS(-2.0_r8*pi*DDf*junk4(t)+Z_bw+cff2)
                    END DO
                  END DO
# else
                  junk1=junk1+A3*COS(-2.0_r8*pi*DDf*junk4+Z_bw)
!                 junk1=junk1+A3*COS(-2.0_r8*pi*DDf*WAVES(ng)%idt*REAL(t,r8)+Z_bw)
# endif
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
# ifdef VARY_ACBC
          WAVES(ng)%bndwave_z(j,:)=WAVES(ng)%bndwave_z(j,:)+            &
     &                             junk1(j,:)
          WAVES(ng)%bndwave_u(j,:)=WAVES(ng)%bndwave_u(j,:)+            &
     &                             junk1(J,:)*oh0*sqrt(g*h0)
          WAVES(ng)%bndwave_v(j,:)=WAVES(ng)%bndwave_v(j,:)+            &
     &                             junk1(j,:)*oh0*sqrt(g*h0)
# else
          WAVES(ng)%bndwave_z(:)=WAVES(ng)%bndwave_z(:)+junk1(:)
          WAVES(ng)%bndwave_u(:)=WAVES(ng)%bndwave_u(:)+                &
     &                           junk1(:)*oh0*sqrt(g*h0)
          WAVES(ng)%bndwave_v(:)=WAVES(ng)%bndwave_v(:)+                &
     &                           junk1(:)*oh0*sqrt(g*h0)
# endif
        END DO
!
        deallocate(dist_x,dist_y)
        deallocate (k, ocoshkh)
        deallocate (junk1)

      END IF

      RETURN

      END SUBROUTINE boundwave_tile

#endif

      END MODULE mod_inwave_swan
