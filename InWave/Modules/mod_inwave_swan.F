#include "cppdefs.h"

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
!  FNYQ        Nyquist frequency                                       !
!  FMAX        Maximum frequency to consider when recontructing the    !
!              free surface signal                                     !
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
!  FP          Peak frequency                                          !
!  DF_FFT      Frequency increment for the positive Fourier Components !
!  FSW         Frequency array                                         !
!  F           Frequency array                                         !
!  DIR         Directional array read from swan                        !
!  THETA       Directional array                                       !
!  SF          Spectral density function read from swan                !
!  SDSW        Directional spectral density function read from swan    !
!  SD          Directional spectral density function                   !
!  SDD         Spectral density function integrated over the frequencies!
!  PHASE       Random phase for each frequency- direction component    !
!  AMP         Amplitude of the Fourier components                     !
!  POS_F       Positive frequencies of the Fourier Components          !
!  ZETA        Free surface elevation for each directional bin         !
!  AMPZETA     Amplitude of the free surface elevation for the         !
!              free surface elevation                                  !
!  AMPZETA_TOT Amplitude of the free surface elevation for the         !
!              free surface elevation                                  !
!  POSITION    positive frequencies of the Fourier Components          !
!  CompFn      Fourier components                                      !
!  Comptemp    Fourier components                                      !
!=======================================================================
!
!
!=======================================================================
!                                                                      !
!  CONTAINS THE FOLLOWING SUBROUTINES:                                 !
!                                                                      !
!     inwave_swan: this is the main driver of boundary conditions      !
!                  computations from swan output spectra               !
!     array_gen  : this subroutine generates the arrays necesary to    !
!                  compute the free surface elevation time series      !
!                  from the directional wave spectra derived from swan !
!     random_phase  : it assignes random phases to each                !
!                     frequency component                              !
!     tpdcalc    : This subroutine computes the respresentative period !
!     amplitudes : This subroutine computes the amplitude for the      !
!                  Fourier components for each frequency in the spectra!
!     FourierComp : Computes the Fourier components                    !
!=======================================================================
!

#ifdef INWAVE_SWAN_COUPLING

        USE mod_kinds
        USE mod_inwave_params
        USE math_tools
        USE mod_iounits

        implicit none

        TYPE SHORT_WAVE

         integer             :: ndir,nfreq, mainangi
         integer             :: ndirsw,nfreqsw
         integer             :: Insteps, Swancircle
         real(r8)            :: dur,dt
         real(r8)            :: fnyq,fmax,fmin,df
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
         real(r8), pointer :: amp(:,:)
         real(r8), pointer :: pos_f(:)
         real(r8), pointer :: zeta(:,:)
# ifdef VARY_ACBC
         real(r8), pointer :: Ampzeta(:,:,:)
# else
         real(r8), pointer :: Ampzeta(:,:)
         real(r8), pointer :: Ampzeta_tot(:)
# endif
         real(r8), pointer :: position(:)

# ifdef VARY_ACBC
         complex(fftkind),dimension(:,:),allocatable :: CompFn
         complex(fftkind),dimension(:),allocatable   :: Ampzeta_tot
         real(r8), pointer :: bndwave_z(:,:)
         real(r8), pointer :: bndwave_u(:,:)
         real(r8), pointer :: bndwave_v(:,:)
# else
         complex(fftkind),dimension(:),allocatable   :: CompFn
         real(r8), pointer :: bndwave_z(:)
         real(r8), pointer :: bndwave_u(:)
         real(r8), pointer :: bndwave_v(:)
# endif

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

!  Local variable declarations.

!      integer :: i, j, dir_index, numavg, numavgh, start, tidx, dindex
!      integer, allocatable :: peaks(:)
!      real(r8) :: cff
!      real(r8), allocatable :: zeta_filt(:), zeta_filt2(:)
!
!-----------------------------------------------------------------------
! Read swan output spectra
!-----------------------------------------------------------------------
!
      call swan_reader (ng, first, tile)
!
!-----------------------------------------------------------------------
! Generate the arrays for the frequency components and Fourier components
!-----------------------------------------------------------------------
!
      call array_gen (ng, first)
!
!-----------------------------------------------------------------------
! Generate random phases for each frequency component
!-----------------------------------------------------------------------
!
      call random_phase (ng)
!
!-----------------------------------------------------------------------
! Compute the amplitudes for each frequency component
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
! Compute the Hilbert transform
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
      integer             :: i, t, d, j, k, p, n_pos_f, offset
      real(r8), parameter :: IJspv = 0.0_r8
      real(r8)            :: Fmin, Fmax, cff
      real(r8)            :: my_min, my_max
      real(r8), allocatable :: angle(:,:), Iout(:,:), Jout(:,:)
      real(r8), allocatable :: fsw_2d(:,:), theta_2d(:,:)
      real(r8), allocatable :: f_2d(:,:), wd_2d(:,:), SDSW_circle(:,:)
!
!-----------------------------------------------------------------------
! Define the Nyquist frequency, the maximum and minimum frequencies and 
! the number of directional bins
!-----------------------------------------------------------------------
!
      IF (first.eq.1) THEN
        WAVES(ng)%nfreq=100   !jcw 200 ! based on Van Dongeren et al. (2004)
        WAVES(ng)%fmax=0.5_r8  !! Tmin=2 sec
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
        allocate (WAVES(ng)%SD(WAVES(ng)%nfreq,ND))
        allocate (WAVES(ng)%SDD(ND))
        DO j=1,WAVES(ng)%ndir
          DO i=1,WAVES(ng)%nfreq
            WAVES(ng)%SD(i,j)=0.0_r8
          ENDDO
          WAVES(ng)%SDD(j)=0.0_r8
        ENDDO
!
        DO i=1,WAVES(ng)%nfreq
          WAVES(ng)%f(i)=REAL(i-1,r8)*WAVES(ng)%df+WAVES(ng)%fmin
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
      DO i=1,WAVES(ng)%nfreqsw
        DO j=1,WAVES(ng)%ndirsw*offset
          angle(i,j)=0.0_r8
        END DO
      END DO
      DO i=1,WAVES(ng)%nfreqsw
        DO j=1,WAVES(ng)%ndirsw*offset
          fsw_2d(i,j)=WAVES(ng)%fsw(i)
          IF (j.le.WAVES(ng)%ndirsw) THEN
            k=j
            theta_2d(i,j)=WAVES(ng)%theta(k)-360.0_r8*pi/180.0_r8
          ELSEIF (j.le.WAVES(ng)%ndirsw*2) THEN
            k=j-WAVES(ng)%ndirsw
            theta_2d(i,j)=WAVES(ng)%theta(k)
          ELSE
            k=j-WAVES(ng)%ndirsw*2
            theta_2d(i,j)=WAVES(ng)%theta(k)+360.0_r8*pi/180.0_r8
          END IF
          SDSW_circle(i,j)=WAVES(ng)%SDSW(i,k)
        END DO
      END DO
!
!  Set up 2d gridded freq and dir arrays for user defined computation grid.
!
      allocate (Iout(1:WAVES(ng)%nfreq,1:ND))
      allocate (Jout(1:WAVES(ng)%nfreq,1:ND))
      allocate (f_2d(1:WAVES(ng)%nfreq,1:ND))
      allocate (wd_2d(1:WAVES(ng)%nfreq,1:ND))
      DO i=1,WAVES(ng)%nfreq
        DO j=1,ND
          Iout(i,j)=0.0_r8
          Jout(i,j)=0.0_r8
        END DO
      END DO
      DO i=1,WAVES(ng)%nfreq
        DO j=1,ND
          f_2d(i,j)=WAVES(ng)%f(i)
          wd_2d(i,j)=WAVEG(ng)%wd(j)
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
!
!  Sum up all wave E for each dir.
!
      DO j=1,WAVES(ng)%ndir
        cff=0.0_r8
        DO i=1,WAVES(ng)%nfreq
          cff=cff+WAVES(ng)%SD(i,j)
        ENDDO
        WAVES(ng)%SDD(j)=cff
      ENDDO
!
!  Sum up E over all freqs
!
      WAVES(ng)%int=0.0_r8
      DO i=1,WAVES(ng)%ndir
        WAVES(ng)%int=WAVES(ng)%int+WAVES(ng)%SDD(i)
      ENDDO
!
      IF (first.eq.1) THEN
!
!  Set bound directional arrays to be equal to the computational dirs.
!
        WAVEB(ng)%ND_bnd=ND
        allocate(WAVEB(ng)%WD_bnd(ND))
        DO i=1,ND
          WAVEB(ng)%WD_BND(i)=WAVEG(ng)%wd(i)
        END DO
!
!  Determine the time dimensions of duration and num steps.
!  The Ampzeta time series will repeat every 1/df time steps.
!
         WAVES(ng)%dur=1./WAVES(ng)%df
!        WAVES(ng)%dur=1800.0_r8
!        WAVES(ng)%dur=900.0_r8
!  Force the number of steps to be 2000 to Hilbert.
        WAVES(ng)%Insteps=nint(WAVES(ng)%dur/dt(ng))
!       WAVES(ng)%Insteps=nint(WAVES(ng)%dur/1.0_r8)
!
!  Allocate and init the computational arrays.
!
# if !defined VARY_ACBC
        allocate(WAVES(ng)%CompFn(WAVES(ng)%Insteps))
        allocate(WAVES(ng)%Ampzeta(WAVES(ng)%Insteps,WAVES(ng)%ndir))
# endif
        allocate(WAVES(ng)%Ampzeta_tot(WAVES(ng)%Insteps))
!       allocate(WAVES(ng)%bndwave(WAVES(ng)%Insteps))
!
!-----------------------------------------------------------------------
! Create frequency and directional arrays for the spectra
!-----------------------------------------------------------------------
!
        allocate (WAVES(ng)%position(WAVES(ng)%nfreq))
        allocate (WAVES(ng)%phase(WAVES(ng)%nfreq,WAVES(ng)%ndir))
        allocate (WAVES(ng)%amp(WAVES(ng)%nfreq,WAVES(ng)%ndir))
      END IF
!
# if !defined VARY_ACBC
      DO j=1,WAVES(ng)%ndir
        DO i=1,WAVES(ng)%Insteps
          WAVES(ng)%Ampzeta(i,j)=0.0_r8
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

      integer :: i, j, k, Npts, MyError
      real(r8) :: twopi
# ifdef DISTRIBUTE
      real(r8), allocatable :: wrk(:)
# endif

      call random_number(WAVES(ng)%phase)

      twopi=2.0_r8*pi

# ifdef DISTRIBUTE
      IF (Master) THEN
# endif
        DO i=1,WAVES(ng)%nfreq
          DO j=1,WAVES(ng)%ndir
            WAVES(ng)%phase(i,j)=WAVES(ng)%phase(i,j)*twopi
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
        DO i=1,WAVES(ng)%nfreq
          DO j=1,WAVES(ng)%ndir
            k=k+1
            wrk(k)=WAVES(ng)%phase(i,j)
          END DO
        END DO
      END IF
      CALL MPI_BCAST(wrk, Npts, MP_FLOAT, 0,                            &
     &               OCN_COMM_WORLD, MyError)
      k=0
      DO i=1,WAVES(ng)%nfreq
        DO j=1,WAVES(ng)%ndir
          k=k+1
          WAVES(ng)%phase(i,j)=wrk(k)
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
      integer                         :: i, j
      integer                         :: dir_index
      real(r8)                        :: cff1, cff2
!
!-----------------------------------------------------------------------
! Compute the amplitude for the Fourier components
! For each frequency in the spectra there is one amplitude
!-----------------------------------------------------------------------
!
      cff1=WAVES(ng)%df
      cff2=WAVEG(ng)%pd
      DO j=1,WAVES(ng)%ndir
        DO i=1,WAVES(ng)%nfreq
          WAVES(ng)%amp(i,j)= sqrt(2.0_r8*WAVES(ng)%SD(i,j)*cff1*cff2)
        END DO
      END DO

! Assing a position in the spectral frequency array to each Fourier component

      DO i=1,WAVES(ng)%nfreq
        WAVES(ng)%position(i)=i
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
      integer                         :: i, j, m, dindex, tidx
      integer                         :: SstrR, SendR, is, jcenter
      integer                         :: MyError, Npts, p1
      real(r8)                        :: cff, cff1, cff2, DIR3, k3_x, k3_y
      real(r8)                        :: twopi, otwopi, x
      real(r8)                        :: error, kedge
      real(r8)                        :: F, FDER, tanhkh
      real(r8)                        :: L0, k0, k1, kh, wr
      real(r8)                        :: h0, oh0, fmin, fmax
      real(r8)                        :: x_mid, y_mid, dx_min
      real(r8), allocatable           :: k(:), angle_edge(:)
      real(r8), allocatable           :: dist_local(:), dist_x(:), dist_y(:)
      real(r8), allocatable           :: dx_local(:)
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

      fmin=1.0_r8/15.0_r8
      fmax=1.0_r8/3.0_r8
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
# ifdef DISTRIBUTE
          IF (Master) THEN
# endif
          write(*,*) 'Computing western AC forcing'
# ifdef DISTRIBUTE
          END IF
# endif
           is=0
        ELSE
# ifdef DISTRIBUTE
          IF (Master) THEN
# endif
          write(*,*) 'Computing eastern AC forcing'
# ifdef DISTRIBUTE
          END IF
# endif
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
# ifdef DISTRIBUTE
          IF (Master) THEN
# endif
            write(*,*) 'Mean offshore water depth for AC is ', h0
# ifdef DISTRIBUTE
          END IF
# endif
        
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
# ifdef DISTRIBUTE
          IF (Master) THEN
# endif
           write(*,*) 'Computing Southern AC forcing'
# ifdef DISTRIBUTE
          END IF
# endif
           is=0
        ELSE
# ifdef DISTRIBUTE
          IF (Master) THEN
# endif
          write(*,*) 'Computing Northern AC forcing'
# ifdef DISTRIBUTE
          END IF
# endif
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
# ifdef DISTRIBUTE
        IF (Master) THEN
# endif
          write(*,*) 'Mean offshore water depth for AC is ', h0
# ifdef DISTRIBUTE
        END IF
# endif

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
      DO i=1,WAVES(ng)%nfreq
        L0=g*otwopi*(1.0_r8/WAVES(ng)%f(i))**2
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
!        write(*,*) 'wave length ', i,2.0_r8*pi/k(i),1.0_r8/WAVES(ng)%f(i)
      END DO
!
      IF (first.eq.1) THEN
        allocate(WAVES(ng)%CompFn(SstrR:SendR,WAVES(ng)%Insteps))
        allocate(WAVES(ng)%Ampzeta(SstrR:SendR,WAVES(ng)%Insteps,WAVES(ng)%ndir))
      END IF
!
! COMPUTES THE FOURIER COMPONENTS
!
      DO tidx=1,WAVES(ng)%Insteps
        DO i=SstrR,SendR
          WAVES(ng)%CompFn(i,tidx)=0.0_r8
          DO dindex=1,WAVES(ng)%ndir
            WAVES(ng)%Ampzeta(i,tidx,dindex)=0.0_r8
          END DO
        END DO
      END DO
      DO dindex=1,ND
!
! check if this dir has at least 10 percent of total energy
!
        cff2=WAVES(ng)%SDD(dindex)/WAVES(ng)%int
        IF (cff2.gt.0.10_r8) THEN
          DO i=SstrR,SendR
            DO tidx=1,WAVES(ng)%Insteps
              cff2=0.0_r8
              cff=REAL((tidx-1),r8)*dt(ng)
              DO j=1,WAVES(ng)%nfreq
                IF ((WAVES(ng)%f(j).ge.fmin).and.(WAVES(ng)%f(j).le.fmax)) THEN
                  DIR3=1.5_r8*pi-WAVEG(ng)%wd(dindex)-angle_edge(i)

!
!  Here we use an avg wavelength for each group, based on min grid size.
!
                  if (twopi/k(j).ge.(8.0_r8*dx_min)) then
                    k(j)=twopi/(8.0_r8*dx_min)
                  elseif (twopi/k(j).ge.(6.0_r8*dx_min)) then
                    k(j)=twopi/(6.0_r8*dx_min)
                  elseif (twopi/k(j).ge.(4.0_r8*dx_min)) then
                    k(j)=twopi/(4.0_r8*dx_min)
                  elseif (twopi/k(j).ge.(2.0_r8*dx_min)) then
!                    k(j)=twopi/(2.0_r8*dx_min)
!                  else
!!                  k(j)=twopi/(2.0_r8*dx_min)
                   k(j)=0.0_r8
                  end if
!
                  k3_x=k(j)*cos(DIR3)
                  k3_y=k(j)*sin(DIR3)
                  kedge=k3_x*dist_x(i)+k3_y*dist_y(i)
                  cff2=cff2+WAVES(ng)%amp(j,dindex)*                    &
     &                      COS(twopi*WAVES(ng)%f(j)*cff+               &
     &                      WAVES(ng)%phase(j,dindex)-                  &
     &                      kedge)
                END IF
              END DO
              WAVES(ng)%CompFn(i,tidx)=WAVES(ng)%CompFn(i,tidx)+cff2            
            END DO
          END DO
        END IF
      END DO
!
!  Sum over all the directions and then partition out.
!
        DO i=SstrR,SendR
          DO tidx=1,WAVES(ng)%Insteps
            WAVES(ng)%Ampzeta_tot(tidx)=WAVES(ng)%CompFn(i,tidx)
          END DO
          call hilbert(WAVES(ng)%Ampzeta_tot,                           &
     &                 size(WAVES(ng)%Ampzeta_tot))
          DO tidx=1,WAVES(ng)%Insteps
            WAVES(ng)%CompFn(i,tidx)=WAVES(ng)%Ampzeta_tot(tidx)
          END DO
        END DO

        DO i=SstrR,SendR
          DO tidx=2,WAVES(ng)%Insteps-1
            WAVES(ng)%CompFn(i,tidx)=0.3333_r8*(WAVES(ng)%CompFn(i,tidx-1)+ &
     &                                           WAVES(ng)%CompFn(i,tidx)+ &
     &                                           WAVES(ng)%CompFn(i,tidx+1))
          END DO
        END DO

        DO i=SstrR,SendR
          DO dindex=1,ND
            DO tidx=1,WAVES(ng)%Insteps
!              WAVES(ng)%Ampzeta(i,tidx,dindex)=0.5_r8*g*rho0*           &
!     &                    ((ABS(WAVES(ng)%Ampzeta_tot(tidx)))**2)*      &
!     &                                   WAVES(ng)%SDD(dindex)/         &
!     &                                   WAVES(ng)%int
              WAVES(ng)%Ampzeta(i,tidx,dindex)=0.5_r8*g*rho0*           &
     &                    ((ABS(WAVES(ng)%CompFn(i,tidx)))**2)*      &
     &                                   WAVES(ng)%SDD(dindex)/         &
     &                                   WAVES(ng)%int
            END DO
          END DO
        END DO
        deallocate (k, angle_edge, dist_x, dist_y)
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
      integer                         :: i, j, dindex
      real(r8)                        :: cff
      real(r8)                        :: twopi, fmin, fmax

      twopi=2.0_r8*pi
      fmin=1.0_r8/400.0_r8
      fmax=1.0_r8/3.0_r8

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
      DO dindex=1,ND
        DO j=1,WAVES(ng)%nfreq
          IF ((WAVES(ng)%f(j).ge.fmin).and.(WAVES(ng)%f(j).le.fmax)) THEN
          DO i=1,WAVES(ng)%Insteps
            cff=REAL((i-1),r8)*dt(ng)
            WAVES(ng)%CompFn(i)=WAVES(ng)%CompFn(i)+                    &
     &                          WAVES(ng)%amp(j,dindex)*                &
     &                          COS(twopi*WAVES(ng)%f(j)*cff+           &
     &                          WAVES(ng)%phase(j,dindex))
          END DO
          END IF
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
      DO dindex=1,ND
        DO i=1,WAVES(ng)%Insteps
          WAVES(ng)%Ampzeta(i,dindex)=0.5_r8*g*rho0*                     &
     &                                (WAVES(ng)%Ampzeta_tot(i)**2)*     &
     &                                 WAVES(ng)%SDD(dindex)/            &
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
      integer                         :: i, j, t, is, p1, p2, f1, f2, dum
      integer                         :: n1, n2, SstrR, SendR
      integer                         :: Npts, MyError, jcenter
      real(r8)                        :: cff, cff1, cff2, cff3, cff4
      real(r8)                        :: twopi, otwopi
      real(r8)                        :: error, Tr_min, fmin, fmax
      real(r8)                        :: fac1, fac2, fac3, fac4
      real(r8)                        :: L0, k0, k1, kh, wr
      real(r8)                        :: F, FDER, tanhkh
      real(r8)                        :: DDf, A3, Z_bw
      real(r8)                        :: DDtheta, k3, k3_x, k3_y
      real(r8)                        :: D1, D2, D3, D4, D1a, D3a, DTOT
      real(r8)                        :: h0, oh0, E3, DIR3, fac5
      real(r8)                        :: x_mid, y_mid
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: Awrk
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: h_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: zeta_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: xr_local
      real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: yr_local

!      real(r8), allocatable           :: E1d(:), D1d(:), P1d(:)
      real(r8), allocatable           :: k(:)
      real(r8), allocatable           :: dist_local(:), dist_x(:), dist_y(:)
      real(r8), allocatable           :: junk1(:), junk2(:), junk3(:), junk4(:)
      real(r8), allocatable           :: junk(:,:)

# ifdef DISTRIBUTE
      real(r8), allocatable :: wrk(:)
# endif
      real(r8), parameter :: maxErr = 0.1_r8
      real(r8), parameter :: eps = 1.0E-10_r8

# include "set_bounds.h"



        twopi=2.0_r8*pi
        otwopi=1.0_r8/(2.0_r8*pi)
        Tr_min=1.0_r8
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
# ifdef DISTRIBUTE
           IF (Master) THEN
# endif
            write(*,*) 'Computing western bound wave'
# ifdef DISTRIBUTE
           END IF
# endif
             is=0
          ELSE
# ifdef DISTRIBUTE
           IF (Master) THEN
# endif
            write(*,*) 'Computing eastern bound wave'
# ifdef DISTRIBUTE
           END IF
# endif
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
# ifdef DISTRIBUTE
        IF (Master) THEN
# endif
          write(*,*) 'Mean offshore water depth for bndwave is ', h0
# ifdef DISTRIBUTE
        END IF
# endif
        
          deallocate(dist_local)

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
# ifdef DISTRIBUTE
           IF (Master) THEN
# endif
             write(*,*) 'Computing Southern bound wave'
# ifdef DISTRIBUTE
           END IF
# endif
             is=0
          ELSE
# ifdef DISTRIBUTE
           IF (Master) THEN
# endif
            write(*,*) 'Computing Northern bound wave'
# ifdef DISTRIBUTE
           END IF
# endif
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
# ifdef DISTRIBUTE
        IF (Master) THEN
# endif
          write(*,*) 'Mean offshore water depth for bndwave is ', h0
# ifdef DISTRIBUTE
        END IF
# endif
        
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
! Compute the wave number or each freq bin
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
! Compute the wave number or each freq bin
!
! Compute the energy transfer for each pair of frequency bins.
!
        allocate( junk1(WAVES(ng)%Insteps))
        allocate( junk2(WAVES(ng)%Insteps))
        allocate( junk3(WAVES(ng)%Insteps))
        allocate( junk4(WAVES(ng)%Insteps))
        DO t=1,WAVES(ng)%Insteps
          junk1(t)=0.0_r8
          junk2(t)=0.0_r8
          junk3(t)=0.0_r8
          junk4(t)=REAL(t,r8)
        END DO
        fac4=WAVES(ng)%df
        fac5=WAVEG(ng)%pd
!
        DO f1=1,WAVES(ng)%nfreq-1
!
!  Do the double summation loop.
!
          DO f2=f1,WAVES(ng)%nfreq

# ifdef DISTRIBUTE
        IF (Master) THEN
# endif
!          write(*,*) 'f1 f2 are  ', f1,f2
# ifdef DISTRIBUTE
        END IF
# endif

            DDf=WAVES(ng)%f(f2)-WAVES(ng)%f(f1)
!
!  Limit the freq diff range from 1/400 to 1/30, that is
!  df=1/400=.0025 to df=1/30=0.033
!
            IF ((DDf.ge.fmin).and.(DDf.le.fmax)) THEN
  
              DO n1=1,WAVES(ng)%ndir
                DO n2=1,WAVES(ng)%ndir
!
! Compute the freq diff, dir diff, wave number, and phase of bound wave.
!
                  DDtheta=abs(WAVEG(ng)%wd(n1)-WAVEG(ng)%wd(n2))
                  k3=sqrt(k(f1)**2+k(f2)**2-                            &
     &               2.0_r8*k(f1)*k(f2)*COS(DDtheta))
                  Z_bw=Waves(ng)%phase(f2,n2)-Waves(ng)%phase(f1,n1)+pi

                  fac1=WAVES(ng)%f(f1)*WAVES(ng)%f(f2)
                  fac2=1.0_r8/cosh(k(f1)*h0)
                  fac3=1.0_r8/cosh(k(f2)*h0)
            
                  D1a=g*k(f1)*k(f2)/(8.0_r8*pi**2*fac1)*                &
     &                fac2*fac3
                  D3a=((2.0_r8*pi)**4*(fac1)**2)/(g**2)
                  D4=-0.5*(-WAVES(ng)%f(f1)*k(f2)**2*fac3**2+           &
     &                      WAVES(ng)%f(f2)*k(f1)**2*fac2**2)
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
                  k3_x=k3*cos(1.5_r8*pi-DIR3)   !DIR3 is in nautical coordinates
                  k3_y=k3*sin(1.5_r8*pi-DIR3)
!
!  Compute this contribution to the bndwave (water level) and associated 
!  currents; currents depend on direction
!  Here we need to include the phase difference along the boundary; this
!  needs to be done; bndwave needs to depend on the boudary (and so does
!  the free surface elevation)
!
# ifdef VARY_ACBC
                  DO j=SstrR,SendR
                  cff2=k3_x*dist_x(j)+k3_y*dist_y(j)
# else
                  cff2=0.0_r8
# endif
!
                  junk2(:)=-2.0_r8*pi*DDf*junk4(:)*dt(ng)+Z_bw+cff2
                  junk1(:)=A3*COS(junk2(:))

!  The currents need to be rotated to the roms grid plus they are
!  defined in the rho points
!  Here we assume these waves are shallow progressive water waves
# ifdef VARY_ACBC
                  WAVES(ng)%bndwave_z(j,:)=WAVES(ng)%bndwave_z(j,:)+    &
     &                                     junk1(:)
                  WAVES(ng)%bndwave_u(j,:)=WAVES(ng)%bndwave_u(j,:)+    &
     &                                     junk1(:)*oh0*sqrt(g*h0)
                  WAVES(ng)%bndwave_v(j,:)=WAVES(ng)%bndwave_v(j,:)+    &
     &                                     junk1(:)*oh0*sqrt(g*h0)
# else
                  WAVES(ng)%bndwave_z(:)=WAVES(ng)%bndwave_z(:)+junk1(:)
                  WAVES(ng)%bndwave_u(:)=WAVES(ng)%bndwave_u(:)+        &
     &                                   junk1(:)*oh0*sqrt(g*h0)
                  WAVES(ng)%bndwave_v(:)=WAVES(ng)%bndwave_v(:)+        &
     &                                   junk1(:)*oh0*sqrt(g*h0)
# endif
# ifdef VARY_ACBC
                END DO
# endif
                END DO
              END DO
            END IF
          END DO
        END DO
!
        deallocate(dist_x,dist_y)
        deallocate (k)
        deallocate( dist_local, junk1, junk2, junk3, junk4 )
      END IF

      RETURN

      END SUBROUTINE boundwave_tile

#endif

      END MODULE mod_inwave_swan
