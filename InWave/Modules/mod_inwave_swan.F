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
!  Comptemp_tot    Fourier components                                  !
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

         integer             :: ndir,nfreq
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
         real(r8), pointer :: Ampzeta(:,:)
         real(r8), pointer :: Ampzeta_tot(:)
         real(r8), pointer :: position(:)
         real(r8), pointer :: bndwave(:)

         complex(fftkind),dimension(:),allocatable   :: CompFn, Comptemp
         complex(fftkind),dimension(:),allocatable   :: Comptemp_tot

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

      integer :: i, dir_index, numavg, numavgh
      real(r8) :: cff
      real(r8), allocatable :: zeta_filt(:)
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
! Generate random phases fro each frequency component
!-----------------------------------------------------------------------
!
      call random_phase (ng)
!
!-----------------------------------------------------------------------
! Compute the amplitudes fro each frequency component
!-----------------------------------------------------------------------
!
      call amplitudes (ng)
!
!-----------------------------------------------------------------------
! Compute the Fourier Components
!-----------------------------------------------------------------------
!
      DO dir_index=1,ND
        call FourierComp (dir_index, ng, first)
      ENDDO
!
!-----------------------------------------------------------------------
! Compute the boundwave
!-----------------------------------------------------------------------
!
! jcw     call boundwave (ng, first)
!
!-----------------------------------------------------------------------
! Compute the Hilbert transform
!-----------------------------------------------------------------------
!
      call hilbert(WAVES(ng)%Comptemp_tot,size(WAVES(ng)%Comptemp_tot))

      WAVES(ng)%Ampzeta_tot(:)=abs(WAVES(ng)%Comptemp_tot)

!
!-----------------------------------------------------------------------
! Compute the wave energy, scale for each direction
!-----------------------------------------------------------------------
!
      allocate (zeta_filt(WAVES(ng)%Insteps))
      DO dir_index=1,ND
        DO i=1,WAVES(ng)%Insteps
          WAVES(ng)%Ampzeta(i,dir_index)=0.5_r8*g*rho0*                 &
     &                                   (WAVES(ng)%Ampzeta_tot(i)*     &
     &                                   WAVES(ng)%SDD(dir_index)/      &
     &                                   WAVES(ng)%int)**2.0_r8
          zeta_filt(i)=WAVES(ng)%Ampzeta(i,dir_index)
        ENDDO
! do a filter on ampzeta
        cff=0.0_r8
        numavg=INT(5.0_r8/dt(ng))
        numavg=numavg-1*(1-MOD(numavg,2)) !force odd
        DO i=1,numavg
          cff=cff+zeta_filt(i)
        ENDDO
        numavgh=(numavg-1)/2
        DO i=numavgh+1,WAVES(ng)%Insteps-numavgh
          WAVES(ng)%Ampzeta(i,dir_index)=cff/REAL(numavg,r8)
          cff=cff-zeta_filt(i-numavgh)+zeta_filt(i+numavgh)
        ENDDO
      ENDDO
      deallocate (zeta_filt)

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
      integer             :: i, j, k, p, n_pos_f, offset
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
        WAVES(ng)%nfreq=1000 ! because Mai wanted it that way.
        WAVES(ng)%fnyq=3.0_r8*WAVES(ng)%fp
!       WAVES(ng)%fmax=WAVES(ng)%fnyq
!       WAVES(ng)%df=(WAVES(ng)%fmax-WAVES(ng)%fmin)/                   &
!   &                REAL((WAVES(ng)%nfreq-1),r8)
        WAVES(ng)%fmax=0.5_r8
        WAVES(ng)%df=WAVES(ng)%fmax/WAVES(ng)%nfreq
        WAVES(ng)%ndir=ND
        WAVES(ng)%fmin=WAVES(ng)%df
        IF (Master) THEN
          WRITE(stdout,*) 'Computing AC boundary forcing'
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
        WAVES(ng)%Insteps=nint(WAVES(ng)%dur/dt(ng))
!
!  Allocate and init the computational arrays.
!
        allocate(WAVES(ng)%CompFn(WAVES(ng)%Insteps))
        allocate(WAVES(ng)%zeta(WAVES(ng)%Insteps,WAVES(ng)%ndir))
        allocate(WAVES(ng)%Ampzeta(WAVES(ng)%Insteps,WAVES(ng)%ndir))
        allocate(WAVES(ng)%Ampzeta_tot(WAVES(ng)%Insteps))
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
      DO j=1,WAVES(ng)%ndir
        DO i=1,WAVES(ng)%Insteps
          WAVES(ng)%zeta(i,j)=0.0_r8
          WAVES(ng)%Ampzeta(i,j)=0.0_r8
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

!
!***********************************************************************
      SUBROUTINE FourierComp (dir_index, ng, first)
!***********************************************************************
!
      USE mod_scalars

      implicit none

!  Imported variable declarations.
!
      integer, intent(in) :: dir_index, ng, first

!  Local variable declarations.
!
      integer                         :: i, j
      real(r8)                        :: cff
      real(r8)                        :: twopi
      complex(kind(0.0d0))            :: compi

      twopi=2.0_r8*pi

      compi = (0.0_r8,1.0_r8)

! COMPUTES THE FOURIER COMPONENTS
      DO i=1,WAVES(ng)%Insteps
        WAVES(ng)%CompFn(i)=0.0_r8
        DO j=1,WAVES(ng)%nfreq
          cff=REAL((i-1),r8)*dt(ng)
          WAVES(ng)%CompFn(i)=WAVES(ng)%CompFn(i)+                      &
     &                        WAVES(ng)%amp(j,dir_index)*               &
     &                        COS(twopi*WAVES(ng)%f(j)*cff+             &
     &                            WAVES(ng)%phase(j,dir_index))
        END DO
      END DO

      IF ((dir_index.eq.1).and.(first.eq.1)) THEN
        allocate(WAVES(ng)%Comptemp_tot(WAVES(ng)%Insteps))
      ENDIF

      IF (dir_index.eq.1) THEN
        DO i=1,WAVES(ng)%Insteps
          WAVES(ng)%Comptemp_tot(i)=0.0_r8
        END DO
      ENDIF

      DO i=1,WAVES(ng)%Insteps
        WAVES(ng)%Comptemp_tot(i)=WAVES(ng)%Comptemp_tot(i)+            &
     &                            WAVES(ng)%CompFn(i)
      END DO

      RETURN

      END SUBROUTINE FourierComp
!
!***********************************************************************
      SUBROUTINE boundwave (ng, first)
!                                                                       !
!                compute the                                            !
!                wave envelope and the associated bound wave (Hasselman,!
!                1962; Herbers et al., 1994; Van Dongeren et al., 2003).!
!                A double summation technique                           !
!***********************************************************************
!
      USE mod_scalars
      USE mod_inwave_vars

      implicit none

!  Imported variable declarations.
!
      integer, intent(in) :: ng, first

!  Local variable declarations.
!
      integer                         :: i, j, p, f1, f2, dum
      real(r8)                        :: cff, cff1, cff2, twopi, otwopi
      real(r8)                        :: error, Tr_min
      real(r8)                        :: fac1, fac2, fac3
      real(r8)                        :: L0, k0, k1, kh, wr
      real(r8)                        :: F, FDER, tanhkh
      real(r8)                        :: DDf, DDtheta, k3, Z_bw
      real(r8)                        :: D1, D2, D3, D4, DTOT, E3, A3
      real(r8)                        :: h0

      real(r8), allocatable           :: k(:)

      real(r8), parameter :: maxErr = 0.1_r8
!     real(r8), parameter :: kwc_max = 10.0_r8

      twopi=2.0_r8*pi
      otwopi=1.0_r8/(2.0_r8*pi)
      Tr_min=1.0_r8
      h0=15.0_r8
!
!  Allocate bouundwave array.
!
      allocate (k(WAVES(ng)%nfreq))
      IF (first.eq.1) THEN
        allocate(WAVES(ng)%bndwave(WAVES(ng)%Insteps))
      ENDIF
      DO i=1,WAVES(ng)%Insteps
        WAVES(ng)%bndwave(i)=0.0_r8
      END DO
!
! Compute the wave number or each freq bin
!
      DO i=1,WAVES(ng)%nfreq
        L0=g*otwopi*WAVEG(ng)%Trep**2.0_r8
        k0=twopi/L0
        error=100.0_r8
        wr=twopi/MAX(Tr_min,WAVEG(ng)%Trep)
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
! Compute the energy transfer for each pair of frequency bins
!
!     dum=find(f<fcutoff);
!     DO f1=dum(end):nf-1

      DO f1=1,WAVES(ng)%nfreq-1
        DO f2=f1,WAVES(ng)%nfreq
          fac1=WAVES(ng)%f(f1)*WAVES(ng)%f(f2)
          fac2=1.0_r8/cosh(k(f1)*h0)
          fac3=1.0_r8/cosh(k(f2)*h0)
          DO p=1,ND
          
          write(*,*) 'doing f1 f2 p ', f1,f2,p,WAVES(ng)%nfreq,ND
!
! Compute the freq diff, dir diff, wave number, and phase of bound wave.
!
            DDf=WAVES(ng)%f(f2)-WAVES(ng)%f(f1)
!           DDtheta=(WAVES(ng)%Sdir(f2)-WAVES(ng)%Sdir(f1))*deg2rad
            DDtheta=WAVEG(ng)%pd*(p-1)
            k3=sqrt(k(f2)**2+k(f1)**2-2.0_r8*k(f2)*k(f1)*COS(DDtheta))
            Z_bw=WAVES(ng)%phase(f1,p)-WAVES(ng)%phase(f2,p)+pi
!
!  Compute the double summation.
!
            D1=g*k(f1)*k(f2)*COS(DDtheta+pi)/                           &
     &         (8.0_r8*pi**2*fac1)*                                     &
     &         cosh(k3*h0)*fac2*fac3
            D2=-g*(DDf)/((g*k3*tanh(k3*h0)-(2.0_r8*pi)**2*DDf**2)*      &
     &                   fac1)
            D3=DDf*(((2.0_r8*pi)**4*(fac1)**2)/                         &
     &         (g**2)-k(f1)*k(f2)*COS(DDtheta+pi))
            D4=-0.5*(-WAVES(ng)%f(f1)*k(f2)**2*fac3**2+                 &
     &                WAVES(ng)%f(f2)*k(f1)**2*fac2**2)
            DTOT=D1+D2*(D3+D4)

!            if DDf>fcutoff
!              DTOT(f1,f2)=0;  % do not consider frequencies higher than fcutoff
!            end   
!        
!            if DDf<fmin
!              DTOT(f1,f2)=0;  % do not consider frequencies lower than fmin
!            end   
!            if (isnan(DTOT(f1,f2))==1) 
!              DTOT(f1,f2)=0;
!            end

!
!  Compute energy of bound wave
!
!           cff=2.0_r8*DTOT**2* Si(f1) * Si(f2) * df
            cff=2.0_r8*DTOT**2*WAVES(ng)%SD(f1,p)*WAVES(ng)%SD(f2,p)*   &
     &          WAVEG(ng)%pd*WAVEG(ng)%pd*WAVES(ng)%df
            A3=SQRT(2.0_r8*cff*WAVES(ng)%df);
!
!            IF (DDf.le.fcutoff) THEN
             DO i=1,WAVES(ng)%Insteps
               cff=REAL((i-1),r8)*dt(ng)
               WAVES(ng)%bndwave(i)=WAVES(ng)%bndwave(i)+0.5_r8*A3*     &
     &                              COS(2.0_r8*pi*DDf*cff+Z_bw)
             END DO
!            ENDIF
          END DO
        END DO
      END DO
      deallocate (k)

      RETURN

      END SUBROUTINE boundwave

#endif

      END MODULE mod_inwave_swan
