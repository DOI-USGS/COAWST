#include "cppdefs.h"
      MODULE mod_inwave_swan
!
!svn $Id: swan_reader.F 1336 2008-01-24 02:45:56Z jcwarner $
! LAST CHANGE: mai 12/28/2010

!=======================================================================
!                                                                      !
!  NDIR        number of directional dimensions in Fourier components  !
!  NFREQ       number of frequency dimensions in Fourier components    !
!  NDIRSW      number of directional dimensions in SWAN spectra        !
!  NFREQSW     number of frequency dimensions in SWAN spectra          !
!  N           number of time steps                                    !
!  DUR         duration of the time series to generate                 !
!  N           number of time steps                                    !
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
!  Currently,  this routine is only executed in serial mode by the     !
!  main thread.                                                        !
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

# ifdef INWAVE_SWAN_COUPLING

        USE mod_kinds
        USE mod_inwave_params
        USE math_tools
        USE mod_iounits

        implicit none

        TYPE SHORT_WAVE

         integer             :: ndir,nfreq
         integer             :: ndirsw,nfreqsw
         integer             :: N
         real(r8)            :: dur,dt
         real(r8)            :: fnyq,fmax,fmin,df,Trep
         real(r8)            :: ddir,fp
         real(r8)            :: mainang, hm0gew, dang,int
         real(r8)            :: df_fft

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
      SUBROUTINE inwave_swan_run (ng)
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
!      USE mod_inwave_spectra

!  Imported variable declarations.
!
      integer, intent(in) :: ng

!  Local variable declarations.

      integer             :: dir_index,i

!
!-----------------------------------------------------------------------
! Read swan output spectra
!-----------------------------------------------------------------------
!
      call swan_reader (ng)
!
!-----------------------------------------------------------------------
! Generate the arrays for the frequency components and Fourier components
!-----------------------------------------------------------------------
!
      call array_gen (ng)
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
      DO dir_index=1,nd(ng)
        call FourierComp (dir_index,ng)
      ENDDO
!
!-----------------------------------------------------------------------
! Compute the Hilbert transform
!-----------------------------------------------------------------------
!

      call hilbert(WAVES(ng)%Comptemp_tot,size(WAVES(ng)%Comptemp_tot))

      WAVES(ng)%Ampzeta_tot(:)=abs(WAVES(ng)%Comptemp_tot)

!
!-----------------------------------------------------------------------
! Compute the wave energy
!-----------------------------------------------------------------------
!
      DO dir_index=1,nd(ng)
        DO i=1,WAVES(ng)%N
          WAVES(ng)%Ampzeta(i,dir_index)=0.5_r8*g*rho0*                 &
     &                                   (WAVES(ng)%Ampzeta_tot(i)*     &
     &                                   WAVES(ng)%SDD(dir_index)/      &
     &                                   WAVES(ng)%int)**2.0_r8
        ENDDO
      ENDDO

      RETURN

      END SUBROUTINE inwave_swan_run

!
!***********************************************************************
      SUBROUTINE array_gen (ng)
!***********************************************************************
!
!=======================================================================
!                                                                      !
!  Generates the arrays necesary to compute the free surface           !
!  elevation time series from the swan directional spectra             !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
      USE mod_inwave_vars
      USE mod_inwave_params
!      USE mod_inwave_spectra

      implicit none

!  Imported variable declarations.
!
      integer, intent(in) :: ng

!  Local variable declarations.

      integer :: i,j,p,n_pos_f
      real    :: N_ini

      real(r8), pointer     :: temp(:)   ! scratch array
      real(r8)              :: Fmin, Fmax
!
!-----------------------------------------------------------------------
! Define the Nyquist frequency, the maximum and minimum frequencies and 
! the number of directional bins
!-----------------------------------------------------------------------
!
      WAVES(ng)%fnyq=3.0_r8*WAVES(ng)%fp
      WAVES(ng)%fmax=WAVES(ng)%fnyq
      WAVES(ng)%fmin=WAVES(ng)%fp/20.0_r8
      WAVES(ng)%ndir=nd(ng)
!
!-----------------------------------------------------------------------
! Create the frequency and directional arrays for the fft
!-----------------------------------------------------------------------
!
      WAVES(ng)%nfreq=1000.0_r8    !!!! I choose this because I wanted

      allocate (temp(WAVES(ng)%nfreq))
      allocate (WAVES(ng)%f(WAVES(ng)%nfreq))
      allocate (WAVES(ng)%SD(WAVES(ng)%nfreq,nd (ng)))
      allocate (WAVES(ng)%SDD(nd(ng)))

      WAVES(ng)%df=(WAVES(ng)%fmax-WAVES(ng)%fmin)/(WAVES(ng)%nfreq-1)

      DO i=1,WAVES(ng)%nfreq
        WAVES(ng)%f(i)=REAL(i+1,r8)*WAVES(ng)%df+WAVES(ng)%fmin
      END DO
!
!-----------------------------------------------------------------------
! Interpolate from the SWAN 2D spectral grid to the 2D spectral grid that
! we predefine
!-----------------------------------------------------------------------
!
      call interp2d ( 1, WAVES(ng)%nfreqsw, 1, WAVES(ng)%ndirsw,        &
     &               WAVES(ng)%fsw, WAVES(ng)%theta, WAVES(ng)%SDSW,    &
     &               1, WAVES(ng)%nfreq, 1, nd (ng),                    &
     &               1, WAVES(ng)%nfreq, 1, nd (ng),                    &
     &               WAVES(ng)%f, WAVEG(ng)%wd,                         &
     &               WAVES(ng)%SD)

      WAVES(ng)%SDD = sum(WAVES(ng)%SD, DIM = 1)*WAVES(ng)%df

      WAVES(ng)%int=0.0_r8

      DO i=1,nd(ng)
        WAVES(ng)%int =WAVES(ng)%int+WAVES(ng)%SDD(i)*WAVEG(ng)%pd(i)
      ENDDO

      WAVES(ng)%dur=1./WAVES(ng)%df
      WAVES(ng)%N=nint(WAVES(ng)%dur/dt(ng))

      IF (mod(WAVES(ng)%N,2)/=0)THEN
        WAVES(ng)%N=WAVES(ng)%N+1
      ENDIF

!     allocate(WAVES(ng)%CompFn(WAVES(ng)%N))
      allocate(WAVES(ng)%CompFn(WAVES(ng)%nfreq))
      allocate(WAVES(ng)%zeta(WAVES(ng)%N,WAVES(ng)%ndir))
      allocate(WAVES(ng)%Ampzeta(WAVES(ng)%N,WAVES(ng)%ndir))
      allocate(WAVES(ng)%Ampzeta_tot(WAVES(ng)%N))

      DO j=1,WAVES(ng)%ndir
        DO i=1,WAVES(ng)%N    
          WAVES(ng)%zeta(i,j)=0.0_r8
          WAVES(ng)%Ampzeta(i,j)=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
! Create frequency and directional arrays for the spectra
!-----------------------------------------------------------------------
!
      allocate (WAVES(ng)%position(size(WAVES(ng)%f)))
      allocate (WAVES(ng)%phase(size(WAVES(ng)%f),WAVES(ng)%ndir))
      allocate (WAVES(ng)%amp(size(WAVES(ng)%f),WAVES(ng)%ndir))

      RETURN

      END SUBROUTINE array_gen

!
!***********************************************************************
      SUBROUTINE random_phase (ng)
!***********************************************************************
!
      USE mod_scalars
!      USE mod_inwave_spectra

      implicit none

!  Imported variable declarations.
!
      integer, intent(in) :: ng
      
!  Local variable declarations.

      integer               :: i,j
      real(r8)              :: twopi

      call random_number(WAVES(ng)%phase)

      twopi= 2.0_r8*pi

      DO j=1,WAVES(ng)%ndir
        DO i=1,size(WAVES(ng)%f)
          WAVES(ng)%phase(i,j)=WAVES(ng)%phase(i,j)*twopi
        END DO
      END DO

      RETURN

      END SUBROUTINE random_phase

!
!***********************************************************************
      SUBROUTINE tpdcalc(Sf,f,Trep)
!***********************************************************************
!
      USE mod_inwave_bound

      implicit none

      real(r8), intent(in)  :: Sf(:), f(:)
      real(r8), pointer     :: temp(:)
      real(r8)              :: Trep


      allocate(temp(size(Sf)))
      temp=0.0_r8
      where (Sf>0.8_r8*maxval(Sf))
       temp=1.0_r8
      end where

      Trep=sum(temp*Sf*f)/sum(temp*Sf)

      RETURN

      END SUBROUTINE tpdcalc

!
!***********************************************************************
      SUBROUTINE amplitudes (ng)
!***********************************************************************
!

      use mod_scalars
      use mod_inwave_params
!      USE mod_inwave_spectra

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
      cff2=WAVES(ng)%dang
      DO j=1,WAVES(ng)%ndir
        DO i=1,WAVES(ng)%nfreq
          WAVES(ng)%amp(i,j)= sqrt(2.0_r8*WAVES(ng)%SD(i,j)*cff1*cff2)
        END DO
      END DO

! Assing a position in the spectral frequency arrat to each Fourier component

      DO i=1,WAVES(ng)%nfreq
        WAVES(ng)%position(i)=i
      END DO

      RETURN

      END SUBROUTINE amplitudes

!
!***********************************************************************
      SUBROUTINE FourierComp (dir_index,ng)
!***********************************************************************
!

      USE mod_scalars
!      USE mod_inwave_spectra

      implicit none

!  Imported variable declarations.
!
      integer, intent(in) :: ng

!  Local variable declarations.
!

      integer                         :: i,F2
      integer                         :: dir_index,nn
      complex(kind(0.0d0))            :: compi

      WAVES(ng)%CompFn=0.0_r8

      F2=0.0_r8

      compi = (0.0_r8,1.0_r8)

! COMPUTES THE FOURIER COMPONENTS

      DO i=1,WAVES(ng)%nfreq
!     DO i=1,WAVES(ng)%N
        WAVES(ng)%CompFn(i)=0.5_r8*WAVES(ng)%amp(i,dir_index)*          &
     &                      exp(compi*WAVES(ng)%phase(i,dir_index))
      END DO

      nn=size(WAVES(ng)%CompFn(WAVES(ng)%N/2+2:WAVES(ng)%N))

      allocate(WAVES(ng)%Comptemp(nn))

      WAVES(ng)%Comptemp=conjg(WAVES(ng)%CompFn(2:WAVES(ng)%N/2))

      call flipiv(WAVES(ng)%Comptemp,size(WAVES(ng)%Comptemp))

      WAVES(ng)%CompFn(WAVES(ng)%N/2+2:WAVES(ng)%N)=WAVES(ng)%Comptemp

      deallocate(WAVES(ng)%Comptemp)

      allocate(WAVES(ng)%Comptemp(size(WAVES(ng)%CompFn)))

      IF (dir_index.eq.1)THEN
        allocate(WAVES(ng)%Comptemp_tot(size(WAVES(ng)%CompFn)))
      ENDIF

      WAVES(ng)%Comptemp=fft(WAVES(ng)%CompFn,inv=.true.,stat=F2)

! SCALE FACTOR

      WAVES(ng)%Comptemp=WAVES(ng)%Comptemp/                            &
     &                   sqrt(real(size(WAVES(ng)%Comptemp)))

      WAVES(ng)%zeta(:,dir_index)=dble(WAVES(ng)%Comptemp*WAVES(ng)%N)

      WAVES(ng)%Comptemp=WAVES(ng)%zeta(:,dir_index)

      WAVES(ng)%Comptemp_tot=WAVES(ng)%Comptemp_tot+WAVES(ng)%Comptemp

! COMPUTES THE HILBERT TRANSFORM

      CALL hilbert(WAVES(ng)%Comptemp,size(WAVES(ng)%Comptemp))

      WAVES(ng)%Ampzeta(:,dir_index)=abs(WAVES(ng)%Comptemp)

      deallocate(WAVES(ng)%Comptemp)

      RETURN

      END SUBROUTINE FourierComp

# endif

      END MODULE mod_inwave_swan