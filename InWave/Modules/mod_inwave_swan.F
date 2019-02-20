#include "cppdefs.h"
#define AMPZETAFILT

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

      integer :: i, j, dir_index, numavg, numavgh, start
      integer, allocatable :: peaks(:)
      real(r8) :: cff
      real(r8), allocatable :: zeta_filt(:), zeta_filt2(:)
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
      call boundwave (ng, tile, first)
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
#ifdef AMPZETAFILT
      allocate (zeta_filt2(WAVES(ng)%Insteps))
!  find peaks of
      allocate (peaks(WAVES(ng)%Insteps))
      DO i=1,WAVES(ng)%Insteps
        peaks(i)=0
        zeta_filt(i)=WAVES(ng)%Ampzeta_tot(i)
      END DO
! 
      peaks(1)=1
      peaks(2)=0
      peaks(WAVES(ng)%Insteps-1)=0
      peaks(WAVES(ng)%Insteps  )=1
      DO i=3,WAVES(ng)%Insteps-2
        IF ((WAVES(ng)%Ampzeta_tot(i).gt.WAVES(ng)%Ampzeta_tot(i-1))    &
     &    .and.(WAVES(ng)%Ampzeta_tot(i).gt.WAVES(ng)%Ampzeta_tot(i+1)))&
     &    THEN
          peaks(i)=1
        ELSE
          peaks(i)=0
        END IF
      END DO
!
      start=2
      DO i=3,WAVES(ng)%Insteps-2
        IF (peaks(i).eq.1) THEN
          DO j=start,i-1
            zeta_filt(j)=(WAVES(ng)%Ampzeta_tot(i)-                     &
     &                    WAVES(ng)%Ampzeta_tot(start))/                &
     &                    REAL((i-start),r8)*REAL(j-start,r8)+          &
     &                    WAVES(ng)%Ampzeta_tot(start)
          END DO
          start=i
        ELSE
        END IF
      END DO
      DO i=1,WAVES(ng)%Insteps
        zeta_filt2(i)=zeta_filt(i)
      END DO
!
! do a filter on zeta_filt
!
      cff=0.0_r8
      numavg=INT(30.0_r8/dt(ng))
      numavg=numavg-1*(1-MOD(numavg,2)) !force odd
      DO i=1,numavg
        cff=cff+zeta_filt(i)
      ENDDO
      numavgh=(numavg-1)/2
      DO i=numavgh+1,WAVES(ng)%Insteps-numavgh
        zeta_filt2(i)=cff/REAL(numavg,r8)
        cff=cff-zeta_filt(i-numavgh)+zeta_filt(i+numavgh)
      ENDDO
!
      DO dir_index=1,ND
        DO i=1,WAVES(ng)%Insteps
          WAVES(ng)%Ampzeta(i,dir_index)=0.5_r8*g*rho0*                 &
     &                                   (zeta_filt2(i)*                &
     &                                   WAVES(ng)%SDD(dir_index)/      &
     &                                   WAVES(ng)%int)**2.0_r8
        ENDDO
      ENDDO
      deallocate (peaks, zeta_filt2)
#else
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
#endif
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
        allocate(WAVES(ng)%bndwave(WAVES(ng)%Insteps))
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

      twopi=2.0_r8*pi

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
     &                     WAVES(ng)%Insteps, WAVES(ng)%bndwave)
      RETURN
      END SUBROUTINE boundwave
!
!***********************************************************************
      SUBROUTINE boundwave_tile (ng, tile, first,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Insteps, bndwave)
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
      USE distribute_mod, ONLY : mp_bcasti, mp_gather2d
# endif
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, first, Insteps
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: bndwave(:)
#else
      real(r8), intent(inout) :: bndwave(1:Insteps)
#endif
!
!  Local variable declarations.
!
      integer                         :: i, j, p1, p2, f1, f2, dum
      integer                         :: Npts, MyError
      real(r8)                        :: cff, cff1, cff2, cff3, cff4
      real(r8)                        :: twopi, otwopi
      real(r8)                        :: error, Tr_min, fmin, fmax
      real(r8)                        :: fac1, fac2, fac3, fac4
      real(r8)                        :: L0, k0, k1, kh, wr
      real(r8)                        :: F, FDER, tanhkh
      real(r8)                        :: DDf, A3, Z_bw
      real(r8)                        :: DDtheta, k3
      real(r8)                        :: D1, D2, D3, D4, D1a, D3a, DTOT
      real(r8)                        :: h0, E3
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: Awrk
      real(r8), dimension((Lm(ng)+2),(Mm(ng)+2)) :: h_local
      real(r8), dimension((Lm(ng)+2),(Mm(ng)+2)) :: zeta_local

      real(r8), allocatable           :: E1d(:), D1d(:), P1d(:)
      real(r8), allocatable           :: k(:)
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
      DO j=1,Mm(ng)+2
        DO i=1,Lm(ng)+2
          h_local(i,j)=0.0_r8
          zeta_local(i,j)=0.0_r8
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_gather2d (ng, iNLM, LBi, UBi, LBj, UBj,                   &
     &                  0, r2dvar, 1.0_r8,                              &
     &                  GRID(ng) % rmask,                               &
     &                  GRID(ng)%h(:,:), Npts, Awrk)
      p1=0
      DO j=1,Mm(ng)+2
        DO i=1,Lm(ng)+2
          p1=p1+1
          h_local(i,j)=Awrk(p1)
        END DO
      END DO
!
      CALL mp_gather2d (ng, iNLM, LBi, UBi, LBj, UBj,                   &
     &                  0, r2dvar, 1.0_r8,                              &
     &                  GRID(ng) % rmask,                               &
     &                  OCEAN(ng)%zeta(:,:,1), Npts, Awrk)
      p1=0
      DO j=1,Mm(ng)+2
        DO i=1,Lm(ng)+2
          p1=p1+1
          zeta_local(i,j)=Awrk(p1)
        END DO
      END DO
# else
      DO j=1,Mm(ng)+2
        DO i=1,Lm(ng)+2
          h_local(i,j)=GRID(ng)%h(i-1,j-1)
          zeta_local(i,j)=OCEAN(ng)%zeta(i-1,j-1,1)
        END DO
      END DO
# endif
!
!  Now set h0 based on the middle boundary total depth.
!
# ifdef DISTRIBUTE
      IF (Master) THEN
# endif
      IF (LBC(iwest,isAC3d,ng)%acquire) THEN
        j=INT(Mm(ng)/2)
        h0=h_local(1,j)+zeta_local(1,j)
      ELSE IF (LBC(ieast,isAC3d,ng)%acquire) THEN
        j=INT(Mm(ng)/2)
        h0=h_local(Lm(ng)+2,j)+zeta_local(Lm(ng)+2,j)
      ELSE IF (LBC(inorth,isAC3d,ng)%acquire) THEN
        i=INT(Lm(ng)/2)
        h0=h_local(i,Mm(ng)+2)+zeta_local(i,Mm(ng)+2)
      ELSE IF (LBC(isouth,isAC3d,ng)%acquire) THEN
        i=INT(Lm(ng)/2)
        h0=h_local(i,1)+zeta_local(i,1)
      END IF        
# ifdef DISTRIBUTE
      END IF
      CALL MPI_BCAST(h0, 1, MP_FLOAT, 0,                                &
     &               OCN_COMM_WORLD, MyError)
# endif
# ifdef DISTRIBUTE
      IF (Master) THEN
# endif
      write(*,*) 'Mean offshore water depth for bndwave is ', h0
# ifdef DISTRIBUTE
      END IF
# endif

!
!  Allocate bouundwave and local wavenumber arrays.
!  Determine length along the active side. This is used to 
!  phase the bndwave along the boundary.
!
      allocate (k(WAVES(ng)%nfreq))
      DO i=1,WAVES(ng)%Insteps
          bndwave(i)=0.0_r8
      END DO
!
! Compute the wave number or each freq bin
! dont we already have this...
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
!  Make it a 1D spectrum for now to compute boudnwave.
!
      allocate (E1d(WAVES(ng)%nfreq))
      allocate (D1d(WAVES(ng)%nfreq))
      allocate (P1d(WAVES(ng)%nfreq))
      DO i=1,WAVES(ng)%nfreq
        E1d(i)=0.0_r8
        D1d(i)=0.0_r8
        P1d(i)=0.0_r8
      ENDDO
!
! Compute mean direction.
!
      DO i=1,WAVES(ng)%nfreq
        cff1=0.0_r8
        cff2=0.0_r8
        DO j=1,ND
          cff1=cff1+WAVES(ng)%SD(i,j)*WAVEG(ng)%wd(j)
          cff2=cff2+WAVES(ng)%SD(i,j)
        END DO
        D1d(i)=cff1/(cff2+eps)
      END DO
!
! Compute 1D Energy.
!
      DO i=1,WAVES(ng)%nfreq
        DO j=1,ND
          E1d(i)=E1d(i)+WAVES(ng)%SD(i,j)*WAVEG(ng)%pd
        END DO
      END DO
!
! Compute Random number for phase.
!
        call random_number(P1d)
# ifdef DISTRIBUTE
      IF (Master) THEN
# endif
        DO i=1,WAVES(ng)%nfreq
          P1d(i)=P1d(i)*twopi
        END DO
# ifdef DISTRIBUTE
      END IF
# endif
# ifdef DISTRIBUTE
!
!  Scatter phase to all the nodes.
!
      Npts=WAVES(ng)%nfreq
      allocate(wrk(Npts))
      IF (Master) THEN
        j=0
        DO i=1,WAVES(ng)%nfreq
            j=j+1
            wrk(j)=P1d(i)
        END DO
      END IF
      CALL MPI_BCAST(wrk, Npts, MP_FLOAT, 0,                            &
     &               OCN_COMM_WORLD, MyError)
      j=0
      DO i=1,WAVES(ng)%nfreq
          j=j+1
          P1d(i)=wrk(j)
      END DO
!
      deallocate(wrk)
# endif
!
! Compute the energy transfer for each pair of frequency bins.
!
# ifdef DISTRIBUTE
      IF (Master) THEN
# endif
      write(*,*) 'Computing bound wave'
# ifdef DISTRIBUTE
      END IF
# endif
      fac4=WAVES(ng)%df
      DO f1=1,WAVES(ng)%nfreq-1
!
!  Do the double summation loop.
!
        DO f2=f1,WAVES(ng)%nfreq
          fac1=WAVES(ng)%f(f1)*WAVES(ng)%f(f2)
          fac2=1.0_r8/cosh(k(f1)*h0)
          fac3=1.0_r8/cosh(k(f2)*h0)
          DDf=WAVES(ng)%f(f2)-WAVES(ng)%f(f1)
!
!  Limit the freq diff range from 1/400 to 1/30, that is
!  df=1/400=.0025 to df=1/30=0.033
!
          IF ((DDf.ge.fmin).and.(DDf.le.fmax)) THEN
!
!  Do part of the D summations here that are not funcions of dtheta.
!
            D1a=g*k(f1)*k(f2)/(8.0_r8*pi**2*fac1)*                      &
     &          fac2*fac3
            D3a=((2.0_r8*pi)**4*(fac1)**2)/(g**2)
            D4=-0.5*(-WAVES(ng)%f(f1)*k(f2)**2*fac3**2+                 &
     &                WAVES(ng)%f(f2)*k(f1)**2*fac2**2)
!
! Compute the freq diff, dir diff, wave number, and phase of bound wave.
!
            DDtheta=D1d(f2)-D1d(f1)
            k3=sqrt(k(f1)**2+k(f2)**2-2.0_r8*k(f1)*k(f2)*COS(DDtheta))
            Z_bw=P1d(f1)-P1d(f2)+pi
!
!  Compute the double summation.
!
            D1=D1a*COS(DDtheta+pi)*MIN(cosh(k3*h0),1.0E10_r8)
            D2=-g*(DDf)/((g*k3*tanh(k3*h0)-(2.0_r8*pi)**2*DDf**2)*      &
     &             fac1+eps)
            D3=DDf*(D3a-k(f1)*k(f2)*COS(DDtheta+pi))
            DTOT=D1+D2*(D3+D4)
!
!  Compute energy of bound wave
!
            E3=2.0_r8*DTOT**2*                                          &
     &         E1d(f1)*E1d(f2)*fac4
            A3=SQRT(2.0_r8*E3*fac4)
!
!  Compute this contribution to the bndwave.
!  NOT sure about the 0.5 here
!
            DO i=1,WAVES(ng)%Insteps
              cff=REAL((i-1),r8)*dt(ng)
              cff1=2.0_r8*pi*DDf*cff+Z_bw
              bndwave(i)=bndwave(i)+0.5_r8*A3*COS(cff1)
            END DO
          END IF
        END DO
      END DO
!
      deallocate (k, E1D, D1d, P1d)

      RETURN

      END SUBROUTINE boundwave_tile

#endif

      END MODULE mod_inwave_swan
