!
!svn $Id: sedflocs_mod.h 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group        John C. Warner   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Sediment Floc Model Kernel Variables:                               !
!                                                                      !
#if defined SED_FLOCS
!  bottom         Exposed sediment layer properties:                   !
!                   bottom(:,:,isd50) => mean grain diameter           !
!                   bottom(:,:,idens) => mean grain density            !
!                   bottom(:,:,iwsed) => mean settling velocity        !
!                   bottom(:,:,itauc) => mean critical erosion stress  !
#endif
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
        logical :: l_ASH
        logical :: l_ADS
        logical :: l_COLLFRAG
        logical :: l_testcase
        INTEGER  :: f_ero_iv
        real(r8) :: f_dp0,f_alpha,f_beta,f_nb_frag
        real(r8) :: f_dmax,f_ater,f_clim
        real(r8) :: f_ero_frac,f_ero_nbfrag
        real(r8) :: f_nf
        real(r8) :: f_frag
        real(r8) :: f_fter
        real(r8) :: f_collfragparam
        real(r8), parameter :: rhoref = 1030.0_r8
!
      TYPE T_SEDFLOCS
!
#if defined SEDIMENT
        real(r8), pointer :: f_diam(:)
        real(r8), pointer :: f_vol(:)
        real(r8), pointer :: f_rho(:)
        real(r8), pointer :: f_cv(:)
        real(r8), pointer :: f_l3(:)
        real(r8), pointer :: f_mass(:)
        real(r8), pointer :: f_coll_prob_sh(:,:)
        real(r8), pointer :: f_coll_prob_ds(:,:)
        real(r8), pointer :: f_l1_sh(:,:)
        real(r8), pointer :: f_l1_ds(:,:)
        real(r8), pointer :: f_g3(:,:)
        real(r8), pointer :: f_l4(:,:)
        real(r8), pointer :: f_g1_sh(:,:,:)
        real(r8), pointer :: f_g1_ds(:,:,:)
        real(r8), pointer :: f_g4(:,:,:)
#endif

      END TYPE T_SEDFLOCS

      TYPE (T_SEDFLOCS), allocatable :: SEDFLOCS(:)

      CONTAINS

      SUBROUTINE allocate_sedflocs (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_sediment
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate structure variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( SEDFLOCS(Ngrids) )
!
!  Nonlinear model state.
!
#if defined SEDIMENT
      allocate ( SEDFLOCS(ng) % f_diam(NCS) )
      allocate ( SEDFLOCS(ng) % f_vol(NCS) )
      allocate ( SEDFLOCS(ng) % f_rho(NCS) )
      allocate ( SEDFLOCS(ng) % f_cv(NCS) )
      allocate ( SEDFLOCS(ng) % f_l3(NCS) )
      allocate ( SEDFLOCS(ng) % f_mass(0:NCS+1) )
      allocate ( SEDFLOCS(ng) % f_coll_prob_sh(NCS,NCS) )
      allocate ( SEDFLOCS(ng) % f_coll_prob_ds(NCS,NCS) )
      allocate ( SEDFLOCS(ng) % f_l1_sh(NCS,NCS) )
      allocate ( SEDFLOCS(ng) % f_l1_ds(NCS,NCS) )
      allocate ( SEDFLOCS(ng) % f_g3(NCS,NCS) )
      allocate ( SEDFLOCS(ng) % f_l4(NCS,NCS) )
      allocate ( SEDFLOCS(ng) % f_g1_sh(NCS,NCS,NCS) )
      allocate ( SEDFLOCS(ng) % f_g1_ds(NCS,NCS,NCS) )
      allocate ( SEDFLOCS(ng) % f_g4(NCS,NCS,NCS) )
#endif


      RETURN
      END SUBROUTINE allocate_sedflocs

      SUBROUTINE initialize_sedflocs (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize structure variables in the module using     !
!  first touch distribution policy. In shared-memory configuration,    !
!  this operation actually performs the propagation of the "shared     !
!  arrays" across the cluster,  unless another policy is specified     !
!  to  override the default.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_sediment
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, itrc, j, k

      real(r8), parameter :: IniVal = 0.0_r8

#include "set_bounds.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (WESTERN_EDGE) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#endif
!
!-----------------------------------------------------------------------
!  Initialize sediment structure variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        CALL initialize_sedflocs_param (ng, tile,                       &
     &                     SEDFLOCS(ng) % f_mass,                       &
     &                     SEDFLOCS(ng) % f_diam,                       &
     &                     SEDFLOCS(ng) % f_g1_sh,                      &
     &                     SEDFLOCS(ng) % f_g1_ds,                      &
     &                     SEDFLOCS(ng) % f_g3,                         &
     &                     SEDFLOCS(ng) % f_l1_sh,                      &
     &                     SEDFLOCS(ng) % f_l1_ds,                      &
     &                     SEDFLOCS(ng) % f_coll_prob_sh,               &
     &                     SEDFLOCS(ng) % f_coll_prob_ds,               &
     &                     SEDFLOCS(ng) % f_l3)
!
      END IF
!
      RETURN
      END SUBROUTINE initialize_sedflocs
!
!***********************************************************************
      SUBROUTINE initialize_sedflocs_param (ng, tile,                   &
     &                              f_mass,f_diam,f_g1_sh,              &
     &                              f_g1_ds,f_g3,f_l1_sh,f_l1_ds,       &
     &                              f_coll_prob_sh,f_coll_prob_ds,      &
     &                              f_l3)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
!
      implicit none 
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      real(r8), intent(inout) :: f_mass(0:NCS+1)
      real(r8), intent(inout) :: f_diam(NCS)
      real(r8), intent(inout) :: f_g1_sh(NCS,NCS,NCS)
      real(r8), intent(inout) :: f_g1_ds(NCS,NCS,NCS)
      real(r8), intent(inout) :: f_g3(NCS,NCS)
      real(r8), intent(inout) :: f_l1_sh(NCS,NCS)
      real(r8), intent(inout) :: f_l1_ds(NCS,NCS)
      real(r8), intent(inout) :: f_coll_prob_sh(NCS,NCS)
      real(r8), intent(inout) :: f_coll_prob_ds(NCS,NCS)
      real(r8), intent(inout) :: f_l3(NCS)
!
!  Local variable declarations.
!
      logical  :: f_test
      real(r8) :: f_weight,mult,dfragmax
      integer  :: iv1,iv2,iv3,iv,itrc
      real(r8) :: f_vol(NCS),f_rho(NCS)
      real(r8), parameter :: mu = 0.001_r8
      real(r8) :: eps 
      eps = epsilon(1.0)
!
! ALA the rest of the initialization
!
!      l_ADS=.false.
!      l_ASH=.true.
!      l_COLLFRAG=.false.
!      f_dp0=0.000004_r8
!      f_nf=1.9_r8
!      f_dmax=0.001500_r8
!      f_nb_frag=2.0_r8
!      f_alpha=0.35_r8
!      f_beta=0.15_r8
!      f_ater=0.0_r8
!      f_ero_frac=0.0_r8
!      f_ero_nbfrag=2.0_r8
!      f_ero_iv=1
!      f_collfragparam=0.01_r8


!!--------------------------------------------------
!! floc characteristics
      DO itrc=1,NCS
         f_diam(itrc)=Sd50(itrc,ng) 
         f_vol(itrc)=pi/6.0_r8*(f_diam(itrc))**3.0_r8
         f_rho(itrc)=rhoref+(2650.0_r8-rhoref)*                          &
     &     (f_dp0/f_diam(itrc))**(3.0_r8-f_nf)
         f_mass(itrc)=f_vol(itrc)*(f_rho(itrc)-rhoref)
      ENDDO
      f_mass(NCS+1)=f_mass(NCS)*2.0_r8+1.0_r8  
      IF (f_diam(1).eq.f_dp0)  THEN
          f_mass(1)=f_vol(1)*Srho(1,ng)
      ENDIF
!TODO - This wont parallelize 
      WRITE(*,*) ' '
      WRITE(*,*) '*** FLOCMOD INIT *** '
      write(*,*) 'NAT, NPT, NCS, NNS:', NAT,NPT,NCS,NNS   
      WRITE(*,*) 'class diameter (um)  volume (m3)  density (kg/m3)  mass (kg)'
      DO itrc=1,NCS
         WRITE(*,*) itrc,f_diam(itrc)*1e6,f_vol(itrc),f_rho(itrc),f_mass(itrc)
      ENDDO
      write(*,*) 'f_mass(0) and f_mass(NCS+1): ',f_mass(0),f_mass(NCS+1)
      WRITE(*,*) ' '
      WRITE(*,*) ' *** PARAMETERS ***'
      WRITE(*,*) 'Primary particle size (f_dp0)                                : ',f_dp0
      WRITE(*,*) 'Fractal dimension (f_nf)                                     : ',f_nf
      WRITE(*,*) 'Flocculation efficiency (f_alpha)                            : ',f_alpha
      WRITE(*,*) 'Floc break up parameter (f_beta)                             : ',f_beta
      WRITE(*,*) 'Nb of fragments (f_nb_frag)                                  : ',f_nb_frag
      WRITE(*,*) 'Ternary fragmentation (f_ater)                               : ',f_ater
      WRITE(*,*) 'Floc erosion (% of mass) (f_ero_frac)                        : ',f_ero_frac
      WRITE(*,*) 'Nb of fragments by erosion (f_ero_nbfrag)                    : ',f_ero_nbfrag
      WRITE(*,*) 'fragment class (f_ero_iv)                                    : ',f_ero_iv
!      WRITE(*,*) 'negative mass tolerated before redistribution (f_mneg_param) : ',f_mneg_param
      WRITE(*,*) 'Boolean for differential settling aggregation (L_ADS)        : ',l_ADS
      WRITE(*,*) 'Boolean for shear aggregation (L_ASH)                        : ',l_ASH
      WRITE(*,*) 'Boolean for collision fragmenation (L_COLLFRAG)              : ',l_COLLFRAG
      WRITE(*,*) 'Collision fragmentation parameter (f_collfragparam)          : ',f_collfragparam
      WRITE(*,*) ' '
      WRITE(*,*) 'Value of eps                                                 : ',eps
      WRITE(*,*) ' '

! kernels computation former SUBROUTINE flocmod_kernels


!!--------------------------------------------------------------------------
!! * Executable part
      f_test=.true.
      dfragmax=0.00003_r8
! compute collision probability former SUBROUTINE flocmod_agregation_statistics

      DO iv1=1,NCS
        DO iv2=1,NCS

          f_coll_prob_sh(iv1,iv2)=1.0_r8/6.0_r8*(f_diam(iv1)+           &
     &            f_diam(iv2))**3.0_r8

          f_coll_prob_ds(iv1,iv2)=0.25_r8*pi*(f_diam(iv1)+              &
     &            f_diam(iv2))**2.0_r8*g/mu*abs((f_rho(iv1)-         &
     &            rhoref)*f_diam(iv1)**2.0_r8-(f_rho(iv2)-rhoref)*      &
     &            f_diam(iv2)**2.0_r8)

        ENDDO
      ENDDO

  !********************************************************************************
  ! agregation : GAIN : f_g1_sh and f_g1_ds
  !********************************************************************************
      DO iv1=1,NCS
       DO iv2=1,NCS
        DO iv3=iv2,NCS
           IF((f_mass(iv2)+f_mass(iv3)) .gt. f_mass(iv1-1)    .and.     &
     &            ((f_mass(iv2)+f_mass(iv3)) .le. f_mass(iv1))) THEN

              f_weight=(f_mass(iv2)+f_mass(iv3)-f_mass(iv1-1))/         &
     &            (f_mass(iv1)-f_mass(iv1-1))

           ELSEIF ((f_mass(iv2)+f_mass(iv3)) .gt. f_mass(iv1) .and.     &
     &            ((f_mass(iv2)+f_mass(iv3)) .lt. f_mass(iv1+1))) THEN

              IF (iv1 .eq. NCS) THEN
                 f_weight=1.0_r8
              ELSE
                 f_weight=1.0_r8-(f_mass(iv2)+f_mass(iv3)-              &
     &                 f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))
              ENDIF

           ELSE
              f_weight=0.0_r8
           ENDIF

           f_g1_sh(iv2,iv3,iv1)=f_weight*f_alpha*                       &
     &            f_coll_prob_sh(iv2,iv3)*(f_mass(iv2)+                 &
     &            f_mass(iv3))/f_mass(iv1)
           f_g1_ds(iv2,iv3,iv1)=f_weight*f_alpha*                       &
     &            f_coll_prob_ds(iv2,iv3)*(f_mass(iv2)+                 &
     &            f_mass(iv3))/f_mass(iv1)

        ENDDO
       ENDDO
      ENDDO

  !********************************************************************************
  ! Shear fragmentation : GAIN : f_g3
  !********************************************************************************

      DO iv1=1,NCS
       DO iv2=iv1,NCS

        IF (f_diam(iv2).gt.dfragmax) THEN
           ! binary fragmentation

           IF (f_mass(iv2)/f_nb_frag .gt. f_mass(iv1-1) &
                .and. f_mass(iv2)/f_nb_frag .le. f_mass(iv1)) THEN

              IF (iv1 .eq. 1) THEN 
                 f_weight=1.0_r8
              ELSE
                 f_weight=(f_mass(iv2)/f_nb_frag-f_mass(iv1-1))/        &
     &                (f_mass(iv1)-f_mass(iv1-1))
              ENDIF

           ELSEIF (f_mass(iv2)/f_nb_frag .gt. f_mass(iv1)               &
     &            .and. f_mass(iv2)/f_nb_frag .lt. f_mass(iv1+1)) THEN

              f_weight=1.0_r8-(f_mass(iv2)/f_nb_frag-f_mass(iv1))/      &
     &              (f_mass(iv1+1)-f_mass(iv1))

           ELSE

              f_weight=0.0_r8

           ENDIF

        ELSE
           f_weight=0.0_r8
        ENDIF

        f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1.0_r8-f_ero_frac)*(1.0_r8-        &
     &            f_ater)*f_weight*f_beta*f_diam(iv2)*((f_diam(iv2)-    &
     &            f_dp0)/f_dp0)**(3.0_r8-f_nf)*f_mass(iv2)/             &
     &            f_mass(iv1)

        ! ternary fragmentation
        IF (f_diam(iv2).gt.dfragmax) THEN
           IF (f_mass(iv2)/(2.0_r8*f_nb_frag) .gt. f_mass(iv1-1) .and.  &
     &          f_mass(iv2)/(2.0_r8*f_nb_frag) .le. f_mass(iv1)) THEN

              IF (iv1 .eq. 1) THEN 
                 f_weight=1.0_r8
              ELSE
                 f_weight=(f_mass(iv2)/(2.0_r8*f_nb_frag)-              &
     &                f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))
              ENDIF

           ELSEIF (f_mass(iv2)/(2.0_r8*f_nb_frag) .gt. f_mass(iv1)      &
     &           .and. f_mass(iv2)/(2.0_r8*f_nb_frag) .lt.              &
     &            f_mass(iv1+1)) THEN

              f_weight=1.0_r8-(f_mass(iv2)/(2.0_r8*f_nb_frag)-          &
     &                 f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))

           ELSE
              f_weight=0.0_r8

           ENDIF
           ! update for ternary fragments
           f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1.0_r8-f_ero_frac)*(f_ater)*    &
     &            f_weight*f_beta*f_diam(iv2)*((f_diam(iv2)-f_dp0)/     &
     &            f_dp0)**(3.0_r8-f_nf)*f_mass(iv2)/f_mass(iv1)   

           ! Floc erosion

           IF ((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .gt.         &
     &            f_mass(f_ero_iv)) THEN

              IF (((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .gt.     &
     &            f_mass(iv1-1)) .and. (f_mass(iv2)-f_mass(f_ero_iv)*   &
     &            f_ero_nbfrag) .le. f_mass(iv1)) THEN

                 IF (iv1 .eq. 1) THEN
                    f_weight=1.0_r8
                 ELSE
                    f_weight=(f_mass(iv2)-f_mass(f_ero_iv)*             &
     &                      f_ero_nbfrag-f_mass(iv1-1))/(f_mass(iv1)-   &
     &                      f_mass(iv1-1))
                 ENDIF

              ELSEIF ((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .gt.  &
     &            f_mass(iv1) .and. (f_mass(iv2)-f_mass(f_ero_iv)*      &
     &            f_ero_nbfrag) .lt. f_mass(iv1+1)) THEN

                 f_weight=1.0_r8-(f_mass(iv2)-f_mass(f_ero_iv)*         &
     &                    f_ero_nbfrag-f_mass(iv1))/(f_mass(iv1+1)-     &
     &                    f_mass(iv1))

              ELSE
                 f_weight=0.0_r8
              ENDIF

              ! update for eroded floc masses 

              f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*f_weight*f_beta*   &
     &            f_diam(iv2)*(max(eps,(f_diam(iv2)-f_dp0))/f_dp0)**    &
     &            (3.0_r8-f_nf)*(f_mass(iv2)-f_mass(f_ero_iv)*          &
     &            f_ero_nbfrag)/f_mass(iv1)

              IF (iv1 .eq. f_ero_iv) THEN

                 f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*f_beta*         &
     &              f_diam(iv2)*(max(eps,(f_diam(iv2)-f_dp0))/f_dp0)**  &
     &              (3.0_r8-f_nf)*f_ero_nbfrag*f_mass(f_ero_iv)/        &
     &              f_mass(iv1)
              ENDIF
           ENDIF
        ENDIF ! condition on dfragmax
       ENDDO
      ENDDO

  !********************************************************************************
  !  Shear agregation : LOSS : f_l1
  !********************************************************************************

      DO iv1=1,NCS
       DO iv2=1,NCS

        IF(iv2 .eq. iv1) THEN
           mult=2.0_r8
        ELSE
           mult=1.0_r8
        ENDIF

        f_l1_sh(iv2,iv1)=mult*f_alpha*f_coll_prob_sh(iv2,iv1) 
        f_l1_ds(iv2,iv1)=mult*f_alpha*f_coll_prob_ds(iv2,iv1) 

       ENDDO
      ENDDO

  !********************************************************************************
  !  Shear fragmentation : LOSS : f_l2
  !********************************************************************************

      DO iv1=1,NCS
       f_l3(iv1)=0.0_r8
       IF (f_diam(iv1).gt.dfragmax) THEN
        ! shear fragmentation
           f_l3(iv1)=f_l3(iv1)+(1.0_r8-f_ero_frac)*f_beta*f_diam(iv1)*  &
     &            ((f_diam(iv1)-f_dp0)/f_dp0)**(3.0_r8-f_nf)

        ! shear erosion
        IF ((f_mass(iv1)-f_mass(f_ero_iv)*f_ero_nbfrag) .gt.            &
     &            f_mass(f_ero_iv)) THEN
           f_l3(iv1)=f_l3(iv1)+f_ero_frac*f_beta*f_diam(iv1)*           &
     &            ((f_diam(iv1)-f_dp0)/f_dp0)**(3.0_r8-f_nf)
        ENDIF
       ENDIF
      ENDDO

      WRITE(*,*) ' '
      write(*,*) 'Sum of kernal coefficients:'
      write(*,*) 'f_coll_prob_sh',sum(f_coll_prob_sh)
      write(*,*) 'f_coll_prob_ds',sum(f_coll_prob_ds)
      write(*,*) 'f_g1_sh',sum(f_g1_sh)
      write(*,*) 'f_g1_ds',sum(f_g1_ds)
      write(*,*) 'f_l1_sh',sum(f_l1_sh)
      write(*,*) 'f_l1_ds',sum(f_l1_ds)
      write(*,*) 'f_g3',sum(f_g3)
      write(*,*) 'f_l3',sum(f_l3)
      WRITE(*,*) ' '
      WRITE(*,*) '*** END FLOCMOD INIT *** '    


!      END FORMER SUBROUTINE flocmod_kernels

      RETURN
      END SUBROUTINE initialize_sedflocs_param
