!
!svn $Id: sedbed_mod.h 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group        John C. Warner   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Sediment Model Kernel Variables:                                    !
!                                                                      !
#if defined BEDLOAD     && \
    defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
! avgbedldu       Time-averaged Bed load u-transport (kg/m/s).         !
! avgbedldv       Time-averaged Bed load v-transport (kg/m/s).         !
#endif
#ifdef SEDIMENT
!  bed            Sediment properties in each bed layer:               !
!                   bed(:,:,:,ithck) => layer thickness                !
!                   bed(:,:,:,iaged) => layer age                      !
!                   bed(:,:,:,iporo) => layer porosity                 !
!                   bed(:,:,:,idiff) => layer bio-diffusivity          !
!  bed_frac       Sediment fraction of each size class in each bed     !
!                   layer(nondimensional: 0-1.0).  Sum of              !
!                   bed_frac = 1.0.                                    !
!  bed_mass       Sediment mass of each size class in each bed layer   !
!                   (kg/m2).
#endif
#if defined SEDIMENT && defined SED_MORPH
!  bed_thick0     Sum all initial bed layer thicknesses (m).           !
!  bed_thick      Instantaneous total bed thickness (m).               !
#endif
#ifdef BEDLOAD
!  bedldu         Bed load u-transport (kg/m/s).                       !
!  bedldv         Bed load v-transport (kg/m/s).                       !
#endif
!  bottom         Exposed sediment layer properties:                   !
!                   bottom(:,:,isd50) => mean grain diameter           !
!                   bottom(:,:,idens) => mean grain density            !
!                   bottom(:,:,iwsed) => mean settling velocity        !
!                   bottom(:,:,itauc) => mean critical erosion stress  !
!                   bottom(:,:,irlen) => ripple length                 !
!                   bottom(:,:,irhgt) => ripple height                 !
!                   bottom(:,:,ibwav) => bed wave excursion amplitude  !
!                   bottom(:,:,izNik) => Nikuradse bottom roughness    !
!                   bottom(:,:,izbio) => biological bottom roughness   !
!                   bottom(:,:,izbfm) => bed form bottom roughness     !
!                   bottom(:,:,izbld) => bed load bottom roughness     !
!                   bottom(:,:,izapp) => apparent bottom roughness     !
!                   bottom(:,:,izwbl) => wave bottom roughness         !
!                   bottom(:,:,izdef) => default bottom roughness      !
!                   bottom(:,:,iactv) => active layer thickness        !
!                   bottom(:,:,ishgt) => saltation height              !
#if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
!                   bottom(:,:,idoff) => tau critical offset           !
!                   bottom(:,:,idslp) => tau critical slope            !
!                   bottom(:,:,idtim) => erodibility time scale        !
!                   bottom(:,:,idbmx) => diffusivity db_max            !
!                   bottom(:,:,idbmm) => diffusivity db_m              !
!                   bottom(:,:,idbzs) => diffusivity db_zs             !
!                   bottom(:,:,idbzm) => diffusivity db_zm             !
!                   bottom(:,:,idbzp) => diffusivity db_zphi           !
#endif
#if defined MIXED_BED
!                   bottom(:,:,idprp) => cohesive behavior             !
#endif
#if defined SEAGRASS_BOTTOM
!                   bottom(:,:,isgrH) => Seagrass height               !
!                   bottom(:,:,isgrD) => Seagrass shoot density        !
#endif
#if defined SEDIMENT && defined SUSPLOAD
!  ero_flux       Flux from erosion.                                   !
!  settling_flux  Flux from settling.                                  !
#endif
#if defined COHESIVE_BED || defined MIXED_BED
!  tcr_min         minimum shear for erosion
!  tcr_max         maximum shear for erosion
!  tcr_slp         Tau_crit profile slope
!  tcr_off         Tau_crit profile offset
!  tcr_tim         Tau_crit consolidation rate
#endif
#if defined MIXED_BED
!  transC          cohesive transition
!  transN          noncohesive transition
#endif
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      TYPE T_SEDBED
!
!  Nonlinear model state.
!
#if defined BEDLOAD     && \
    defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
        real(r8), pointer :: avgbedldu(:,:,:)
        real(r8), pointer :: avgbedldv(:,:,:)
#endif
#ifdef SEDIMENT
        real(r8), pointer :: bed(:,:,:,:)
        real(r8), pointer :: bed_frac(:,:,:,:)
        real(r8), pointer :: bed_mass(:,:,:,:,:)
#endif
#if defined SEDIMENT && defined SED_MORPH
        real(r8), pointer :: bed_thick0(:,:)
        real(r8), pointer :: bed_thick(:,:,:)
#endif
#ifdef BEDLOAD
        real(r8), pointer :: bedldu(:,:,:)
        real(r8), pointer :: bedldv(:,:,:)
#endif
        real(r8), pointer :: bottom(:,:,:)
#if defined SEDIMENT && defined SUSPLOAD
        real(r8), pointer :: ero_flux(:,:,:)
        real(r8), pointer :: settling_flux(:,:,:)
#endif
#if defined SEDIMENT && defined SED_BIOMASS
        real(r8), pointer :: Dstp_max(:,:,:)
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
# ifdef SEDIMENT
        real(r8), pointer :: tl_bed(:,:,:,:)
        real(r8), pointer :: tl_bed_frac(:,:,:,:)
        real(r8), pointer :: tl_bed_mass(:,:,:,:,:)
# endif
# if defined SEDIMENT && defined SED_MORPH
        real(r8), pointer :: tl_bed_thick0(:,:)
        real(r8), pointer :: tl_bed_thick(:,:,:)
# endif
# ifdef BEDLOAD
        real(r8), pointer :: tl_bedldu(:,:,:)
        real(r8), pointer :: tl_bedldv(:,:,:)
# endif
        real(r8), pointer :: tl_bottom(:,:,:)
# if defined SEDIMENT && defined SUSPLOAD
        real(r8), pointer :: tl_ero_flux(:,:,:)
        real(r8), pointer :: tl_settling_flux(:,:,:)
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
# ifdef SEDIMENT
        real(r8), pointer :: ad_bed(:,:,:,:)
        real(r8), pointer :: ad_bed_frac(:,:,:,:)
        real(r8), pointer :: ad_bed_mass(:,:,:,:,:)
# endif
# if defined SEDIMENT && defined SED_MORPH
        real(r8), pointer :: ad_bed_thick0(:,:)
        real(r8), pointer :: ad_bed_thick(:,:,:)
# endif
# ifdef BEDLOAD
        real(r8), pointer :: ad_bedldu(:,:,:)
        real(r8), pointer :: ad_bedldv(:,:,:)
# endif
        real(r8), pointer :: ad_bottom(:,:,:)
# if defined SEDIMENT && defined SUSPLOAD
        real(r8), pointer :: ad_ero_flux(:,:,:)
        real(r8), pointer :: ad_settling_flux(:,:,:)
# endif
#endif

      END TYPE T_SEDBED

      TYPE (T_SEDBED), allocatable :: SEDBED(:)

      CONTAINS

      SUBROUTINE allocate_sedbed (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
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
      IF (ng.eq.1) allocate ( SEDBED(Ngrids) )
!
!  Nonlinear model state.
!
#if defined BEDLOAD     && \
    defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
      IF (ANY(Aout(idUbld(:),ng))) THEN
        allocate ( SEDBED(ng) % avgbedldu(LBi:UBi,LBj:UBj,NST) )
      END IF
      IF (ANY(Aout(idVbld(:),ng))) THEN
        allocate ( SEDBED(ng) % avgbedldv(LBi:UBi,LBj:UBj,NST) )
      END IF
#endif
#if defined SEDIMENT
      allocate ( SEDBED(ng) % bed(LBi:UBi,LBj:UBj,Nbed,MBEDP) )
      allocate ( SEDBED(ng) % bed_frac(LBi:UBi,LBj:UBj,Nbed,NST) )
      allocate ( SEDBED(ng) % bed_mass(LBi:UBi,LBj:UBj,Nbed,2,NST) )
#endif
#if defined SEDIMENT && defined SED_MORPH
      allocate ( SEDBED(ng) % bed_thick0(LBi:UBi,LBj:UBj) )
      allocate ( SEDBED(ng) % bed_thick(LBi:UBi,LBj:UBj,1:2) )
#endif
#ifdef BEDLOAD
      allocate ( SEDBED(ng) % bedldu(LBi:UBi,LBj:UBj,NST) )
      allocate ( SEDBED(ng) % bedldv(LBi:UBi,LBj:UBj,NST) )
#endif
      allocate ( SEDBED(ng) % bottom(LBi:UBi,LBj:UBj,MBOTP) )
#if defined SEDIMENT && defined SUSPLOAD
      allocate ( SEDBED(ng) % ero_flux(LBi:UBi,LBj:UBj,NST) )
      allocate ( SEDBED(ng) % settling_flux(LBi:UBi,LBj:UBj,NST) )
#endif
#if defined SEDIMENT && defined SED_BIOMASS
      allocate ( SEDBED(ng) % Dstp_max(LBi:UBi,LBj:UBj,24) )
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
# if defined SEDIMENT
      allocate ( SEDBED(ng) % tl_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP) )
      allocate ( SEDBED(ng) % tl_bed_frac(LBi:UBi,LBj:UBj,Nbed,NST) )
      allocate ( SEDBED(ng) % tl_bed_mass(LBi:UBi,LBj:UBj,Nbed,2,NST) )
# endif
# if defined SEDIMENT && defined SED_MORPH
      allocate ( SEDBED(ng) % tl_bed_thick0(LBi:UBi,LBj:UBj) )
      allocate ( SEDBED(ng) % tl_bed_thick(LBi:UBi,LBj:UBj,1:2) )
# endif
# ifdef BEDLOAD
      allocate ( SEDBED(ng) % tl_bedldu(LBi:UBi,LBj:UBj,NST) )
      allocate ( SEDBED(ng) % tl_bedldv(LBi:UBi,LBj:UBj,NST) )
# endif
      allocate ( SEDBED(ng) % tl_bottom(LBi:UBi,LBj:UBj,MBOTP) )
# if defined SEDIMENT && defined SUSPLOAD
      allocate ( SEDBED(ng) % tl_ero_flux(LBi:UBi,LBj:UBj,NST) )
      allocate ( SEDBED(ng) % tl_settling_flux(LBi:UBi,LBj:UBj,NST) )
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
# if defined SEDIMENT
      allocate ( SEDBED(ng) % ad_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP) )
      allocate ( SEDBED(ng) % ad_bed_frac(LBi:UBi,LBj:UBj,Nbed,NST) )
      allocate ( SEDBED(ng) % ad_bed_mass(LBi:UBi,LBj:UBj,Nbed,2,NST) )
# endif
# if defined SEDIMENT && defined SED_MORPH
      allocate ( SEDBED(ng) % ad_bed_thick0(LBi:UBi,LBj:UBj) )
      allocate ( SEDBED(ng) % ad_bed_thick(LBi:UBi,LBj:UBj,1:2) )
# endif
# ifdef BEDLOAD
      allocate ( SEDBED(ng) % ad_bedldu(LBi:UBi,LBj:UBj,NST) )
      allocate ( SEDBED(ng) % ad_bedldv(LBi:UBi,LBj:UBj,NST) )
# endif
      allocate ( SEDBED(ng) % ad_bottom(LBi:UBi,LBj:UBj,MBOTP) )
# if defined SEDIMENT && defined SUSPLOAD
      allocate ( SEDBED(ng) % ad_ero_flux(LBi:UBi,LBj:UBj,NST) )
      allocate ( SEDBED(ng) % ad_settling_flux(LBi:UBi,LBj:UBj,NST) )
# endif
#endif

      RETURN
      END SUBROUTINE allocate_sedbed

      SUBROUTINE initialize_sedbed (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize structure variables in the module using     !
!  first touch distribution policy. In shared-memory configuration,    !
!  this operation actually performs the propagation of the  shared     !
!  arrays  across the cluster,  unless another policy is specified     !
!  to  override the default.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
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
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
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

#if defined BEDLOAD     && \
    defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
        IF (ANY(Aout(idUbld(:),ng))) THEN
          DO itrc=1,NST
            DO j=Jmin, Jmax
              DO i=Imin,Imax
                SEDBED(ng) % avgbedldu(i,j,itrc) = IniVal
              END DO
            END DO
          END DO
        END IF
        IF (ANY(Aout(idVbld(:),ng))) THEN
          DO itrc=1,NST
            DO j=Jmin, Jmax
              DO i=Imin,Imax
                SEDBED(ng) % avgbedldv(i,j,itrc) = IniVal
              END DO
            END DO
          END DO
        END IF
#endif
        DO j=Jmin,Jmax
#ifdef SEDIMENT
          DO itrc=1,MBEDP
            DO k=1,Nbed
              DO i=Imin,Imax
                SEDBED(ng) % bed(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
          DO itrc=1,NST
            DO k=1,Nbed
              DO i=Imin,Imax
                SEDBED(ng) % bed_frac(i,j,k,itrc) = IniVal
                SEDBED(ng) % bed_mass(i,j,k,1,itrc) = IniVal
                SEDBED(ng) % bed_mass(i,j,k,2,itrc) = IniVal
              END DO
            END DO
          END DO
#endif
#if defined SEDIMENT && defined SED_MORPH
          DO i=Imin,Imax
            SEDBED(ng) % bed_thick0(i,j) = IniVal
            SEDBED(ng) % bed_thick(i,j,1) = IniVal
            SEDBED(ng) % bed_thick(i,j,2) = IniVal
          END DO
#endif
#ifdef BEDLOAD
          DO itrc=1,NST
            DO i=Imin,Imax
              SEDBED(ng) % bedldu(i,j,itrc) = IniVal
              SEDBED(ng) % bedldv(i,j,itrc) = IniVal
            END DO
          END DO
#endif
          DO itrc=1,MBOTP
            DO i=Imin,Imax
              SEDBED(ng) % bottom(i,j,itrc) = IniVal
            END DO
          END DO
#if defined SEDIMENT && defined SUSPLOAD
          DO itrc=1,NST
            DO i=Imin,Imax
              SEDBED(ng) % ero_flux(i,j,itrc) = IniVal
              SEDBED(ng) % settling_flux(i,j,itrc) = IniVal
            END DO
          END DO
#endif
#if defined SEDIMENT && defined SED_BIOMASS
          DO itrc=1,24
            DO i=Imin,Imax
              SEDBED(ng) % Dstp_max(i,j,itrc) = 0.1_r8
            END DO
          END DO
#endif
        END DO
      END IF

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      IF ((model.eq.0).or.(model.eq.iTLM).or.(model.eq.iRPM)) THEN
        DO j=Jmin,Jmax
# ifdef SEDIMENT
          DO itrc=1,MBEDP
            DO k=1,Nbed
              DO i=Imin,Imax
                SEDBED(ng) % tl_bed(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
          DO itrc=1,NST
            DO k=1,Nbed
              DO i=Imin,Imax
                SEDBED(ng) % tl_bed_frac(i,j,k,itrc) = IniVal
                SEDBED(ng) % tl_bed_mass(i,j,k,1,itrc) = IniVal
                SEDBED(ng) % tl_bed_mass(i,j,k,2,itrc) = IniVal
              END DO
            END DO
          END DO
# endif
# if defined SEDIMENT && defined SED_MORPH
          DO i=Imin,Imax
            SEDBED(ng) % tl_bed_thick0(i,j) = IniVal
            SEDBED(ng) % tl_bed_thick(i,j,1) = IniVal
            SEDBED(ng) % tl_bed_thick(i,j,2) = IniVal
          END DO
# endif
# ifdef BEDLOAD
          DO itrc=1,NST
            DO i=Imin,Imax
              SEDBED(ng) % tl_bedldu(i,j,itrc) = IniVal
              SEDBED(ng) % tl_bedldv(i,j,itrc) = IniVal
            END DO
          END DO
# endif
          DO itrc=1,MBOTP
            DO i=Imin,Imax
              SEDBED(ng) % tl_bottom(i,j,itrc) = IniVal
            END DO
          END DO
# if defined SEDIMENT && defined SUSPLOAD
          DO itrc=1,NST
            DO i=Imin,Imax
              SEDBED(ng) % tl_ero_flux(i,j,itrc) = IniVal
              SEDBED(ng) % tl_settling_flux(i,j,itrc) = IniVal
            END DO
          END DO
# endif
        END DO
      END IF
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
      IF ((model.eq.0).or.(model.eq.iADM)) THEN
        DO j=Jmin,Jmax
# ifdef SEDIMENT
          DO itrc=1,MBEDP
            DO k=1,Nbed
              DO i=Imin,Imax
                SEDBED(ng) % ad_bed(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
          DO itrc=1,NST
            DO k=1,Nbed
              DO i=Imin,Imax
                SEDBED(ng) % ad_bed_frac(i,j,k,itrc) = IniVal
                SEDBED(ng) % ad_bed_mass(i,j,k,1,itrc) = IniVal
                SEDBED(ng) % ad_bed_mass(i,j,k,2,itrc) = IniVal
              END DO
            END DO
          END DO
# endif
# if defined SEDIMENT && defined SED_MORPH
          DO i=Imin,Imax
            SEDBED(ng) % ad_bed_thick0(i,j) = IniVal
            SEDBED(ng) % ad_bed_thick(i,j,1) = IniVal
            SEDBED(ng) % ad_bed_thick(i,j,2) = IniVal
          END DO
# endif
# ifdef BEDLOAD
          DO itrc=1,NST
            DO i=Imin,Imax
              SEDBED(ng) % ad_bedldu(i,j,itrc) = IniVal
              SEDBED(ng) % ad_bedldv(i,j,itrc) = IniVal
            END DO
          END DO
# endif
          DO itrc=1,MBOTP
            DO i=Imin,Imax
              SEDBED(ng) % ad_bottom(i,j,itrc) = IniVal
            END DO
          END DO
# if defined SEDIMENT && defined SUSPLOAD
          DO itrc=1,NST
            DO i=Imin,Imax
              SEDBED(ng) % ad_ero_flux(i,j,itrc) = IniVal
              SEDBED(ng) % ad_settling_flux(i,j,itrc) = IniVal
            END DO
          END DO
# endif
        END DO
      END IF
#endif

      RETURN
      END SUBROUTINE initialize_sedbed
