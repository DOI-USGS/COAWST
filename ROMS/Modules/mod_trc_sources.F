#include "cppdefs.h"
      MODULE mod_trc_sources
#if defined TRC_PSOURCE
!=======================================================================
!
        USE mod_kinds
        USE mod_param

        implicit none

        integer, allocatable :: Msrcpt(:)
        integer, allocatable :: Nsrcpt(:)

        TYPE T_TRC_SOURCES

          logical, pointer :: Lsrcpt(:,:)
          integer, pointer :: Isrcpt(:)
          integer, pointer :: Jsrcpt(:)
          real(r8), pointer :: Tsrcpt(:,:,:)

        END TYPE T_TRC_SOURCES

        TYPE (T_TRC_SOURCES), allocatable :: TRC_SOURCES(:)

      CONTAINS

      SUBROUTINE allocate_trc_sources (ng)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes all variables in the module  !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!  Imported variable declarations.
!
      integer :: ng
!
!  Local variable declarations.
!
      integer :: is, itrc, k, mg

      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
!  Set number of point sources to maximum number of analytical sources.
!  Notice that a maximum of 200 analytical sources are set-up here.
!
      allocate ( Msrcpt(Ngrids) )
      allocate ( Nsrcpt(Ngrids) )
      Msrcpt(ng)=200
      Nsrcpt(ng)=Msrcpt(ng)
!
!  Allocate structure.
!
      IF (ng.eq.1) allocate ( TRC_SOURCES(Ngrids) )
!
!  Allocate point Sources/Sinks variables.
!
      allocate ( TRC_SOURCES(ng) % Isrcpt(Nsrcpt(ng)) )
      allocate ( TRC_SOURCES(ng) % Jsrcpt(Nsrcpt(ng)) )
      allocate ( TRC_SOURCES(ng) % Lsrcpt(Nsrcpt(ng),NT(ng)) )
      allocate ( TRC_SOURCES(ng) % Tsrcpt(Nsrcpt(ng),N(ng),NT(ng)) )
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO is=1,Nsrcpt(ng)
        TRC_SOURCES(ng) % Isrcpt(is) = IniVal
        TRC_SOURCES(ng) % Jsrcpt(is) = IniVal
      END DO
      DO itrc=1,NT(ng)
        DO is=1,Nsrcpt(ng)
          TRC_SOURCES(ng) % Lsrcpt(is,itrc) = .FALSE.
        END DO
        DO k=1,N(ng)
          DO is=1,Nsrcpt(ng)
            TRC_SOURCES(ng) % Tsrcpt(is,k,itrc) = IniVal
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE allocate_trc_sources
#endif
      END MODULE mod_trc_sources
