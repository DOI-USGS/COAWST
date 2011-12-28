#include "cppdefs.h"
      MODULE mod_nesting
#ifdef COMPOSED_GRID_NOUSE
!
!svn $Id: mod_nesting.F 737 2008-09-07 02:06:44Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module defines generic structures used for nesting, composed,  !
!  and mosaic grids.                                                   !
!                                                                      !
!=======================================================================
!                                                 
        USE mod_kinds

        implicit none
!
! Generic structure for 2D state fields.
!
        TYPE T_NEST2D
          integer, pointer :: id             ! state field id
          integer, pointer :: gtype          ! staggered grid type
          integer, pointer :: tile(:)        ! tile partition

          integer, pointer :: LBi(:)         ! lower bound I-dimension
          integer, pointer :: UBi(:)         ! upper bound I-dimension
          integer, pointer :: LBj(:)         ! lower bound J-dimension
          integer, pointer :: UBj(:)         ! upper bound J-dimension

          integer, pointer :: Istr(:)        ! starting tile I-direction
          integer, pointer :: Iend(:)        ! ending   tile I-direction
          integer, pointer :: Jstr(:)        ! starting tile J-direction
          integer, pointer :: Jend(:)        ! ending   tile J-direction

          integer, pointer :: IstrR(:)       ! starting RHO  I-direction
          integer, pointer :: IstrU(:)       ! starting U    I-direction
          integer, pointer :: IendR(:)       ! ending   RHO  I-direction

          integer, pointer :: JstrR(:)       ! starting RHO  J-direction
          integer, pointer :: JstrV(:)       ! starting V    J-direction
          integer, pointer :: JendR(:)       ! ending   RHO  J-direction

          real (r8), pointer :: x(:,:)       ! X-positions
          real (r8), pointer :: y(:,:)       ! Y-positions
          real (r8), pointer :: angler(:,:)  ! grid angle
          real (r8), pointer :: ipos(:,:)    ! fractioanl X-positions
          real (r8), pointer :: jpos(:,:)    ! fractioanl X-positions
# ifdef MASKING
          real (r8), pointer :: mask(:,:)    ! land-sea masking
# endif
          real (r8), pointer :: s(:,:)       ! state array(i,j)
        END TYPE T_NEST2D

# ifdef SOLVE3D
!
! Generic structure for 3D state fields.
!
        TYPE T_NEST3D
          integer, pointer :: id             ! state field id
          integer, pointer :: gtype          ! staggered grid type
          integer, pointer :: tile(:)        ! tile partition

          integer, pointer :: LBi(:)         ! lower bound I-dimension
          integer, pointer :: UBi(:)         ! upper bound I-dimension
          integer, pointer :: LBj(:)         ! lower bound J-dimension
          integer, pointer :: UBj(:)         ! upper bound J-dimension
          integer, pointer :: LBk            ! lower bound K-dimension
          integer, pointer :: UBk            ! upper bound K-dimension

          integer, pointer :: Istr(:)        ! starting tile I-direction
          integer, pointer :: Iend(:)        ! ending   tile I-direction
          integer, pointer :: Jstr(:)        ! starting tile J-direction
          integer, pointer :: Jend(:)        ! ending   tile J-direction

          integer, pointer :: IstrR(:)       ! starting RHO  I-direction
          integer, pointer :: IstrU(:)       ! starting U    I-direction
          integer, pointer :: IendR(:)       ! ending   RHO  I-direction

          integer, pointer :: JstrR(:)       ! starting RHO  J-direction
          integer, pointer :: JstrV(:)       ! starting V    J-direction
          integer, pointer :: JendR(:)       ! ending   RHO  J-direction

          integer, pointer :: Kstr           ! starting K-index
          integer, pointer :: Kend           ! ending   K-index

          real (r8), pointer :: x(:,:)       ! X-positions
          real (r8), pointer :: y(:,:)       ! Y-positions
          real (r8), pointer :: angler(:,:)  ! grid angle
          real (r8), pointer :: ipos(:,:)    ! fractioanl X-positions
          real (r8), pointer :: jpos(:,:)    ! fractioanl X-positions
#  ifdef MASKING
          real (r8), pointer :: mask(:,:)    ! land-sea masking
#  endif
          real (r8), pointer :: s(:,:,:)     ! state array(i,j,k)
        END TYPE T_NEST3D
!
! Generic structure for 4D state fields.
!
        TYPE T_NEST4D
          integer, pointer :: id             ! state field id
          integer, pointer :: gtype          ! staggered grid type
          integer, pointer :: tile(:)        ! tile partition

          integer, pointer :: LBi(:)         ! lower bound I-dimension
          integer, pointer :: UBi(:)         ! upper bound I-dimension
          integer, pointer :: LBj(:)         ! lower bound J-dimension
          integer, pointer :: UBj(:)         ! upper bound J-dimension
          integer, pointer :: LBk            ! lower bound K-dimension
          integer, pointer :: UBk            ! upper bound K-dimension
          integer, pointer :: LBl            ! lower bound L-dimension
          integer, pointer :: UBl            ! upper bound L-dimension

          integer, pointer :: Istr(:)        ! starting tile I-direction
          integer, pointer :: Iend(:)        ! ending   tile I-direction
          integer, pointer :: Jstr(:)        ! starting tile J-direction
          integer, pointer :: Jend(:)        ! ending   tile J-direction

          integer, pointer :: IstrR(:)       ! starting RHO  I-direction
          integer, pointer :: IstrU(:)       ! starting U    I-direction
          integer, pointer :: IendR(:)       ! ending   RHO  I-direction

          integer, pointer :: JstrR(:)       ! starting RHO  J-direction
          integer, pointer :: JstrV(:)       ! starting V    J-direction
          integer, pointer :: JendR(:)       ! ending   RHO  J-direction

          integer, pointer :: Kstr           ! starting K-index
          integer, pointer :: Kend           ! ending   K-index
          integer, pointer :: Lstr           ! starting L-index
          integer, pointer :: Lend           ! ending   L-index

          real (r8), pointer :: x(:,:)       ! X-positions
          real (r8), pointer :: y(:,:)       ! Y-positions
          real (r8), pointer :: angler(:,:)  ! grid angle
          real (r8), pointer :: ipos(:,:)    ! fractioanl X-positions
          real (r8), pointer :: jpos(:,:)    ! fractioanl X-positions
#  ifdef MASKING
          real (r8), pointer :: mask(:,:)    ! land-sea masking
#  endif
          real (r8), pointer :: s(:,:,:,:)   ! state array(i,j,k,l)
        END TYPE T_NEST4D
# endif
#endif
      END MODULE mod_nesting

