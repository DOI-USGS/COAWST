#include "cppdefs.h"
      MODULE state_assign_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module associates/disassociates temporary work state variables !
!  "my_var" from ROMS kernel variables for efficient use of computer   !
!  memory.                                                             !
!                                                                      !
!  In the 4D-Var algorithms, several state/control vector copies are   !
!  needed to store ROMS trajectories, background, error covariance     !
!  standard deviations, error covariance normalization coefficients,   !
!  impulse forcing, eigenvectors, conjugate directions, and other      !
!  transformations. They are also used to read/write state vectors     !
!  into NetCDF files. The list below documents how these state vectors !
!  are used to help manage memory better in some parts of code.        !
!                                                                      !
!  STATE     ID    INDEX                                               !
!  --------  ----  -----                                               !
!  nl_state  iNLM  1,2     basic nonlinear kernel/prior state          !
!  ad_state  iADM  1,2     basic adjoint kernel/control state          !
!  tl_state  iTLM  1,2     basic tangent linear kernel/control state   !
!            iRPM  1,2     basic representer kernel/control state      !
!  b_state   14    1       initial conditions B normalization          !
!            15    2       model error B normalization                 !
!  e_state   10    1       initial conditions B standard deviations    !
!            11    1       model error B standard deviation            !
!  d_state   18    N/A     conjugate direction, eigenvectors           !
!  f_state   7     N/A     time interpolated impulse forcing           !
!  fG_state  19    1,2     impulse forcing snapshots                   !
!  fS_state  N/A   1,2     time convolutions snapshots (not used here) !
!                                                                      !
!  By default, the pointer sections remap bounds begin at index 1.     !
!  Therefore, it is necessary to explicitly specify all bounds,        !
!  especially for tiled array dimensions.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ocean
!
      implicit none
!
      PUBLIC :: state_assign
      PUBLIC :: state_unassign
!
!  Declare pointers for temporary work state variables.
!
      real (r8), pointer :: my_zeta(:,:)  => NULL()
      real (r8), pointer :: my_ubar(:,:)  => NULL()
      real (r8), pointer :: my_vbar(:,:)  => NULL()
#ifdef SOLVE3D
      real (r8), pointer :: my_t(:,:,:,:) => NULL()
      real (r8), pointer :: my_u(:,:,:)   => NULL()
      real (r8), pointer :: my_v(:,:,:)   => NULL()
#endif 
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CONTAINS
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     
!***********************************************************************
      SUBROUTINE state_assign (ng, tile, Istate, TimeIndex)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Istate, TimeIndex
!
!  Local variable declarations.
!
      integer :: is, ie, js, je, ks, ke, ls, le
!
!-----------------------------------------------------------------------
!  Assign temporarily work state variables to requested ROMS state
!  vector.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      is=BOUNDS(ng)%LBi(tile)
      ie=BOUNDS(ng)%UBi(tile)
      js=BOUNDS(ng)%LBj(tile)
      je=BOUNDS(ng)%UBj(tile)
      ks=1
      ke=N(ng)
      ls=1
      le=NT(ng)
!
!  Associate to required state vector variables.
!
      SELECT CASE (Istate)
        CASE (iNLM)
          IF (associated(OCEAN(ng)%zeta)) THEN
            my_zeta(is:,js:)  => OCEAN(ng)%zeta(is:ie,js:je,            &
     &                                          TimeIndex)
            my_ubar(is:,js:)  => OCEAN(ng)%ubar(is:ie,js:je,            &
     &                                          TimeIndex)
            my_vbar(is:,js:)  => OCEAN(ng)%vbar(is:ie,js:je,            &
     &                                          TimeIndex)
#ifdef SOLVE3D
            my_t(is:,js:,ks:,ls:) => OCEAN(ng)%t(is:ie,js:je,ks:ke,     &
     &                                           TimeIndex,ls:le)
            my_u(is:,js:,ks:)     => OCEAN(ng)%u(is:ie,js:je,ks:ke,     &
     &                                           TimeIndex)
            my_v(is:,js:,ks:)     => OCEAN(ng)%v(is:ie,js:je,ks:ke,     &
     &                                           TimeIndex)
#endif
            IF (Master) WRITE (stdout,10) 'Nonlinear', TimeIndex
          ELSE
            IF (Master) WRITE (stdout,30) 'Nonlinear'
          END IF

#if defined TANGENT || defined TL_IOMS
        CASE (iTLM, iRPM)
          IF (associated(OCEAN(ng)%tl_zeta)) THEN
            my_zeta(is:,js:)  => OCEAN(ng)%tl_zeta(is:ie,js:je,         &
     &                                             TimeIndex)
            my_ubar(is:,js:)  => OCEAN(ng)%tl_ubar(is:ie,js:je,         &
     &                                             TimeIndex)
            my_vbar(is:,js:)  => OCEAN(ng)%tl_vbar(is:ie,js:je,         &
     &                                             TimeIndex)
# ifdef SOLVE3D
            my_t(is:,js:,ks:,ls:) => OCEAN(ng)%tl_t(is:ie,js:je,ks:ke,  &
     &                                              TimeIndex,ls:le)
            my_u(is:,js:,ks:)     => OCEAN(ng)%tl_u(is:ie,js:je,ks:ke,  &
     &                                              TimeIndex)
            my_v(is:,js:,ks:)     => OCEAN(ng)%tl_v(is:ie,js:je,ks:ke,  &
     &                                              TimeIndex)
# endif
            IF (Master) WRITE (stdout,10) 'Tangent Linear', TimeIndex
          ELSE
            IF (Master) WRITE (stdout,30) 'Tangent Linear'
          END IF
#endif

#if defined ADJOINT
        CASE (iADM)
          IF (associated(OCEAN(ng)%ad_zeta)) THEN
            my_zeta(is:,js:)  => OCEAN(ng)%ad_zeta(is:ie,js:je,         &
     &                                             TimeIndex)
            my_ubar(is:,js:)  => OCEAN(ng)%ad_ubar(is:ie,js:je,         &
     &                                             TimeIndex)
            my_vbar(is:,js:)  => OCEAN(ng)%ad_vbar(is:ie,js:je,         &
     &                                             TimeIndex)
# ifdef SOLVE3D
            my_t(is:,js:,ks:,ls:) => OCEAN(ng)%ad_t(is:ie,js:je,ks:ke,  &
     &                                              TimeIndex,ls:le)
            my_u(is:,js:,ks:)     => OCEAN(ng)%ad_u(is:ie,js:je,ks:ke,  &
     &                                              TimeIndex)
            my_v(is:,js:,ks:)     => OCEAN(ng)%ad_v(is:ie,js:je,ks:ke,  &
     &                                              TimeIndex)
# endif
            IF (Master) WRITE (stdout,10) 'Adjoint', TimeIndex
          ELSE
            IF (Master) WRITE (stdout,30) 'Adjoint'
          END IF
#endif

#if defined FOUR_DVAR  || \
   (defined HESSIAN_SV && defined BNORM)
        CASE (14, 15)
          IF (associated(OCEAN(ng)%b_zeta)) THEN
            my_zeta(is:,js:)  => OCEAN(ng)%b_zeta(is:ie,js:je,        &
     &                                            TimeIndex)
            my_ubar(is:,js:)  => OCEAN(ng)%b_ubar(is:ie,js:je,        &
     &                                            TimeIndex)
            my_vbar(is:,js:)  => OCEAN(ng)%b_vbar(is:ie,js:je,        &
     &                                            TimeIndex)
# ifdef SOLVE3D
            my_t(is:,js:,ks:,ls:) => OCEAN(ng)%b_t(is:ie,js:je,ks:ke, &
     &                                             TimeIndex,ls:le)
            my_u(is:,js:,ks:)     => OCEAN(ng)%b_u(is:ie,js:je,ks:ke, &
     &                                             TimeIndex)
            my_v(is:,js:,ks:)     => OCEAN(ng)%b_v(is:ie,js:je,ks:ke, &
     &                                             TimeIndex)
# endif
            IF (Master) WRITE (stdout,10) 'b_state', TimeIndex
          ELSE
            IF (Master) WRITE (stdout,30) 'b_state'
          END IF
#endif

#if defined FOUR_DVAR  || \
   (defined HESSIAN_SV && defined BNORM)
        CASE (18)
          IF (associated(OCEAN(ng)%d_zeta)) THEN
            my_zeta => OCEAN(ng) % d_zeta
            my_ubar => OCEAN(ng) % d_ubar
            my_vbar => OCEAN(ng) % d_vbar
# ifdef SOLVE3D
            my_t    => OCEAN(ng) % d_t
            my_u    => OCEAN(ng) % d_u
            my_v    => OCEAN(ng) % d_v
# endif
            IF (Master) WRITE (stdout,20) 'd_state'
          ELSE
            IF (Master) WRITE (stdout,30) 'd_state'
          END IF
#endif

#if defined FOUR_DVAR  || \
   (defined HESSIAN_SV && defined BNORM)
        CASE (10, 11)
          IF (associated(OCEAN(ng)%e_zeta)) THEN
            my_zeta(is:,js:)  => OCEAN(ng) % e_zeta(is:ie,js:je,        &
     &                                              TimeIndex)
            my_ubar(is:,js:)  => OCEAN(ng) % e_ubar(is:ie,js:je,        &
     &                                              TimeIndex)
            my_vbar(is:,js:)  => OCEAN(ng) % e_vbar(is:ie,js:je,        &
     &                                              TimeIndex)
# ifdef SOLVE3D
            my_t(is:,js:,ks:,ls:) => OCEAN(ng) % e_t(is:ie,js:je,ks:ke, &
     &                                               TimeIndex,ls:le)
            my_u(is:,js:,ks:)     => OCEAN(ng) % e_u(is:ie,js:je,ks:ke, &
     &                                               TimeIndex)
            my_v(is:,js:,ks:)     => OCEAN(ng) % e_v(is:ie,js:je,ks:ke, &
     &                                               TimeIndex)
# endif
            IF (Master) WRITE (stdout,10) 'e_state', TimeIndex
          ELSE
            IF (Master) WRITE (stdout,30) 'e_state'
          END IF
#endif

#if defined WEAK_CONSTRAINT || defined IMPULSE || \
    defined FORCING_SV      || defined HESSIAN_FSV
        CASE (7)
          IF (associated(OCEAN(ng)%f_zeta)) THEN
            my_zeta => OCEAN(ng)%f_zeta
            my_ubar => OCEAN(ng)%f_ubar
            my_vbar => OCEAN(ng)%f_vbar
# ifdef SOLVE3D
            my_t    => OCEAN(ng)%f_t
            my_u    => OCEAN(ng)%f_u
            my_v    => OCEAN(ng)%f_v
# endif
            IF (Master) WRITE (stdout,30) 'f_state'
          ELSE
            IF (Master) WRITE (stdout,20) 'f_state'
          END IF
#endif

#if defined WEAK_CONSTRAINT || defined IMPULSE || \
    defined FORCING_SV      || defined HESSIAN_FSV
        CASE (19)
          IF (associated(OCEAN(ng)%f_zetaG)) THEN
            my_zeta(is:,js:)  => OCEAN(ng)%f_zetaG(is:ie,js:je,         &
     &                                             TimeIndex)
            my_ubar(is:,js:)  => OCEAN(ng)%f_ubarG(is:ie,js:je,         &
     &                                             TimeIndex)
            my_vbar(is:,js:)  => OCEAN(ng)%f_vbarG(is:ie,js:je,         &
     &                                             TimeIndex)
# ifdef SOLVE3D
            my_t(is:,js:,ks:,ls:) => OCEAN(ng)%f_tG(is:ie,js:je,ks:ke,  &
     &                                              TimeIndex,ls:le)
            my_u(is:,js:,ks:)     => OCEAN(ng)%f_uG(is:ie,js:je,ks:ke,  &
     &                                              TimeIndex)
            my_v(is:,js:,ks:)     => OCEAN(ng)%f_vG(is:ie,js:je,ks:ke,  &
     &                                              TimeIndex)
# endif
            IF (Master) WRITE (stdout,10) 'fG_state', TimeIndex
          ELSE
            IF (Master) WRITE (stdout,30) 'fG_state'
          END IF
#endif
      END SELECT
!
 10   FORMAT (/,2x,'STATE_ASSIGN     - Work state variables ',          &
     &        'pointers are associated with => ',a,' variables, ',      &
     &        'TimeIndex = ',i0)
 20   FORMAT (/,2x,'STATE_ASSIGN     - Work state variables ',          &
     &        'pointers are associated with => ',a,' variables.')
 30   FORMAT (/,2x,'STATE_ASSIGN     - unallocated ',a,                 &
     &        'state variables.')
!
      RETURN
      END SUBROUTINE state_assign
!     
!***********************************************************************
      SUBROUTINE state_unassign
!***********************************************************************
!
!-----------------------------------------------------------------------
!  Nullify temporary work state variables.
!-----------------------------------------------------------------------
!
      IF (associated(my_zeta)) nullify (my_zeta)
      IF (associated(my_ubar)) nullify (my_ubar)
      IF (associated(my_vbar)) nullify (my_vbar)
#ifdef SOLVE3D
      IF (associated(my_t))    nullify (my_t)
      IF (associated(my_u))    nullify (my_u)
      IF (associated(my_v))    nullify (my_v)
#endif
      IF (Master) WRITE(stdout,10)
 10   FORMAT (/,2x,'STATE_UNASSIGN   - Work state variables ',          &
     &        'pointers nullified.',/)
!
      RETURN
      END SUBROUTINE state_unassign
!
      END MODULE state_assign_mod
