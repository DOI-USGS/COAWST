#include "cppdefs.h"
      SUBROUTINE checkadj
!
!svn $Id: checkadj.F 791 2016-05-05 22:39:42Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine checks activated C-preprocessing options for        !
!  consistency with all available algorithms.                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_strings
!
      USE strings_mod, ONLY : uppercase
!
      implicit none
!
!  Local variable declarations.
!
      integer :: ic = 0
      integer :: ifound

      character (len=40) :: string
!
!-----------------------------------------------------------------------
!  Report issues with various C-preprocessing options.
!-----------------------------------------------------------------------
!
#ifdef BIO_FASHAM
      ic=ic+1
      string=uppercase('bio_fasham')
      IF (Master) WRITE(stdout,10) TRIM(string),                        &
     &            'CPP option renamed to '//uppercase('bio_fennel')
#endif

      string=uppercase('ts_smagorinsky')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        IF (Master) WRITE(stdout,10) TRIM(string),                      &
     &                               'stability problems, WARNING'
      END IF

      string=uppercase('ts_u3adv_split')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        IF (Master) WRITE(stdout,10) TRIM(string),                      &
     &                               'stability problems, WARNING'
      END IF

      string=uppercase('uv_smagorinsky')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        IF (Master) WRITE(stdout,10) TRIM(string),                      &
     &                               'stability problems, WARNING'
      END IF

      string=uppercase('uv_u3adv_split')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        IF (Master) WRITE(stdout,10) TRIM(string),                      &
     &                               'stability problems, WARNING'
      END IF

#if defined TANGENT || defined TL_IOMS || defined ADJOINT
!
!-----------------------------------------------------------------------
!  Stop if unsupported C-preprocessing options are activated for the
!  adjoint-based algorithms.
!-----------------------------------------------------------------------
!
      string=uppercase('atm_press')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not tested, FATAL ERROR'
      END IF

      string=uppercase('bedload_mpm')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('bedload_soulsby')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('bio_fennel')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('bio_umaine')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

# ifndef NL_BULK_FLUXES
      string=uppercase('bulk_fluxes')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
!!      ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not working, WARNING'
      END IF
# endif

# ifndef FORWARD_MIXING
      string=uppercase('bvf_mixing')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF
# endif

      string=uppercase('clima_ts_mix')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not tested, WARNING'
      END IF

      IF (ANY(ad_LBC(ieast,isFsur,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(ieast,isFsur)',            &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(iwest,isFsur,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(iwest,isFsur)',            &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(isouth,isFsur,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(isouth,isFsur)',           &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(inorth,isFsur,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(inorth,isFsur)',           &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(ieast,isUbar,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(ieast,isUbar)',            &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(iwest,isUbar,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(iwest,isUbar)',            &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(isouth,isUbar,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(isouth,isUbar)',           &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(inorth,isUbar,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(inorth,isUbar)',           &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(ieast,isVbar,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(ieast,isVbar)',            &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(iwest,isVbar,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(iwest,isVbar)',            &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(isouth,isVbar,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(isouth,isVbar)',           &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(inorth,isVbar,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(inorth,isVbar)',           &
     &                               'not finished, FATAL ERROR'
      END IF

# ifdef SOLVE3D
      IF (ANY(ad_LBC(ieast,isUvel,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(ieast,isUvel)',            &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(iwest,isUvel,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(iwest,isUvel)',            &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(isouth,isUvel,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(isouth,isUvel)',           &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(inorth,isUvel,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(inorth,isUvel)',           &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(ieast,isVvel,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(ieast,isVvel)',            &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(iwest,isVvel,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(iwest,isVvel)',            &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(isouth,isVvel,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(isouth,isVvel)',           &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(inorth,isVvel,:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(inorth,isVvel)',           &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(ieast,isTvar(:),:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(ieast,isTvar(:))',         &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(iwest,isTvar(:),:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(iwest,isTvar(:))',         &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(isouth,isTvar(:),:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(isouth,isTvar(:))',        &
     &                               'not finished, FATAL ERROR'
      END IF

      IF (ANY(ad_LBC(inorth,isTvar(:),:)%radiation)) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) 'ad_LBC(inorth,isTvar(:))',        &
     &                               'not finished, FATAL ERROR'
      END IF
# endif

      string=uppercase('ecosim')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

# ifndef FORWARD_MIXING
      string=uppercase('gls_mixing')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
!!      ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not differentiable, WARNING'
      END IF
# endif

# ifndef FORWARD_MIXING
      string=uppercase('lmd_mixing')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
!!      ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not differentiable, WARNING'
      END IF
# endif

      string=uppercase('mb_bbl')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

# ifndef FORWARD_MIXING
      string=uppercase('my25_mixing')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
!!      ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not differentiable, WARNING'
      END IF
# endif

      string=uppercase('nearshore_mellor')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('nemuro')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('nesting')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('npzd_franks')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
!!      ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not working, FATAL ERROR'
      END IF

      string=uppercase('pj_gradpq2')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('pj_gradpq4')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('red_tide')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('refdif_coupling')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not allowed, FATAL ERROR'
      END IF

      string=uppercase('sediment')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('sed_dens')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not tested, FATAL ERROR'
      END IF

      string=uppercase('sed_morph')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('sg_bbl')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('ssw_bbl')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('suspload')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('swan_coupling')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not allowed, FATAL ERROR'
      END IF

# if !(defined TS_HADVECTION_TL && defined TS_VADVECTION_TL)
      string=uppercase('ts_mpdata')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF
# endif

      string=uppercase('ts_smagorinsky')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

# if !defined TS_HADVECTION_TL
      string=uppercase('ts_u3adv_split')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF
# endif

      string=uppercase('uv_smagorinsky')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('uv_u3adv_split')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF

      string=uppercase('wet_dry')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20)  TRIM(string),                     &
     &                                'not coded, FATAL ERROR'
      END IF

      string=uppercase('wrf_coupling')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),                      &
     &                               'not coded, FATAL ERROR'
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Set execution error flag to stop execution.
!-----------------------------------------------------------------------
!
      IF (ic.gt.0) THEN
        exit_flag=5
      END IF
!
 10   FORMAT (/,' CHECKADJ - use caution when activating: ', a,/,12x,   &
     &        'REASON: ',a,'.')
#if defined TANGENT || defined TL_IOMS || defined ADJOINT
 20   FORMAT (/,' CHECKADJ - unsupported option in adjoint-based',      &
     &        ' algorithms: ',a,/,12x,'REASON: ',a,'.')
#endif

      RETURN
      END SUBROUTINE checkadj
