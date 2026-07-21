#include "cppdefs.h"
      MODULE state_regrid_mod

#if defined SPLIT_4DVAR  && \
   (defined GRID_EXTRACT && defined SPLIT_EXECUTABLE)
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module regrids the state vector from source to target grid     !
!  resolution via interpolation.                                       !
!                                                                      !
!  Currently, it is used to remap the mixed resolution, split scheme   !
!  4D-Var coarser increment to the finer resolution of the outer loop  !
!  state vector in the 'analysis' phase.                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     tile       Domain partition (integer)                            !
!     model      Calling model identifier (integer)                    !
!     Linp       Coarse NetCDF time record to process (integer)        !
!     Lout       State array time index to update (integer)            !
!     ncname     Coarse resolution filename to process (string)        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_forces
      USE mod_iounits
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE state_assign_mod
!
      USE roms_interp_mod,  ONLY : roms2roms,                           &
     &                             roms2roms_create,                    &
     &                             roms_interp_regrid2d,                &
     &                             roms_interp_regrid_bylevel       
      USE strings_mod,      ONLY : FoundError
!
      implicit none
!
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      PUBLIC  :: frc_regrid
# endif
      PUBLIC  :: state_regrid
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      PRIVATE :: frc_regrid_tile
# endif
      PRIVATE :: state_regrid_tile
!
      PRIVATE
!
!-----------------------------------------------------------------------
!  Module local variables.
!-----------------------------------------------------------------------
!
!  Set ROMS-to-ROMS interpolation CLASS object, which it is used for
!  regridging.
!
      logical :: LcreatedC2F = .FALSE.    ! coarse-to-fine switch
!
      TYPE (roms2roms) :: C2F             ! coarse-to-fine interpolation
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CONTAINS
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS

!***********************************************************************
      SUBROUTINE frc_regrid (ng, tile, model, Linp, Lout, ncname)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: Linp, Lout
!
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, model, 88, __LINE__, MyFile)
#  endif
      CALL frc_regrid_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Linp, Lout,                                 &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
#  endif
#  ifdef ADJUST_WSTRESS
     &                      FORCES(ng) % ustr,                          &
     &                      FORCES(ng) % vstr,                          &
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr,                         &
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
     &                      FORCES(ng) % tflux,                         &
     &                      FORCES(ng) % stflx,                         &
#  endif
     &                      ncname)
#  ifdef PROFILE
      CALL wclock_off (ng, model, 88, __LINE__, MyFile)
#  endif
!
      RETURN
      END SUBROUTINE frc_regrid
!
!***********************************************************************
      SUBROUTINE frc_regrid_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Linp, Lout,                           &
#  ifdef MASKING
     &                            rmask, umask, vmask,                  &
#  endif
#  ifdef ADJUST_WSTRESS
     &                            ustr, vstr,                           &
     &                            nl_sustr, nl_svstr,                   &
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
     &                            tflux,                                &
     &                            nl_stflx,                             &
#  endif
     &                            ncname)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Linp, Lout
!
      character (len=*), intent(in) :: ncname
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in   ) :: rmask(LBi:,LBj:)
      real(r8), intent(in   ) :: umask(LBi:,LBj:)
      real(r8), intent(in   ) :: vmask(LBi:,LBj:)
#   endif
#   ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: vstr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: nl_svstr(LBi:,LBj:)
#   endif
#   if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tflux(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: nl_stflx(LBi:,LBj:,:)
#   endif

#  else

#   ifdef MASKING
      real(r8), intent(in   ) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: nl_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: nl_svstr(LBi:UBi,LBj:UBj)
#   endif
#   if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tflux(LBi:UBi,LBj:UBj,                    &
     &                                 Nfrec(ng),2,NT(ng))
      real(r8), intent(inout) :: nl_stflx(LBi:UBi,LBj:UBj,NT(ng))
#   endif
#  endif
!
!  Local variable declarations.
!
      logical, save :: FirstPass = .TRUE.
!
      integer :: gtype, i, it1, it2, j
      integer :: method
#  if defined ADJUST_STFLUX && defined SOLVE3D
      integer :: ifield, itrc
#  endif
!
      real(r8) :: fac, fac1, fac2
!
      character (len=1) :: Ctype(3)           ! C-grid locations: r,u,v
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", frc_regrid_tile"

#  include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Read in coarse grid surface adjustment increments from input NetCDF
!  and interpolate to finer resolution model grid. Then, adjust finer
!  grid surface focing.
!-----------------------------------------------------------------------
!
!  Intepolation scheme. 
!
      method = linear                        ! bi-linear
!
!  If first pass, interpolate surface adjutment increments for all
!  "Nfrc" time records snapshots.
!
      IF (FirstPass) THEN
        FirstPass=.FALSE.
!
!  If applicable, create ROMS-to-ROMS interpotation object: Coarse-to-
!  Fine state regridding. It is very likely, C2F object has been
!  already created.
!
      IF (.not.LcreatedC2F) THEN
        LcreatedC2F=.TRUE.
        Ctype(1)='r'
        Ctype(2)='u'
        Ctype(3)='v'
        CALL roms2roms_create (C2F, ng, tile, model, Ctype, ncname)
      END IF

#  ifdef ADJUST_WSTRESS
!
!  Regrid surface stress adjustments for all time levels.
!
      CALL roms_interp_regrid_bylevel (C2F%Cvars(2),                    &
     &                                 ng, tile, model, method, idUsms, &
     &                                 Linp, ncname, Vname(1,idUsms),   &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 1, Nfrec(ng),                    &
     &                                 Istr, IendR, JstrR, JendR,       &
#   ifdef MASKING
     &                                 umask,                           &
#   endif
     &                                 ustr(:,:,:,Lout))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL roms_interp_regrid_bylevel (C2F%Cvars(3),                    &
     &                                 ng, tile, model, method, idVsms, &
     &                                 Linp, ncname, Vname(1,idVsms),   &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 1, Nfrec(ng),                    &
     &                                 Istr, IendR, JstrR, JendR,       &
#   ifdef MASKING
     &                                 vmask,                           &
#   endif
     &                                 vstr(:,:,:,Linp))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Regrid surface fluxes adjustments for all time levels.
!
        gtype=r2dvar
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            ifield=idTsur(itrc)
            CALL roms_interp_regrid_bylevel (C2F%Cvars(1),              &
     &                                       ng, tile, model, method,   &
     &                                       ifield, Linp, ncname,      &
     &                                       Vname(1,ifield),           &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       1, Nfrec(ng),              &
     &                                       IstrR, IendR,              &
     &                                       JstrR, JendR,              &
#   ifdef MASKING
     &                                       rmask,                     &
#   endif
     &                                       tflux(:,:,:,Linp,itrc))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Set time records and interpolation factors, if any.
!-----------------------------------------------------------------------
!
      IF (Nfrec(ng).eq.1) THEN
        it1=1
        it2=1
        fac1=1.0_r8
        fac2=0.0_r8
      ELSE
#  ifdef GENERIC_DSTART
        it1=MAX(0,(iic(ng)-ntstart(ng))/nSFF(ng))+1
#  else
        it1=MAX(0,(iic(ng)-1)/nSFF(ng))+1
#  endif
        it2=MIN(it1+1,Nfrec(ng))
        fac1=SF_time(it2,ng)-(time(ng)+dt(ng))
        fac2=(time(ng)+dt(ng))-SF_time(it1,ng)
        fac=1.0_r8/(fac1+fac2)
        fac1=fac*fac1
        fac2=fac*fac2
      END IF

#  ifdef ADJUST_WSTRESS
!
!  Adjust surface momentum stress. Interpolate between surface forcing
!  increments, if appropriate.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          nl_sustr(i,j)=nl_sustr(i,j)+                                  &
     &                  fac1*ustr(i,j,it1,Lout)+                        &
     &                  fac2*ustr(i,j,it2,Lout)
        END DO
      END DO
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          nl_svstr(i,j)=nl_svstr(i,j)+                                  &
     &                  fac1*vstr(i,j,it1,Lout)+                        &
     &                  fac2*vstr(i,j,it2,Lout)
        END DO
      END DO
#  endif

#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Adjust surface tracer fluxes. Interpolate between surface forcing
!  increments, if appropriate.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              nl_stflx(i,j,itrc)=nl_stflx(i,j,itrc)+                    &
     &                           fac1*tflux(i,j,it1,Lout,itrc)+         &
     &                           fac2*tflux(i,j,it2,Lout,itrc)
            END DO
          END DO
        END IF
      END DO
#  endif
!
      RETURN
      END SUBROUTINE frc_regrid_tile
# endif
!
!***********************************************************************
      SUBROUTINE state_regrid (ng, tile, model, Linp, Lout, ncname)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: Linp, Lout
!
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 88, __LINE__, MyFile)
# endif
!
!  If not available, set finer and coarser grid structures.
!
      IF (SetGridConfig(ng)) THEN
        CALL set_grid (ng, model)
        SetGridConfig(ng)=.FALSE.
      END IF
!
!  Associate temporarily work state variables, my_var.
!
      CALL state_assign (ng, tile, model, Lout)
!
!  Regrid nonlinear model initial conditions,
!
      CALL state_regrid_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Linp, Lout, ncname,                       &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % vmask,                         &
# endif
# ifdef SOLVE3D
     &                        my_t,                                     &
     &                        my_u,                                     &
     &                        my_v,                                     &
# else
     &                        my_ubar,                                  &
     &                        my_vbar,                                  &
# endif
     &                        my_zeta)
!
!  Diassociate temporary work state variables, my_var.
!
      CALL state_unassign

# ifdef PROFILE
      CALL wclock_off (ng, model, 88, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE state_regrid
!
!***********************************************************************
      SUBROUTINE state_regrid_tile (ng, tile, model,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Linp, Lout, ncname,                 &
# ifdef MASKING
     &                              rmask, umask, vmask,                &
# endif
# ifdef SOLVE3D
     &                              t_regrid, u_regrid, v_regrid,       &
# else
     &                              ubar_regrid, vbar_regrid,           &
# endif
     &                              zeta_regrid)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Linp, Lout
!
      character (len=*), intent(in) :: ncname
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in   ) :: rmask(LBi:,LBj:)
      real(r8), intent(in   ) :: umask(LBi:,LBj:)
      real(r8), intent(in   ) :: vmask(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: t_regrid(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: u_regrid(LBi:,LBj:,:)
      real(r8), intent(inout) :: v_regrid(LBi:,LBj:,:)
#  else
      real(r8), intent(inout) :: ubar_regid(LBi:,LBj:)
      real(r8), intent(inout) :: vbar_regrd(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: zeta_regrid(LBi:,LBj:)

# else

#  ifdef MASKING
      real(r8), intent(in   ) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: t_regrid(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: u_regrid(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: v_regrid(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(inout) :: ubar_regrid(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vbar_regrid(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: zeta_regrid(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: gtype, i, j
# ifdef SOLVE3D
      integer :: ifield, itrc, k
# endif
      integer :: method
!
      character (len=1) :: Ctype(3)           ! C-grid locations: r,u,v
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ini_regrid_tile"

# include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Read in coarse grid state increment from input NetCDF and interpolate
!  to finer resolution model grid. Then, adjust finer grid initial
!  conditions for 4D-Var analysis phase.
!-----------------------------------------------------------------------
!
!  Initialize. 
!
      method = linear                        ! interpolation scheme
!
      IF (Master) WRITE (stdout,10) ng, method, Lout, Linp, TRIM(ncname)
!
!  Create ROMS-to-ROMS interpotation object: Coarse-to-Fine state
!  regridding.
!
      IF (.not.LcreatedC2F) THEN
        LcreatedC2F=.TRUE.
        Ctype(1)='r'
        Ctype(2)='u'
        Ctype(3)='v'
        CALL roms2roms_create (C2F, ng, tile, model, Ctype, ncname)
      END IF
!
!  Process free-surface from coarse to finer grid.
!
      zeta_regrid=0.0_r8
      CALL roms_interp_regrid2d (C2F%Cvars(1),                          &
     &                           ng, tile, model, method, idFsur,       &
     &                           Linp, ncname, Vname(1,idFsur),         &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IstrR, IendR, JstrR, JendR,            &
# ifdef MASKING
     &                           rmask,                                 &
# endif
     &                           zeta_regrid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifndef SOLVE3D
!
!  Process 2D U-momentum coarse increment and add it to finner prior
!  initial conditons.
!
      ubar_regrid=0.0_r8
      CALL roms_interp_regrid2d (C2F%Cvars(2),                          &
     &                           ng, tile, model, method, idUbar,       &
     &                           Linp, ncname, Vname(1,idUbar),         &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Istr, IendR, JstrR, JendR,             &
#  ifdef MASKING
     &                           umask,                                 &
#  endif
     &                           ubar_regrid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Process 2D V-momentum coarse increment and add it to finer  prior
!  initial conditons.
!
      vbar_regrid=0.0_r8
      CALL roms_interp_regrid2d (C2F%Cvars(3),                          &
     &                           ng, tile, model, method, idVbar,       &
     &                           Linp, ncname, Vname(1,idVbar),         &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IstrR, IendR, Jstr, JendR,             &
#  ifdef MASKING
     &                           vmask,                                 &
#  endif
     &                           vbar_rgrid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SOLVE3D
!
!  Process 3D U-momentum coarse increment and add it to finer prior
!  initial conditons. The regridding is done level-by-level.
!
      u_regrid=0.0_r8
      CALL roms_interp_regrid_bylevel (C2F%Cvars(2),                    &
     &                                 ng, tile, model, method, idUvel, &
     &                                 Linp, ncname, Vname(1,idUvel),   &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 Istr, IendR, JstrR, JendR,       &
#  ifdef MASKING
     &                                 umask,                           &
#  endif
     &                                 u_regrid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Process 3D V-momentum coarse increment and add it to finer prior
!  initial conditons. The regridding is done level-by-level.
!
      v_regrid=0.0_r8
      CALL roms_interp_regrid_bylevel (C2F%Cvars(3),                    &
     &                                 ng, tile, model, method, idVvel, &
     &                                 Linp, ncname, Vname(1,idVvel),   &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 IstrR, IendR, Jstr, JendR,       &
#  ifdef MASKING
     &                                 vmask,                           &
#  endif
     &                                 v_regrid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Process tracer-type variables coarse increment and add it to finer
!  prior initial conditions. The regridding is done level-by-level.
!
      t_regrid=0.0_r8
      DO itrc=1,NT(ng)
        ifield=idTvar(itrc)
        CALL roms_interp_regrid_bylevel (C2F%Cvars(1), ng,              &
     &                                   tile, model, method, ifield,   &
     &                                   Linp, ncname, Vname(1,ifield), &
     &                                   LBi, UBi, LBj, UBj, 1, N(ng),  &
     &                                   IstrR, IendR, JstrR, JendR,    &
#  ifdef MASKING
     &                                   rmask,                         &
#  endif
     &                                   t_regrid(:,:,:,itrc))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
# endif
!
 10   FORMAT (/,2x,'STATE_REGRID     - ROMS-to-ROMS State Vector ',     &
     &        'Regridding',/,22x,'(Grid ',i2.2,', Method=',i0,          &
     &        ', Index= ',i0,', Rec= ',i4.4,', Source: ',a,')')
!
      RETURN
      END SUBROUTINE state_regrid_tile
#endif
      END MODULE state_regrid_mod

