      SUBROUTINE ana_psource (ng, tile, model)
!
!! svn $Id: ana_psource.h 429 2009-12-20 17:30:26Z arango $
!!======================================================================
!! Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets analytical tracer and mass point Sources       !
!  and/or Sinks.  River runoff can be consider as a point source.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_sources
      USE mod_stepping
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_psource_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nnew(ng), knew(ng), Msrc(ng), Nsrc(ng),    &
     &                       OCEAN(ng) % zeta,                          &
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
#ifdef SOLVE3D
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       GRID(ng) % z_w,                            &
#endif
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % on_u,                           &
     &                       GRID(ng) % om_v,                           &
     &                       SOURCES(ng) % Isrc,                        &
     &                       SOURCES(ng) % Jsrc,                        &
     &                       SOURCES(ng) % Dsrc,                        &
#ifdef SOLVE3D
# if defined UV_PSOURCE || defined Q_PSOURCE
     &                       SOURCES(ng) % Qshape,                      &
     &                       SOURCES(ng) % Qsrc,                        &
# endif
# ifdef TS_PSOURCE
     &                       SOURCES(ng) % Tsrc,                        &
# endif
#endif
     &                       SOURCES(ng) % Qbar)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(20)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_psource
!
!***********************************************************************
      SUBROUTINE ana_psource_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nnew, knew, Msrc, Nsrc,              &
     &                             zeta, ubar, vbar,                    &
#ifdef SOLVE3D
     &                             u, v, z_w,                           &
#endif
     &                             h, on_u, om_v,                       &
     &                             Isrc, Jsrc, Dsrc,                    &
#ifdef SOLVE3D
# if defined UV_PSOURCE || defined Q_PSOURCE
     &                             Qshape, Qsrc,                        &
# endif
# ifdef TS_PSOURCE
     &                             Tsrc,                                &
# endif
#endif
     &                             Qbar)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
#ifdef SEDIMENT
      USE mod_sediment
#endif
#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
      USE distribute_mod, ONLY : mp_collect, mp_reduce
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nnew, knew
      integer, intent(in) :: Msrc

      integer, intent(out) :: Nsrc
!
#ifdef ASSUMED_SHAPE
      integer, intent(inout) :: Isrc(:)
      integer, intent(inout) :: Jsrc(:)

      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
# ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
# endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)

      real(r8), intent(inout) :: Dsrc(:)
      real(r8), intent(inout) :: Qbar(:)
# ifdef SOLVE3D
#  if defined UV_PSOURCE || defined Q_PSOURCE
      real(r8), intent(inout) :: Qshape(:,:)
      real(r8), intent(inout) :: Qsrc(:,:)
#  endif
#  ifdef TS_PSOURCE
      real(r8), intent(inout) :: Tsrc(:,:,:)
#  endif
# endif
#else
      integer, intent(inout) :: Isrc(Msrc)
      integer, intent(inout) :: Jsrc(Msrc)

      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
# ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
# endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: Dsrc(Msrc)
      real(r8), intent(inout) :: Qbar(Msrc)
# ifdef SOLVE3D
#  if defined UV_PSOURCE || defined Q_PSOURCE
      real(r8), intent(inout) :: Qshape(Msrc,N(ng))
      real(r8), intent(inout) :: Qsrc(Msrc,N(ng))
#  endif
#  ifdef TS_PSOURCE
      real(r8), intent(inout) :: Tsrc(Msrc,N(ng),NT(ng))
#  endif
# endif
#endif
!
!  Local variable declarations.
!
      integer :: Npts, NSUB, is, i, j, k, ised

      real(r8) :: Pspv = 0.0_r8
      real(r8), save :: area_east, area_west
      real(r8) :: cff, fac, my_area_east, my_area_west

#if defined DISTRIBUTE && defined SOLVE3D
      real(r8), dimension(Msrc*N(ng)) :: Pwrk
#endif
#if defined DISTRIBUTE
      real(r8), dimension(2) :: buffer

      character (len=3), dimension(2) :: io_handle
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set tracer and/or mass point sources and/or sink.
!-----------------------------------------------------------------------
!
      IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Set-up point Sources/Sink number (Nsrc), direction (Dsrc), I- and
!  J-grid locations (Isrc,Jsrc), and logical switch for type of tracer
!  to apply (LtracerSrc).  Currently, the direction can be along
!  XI-direction (Dsrc = 0) or along ETA-direction (Dsrc > 0).  The
!  mass sources are located at U- or V-points so the grid locations
!  should range from 1 =< Isrc =< L  and  1 =< Jsrc =< M.
!
#if defined RIVERPLUME1
        IF (Master.and.SOUTH_WEST_TEST) THEN
          Nsrc=1
          Dsrc(Nsrc)=0.0_r8
          Isrc(Nsrc)=1
          Jsrc(Nsrc)=50
          LtracerSrc(itemp,ng)=.TRUE.
          LtracerSrc(isalt,ng)=.TRUE.
        END IF
#elif defined RIVERPLUME2
        IF (Master.and.SOUTH_WEST_TEST) THEN
          Nsrc=1+Lm(ng)*2
          LtracerSrc(itemp,ng)=.TRUE.
          LtracerSrc(isalt,ng)=.TRUE.
          DO is=1,(Nsrc-1)/2
            Dsrc(is)=1.0_r8
            Isrc(is)=is
            Jsrc(is)=1
          END DO
          DO is=(Nsrc-1)/2+1,Nsrc-1
            Dsrc(is)=1.0_r8
            Isrc(is)=is-Lm(ng)
            Jsrc(is)=Mm(ng)+1
          END DO
          Dsrc(Nsrc)=0.0_r8
          Isrc(Nsrc)=1
          Jsrc(Nsrc)=60
        END IF
#elif defined SED_TEST1
        IF (Master.and.SOUTH_WEST_TEST) THEN
          Nsrc=Mm(ng)*2
          LtracerSrc(itemp,ng)=.TRUE.
          LtracerSrc(isalt,ng)=.TRUE.
          DO is=1,Nsrc/2
            Dsrc(is)=0.0_r8
            Isrc(is)=1
            Jsrc(is)=is
          END DO
          DO is=Nsrc/2+1,Nsrc
            Dsrc(is)=0.0_r8
            Isrc(is)=Lm(ng)+1
            Jsrc(is)=is-Mm(ng)
          END DO
        END IF
#else
        ana_psource.h: No values provided for LtracerSrc, Nsrc, Dsrc,
                                              Isrc, Jsrc.
#endif
#ifdef DISTRIBUTE
!
!  Broadcast point sources/sinks information to all nodes.
!
        CALL mp_bcasti (ng, iNLM, Nsrc)
        CALL mp_bcasti (ng, iNLM, Isrc)
        CALL mp_bcasti (ng, iNLM, Jsrc)
        CALL mp_bcastf (ng, iNLM, Dsrc)
#endif
      END IF

#if defined UV_PSOURCE || defined Q_PSOURCE
# ifdef SOLVE3D
!
!  If appropriate, set-up nondimensional shape function to distribute
!  mass point sources/sinks vertically.  It must add to unity!!.
!
#  ifdef DISTRIBUTE
      Qshape=Pspv
#  endif
      Npts=Msrc*N(ng)

!$OMP BARRIER

#  if defined SED_TEST1
      DO k=1,N(ng)
        DO is=1,Nsrc
          i=Isrc(is)
          j=Jsrc(is)
          IF (((IstrR.le.i).and.(i.le.IendR)).and.                      &
     &        ((JstrR.le.j).and.(j.le.JendR))) THEN
            IF (ubar(i,j,knew).ne.0.0_r8) THEN
              cff=ABS(u(i,j,k,nnew)/ubar(i,j,knew))
            ELSE
              cff=1.0_r8
            END IF
            Qshape(is,k)=cff*                                           &
     &                   (z_w(i-1,j,k    )-z_w(i-1,j,k-1)+              &
     &                    z_w(i  ,j,k    )-z_w(i  ,j,k-1))/             &
     &                   (z_w(i-1,j,N(ng))-z_w(i-1,j,0  )+              &
     &                    z_w(i  ,j,N(ng))-z_w(i  ,j,0  ))
          END IF
        END DO
      END DO
#   ifdef DISTRIBUTE
      Pwrk=RESHAPE(Qshape,(/Npts/))
      CALL mp_collect (ng, iNLM, Npts, Pspv, Pwrk)
      Qshape=RESHAPE(Pwrk,(/Msrc,N(ng)/))
#   endif
#  elif defined RIVERPLUME2
      DO k=1,N(ng)
        DO is=1,Nsrc-1
          i=Isrc(is)
          j=Jsrc(is)
          IF (((IstrR.le.i).and.(i.le.IendR)).and.                      &
     &        ((JstrR.le.j).and.(j.le.JendR))) THEN
            IF (vbar(i,j,knew).ne.0.0_r8) THEN
              cff=ABS(v(i,j,k,nnew)/vbar(i,j,knew))
            ELSE
              cff=1.0_r8
            END IF
            Qshape(is,k)=cff*                                           &
     &                   (z_w(i,j-1,k)-z_w(i,j-1,k-1)+                  &
     &                    z_w(i,j  ,k)-z_w(i,j  ,k-1))/                 &
     &                   (z_w(i,j-1,N(ng))-z_w(i,j-1,0)+                &
     &                    z_w(i,j  ,N(ng))-z_w(i,j  ,0))
          END IF
        END DO
      END DO
      IF (Master.and.SOUTH_WEST_TEST) THEN
        DO k=1,N(ng)
          Qshape(Nsrc,k)=1.0_r8/REAL(N(ng),r8)
        END DO
      END IF
#   ifdef DISTRIBUTE
      Pwrk=RESHAPE(Qshape,(/Npts/))
      CALL mp_collect (ng, iNLM, Npts, Pspv, Pwrk)
      Qshape=RESHAPE(Pwrk,(/Msrc,N(ng)/))
#   endif
#  else
!!
!!  Notice that there is not need for distributed-memory communications
!!  here since the computation below does not have a parallel tile
!!  dependency. All the nodes are computing this simple statement.
!!
      IF (NORTH_EAST_TEST) THEN
        DO k=1,N(ng)
          DO is=1,Nsrc
            Qshape(is,k)=1.0_r8/REAL(N(ng),r8)
          END DO
        END DO
      END IF
#  endif
# endif
!
!  Set-up vertically integrated mass transport (m3/s) of point
!  Sources/Sinks (positive in the positive U- or V-direction and
!  viceversa).
!
# ifdef DISTRIBUTE
      Qbar=Pspv
# endif

!$OMP BARRIER

# if defined RIVERPLUME1
      IF ((tdays(ng)-dstart).lt.0.5_r8) THEN
        fac=1.0_r8+TANH((time(ng)-43200.0_r8)/43200.0_r8)
      ELSE
        fac=1.0_r8
      END IF
      DO is=1,Nsrc
        Qbar(is)=fac*1500.0_r8
      END DO
# elif defined RIVERPLUME2
      DO is=1,(Nsrc-1)/2                     ! North end
        i=Isrc(is)
        j=Jsrc(is)
        IF (((IstrR.le.i).and.(i.le.IendR)).and.                        &
     &      ((JstrR.le.j).and.(j.le.JendR))) THEN
          Qbar(is)=-0.05_r8*om_v(i,j)*                                  &
     &             (0.5_r8*(zeta(i,j-1,knew)+h(i,j-1)+                  &
     &                      zeta(i  ,j,knew)+h(i  ,j)))
        END IF
      END DO
      DO is=(Nsrc-1)/2+1,Nsrc-1              ! South end
        i=Isrc(is)
        j=Jsrc(is)
        IF (((IstrR.le.i).and.(i.le.IendR)).and.                        &
     &      ((JstrR.le.j).and.(j.le.JendR))) THEN
          Qbar(is)=-0.05_r8*om_v(i,j)*                                  &
     &             (0.5_r8*(zeta(i,j-1,knew)+h(i,j-1)+                  &
     &                      zeta(i  ,j,knew)+h(i  ,j)))
        END IF
      END DO
      IF (Master.and.SOUTH_WEST_TEST) THEN
        Qbar(Nsrc)=1500.0_r8                 ! West wall
      END IF
#  ifdef DISTRIBUTE
      CALL mp_collect (ng, iNLM, Msrc, Pspv, Qbar)
#  endif
# elif defined SED_TEST1
      my_area_west=0.0_r8                    ! West end
      fac=-36.0_r8*10.0_r8*1.0_r8
      DO is=1,Nsrc/2
        i=Isrc(is)
        j=Jsrc(is)
        IF (((IstrR.le.i).and.(i.le.IendR)).and.                        &
     &      ((JstrR.le.j).and.(j.le.JendR))) THEN
          cff=0.5_r8*(zeta(i-1,j,knew)+h(i-1,j)+                        &
     &                zeta(i  ,j,knew)+h(i  ,j))*on_u(i,j)
          Qbar(is)=fac*cff
          my_area_west=my_area_west+cff
        END IF
      END DO
!
      my_area_east=0.0_r8                    ! East end
      fac=-36.0_r8*10.0_r8*1.0_r8
      DO is=Nsrc/2+1,Nsrc
        i=Isrc(is)
        j=Jsrc(is)
        IF (((IstrR.le.i).and.(i.le.IendR)).and.                        &
     &      ((JstrR.le.j).and.(j.le.JendR))) THEN
          cff=0.5_r8*(zeta(i-1,j,knew)+h(i-1,j)+                        &
     &                zeta(i  ,j,knew)+h(i  ,j))*on_u(i,j)
          Qbar(is)=fac*cff
          my_area_east=my_area_east+cff
        END IF
      END DO
!
      IF (SOUTH_WEST_CORNER.and.                                        &
     &    NORTH_EAST_CORNER) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
!$OMP CRITICAL (PSOURCE)
      IF (tile_count.eq.0) THEN
        area_west=0.0_r8
        area_east=0.0_r8
      END IF
      area_west=area_west+my_area_west
      area_east=area_east+my_area_east
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#  ifdef DISTRIBUTE
        buffer(1)=area_west
        buffer(2)=area_east
        io_handle(1)='SUM'
        io_handle(2)='SUM'
        CALL mp_reduce (ng, iNLM, 2, buffer, io_handle)
        area_west=buffer(1)
        area_east=buffer(1)
#  endif
        DO is=1,Nsrc/2
          Qbar(is)=Qbar(is)/area_west
        END DO
        DO is=Nsrc/2+1,Nsrc
          Qbar(is)=Qbar(is)/area_east
        END DO
      END IF
!$OMP END CRITICAL (PSOURCE)

#  ifdef DISTRIBUTE
      CALL mp_collect (ng, iNLM, Msrc, Pspv, Qbar)
#  endif
# else
      ana_psource.h: No values provided for Qbar.
# endif

# ifdef SOLVE3D
!
!  Set-up mass transport profile (m3/s) of point Sources/Sinks.
!
!$OMP BARRIER

      IF (NORTH_EAST_TEST) THEN
        DO k=1,N(ng)
          DO is=1,Nsrc
            Qsrc(is,k)=Qbar(is)*Qshape(is,k)
          END DO
        END DO
      END IF
# endif
#endif

#if defined TS_PSOURCE && defined SOLVE3D
!
!  Set-up tracer (tracer units) point Sources/Sinks.
!
# if defined RIVERPLUME1
      IF (NORTH_EAST_TEST) THEN
        DO k=1,N(ng)
          DO is=1,Nsrc
            Tsrc(is,k,itemp)=T0(ng)
            Tsrc(is,k,isalt)=0.0_r8
#  ifdef SEDIMENT
            DO ised=1,NST
              Tsrc(is,k,ised+2)=0.0_r8
            END DO
#  endif
          END DO
        END DO
      END IF
# elif defined RIVERPLUME2
      IF (NORTH_EAST_TEST) THEN
        DO k=1,N(ng)
          DO is=1,Nsrc-1
            Tsrc(is,k,itemp)=T0(ng)
            Tsrc(is,k,isalt)=S0(ng)
          END DO
          Tsrc(Nsrc,k,itemp)=T0(ng)
          Tsrc(Nsrc,k,isalt)=0.0_r8
        END DO
      END IF
# else
      ana_psource.h: No values provided for Tsrc.
# endif
#endif

      RETURN
      END SUBROUTINE ana_psource_tile
