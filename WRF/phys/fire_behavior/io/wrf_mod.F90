  module wrf_mod

    use constants_mod, only : CP, XLV
    use emis_mod, only : Calc_smoke_aod
    use interp_mod, only : Interp_profile !, VINTERP_WINDS_FROM_3D_WINDS, VINTERP_WINDS_FROM_10M_WINDS
    use namelist_mod, only : namelist_t
    use proj_lc_mod, only : proj_lc_t
    use stderrout_mod, only : Stop_simulation

    implicit none

    private

    public :: Interp_wrf2dvar_to_cfbm, Interp_wrfwinds_to_cfbm, Provide_atm_feedback

  contains

    subroutine Interp_wrf2dvar_to_cfbm (wrfatm2dvar, ims, ime, jms, jme, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, &
        lats_in, lons_in, proj, vals_out)

      implicit none

      integer, intent (in) :: ims, ime, jms, jme, &
                              ifms, ifme, jfms, jfme, &
                              ifps, ifpe, jfps, jfpe
      real, dimension(ims:ime, jms:jme), intent (in) :: wrfatm2dvar
      real, dimension(ifms:ifme, jfms:jfme), intent (in) :: lats_in, lons_in
      type (proj_lc_t), intent (in) :: proj
      real, dimension(ifms:ifme, jfms:jfme), intent (in out) :: vals_out

      integer :: i, j, i_wrf, j_wrf
      real :: i_real, j_real


      do j = jfps, jfpe
        do i = ifps, ifpe
          call proj%Calc_ij (lats_in(i, j), lons_in(i, j), i_real, j_real)
          i_wrf = min (max (ims, nint (i_real)), ime)
          j_wrf = min (max (jms, nint (j_real)), jme)
          vals_out(i, j) = wrfatm2dvar(i_wrf, j_wrf)
        end do
      end do

    end subroutine Interp_wrf2dvar_to_cfbm

    subroutine Interp_wrfwinds_to_cfbm (u_phy, v_phy, z_at_w, ims, ime, kms, kme, jms, jme, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, &
        kfds, kfde, lats_in, lons_in, proj, z0f, fire_lsm_zcoupling, fire_lsm_zcoupling_ref, fire_wind_height, u_out, v_out)

      implicit none

      integer, intent (in) :: ims, ime, kms, kme, jms, jme, &
                              ifms, ifme, jfms, jfme, &
                              ifps, ifpe, jfps, jfpe, &
                              kfds, kfde
      logical, intent (in) :: fire_lsm_zcoupling
      real, dimension(ims:ime, kms:kme, jms:jme), intent (in) :: u_phy, v_phy, z_at_w
      real, dimension(ifms:ifme, jfms:jfme), intent (in) :: lats_in, lons_in, z0f
      type (proj_lc_t), intent (in) :: proj
      real, intent (in) :: fire_wind_height, fire_lsm_zcoupling_ref
      real, dimension(ifms:ifme, jfms:jfme), intent (in out) :: u_out, v_out

      integer :: i, j, i_wrf, j_wrf
      real :: i_real, j_real, uout, vout


      do j = jfps, jfpe
        do i = ifps, ifpe
          call proj%Calc_ij (lats_in(i, j), lons_in(i, j), i_real, j_real)
          i_wrf = min (max (ims, nint (i_real)), ime)
          j_wrf = min (max (jms, nint (j_real)), jme)
          call Interp_profile (fire_lsm_zcoupling, fire_lsm_zcoupling_ref, fire_wind_height, kfds, kfde, &
              u_phy(i_wrf, :, j_wrf), v_phy(i_wrf, :, j_wrf), z_at_w(i_wrf, :, j_wrf), z0f(i, j), &
              uout, vout)
          u_out(i, j) = uout
          v_out(i, j) = vout
        end do
      end do

    end subroutine  Interp_wrfwinds_to_cfbm

    subroutine Provide_atm_feedback (config_flags, &
            ifms, ifme, jfms, jfme,                &
            ifts, ifte, jfts, jfte,                &
            ifps, ifpe, jfps, jfpe,                &
            ids, ide, kds, kde, jds, jde,          &
            ims, ime, kms, kme, jms, jme,          &
            its, ite, kts, kte, jts, jte,          &
            sr_x, sr_y,                            &
            emis_smoke, smoke_tracer, tracer_opt,  &
            p_phy, t_phy, qv,                      &
            aod5502d_smoke,                        &
            fgrnhfx, fgrnqfx,                      &
            grnhfx, grnqfx, canhfx, canqfx,        &
            grnsmk,                                &
            alfg, alfc, z1can,                     &
            rho, dz8w, z_at_w,                     &
            mu, c1h, c2h,                          &
            rthfrten, rqvfrten)

      implicit none

      type (namelist_t), intent (in) :: config_flags
      integer, intent (in) :: ifms, ifme, jfms, jfme,       &
                              ifps, ifpe, jfps, jfpe,       &
                              ids, ide, kds, kde, jds, jde, &
                              ims, ime, kms, kme, jms, jme, &
                              its, ite, kts, kte, jts, jte, &
                              ifts, ifte, jfts, jfte, sr_x, sr_y

      real, dimension(ifms:ifme, jfms:jfme), intent (in) :: emis_smoke
      real, dimension(ims:ime, kms:kme, jms:jme), intent (in out), optional :: smoke_tracer
      integer, intent (in) :: tracer_opt
      real, dimension(ifms:ifme, jfms:jfme), intent (in) :: fgrnhfx, fgrnqfx
      real, dimension(ims:ime, kms:kme, jms:jme), intent (in) :: rho, dz8w, z_at_w, p_phy, t_phy, qv
      real, intent(in), dimension(ims:ime, jms:jme) :: mu   ! dry air mass (pa)
      real, intent(in), dimension(kms:kme) :: c1h, c2h      ! hybrid coordinate weights
      real, intent(in) :: alfg                              ! extinction depth surface fire heat (m)
      real, intent(in) :: alfc                              ! extinction depth crown  fire heat (m)
      real, intent(in) :: z1can                             ! height of crown fire heat release (m)
      real, dimension(ims:ime, jms:jme), intent (out) :: grnhfx, grnqfx, canhfx, canqfx, grnsmk, aod5502d_smoke
      real, intent(out), dimension(ims:ime, kms:kme, jms:jme) ::   &
           rthfrten, & ! theta tendency from fire (in mass units)
           rqvfrten    ! Qv tendency from fire (in mass units)

      logical, parameter :: DEBUG_LOCAL = .true.
      integer :: i, j, ibase, jbase, i_f, j_f, ioff, joff
      real :: avgw, convert_kg_m2_to_g_kg


      if (DEBUG_LOCAL) call Check_dims (its, ite, jts, jte, ifts, ifte, jfts, jfte, sr_x, sr_y)

      avgw = 1.0 / (sr_x * sr_y)
      do j = max (jds + 1, jts), min (jte, jde - 2)
        jbase = jfts + sr_y * (j - jts)
        do i = max (ids + 1, its), min (ite, ide - 2)
          ibase = ifts + sr_x * (i - its)
          canqfx(i, j) = 0.0
          canhfx(i, j) = 0.0
          grnsmk(i, j) = 0.0
          grnhfx(i, j) = 0.0
          grnqfx(i, j) = 0.0
          convert_kg_m2_to_g_kg = 1000.0 / (rho(i, kts, j) * dz8w(i, kts, j))
          do joff = 0, sr_y - 1
            j_f = joff + jbase
            do ioff = 0, sr_x - 1
              i_f = ioff + ibase
              grnsmk(i, j) = grnsmk(i, j) + emis_smoke(i_f, j_f)
              grnhfx(i, j) = grnhfx(i, j) + fgrnhfx(i_f, j_f) ! * config_flags%fire_atm_feedback
              grnqfx(i, j) = grnqfx(i, j) + fgrnqfx(i_f, j_f) ! * config_flags%fire_atm_feedback
            end do
          end do
          grnhfx(i, j) = grnhfx(i, j) * avgw
          grnqfx(i, j) = grnqfx(i, j) * avgw
          grnsmk(i, j) = grnsmk(i, j) * convert_kg_m2_to_g_kg * avgw
          if (tracer_opt == 3) smoke_tracer(i, kts, j) = smoke_tracer(i, kts, j) + grnsmk(i, j)
        end do
      end do

      call Fire_tendency (               &
            ids,ide - 1,kds,kde,jds,jde - 1,     & ! dimensions
            ims,ime,kms,kme,jms,jme,     &
            its,min (ite, ide-1),kts,kte,jts,min (jte, jde - 1),     &
            grnhfx,grnqfx,canhfx,canqfx, & ! heat fluxes summed up to  atm grid
            alfg,alfc,z1can,             & ! coeffients, properties, geometry
            z_at_w,dz8w,mu,c1h,c2h,rho,  &
            config_flags%fire_atm_feedback, &
            rthfrten,rqvfrten)             ! theta and Qv tendencies

      if (tracer_opt == 3) call Calc_smoke_aod (dz8w, p_phy, t_phy, qv, rho, smoke_tracer, aod5502d_smoke, &
           ids, ide, kds, kde, jds, jde,          &
           ims, ime, kms, kme, jms, jme,          &
           its, ite, kts, kte, jts, jte)

    contains

      subroutine Check_dims (its, ite, jts, jte, ifts, ifte, jfts, jfte, sr_x, sr_y)

        use, intrinsic :: iso_fortran_env, only : OUTPUT_UNIT

        implicit none

        integer, intent (in) :: its, ite, jts, jte, ifts, ifte, jfts, jfte, sr_x, sr_y

        integer :: isz1, jsz1, isz2, jsz2, ir, jr
        logical, parameter :: DEBUG_LOCAL = .false.


        isz1 = ite - its + 1
        jsz1 = jte - jts + 1
        isz2 = ifte - ifts + 1
        jsz2 = jfte - jfts + 1
        ir = isz2 / isz1
        jr = jsz2 / jsz1

        if (DEBUG_LOCAL) write (OUTPUT_UNIT, *) 'its, ite, jts, jte =', its, ite, jts, jte
        if (DEBUG_LOCAL) write (OUTPUT_UNIT, *) 'ifts, ifte, jfts, jfte =', ifts, ifte, jfts, jfte 
        if (DEBUG_LOCAL) write (OUTPUT_UNIT, *) 'isz1, jsz1, isz2, jsz2 =', isz1, jsz1, isz2, jsz2
        if (DEBUG_LOCAL) write (OUTPUT_UNIT, *) 'ir, jz =', ir, jr

        if (ir /= sr_x .or. jr /= sr_y) call Stop_simulation ('Tile dims do not preserve fire/atm ratio')

      end subroutine Check_dims

    end subroutine Provide_atm_feedback

    subroutine Fire_tendency(   &
        ids,ide, kds,kde, jds,jde,   & ! dimensions
        ims,ime, kms,kme, jms,jme,   &
        its,ite, kts,kte, jts,jte,   &
        grnhfx,grnqfx,canhfx,canqfx, & ! heat fluxes summed up to  atm grid
        alfg,alfc,z1can,             & ! coeffients, properties, geometry
        z_at_w,dz8w,mu,c1h,c2h,rho,  &
        fire_atm_feedback,           &
        rthfrten,rqvfrten)             ! theta and Qv tendencies

    ! This routine is atmospheric physics

    ! --- this routine takes fire generated heat and moisture fluxes and
    !     calculates their influence on the theta and water vapor
    ! --- note that these tendencies are valid at the Arakawa-A location

      use, intrinsic :: iso_fortran_env, only : OUTPUT_UNIT
      implicit none

    ! --- incoming variables

      integer, intent (in) :: ids, ide, kds, kde, jds, jde, &
                              ims, ime, kms, kme, jms, jme, &
                              its, ite, kts, kte, jts, jte

      real, intent(in), dimension(ims:ime, jms:jme) :: grnhfx,grnqfx  ! w/m^2
      real, intent(in), dimension(ims:ime, jms:jme) :: canhfx,canqfx  ! w/m^2
      real, intent(in), dimension(ims:ime, jms:jme) :: mu             ! dry air mass (pa)
      real, intent(in), dimension(kms:kme) :: c1h, c2h       ! hybrid coordinate weights

      real, intent(in), dimension(ims:ime, kms:kme, jms:jme) :: z_at_w ! m abv sealvl
      real, intent(in), dimension(ims:ime, kms:kme, jms:jme) :: dz8w   ! dz across w-lvl
      real, intent(in), dimension(ims:ime, kms:kme, jms:jme) :: rho    ! density

      real, intent(in) :: alfg     ! extinction depth surface fire heat (m)
      real, intent(in) :: alfc     ! extinction depth crown  fire heat (m)
      real, intent(in) :: z1can    ! height of crown fire heat release (m)

      real, intent(in) :: fire_atm_feedback

    ! --- outgoing variables

      real, intent(out), dimension(ims:ime, kms:kme, jms:jme) ::   &
           rthfrten, & ! theta tendency from fire (in mass units)
           rqvfrten    ! Qv tendency from fire (in mass units)
    ! --- local variables

      integer :: i,j,k
      integer :: i_st,i_en, j_st,j_en, k_st,k_en

      real :: cp_i
      real :: rho_i
      real :: xlv_i
      real :: z_w
      real :: fact_g, fact_c
      real :: alfg_i, alfc_i

      real, dimension( its:ite,kts:kte,jts:jte ) :: hfx,qfx


!      write (OUTPUT_UNIT, *) 'pajm: its, ite, jts, jte, kts, kte = ', its, ite, jts, jte, kts, kte
!      write (OUTPUT_UNIT, *) 'pajm: ids, ide, jds, jde, kds, kde = ', ids, ide, jds, jde, kds, kde
!      write (OUTPUT_UNIT, *) 'pajm: alfg,alfc,z1can = ', alfg,alfc,z1can

      do j=jts,jte
        do k=kts,min(kte+1,kde)
          do i=its,ite
            rthfrten(i,k,j)=0.
            rqvfrten(i,k,j)=0.
          enddo
        enddo
      enddo

      if (fire_atm_feedback <= 0.0) return

    ! --- set some local constants

      cp_i = 1.0 / CP     ! inverse of specific heat
      xlv_i = 1.0 / XLV   ! inverse of latent heat
      alfg_i = 1./alfg
      alfc_i = 1./alfc

    ! --- set loop indicies : note that

      i_st = MAX(its,ids+1)
      i_en = MIN(ite,ide-1)
      k_st = kts
      k_en = MIN(kte,kde-1)
      j_st = MAX(jts,jds+1)
      j_en = MIN(jte,jde-1)
    ! --- distribute fluxes

      do j = j_st,j_en
        do k = k_st,k_en
          do i = i_st,i_en

            ! --- set z (in meters above ground)

            z_w = z_at_w(i,k,j) - z_at_w(i, 1, j)

            ! --- heat flux

            fact_g = cp_i * EXP( - alfg_i * z_w )
            if ( z_w < z1can ) then
                   fact_c = cp_i
            else
                   fact_c = cp_i * EXP( - alfc_i * (z_w - z1can) )
            end if
            hfx(i,k,j) = fact_g * grnhfx(i,j) * fire_atm_feedback+ fact_c * canhfx(i,j)

            ! --- vapor flux

            fact_g = xlv_i * EXP( - alfg_i * z_w )
            if (z_w < z1can) then
                   fact_c = xlv_i
            else
                   fact_c = xlv_i * EXP( - alfc_i * (z_w - z1can) )
            end if
            qfx(i,k,j) = fact_g * grnqfx(i,j) * fire_atm_feedback + fact_c * canqfx(i,j)

!            if ((grnhfx(i,j) * fire_atm_feedback >0.) .and. (k == 1)) then
!              write (OUTPUT_UNIT, *) 'masih: grnhfx, grnqfx', grnhfx(i,j), grnqfx(i,j)
!              write (OUTPUT_UNIT, *) 'masih: hfx, qfx', hfx(i,1,j), qfx(i,1,j), hfx(i,1,j), qfx(i,1,j)
!            end if

          end do
        end do
      end do
    ! --- add flux divergence to tendencies
    !
    !   multiply by dry air mass (mu) to eliminate the need to
    !   call sr. calculate_phy_tend (in dyn_em/module_em.F)

      do j = j_st,j_en
        do k = k_st,k_en-1
          do i = i_st,i_en

            rho_i = 1./rho(i,k,j)

            rthfrten(i,k,j) = - (c1h(k)*mu(i,j)+c2h(k)) * rho_i * (hfx(i,k+1,j)-hfx(i,k,j)) / dz8w(i,k,j)
            rqvfrten(i,k,j) = - (c1h(k)*mu(i,j)+c2h(k)) * rho_i * (qfx(i,k+1,j)-qfx(i,k,j)) / dz8w(i,k,j)

          end do
        end do
      end do

      return

    end subroutine Fire_tendency

  end module wrf_mod
