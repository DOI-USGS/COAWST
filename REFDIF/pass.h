
!C ------------------------------------------------------------------
!c   Common blocks for interface of Nearshore Community Model
!c   It is used in master program, wave module, circulation module,
!c   and sediment module.
!c      Fyshi 01/21/2002   
!c ------------------------------------------------------------------

       integer Nx_Max, Ny_Max, Nx_Circ, Ny_Circ, Nx_Wave,Ny_Wave,       &
     &         Nx_Mast, Ny_Mast, Nx_Sedi, Ny_Sedi
  
       parameter (Nx_Max = 300, Ny_Max = 300) 

!c -- wave module:

       real Pass_Sxx(Nx_Max,Ny_Max), Pass_Sxy(Nx_Max,Ny_Max),           &
     &      Pass_Syy(Nx_Max,Ny_Max),                                    &
     &      Pass_Sxx_body(Nx_Max,Ny_Max),Pass_Sxy_body(Nx_Max,Ny_Max),  &
     &      Pass_Syy_body(Nx_Max,Ny_Max),                               &
     &      Pass_Sxx_surf(Nx_Max,Ny_Max),                               &
     &      Pass_Sxy_surf(Nx_Max,Ny_Max),                               &
     &      Pass_Syy_surf(Nx_Max,Ny_Max),                               &
     &      Pass_Wave_Fx(Nx_Max,Ny_Max),Pass_Wave_Fy(Nx_Max,Ny_Max),    &
     &      Pass_MassFluxU(Nx_Max,Ny_Max),                              &
     &      Pass_MassFluxV(Nx_Max,Ny_Max),                              &
     &      Pass_MassFlux(Nx_Max,Ny_Max),                               &
     &      Pass_Diss(Nx_Max,Ny_Max),                                   &
     &      Pass_WaveNum(Nx_Max,Ny_Max), Pass_Theta(Nx_Max,Ny_Max),     &
     &      Pass_ubott(Nx_Max,Ny_Max), Pass_Height(Nx_Max,Ny_Max),      &
     &      Pass_Cg(Nx_Max,Ny_Max),                                     &
     &      Pass_C(Nx_Max,Ny_Max),pass_area(nx_max,ny_max),             &
     &      Pass_period,                                                &
     &      Intp_U_Wave(Nx_Max,Ny_Max), Intp_V_Wave(Nx_Max,Ny_Max),     &
     &      Intp_eta_Wave(Nx_Max,Ny_Max),Pass_uw(nx_max,ny_max,100)
 
       real Pass_ibrk(Nx_Max,Ny_Max)

!c -- circulation module:

       real Pass_U(Nx_Max,Ny_Max),Pass_V(Nx_Max,Ny_Max),                &
     &      Pass_Ub(Nx_Max,Ny_Max),Pass_Vb(Nx_Max,Ny_Max),              &
     &      Pass_eta(Nx_Max,Ny_Max),                                    &
     &      Pass_d11(Nx_Max,Ny_Max), Pass_d12(Nx_Max,Ny_Max),           &
     &      Pass_e11(Nx_Max,Ny_Max), Pass_e12(Nx_Max,Ny_Max),           &
     &      Pass_f11(Nx_Max,Ny_Max), Pass_f12(Nx_Max,Ny_Max),           &
     &      Pass_fw(Nx_Max,Ny_Max),Pass_vt(Nx_Max,Ny_Max),              &
     &      Intp_Fx_Circ(Nx_Max,Ny_Max),Intp_Fy_Circ(Nx_Max,Ny_Max),    &
     &      Intp_ubott_Circ(Nx_Max,Ny_Max),                             &
     &      Intp_Theta_Circ(Nx_Max,Ny_Max),                             &
     &      Intp_Height_Circ(Nx_Max,Ny_Max),                            &
     &      Intp_C_Circ(Nx_Max,Ny_Max),                                 &
     &      Intp_Cg_Circ(Nx_Max,Ny_Max),                                &
     &      Intp_WaveNum_Circ(Nx_Max,Ny_Max),                           &
     &      Intp_Diss_Circ(Nx_Max,Ny_Max),                              &
     &      Intp_ibrk_Circ(Nx_Max,Ny_Max),                              &
     &      Intp_Sxx_Circ(Nx_Max,Ny_Max),Intp_Sxy_Circ(Nx_Max,Ny_Max),  &
     &      Intp_Syy_Circ(Nx_Max,Ny_Max),                               &
     &      Intp_Sxx_Surf(Nx_Max,Ny_Max),Intp_Sxy_Surf(Nx_Max,Ny_Max),  &
     &      Intp_Syy_Surf(Nx_Max,Ny_Max),                               &
     &      Intp_Sxx_Body(Nx_Max,Ny_Max),Intp_Sxy_Body(Nx_Max,Ny_Max),  &
     &      Intp_Syy_Body(Nx_Max,Ny_Max),                               &
     &      Intp_MassFluxU_Circ(Nx_Max,Ny_Max),                         &
     &      Intp_MassFluxV_Circ(Nx_Max,Ny_Max)


!c -- sediment module:

      real Pass_Dupdated(Nx_Max,Ny_Max),                                &
     &      Pass_sedflux_x(Nx_Max,Ny_Max),                              &
     &      Pass_sedflux_y(Nx_Max,Ny_Max),                              &
     &      Pass_sedfluxcum_x(Nx_Max,Ny_Max),                           &
     &      Pass_sedfluxcum_y(Nx_Max,Ny_Max),                           &
     &      Intp_U_Sedi(Nx_Max,Ny_Max),                                 &
     &      Intp_V_Sedi(Nx_Max,Ny_Max),                                 &
     &      Intp_Ub_Sedi(Nx_Max,Ny_Max),                                &
     &      Intp_Vb_Sedi(Nx_Max,Ny_Max),                                &
     &      Intp_ubott_Sedi(Nx_Max,Ny_Max),                             &
     &      Intp_eta_Sedi(Nx_Max,Ny_Max),                               &
     &      Intp_fw_Sedi(Nx_Max,Ny_Max),                                &
     &      Intp_vt_Sedi(Nx_Max,Ny_Max),                                &
     &      Intp_Theta_Sedi(Nx_Max,Ny_Max),                             &
     &      Intp_Height_Sedi(Nx_Max,Ny_Max),                            &
     &      Intp_ibrk_Sedi(Nx_Max,Ny_Max)

!c -- coordinate systems, depth and ...

       real Depth_Circ(Nx_Max,Ny_Max),Depth_Wave(Nx_Max,Ny_Max),        &
     &      X_Wave(Nx_Max,Ny_Max),Y_Wave(Nx_Max,Ny_Max),                &
     &      X_Circ(Nx_Max,Ny_Max),                                      &
     &      Y_Circ(Nx_Max,Ny_Max), Pass_tide,                           &
     &      U_wind_Mast(Nx_Max,Ny_Max),V_wind_Mast(Nx_Max,Ny_Max),      &
     &      U_wind_Circ(Nx_Max,Ny_Max),V_wind_CIrc(Nx_Max,Ny_Max),      &
     &      U_wind_Wave(Nx_Max,Ny_Max),V_wind_Wave(Nx_Max,Ny_Max),      &
     &      Depth_Sedi(Nx_max,Ny_Max),X_Sedi(Nx_Max,Ny_Max),            &
     &      Y_Sedi(Nx_Max,Ny_Max),                                      &
     &      X_Mast(Nx_Max,Ny_Max),Y_Mast(Nx_Max,Ny_Max),                &
     &      Depth_Mast(Nx_Max,Ny_Max)

!c -- vector rotate
       real Circ_Rotate_Angle, Wave_Rotate_Angle,                       &
     &      Sedi_Rotate_Angle       
 
!c -- control parameters

       integer Master_Start,nWave,nCirc,nSedi,nOut

       integer Wave_Stag_huv(3), Circ_Stag_huv(3),Sedi_Stag_huv(3)

       real N_Interval_CallWave,N_Interval_CallCirc,N_Interval_CallSedi,&
     &     N_Delay_CallSedi,N_Interval_Output   
       real Total_Time, Master_dt,Time_Master
       logical Grid_Mast_Wave_Same, Grid_Mast_Circ_Same,                &
     &         Grid_Mast_Sedi_Same, Grid_Wave_Circ_Same,                &
     &         Grid_Wave_Sedi_Same, Grid_Circ_Sedi_Same,                &
     &         Wave_Staggered, Circ_Staggered,Sedi_Staggered,           &
     &         Wave_Structured, Circ_Structured, Sedi_Structured,       &
     &         Grid_Extrapolation,                                      &
     &         Wave_Curr_Interact,                                      &
     &         Wave_Bed_Interact,                                       &
     &         Curr_Bed_Interact,                                       &
     &         Wave_To_Circ_Height,                                     &
     &         Wave_To_Circ_Angle,                                      &
     &         Wave_To_Circ_WaveNum,                                    &
     &         Wave_To_Circ_C,                                          &
     &         Wave_To_Circ_Cg,                                         &
     &         Wave_To_Circ_Radiation,                                  &
     &         Wave_To_Circ_Rad_Surf,                                   &
     &         Wave_To_Circ_Rad_Body,                                   &
     &         Wave_To_Circ_Forcing,                                    &
     &         Wave_To_Circ_MassFlux,                                   &
     &         Wave_To_Circ_Dissipation,                                &
     &         Wave_To_Circ_BottomUV,                                   &
     &         Wave_To_Circ_Brkindex,                                   &
     &         Circ_To_Wave_UV,                                         &
     &         Circ_To_Wave_eta,                                        &
     &         Wave_To_Sedi_Height,                                     &
     &         Wave_To_Sedi_Angle,                                      &
     &         Wave_To_Sedi_BottomUV,                                   &
     &         Circ_To_Sedi_UV,                                         &
     &         Circ_To_Sedi_UVb,                                        &
     &         Circ_To_Sedi_eta,                                        &
     &         Circ_To_Sedi_UV3D,                                       &
     &         Circ_To_Sedi_fw,                                         &
     &         Circ_To_Sedi_vt,                                         &
     &         Circ_To_Sedi_UVquasi3D,                                  &
     &         Sedi_To_Wave_Depth,                                      &
     &         Sedi_To_Circ_Depth,                                      &
     &         Circ_POM,                                                &
     &         Circ_SC,                                                 &
     &         Waveupdat_for_Circ,                                      &
     &         Waveupdat_for_Sedi,                                      &
     &         Circupdat_for_Wave,                                      &
     &         Circupdat_for_Sedi,                                      &
     &         Sediupdat_for_Wave,                                      &
     &         Sediupdat_for_Circ

!c -- file names

      character*255 f_depth,f_xymast,f_xywave,f_xycirc,f_xysedi,        &
     &              f_name6,f_name7,                                    &
     &              f_name8,                                            &
     &              f_name9,f_name10,f_name11,f_name12,                 &
     &              f_name13,f_name14,                                  &
     &              f_name15,f_name16

!c -- common blocks

       common/mast_waves/Pass_Sxx,Pass_Sxy,Pass_Syy,                    &
     &             Pass_Sxx_body,Pass_Sxy_body,Pass_Syy_body,           &
     &             Pass_Sxx_surf,Pass_Sxy_surf,Pass_Syy_surf,           &
     &             Pass_Wave_Fx,Pass_Wave_Fy,                           &
     &             Pass_MassFluxU,                                      &
     &             Pass_MassFluxV,Pass_MassFlux,                        &
     &             Pass_Diss, Pass_WaveNum, Pass_Theta,                 &
     &             Pass_ubott,Pass_Height, Pass_tide,                   &
     &             Depth_Wave,                                          &
     &             X_Wave,Y_Wave,                                       &
     &             Pass_C,Pass_Cg,pass_area,                            &
     &             Pass_ibrk,Pass_period,                               &
     &             Intp_U_Wave, Intp_V_Wave,                            &
     &             Intp_eta_Wave,pass_uw,                               &
     &             Nx_Wave,Ny_Wave
        common/mast_circ/Depth_Circ,Pass_U,Pass_V,                      &
     &              Pass_Ub,Pass_Vb,                                    &
     &              X_Circ,Y_Circ,                                      &
     &              Pass_eta,                                           &
     &              Pass_d11,Pass_d12,Pass_e11,Pass_e12,                &
     &              Pass_f11,Pass_f12,                                  &
     &              Pass_fw,pass_vt,                                    &
     &              Intp_Fx_Circ,Intp_Fy_Circ,                          &
     &              Intp_ubott_Circ,                                    &
     &              Intp_Theta_Circ,                                    &
     &              Intp_Height_Circ,                                   &
     &              Intp_C_Circ,                                        &
     &              Intp_Cg_Circ,                                       &
     &              Intp_WaveNum_Circ,                                  &
     &              Intp_ibrk_Circ,                                     &
     &              Intp_Sxx_Circ,Intp_Sxy_Circ,Intp_Syy_Circ,          &
     &              Intp_Sxx_Surf,Intp_Sxy_Surf,Intp_Syy_Surf,          &
     &              Intp_Sxx_Body,Intp_Sxy_Body,Intp_Syy_Body,          &
     &              Intp_MassFluxU_Circ,                                &
     &              Intp_MassFluxV_Circ,                                &
     &              Intp_Diss_Circ,                                     &
     &              Nx_Circ,Ny_Circ

        common/mast_sedi/Pass_Dupdated,                                 &
     &              pass_sedflux_x,pass_sedflux_y,                      &
     &              pass_sedfluxcum_x,pass_sedfluxcum_y,                &
     &              Depth_Sedi,X_Sedi,Y_Sedi,                           &
     &              Intp_U_Sedi,Intp_V_Sedi,                            &
     &              Intp_Ub_Sedi,Intp_Vb_Sedi,                          &
     &              Intp_ubott_Sedi,                                    &
     &              Intp_eta_Sedi,                                      &
     &              Intp_fw_Sedi,                                       &
     &              Intp_vt_Sedi,                                       &
     &              Intp_Theta_Sedi,                                    &
     &              Intp_ibrk_Sedi,                                     &
     &              Intp_Height_Sedi,                                   &
     &              Nx_Sedi,Ny_Sedi

      common/mast_wind/U_Wind_Mast, V_Wind_Mast,                        &
     &              U_Wind_Circ, V_Wind_Circ,                           &
     &              U_Wind_Wave, V_Wind_Wave

      common/mast_rotate/Circ_Rotate_Angle, Wave_Rotate_Angle,          &
     &                Sedi_Rotate_Angle 

      common/mast_contr/X_Mast,Y_Mast,                                  &
     &              Depth_Mast,                                         &
     &              Total_Time, Master_dt, Time_Master,                 &
     &              N_Interval_CallCirc, N_Interval_CallWave,           &
     &              N_Interval_CallSedi, N_Delay_CallSedi,              &
     &              N_Interval_Output,                                  &
     &              Master_Start,                                       &
     &              Nx_Mast, Ny_Mast,                                   &
     &              nWave,nCirc,nSedi,nOut,                             &
     &              Wave_Stag_huv,Circ_Stag_huv,Sedi_Stag_huv,          &
     &              Grid_Mast_Wave_Same, Grid_Mast_Circ_Same,           &
     &              Grid_Mast_Sedi_Same,                                &
     &              Grid_Wave_Circ_Same,                                &
     &              Grid_Wave_Sedi_Same, Grid_Circ_Sedi_Same,           &
     &              Wave_Staggered, Circ_Staggered,Sedi_Staggered,      &
     &              Wave_structured, Circ_Structured, Sedi_Structured,  &
     &              Grid_Extrapolation,                                 &
     &              Wave_Curr_Interact,                                 &
     &              Wave_Bed_Interact,                                  &
     &              Curr_Bed_Interact,                                  &
     &              Wave_To_Circ_Height,                                &
     &              Wave_To_Circ_Angle,                                 &
     &              Wave_To_Circ_WaveNum,                               &
     &              Wave_To_Circ_C,                                     &
     &              Wave_To_Circ_Cg,                                    &
     &              Wave_To_Circ_Radiation,                             &
     &              Wave_To_Circ_Rad_Surf,                              &
     &              Wave_To_Circ_Rad_Body,                              &
     &              Wave_To_Circ_Forcing,                               &
     &              Wave_To_Circ_MassFlux,                              &
     &              Wave_To_Circ_Dissipation,                           &
     &              Wave_To_Circ_BottomUV,                              &
     &              Wave_To_Circ_Brkindex,                              &
     &              Circ_To_Wave_UV,                                    &
     &              Circ_To_Wave_eta,                                   &
     &              Wave_To_Sedi_Height,                                &
     &              Wave_To_Sedi_Angle,                                 &
     &              Wave_To_Sedi_BottomUV,                              &
     &              Circ_To_Sedi_UV,                                    &
     &              Circ_To_Sedi_UVb,                                   &
     &              Circ_To_Sedi_eta,                                   &
     &              Circ_To_Sedi_UV3D,                                  &
     &              Circ_To_Sedi_UVquasi3D,                             &
     &              Sedi_To_Wave_Depth,                                 &
     &              Sedi_To_Circ_Depth,                                 &
     &              Circ_POM,                                           &
     &              Circ_SC,                                            &
     &              Waveupdat_for_Circ,                                 &
     &              Waveupdat_for_Sedi,                                 &
     &              Circupdat_for_Wave,                                 &
     &              Circupdat_for_Sedi,                                 &
     &              Sediupdat_for_Wave,                                 &
     &              Sediupdat_for_Circ


      common/mast_filename/f_depth,f_xymast,f_xywave,f_xycirc,          &
     &        f_xysedi,f_name6,f_name7,                                 &
     &        f_name8,                                                  &
     &        f_name9,f_name10,f_name11,f_name12,f_name13,f_name14,     &
     &        f_name15,f_name16 
