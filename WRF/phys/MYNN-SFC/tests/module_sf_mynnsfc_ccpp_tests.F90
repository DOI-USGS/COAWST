! Module for MYNN SFC scheme tests
module module_sf_mynnsfc_ccpp_tests
  implicit none

  contains
    !=================================================================================================================    

    subroutine init_mynn_sfc_flags_for_test_all_true()
      write(*,*) '--- calling  init_mynn_sfc_flags_for_test_all_true()'
      ! for future use
    end subroutine init_mynn_sfc_flags_for_test_all_true
    !=================================================================================================================    
    subroutine init_input_data_for_test(saveoutput)
      logical, intent(in) :: saveoutput
      integer :: iostat, line_num
      character(len=2000) :: input_line
      integer, parameter :: input_unit = 10
      integer, parameter :: output_unit = 20
      
      write(*,*) '--- opening data files ---'
      ! Open input file
      open(unit=input_unit, file='./data/ccpp_input_lnd.txt', status='old', action='read', iostat=iostat)

      if (iostat /= 0) then
          print *, 'Error opening input file'
          stop
      end if

      ! Open output file
      if (saveoutput) then
          open(unit=output_unit, file='./data/ccpp_output_lnd.txt', status='replace', action='write', iostat=iostat)
          write(output_unit,'(A5, A5, A10, A10, A10, A10, A10, A10, A10, A10, A10)')                &
                'itimestep', 'iter', 'T2', 'Q2', 'TH2', 'U10', 'V10', 'HFX', 'LH', 'UST_lnd','PBLH'
          if (iostat /= 0) then
              print *, 'Error opening output file'
              close(output_unit)
              stop
          end if
      end if

    end subroutine init_input_data_for_test

    !===============================================================================
    ! Subroutine to process each line
    !===============================================================================
    subroutine process_line(line, out_unit, line_number, flag_iter,  &                         
           U1D,V1D,T1D,QV1D,P1D,dz8w1d,                              &
           U1D2,V1D2,dz2w1d,                                         &
           PSFCPA,PBLH,MAVAIL,XLAND,DX,                              &
           ISFFLX,isftcflx,iz0tlnd,psi_opt,                          &
           compute_flux,compute_diag,                                &
           sigmaf,vegtype,shdmax,ivegsrc,                            & 
           z0pert,ztpert,                                            &
           redrag,sfc_z0_type,                                       &
           itimestep,iter,flag_restart,lsm,lsm_ruc,                  &
           wet,          dry,          icy,                          &
           tskin_wat,    tskin_lnd,    tskin_ice,                    &
           tsurf_wat,    tsurf_lnd,    tsurf_ice,                    &
           qsfc_wat,     qsfc_lnd,     qsfc_ice,                     &
           snowh_wat,    snowh_lnd,    snowh_ice,                    &
           ZNT_wat,      ZNT_lnd,      ZNT_ice,                      &
           UST_wat,      UST_lnd,      UST_ice,                      &
           cm_wat,       cm_lnd,       cm_ice,                       &
           ch_wat,       ch_lnd,       ch_ice,                       &
           rb_wat,       rb_lnd,       rb_ice,                       &
           stress_wat,   stress_lnd,   stress_ice,                   &
           psix_wat,     psix_lnd,     psix_ice,                     &
           psit_wat,     psit_lnd,     psit_ice,                     &
           psix10_wat,   psix10_lnd,   psix10_ice,                   &
           psit2_wat,    psit2_lnd,    psit2_ice,                    &
           HFLX_wat,     HFLX_lnd,     HFLX_ice,                     &
           QFLX_wat,     QFLX_lnd,     QFLX_ice,                     &
           ch,CHS,CHS2,CQS2,CPM,                                     &
           ZNT,USTM,ZOL,MOL,RMOL,                                    &
           PSIM,PSIH,                                                &
           HFLX,HFX,QFLX,QFX,LH,FLHC,FLQC,                           &
           QGH,QSFC,                                                 &
           U10,V10,TH2,T2,Q2,                                        &
           GZ1OZ0,WSPD,wstar,qstar,                                  &
           spp_sfc,  rstoch1D )         
        
        implicit none
        
        character(len=2000), intent(in) :: line
        integer, intent(in) :: out_unit
        integer, intent(in) :: line_number
        
        ! Variables to store parsed values
        logical, intent(out) :: flag_iter, compute_flux, compute_diag, redrag, flag_restart
        logical, intent(out) :: wet, dry, icy
        real, intent(out) :: U1D, V1D, T1D, QV1D, P1D, dz8w1d, U1D2, V1D2, dz2w1d, &
                PSFCPA, PBLH, MAVAIL, XLAND, DX, sigmaf, shdmax, z0pert,           &
                ztpert
        real, intent(out) :: tskin_wat, tskin_lnd, tskin_ice, tsurf_wat, tsurf_lnd, tsurf_ice, &
                qsfc_wat, qsfc_lnd, qsfc_ice, snowh_wat, snowh_lnd, snowh_ice,                 &
                ZNT_wat, ZNT_lnd, ZNT_ice, UST_wat, UST_lnd, UST_ice,                          &
                cm_wat, cm_lnd, cm_ice, ch_wat, ch_lnd, ch_ice,                                &
                rb_wat, rb_lnd, rb_ice, stress_wat, stress_lnd, stress_ice,                    &                   
                psix_wat,     psix_lnd,     psix_ice,                                          &
                psit_wat,     psit_lnd,     psit_ice,                                          &
                psix10_wat,   psix10_lnd,   psix10_ice,                                        &
                psit2_wat,    psit2_lnd,    psit2_ice,                                         &
                HFLX_wat,     HFLX_lnd,     HFLX_ice,                                          &
                QFLX_wat,     QFLX_lnd,     QFLX_ice,                                          &
                ch,CHS,CHS2,CQS2,CPM,                                                          &
                ZNT,USTM,ZOL,MOL,RMOL,                                                         &
                PSIM,PSIH,                                                                     &
                HFLX,HFX,QFLX,QFX,LH,FLHC,FLQC,                                                &
                QGH,QSFC,                                                                      &
                U10,V10,TH2,T2,Q2,                                                             &
                GZ1OZ0,WSPD,wstar,qstar,                                                       &
                rstoch1D  
        integer :: read_stat           
        integer, intent(out) :: ISFFLX, isftcflx, iz0tlnd, psi_opt, vegtype,                   &
                   ivegsrc, sfc_z0_type, itimestep, iter, lsm, lsm_ruc, spp_sfc
        
        ! Read values from the line with specified format
        read (line, '(L5,E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,'          // &
                     'E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,' // &
                     'I5,I5,I5,I5,'                                     // &
                     'L5,L5,'                                           // &
                     'E15.3,I5,E15.3,I5,'                               // &
                     'E15.3,E15.3,'                                     // &
                     'L5,I5,'                                           // &
                     'I5,I5,L5,I5,I5,'                                  // &
                     'L5,L5,L5,'                                        // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // & !stress_* values
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,E15.3,E15.3,'                   // &
                     'E15.3,E15.3,E15.3,E15.3,E15.3,'                   // &
                     'E15.3,E15.3,'                                     // &
                     'E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,'       // &
                     'E15.3,E15.3,'                                     // &
                     'E15.3,E15.3,E15.3,E15.3,E15.3,'                   // &
                     'E15.3,E15.3,E15.3,E15.3,'                         // &
                     'I5,E15.3)', iostat=read_stat)                        &
                      flag_iter,                                           &
                      U1D,   V1D,   T1D,   QV1D,   P1D,   dz8w1d,          &
                      U1D2,   V1D2,   dz2w1d,                              &
                      PSFCPA,   PBLH,   MAVAIL,   XLAND,   DX,             &
                      ISFFLX,isftcflx,iz0tlnd,psi_opt,                     &
                      compute_flux,compute_diag,                           &
                      sigmaf,   vegtype,   shdmax,   ivegsrc,              &
                      z0pert,   ztpert,                                    &
                      redrag,sfc_z0_type,                                  &
                      itimestep,iter,flag_restart,lsm,lsm_ruc,             &
                            wet,             dry,             icy,         &
                      tskin_wat,       tskin_lnd,       tskin_ice,         &
                      tsurf_wat,       tsurf_lnd,       tsurf_ice,         &
                       qsfc_wat,        qsfc_lnd,        qsfc_ice,         &
                      snowh_wat,       snowh_lnd,       snowh_ice,         &
                        ZNT_wat,         ZNT_lnd,         ZNT_ice,         &
                        UST_wat,         UST_lnd,         UST_ice,         &
                         cm_wat,          cm_lnd,          cm_ice,         &
                         ch_wat,          ch_lnd,          ch_ice,         &
                         rb_wat,          rb_lnd,          rb_ice,         &
                     stress_wat,      stress_lnd,      stress_ice,         &
                     psix_wat,     psix_lnd,     psix_ice,                 & 
                     psit_wat,     psit_lnd,     psit_ice,                 & 
                     psix10_wat,   psix10_lnd,   psix10_ice,               &
                     psit2_wat,    psit2_lnd,    psit2_ice,                &
                     HFLX_wat,     HFLX_lnd,     HFLX_ice,                 &
                     QFLX_wat,     QFLX_lnd,     QFLX_ice,                 &
                     ch,CHS,CHS2,CQS2,CPM,                                 &
                     ZNT,USTM,ZOL,MOL,RMOL,                                &
                     PSIM,PSIH,                                            &
                     HFLX,HFX,QFLX,QFX,LH,FLHC,FLQC,                       &
                     QGH,QSFC,                                             &
                     U10,V10,TH2,T2,Q2,                                    &
                     GZ1OZ0,WSPD,wstar,qstar,                              &
                     spp_sfc,rstoch1D                                              
              
      ! Debug: print the line being read
      write(0,*) "Line length:", len_trim(line)
      write(0,*) "U1D=",U1D !, "flag_iter=",flag_iter, "V1D=", V1D," ZNT", ZNT_lnd, "WSPD=",WSPD
        
    end subroutine process_line

    !=================================================================================================================
    subroutine ccpp_test(saveoutput)

      !use module_sf_mynnsfc_driver, only : SFCLAY1D_mynn
      use module_sf_mynnsfc, only : SFCLAY1D_mynn

      logical, intent(in) :: saveoutput
      integer :: iostat, line_num
      integer, parameter :: ids=1, ide=1, jds=1, jde=1, kds=1,kde=1
      integer, parameter :: ims=1, ime=1, jms=1, jme=1, kms=1,kme=1
      integer, parameter :: its=1, ite=1, jts=1, jte=1, kts=1,kte=1

      character(len=2000) :: input_line
      integer, parameter :: input_unit = 10
      integer, parameter :: output_unit = 20
      integer, parameter :: n = 1  ! Number of points

      logical, dimension(n)  :: flag_iter 
      logical :: compute_flux, compute_diag, redrag, flag_restart
      logical, dimension(n)  :: wet, dry, icy
      real, dimension(n) ::  U1D, V1D, T1D, QV1D, P1D, dz8w1d, U1D2, V1D2, dz2w1d, &
              PSFCPA, PBLH, MAVAIL, XLAND, DX, sigmaf, shdmax, z0pert,             &
              ztpert
      integer :: J,spp_sfc
      real, dimension(n) :: tskin_wat, tskin_lnd, tskin_ice, tsurf_wat, tsurf_lnd, tsurf_ice, &
              qsfc_wat, qsfc_lnd, qsfc_ice, snowh_wat, snowh_lnd, snowh_ice,                  &
              ZNT_wat, ZNT_lnd, ZNT_ice, UST_wat, UST_lnd, UST_ice,                           &
              cm_wat, cm_lnd, cm_ice, ch_wat, ch_lnd, ch_ice,                                 &
              rb_wat, rb_lnd, rb_ice, stress_wat, stress_lnd, stress_ice,                     &
                   psix_wat,     psix_lnd,     psix_ice,                                      &
                   psit_wat,     psit_lnd,     psit_ice,                                      &
                 psix10_wat,   psix10_lnd,   psix10_ice,                                      &
                  psit2_wat,    psit2_lnd,    psit2_ice,                                      &
                   HFLX_wat,     HFLX_lnd,     HFLX_ice,                                      &
                   QFLX_wat,     QFLX_lnd,     QFLX_ice,                                      &
                 ch,CHS,CHS2,CQS2,CPM,                                                        &
                 ZNT,USTM,ZOL,MOL,RMOL,                                                       &
                 PSIM,PSIH,                                                                   &
                 HFLX,HFX,QFLX,QFX,LH,FLHC,FLQC,                                              &
                 QGH,QSFC,                                                                    &
                 U10,V10,TH2,T2,Q2,                                                           &
                 GZ1OZ0,WSPD,wstar,qstar,                                                     &
                 rstoch1D     
      integer, dimension(n) :: vegtype
                  
      integer :: read_stat,ISFFLX,isftcflx,iz0tlnd,psi_opt,ivegsrc,sfc_z0_type,itimestep, iter,lsm, lsm_ruc
      
      ! Variables for error handling
      character(len=512) :: errmsg
      integer :: errflg
      ! Initialize error variables
      errmsg = ''
      errflg = 0

      write(*,*) '--- entering ccpp_test subroutine'    
      ! Initialize input data for tests
      call init_input_data_for_test(saveoutput)

      ! Read header
      read(input_unit, '(A)', iostat=iostat) input_line

      ! Process each line
      line_num = 0
      do
          read(input_unit, '(A)', iostat=iostat) input_line
          write(0,*) input_line
          
          ! Check for end of file or error
          if (iostat < 0) exit  ! End of file
          if (iostat > 0) then
              print *, 'Error reading line', line_num + 1
              exit
          end if
          
          line_num = line_num + 1
          
          ! Call subroutine to process the line
          call process_line(input_line, output_unit, line_num, flag_iter(1),   &
           U1D(1),  V1D(1),  T1D(1),  QV1D(1),  P1D(1),  dz8w1d(1),            &
           U1D2(1),  V1D2(1),  dz2w1d(1),                                      &
           PSFCPA(1),  PBLH(1),  MAVAIL(1),  XLAND(1),  DX(1),                 &
           ISFFLX,  isftcflx,  iz0tlnd,  psi_opt,                              &
           compute_flux,  compute_diag,                                        &
           sigmaf(1),  vegtype(1),  shdmax(1),  ivegsrc,                       &
           z0pert(1),  ztpert(1),                                              &
           redrag,  sfc_z0_type,                                               &
           itimestep,  iter,  flag_restart,  lsm,  lsm_ruc,                    &
                  wet(1),            dry(1),            icy(1),                &
            tskin_wat(1),      tskin_lnd(1),      tskin_ice(1),                &
            tsurf_wat(1),      tsurf_lnd(1),      tsurf_ice(1),                &
             qsfc_wat(1),       qsfc_lnd(1),       qsfc_ice(1),                &
            snowh_wat(1),      snowh_lnd(1),      snowh_ice(1),                &
              ZNT_wat(1),        ZNT_lnd(1),        ZNT_ice(1),                &
              UST_wat(1),        UST_lnd(1),        UST_ice(1),                &
               cm_wat(1),         cm_lnd(1),         cm_ice(1),                &
               ch_wat(1),         ch_lnd(1),         ch_ice(1),                &
               rb_wat(1),         rb_lnd(1),         rb_ice(1),                &
           stress_wat(1),     stress_lnd(1),     stress_ice(1),                &
             psix_wat(1),       psix_lnd(1),       psix_ice(1),                &
             psit_wat(1),       psit_lnd(1),       psit_ice(1),                &
             psix10_wat(1),     psix10_lnd(1),     psix10_ice(1),              &
             psit2_wat(1),      psit2_lnd(1),      psit2_ice(1),               &
             HFLX_wat(1),       HFLX_lnd(1),       HFLX_ice(1),                &
             QFLX_wat(1),       QFLX_lnd(1),       QFLX_ice(1),                &
             ch(1),  CHS(1),  CHS2(1),  CQS2(1),  CPM(1),                      &
             ZNT(1),  USTM(1),  ZOL(1),  MOL(1),  RMOL(1),                     &
             PSIM(1),  PSIH(1),                                                &
             HFLX(1),  HFX(1),  QFLX(1),  QFX(1),  LH(1),  FLHC(1),  FLQC(1),  &
             QGH(1),  QSFC(1),                                                 &
             U10(1),  V10(1),  TH2(1),  T2(1),  Q2(1),                         &
             GZ1OZ0(1),  WSPD(1),  wstar(1),  qstar(1),                        &
             spp_sfc,  rstoch1D(1))

         ! Initialize MYNN SFC
          write(*,*) '--- calling  SFCLAY1D_mynn()'

          call SFCLAY1D_mynn(flag_iter,                                                                             &
             J=99999,U1D=U1D,V1D=V1D,T1D=T1D,QV1D=QV1D,P1D=P1D,dz8w1d=dz8w1d,U1D2=U1D2,V1D2=V1D2,dz2w1d=dz2w1d,     &
             PSFCPA=PSFCPA,PBLH=PBLH,MAVAIL=MAVAIL,XLAND=XLAND,DX=DX,                                               &
             ISFFLX=ISFFLX,isftcflx=isftcflx,iz0tlnd=iz0tlnd,psi_opt=psi_opt,                                       &
             compute_flux=compute_flux,compute_diag=compute_diag,                                                   &
             sigmaf=sigmaf,vegtype=vegtype,shdmax=shdmax,ivegsrc=ivegsrc,                                           &  !intent(in)
             z0pert=z0pert,ztpert=ztpert,                                                                           &  !intent(in)
             redrag=redrag,sfc_z0_type=sfc_z0_type,                                                                 &  !intent(in)
             itimestep=itimestep,iter=iter,flag_restart=flag_restart,lsm=lsm,lsm_ruc=lsm_ruc,                       &
                    wet=wet,          dry=dry,          icy=icy,                                                    &  !intent(in)
              tskin_wat=tskin_wat,    tskin_lnd=tskin_lnd,    tskin_ice=tskin_ice,                                  &  !intent(in)
              tsurf_wat=tsurf_wat,    tsurf_lnd=tsurf_lnd,    tsurf_ice=tskin_ice,                                  &  !intent(in)
               qsfc_wat=qsfc_wat,     qsfc_lnd=qsfc_lnd,     qsfc_ice=qsfc_ice,                                     &  !intent(in)
              snowh_wat=snowh_wat,    snowh_lnd=snowh_lnd,    snowh_ice=snowh_ice,                                  &  !intent(in)
                ZNT_wat=ZNT_wat,      ZNT_lnd=ZNT_lnd,      ZNT_ice=ZNT_ice,                                        &  !intent(inout)
                UST_wat=UST_wat,      UST_lnd=UST_lnd,      UST_ice=UST_ice,                                        &  !intent(inout)
                 cm_wat=cm_wat,       cm_lnd=cm_lnd,       cm_ice=cm_ice,                                           &  !intent(inout)
                 ch_wat=ch_wat,       ch_lnd=ch_lnd,       ch_ice=ch_ice,                                           &  !intent(inout)
                 rb_wat=rb_wat,       rb_lnd=rb_lnd,       rb_ice=rb_ice,                                           &  !intent(inout)
             stress_wat=stress_wat,   stress_lnd=stress_lnd,   stress_ice=stress_ice,                               &  !intent(inout)
                psix_wat=psix_wat,     psix_lnd=psix_lnd,     psix_ice=psix_ice,                                    &  !=fm, intent(inout)
               psit_wat=psit_wat,     psit_lnd=psit_lnd,     psit_ice=psix_ice,                                     &  !=fh, intent(inout)
             psix10_wat=psix10_wat,   psix10_lnd=psix10_lnd,   psix10_ice=psix10_ice,                               &  !=fm10, intent(inout)
              psit2_wat=psit2_wat,    psit2_lnd=psit2_lnd,    psit2_ice=psit2_ice,                                  &  !=fh2, intent(inout)
               HFLX_wat=HFLX_wat,     HFLX_lnd=HFLX_lnd,     HFLX_ice=HFLX_ice,                                     &
               QFLX_wat=QFLX_wat,     QFLX_lnd=QFLX_lnd,     QFLX_ice=QFLX_ice,                                     &
             ch=ch,CHS=CHS,CHS2=CHS2,CQS2=CQS2,CPM=CPM,                                                             &
             ZNT=ZNT,USTM=USTM,ZOL=ZOL,MOL=MOL,RMOL=RMOL,                                                           &
             PSIM=PSIM,PSIH=PSIH,                                                                                   &
             HFLX=HFLX,HFX=HFX,QFLX=QFLX,QFX=QFX,LH=LH,FLHC=FLHC,FLQC=FLQC,                                         &
             QGH=QGH,QSFC=QSFC,                                                                                     &
             U10=U10,V10=V10,TH2=TH2,T2=T2,Q2=Q2,                                                                   &
             GZ1OZ0=GZ1OZ0,WSPD=WSPD,wstar=wstar,qstar=qstar,                                                       &
             spp_sfc=spp_sfc,rstoch1D=rstoch1D,                                                                     &
             ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde,                                                     &
             ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme,                                                     &
             its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte,                                                     &
             errmsg= errmsg, errflg=errflg                                     )
         
         write(0,*) "T2=",T2
         write(0,*) "Read status:", read_stat
         if (saveoutput) then
             open(output_unit, file = './data/ccpp_output_lnd.txt')
             write(output_unit,'(I5, I5, F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F15.2,F10.2)') itimestep, &
                  iter,T2,Q2,TH2,U10,V10,HFX,LH,UST_lnd,PBLH
         end if

      end do
    end subroutine ccpp_test
    !=================================================================================================================

end module module_sf_mynnsfc_ccpp_tests           
