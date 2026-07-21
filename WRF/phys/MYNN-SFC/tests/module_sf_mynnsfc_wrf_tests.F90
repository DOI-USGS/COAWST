! Module for MYNN SFC scheme tests
module module_sf_mynnsfc_wrf_tests
    use module_sf_mynnsfc_driver
    ! public
    !=================================================================================================================    
    implicit none
    logical :: cycling,restart,flag_iter
    integer :: initflag, spp_pbl, isfflx, flag_lsm,lsm
    real,dimension(1) :: pattern_spp_pbl
    contains

    subroutine init_mynn_sfc_flags()
      write(*,*) '--- calling  init_mynn_sfc_flags() ---'
       cycling=.false.
       initflag=1
       spp_pbl=1
       restart=.false.
       pattern_spp_pbl=0.0
       isfflx =1
       flag_lsm=3
       lsm=3
       flag_iter = .true.

    end subroutine init_mynn_sfc_flags
    !=================================================================================================================    
    subroutine init_input_data_for_test(case,saveoutput)
      integer :: iostat, line_num
      character(len=*),intent(in) :: case
      logical,intent(in) :: saveoutput
      character(len=2000) :: input_line
      integer, parameter :: input_unit = 10
      integer, parameter :: output_unit = 20
      
      write(*,*) '--- opening data files ---'
      ! Open input file
      close(unit=input_unit)
      open(unit=input_unit, file='./data/input_'//trim(case)//'.txt', status='old', action='read', iostat=iostat)

      if (iostat /= 0) then
          print *, 'Error opening input file'
          stop
      end if

      ! Open output file
      if (saveoutput) then 
          close(unit=output_unit)
          open(unit=output_unit, file='./data/output_'//trim(case)//'.txt', status='replace', action='write', iostat=iostat)
          write(output_unit,'(A5, A10, A10, A10, A10, A10, A10, A10, A10, A10)')                &
                'itimestep',  'T2', 'Q2', 'TH2', 'U10', 'V10', 'HFX', 'LH', 'UST','PBLH'
          if (iostat /= 0) then
              print *, 'Error opening output file'
              close(output_unit)
              stop
          end if    
      end if

    end subroutine init_input_data_for_test


    !===============================================================================
    ! Subroutine to process each line from wrf v4.5.1 processed input
    !===============================================================================
    subroutine process_line_wrf(line, out_unit, line_number,    &                         
            xlat,xlong,itimestep,                               &
            dx, xland,                                          & 
            U1D,V1D,T1D,QV1D,                                   & 
            P1D, dz8w1d,RHO1D, U1D2,                            &
            V1D2, dz2w1d,                                       &
            pblh, znt,PSFCPA,MAVAIL,                            &
            TSK,snowh,                                          &
            chs,chs2,cqs2,                                      &
            rmol,zol,mol,                                       &
            hfx,qfx,                                            &
            ust,ustm,qsfc,                                      &
            u10,v10,th2,                                        &
            t2,q2,flhc,flqc,                                    &
            lh,gz1oz0,WSPD,br,                                  &
            cpm,ch,rstoch1D,                                    &
            wstar,qstar,                                        &
            ck,cka,cd,cda,                                      &
            psix,psit,psix10,psit2)


        
        character(len=2000), intent(in) :: line
        integer, intent(in) :: out_unit
        integer, intent(in) :: line_number
        
        ! Case-specific locations
        real, intent(out) ::             xlat,xlong
        
        ! Variables to store parsed values
        real, intent(out) ::    dx, xland,          & 
            U1D,V1D,T1D,QV1D,                       & 
            P1D, dz8w1d,RHO1D, U1D2,                &
            V1D2, dz2w1d,                           &
            pblh, znt,PSFCPA,MAVAIL,                &
            TSK,snowh,                              &
            chs,chs2,cqs2,                          &
            ust,ustm,qsfc,                          &
            rmol,zol,mol,                           &
            hfx,qfx,                                &
            u10,v10,th2,                            &
            t2,q2,flhc,flqc,                        &
            lh,gz1oz0,WSPD,br,                      &
            cpm,ch,rstoch1D,                        &
            wstar,qstar,                            &
            ck,cka,cd,cda,                          &
            psix,psit,psix10,psit2      

        integer :: read_stat           
        integer, intent(out) :: itimestep
        
        ! Read values from the line with specified format

        read (line, '(F10.2,F10.2,I10,'                              //&
                    'F10.2,F10.2,'                                   //&
                    'F10.2,F10.2,F10.2,F10.4,'                       //&
                    'F10.2,F10.2,F10.2,F10.2,'                       //&
                    'F10.2,F10.2,'                                   //&
                    'F10.2,F10.2,F10.2,F10.2,'                       //&
                    'F10.2,F10.2,'                                   //&
                    'F10.2,F10.2,F10.2,'                             //&
                    'F10.2,F10.2,F10.2,'                             //&
                    'F10.2,F10.2,F10.2,'                             //&
                    'F10.2,F10.2,'                                   //&
                    'F10.2,F10.2,F10.2,'                             //&
                    'F10.2,F10.2,F10.2,F10.2,'                       //&
                    'F10.2,F10.2,F10.2,F10.2,'                       //&
                    'F10.2,F10.2,F10.2,'                             //&
                    'F10.2,F10.2,'                                   //&
                    'F10.2,F10.2,F10.2,F10.2,'                       //&
                    'F10.2,F10.2,F10.2,F10.2)', iostat=read_stat)      &
                    xlat,xlong,itimestep,                              &
                    dx, xland,                                         & 
                    U1D,V1D,T1D,QV1D,                                  & 
                    P1D, dz8w1d,RHO1D, U1D2,                           &
                    V1D2, dz2w1d,                                      &
                    pblh, znt,PSFCPA,MAVAIL,                           &
                    TSK,snowh,                                         &
                    chs,chs2,cqs2,                                     &
                    ust,ustm,qsfc,                                     &
                    rmol,zol,mol,                                      &
                    hfx,qfx,                                           &
                    u10,v10,th2,                                       &
                    t2,q2,flhc,flqc,                                   &
                    lh,gz1oz0,WSPD,br,                                 &
                    cpm,ch,rstoch1D,                                   &
                    wstar,qstar,                                       &
                    ck,cka,cd,cda,                                     &
                    psix,psit,psix10,psit2                                      
                                      
              
      ! Debug: print the line being read
      write(0,*) "Line length:", len_trim(line)
      write(0,*) "U1D=",U1D, "V1D=", V1D," ZNT", ZNT, "WSPD=",WSPD, "UST=",UST, &
                "HFX=",HFX, "LH=",LH, "tsk=", tsk,"t2=",t2, "xland=",xland
        
    end subroutine process_line_wrf
    !=================================================================================================================
 
    subroutine wrf_test (case, sf_mynn_sfcflux_water, sf_mynn_sfcflux_land, saveoutput)

      character(len=*), intent(in) :: case
      integer, intent(in) :: sf_mynn_sfcflux_water, sf_mynn_sfcflux_land
      logical, intent(in) :: saveoutput
      integer :: iostat, line_num
      integer, parameter :: i=1, j=1
      real, parameter :: xice_threshold=0.5
      character(len=2000) :: input_line
      integer, parameter :: input_unit = 10
      integer, parameter :: output_unit = 20
       
      logical :: compute_flux, compute_diag, redrag, flag_restart
      real, dimension(1) ::                      &
            dx, xland, xice,                     & 
            U1D,V1D,T1D,QV1D,                    & 
            P1D, dz8w1d,RHO1D, U1D2,             &
            V1D2, dz2w1d,                        &
            pblh, znt,PSFCPA,MAVAIL,             &
            TSK,snowh,                           &
            chs,chs2,cqs2,                       &
            rmol,zol,mol,                        &
            ustm,qsfc,qgh,                       &
            u10,v10,th2,                         &
            t2,q2,flhc,flqc,                     &
            gz1oz0,WSPD,br,                      &
            cpm,ch,rstoch1D,                     &
            wstar,qstar,                         &
            ck,cka,cd,cda,                       &
            psix,psit,psix10,psit2,              &
            PSIM,PSIH                                     

      !GFS-related input
      real, dimension(1) ::   sigmaf, shdmax, z0pert, ztpert
      integer :: sfc_z0_type, ivegsrc
      integer, dimension(1) :: vegtype

      real, dimension(1) :: HFLX_mod, HFX_mod , QFLX_mod, QFX_mod, LH_mod, ust_mod, &
                            u10_mod, v10_mod,t2_mod,th2_mod,q2_mod        
      real, dimension(1) :: HFLX, HFX, QFLX, QFX, LH, ust,stress
                  
      integer :: read_stat,itimestep
      real :: xlat,xlong

      ! Variables for error handling
      character(len=512) :: errmsg
      integer :: errflg
      ! Initialize 3D variables
      real, dimension(64) ::  u3d, v3d, t3d, qv3d, p3d, dz8w, th3d, rho3d  ! 64 is random number
      u3d=0.0
      v3d=0.0
      t3d=0.0
      qv3d=0.0
      p3d=0.0
      dz8w=0.0 
      th3d=0.0
      rho3d=0.0

      xice=0.0
      
      ! Initialize error variables
      errmsg = ''
      errflg = 0

      write(*,*) '--- entering wrf_test subroutine ---'   
      call init_mynn_sfc_flags ()

      ! Initialize input data for tests
      call init_input_data_for_test(case,saveoutput)

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
 
          call process_line_wrf(input_line, output_unit, line_num,     &
            xlat,xlong,itimestep,                                      &
            dx(1), xland(1),                                           & 
            U1D(1),V1D(1),T1D(1),QV1D(1),                              & 
            P1D(1), dz8w1d(1),RHO1D(1), U1D2(1),                       &
            V1D2(1), dz2w1d(1),                                        &
            pblh(1), znt(1),PSFCPA(1),MAVAIL(1),                       &
            TSK(1),snowh(1),                                           &
            chs(1),chs2(1),cqs2(1),                                    &
            rmol(1),zol(1),mol(1),                                     &
            hfx(1),qfx(1),                                             &
            ust(1),ustm(1),qsfc(1),                                    &
            u10(1),v10(1),th2(1),                                      &
            t2(1),q2(1),flhc(1),flqc(1),                               &
            lh(1),gz1oz0(1),WSPD(1),br(1),                             &
            cpm(1),ch(1),rstoch1D(1),                                  &
            wstar(1),qstar(1),                                         &
            ck(1),cka(1),cd(1),cda(1),                                 &
            psix(1),psit(1),psix10(1),psit2(1) )                                     
                                      
         ! Construct 3D arrays
          u3d(1)=U1D(1)
          u3d(2)=U1D2(1)
          v3d(1)=V1D(1)
          v3d(2)=V1D2(1)
          t3d(1)=T1D(1)
          qv3d(1)=QV1D(1)
          p3d(1)=P1D(1)
          dz8w(1)=dz8w1d(1)
          dz8w(2)=dz2w1d(1)

         ! Initialize MYNN SFC
          write(*,*) '--- calling  mynnsfc_driver() ---'

          call mynnsfc_init(allowed_to_read=.true.,errmsg=errmsg,errflg=errflg)
          !write(*,*) 'psih_stab_in_wrftests=',psih_stab

          call mynnsfc_driver(u3d=u3d , v3d =v3d , t3d=t3d , qv3d =qv3d, p3d =p3d, dz8w=dz8w,        &
              th3d=th3d, rho3d=rho3d,                                                                &
              !GFS-related input
              sigmaf=sigmaf, vegtype=vegtype, shdmax=shdmax, ivegsrc=ivegsrc,                        &  !intent(in)
              z0pert=z0pert, ztpert=ztpert, redrag=redrag, sfc_z0_type=sfc_z0_type,                  &  !intent(in)
              !2d variables
              psfcpa=psfcpa , chs=CHS, chs2=CHS2, cqs=chs , cqs2=cqs2, cpm=CPM,                      &
              znt=ZNT, ust=ust_mod, ustm=USTM, pblh=pblh, mavail=mavail, zol=ZOL,                    &
              mol=MOL, rmol=RMOL, psim=PSIM , psih=PSIH, xland=XLAND, qgh=QGH,                       &
              hfx=HFX_mod, qfx=QFX_mod, lh=LH_mod, tsk=tsk, flhc=FLHC, flqc=FLQC,                    &
              qsfc=QSFC, u10=u10_mod, v10=v10_mod, th2 =th2_mod, t2=t2_mod ,                         &
              q2=q2_mod, snowh=snowh, gz1oz0=GZ1OZ0, wspd=WSPD, br=br, dx=DX,                        &
              ch=ch, ck=ck, cka=cka, cd=cd, cda=cda ,                                                &
              xice=xice, xice_threshold=xice_threshold,                                              &
              stress=stress , hflx=hflx, qflx=qflx, fm=psix, fh=psit,                                &
              fm10=psix10, fh2=psit2,                                                                &
              tsurf=tsk    ,                                                                         &
              !configuration options
              spp_pbl=spp_pbl, pattern_spp_pbl=pattern_spp_pbl,                                      &
              sf_mynn_sfcflux_water=sf_mynn_sfcflux_water          ,                                 &
              sf_mynn_sfcflux_land=sf_mynn_sfcflux_land           ,                                  &
              isfflx=ISFFLX, restart=restart, cycling=cycling, initflag=1,                           &
              flag_iter=flag_iter, flag_lsm=lsm,                                                     &
              !model information
              itimestep=itimestep,                                                                   &
              ids=1     , ide=1      , jds=1    , jde=1     , kds=1     , kde=1     ,                &
              ims=1     , ime=1      , jms=1    , jme=1     , kms=1     , kme=1     ,                &
              its=1     , ite=1      , jts=1    , jte=1     , kts=1     , kte=1     ,                &
              errmsg=errmsg , errflg=errflg                                                          &
              )


         write(0,*) "flag_iter",flag_iter,"T2=",t2,'chs=',ch,'ust=',UST,'hfx=',hfx_mod,'lh=',lh_mod,'xland=',xland, 'tsk', tsk
         write(0,*) "Read status:", read_stat
         if (saveoutput) then
             open(output_unit, file = './data/output_'//trim(case)//'.txt')
             write(output_unit,'(I5, F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)')  &
                  itimestep,t2_mod,q2_mod,th2_mod,u10_mod,v10_mod,hfx_mod,lh_mod,ust_mod,pblh
          end if

      end do
    end subroutine wrf_test

end module module_sf_mynnsfc_wrf_tests           
