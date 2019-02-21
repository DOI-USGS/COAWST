program da_verif_grid 
   !---------------------------------------------------------------------------
   ! History:
   !
   !  Abstract: Program to calculate statistics for multiple experiments for 
   !            verification analayis/forecasts against any desired analysis
   !
   !  Author:   Syed RH Rizvi     NCAR/MMM         05/25/2006
   !  Udates:
   !            Syed RH Rizvi     NCAR/MMM         06/05/2009
   !            a) Added verification for wind vector & geopotentials 
   !            b) For vertical profiles, added the choice to display the desired
   !               number of pressure levels 
   !---------------------------------------------------------------------------

   use da_verif_grid_control, only : control_main, control_times, control_vars, sub_domain, &
       max_3d_variables, max_2d_variables,num_vert_levels,verification_file_string,&
       missing,namelist_unit,time_series_unit, time_series_2d, profile_time_series_3d,&
       filename, stime, etime, hstart, hend, hdate, date, pdate, vert_levels, &
       nx, ny, nz, io_status,  debug1, debug2, verify_its_own_analysis, &
       num_verifying_experiments, verify_forecast_hour, domain, control_exp_dir, verif_dirs, &
       out_dirs,start_year, end_year, start_month, end_month, start_day, end_day, &
       start_hour, end_hour,start_minutes, end_minutes, start_seconds, end_seconds,interval_hour, &
       num3dvar, num2dvar, var3d, var2d, num_scores, score_names, vertical_type, &
       istart, iend, jstart, jend, unit_all, unit_land, unit_water

   use da_netcdf_interface, only : da_get_dims_cdf, da_get_gl_att_int_cdf, da_get_gl_att_real_cdf, &
      da_get_var_3d_real_cdf, da_get_var_2d_real_cdf, da_get_var_2d_int_cdf

   implicit none

   character (len=512) :: control_file, verif_file
   character (len=512) :: out_dir, outname_all, outname_land, outname_water

   integer, parameter                   :: imiss = -99
   real, parameter                      :: rmiss = -99.99
   integer                              :: time_loop_count
   integer                              :: time(6), ptime(6)
   integer                              :: nx1, ny1, nz1
   integer                              :: nx2, ny2, nz2
   integer                              :: i,k
   integer                              :: ivar, iexp, iscore
   character (len=10)                   :: sdate
   character (len=20)                   :: file_string, domain_string, out_hr
   logical, allocatable,dimension(:)    :: first_score
   logical                              :: exist_file1, exist_file2

   real, allocatable, dimension(:,:,:)  :: data_out1, data_out2
   real, allocatable, dimension(:,:,:)  :: u1, u2, v1, v2       
 
   real, allocatable, dimension(:,:,:)  :: sum3d, asum3d, sqr3d, diff, absdiff, sqdiff
   real, allocatable, dimension(:,:)    :: score_avg_prof
   real, allocatable, dimension(:)      :: avg_prof
   real, allocatable, dimension(:,:)    :: landmask  ! 1 for land, 0 for water
   integer, parameter                   :: island = 1, iswater = 0
   integer, allocatable, dimension(:)   :: num_counter
   integer, allocatable, dimension(:,:,:)  :: mask_count
   real, allocatable, dimension(:,:,:)  :: mask_stats ! 1st dimension: 1:bias, 2:rmse, 3:abias
                                                      ! 2nd dimension: 1:all, 2:land, 3:water
                                                      ! 3rd dimension: number of vertical levels
   
!----------------------------------------------------------------------------
   verify_forecast_hour = 0
!----------------------------------------------------------------------------
   debug1 = .false.
   debug2 = .false.
   vertical_type = 'p'

!--3D need update
   num3dvar=max_3d_variables
   var3d(1)='U'
   var3d(2)='V'
   var3d(3)='T'
   var3d(4)='QVAPOR'
   var3d(5)='Z'
   var3d(6)='WV'

!--2D need update
   num2dvar=max_2d_variables
   var2d(1)='SLP'
   var2d(2)='PSFC'
   var2d(3)='U10M'
   var2d(4)='V10M'
   var2d(5)='T2M'
   var2d(6)='Q2M'
   var2d(7)='MU'

!--Score names
   num_scores = 3
   score_names(1) = 'BIAS'
   score_names(2) = 'RMSE'
   score_names(3) = 'ABIAS'
                                     
   domain = 1
   verification_file_string = 'wrfout'
   out_hr ='_00'
!---------------------------------------------------------------------
!  Read namelist
!----------------------------------------------------------------------------
   io_status = 0

                                                                      
   open(unit = namelist_unit, file = 'namelist.in', &
          status = 'old' , access = 'sequential', &
          form   = 'formatted', action = 'read', &
          iostat = io_status )

                                                                   
   if(io_status /= 0) then
      print *, 'Error to open namelist.in file: '
   else
      read(unit=namelist_unit, nml = control_main , iostat = io_status)

      if(io_status /= 0) then
         print *, 'Error to read control_main. Stopped.'
         stop
      endif

                                                                   
      read(unit=namelist_unit, nml = control_times , iostat = io_status )
      
      if(io_status /= 0) then
         print *, 'Error to read control_times Stopped.'
         stop
      endif
   
      read(unit=namelist_unit, nml = control_vars , iostat = io_status )
     
      if(io_status /= 0) then
         print *, 'Error to read control_vars Stopped.'
         stop
      endif
   
      read(unit=namelist_unit, nml = sub_domain , iostat = io_status )
     
      if(io_status /= 0) then
         print *, 'Error reading sub_domain'
         istart = 1
         iend   = 10000
         jstart = 1
         jend   = 10000
         !stop
      endif
   
      close(unit=namelist_unit)
   endif
!----------------------------------------------------------------------------
!---------------------------------------------------------------------
!
    write(domain_string, fmt ='("_d",i2.2,"_")') domain
    file_string = trim(verification_file_string)//trim(domain_string) 
    write(out_hr, fmt ='("_",i2.2)') verify_forecast_hour
    allocate(first_score(num_scores))
!---------------------------------------------------------------------

   stime(1) = start_year
   stime(2) = start_month
   stime(3) = start_day  
   stime(4) = start_hour 
   stime(5) = start_minutes
   stime(6) = start_seconds


   call build_hdate(hstart, stime )

   etime(1) = end_year
   etime(2) = end_month
   etime(3) = end_day  
   etime(4) = end_hour 
   etime(5) = end_minutes
   etime(6) = end_seconds

   call build_hdate(hend, etime )

                                                                   
   if ( hend < hstart ) then
      print*, '****************************************************************'
      print*, 'End time is before the start time'
      print*, ' Start time = ', hstart ,'  & End time = ', hend
      print*, '****************************************************************'
      stop
   endif

   hdate = hstart 
   call build_hdate(sdate, stime )

   filename = trim(file_string)//hdate
   
!---------------------------------------------------------------------
   loop_verif_exp : do iexp = 1, num_verifying_experiments

      first_score  = .true.
      out_dir = trim(out_dirs(iexp))//'/'
!---------------------------------------------------------------------
!--For 3D variables
!---------------------------------------------------------------------
      allocate( avg_prof(num_vert_levels) )
      allocate( num_counter(num_vert_levels) )
      allocate( score_avg_prof(num_scores, num_vert_levels) )

      loop_3d : do ivar=1,num3dvar

!  Open profile units

        profile_time_series_3d = trim(out_dir)//trim(var3d(ivar))//'_time_series'//trim(out_hr) 
        open(time_series_unit, file=profile_time_series_3d,form='formatted', &
                               status='unknown')
       
        outname_all   = trim(profile_time_series_3d)//'_sub'
        outname_land  = trim(profile_time_series_3d)//'_sub_land'
        outname_water = trim(profile_time_series_3d)//'_sub_water'
        open(unit_all,   file=trim(outname_all),  form='formatted', status='unknown')
        open(unit_land,  file=trim(outname_land), form='formatted', status='unknown')
        open(unit_water, file=trim(outname_water),form='formatted', status='unknown')
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
        time_loop_count = 0
        hdate = hstart
        time = stime
     
        time_loop_3d : do 
              
!----------------------------------------------------------------------------
           call build_hdate(hdate, time)
           print*,' processing exp: ',iexp,' 3d var: ',trim(var3d(ivar)),' for time: ',hdate

           if ( hdate > hend ) exit time_loop_3d

           call build_hdate(date, time )
           ptime = time
           call advance_date(ptime,-verify_forecast_hour) 
           call build_hdate(pdate, ptime)


           filename = trim(file_string)//hdate

           if( verify_its_own_analysis) then
           control_file = trim(verif_dirs(iexp))//'/'//date//'/'//trim(filename) 
           else
           control_file = trim(control_exp_dir)//'/'//date//'/'//trim(filename) 
           endif

           verif_file = trim(verif_dirs(iexp))//'/'//pdate//'/'//trim(filename)

           inquire(file=trim(control_file), exist=exist_file1)
           inquire(file=trim(verif_file), exist=exist_file2)
           if ( (.not. exist_file1) .or. (.not. exist_file2) ) then
              avg_prof       = rmiss
              num_counter    = imiss
              score_avg_prof = rmiss
              call write_profile(date, score_avg_prof, num_counter, num_vert_levels, num_scores, &
                                 time_series_unit) 
              call write_profile(date, score_avg_prof, num_counter, num_vert_levels, num_scores, unit_all)
              call write_profile(date, score_avg_prof, num_counter, num_vert_levels, num_scores, unit_land)
              call write_profile(date, score_avg_prof, num_counter, num_vert_levels, num_scores, unit_water)
              call advance_date(time,interval_hour) 
              if ( .not. exist_file1 ) print*, trim(control_file), ' file not found'
              if ( .not. exist_file2 ) print*, trim(verif_file), ' file not found'
              print*, 'skipping date ', date
              cycle time_loop_3d
           end if

!          print*,' control_file  : ', trim(control_file)
!          print*,' verif_file    : ', trim(verif_file)

!  Get the dimensions of the both files and check
           call get_dimensions(control_file,nx1,ny1,nz1)
           call get_dimensions(verif_file,nx2,ny2,nz2)

           if ( nx1 /= nx2 .or. ny1 /= ny2 .or. nz1 /= nz2 ) then
              print*, '********************************************************'
              print*, 'Dimension mismatch between files of the experiments ....' 
              print*, '********************************************************'
              stop
           else
             nx = nx1
             ny = ny1
             nz = nz1
             if (time_loop_count == 0 ) then
               allocate( sum3d(nx,ny,num_vert_levels))
               allocate( asum3d(nx,ny,num_vert_levels))
               allocate( sqr3d(nx,ny,num_vert_levels))
               sum3d = 0.0
               asum3d = 0.0
               sqr3d = 0.0

               if ( ivar == 1 ) then
                  allocate(landmask(nx,ny))
                  call da_get_var_2d_real_cdf(verif_file,"LANDMASK", landmask, nx, ny, 1, .false.) 
                  istart = max(1, istart)
                  iend   = min(nx, iend)
                  jstart = max(1, jstart)
                  jend   = min(ny, jend)
               end if

             endif
           endif
           allocate(diff(nx,ny,num_vert_levels))
           allocate(absdiff(nx,ny,num_vert_levels))
           allocate(sqdiff(nx,ny,num_vert_levels))
!  first, get control data
           if( trim(var3d(ivar)) .eq. "WV" ) then
           allocate ( u1 (nx, ny, num_vert_levels) )
           allocate ( u2 (nx, ny, num_vert_levels) )
           allocate ( v1 (nx, ny, num_vert_levels) )
           allocate ( v2 (nx, ny, num_vert_levels) )
           call compute_wind_3d( control_file, nx, ny, nz, num_vert_levels, 1, &
                                 vert_levels, vertical_type, missing, u1, v1, debug1 )
           call compute_wind_3d( verif_file, nx, ny, nz, num_vert_levels, 1, &
                                 vert_levels, vertical_type, missing, u2, v2, debug1 )
           call get_diffs_wind(u1, v1, u2, v2, diff, absdiff, sqdiff, nx, ny, &
                          num_vert_levels, missing)
           deallocate(u1, v1, u2, v2)
           else
           allocate ( data_out1 (nx, ny, num_vert_levels) )
           allocate ( data_out2 (nx, ny, num_vert_levels) )
           call compute_data_3d( control_file, trim(var3d(ivar)), len_trim(var3d(ivar)),    &
                                nx, ny, nz, num_vert_levels, 1, vert_levels, &
                                vertical_type, missing, data_out1, debug1 )
!  second, get verifying data
           call compute_data_3d( verif_file, trim(var3d(ivar)), len_trim(var3d(ivar)),    &
                                nx, ny, nz, num_vert_levels, 1, vert_levels, &
                                vertical_type, missing, data_out2, debug2 )

           call get_diffs(data_out1, data_out2, diff, absdiff, sqdiff, nx, ny,  &
                          num_vert_levels, missing)
           deallocate(data_out1)
           deallocate(data_out2)
           end if

           do iscore = 1, num_scores
             if ( trim(score_names(iscore)) == 'BIAS' ) then
                call domain_average( diff, avg_prof, num_counter, nx, ny, num_vert_levels, missing,0)
             elseif ( trim(score_names(iscore)) == 'RMSE' ) then
                call domain_average( sqdiff, avg_prof, num_counter, nx, ny, num_vert_levels, missing,1)
             elseif ( trim(score_names(iscore)) == 'ABIAS' ) then
                call domain_average( absdiff, avg_prof, num_counter, nx, ny, num_vert_levels, missing,0)
             endif
             score_avg_prof(iscore,:) = avg_prof(:)
           enddo
           call write_profile(date, score_avg_prof, num_counter, num_vert_levels, num_scores, &
                              time_series_unit) 

           allocate(mask_stats(3,3,nz))
           allocate(mask_count(3,3,nz))
           call mask_domain_average( diff,    mask_stats(1,:,:), mask_count(1,:,:), nx, ny, num_vert_levels, missing,0)
           call mask_domain_average( sqdiff,  mask_stats(2,:,:), mask_count(2,:,:), nx, ny, num_vert_levels, missing,1)
           call mask_domain_average( absdiff, mask_stats(3,:,:), mask_count(3,:,:), nx, ny, num_vert_levels, missing,0)
           call write_profile(date, mask_stats(:,1,:), mask_count(1,1,:), num_vert_levels, 3, unit_all) 
           call write_profile(date, mask_stats(:,2,:), mask_count(1,2,:), num_vert_levels, 3, unit_land) 
           call write_profile(date, mask_stats(:,3,:), mask_count(1,3,:), num_vert_levels, 3, unit_water) 
           deallocate(mask_stats)
           deallocate(mask_count)
           
           call get_sum(sum3d,diff,nx,ny,num_vert_levels,missing)
           call get_sum(asum3d,absdiff,nx,ny,num_vert_levels,missing)
           call get_sum(sqr3d,sqdiff,nx,ny,num_vert_levels,missing)

           deallocate(diff)
           deallocate(absdiff)
           deallocate(sqdiff)

           time_loop_count = time_loop_count + 1

           call advance_date(time,interval_hour) 

        enddo  time_loop_3d       ! time loop over

        close(time_series_unit)
        close(unit_all)
        close(unit_land)
        close(unit_water)
        deallocate(sum3d)
        deallocate(asum3d)
        deallocate(sqr3d)

     enddo loop_3d
     print*, ' successful completion of loop_3d '

     deallocate( avg_prof )
     deallocate( num_counter )
     deallocate( score_avg_prof )

!--------------------------------------------------------------------------------
!--Loop For 2D variables
!--------------------------------------------------------------------------------
     allocate( avg_prof(1) )
     allocate( num_counter(1) )
     allocate( score_avg_prof(num_scores, 1) )

     loop_2d : do ivar = 1, num2dvar

!  Open profile units

        time_series_2d  = trim(out_dir)//trim(var2d(ivar))//'_time_series'//trim(out_hr) 
        open(time_series_unit, file=time_series_2d,form='formatted', &
                               status='unknown')
       
        outname_all   = trim(time_series_2d)//'_sub'
        outname_land  = trim(time_series_2d)//'_sub_land'
        outname_water = trim(time_series_2d)//'_sub_water'
        open(unit_all,   file=trim(outname_all),  form='formatted', status='unknown')
        open(unit_land,  file=trim(outname_land), form='formatted', status='unknown')
        open(unit_water, file=trim(outname_water),form='formatted', status='unknown')
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
        
        time_loop_count = 0
        hdate = hstart
        time = stime
     
        time_loop_2d : do 
              
!----------------------------------------------------------------------------
           call build_hdate(hdate, time )
           print*,' processing exp: ',iexp,' 2d var: ',trim(var2d(ivar)),' for time: ',hdate

           if ( hdate > hend ) exit time_loop_2d

           call build_hdate(date, time )
           ptime = time
           call advance_date(ptime,-verify_forecast_hour) 
           call build_hdate(pdate, ptime)

           filename = trim(file_string)//hdate
           if( verify_its_own_analysis) then
           control_file = trim(verif_dirs(iexp))//'/'//date//'/'//trim(filename) 
           else
           control_file = trim(control_exp_dir)//'/'//date//'/'//trim(filename) 
           endif

           verif_file = trim(verif_dirs(iexp))//'/'//pdate//'/'//trim(filename)

           inquire(file=trim(control_file), exist=exist_file1)
           inquire(file=trim(verif_file), exist=exist_file2)
           if ( (.not. exist_file1) .or. (.not. exist_file2) ) then
              avg_prof       = rmiss
              num_counter    = imiss
              score_avg_prof = rmiss
              call write_profile(date, score_avg_prof, num_counter, 1, num_scores, &
                                 time_series_unit) 
              call write_profile(date, score_avg_prof, num_counter, 1, num_scores, unit_all)
              call write_profile(date, score_avg_prof, num_counter, 1, num_scores, unit_land)
              call write_profile(date, score_avg_prof, num_counter, 1, num_scores, unit_water)
              call advance_date(time,interval_hour) 
              if ( .not. exist_file1 ) print*, trim(control_file), ' file not found'
              if ( .not. exist_file2 ) print*, trim(verif_file), ' file not found'
              print*, 'skipping date ', date
              cycle time_loop_2d
           end if

!
!  Get the dimensions of the both files and check
           call get_dimensions(control_file,nx1,ny1,nz1)
           call get_dimensions(verif_file,nx2,ny2,nz2)

           if ( nx1 /= nx2 .or. ny1 /= ny2 .or. nz1 /= nz2 ) then
              print*, '********************************************************'
              print*, 'Dimension mismatch between files of the experiments ....' 
              print*, '********************************************************'
              stop
           else
             nx = nx1
             ny = ny1
             nz = nz1
             if (time_loop_count == 0 ) then
               allocate( sum3d(nx,ny,1))
               allocate( asum3d(nx,ny,1))
               allocate( sqr3d(nx,ny,1))
               sum3d = 0.0
               asum3d = 0.0
               sqr3d = 0.0
             endif
           endif

           allocate(data_out1(nx, ny, 1))
           allocate(data_out2(nx, ny, 1))
 
           call g_output_2d (control_file, 1, trim(var2d(ivar)), len_trim(var2d(ivar)), &
                             nx, ny, nz, data_out1, debug1)

           call g_output_2d (verif_file, 1, trim(var2d(ivar)), len_trim(var2d(ivar)), &
                             nx, ny, nz, data_out2, debug2)

           allocate(diff(nx,ny,1))
           allocate(absdiff(nx,ny,1))
           allocate(sqdiff(nx,ny,1))
           call get_diffs(data_out1, data_out2, diff, absdiff, sqdiff, nx, ny, 1, missing)
           deallocate(data_out1)
           deallocate(data_out2)

           do iscore = 1, num_scores
             if ( trim(score_names(iscore)) == 'BIAS' ) then
                call domain_average( diff, avg_prof, num_counter, nx, ny, 1, missing,0)
             elseif ( trim(score_names(iscore)) == 'RMSE' ) then
                call domain_average( sqdiff, avg_prof, num_counter, nx, ny, 1, missing,1)
             elseif ( trim(score_names(iscore)) == 'ABIAS' ) then
                call domain_average( absdiff, avg_prof, num_counter, nx, ny, 1, missing,0)
             endif
             score_avg_prof(iscore,:) = avg_prof(:)
           enddo
           call write_profile(date, score_avg_prof, num_counter, 1, num_scores, &
                              time_series_unit)

           allocate(mask_stats(3,3,1))
           allocate(mask_count(3,3,1))
           call mask_domain_average( diff,    mask_stats(1,:,:), mask_count(1,:,:), nx, ny, 1, missing,0)
           call mask_domain_average( sqdiff,  mask_stats(2,:,:), mask_count(2,:,:), nx, ny, 1, missing,1)
           call mask_domain_average( absdiff, mask_stats(3,:,:), mask_count(3,:,:), nx, ny, 1, missing,0)
           call write_profile(date, mask_stats(:,1,:), mask_count(1,1,:), 1, 3, unit_all) 
           call write_profile(date, mask_stats(:,2,:), mask_count(1,2,:), 1, 3, unit_land) 
           call write_profile(date, mask_stats(:,3,:), mask_count(1,3,:), 1, 3, unit_water) 
           deallocate(mask_stats)
           deallocate(mask_count)
           
           call get_sum(sum3d,diff,nx,ny,1,missing)
           call get_sum(asum3d,absdiff,nx,ny,1,missing)
           call get_sum(sqr3d,sqdiff,nx,ny,1,missing)

           deallocate(diff)
           deallocate(absdiff)
           deallocate(sqdiff)
           time_loop_count = time_loop_count + 1

           call advance_date(time,interval_hour) 


        enddo  time_loop_2d
        
        close(time_series_unit)
        close(unit_all)
        close(unit_land)
        close(unit_water)
        deallocate(sum3d)
        deallocate(asum3d)
        deallocate(sqr3d)

     enddo loop_2d

     deallocate( avg_prof )
     deallocate( num_counter )
     deallocate( score_avg_prof )
     deallocate( landmask )

     print*, ' successful completion of loop_2d '
   print*,' Finished Experiment : ', trim(verif_dirs(iexp))
   enddo loop_verif_exp

!-----------------------------------------------------

   contains  
!-----------------------------------------------------
   subroutine advance_date( time, delta )

      implicit none

      integer, intent(inout) :: time(6)         
      integer, intent(in)    :: delta

      integer                :: ccyy, mm, dd, hh
      integer, dimension(12) :: mmday
!      character(len=10) :: ccyymmddhh

!-----------------------------------------------------
      mmday = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      mmday(2) = 28
!-----------------------------------------------------
    ccyy = time(1)
    mm   = time(2)
    dd   = time(3)
    hh   = time(4)
!-----------------------------------------------------

    hh = hh + delta

    do while (hh < 0)
       hh = hh + 24
       call change_date ( ccyy, mm, dd, -1 )
    end do

    do while (hh > 23)
       hh = hh - 24
       call change_date ( ccyy, mm, dd, 1 )
    end do

!    write(ccyymmddhh(1:10), fmt='(i4, 3i2.2)')  ccyy, mm, dd, hh
    time(1) = ccyy
    time(2) = mm
    time(3) = dd
    time(4) = hh

    end subroutine advance_date
!---------------------------------------------------------------------------
    subroutine change_date ( ccyy, mm, dd, delta)
      integer, intent(inout) :: ccyy, mm, dd
      integer, intent(in) :: delta
!
      integer, dimension(12) :: mmday
      mmday = (/31,28,31,30,31,30,31,31,30,31,30,31/)

      mmday(2) = 28

      if (mod(ccyy,4) == 0) then
         mmday(2) = 29

         if ( mod(ccyy,100) == 0) then
            mmday(2) = 28
         endif

         if(mod(ccyy,400) == 0) then
            mmday(2) = 29
         end if
      endif

      dd = dd + delta

      if(dd == 0) then
         mm = mm - 1

         if(mm == 0) then
            mm = 12         
            ccyy = ccyy - 1
         endif

         dd = mmday(mm)
      elseif ( dd .gt. mmday(mm) ) then
         dd = 1
         mm = mm + 1
         if(mm > 12 ) then
            mm = 1
            ccyy = ccyy + 1
         end if
      end if
   end subroutine change_date

   subroutine build_hdate(hdate, time)

! PURPOSE:
!      From the Year, Month, Day, Hour, Minute, and Second values,
!      creates a 19-character string representing the date, in the
!      format:  "YYYY-MM-DD hh:mm:ss"

! INPUT:
      integer, intent(in) :: time(6) ! all time specs in it
! OUTPUT:
      character*(*), intent(out) :: hdate ! 'YYYY-MM-DD hh:mm:ss'

! LOCAL:
      integer iyr     ! year (e.g., 1997, 2001)
      integer imo     ! month (01 - 12)
      integer idy     ! day of the month (01 - 31)
      integer ihr     ! hour (00-23)
      integer imi     ! minute (00-59)
      integer isc     ! second (00-59)
!
!      integer i  ! Loop counter.
      integer hlen ! Length of hdate string

      hlen = len(hdate)
      iyr = time(1)
      imo = time(2)
      idy = time(3)
      ihr = time(4)
      imi = time(5)
      isc = time(6)

      if (hlen.eq.19) then
         write(hdate,19) iyr, imo, idy, ihr, imi, isc
 19      format(i4,'-',i2.2,'-',i2.2,'_',i2.2,':',i2.2,':',i2.2)

      elseif (hlen.eq.16) then
         write(hdate,16) iyr, imo, idy, ihr, imi
 16      format(i4,'-',i2.2,'-',i2.2,'_',i2.2,':',i2.2)

      elseif (hlen.eq.13) then
         write(hdate,13) iyr, imo, idy, ihr
 13      format(i4,'-',i2.2,'-',i2.2,'_',i2.2)

      elseif (hlen.eq.10) then
         write(hdate,10) iyr, imo, idy, ihr
 10      format(i4,i2.2,i2.2,i2.2)
      endif

      return
      end subroutine build_hdate

!----------------------------------------------------------------------------------
  subroutine write_profile(date, profile, counter, nlevel, nscore, out_unit)
  
  integer, intent(in)                 :: nlevel, nscore, out_unit
  real, intent(in), dimension(:,:)    :: profile
  integer, intent(in), dimension(:)   :: counter
  character (len=10), intent(in)      :: date
  write(out_unit,fmt='(a10,1x,100(i8,1x,3(f14.8,1x)))') date,  &
                              (counter(k), (profile(i,k),i=1,nscore),k=1,nlevel)

  end subroutine write_profile
  
!---------------------------------------------------------------------------------
  subroutine time_calc( time, timestamp, datestamp, debug , tdef,it)

  implicit none

  character (len=19), intent(in) :: time
  character (len=35), intent(inout)             :: tdef
  integer, intent(out)           :: timestamp, datestamp
  logical, intent(in)            :: debug

   integer, intent(in) :: it
  integer :: hours, minutes, seconds, year, month, day,hour1,hourint
  integer :: mins1,minsint

  save hourint 
  save minsint

  read(time(18:19),*) seconds
  read(time(15:16),*) minutes
  read(time(12:13),*) hours
  read(time(1:4),*)   year
  read(time(6:7),*)   month
  read(time(9:10),*)  day

  if(debug) write(6,*) ' day, month, year, hours, minutes, seconds '
  if(debug) write(6,*) day, month, year, hours, minutes, seconds 

  if ( it == 1) then
    write (tdef(19:20),'(i2)') hours
    if ( day < 10 ) then
      write (tdef(23:23),'(i1)') day
    else
      write (tdef(22:23),'(i2)') day
    endif
    write (tdef(27:30),'(i4)') year
    if (month == 1) write (tdef(24:26),'(a3)') 'jan'
    if (month == 2) write (tdef(24:26),'(a3)') 'feb'
    if (month == 3) write (tdef(24:26),'(a3)') 'mar'
    if (month == 4) write (tdef(24:26),'(a3)') 'apr'
    if (month == 5) write (tdef(24:26),'(a3)') 'may'
    if (month == 6) write (tdef(24:26),'(a3)') 'jun'
    if (month == 7) write (tdef(24:26),'(a3)') 'jul'
    if (month == 8) write (tdef(24:26),'(a3)') 'aug'
    if (month == 9) write (tdef(24:26),'(a3)') 'sep'
    if (month ==10) write (tdef(24:26),'(a3)') 'oct'
    if (month ==11) write (tdef(24:26),'(a3)') 'nov'
    if (month ==12) write (tdef(24:26),'(a3)') 'dec'
    hour1=hours
    mins1=minutes
  elseif ( it == 2) then
    hourint = abs(hours-hour1)
    minsint = abs(minutes-mins1)
    if (hourint == 0 ) then
      if (minsint == 0 ) minsint = 1
      if(debug) write(6,*) "interval is",minsint
      write (tdef(34:35),'(a2)') "mn"
      write (tdef(32:33),'(i2)') minsint
      if(debug) write(6,*) "TDEF is",tdef
    else
      if(debug) write(6,*) "Interval is",hourint
      write (tdef(32:33),'(i2)') hourint
      if(debug) write(6,*) "TDEF is",tdef
    endif
  endif

  timestamp = seconds+100*minutes+10000*hours

  if((year > 1800) .and. (year < 2000)) year = year-1900
  if((year >= 2000)) year = year-2000

  if(month >= 2) day = day+31  ! add january
  if(month >= 3) day = day+28  ! add february
  if(month >= 4) day = day+31  ! add march
  if(month >= 5) day = day+30  ! add april
  if(month >= 6) day = day+31  ! add may
  if(month >= 7) day = day+30  ! add june
  if(month >= 8) day = day+31  ! add july
  if(month >= 9) day = day+31  ! add august
  if(month >= 10) day = day+30 ! add september
  if(month >= 11) day = day+31 ! add october
  if(month >= 12) day = day+30 ! add november
  if((month > 2) .and. (mod(year,4) == 0)) day = day+1  ! get leap year day

  datestamp = (year)*1000 + day
!  datestamp = (year+2100)*1000 + day

  if(debug) then
    write(6,*) ' time, timestamp, datestamp ',time(1:19),timestamp,datestamp
  endif

  end subroutine time_calc

!-------------------------------------------------------------------------

  subroutine g_output_3d (file, file_time_index, var, length_var,            &
                          nx, ny, nz, data_out, debug)
  implicit none

  character (len=*), intent(in)             ::   file
  integer, intent(in)                       ::   file_time_index
  character (len=*), intent(in)             ::   var
  integer, intent(in)                       ::   length_var
  integer , intent(in)                      ::   nx, ny, nz           
  real, intent(out), dimension(:,:,:)       ::   data_out
  logical, intent(in)                                   ::   debug
  real,    allocatable, dimension(:,:,:)    ::   data_tmp, data_tmp2
  real,    allocatable, dimension(:,:,:)    ::   u, v
  real,    allocatable, dimension(:,:)      ::   xlat, xlon
!  real,    allocatable, dimension(:,:,:)    ::   z 
  real,    allocatable, dimension(:,:,:)    ::   ph, phb  
  real,    allocatable, dimension(:,:,:)    ::   p, pb  
  real,    allocatable, dimension(:,:,:)    ::   t, qv 
  integer                                   ::   map_proj
  real                                      ::   cen_lon, truelat1, truelat2


   REAL    , PARAMETER :: g            = 9.81  ! acceleration due to gravity (m {s}^-2)
   REAL    , PARAMETER :: r_d          = 287.
   REAL    , PARAMETER :: r_v          = 461.6
   REAL    , PARAMETER :: cp           = 7.*r_d/2.
   REAL    , PARAMETER :: cv           = cp-r_d
   REAL    , PARAMETER :: cliq         = 4190.
   REAL    , PARAMETER :: cice         = 2106.
   REAL    , PARAMETER :: psat         = 610.78
   REAL    , PARAMETER :: rcv          = r_d/cv
   REAL    , PARAMETER :: rcp          = r_d/cp
   REAL    , PARAMETER :: c2           = cp * rcv
   REAL    , PARAMETER :: T0           = 273.16

   REAL    , PARAMETER :: p1000mb      = 100000.
   REAL    , PARAMETER :: cpovcv       = cp/(cp-r_d)
   REAL    , PARAMETER :: cvovcp       = 1./cpovcv
!   REAL   :: pp

  if(debug) then
    write(6,*) ' calculations for variable ',var
  end if

       if(var == 'U' ) then

          allocate ( data_tmp(nx+1,ny,nz) )
          call da_get_var_3d_real_cdf( file,"U", data_tmp, nx+1, ny, nz,            &
                                file_time_index, debug  )
          data_out = 0.5*(data_tmp(1:nx,:,:)+data_tmp(2:nx+1,:,:))

          deallocate ( data_tmp )

  else if(var == 'V' ) then

          allocate ( data_tmp(nx,ny+1,nz) )
          call da_get_var_3d_real_cdf( file,"V", data_tmp, nx, ny+1, nz,            &
                                file_time_index, debug  )
          data_out = 0.5*(data_tmp(:,1:ny,:)+data_tmp(:,2:ny+1,:))
          deallocate ( data_tmp )

  else if(var == 'UMET' ) then

          call da_get_gl_att_int_cdf ( file, 'MAP_PROJ', map_proj, debug )

          IF ( map_proj == 1  .OR.  map_proj == 2 ) THEN

              allocate (        u(nx,ny,nz)   )
              allocate (        v(nx,ny,nz)   )
              allocate (     xlat(nx,ny)      )             
              allocate (     xlon(nx,ny)      )             
    
              allocate ( data_tmp(nx+1,ny,nz) )
              call da_get_var_3d_real_cdf( file,"U", data_tmp, nx+1, ny, nz,        &
                                    file_time_index, debug  )
              u = 0.5*(data_tmp(1:nx,:,:)+data_tmp(2:nx+1,:,:))
              deallocate ( data_tmp )
    
              allocate ( data_tmp(nx,ny+1,nz) )
              call da_get_var_3d_real_cdf( file,"V", data_tmp, nx, ny+1, nz,        &
                                    file_time_index, debug  )
              v = 0.5*(data_tmp(:,1:ny,:)+data_tmp(:,2:ny+1,:))
              deallocate ( data_tmp )
 
              call da_get_gl_att_real_cdf( file, 'STAND_LON', cen_lon, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT1', truelat1, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT2', truelat2, debug )
              call da_get_var_2d_real_cdf( file, 'XLAT', xlat,nx,ny, 1,debug )
              call da_get_var_2d_real_cdf( file, 'XLONG',xlon,nx,ny, 1,debug )

              call rotate_wind (u,v,nx,ny,nz,var,                                &
                                map_proj,cen_lon,xlat,xlon,                      &
                                truelat1,truelat2,data_out)

              deallocate ( xlat )             
              deallocate ( xlon )             
              deallocate ( u    )
              deallocate ( v    )

          ELSE

              allocate ( data_tmp(nx+1,ny,nz) )
              call da_get_var_3d_real_cdf( file,"U", data_tmp, nx+1, ny, nz,        &
                                    file_time_index, debug  )
              data_out = 0.5*(data_tmp(1:nx,:,:)+data_tmp(2:nx+1,:,:))
              deallocate ( data_tmp )
    
          ENDIF

  else if(var == 'VMET' ) then

          call da_get_gl_att_int_cdf ( file, 'MAP_PROJ', map_proj, debug )

          IF ( map_proj == 1  .OR.  map_proj == 2 ) THEN

              allocate (        u(nx,ny,nz)   )
              allocate (        v(nx,ny,nz)   )
              allocate (     xlat(nx,ny)      )             
              allocate (     xlon(nx,ny)      )             
    
              allocate ( data_tmp(nx+1,ny,nz) )
              call da_get_var_3d_real_cdf( file,"U", data_tmp, nx+1, ny, nz,        &
                                    file_time_index, debug  )
              u = 0.5*(data_tmp(1:nx,:,:)+data_tmp(2:nx+1,:,:))
              deallocate ( data_tmp )
    
              allocate ( data_tmp(nx,ny+1,nz) )
              call da_get_var_3d_real_cdf( file,"V", data_tmp, nx, ny+1, nz,        &
                                    file_time_index, debug  )
              v = 0.5*(data_tmp(:,1:ny,:)+data_tmp(:,2:ny+1,:))
              deallocate ( data_tmp )
 
              call da_get_gl_att_real_cdf( file, 'STAND_LON', cen_lon, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT1', truelat1, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT2', truelat2, debug )
              call da_get_var_2d_real_cdf( file, 'XLAT', xlat,nx,ny, 1,debug )
              call da_get_var_2d_real_cdf( file, 'XLONG',xlon,nx,ny, 1,debug )

              call rotate_wind (u,v,nx,ny,nz,var,                                &
                                map_proj,cen_lon,xlat,xlon,                      &
                                truelat1,truelat2,data_out)

              deallocate ( xlat )             
              deallocate ( xlon )             
              deallocate ( u    )
              deallocate ( v    )

          ELSE

              allocate ( data_tmp(nx,ny+1,nz) )
              call da_get_var_3d_real_cdf( file,"V", data_tmp, nx, ny+1, nz,        &
                                    file_time_index, debug  )
              data_out = 0.5*(data_tmp(:,1:ny,:)+data_tmp(:,2:ny+1,:))
              deallocate ( data_tmp )
 
          ENDIF

  else if(var == 'W' ) then

          allocate ( data_tmp(nx,ny,nz+1) )
          call da_get_var_3d_real_cdf( file,"W", data_tmp, nx, ny, nz+1,            &
                                file_time_index, debug  )
          data_out = 0.5*(data_tmp(:,:,1:nz)+data_tmp(:,:,2:nz+1))
          deallocate ( data_tmp )
 
  else if(var == 'P' ) then

          allocate (  p(nx,ny,nz) )
          allocate ( pb(nx,ny,nz) )             

          call da_get_var_3d_real_cdf( file,"P", p, nx, ny, nz,                     &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PB", pb, nx, ny, nz,                   &
                                file_time_index, debug  ) 
          data_out = (p+pb)*.01

          deallocate (  p )
          deallocate ( pb )             
 
  else if(var == 'Z' ) then

          allocate (  ph(nx,ny,nz+1) )
          allocate ( phb(nx,ny,nz+1) )             

          call da_get_var_3d_real_cdf( file,"PH", ph, nx, ny, nz+1,                 &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PHB", phb, nx, ny, nz+1,               &
                                file_time_index, debug  ) 
          ph = (ph+phb)/9.81
          data_out = 0.5*(ph(:,:,1:nz)+ph(:,:,2:nz+1))

          deallocate (  ph )
          deallocate ( phb )             

  else if(var == 'THETA' ) then

          call da_get_var_3d_real_cdf( file,"T", data_out, nx, ny, nz,              &
                                file_time_index, debug  )
          data_out = data_out + 300.

  else if(var == 'T' ) then

          allocate (        p(nx,ny,nz) )
          allocate (       pb(nx,ny,nz) )             
          allocate ( data_tmp(nx,ny,nz) )             

          call da_get_var_3d_real_cdf( file,"P", p, nx, ny, nz,                     &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PB", pb, nx, ny, nz,                   &
                                file_time_index, debug  ) 
          p = p+pb

          call da_get_var_3d_real_cdf( file,"T", data_tmp, nx, ny, nz,              &
                                file_time_index, debug  )
          data_out = (data_tmp+300.)*(p/p1000mb)**rcp

          deallocate (  p       )
          deallocate ( pb       )             
          deallocate ( data_tmp )

  else if(var == 'TC' ) then

          allocate (        p(nx,ny,nz) )
          allocate (       pb(nx,ny,nz) )             
          allocate ( data_tmp(nx,ny,nz) )             

          call da_get_var_3d_real_cdf( file,"P", p, nx, ny, nz,                     &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PB", pb, nx, ny, nz,                   &
                                file_time_index, debug  ) 
          p = p+pb

          call da_get_var_3d_real_cdf( file,"T", data_tmp, nx, ny, nz,              &
                                file_time_index, debug  )
          data_out = (data_tmp+300.)*(p/p1000mb)**rcp -T0

          deallocate (  p       )
          deallocate ( pb       )             
          deallocate ( data_tmp )

  else if(var == 'TD' ) then

          allocate (        p(nx,ny,nz) )
          allocate (       pb(nx,ny,nz) )             
          allocate (       qv(nx,ny,nz) )             
          allocate ( data_tmp(nx,ny,nz) )             

          call da_get_var_3d_real_cdf( file,"P", p, nx, ny, nz,                     &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PB", pb, nx, ny, nz,                   &
                                file_time_index, debug  ) 
          p = p+pb

          call da_get_var_3d_real_cdf( file,"QVAPOR", qv, nx, ny, nz,               &
                                file_time_index, debug  )

          data_tmp = qv*(p/100.)/(0.622+qv)
          data_tmp = AMAX1(data_tmp,0.001)
          data_out = (243.5*log(data_tmp)-440.8)/(19.48-log(data_tmp))

          deallocate (  p       )
          deallocate ( pb       )             
          deallocate ( qv       )             
          deallocate ( data_tmp )

  else if(var == 'RH' ) then

          allocate (         p(nx,ny,nz) )
          allocate (        pb(nx,ny,nz) )             
          allocate (        qv(nx,ny,nz) )             
          allocate (         t(nx,ny,nz) )             
          allocate (  data_tmp(nx,ny,nz) )             
          allocate ( data_tmp2(nx,ny,nz) )             

          call da_get_var_3d_real_cdf( file,"P", p, nx, ny, nz,                     &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PB", pb, nx, ny, nz,                   &
                                file_time_index, debug  ) 
          p = p+pb

          call da_get_var_3d_real_cdf( file,"T", t, nx, ny, nz,                     &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"QVAPOR", qv, nx, ny, nz,               &
                                file_time_index, debug  )

          t = (t+300.)*(p/p1000mb)**rcp
          data_tmp2 = 10.*0.6112*exp(17.67*(t-T0)/(t-29.65))
          data_tmp  = 0.622*data_tmp2/(0.01 * p -  (1.-0.622)*data_tmp2)
          data_out  = 100.*AMAX1(AMIN1(qv/data_tmp,1.0),0.0)


          deallocate (  p        )
          deallocate ( pb        )             
          deallocate ( qv        )             
          deallocate ( t         )             
          deallocate ( data_tmp  )
          deallocate ( data_tmp2 )

  else 
          call da_get_var_3d_real_cdf( file,var(1:length_var),                      &
                                    data_out, nx,ny,nz,                          &
                                    file_time_index, debug  )
  endif


  end subroutine g_output_3d

!-------------------------------------------------------------------------

  subroutine g_output_2d (file, file_time_index, var, length_var,            &
                          nx, ny, nz, data_out, debug)
  implicit none

  character (len=*), intent(in)             ::   file
  integer, intent(in)                       ::   file_time_index
  character (len=*), intent(in)             ::   var
  integer, intent(in)                       ::   length_var
  integer, intent(in)                       ::   nx, ny, nz           
  real, intent(out), dimension(:,:,:)       ::   data_out
  logical, intent(in)                                   ::   debug
  integer, allocatable, dimension(:,:,:)    ::   data_int
  real,    allocatable, dimension(:,:,:)    ::   u10, v10
  real,    allocatable, dimension(:,:)      ::   psfc,t2m,q2m,mu
  real,    allocatable, dimension(:,:)      ::   xlat, xlon
  real,    allocatable, dimension(:,:,:)    ::   z,ph,phb  
  real,    allocatable, dimension(:,:,:)    ::   p,pb  
  real,    allocatable, dimension(:,:,:)    ::   ts,qv 
  integer                                   ::   map_proj
  real                                      ::   cen_lon, truelat1, truelat2

  if(debug) then
     write(6,*) ' calculations for variable ',var
  end if

       if(var == 'SLP') then

          allocate (   z(nx,ny,nz)   )
          allocate (  ph(nx,ny,nz+1) )
          allocate ( phb(nx,ny,nz+1) )             
          allocate (   p(nx,ny,nz)   )             
          allocate (  pb(nx,ny,nz)   )             
          allocate (  ts(nx,ny,nz)   )             
          allocate (  qv(nx,ny,nz)   )             

          call da_get_var_3d_real_cdf( file,"PH", ph, nx, ny,nz+1,                  &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PHB", phb, nx, ny,nz+1,                &
                                file_time_index, debug  ) 
          ph = (ph+phb)/9.81
          z = 0.5*(ph(:,:,1:nz)+ph(:,:,2:nz+1))

          call da_get_var_3d_real_cdf( file,"P", p, nx, ny,nz,                      &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PB", pb, nx, ny,nz,                    &
                                file_time_index, debug  ) 
          p = p+pb

          call da_get_var_3d_real_cdf( file,"T", ts, nx, ny,nz,                     &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"QVAPOR", qv, nx, ny,nz,                &
                                file_time_index, debug  )

          call compute_seaprs (nx, ny, nz, z, ts, p, qv, data_out, debug)


          deallocate (   z )              
          deallocate (  ph )              
          deallocate ( phb )                         
          deallocate (   p )                       
          deallocate (  pb )                       
          deallocate (  ts )                       
          deallocate (  qv )                       

  else if(var == 'MU' ) then
          allocate ( mu(nx,ny) )
          call da_get_var_2d_real_cdf( file,"MU", mu, nx, ny,                 &
                                file_time_index, debug  )
          data_out(:,:,1) = mu(:,:)
          deallocate ( mu )

  else if(var == 'PSFC' ) then
          allocate ( psfc(nx,ny) )
          call da_get_var_2d_real_cdf( file,"PSFC", psfc, nx, ny,                 &
                                file_time_index, debug  )
          data_out(:,:,1) = psfc(:,:)
          deallocate ( psfc )

  else if(var == 'T2M'  ) then
          allocate ( t2m(nx,ny) )
          call da_get_var_2d_real_cdf( file,"T2", t2m, nx, ny,                  &
                                file_time_index, debug  )
          data_out(:,:,1) = t2m(:,:)
          deallocate ( t2m )

  else if(var == 'Q2M'  ) then
          allocate ( q2m(nx,ny) )
          call da_get_var_2d_real_cdf( file,"Q2", q2m, nx, ny,                  &
                                file_time_index, debug  )
          data_out(:,:,1) = q2m(:,:)
          deallocate ( q2m )

  else if(var == 'U10M' ) then

          call da_get_gl_att_int_cdf ( file, 'MAP_PROJ', map_proj, debug )

          IF ( map_proj == 1  .OR.  map_proj == 2 ) THEN

              allocate ( u10(nx,ny,1) )
              allocate ( v10(nx,ny,1) )
              allocate ( xlat(nx, ny) )             
              allocate ( xlon(nx, ny) )             
              call da_get_var_2d_real_cdf( file,"U10", u10, nx, ny,                 &
                                    file_time_index, debug  )
              call da_get_var_2d_real_cdf( file,"V10", v10, nx, ny,                 &
                                    file_time_index, debug  )
 
              call da_get_gl_att_real_cdf( file, 'STAND_LON', cen_lon, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT1', truelat1, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT2', truelat2, debug )
              call da_get_var_2d_real_cdf( file, 'XLAT', xlat,nx,ny, 1,debug )
              call da_get_var_2d_real_cdf( file, 'XLONG',xlon,nx,ny, 1,debug )

              call rotate_wind (u10,v10,nx,ny,1,var,                             &
                                map_proj,cen_lon,xlat,xlon,                      &
                                truelat1,truelat2,data_out)

              deallocate ( xlat )             
              deallocate ( xlon )             
              deallocate ( u10  )
              deallocate ( v10  )

          ELSE

              call da_get_var_2d_real_cdf( file,"U10", data_out, nx, ny,            &
                                    file_time_index, debug  )

          ENDIF

  else if(var == 'V10M' ) then

          call da_get_gl_att_int_cdf ( file, 'MAP_PROJ', map_proj, debug )

          IF ( map_proj == 1  .OR.  map_proj == 2 ) THEN

              allocate ( u10(nx,ny,1) )
              allocate ( v10(nx,ny,1) )
              allocate ( xlat(nx, ny) )             
              allocate ( xlon(nx, ny) )             
              call da_get_var_2d_real_cdf( file,"U10", u10, nx, ny,                 &
                                    file_time_index, debug  )
              call da_get_var_2d_real_cdf( file,"V10", v10, nx, ny,                 &
                                    file_time_index, debug  )
 
              call da_get_gl_att_real_cdf( file, 'STAND_LON', cen_lon, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT1', truelat1, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT2', truelat2, debug )
              call da_get_var_2d_real_cdf( file, 'XLAT', xlat,nx,ny, 1,debug )
              call da_get_var_2d_real_cdf( file, 'XLONG',xlon,nx,ny, 1,debug )

              call rotate_wind (u10,v10,nx,ny,1,var,                             &
                                map_proj,cen_lon,xlat,xlon,                      &
                                truelat1,truelat2,data_out)

              deallocate ( xlat )             
              deallocate ( xlon )             
              deallocate ( u10  )
              deallocate ( v10  )

          ELSE

              call da_get_var_2d_real_cdf( file,"V10", data_out, nx, ny,            &
                                    file_time_index, debug  )

          ENDIF

  else if(var == 'XLONG' ) then
          call da_get_var_2d_real_cdf( file,var(1:length_var),                      &
                                    data_out, nx,ny,                             &
                                    file_time_index, debug  )
          WHERE ( data_out < 0.0 )
             data_out = data_out + 360.0
          ENDWHERE

  else if(var == 'IVGTYP' .or. var == 'ISLTYP') then

          allocate (data_int(nx,ny,1))
          call da_get_var_2d_int_cdf( file,var(1:length_var),                       &
                                    data_int, nx,ny,                             &
                                    file_time_index, debug  )
          data_out = data_int
          deallocate (data_int)

  else 
          call da_get_var_2d_real_cdf( file,var(1:length_var),                      &
                                    data_out, nx,ny,                             &
                                    file_time_index, debug  )
  endif


  end subroutine g_output_2d

!------------------------------------------------------------------

  subroutine interp_to_z( data_in , nx_in , ny_in , nz_in , &
                              data_out, nx_out, ny_out, nz_out, &
                              z_in, z_out, missing_value, &
                              vertical_type, debug   )
  implicit none
  integer, intent(in)                                  :: nx_in , ny_in , nz_in 
  integer, intent(in)                                  :: nx_out, ny_out, nz_out
  real, intent(in)                                     :: missing_value
  real, dimension(nx_in , ny_in , nz_in ), intent(in ) :: data_in, z_in
  real, dimension(nx_out, ny_out, nz_out), intent(out) :: data_out
  real, dimension(nz_out), intent(in)                  :: z_out
  logical, intent(in)                                  :: debug
  character (len=1)                , intent(in)                     :: vertical_type

  real, dimension(nz_in)                               :: data_in_z, zz_in
  real, dimension(nz_out)                              :: data_out_z

  integer :: i,j,k

    do i=1,nx_in
    do j=1,ny_in

      do k=1,nz_in
        data_in_z(k) = data_in(i,j,k)
        zz_in(k) = z_in(i,j,k)
      enddo

!Hui      do k=1,nz_out
!Hui        data_out_z(k) = data_out(i,j,k)
!Hui      enddo

      call interp_1d( data_in_z, zz_in, nz_in, &
                      data_out_z, z_out, nz_out, &
                      vertical_type, missing_value )

      do k=1,nz_out
        data_out(i,j,k) = data_out_z(k)
      enddo


    enddo
    enddo

  end subroutine interp_to_z

!----------------------------------------------

  subroutine interp_1d( a, xa, na, &
                        b, xb, nb, vertical_type, missing_value )
  implicit none
  integer, intent(in)              ::  na, nb
  real, intent(in), dimension(na)  :: a, xa
  real, intent(in), dimension(nb)  :: xb
  real, intent(out), dimension(nb) :: b
  real, intent(in)                 :: missing_value

  integer                          :: n_in, n_out
  logical                          :: interp
  real                             :: w1, w2
  character (len=1) ,intent(in)               :: vertical_type


  if ( vertical_type == 'p' ) then

  do n_out = 1, nb

    b(n_out) = missing_value
    interp = .false.
    n_in = 1

    do while ( (.not.interp) .and. (n_in < na) )

      if( (xa(n_in)   >= xb(n_out)) .and. &
          (xa(n_in+1) <= xb(n_out))        ) then
        interp = .true.
        w1 = (xa(n_in+1)-xb(n_out))/(xa(n_in+1)-xa(n_in))
        w2 = 1. - w1
        b(n_out) = w1*a(n_in) + w2*a(n_in+1)
      end if
      n_in = n_in +1

    enddo

  enddo
  
  else

  do n_out = 1, nb

    b(n_out) = missing_value
    interp = .false.
    n_in = 1

    do while ( (.not.interp) .and. (n_in < na) )

      if( (xa(n_in)   <= xb(n_out)) .and. &
          (xa(n_in+1) >= xb(n_out))        ) then
        interp = .true.
        w1 = (xa(n_in+1)-xb(n_out))/(xa(n_in+1)-xa(n_in))
        w2 = 1. - w1
        b(n_out) = w1*a(n_in) + w2*a(n_in+1)
      end if
      n_in = n_in +1

    enddo

  enddo
  
  endif

  end subroutine interp_1d

!-------------------------------------------------------------------------
!
! This routines has been taken "as is" from wrf_user_fortran_util_0.f
!
! This routine assumes
!    index order is (i,j,k)
!    wrf staggering
!    units: pressure (Pa), temperature(K), height (m), mixing ratio (kg kg{-1})
!    availability of 3d p, t, and qv; 2d terrain; 1d half-level zeta string
!    output units of SLP are Pa, but you should divide that by 100 for the
!          weather weenies.
!    virtual effects are included
!
! Dave

  subroutine compute_seaprs ( nx , ny , nz  ,         &
                                  z, t , p , q ,          &
                                  sea_level_pressure,debug)
!     &                            t_sea_level, t_surf, level )
      IMPLICIT NONE
!     Estimate sea level pressure.
      INTEGER, intent(in) :: nx , ny , nz
      REAL, intent(in) ::    z(nx,ny,nz)
      REAL, intent(in)  ::  p(nx,ny,nz) , q(nx,ny,nz)
      REAL, intent(inout)  ::  t(nx,ny,nz)
!     The output is the 2d sea level pressure.
      REAL, intent(out)   :: sea_level_pressure(nx,ny)
      INTEGER level(nx,ny)
      REAL t_surf(nx,ny) , t_sea_level(nx,ny)
      LOGICAL, intent(in) :: debug

!     Some required physical constants:

      REAL R, G, GAMMA
      PARAMETER (R=287.04, G=9.81, GAMMA=0.0065)

!     Specific constants for assumptions made in this routine:

      REAL    TC, PCONST
      PARAMETER (TC=273.16+17.5, PCONST = 10000)
      LOGICAL ridiculous_mm5_test
      PARAMETER (ridiculous_mm5_test = .TRUE.)
!      PARAMETER (ridiculous_mm5_test = .false.)

!     Local variables:

      INTEGER i , j , k
      INTEGER klo , khi


      REAL plo , phi , tlo, thi , zlo , zhi
      REAL p_at_pconst , t_at_pconst , z_at_pconst
      REAL z_half_lowest

      REAL    , PARAMETER :: cp           = 7.*R/2.
      REAL    , PARAMETER :: rcp          = R/cp
      REAL    , PARAMETER :: p1000mb      = 100000.

      LOGICAL  l1 , l2 , l3, found

!     Find least zeta level that is PCONST Pa above the surface.  We later use this
!     level to extrapolate a surface pressure and temperature, which is supposed
!     to reduce the effect of the diurnal heating cycle in the pressure field.

      t(:,:,:) = (t(:,:,:)+300.)*(p(:,:,:)/p1000mb)**rcp

      DO j = 1 , ny
         DO i = 1 , nx
            level(i,j) = -1

            k = 1
            found = .false.
            do while( (.not. found) .and. (k.le.nz))
               IF ( p(i,j,k) .LT. p(i,j,1)-PCONST ) THEN
                  level(i,j) = k
                  found = .true.
               END IF
               k = k+1
            END DO

            IF ( level(i,j) .EQ. -1 ) THEN
            PRINT '(A,I4,A)','Troubles finding level ',   &
                        NINT(PCONST)/100,' above ground.'
            PRINT '(A,I4,A,I4,A)',                        &
                  'Problems first occur at (',i,',',j,')'
            PRINT '(A,F6.1,A)',                           &
                  'Surface pressure = ',p(i,j,1)/100,' hPa.'
            STOP 'Error_in_finding_100_hPa_up'
         END IF


         END DO
      END DO

!     Get temperature PCONST Pa above surface.  Use this to extrapolate
!     the temperature at the surface and down to sea level.

      DO j = 1 , ny
         DO i = 1 , nx

            klo = MAX ( level(i,j) - 1 , 1      )
            khi = MIN ( klo + 1        , nz - 1 )

            IF ( klo .EQ. khi ) THEN
               PRINT '(A)','Trapping levels are weird.'
               PRINT '(A,I3,A,I3,A)','klo = ',klo,', khi = ',khi, &
                            ': and they should not be equal.'
               STOP 'Error_trapping_levels'
            END IF

         plo = p(i,j,klo)
         phi = p(i,j,khi)
         tlo = t(i,j,klo)*(1. + 0.608 * q(i,j,klo) )
         thi = t(i,j,khi)*(1. + 0.608 * q(i,j,khi) )
!         zlo = zetahalf(klo)/ztop*(ztop-terrain(i,j))+terrain(i,j)
!         zhi = zetahalf(khi)/ztop*(ztop-terrain(i,j))+terrain(i,j)
         zlo = z(i,j,klo)
         zhi = z(i,j,khi)

         p_at_pconst = p(i,j,1) - pconst
         t_at_pconst = thi-(thi-tlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)
         z_at_pconst = zhi-(zhi-zlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)

         t_surf(i,j) = t_at_pconst*(p(i,j,1)/p_at_pconst)**(gamma*R/g)
         t_sea_level(i,j) = t_at_pconst+gamma*z_at_pconst

         END DO
      END DO

!     If we follow a traditional computation, there is a correction to the sea level
!     temperature if both the surface and sea level temnperatures are *too* hot.

      IF ( ridiculous_mm5_test ) THEN
         DO j = 1 , ny
            DO i = 1 , nx
               l1 = t_sea_level(i,j) .LT. TC
               l2 = t_surf     (i,j) .LE. TC
               l3 = .NOT. l1
               IF ( l2 .AND. l3 ) THEN
                  t_sea_level(i,j) = TC
               ELSE
                  t_sea_level(i,j) = TC - 0.005*(t_surf(i,j)-TC)**2
               END IF
            END DO
         END DO
      END IF

!     The grand finale: ta da!

      DO j = 1 , ny
      DO i = 1 , nx
!         z_half_lowest=zetahalf(1)/ztop*(ztop-terrain(i,j))+terrain(i,j)
         z_half_lowest=z(i,j,1)
         sea_level_pressure(i,j) = p(i,j,1) *              &
                               EXP((2.*g*z_half_lowest)/   &
                                   (R*(t_sea_level(i,j)+t_surf(i,j))))
      END DO
      END DO

    if (debug) then
      print *,'sea pres input at weird location i=20,j=1,k=1'
      print *,'t=',t(20,1,1),t(20,2,1),t(20,3,1)
      print *,'z=',z(20,1,1),z(20,2,1),z(20,3,1)
      print *,'p=',p(20,1,1),p(20,2,1),p(20,3,1)
      print *,'slp=',sea_level_pressure(20,1),     &
               sea_level_pressure(20,2),sea_level_pressure(20,3)
    endif
!      print *,'t=',t(10:15,10:15,1),t(10:15,2,1),t(10:15,3,1)
!      print *,'z=',z(10:15,1,1),z(10:15,2,1),z(10:15,3,1)
!      print *,'p=',p(10:15,1,1),p(10:15,2,1),p(10:15,3,1)
!      print *,'slp=',sea_level_pressure(10:15,10:15),     &
!         sea_level_pressure(10:15,10:15),sea_level_pressure(20,10:15)

  end subroutine compute_seaprs

!------------------------------------------------------------------


  subroutine rotate_wind (u,v,d1,d2,d3,var,                 &
                          map_proj,cen_lon,xlat,xlon,       &
                          truelat1,truelat2,data_out)

  implicit none

  integer, intent(in)            ::  d1, d2, d3

  real, dimension(d1,d2,d3), intent(out)      :: data_out
  integer, intent(in)                        :: map_proj
  integer                        ::i,j,k
  real, intent(in)                           :: cen_lon, truelat1, truelat2
  real                          :: cone
  real, dimension(d1,d2,d3), intent(in)      :: u,v 
  real, dimension(d1,d2), intent(in)         :: xlat, xlon
  real, dimension(d1,d2)         :: diff, alpha

  character (len=10), intent(in) :: var

   REAL    , PARAMETER           :: pii = 3.14159265
   REAL    , PARAMETER           :: radians_per_degree = pii/180.




       cone = 1.                                          !  PS
       if( map_proj .eq. 1) then                          !  Lambert Conformal mapping
         IF (ABS(truelat1-truelat2) .GT. 0.1) THEN
            cone=(ALOG(COS(truelat1*radians_per_degree))-            &
                  ALOG(COS(truelat2*radians_per_degree))) /          &
            (ALOG(TAN((90.-ABS(truelat1))*radians_per_degree*0.5 ))- &
             ALOG(TAN((90.-ABS(truelat2))*radians_per_degree*0.5 )) )
         ELSE
            cone = SIN(ABS(truelat1)*radians_per_degree )  
         ENDIF
       end if


       diff = xlon - cen_lon

       do i = 1, d1
       do j = 1, d2
         if(diff(i,j) .gt. 180.) then
           diff(i,j) = diff(i,j) - 360.
         end if
         if(diff(i,j) .lt. -180.) then
           diff(i,j) = diff(i,j) + 360.
         end if
       end do
       end do

       do i = 1, d1
       do j = 1, d2
          if(xlat(i,j) .lt. 0.) then
            alpha(i,j) = - diff(i,j) * cone * radians_per_degree
          else
            alpha(i,j) = diff(i,j) * cone * radians_per_degree
          end if
       end do
       end do


       if(var(1:1) .eq. "U") then
         do k=1,d3
           data_out(:,:,k) = v(:,:,k)*sin(alpha) + u(:,:,k)*cos(alpha)
         end do
       else if(var(1:1) .eq. "V") then    
         do k=1,d3
           data_out(:,:,k) = v(:,:,k)*cos(alpha) - u(:,:,k)*sin(alpha)
         end do
       end if


  end subroutine rotate_wind


!------------------------------------------------------------------
  subroutine handle_err(rmarker,nf_status)

#include "netcdf.inc"
      integer, intent(in) :: nf_status
      character*(*), intent(in)        :: rmarker
      if (nf_status .ne. nf_noerr) then
         write(*,*)  'NetCDF error : ',rmarker
         write(*,*)  '  ',nf_strerror(nf_status)
         stop 
      endif
      
  end subroutine handle_err

!------------------------------------------------------------------
  subroutine get_dimensions(infile,nx,ny,nz)
    
#include "netcdf.inc"
     character (len=512), intent(in)                          :: infile
     integer                                                  :: ncid, dimid, nf_status
     integer, intent(inout)                                   :: nx, ny, nz
     integer                                                  :: nlgen

!  need to pull out some data to set up dimensions, etc.
      nf_status = nf_open (infile, nf_nowrite, ncid)
      call handle_err('Error opening file: '//trim(infile),nf_status)
!
        nf_status = nf_inq_dimid (ncid, 'west_east_stag', dimid)
        call handle_err('west_east_stag',nf_status)
        nf_status = nf_inq_dimlen (ncid, dimid, nx)
        call handle_err('Get NX',nf_status)
        nx = nx-1
!
        nf_status = nf_inq_dimid (ncid, 'south_north_stag', dimid)
        call handle_err('south_north_stag',nf_status)
        nf_status = nf_inq_dimlen (ncid, dimid, ny)
        call handle_err('Get NY',nf_status)
        ny = ny-1
!
        nf_status = nf_inq_dimid (ncid, 'bottom_top', dimid)
        call handle_err('bottom_top',nf_status)
        nf_status = nf_inq_dimlen (ncid, dimid, nz)
        call handle_err('Get NZ',nf_status)
        nlgen = nz
!
      nf_status = nf_close (ncid)
      call handle_err('Error closing file',nf_status)

  end subroutine get_dimensions
!------------------------------------------------------------------
  
  subroutine get_diffs(var1, var2, diff, absdiff, sqdiff, nx, ny, nz, missing)
    
    real, intent(in), dimension(:,:,:)             :: var1, var2
    real, intent(out), dimension(:,:,:)            :: diff, absdiff, sqdiff
    integer, intent(in)                            :: nx, ny, nz
    integer                                        :: i,j,k
    real, intent(in)                               :: missing

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
       if ( var1(i,j,k) /= missing .and. var2(i,j,k) /= missing ) then
         diff(i,j,k)   = var2(i,j,k) - var1(i,j,k) 
         absdiff(i,j,k)   = abs(var2(i,j,k) - var1(i,j,k)) 
         sqdiff(i,j,k) = (var2(i,j,k) - var1(i,j,k) )*(var2(i,j,k) - var1(i,j,k) ) 
       else
         diff(i,j,k)   = missing 
         absdiff(i,j,k)   = missing 
         sqdiff(i,j,k) = missing
       endif
    enddo
    enddo
    enddo

  end subroutine get_diffs 

!------------------------------------------------------------------
  subroutine get_diffs_wind(u1, v1, u2,v2,diff, absdiff, sqdiff, nx, ny, nz, missing)
    
    real, intent(in), dimension(:,:,:)             :: u1, u2, v1, v2
    real, intent(out), dimension(:,:,:)            :: diff, absdiff, sqdiff
    integer, intent(in)                            :: nx, ny, nz
    integer                                        :: i,j,k
    real, intent(in)                               :: missing
    real                                           :: u, v

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
       if ( u1(i,j,k) /= missing .and. v1(i,j,k) /= missing .and. &
            u2(i,j,k) /= missing .and. v2(i,j,k) /= missing ) then
         u  = u1(i,j,k) - u2(i,j,k) 
         v  = v1(i,j,k) - v2(i,j,k) 
         diff(i,j,k)   = u + v
         absdiff(i,j,k)   = abs(u) + abs (v)                
         sqdiff(i,j,k) = u*u + v*v
       else
         diff(i,j,k)   = missing 
         absdiff(i,j,k)   = missing 
         sqdiff(i,j,k) = missing
       endif
    enddo
    enddo
    enddo

  end subroutine get_diffs_wind 

!------------------------------------------------------------------
  subroutine domain_average(var, avg_prof, counter, nx, ny, nz, missing,isq)

  integer, intent(in)                                 :: nx,ny,nz,isq
  real, intent(in), dimension(:,:,:)                  :: var
  integer, intent(out), dimension(:)                  :: counter
  real, intent(out), dimension(:)                     :: avg_prof
  real, intent(in)                                    :: missing

  integer                                             :: i,j,k,icount
  real                                                :: dsum, dmiss
  integer                                             :: imiss

!9999
!Hui: set dmiss value consistent with plot script
!  dmiss = -9999.999
  dmiss = -99.99
  imiss = -99
  do k = 1, nz
     icount = 0
     dsum = 0
     do j = 1, ny
     do i = 1, nx
         if ( var(i,j,k) /= missing ) then
            icount = icount + 1
            dsum = dsum + var(i,j,k)
         endif
     enddo
     enddo   
     avg_prof(k) = dmiss
!Hui     counter(k) = 0
     counter(k) = imiss
     if (icount /= 0 ) then
        counter(k) = icount
        if ( isq .eq. 0 ) then
          avg_prof(k) = dsum /float(icount)
        else
          avg_prof(k) = sqrt(dsum /float(icount))
        endif
      endif
  enddo

  end subroutine domain_average
!------------------------------------------------------------------
  
  subroutine mask_domain_average(var, avg_prof, counter, nx, ny, nz, missing,isq)

  integer, intent(in)                                 :: nx,ny,nz,isq
  real, intent(in), dimension(nx,ny,nz)               :: var
  integer, intent(out), dimension(3,nz)               :: counter
  real, intent(out), dimension(3,nz)                  :: avg_prof
  real, intent(in)                                    :: missing

  integer                                             :: i,j,k,imask
  integer                                             :: icount(3)
  real                                                :: dsum(3)
  real                                                :: dmiss
  integer                                             :: imiss

  ! 1: all, 2:land, 3:water

  dmiss = -99.99
  imiss = -99
  do k = 1, nz
     icount = 0
     dsum = 0.0
     do j = 1, ny
        do i = 1, nx
           if ( var(i,j,k) /= missing ) then
              if ( i>=istart .and. i<=iend .and. &
                   j>=jstart .and. j<=jend ) then
                  icount(1) = icount(1) + 1
                  dsum(1) = dsum(1) + var(i,j,k)
                  if ( int(landmask(i,j)) == island ) then
                     icount(2) = icount(2) + 1
                     dsum(2) = dsum(2) + var(i,j,k)
                  else if ( int(landmask(i,j)) == iswater ) then
                     icount(3) = icount(3) + 1
                     dsum(3) = dsum(3) + var(i,j,k)
                  end if
              end if
           end if
        end do
     end do
     avg_prof(:,k) = dmiss
     counter(:,k) = imiss
     do imask = 1, 3
        if ( icount(imask) /= 0 ) then
           counter(imask,k) = icount(imask)
           if ( isq .eq. 0 ) then
             avg_prof(imask,k) = dsum(imask)/float(icount(imask))
           else
             avg_prof(imask,k) = sqrt(dsum(imask)/float(icount(imask)))
           end if
        end if
     end do
  enddo

  end subroutine mask_domain_average
!------------------------------------------------------------------

  subroutine get_sum(dsum, dvar, nx, ny, nz, missing)
    
    integer, intent(in)                               :: nx, ny, nz
    real, intent(in)                                  :: missing
    real, intent(in),dimension(:,:,:)                 :: dvar
    real, intent(inout),dimension(:,:,:)              :: dsum

    integer                                           :: i,j,k

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
       if ( dvar(i,j,k) /= missing .and. dsum(i,j,k) /= missing ) then
         dsum(i,j,k)   = dsum(i,j,k) + dvar(i,j,k) 
       else
         dsum(i,j,k) = missing
       endif
    enddo
    enddo
    enddo

  end subroutine get_sum
!------------------------------------------------------------------
  subroutine time_average(dvar, nx, ny, nz, time_count,missing, isqr)
    
    integer, intent(in)                               :: nx, ny, nz,time_count,isqr
    real, intent(in)                                  :: missing
    real, intent(inout), dimension(:,:,:)             :: dvar 

    integer                                           :: i,j,k

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
       if ( dvar(i,j,k) /= missing ) then
         if (isqr .eq. 1 ) then
           dvar(i,j,k) = sqrt(dvar(i,j,k)/float(time_count))
         else
           dvar(i,j,k)   =  dvar(i,j,k)/float(time_count)
         endif
       else
         dvar(i,j,k) = missing
       endif
    enddo
    enddo
    enddo

  end subroutine time_average
!------------------------------------------------------------------
  subroutine compute_data_3d( infile, var, length, nx, ny, nz, levels, time_idx,   &
                 vert_levels, vertical_type, missing, data_out_z, debug )

  integer, intent(in)                      :: time_idx
  integer, intent(in)                      :: nx, ny, nz, levels
  integer, intent(in)                      :: length
  real, intent(in)                         :: missing
  real, intent(in)                         :: vert_levels(:)   
  character (len=1), intent(in)            :: vertical_type
  character (len=*), intent(in)            :: var
  character (len=*), intent(in)            :: infile
  logical, intent(in)                      :: debug

  real, intent(out), dimension(:,:,:)      :: data_out_z
  real, allocatable, dimension(:,:,:)      :: data_out
  real, allocatable, dimension(:,:,:)      :: z, ph, phb
  real, allocatable, dimension(:,:,:)      :: p, pb


  !  first, get some base-state-stuff

      if ( vertical_type == 'p' ) then
        allocate( p (nx, ny, nz)  )
        allocate( pb(nx, ny, nz)  )
        call da_get_var_3d_real_cdf( infile,'PB',pb,nx,ny,nz,time_idx,debug )
        call da_get_var_3d_real_cdf( infile,'P',p, nx, ny, nz, time_idx,debug )
        p = (p+pb)/100.0 ! Convert to hPa
      endif
      if ( vertical_type == 'z' ) then
        allocate( z (nx, ny, nz)  )
        allocate( ph (nx, ny, nz+1)  )
        allocate( phb(nx, ny, nz+1)  )
        call da_get_var_3d_real_cdf( infile,'PHB',phb,nx,ny,nz+1,time_idx,debug )
        call da_get_var_3d_real_cdf( infile,'PH',ph, nx, ny, nz+1, time_idx,debug )
        ph = (ph+phb)/9.81
        z = 0.5*(ph(:,:,1:nz)+ph(:,:,2:nz+1))
        z = z/1000.  ! Convert to kilometers 
      endif

      allocate ( data_out (nx, ny, nz) )

      call g_output_3d (infile, time_idx, var, length, nx, ny, nz, data_out, debug)

      if ( vertical_type == 'p' ) then
         call interp_to_z( data_out, nx, ny, nz, data_out_z, nx, ny, levels,     &
                         p, vert_levels, missing, vertical_type, debug )

      else if ( vertical_type == 'z' ) then
         call interp_to_z( data_out, nx, ny, nz, data_out_z, nx, ny, levels,     &
                         z, vert_levels, missing, vertical_type, debug )
      else
        data_out_z = data_out
      endif
      deallocate ( data_out )
      if ( vertical_type == 'p' ) then
              deallocate( p )
              deallocate( pb )
      endif
      if ( vertical_type == 'z' ) then
              deallocate( z )
              deallocate( ph )
              deallocate( phb )
      endif

  end subroutine compute_data_3d
  
  subroutine compute_wind_3d( infile, nx, ny, nz, levels, time_idx,   &
                 vert_levels, vertical_type, missing, data_out_z1, data_out_z2, debug )

  integer, intent(in)                      :: time_idx
  integer, intent(in)                      :: nx, ny, nz, levels
  real, intent(in)                         :: missing
  real, intent(in)                         :: vert_levels(:)   
  character (len=1), intent(in)            :: vertical_type
  character (len=*), intent(in)            :: infile
  logical, intent(in)                      :: debug

  real, intent(out), dimension(:,:,:)      :: data_out_z1
  real, intent(out), dimension(:,:,:)      :: data_out_z2
  real, allocatable, dimension(:,:,:)      :: data_out
  real, allocatable, dimension(:,:,:)      :: z, ph, phb
  real, allocatable, dimension(:,:,:)      :: p, pb
  character (len=1)                        :: var


  !  first, get some base-state-stuff

      if ( vertical_type == 'p' ) then
        allocate( p (nx, ny, nz)  )
        allocate( pb(nx, ny, nz)  )
        call da_get_var_3d_real_cdf( infile,'PB',pb,nx,ny,nz,time_idx,debug )
        call da_get_var_3d_real_cdf( infile,'P',p, nx, ny, nz, time_idx,debug )
        p = (p+pb)/100.0 ! Convert to hPa
      endif
      if ( vertical_type == 'z' ) then
        allocate( z (nx, ny, nz)  )
        allocate( ph (nx, ny, nz+1)  )
        allocate( phb(nx, ny, nz+1)  )
        call da_get_var_3d_real_cdf( infile,'PHB',phb,nx,ny,nz+1,time_idx,debug )
        call da_get_var_3d_real_cdf( infile,'PH',ph, nx, ny, nz+1, time_idx,debug )
        ph = (ph+phb)/9.81
        z = 0.5*(ph(:,:,1:nz)+ph(:,:,2:nz+1))
        z = z/1000.  ! Convert to kilometers 
      endif

      allocate ( data_out2 (nx, ny, nz) )
      allocate ( data_out1 (nx, ny, nz) )
      var='U'
      call g_output_3d (infile, time_idx, var, 1, nx, ny, nz, data_out1, debug)
      var='V'
      call g_output_3d (infile, time_idx, var, 1, nx, ny, nz, data_out2, debug)

      if ( vertical_type == 'p' ) then
         call interp_to_z( data_out1, nx, ny, nz, data_out_z1, nx, ny, levels,     &
                         p, vert_levels, missing, vertical_type, debug )
         call interp_to_z( data_out2, nx, ny, nz, data_out_z2, nx, ny, levels,     &
                         p, vert_levels, missing, vertical_type, debug )

      else if ( vertical_type == 'z' ) then
         call interp_to_z( data_out1, nx, ny, nz, data_out_z1, nx, ny, levels,     &
                         z, vert_levels, missing, vertical_type, debug )
         call interp_to_z( data_out2, nx, ny, nz, data_out_z2, nx, ny, levels,     &
                         z, vert_levels, missing, vertical_type, debug )
      else
        data_out_z1 = data_out1
        data_out_z2 = data_out2
      endif
      deallocate ( data_out1 )
      deallocate ( data_out2 )
      if ( vertical_type == 'p' ) then
              deallocate( p )
              deallocate( pb )
      endif
      if ( vertical_type == 'z' ) then
              deallocate( z )
              deallocate( ph )
              deallocate( phb )
      endif

  end subroutine compute_wind_3d
  
!---------------------------------------------------------------------
end program da_verif_grid
