% script create_roms_init
%
% Create a netcdf file that contains initialization data for ROMS.
% Initializes temp, salt, u, v, ubar, vbar, and 
% all sediment parameters.
%
% In cppdefs.h you should have 
% #undef ana_initial
% #undef ana_sediment
%
% Users can adapt this file to their own application.
%
% jcw 7-8-2008
%

%!         W-level  RHO-level                                           !
%!                                                                      !
%!            N     _________                                           !
%!                 |         |                                          !
%!                 |    N    |                                          !
%!          N-1    |_________|  S                                       !
%!                 |         |  E                                       !
%!                 |   N-1   |  A                                       !
%!            2    |_________|  W                                       !
%!                 |         |  A                                       !
%!                 |    2    |  T                                       !
%!            1    |_________|  E                                       !
%!                 |         |  R                                       !
%!                 |    1    |                                          !
%!            0    |_________|_____ bathymetry                          !
%!                 |/////////|                                          !
%!                 |    1    |                                          !
%!            1    |_________|  S                                       !
%!                 |         |  E                                       !
%!                 |    2    |  D                                       !
%!            2    |_________|  I                                       !
%!                 |         |  M                                       !
%!                 |  Nbed-1 |  E                                       !
%!        Nbed-1   |_________|  N                                       !
%!                 |         |  T                                       !
%!                 |  Nbed   |                                          !
%!         Nbed    |_________|                                          !
%!                                                                      !

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Begin user input section                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1) Enter name of netcdf initial file to be created.
%   If it already exists it will be overwritten!!.
     init_file='wetdry_dam_break_init.nc';

%2) Enter start time of initial file, in seconds.
%   This time needs to be consistent with model time (ie dstart and time_ref).
%   See *.in files for more detail. 
    init_time=0*24*3600;

%3) Enter number of vertical sigma levels in model.
%   This will be same value as entered in mod_param.F
    N=8;

%4) Enter the values of theta_s, theta_b, and Tcline from your *.in file.
    theta_s = 0.0;  
    theta_b = 0.0;  
    Tcline =  0.0;  

%5) Enter value of h, Lm, and Mm.
%   This info can come from a grid file or user supplied here.
%   
%   Are you entering a grid file name (1 = yes, 0 = no)? 
    get_grid = 1;    %<--- put a 1 or 0 here
  
    if (get_grid)
       grid_file='wetdry_dam_break_grd.nc'    %<-enter name of grid here
%
% Get some grid info, do not change this.
% 
      nc=netcdf(grid_file);
      h=nc{'h'}(:);
      [MP,LP]=size(h);
      close(nc);
%
    else
      Lm=100;       %<--- else put size of grid here, from mod_param.F
      Mm=20;        %<--- else put size of grid here, from mod_param.F
      LP = Lm+2;    %don't change this.
      MP = Mm+2;    %don't change this.

      % enter depth, same as in ana_grid
      for j=1:MP
        for i=1:LP
          h(j,i)=18-16*(Mm-j)/(Mm-1);
        end
      end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc some grid stuff here - do not change this.
% You go on to step 6.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   L  = LP-1;
   M  = MP-1;
   xi_psi  = L;
   xi_rho  = LP;
   xi_u    = L;
   xi_v    = LP;
   eta_psi = M;
   eta_rho = MP;
   eta_u   = MP;
   eta_v   = M;
%
% Don't change this either.  This is from set_scoord.F
% Go to step 6 now.
%
   hmin=0;
   hc=min([hmin,Tcline]);
   if (theta_s~=0.0)
     cff1=1.0/sinh(theta_s);
     cff2=0.5/tanh(0.5*theta_s);
   end
   sc_w(1)=-1.0;
   Cs_w(1)=-1.0;
   cff=1.0/N;
   for k=1:N
     sc_w(k+1)=cff*(k-N);
     sc_r(k)=cff*((k-N)-0.5);
     if (theta_s~=0)
       Cs_w(k+1)=(1.0-theta_b)*cff1*sinh(theta_s*sc_w(k+1))+   ...
                      theta_b*(cff2*tanh(theta_s*(sc_w(k+1)+0.5))-0.5);
       Cs_r(k)  =(1.0-theta_b)*cff1*sinh(theta_s*sc_r(k))+   ...
                      theta_b*(cff2*tanh(theta_s*(sc_r(k)+0.5))-0.5);
     else
       Cs_w(k+1)=sc_w(k+1);
       Cs_r(k)=sc_r(k);
      end
    end

%6) Initialize zeta, salt, temp, u, v, ubar, vbar.
%
% Init values for zeta.
  display('Initializing zeta')
%
    for j=1:eta_rho
      for i=1:xi_rho
        for time=1:length(init_time)
          zeta(time,j,i) = 0.6;
          if (i>26)
            zeta(time,j,i) = 0.0;
          end
        end
      end
    end
%
%  Calc z at rho and w points.
%  Don't change this.
%
    for  k = 1:N+1
      z_w(k,:,:) = squeeze(zeta(1,:,:)).*(1+sc_w(k)*hc./h-hc*Cs_w(k)./h+Cs_w(k))+hc*sc_w(k)+(h - hc)*Cs_w(k);
    end
    for  k = 1:N
      z_r(k,:,:) = squeeze(zeta(1,:,:)).*(1+sc_r(k)*hc./h-hc*Cs_r(k)./h+Cs_r(k))+hc*sc_r(k)+(h - hc)*Cs_r(k);
      Hz(k,:,:) = z_w(k+1,:,:)-z_w(k,:,:);
    end
%
% Init values for u, ubar, v, and vbar.
  display('Initializing u, v, ubar, and vbar')
%
  for j=1:eta_u
    for i=1:xi_u
      for time=1:length(init_time)
        for k=1:N 
          u(time,k,j,i) = 0;
        end
        ubar(time,j,i) = 0;
      end
    end
  end
%
  for j=1:eta_v
    for i=1:xi_v
      for time=1:length(init_time)
        for k=1:N
          v(time,k,j,i)  = 0;
        end
        vbar(time,j,i) = 0;
      end
    end
  end
%
% Init values for temp and salt.
  display('Initializing temp and salt')
%
    NAT=2;                            % Number of active tracers. Usually 2 (temp + salt). Same 
                                      % value as used in mod_param.F
    for j=1:eta_rho
      for i=1:xi_rho
        for time=1:length(init_time)
          for k=1:N
            salt(time,k,j,i) = 30;
            temp(time,k,j,i) = 10;
          end
        end
      end
    end


%7) Enter number of mud sediments (NCS) and number of sand sediments (NNS).
%   These values should be the same as in mod_param.F
    NCS = 0;   %number of cohesive sed classes
    NNS = 1;   %number of non-cohesive sed classes
%
% calc sed parameters. Do not alter.
%
   NST = NCS + NNS;     % total number of sed tracers.
   NT = NAT+NST;        % total number of tracers.

%8) Enter number of bed layers
%   This value should be the same as in mod_param.F
    Nbed = 1;

%9) Sediment class properties (in order, mud first then sand).
%  These values should coincide with your sediment.in file.
  mud_Srho=ones(1,NCS)*2650;        %kg m-3, NCS values
  mud_Sd50=[0.06]/1000;             %m,      NCS values
  mud_Wsed=[1.0]/1000;              %m s-1,  NCS values
  mud_tau_ce=[0.05];                %N m-2,  NCS values
  mud_Erate=[5]*1e-5;               %kg m-2 s-1, NCS values
  sand_Srho=ones(1,NNS)*2650;       %kg m-3, NNS values
  sand_Sd50=[0.2]/1000;             %m,      NNS values
  sand_Wsed=[27]/1000;              %m s-1,  NNS values
  sand_tau_ce=[0.19];               %N m-2,  NNS values
  sand_Erate=[1]*1e-5;              %kg m-2 s-1, NNS values
%
% make some combined arrays.  Do not alter.
%
  Srho=  [mud_Srho,sand_Srho];
  Sd50=  [mud_Sd50,sand_Sd50];
  Wsed=  [mud_Wsed,sand_Wsed];
  tau_ce=[mud_tau_ce,sand_tau_ce];
  Erate= [mud_Erate,sand_Erate];


%9) Provide initial sediment properties in water column.
    display('Initializing suspended sediments.')
%
% mud.
%
for idmud=1:NCS
  count=['0',num2str(idmud)];
  count=count(end-1:end);
  for k=1:N
    for j=1:eta_rho
      for i=1:xi_rho
        for time=1:length(init_time)
          eval(['mud_',count,'(time,k,j,i) = 0;'])               %mud conc in water column
        end
      end
    end
  end
end
%
% sand.
%
for isand=1:NNS
  count=['0',num2str(isand)];
  count=count(end-1:end);
  for k=1:N
    for j=1:eta_rho
      for i=1:xi_rho
        for time=1:length(init_time)
          eval(['sand_',count,'(time,k,j,i) = 0;'])               %sand conc in water column
        end
      end
    end
  end
end

%10) Provide initial sediment properties in bed.
%
% bed properties
  display('Initializing sediment bed.')
%
for k=1:Nbed
  for j=1:eta_rho
    for i=1:xi_rho
      for time=1:length(init_time)
        bed_thickness(time,k,j,i) = 1.0;
        bed_age(time,k,j,i)       = init_time(1);
        bed_porosity(time,k,j,i)  = 0.5;
        bed_biodiff(time,k,j,i)   = 0.0;
      end
    end
  end
end
%
% for mud
%
for idsed=1:NCS
  count=['0',num2str(idsed)];
  count=count(end-1:end);
  for k=1:Nbed
    for j=1:eta_rho
      for i=1:xi_rho
        for time=1:length(init_time)
          eval(['mudfrac_',count,'(time,k,j,i) = 1/NST;'])      %fraction of each sed class in each bed cell
          eval(['mudmass_',count,'(time,k,j,i) = bed_thickness(time,k,j,i)*Srho(idsed)*(1.0-bed_porosity(time,k,j,i))*mudfrac_',count,'(time,k,j,i);'])          %mass of each sed class in each bed cell
        end
      end
    end
  end
end
%
% for sand
%
for idsed=1:NNS
 count=['0',num2str(idsed)];
 count=count(end-1:end);
  for k=1:Nbed
    for j=1:eta_rho
      for i=1:xi_rho
        for time=1:length(init_time)
          eval(['sandfrac_',count,'(time,k,j,i) = 1/NST;'])      %fraction of each sed class in each bed cell
          eval(['sandmass_',count,'(time,k,j,i) = bed_thickness(time,k,j,i)*Srho(idsed)*(1.0-bed_porosity(time,k,j,i))*sandfrac_',count,'(time,k,j,i);'])          %mass of each sed class in each bed cell
        end
      end
    end
  end
end

%11)
%
% set some surface properties
  display('Initializing sediment surface properties.')
%
for j=1:eta_rho
  for i=1:xi_rho
    cff1=1.0;
    cff2=1.0;
    cff3=1.0;
    cff4=1.0;
    for ised=1:NCS
      count=['0',num2str(ised)];
      count=count(end-1:end);
      eval(['cff1=cff1*mud_Sd50(ised)^squeeze(mudfrac_',count,'(1,1,j,i));'])
      eval(['cff2=cff2*mud_Srho(ised)^squeeze(mudfrac_',count,'(1,1,j,i));'])
      eval(['cff3=cff3*mud_Wsed(ised)^squeeze(mudfrac_',count,'(1,1,j,i));'])
      eval(['cff4=cff4*mud_tau_ce(ised)^squeeze(mudfrac_',count,'(1,1,j,i));'])
    end
    for ised=1:NNS
      count=['0',num2str(ised)];
      count=count(end-1:end);
      eval(['cff1=cff1*sand_Sd50(ised)^squeeze(sandfrac_',count,'(1,1,j,i));'])
      eval(['cff2=cff2*sand_Srho(ised)^squeeze(sandfrac_',count,'(1,1,j,i));'])
      eval(['cff3=cff3*sand_Wsed(ised)^squeeze(sandfrac_',count,'(1,1,j,i));'])
      eval(['cff4=cff4*sand_tau_ce(ised)^squeeze(sandfrac_',count,'(1,1,j,i));'])
    end
    grain_diameter(time,j,i)=cff1;
    grain_density(time,j,i)=cff2;
    settling_vel(time,j,i)=cff3;
    erosion_stress(time,j,i)=cff4;
    ripple_length(time,j,i)=0.10;
    ripple_height(time,j,i)=0.01;
    dmix_offset(time,j,i)=0.0;
    dmix_slope(time,j,i)=0.0;
    dmix_time(time,j,i)=0.0;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  END of USER INPUT                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%create init file
nc_init=netcdf(init_file,'clobber');
if isempty(nc_init), return, end
 
%% Global attributes:

disp(' ## Defining Global Attributes...')
nc_init.history = ncchar(['Created by "' mfilename '" on ' datestr(now)]);
nc_init.type = ncchar('Initialization file from create_roms_init.m');

%% Dimensions:

disp(' ## Defining Dimensions...')
 
nc_init('xi_psi') = L;
nc_init('xi_rho') = LP;
nc_init('xi_u') = L;
nc_init('xi_v') = LP;

nc_init('eta_psi') = M;
nc_init('eta_rho') = MP;
nc_init('eta_u') = MP;
nc_init('eta_v') = M;

nc_init('s_rho') = N;
nc_init('s_w') = N+1;
nc_init('tracer') = NT;
nc_init('Nbed') = Nbed;
%nc('boundary') = 4;
nc_init('time')=length(init_time);

nc_init('one') = 1;
nc_init('two') = 2;

 
%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')
 
nc_init{'theta_b'} = ncdouble; %% 1 element.
nc_init{'theta_b'}.long_name = ncchar('S-coordinate bottom control parameter');
nc_init{'theta_b'}.units = ncchar('1');
 
nc_init{'theta_s'} = ncdouble; %% 1 element.
nc_init{'theta_s'}.long_name = ncchar('S-coordinate surface control parameter');
nc_init{'theta_s'}.units = ncchar('1');

nc_init{'Tcline'} = ncdouble; %% 1 element.
nc_init{'Tcline'}.long_name = ncchar('S-coordinate surface/bottom layer width');
nc_init{'Tcline'}.units = ncchar('meter');
 
nc_init{'hc'} = ncdouble; %% 1 element.
nc_init{'hc'}.long_name = ncchar('S-coordinate parameter, critical depth');
nc_init{'hc'}.units = ncchar('meter');

nc_init{'Cs_r'} = ncdouble('s_rho');
nc_init{'Cs_r'}.long_name = ncchar('S-coordinate stretching curves at RHO-points');
nc_init{'Cs_r'}.units = ncchar('1');
nc_init{'Cs_r'}.valid_min = ncdouble(-1);
nc_init{'Cs_r'}.valid_max = ncdouble(0);
nc_init{'Cs_r'}.field = ncchar('Cs_r, scalar');
 
nc_init{'Cs_w'} = ncdouble('s_w');
nc_init{'Cs_w'}.long_name = ncchar('S-coordinate stretching curves at W-points');
nc_init{'Cs_w'}.units = ncchar('1');
nc_init{'Cs_w'}.valid_min = ncdouble(-1);
nc_init{'Cs_w'}.valid_max = ncdouble(0);
nc_init{'Cs_w'}.field = ncchar('Cs_w, scalar');

nc_init{'sc_r'} = ncdouble('s_rho');
nc_init{'sc_r'}.long_name = ncchar('S-coordinate at RHO-points');
nc_init{'sc_r'}.units = ncchar('1');
nc_init{'sc_r'}.valid_min = ncdouble(-1);
nc_init{'sc_r'}.valid_max = ncdouble(0);
nc_init{'sc_r'}.field = ncchar('sc_r, scalar');
 
nc_init{'sc_w'} = ncdouble('s_w');
nc_init{'sc_w'}.long_name = ncchar('S-coordinate at W-points');
nc_init{'sc_w'}.units = ncchar('1');
nc_init{'sc_w'}.valid_min = ncdouble(-1);
nc_init{'sc_w'}.valid_max = ncdouble(0);
nc_init{'sc_w'}.field = ncchar('sc_w, scalar');

nc_init{'ocean_time'} = ncdouble('time'); %% 1 element.
nc_init{'ocean_time'}.long_name = ncchar('time since initialization');
nc_init{'ocean_time'}.units = ncchar('second');
nc_init{'ocean_time'}.field = ncchar('ocean_time, scalar, series');

nc_init{'salt'} = ncfloat('time', 's_rho', 'eta_rho', 'xi_rho');
nc_init{'salt'}.long_name = ncchar('salinity');
nc_init{'salt'}.units = ncchar('PSU');
nc_init{'salt'}.field = ncchar('salinity, scalar, series');
 
nc_init{'temp'} = ncfloat('time', 's_rho', 'eta_rho', 'xi_rho'); 
nc_init{'temp'}.long_name = ncchar('temperature');
nc_init{'temp'}.units = ncchar('C');
nc_init{'temp'}.field = ncchar('temperature, scalar, series');
 
nc_init{'u'} = ncfloat('time', 's_rho', 'eta_u', 'xi_u');
nc_init{'u'}.long_name = ncchar('u-momentum component');
nc_init{'u'}.units = ncchar('meter second-1');
nc_init{'u'}.field = ncchar('u-velocity, scalar, series');

nc_init{'ubar'} = ncfloat('time', 'eta_u', 'xi_u');
nc_init{'ubar'}.long_name = ncchar('vertically integrated u-momentum component');
nc_init{'ubar'}.units = ncchar('meter second-1');
nc_init{'ubar'}.field = ncchar('ubar-velocity, scalar, series');
 
nc_init{'v'} = ncfloat('time', 's_rho', 'eta_v', 'xi_v');
nc_init{'v'}.long_name = ncchar('v-momentum component');
nc_init{'v'}.units = ncchar('meter second-1');
nc_init{'v'}.field = ncchar('v-velocity, scalar, series');
 
nc_init{'vbar'} = ncfloat('time', 'eta_v', 'xi_v');
nc_init{'vbar'}.long_name = ncchar('vertically integrated v-momentum component');
nc_init{'vbar'}.units = ncchar('meter second-1');
nc_init{'vbar'}.field = ncchar('vbar-velocity, scalar, series');
 
nc_init{'zeta'} = ncfloat('time', 'eta_rho', 'xi_rho');
nc_init{'zeta'}.long_name = ncchar('free-surface');
nc_init{'zeta'}.units = ncchar('meter');
nc_init{'zeta'}.field = ncchar('free-surface, scalar, series');
 
for mm=1:NCS
  count=['00',num2str(mm)];
  count=count(end-1:end);

  eval(['nc_init{''mud_',count,'''} = ncdouble(''time'', ''s_rho'', ''eta_rho'', ''xi_rho'');'])
  eval(['nc_init{''mud_',count,'''}.long_name = ncchar(''suspended cohesive sediment, size class ',count,''');'])
  eval(['nc_init{''mud_',count,'''}.units = ncchar(''kilogram meter-3'');'])
  eval(['nc_init{''mud_',count,'''}.time = ncchar(''ocean_time'');'])
  eval(['nc_init{''mud_',count,'''}.field = ncchar(''mud_',count,', scalar, series'');'])

  eval(['nc_init{''mudfrac_',count,'''} = ncdouble(''time'', ''Nbed'', ''eta_rho'', ''xi_rho'');'])
  eval(['nc_init{''mudfrac_',count,'''}.long_name = ncchar(''cohesive sediment fraction, size class ',count,''');'])
  eval(['nc_init{''mudfrac_',count,'''}.units = ncchar(''nondimensional'');'])
  eval(['nc_init{''mudfrac_',count,'''}.time = ncchar(''ocean_time'');'])
  eval(['nc_init{''mudfrac_',count,'''}.field = ncchar(''mudfrac_',count,', scalar, series'');'])

  eval(['nc_init{''mudmass_',count,'''} = ncdouble(''time'', ''Nbed'', ''eta_rho'', ''xi_rho'');'])
  eval(['nc_init{''mudmass_',count,'''}.long_name = ncchar(''cohesive sediment mass, size class ',count,''');'])
  eval(['nc_init{''mudmass_',count,'''}.units = ncchar(''kilogram meter-2'');'])
  eval(['nc_init{''mudmass_',count,'''}.time = ncchar(''ocean_time'');'])
  eval(['nc_init{''mudmass_',count,'''}.field = ncchar(''mudmass_',count,', scalar, series'');'])
end
for mm=1:NNS
  count=['00',num2str(mm)];
  count=count(end-1:end);

  eval(['nc_init{''sand_',count,'''} = ncdouble(''time'', ''s_rho'', ''eta_rho'', ''xi_rho'');'])
  eval(['nc_init{''sand_',count,'''}.long_name = ncchar(''suspended noncohesive sediment, size class ',count,''');'])
  eval(['nc_init{''sand_',count,'''}.units = ncchar(''kilogram meter-3'');'])
  eval(['nc_init{''sand_',count,'''}.time = ncchar(''ocean_time'');'])
  eval(['nc_init{''sand_',count,'''}.field = ncchar(''sand_',count,', scalar, series'');'])

  eval(['nc_init{''sandfrac_',count,'''} = ncdouble(''time'', ''Nbed'', ''eta_rho'', ''xi_rho'');'])
  eval(['nc_init{''sandfrac_',count,'''}.long_name = ncchar(''noncohesive sediment fraction, size class ',count,''');'])
  eval(['nc_init{''sandfrac_',count,'''}.units = ncchar(''nondimensional'');'])
  eval(['nc_init{''sandfrac_',count,'''}.time = ncchar(''ocean_time'');'])
  eval(['nc_init{''sandfrac_',count,'''}.field = ncchar(''sandfrac_',count,', scalar, series'');'])

  eval(['nc_init{''sandmass_',count,'''} = ncdouble(''time'', ''Nbed'', ''eta_rho'', ''xi_rho'');'])
  eval(['nc_init{''sandmass_',count,'''}.long_name = ncchar(''noncohesive sediment mass, size class ',count,''');'])
  eval(['nc_init{''sandmass_',count,'''}.units = ncchar(''kilogram meter-2'');'])
  eval(['nc_init{''sandmass_',count,'''}.time = ncchar(''ocean_time'');'])
  eval(['nc_init{''sandmass_',count,'''}.field = ncchar(''sandmass_',count,', scalar, series'');'])
end

nc_init{'bed_thickness'} = ncdouble ('time', 'Nbed', 'eta_rho', 'xi_rho') ;
nc_init{'bed_thickness'}.long_name = ncchar('sediment layer thickness');
nc_init{'bed_thickness'}.units = ncchar('meter');
nc_init{'bed_thickness'}.time = ncchar('ocean_time');
nc_init{'bed_thickness'}.field = ncchar('bed_thickness, scalar, series');

nc_init{'bed_age'} = ncdouble ('time', 'Nbed', 'eta_rho', 'xi_rho') ;
nc_init{'bed_age'}.long_name = ncchar('sediment layer age');
nc_init{'bed_age'}.units = ncchar('day');
nc_init{'bed_age'}.time = ncchar('ocean_time');
nc_init{'bed_age'}.field = ncchar('bed_age, scalar, series');

nc_init{'bed_porosity'} = ncdouble ('time', 'Nbed', 'eta_rho', 'xi_rho') ;
nc_init{'bed_porosity'}.long_name = ncchar('sediment layer porosity');
nc_init{'bed_porosity'}.units = ncchar('nondimensional');
nc_init{'bed_porosity'}.time = ncchar('ocean_time');
nc_init{'bed_porosity'}.field = ncchar('bed_porosity, scalar, series');

nc_init{'bed_biodiff'} = ncdouble ('time', 'Nbed', 'eta_rho', 'xi_rho') ;
nc_init{'bed_biodiff'}.long_name = ncchar('biodiffusivity at bottom of each layer');
nc_init{'bed_biodiff'}.units = ncchar('meter2 second-1');
nc_init{'bed_biodiff'}.time = ncchar('ocean_time');
nc_init{'bed_biodiff'}.field = ncchar('bed_biodiff, scalar, series');

nc_init{'grain_diameter'} = ncdouble ('time', 'eta_rho', 'xi_rho') ;
nc_init{'grain_diameter'}.long_name = ncchar('sediment median grain diameter size');
nc_init{'grain_diameter'}.units = ncchar('meter');
nc_init{'grain_diameter'}.time = ncchar('ocean_time');
nc_init{'grain_diameter'}.field = ncchar('grain_diameter, scalar, series');

nc_init{'grain_density'} = ncdouble ('time', 'eta_rho', 'xi_rho') ;
nc_init{'grain_density'}.long_name = ncchar('sediment median grain density');
nc_init{'grain_density'}.units = ncchar('kilogram meter-3');
nc_init{'grain_density'}.time = ncchar('ocean_time');
nc_init{'grain_density'}.field = ncchar('grain_density, scalar, series');

nc_init{'settling_vel'} = ncdouble ('time', 'eta_rho', 'xi_rho') ;
nc_init{'settling_vel'}.long_name = ncchar('sediment median grain settling velocity');
nc_init{'settling_vel'}.units = ncchar('meter second-1');
nc_init{'settling_vel'}.time = ncchar('ocean_time');
nc_init{'settling_vel'}.field = ncchar('settling_vel, scalar, series');

nc_init{'erosion_stress'} = ncdouble ('time', 'eta_rho', 'xi_rho') ;
nc_init{'erosion_stress'}.long_name = ncchar('sediment median critical erosion stress');
nc_init{'erosion_stress'}.units = ncchar('meter2 second-2');
nc_init{'erosion_stress'}.time = ncchar('ocean_time');
nc_init{'erosion_stress'}.field = ncchar('erosion_stress, scalar, series');

nc_init{'ripple_length'} = ncdouble ('time', 'eta_rho', 'xi_rho') ;
nc_init{'ripple_length'}.long_name = ncchar('bottom ripple length');
nc_init{'ripple_length'}.units = ncchar('meter');
nc_init{'ripple_length'}.time = ncchar('ocean_time');
nc_init{'ripple_length'}.field = ncchar('ripple length, scalar, series');

nc_init{'ripple_height'} = ncdouble ('time', 'eta_rho', 'xi_rho') ;
nc_init{'ripple_height'}.long_name = ncchar('bottom ripple height');
nc_init{'ripple_height'}.units = ncchar('meter');
nc_init{'ripple_height'}.time = ncchar('ocean_time');
nc_init{'ripple_height'}.field = ncchar('ripple height, scalar, series');

nc_init{'dmix_offset'} = ncdouble ('time', 'eta_rho', 'xi_rho') ;
nc_init{'dmix_offset'}.long_name = ncchar('dmix erodibility profile offset');
nc_init{'dmix_offset'}.units = ncchar('meter');
nc_init{'dmix_offset'}.time = ncchar('ocean_time');
nc_init{'dmix_offset'}.field = ncchar('dmix_offset, scalar, series');

nc_init{'dmix_slope'} = ncdouble ('time', 'eta_rho', 'xi_rho') ;
nc_init{'dmix_slope'}.long_name = ncchar('dmix erodibility profile slope');
nc_init{'dmix_slope'}.units = ncchar('_');
nc_init{'dmix_slope'}.time = ncchar('ocean_time');
nc_init{'dmix_slope'}.field = ncchar('dmix_slope, scalar, series');

nc_init{'dmix_time'} = ncdouble ('time', 'eta_rho', 'xi_rho') ;
nc_init{'dmix_time'}.long_name = ncchar('dmix erodibility profile time scale');
nc_init{'dmix_time'}.units = ncchar('seconds');
nc_init{'dmix_time'}.time = ncchar('ocean_time');
nc_init{'dmix_time'}.field = ncchar('dmix_time, scalar, series');

%now write the data from the arrays to the netcdf file
disp(' ## Filling Variables in netcdf file with data...')

nc_init{'theta_s'}(:) = theta_s;
nc_init{'theta_b'}(:) = theta_b;
nc_init{'Tcline'}(:)  = Tcline;
nc_init{'Cs_r'}(:) = Cs_r;
nc_init{'Cs_w'}(:) = Cs_w;
nc_init{'sc_w'}(:) = sc_w;
nc_init{'sc_r'}(:) = sc_r;
nc_init{'hc'}(:) = hc;

nc_init{'ocean_time'}(:) = init_time;
nc_init{'temp'}(:) = temp;
nc_init{'salt'}(:) = salt;
nc_init{'u'}(:)    = u;
nc_init{'ubar'}(:) = ubar;
nc_init{'v'}(:)    = v;
nc_init{'vbar'}(:) = vbar;
nc_init{'zeta'}(:) = zeta;

for mm=1:NCS
  count=['00',num2str(mm)];
  count=count(end-1:end);
  eval(['nc_init{''mud_',count,'''}(:)     = mud_',count,';'])           %sed conc in water column
  eval(['nc_init{''mudfrac_',count,'''}(:) = mudfrac_',count,';'])       %fraction of each sed class in each bed cell
  eval(['nc_init{''mudmass_',count,'''}(:) = mudmass_',count,';'])       %mass of each sed class in each bed cell
end
for mm=1:NNS
  count=['00',num2str(mm)];
  count=count(end-1:end);
  eval(['nc_init{''sand_',count,'''}(:)     = sand_',count,';'])           %sed conc in water column
  eval(['nc_init{''sandfrac_',count,'''}(:) = sandfrac_',count,';'])       %fraction of each sed class in each bed cell
  eval(['nc_init{''sandmass_',count,'''}(:) = sandmass_',count,';'])       %mass of each sed class in each bed cell
end

nc_init{'bed_thickness'}(:) = bed_thickness;
nc_init{'bed_age'}(:)       = bed_age;
nc_init{'bed_porosity'}(:)  = bed_porosity;
nc_init{'bed_biodiff'}(:)   = bed_biodiff;

nc_init{'grain_diameter'}(:) = grain_diameter;
nc_init{'grain_density'}(:)  = grain_density;
nc_init{'settling_vel'}(:)   = settling_vel;
nc_init{'erosion_stress'}(:) = erosion_stress;
nc_init{'ripple_height'}(:) = ripple_height;
nc_init{'ripple_length'}(:) = ripple_length;
nc_init{'dmix_offset'}(:) = dmix_offset;
nc_init{'dmix_slope'}(:) = dmix_slope;
nc_init{'dmix_time'}(:) = dmix_time;

%close file
close(nc_init)
disp(['created ', init_file])


