function create_roms_netcdf_init(init_file,gn,t_clim,Nbed,NNS,NCS)

%create init file
nc_init=netcdf(init_file,'clobber');
 
%% Global attributes:

disp(' ## Defining Global Attributes...')
nc_init.history = ncchar(['Created by "' mfilename '" on ' datestr(now)]);
nc_init.type = ncchar('Initialization file from create_roms_init.m');

%% Dimensions:

disp(' ## Defining Dimensions...')
 
%get some grid info
  [MP,LP]=size(gn.lon_rho);
  L=LP-1;
  Lm=L-1;
  M=MP-1;
  Mm=M-1;
  L  = Lm+1;
  M  = Mm+1;
  xi_psi  = L;
  xi_rho  = LP;
  xi_u    = L;
  xi_v    = LP;
  eta_psi = M;
  eta_rho = MP;
  eta_u   = MP;
  eta_v   = M;
  N       = gn.N;

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
nc_init('Nbed') = Nbed;
nc_init('time')=1;

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
nc_init{'ocean_time'}.units = ncchar('days');
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

close(nc_init)



