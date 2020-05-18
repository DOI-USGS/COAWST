function Tide = otps2roms(gfile,base_date,pred_date,ofile,model_file)

%
% OTPS2ROMS:  Generates OTPS tidal forcing file for ROMS
%
% Tide = otps2roms(gfile,base_date,pred_date,ofile,model_file)
%
% Generates a ROMS tidal forcing from OTPS dataset. The 'base_date'
% is related to ROMS input parameter the TIDE_START. It is used in
% ROMS to compute the correct phase lag with respect initialization
% time.
%
% On Input:
%
%    gfile       ROMS NetCDF Grid file name (string)
%
%    base_date   Tidal data reference day (serial date number), for
%                  example:
%
%                  base_date = datenum('01-Jan-2000')
%
%    pred_date   Prediction date for nodal corrections (serial date
%                  number)
%
%    ofile       ROMS tidal forcing NetCDF file name (string)
%
%    model_file  OTPS file name (string, full or relative path) for
%                  global or regional dataset. For example,
%
%                  'DATA/Model_tpxo7'          Global file (relative path)
%                  '~/OTPS/DATA/Model_tpxo7'   Global file (full path)
%
%                  'DATA/Model_Mex'            GOM file (relative path)
%                  '~/OTPS/DATA/Model_Mex'     GOM file (relative path)
%
% On Output:
%
%    Tide        Extracted tide data (struct array):
%
%                  Tide.period      Tide angular period (hour)
%                  Tide.names       Tide constituents name (string)
%                  Tide.Ephase      Tide elevation phase angle (degree)
%                  Tide.Eamp        Tide elevation amplitude (m)
%                  Tide.Cmax        Maximum tidal current (m/s)
%                  Tide.Cmin        Minimum tidal current (m/s)
%                  Tide.Cangle      Tide current inclination angle (degree)
%                  Tide.Cphase      Tide current phase angle (degree)
%
% Requirements:
%
%    T_TIDE           Tidal analysis package
%    tidal_ellipse    (Zhigang Xu) package, ap2ep.m
%    extract_HC.f     unpacking program to be compiled as extract_HC
%
% History:
%
% Originally written by John Evans
% Revised by Eli Hunter    03/07/2007
% Revised by Eli Hunter    05/25/2007  Corrected phase lag error)
% Revised by Eli Hunter    08/02/2011  Added 'nodal', to t_vuf and changed
%                                      Tide.period
% Revised by Hernan Arango 08/02/2012
%
  
% svn $Id: otps2roms.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Collect the information necessary to run extract_HC and create setup
% files. 
%--------------------------------------------------------------------------

vars=['z','u','v'];
harmonics_prefix='otps_harmonics_';

% Build grid file.

llfile='ll.dat';
roms2ll(gfile,llfile);
disp(blanks(1));
disp(['Gridfile ' gfile ' lat/lon written to ' llfile]);

if ~(exist(model_file,'file'))
  error(['No Such model file: ' model_file]);
end

disp(blanks(1));
disp(['Mode parameter file used: ' model_file]);

for i=1:length(vars)
  disp(blanks(1));
  disp(['Generating parameter file for ' vars(i)])
  fid=fopen('otps_input','w');
  fprintf(fid,'%s\n',model_file);
  fprintf(fid,'%s\n',llfile);
  fprintf(fid,'%s\n',vars(i));
  fprintf(fid,'%s\n','m2,s2,k1,o1,');
  fprintf(fid,'%s\n','AP');
  fprintf(fid,'%s\n','oce');
  fprintf(fid,'%s\n','0');
  fprintf(fid,'%s\n',[harmonics_prefix vars(i)]);
  fclose(fid);
  
  disp(['Extracting ' vars(i) ' Harmonics'])
  
  [s,w]=unix('~/bin/extract_HC < otps_input');
  disp([' '])
% disp(['*******************************************'])
  disp(w)
% disp(['*******************************************'])
end  

mask_rho = nc_varget ( gfile, 'mask_rho' );
land = find(mask_rho==0);
water = find(mask_rho==1);

fprintf ( 1, 'Reading %s...\n', [harmonics_prefix vars(1)] );
[z_hc,lon,lat] = read_otps_output([harmonics_prefix vars(1)]);
[period, z_amp, z_phase, names] = reshape_to_grid ( z_hc, gfile );
fprintf ( 1, 'Reading %s...\n', [harmonics_prefix vars(2)]);
[u_hc,lon,lat] = read_otps_output([harmonics_prefix vars(2)]);
[period, u_amp, u_phase, names] = reshape_to_grid ( u_hc, gfile );
fprintf ( 1, 'Reading %s...\n', [harmonics_prefix vars(3)]);
[v_hc,lon,lat] = read_otps_output([harmonics_prefix vars(3)]);
[period, v_amp, v_phase, names] = reshape_to_grid ( v_hc, gfile );

cnames=upper(char(names));

%--------------------------------------------------------------------------
% Make sure that the OTPS mask agrees with the ROMS mask.
% Fill in any points that ROMS thinks is water but OTPS thinks is land.
%--------------------------------------------------------------------------

num_constituents = length(period);
lon_rho = nc_varget ( gfile, 'lon_rho' );
lat_rho = nc_varget ( gfile, 'lat_rho' );
mask_rho = nc_varget ( gfile, 'mask_rho' );
a=t_getconsts;

for j = 1:num_constituents
  component = squeeze ( z_amp(j,:,: ) );
  z_amp(j,:,:) = match_roms_mask ( lon_rho, lat_rho, mask_rho, component );

  component = squeeze ( z_phase(j,:,: ) );
  z_phase(j,:,:) = match_roms_mask ( lon_rho, lat_rho, mask_rho, component );

  component = squeeze ( u_amp(j,:,: ) );
  u_amp(j,:,:) = match_roms_mask ( lon_rho, lat_rho, mask_rho, component );

  component = squeeze ( u_phase(j,:,: ) );
  u_phase(j,:,:) = match_roms_mask ( lon_rho, lat_rho, mask_rho, component );

  component = squeeze ( v_amp(j,:,: ) );
  v_amp(j,:,:) = match_roms_mask ( lon_rho, lat_rho, mask_rho, component );

  component = squeeze ( v_phase(j,:,: ) );
  v_phase(j,:,:) = match_roms_mask ( lon_rho, lat_rho, mask_rho, component );

  iconst(j)=strmatch(cnames(j,:), a.name);
  Tide.period(j)=1/a.freq(iconst(j));
end

Tide.names = names;

Ntide = length(z_hc);
rg = roms_get_grid ( gfile );
[Lp,Mp] = size(rg.lon_rho);

z_amp = zero_out_land ( z_amp, land );
z_phase = zero_out_land ( z_phase, land );
	
%--------------------------------------------------------------------------
% Correct phase lag to the specified base_time. Also, the amplitude
% is corrected for nodal adjustment. Use T_TIDE function 't_vuf'.
%--------------------------------------------------------------------------
%
% The reference latitude for 3rd order satellites (degrees) is set to 55.
% You don't need to adjust this to your local latitude. It could also be
% set to NaN as in Xtide, with very little effect. See 't_vuf' for more
% information.

reflat=55;
datestr(base_date)
iconst
reflat

% Correct phase lag.

[V,U,F]=t_vuf('nodal',base_date,iconst,reflat);

% Apply nodal correction. 

[Vp,Up,Fp] = t_vuf('nodal',pred_date,iconst,reflat);

% The tidal currents are returned in cycles, so multiply by 360 to get
% degrees or 2*pi to get radians.

V  = V *360;               % convert phase to degrees
U  = U *360;               % convert phase to degrees
Vp = Vp*360;               % convert phase to degrees
Up = Up*360;               % convert phase to degrees

for k=1:Ntide;
  z_phase(k,:,:) = z_phase(k,:,:) - Up(k)  - V(k);   % degrees
  z_amp(k,:,:) =z_amp(k,:,:) .* Fp(k);

  u_phase(k,:,:) =u_phase(k,:,:) - Up(k) - V(k);   % degrees
  u_amp(k,:,:) = u_amp(k,:,:) .* Fp(k);
    
  v_phase(k,:,:) = v_phase(k,:,:) - Up(k)  - V(k);   % degrees
  v_amp(k,:,:) =v_amp(k,:,:) .* Fp(k);
end

z_phase=mod(z_phase,360);
u_phase=mod(u_phase,360);
v_phase=mod(v_phase,360);

Tide.Ephase    = z_phase(:,:,:);
Tide.Eamp      = z_amp(:,:,:);

%--------------------------------------------------------------------------
%  Convert tidal current amplitude and phase lag parameters to tidal
%  current ellipse parameters: Major axis, ellipticity, inclination,
%  and phase.  Use "tidal_ellipse" (Zhigang Xu) package.
%--------------------------------------------------------------------------

[major,eccentricity,inclination,phase]=ap2ep(u_amp,u_phase,v_amp,v_phase);
	
major = zero_out_land ( major, land );
eccentricity = zero_out_land ( eccentricity, land );
major = major/100;
Tide.Cmax=major;
Tide.Cmin=major.*eccentricity;
Tide.Cangle= zero_out_land ( inclination, land );
Tide.Cphase = zero_out_land ( phase, land );

%--------------------------------------------------------------------------
%  Write out extracted into NetCDF file
%--------------------------------------------------------------------------

write_tides(Tide, gfile, ofile, base_date, model_file);

return
