% ROMS_TRI_TIDES:  Makes ROMS tidal forcing from triangle based tidal data
% This script interpolates tidal elevation and current data from 
% triangle-based tidal data.  This could be output from a Finite Element 
% Model such as ADCIRC or QUODDY or any data run through DELAUNAY triangulation.
%
% User must modify the file names, start time and tidal constituents 
% for the particular application in the user-modified section below.  
%
% Rich Signell (rsignell@usgs.gov)
% John Warner  (jcwarner@usgs.gov)
% 
% 12-Jan-2005
% 17Jul2012 - update to use matlab native netcdf interface
%           - remove tri dependency, use TriScatteredInterp instead
%
% Requires: T_TIDE and TIDAL_ELLIPSE packages 
 
%%%%% BEGINNING of USER-MODIFIED SECTION %%%%%
IPLOT=1;            % 1 to make plots, 0 for no plots
IWRITE=1;           % 1 to write output to netcdf, 0 for no output

% (1) Specify existing Grid File and new tide Forcing Files 
      Gname='C:\work\models\COAWST\Projects\Sandy2\Sandy_roms_grid.nc';
      Fname='tide_forc_Sandy.nc';

% (2) Enter ROMS start time.  This will be used to calculate the proper phase
%     for the tidal constituents.  If you change the start time, this routine
%     must be run again to adjust the tidal phases. 
%            YYYY   MO   DA    HR     MI    SC
      g = [  2012,  10,  28,   0,     0,    0];   % GMT
      disp(['Tidal Start Time =' datestr(g)])

% (3) Select Adcirc or topex/osu tidal constituent data 
      adcirc=1
      osu=0
      if (adcirc)
        load c:\work\models\COAWST_data\tide\adcirc_ec2001v2e_fix.mat tri lon lat u v elev depth periods freq names
        names=names';
      end
      if (osu)
        load ('COAWST_data\tide\tpx_uv.mat')
        load ('COAWST_data\tide\tpx_h.mat')
        names=con;
      end
% (4) Specify which tidal constituents to use.  Look at "names" in the .mat file to
%     see what is available.  Program will halt if these are not found.
      if (adcirc)
        tides_to_use=['K1'
                     'O1'
                     'Q1'
                     'M2'
                     'S2'
                     'N2'
                     'K2'
                     'M4'
                     'M6'];
      end
      if (osu)
        tides_to_use=['m2  '
                      's2  '
                      'n2  '
                      'k2  '
                      'k1  '
                      'o1  '
                      'p1  '
                      'q1  '
                      'mf  '
                      'mm  '
                      'm4  '
                      'ms4 '
                      'mn4 '];
%       periods=[12.4211  12.0000  12.6584 11.9672 23.9344 25.8192 24.0659 ...
%                26.8685 327.8689 661.3392  6.2103 6.1033   6.2692];
        periods=[12.420601221208152  12.000000002711209  12.658348243602962 ...
                 11.967234788246559  23.934469655266863  25.819341665032908 ...
                 24.065890161111373  26.868356672318782 327.8589689145072   ...
                661.3096907770995     6.210300610604076   6.103339299686573 ...
                  6.269173876477840];
      end
%%%%% END of USER-MODIFIED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(adcirc)
  % make sure that tidal data is double 
  tri=double(tri);
  lon=double(lon);
  lat=double(lat);
  u=double(u);
  v=double(v);
  elev=double(elev);
  depth=double(depth);
end

[Ntide,ncol]=size(tides_to_use);
for k=1:Ntide
  imatch=strmatch(tides_to_use(k,:),names');
  if isempty(imatch),
    disp(['Error: Did not find constituent ' tides_to_use(k,:) ' in .mat file']);
    return
  end
  indx(k)=imatch;
end

% create the netcdf tide forc file
h=ncread(Gname,'h');
[Lp,Mp]=size(h);
roms_tide_forc_file(Fname,Lp,Mp,periods(indx));

%---------------------------------------------------------------------
%  Read in application grid data.
%---------------------------------------------------------------------

rlon=ncread(Gname,'lon_rho');
rlat=ncread(Gname,'lat_rho');
[Lp,Mp]=size(rlon);

rmask =ncread(Gname,'mask_rho');
rangle=ncread(Gname,'angle');

% do only water points
% [iwater]=find(~rmask);
% do all points
if (adcirc)
  iwater=[1:size(rmask,1)*size(rmask,2)];
end
% Calculate V,U,F arguments using Rich P's T_Tide.  
% These are the astronomical offsets
% that allow you to go from  amplitude and Greenwich phase to 
% real tidal predictions for a certain time period.  See the
% T_TIDE toolkit and the T_VUF.M program for more info.

a=t_getconsts;

%---------------------------------------------------------------------
%  Interpolate tide data to application grid.
%---------------------------------------------------------------------
%
%  Initialize arrays for interpolated tidal data.

if (adcirc)
  ei=ones(Lp,Mp,Ntide).*NaN;                       % elevation
  ui=ones(Lp,Mp,Ntide).*NaN;                       % u 
  vi=ones(Lp,Mp,Ntide).*NaN;                       % v

  tmp=ones(Lp,Mp).*NaN;

  %  Interpolate elevation and velocity from ADCIRC model data using
  %  "triterp", linear-interpolation from triangulated grid.
  for kk=1:Ntide,
    k=indx(kk);
    disp(['Interpolating ' names(:,k)']);
    format loose

%    tmp(iwater)=triterp(tri,lon,lat,elev(:,k),rlon(iwater),rlat(iwater));
%    ei(:,:,kk)=tmp;
%    tmp(iwater)=triterp(tri,lon,lat,u(:,k),rlon(iwater),rlat(iwater));
%    ui(:,:,kk)=tmp;
%    tmp(iwater)=triterp(tri,lon,lat,v(:,k),rlon(iwater),rlat(iwater));
%    vi(:,:,kk)=tmp;

    F=TriScatteredInterp(lon,lat,elev(:,k));
    zz=F(rlon(iwater),rlat(iwater));
    ei(:,:,kk)=reshape(zz,Lp,Mp);

    F=TriScatteredInterp(lon,lat,u(:,k));
    zz=F(rlon(iwater),rlat(iwater));
    ui(:,:,kk)=reshape(zz,Lp,Mp);
    
    F=TriScatteredInterp(lon,lat,v(:,k));
    zz=F(rlon(iwater),rlat(iwater));
    vi(:,:,kk)=reshape(zz,Lp,Mp);
    
  % find T_TIDE constituent names (a.name) that match 
  % selected constituent names (names(tides_to_use))
    iconst(kk)=strmatch(names(:,k), a.name);
  end
end

if (osu)
  eip=ones(Lp,Mp,Ntide).*NaN;                       % elevation
  eia=ones(Lp,Mp,Ntide).*NaN;                       % elevation
  uip=ones(Lp,Mp,Ntide).*NaN;                       % u 
  uia=ones(Lp,Mp,Ntide).*NaN;                       % u 
  vip=ones(Lp,Mp,Ntide).*NaN;                       % v
  via=ones(Lp,Mp,Ntide).*NaN;                       % v

  %  Interpolate elevation and velocity from TOPEX7.2 data.
  ua=ua/100;
  va=va/100;
 %rlon=rlon+360;
  rlon(rlon<=0) = rlon(rlon<=0) + 360.00;
  for kk=1:Ntide,
    k=indx(kk);
    disp(['Interpolating ' names(:,k)']);
    format loose

    tmp2=squeeze(ha(:,:,k));
    tmp=griddata(lon_z',lat_z',tmp2',rlon,rlat);
    eia(:,:,kk)=tmp;

    tmp2=squeeze(hp(:,:,k));
    tmp=griddata(lon_z',lat_z',tmp2',rlon,rlat);
    eip(:,:,kk)=tmp;

    tmp2=squeeze(ua(:,:,k));
    tmp=griddata(lon_u',lat_u',tmp2',rlon,rlat);
    uia(:,:,kk)=tmp;

    tmp2=squeeze(up(:,:,k));
    tmp=griddata(lon_u',lat_u',tmp2',rlon,rlat);
    uip(:,:,kk)=tmp;

    tmp2=squeeze(va(:,:,k));
    tmp=griddata(lon_v',lat_v',tmp2',rlon,rlat);
    via(:,:,kk)=tmp;

    tmp2=squeeze(vp(:,:,k));
    tmp=griddata(lon_v',lat_v',tmp2',rlon,rlat);
    vip(:,:,kk)=tmp;

  % find T_TIDE constituent names (a.name) that match 
  % selected constituent names (names(tides_to_use))
    iconst(kk)=strmatch(upper(names(:,k)), a.name);
  end
end

clear tmp

%***********************************************************************
% This is the call to t_vuf that
% will correct the phase to be at the user specified time.  Also, the amplitude
% is corrected for nodal adjustment.

start_year=g(1,1);
start_month=g(1,2);
start_day=g(1,3);
start_hour=g(1,4);
start_min=g(1,5);
start_sec=g(1,6);

omega=2*pi.*a.freq(iconst);   %tidal frequencies in radians/hour

% Reference latitude for 3rd order satellites (degrees) is 
% set to 55.  You don't need to adjust this to your local latitude
% It could also be set to NaN as in Xtide, with very little effect.  
% See T_VUF for more info.

reflat=55;  
[vv,uu,ff]=t_vuf('nodal',datenum(start_year,start_month,start_day,start_hour,start_min,start_sec),iconst,reflat)
%vv and uu are returned in cycles, so * by 360 to get degrees or * by 2 pi to get radians

vv=vv*360;  % convert vv to phase in degrees
uu=uu*360;  % convert uu to phase in degrees

%adjust phase and ampl of tide to account for user specified time
ei_phase=ones(Lp,Mp,Ntide).*NaN;                       % elevation
ei_amplitude=ones(Lp,Mp,Ntide).*NaN;                       % elevation
ui_phase=ones(Lp,Mp,Ntide).*NaN;                       % u vel 
ui_amplitude=ones(Lp,Mp,Ntide).*NaN;                       %   u vel
vi_phase=ones(Lp,Mp,Ntide).*NaN;                       % v vel
vi_amplitude=ones(Lp,Mp,Ntide).*NaN;                       %   v vel

if (adcirc)
  for k=1:Ntide;
      ei_phase(:,:,k) = angle(ei(:,:,k)).*180./pi - vv(k) - uu(k);   % degrees
      ei_amplitude(:,:,k) = abs(ei(:,:,k)) .* ff(k);

      ui_phase(:,:,k) = angle(ui(:,:,k)).*180./pi - vv(k) - uu(k);   % degrees
      ui_amplitude(:,:,k) = abs(ui(:,:,k)) .* ff(k);
    
      vi_phase(:,:,k) = angle(vi(:,:,k)).*180./pi - vv(k) - uu(k);   % degrees
      vi_amplitude(:,:,k) = abs(vi(:,:,k)) .* ff(k);
  end
end
if (osu)
  for k=1:Ntide;
      ei_phase(:,:,k) = eip(:,:,k) - vv(k) - uu(k);   % degrees
      ei_amplitude(:,:,k) = eia(:,:,k) .* ff(k);

      ui_phase(:,:,k) = uip(:,:,k) - vv(k) - uu(k);   % degrees
      ui_amplitude(:,:,k) = uia(:,:,k) .* ff(k);
    
      vi_phase(:,:,k) = vip(:,:,k) - vv(k) - uu(k);   % degrees
      vi_amplitude(:,:,k) = via(:,:,k) .* ff(k);
  end
end


Tide.period=periods(indx);
Tide.component=a.name(iconst,:);
Tide.Ephase = ei_phase(:,:,:);
Tide.Eamp   = ei_amplitude(:,:,:);
clear e*
size(Tide.Eamp)

%---------------------------------------------------------------------
%  Convert tidal current amplitude and phase lag parameters to tidal
%  current ellipse parameters: Major axis, ellipticity, inclination,
%  and phase.  Use "tidal_ellipse" package by Zhigang Xu
%---------------------------------------------------------------------

%ui_pha=ui_phase(:,:,:);
%ui_amp=ui_amplitude(:,:,:);

%vi_pha=vi_phase(:,:,:);
%vi_amp=vi_amplitude(:,:,:);

%[major,eccentricity,inclination,phase]=ap2ep(ui_amp,ui_pha,vi_amp,vi_pha);
[major,eccentricity,inclination,phase]=ap2ep(ui_amplitude,ui_phase,vi_amplitude,vi_phase);

Tide.Cmax=major;
Tide.Cmin=major.*eccentricity;
Tide.Cangle=inclination;
Tide.Cphase=phase;

clear ui ui_amp ui_pha vi vi_amp vi_pha ei
clear major eccentricity inclination phase

%---------------------------------------------------------------------
%  Plot results.
%---------------------------------------------------------------------

if (IPLOT),

   for k = 1: Ntide
    figure

    subplot(1,2,1)
    pcolor(rlon,rlat,squeeze(Tide.Eamp(:,:,k)));
    shading('interp'); colorbar; dasp(43); grid on;
    xlabel('Amplitude (m)');
    title(strcat(Tide.component(k,:),' Tidal Component'));

    subplot(1,2,2)
    pcolor(rlon,rlat,squeeze(Tide.Ephase(:,:,k)));
    shading('interp'); colorbar; dasp(43); grid on;
    xlabel('Phase Angle (degree)');
    title(strcat(Tide.component(k,:),' Tidal Component'));
  
  end,

   for k=1:Ntide

    figure

    subplot(2,2,1)
    pcolor(rlon,rlat,squeeze(Tide.Cmax(:,:,k)));
    shading('interp'); colorbar; dasp(43); grid on;
    xlabel('Major Axis amplitude (m/s)');
    title(strcat(Tide.component(k,:),' Tidal Component'));

    subplot(2,2,2)
    pcolor(rlon,rlat,squeeze(Tide.Cmin(:,:,k)));
    shading('interp'); colorbar; dasp(43); grid on;
    xlabel('Minor Axis amplitude (m/s)');
  
    subplot(2,2,3)
    pcolor(rlon,rlat,squeeze(Tide.Cangle(:,:,k)));
    shading('interp'); colorbar; dasp(43); grid on;
    xlabel('Inclination Angle (degrees)');

    subplot(2,2,4)
    pcolor(rlon,rlat,squeeze(Tide.Cphase(:,:,k)));
    shading('interp'); colorbar; dasp(43); grid on;
    xlabel('Phase Angle(degrees)');

 end,

end,

%---------------------------------------------------------------------
%  Write out tide data into FORCING NetCDF file.
%---------------------------------------------------------------------

if (IWRITE),
  Tide.Ephase(isnan(Tide.Ephase))=0;
  Tide.Eamp(isnan(Tide.Eamp))=0;
  Tide.Cphase(isnan(Tide.Cphase))=0;
  Tide.Cangle(isnan(Tide.Cangle))=0;
  Tide.Cmin(isnan(Tide.Cmin))=0;
  Tide.Cmax(isnan(Tide.Cmax))=0;
  ncwrite(Fname,'tide_period',Tide.period);
  ncwrite(Fname,'tide_Ephase',Tide.Ephase);
  ncwrite(Fname,'tide_Eamp',Tide.Eamp);
  ncwrite(Fname,'tide_Cphase',Tide.Cphase);
  ncwrite(Fname,'tide_Cangle',Tide.Cangle);
  ncwrite(Fname,'tide_Cmin',Tide.Cmin);
  ncwrite(Fname,'tide_Cmax',Tide.Cmax);
end 
