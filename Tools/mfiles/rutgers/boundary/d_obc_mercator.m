%
%  D_OBC_MERCATOR:  Driver script to create a ROMS boundary conditions
%
%  This a user modifiable script that can be used to prepare ROMS open
%  boundary conditions NetCDF file from Mercator dataset. It sets-up all
%  the necessary parameters and variables. USERS can use this as a
%  prototype for their application.
%

% svn $Id: d_obc_mercator.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Set input/output NetCDF files.

OPR_Dir = '/home/arango/ocean/toms/repository/Projects/philex/Mercator/OPR';
RTR_Dir = '/home/arango/ocean/toms/repository/Projects/philex/Mercator/RTR';

 GRDname = 'philex_grd_b1.nc';
 BRYname = 'philex_bry_mercator_b1.nc';

CREATE = 1;                      % logical switch to create NetCDF
WRITE  = 1;                      % logical swithc to write out data
report = 1;                      % report vertical grid information

% Get number of grid points.

[Lr,Mr]=size(nc_read(GRDname,'h'));

Lu = Lr-1;   Lv = Lr;
Mu = Mr;     Mv = Mr-1;

%--------------------------------------------------------------------------
% Set full path of Mercator files for boundary conditions.
%--------------------------------------------------------------------------

Tfile1=dir(fullfile(RTR_Dir,'*_gridT_*.nc.gz'));
Ufile1=dir(fullfile(RTR_Dir,'*_gridU_*.nc.gz'));
Vfile1=dir(fullfile(RTR_Dir,'*_gridV_*.nc.gz'));

nfiles1=length(Tfile1);

Tfile2=dir(fullfile(OPR_Dir,'*_gridT_*.nc.gz'));
Ufile2=dir(fullfile(OPR_Dir,'*_gridU_*.nc.gz'));
Vfile2=dir(fullfile(OPR_Dir,'*_gridV_*.nc.gz'));

nfiles2=length(Tfile2);

Tfile=[Tfile1; Tfile2];
Ufile=[Ufile1; Ufile2];
Vfile=[Vfile1; Vfile2];

nfiles=length(Tfile);

%--------------------------------------------------------------------------
%  Set application parameters in structure array, S.
%--------------------------------------------------------------------------

S.ncname      = BRYname;     % output NetCDF file

S.spherical   = 1;           % spherical grid

S.boundary(1) = 1;           % process western  boundary segment (0=no)
S.boundary(2) = 1;           % process eastern  boundary segment (0=no)
S.boundary(3) = 1;           % process southern boundary segment (0=no)
S.boundary(4) = 1;           % process northern boundary segment (0=no)

S.Lm          = Lr-2;        % number of interior RHO-points, X-direction
S.Mm          = Mr-2;        % number of interior RHO-points, Y-direction
S.N           = 42;          % number of vertical levels at RHO-points
S.NT          = 2;           % total number of tracers

S.Vtransform  = 2;           % vertical transfomation equation
S.Vstretching = 2;           % vertical stretching function

S.theta_s     = 7.0;         % S-coordinate surface control parameter
S.theta_b     = 0.1;         % S-coordinate bottom control parameter
S.Tcline      = 150.0;       % S-coordinate surface/bottom stretching width
S.hc          = S.Tcline;    % S-coordinate stretching width

%--------------------------------------------------------------------------
%  Set grid variables.
%--------------------------------------------------------------------------

S.h           = nc_read(GRDname, 'h');            % bathymetry

S.lon_rho     = nc_read(GRDname, 'lon_rho');      % RHO-longitude
S.lat_rho     = nc_read(GRDname, 'lat_rho');      % RHO-latitude

S.lon_u       = nc_read(GRDname, 'lon_u');        % U-longitude
S.lat_u       = nc_read(GRDname, 'lat_u');        % U-latitude

S.lon_v       = nc_read(GRDname, 'lon_v');        % V-longitude
S.lat_v       = nc_read(GRDname, 'lat_v');        % V-latitude

S.mask_rho    = nc_read(GRDname, 'mask_rho');     % RHO-mask
S.mask_u      = nc_read(GRDname, 'mask_u');       % U-mask
S.mask_v      = nc_read(GRDname, 'mask_v');       % V-mask

S.angle       = nc_read(GRDname, 'angle');        % curvilinear angle

%  Set boundary conditions locations.

S.lon_rho_west =S.lon_rho(1,:);
S.lon_u_west   =S.lon_u(1,:);
S.lon_v_west   =S.lon_v(1,:);

S.lat_rho_west =S.lat_rho(1,:);
S.lat_u_west   =S.lat_u(1,:);
S.lat_v_west   =S.lat_v(1,:);

S.lon_rho_east =S.lon_rho(end,:);
S.lon_u_east   =S.lon_u(end,:);
S.lon_v_east   =S.lon_v(end,:);

S.lat_rho_east =S.lat_rho(end,:);
S.lat_u_east   =S.lat_u(end,:);
S.lat_v_east   =S.lat_v(end,:);

S.lon_rho_south=S.lon_rho(:,1);
S.lon_u_south  =S.lon_u(:,1);
S.lon_v_south  =S.lon_v(:,1);

S.lat_rho_south=S.lat_rho(:,1);
S.lat_u_south  =S.lat_u(:,1);
S.lat_v_south  =S.lat_v(:,1);

S.lon_rho_north=S.lon_rho(:,end);
S.lon_u_north  =S.lon_u(:,end);
S.lon_v_north  =S.lon_v(:,end);

S.lat_rho_north=S.lat_rho(:,end);
S.lat_u_north  =S.lat_u(:,end);
S.lat_v_north  =S.lat_v(:,end);

%  Set vertical grid variables.

kgrid=0;                                          % RHO-points

[S.s_rho, S.Cs_r]=stretching(S.Vstretching,                           ...
                             S.theta_s, S.theta_b, S.hc, S.N,         ...
                             kgrid, report);

kgrid=1;                                          % W-points

[S.s_w,   S.Cs_w]=stretching(S.Vstretching, ...
                             S.theta_s, S.theta_b, S.hc, S.N,         ...
                             kgrid, report);

%  Compute ROMS model depths.  Ignore free-sruface contribution
%  so interpolation is bounded below mean sea level.

ssh=zeros(size(S.h));

igrid=1;
[S.z_r]=set_depth(S.Vtransform, S.Vstretching,                        ...
                  S.theta_s, S.theta_b, S.hc, S.N,                    ...
                  igrid, S.h, ssh, report);
      
igrid=3;
[S.z_u]=set_depth(S.Vtransform, S.Vstretching,                        ...
                  S.theta_s, S.theta_b, S.hc, S.N,                    ...
                  igrid, S.h, ssh, report);

igrid=4;
[S.z_v]=set_depth(S.Vtransform, S.Vstretching,                        ...
                  S.theta_s, S.theta_b, S.hc, S.N,                    ...
                  igrid, S.h, ssh, report);

%  Compute ROMS vertical level thicknesses (m).
      
N=S.N;
igrid=5;
[S.z_w]=set_depth(S.Vtransform, S.Vstretching, ...
                  S.theta_s, S.theta_b, S.hc, S.N, ...
                  igrid, S.h, ssh, report);

S.Hz=S.z_w(:,:,2:N+1)-S.z_w(:,:,1:N);

%---------------------------------------------------------------------------
%  Create boundary condition Netcdf file.
%---------------------------------------------------------------------------

if (CREATE),

  [status]=c_boundary(S);

%  Set attributes for "bry_time".

  avalue='seconds since 2007-01-01 00:00:00';
  [status]=nc_attadd(BRYname,'units',avalue,'bry_time');
  
  avalue='gregorian';
  [status]=nc_attadd(BRYname,'calendar',avalue,'bry_time');

%  Set global attribute.

  avalue='Philippine Archipelago Straits, ~5.5 km resolution, Grid b';
  [status]=nc_attadd(BRYname,'title',avalue);

  avalue='Mercator system PSY3V2 daily average, 0.25 degree resolution';
  [status]=nc_attadd(BRYname,'source',avalue);

  [status]=nc_attadd(BRYname,'grd_file',GRDname);
  
%  Write out grid data.
  
  [status]=nc_write(BRYname, 'spherical',   S.spherical);
  [status]=nc_write(BRYname, 'Vtransform',  S.Vtransform);
  [status]=nc_write(BRYname, 'Vstretching', S.Vstretching);
  [status]=nc_write(BRYname, 'theta_s',     S.theta_s);
  [status]=nc_write(BRYname, 'theta_b',     S.theta_b);
  [status]=nc_write(BRYname, 'Tcline',      S.Tcline);
  [status]=nc_write(BRYname, 'hc',          S.hc);
  [status]=nc_write(BRYname, 's_rho',       S.s_rho);
  [status]=nc_write(BRYname, 's_w',         S.s_w);
  [status]=nc_write(BRYname, 'Cs_r',        S.Cs_r);
  [status]=nc_write(BRYname, 'Cs_w',        S.Cs_w);

  if (S.boundary(1)),
    [status]=nc_write(BRYname, 'lon_rho_west' , S.lon_rho_west);
    [status]=nc_write(BRYname, 'lat_rho_west' , S.lat_rho_west);
    [status]=nc_write(BRYname, 'lon_u_west'   , S.lon_u_west);
    [status]=nc_write(BRYname, 'lat_u_west'   , S.lat_u_west);
    [status]=nc_write(BRYname, 'lon_v_west'   , S.lon_v_west);
    [status]=nc_write(BRYname, 'lat_v_west'   , S.lat_v_west);
  end,
  if (S.boundary(2)),
    [status]=nc_write(BRYname, 'lon_rho_east' , S.lon_rho_east);
    [status]=nc_write(BRYname, 'lat_rho_east' , S.lat_rho_east);
    [status]=nc_write(BRYname, 'lon_u_east'   , S.lon_u_east);
    [status]=nc_write(BRYname, 'lat_u_east'   , S.lat_u_east);
    [status]=nc_write(BRYname, 'lon_v_east'   , S.lon_v_east);
    [status]=nc_write(BRYname, 'lat_v_east'   , S.lat_v_east);
  end,
  if (S.boundary(3)),
    [status]=nc_write(BRYname, 'lon_rho_south', S.lon_rho_south);
    [status]=nc_write(BRYname, 'lat_rho_south', S.lat_rho_south);
    [status]=nc_write(BRYname, 'lon_u_south'  , S.lon_u_south);
    [status]=nc_write(BRYname, 'lat_u_south'  , S.lat_u_south);
    [status]=nc_write(BRYname, 'lon_v_south'  , S.lon_v_south);
    [status]=nc_write(BRYname, 'lat_v_south'  , S.lat_v_south);
  end,
  if (S.boundary(4)),
    [status]=nc_write(BRYname, 'lon_rho_north', S.lon_rho_north);
    [status]=nc_write(BRYname, 'lat_rho_north', S.lat_rho_north);
    [status]=nc_write(BRYname, 'lon_u_north'  , S.lon_u_north);
    [status]=nc_write(BRYname, 'lat_u_north'  , S.lat_u_north);
    [status]=nc_write(BRYname, 'lon_v_north'  , S.lon_v_north);
    [status]=nc_write(BRYname, 'lat_v_north'  , S.lat_v_north);
  end,

%  Initialize unlimited dimension record counter.

  BryRec=0;

end,

%---------------------------------------------------------------------------
%  Interpolate boundary conditions from Mercator data to application grid.
%---------------------------------------------------------------------------

disp(' ');
disp([ 'Interpolating from Mercator to ROMS grid ...']);

for n=1:nfiles,

%  Uncompress input Mercator files.

  if (n <= nfiles1),
    Dir=RTR_Dir;
  else,
    Dir=OPR_Dir;
  end,
  
  fileT=fullfile(Dir,Tfile(n).name);  lenT=length(fileT)-3;
  fileU=fullfile(Dir,Ufile(n).name);  lenU=length(fileU)-3;
  fileV=fullfile(Dir,Vfile(n).name);  lenV=length(fileV)-3;
  
  s=unix(['gunzip ',fileT]);
  s=unix(['gunzip ',fileU]);
  s=unix(['gunzip ',fileV]);

%  Read Mercator data has a time coordinate counter (seconds) that
%  starts on 11-Oct-2006.

  time=nc_read(fileT(1:lenT),'time_counter');
  mydate=datestr(datenum('11-Oct-2006')+time/86400-0.5,0);

  disp(' ');
  disp([ '*** Processing: ', mydate]);
  disp(' ');
  
%  Get Mercator grid.

  Tlon=nc_read(fileT(1:lenT),'nav_lon');
  Tlat=nc_read(fileT(1:lenT),'nav_lat');
  Tdepth=nc_read(fileT(1:lenT),'deptht');

  Ulon=nc_read(fileU(1:lenU),'nav_lon');
  Ulat=nc_read(fileU(1:lenU),'nav_lat');
  Udepth=nc_read(fileU(1:lenU),'depthu');

  Vlon=nc_read(fileV(1:lenV),'nav_lon');
  Vlat=nc_read(fileV(1:lenV),'nav_lat');
  Vdepth=nc_read(fileV(1:lenV),'depthv');

%  In the western Pacific, the level 50 (z=5727.9 m) of the Mercator data
%  is all zeros. Our grid needs depth of Zr=-5920 m.  Therefore, the depths
%  are modified in level 49 (z=5274.7 m) to bound the vertical interpolation.

  Tdepth(49)=6000; Udepth(49)=6000; Vdepth(49)=6000;
  Tdepth(50)=6100; Udepth(50)=6100; Vdepth(50)=6100;

%  Read in initial conditions fields.

  Zeta=nc_read(fileT(1:lenT),'sossheig');
  Temp=nc_read(fileT(1:lenT),'votemper');
  Salt=nc_read(fileT(1:lenT),'vosaline');
  Uvel=nc_read(fileU(1:lenU),'vozocrtx');
  Vvel=nc_read(fileV(1:lenV),'vomecrty');

%  Determine Mercator Land/Sea mask.  Since Mercator is a Z-level
%  model, the mask is 3D.

  Rmask3d=ones(size(Temp));
  ind=find(Temp == 0);
  Rmask3d(ind)=0;

  clear ind

%  Compress input Mercator files.

  s=unix(['gzip ',fileT(1:lenT)]);
  s=unix(['gzip ',fileU(1:lenU)]);
  s=unix(['gzip ',fileV(1:lenV)]);

%  Set boundary conditions time (seconds). The time coordinate for
%  this ROMS application is "seconds since 2007-01-01 00:00:00".
%  coefficient here is to account Mecator daily average. So the data
%  is centered at 12:00:00 hours.

  MyTime=time/86400-(datenum('1-Jan-2007')-datenum('11-Oct-2006'))-0.5;

  if (BryRec == 0),
    S.bry_time = 0;       % to bound initilization on 01-Jan-2007 00:00:00
  else,
    S.bry_time = MyTime*86400;
  end,

%  Interpolate free-surface initial conditions.

  zeta=obc_mercator('zeta',S,Zeta,Tlon,Tlat,Rmask3d(:,:,1));

%  Interpolate temperature and salinity.

  temp=obc_mercator('temp',S,Temp,Tlon,Tlat,Rmask3d,Tdepth);
  salt=obc_mercator('salt',S,Salt,Tlon,Tlat,Rmask3d,Tdepth);
  Urho=obc_mercator('u'   ,S,Uvel,Ulon,Ulat,Rmask3d,Udepth);
  Vrho=obc_mercator('v'   ,S,Vvel,Vlon,Vlat,Rmask3d,Vdepth);

%  Process velocity: rotate and/or average to staggered C-grid locations.

  [u,v]=roms_vectors(Urho,Vrho,S.angle,S.mask_u,S.mask_v,S.boundary);

%  Compute barotropic velocities by vertically integrating (u,v).

  [ubar,vbar]=uv_barotropic(u,v,S.Hz,S.boundary);

  disp(' ');

%--------------------------------------------------------------------------
%  Write out boundary conditions.
%--------------------------------------------------------------------------

  if (WRITE),

    BryRec = BryRec+1;

    [status]=nc_write(BRYname, 'bry_time', S.bry_time, BryRec);

    if (S.boundary(1)),
      [status]=nc_write(BRYname, 'zeta_west' , zeta.west, BryRec);
      [status]=nc_write(BRYname, 'ubar_west' , ubar.west, BryRec);
      [status]=nc_write(BRYname, 'vbar_west' , vbar.west, BryRec);
      [status]=nc_write(BRYname, 'u_west',     u.west,    BryRec);
      [status]=nc_write(BRYname, 'v_west',     v.west,    BryRec);
      [status]=nc_write(BRYname, 'temp_west' , temp.west, BryRec);
      [status]=nc_write(BRYname, 'salt_west' , salt.west, BryRec);
    end,
    if (S.boundary(2)),
      [status]=nc_write(BRYname, 'zeta_east' , zeta.east, BryRec);
      [status]=nc_write(BRYname, 'ubar_east' , ubar.east, BryRec);
      [status]=nc_write(BRYname, 'vbar_east' , vbar.east, BryRec);
      [status]=nc_write(BRYname, 'u_east',     u.east,    BryRec);
      [status]=nc_write(BRYname, 'v_east',     v.east,    BryRec);
      [status]=nc_write(BRYname, 'temp_east' , temp.east, BryRec);
      [status]=nc_write(BRYname, 'salt_east' , salt.east, BryRec);
    end,
    if (S.boundary(3)),
      [status]=nc_write(BRYname, 'zeta_south', zeta.south, BryRec);
      [status]=nc_write(BRYname, 'ubar_south', ubar.south, BryRec);
      [status]=nc_write(BRYname, 'vbar_south', vbar.south, BryRec);
      [status]=nc_write(BRYname, 'u_south',    u.south,    BryRec);
      [status]=nc_write(BRYname, 'v_south',    v.south,    BryRec);
      [status]=nc_write(BRYname, 'temp_south', temp.south, BryRec);
      [status]=nc_write(BRYname, 'salt_south', salt.south, BryRec);
    end,
    if (S.boundary(4)),
      [status]=nc_write(BRYname, 'zeta_north', zeta.north, BryRec);
      [status]=nc_write(BRYname, 'ubar_north', ubar.north, BryRec);
      [status]=nc_write(BRYname, 'vbar_north', vbar.north, BryRec);
      [status]=nc_write(BRYname, 'u_north',    u.north,    BryRec);
      [status]=nc_write(BRYname, 'v_north',    v.north,    BryRec);
      [status]=nc_write(BRYname, 'temp_north', temp.north, BryRec);
      [status]=nc_write(BRYname, 'salt_north', salt.north, BryRec);
    end,

  end,

%  Process next boundary record.

end,
