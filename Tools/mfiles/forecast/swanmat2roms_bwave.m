function swanmat2roms_bwave(iname,fname,plot_key);
% SWAN2ROMS_v3 - Convert SWAN output files to ROMS input netCDF files
% Usage: SWAN2ROMS(iname,fname,plot_key);
% 
% Inputs: iname = SWAN input file name (e.g. 'INPUT')
%         fname = ROMS forcing file name (.e.g. 'swan_frc.nc')
%         plot_key = 1 to create plots, 0 otherwise
% Secondary functions: swan_meta_v3.m, netcdf_toolbox

% Function to convert SWAN output data files of interest to ROMS (depth.dat, hsig.dat,
% rtp.dat, ubot.dat, wdir.dat, xp.dat, and yp.dat) to a netCDF file ('fname').
% If plot_key = 0, the program does not plot the variables; if plot_key = 1, the first
% time of each output variable is plotted in a matlab figure window.

% SWAN2ROMS expects to read 'INPUT' in Layout 4 format.  Snippet:

% ! Model specification
% CGRID CURVILINEAR 159 59 EXC 99999. CIRCLE 36 0.05 0.6 26
% ...
% ! Model output
% GROUP 'group1' SUBG 0 159 0 59
% BLOCK 'group1' NOHEAD 'xp.dat' LAY 4 XP
% BLOCK 'group1' NOHEAD 'yp.dat' LAY 4 YP
% BLOCK 'group1' NOHEAD 'depth.dat'  LAY 4 DEPTH
% BLOCK 'group1' NOHEAD 'hsig.dat'   LAY 4 HSIG   OUTPUT  20020912.000000 3 HR
% BLOCK 'group1' NOHEAD 'wdir.dat'   LAY 4 DIR    OUTPUT  20020912.000000 3 HR
% BLOCK 'group1' NOHEAD 'rtp.dat'    LAY 4 RTP    OUTPUT  20020912.000000 3 HR
% BLOCK 'group1' NOHEAD 'tmbot.dat'  LAY 4 TMBOT  OUTPUT  20020912.000000 3 HR
% BLOCK 'group1' NOHEAD 'ubot.dat'   LAY 4 UBOT   OUTPUT  20020912.000000 3 HR
% BLOCK 'group1' NOHEAD 'wind.dat'   LAY 4 WIND   OUTPUT  20020912.000000 3 HR
% BLOCK 'group1' NOHEAD 'wlen.dat'   LAY 4 WLEN   OUTPUT  20020912.000000 3 HR

% ! Wave-current interaction output
% BLOCK 'group1' NOHEAD 'dissip.dat' LAY 4 DISSIP OUTPUT  20020912.000000 3 HR
% BLOCK 'group1' NOHEAD 'force.dat'  LAY 4 FORCE  OUTPUT  20020912.000000 3 HR
% BLOCK 'group1' NOHEAD 'vel.dat'    LAY 4 VEL    OUTPUT  20020912.000000 3 HR
% BLOCK 'group1' NOHEAD 'setup.dat'  LAY 4 SETUP  OUTPUT  20020912.000000 3 HR

% ! Wave-current interaction output
% BLOCK 'COMPGRID' NOHEADER 'dissip.dat' LAY 4 DISSIP 1.
% BLOCK 'COMPGRID' NOHEADER 'force.dat' LAY 4 FORCE 1.
% BLOCK 'COMPGRID' NOHEADER 'setup.dat' LAY 4 SETUP 1.
% BLOCK 'COMPGRID' NOHEADER 'vel.dat' LAY 4 VEL 1.

%Soupy Alexander, 6/13/02
%jcwarner: added internal netcdf creation
%csherwood@usgs.gov: changed documentation a little, and made
%      wave-current interaction terms optional
% rsignell@usgs.gov: Rich Signell : rewrote to use "textread" instead of "fgetl" (25x faster)
%                                 : writes correct time vector in Modified Julian Day 
%                                 : based on info from INPUT 
%
%jcwarner 30Sep2005 - fix a few small plotting bug fixes in var names. Calc nt only once.


% these are required
write_hsig = 1;
write_xp = 1;
write_yp = 1;
write_depth = 0;

% these are standard
% existence of these files is not tested
write_rtp = 0;
write_wdir = 1;
write_ubot = 0;
write_tmbot = 0;

% these are optional and only written if the corresponding SWAN output file
% exists
write_wlen = 0;
write_wind = 0;
write_dissip = 0;
write_force = 0;
write_setup = 0;
write_vel = 0;
write_break = 0;

%jcw enter file names here
hsig_file  = 'hsig.mat';
depth_file ='depth.mat';
dissip_file= 'dissip.mat';
force_file ='force.mat';
setup_file ='setup.mat';
rtp_file   =  'rtp.mat';
tmbot_file ='tmbot.mat';
ubot_file  = 'ubot.mat';
vel_file   =  'vel.mat';
wdir_file  = 'wdir.mat';
wlen_file  = 'wlen.mat';
wind_file  = 'wind.mat';
break_file  = 'qb.mat';

%if(exist('dissip.mat','file')==2),write_dissip=1;,end
%if(exist('hsig.mat','file')==2),write_hsig=1;,end
%if(exist('rtp.mat','file')==2), write_rtp=1; end
%if(exist('tmbot.mat','file')==2),write_tmbot=1; end
%if(exist('ubot.mat','file')==2),write_ubot=1;,end
%if(exist('wdir.mat','file')==2),write_wdir=1;,end
%if(exist('wlen.mat','file')==2),write_wlen=1;,end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read meta data from SWAN input file
s=swan_meta_v3(iname);
xsize=s.xsize;
ysize=s.ysize;
jdmat0=s.jdmat0;
dt=s.dt;    % seconds
nt=s.nt;    % # of records each file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create netcdf file.
nc=netcdf(fname,'clobber');
if isempty(nc), return, end

%% Global attributes:

disp(' ## Defining Global Attributes...')
nc.history = ncchar(['Created by ROMS2SWAN_v3 on ' datestr(now)]);
nc.type = ncchar('forcing file from swan output via swan2roms.m');

%% Dimensions:

disp(' ## Defining Dimensions...')

LP=xsize;
MP=ysize;
L=LP-1;
M=MP-1;
N=0;

nc('xi_psi') = L;
nc('xi_rho') = LP;
nc('xi_u') = L;
nc('xi_v') = LP;

nc('eta_psi') = M;
nc('eta_rho') = MP;
nc('eta_u') = MP;
nc('eta_v') = M;

nc('wave_time')=nt;
nc('one') = 1;

%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')

nc{'wave_time'} = ncdouble('wave_time'); %% 1 element.
nc{'wave_time'}.long_name = ncchar('wave field time');
nc{'wave_time'}.units = ncchar('Julian day');
nc{'wave_time'}.field = ncchar('wave_time, scalar, series');

if(write_depth)
   nc{'swan_depth'} = ncfloat('wave_time','eta_rho', 'xi_rho');
   nc{'swan_depth'}.long_name = ncchar('swan_depth');
   nc{'swan_depth'}.units = ncchar('meter');
   nc{'swan_depth'}.field = ncchar('depth, scalar, series');
end
if(write_dissip)
   nc{'Wave_dissip'} = ncfloat('wave_time', 'eta_rho', 'xi_rho');
   nc{'Wave_dissip'}.long_name = ncchar('Wave_dissip');
   nc{'Wave_dissip'}.units = ncchar('Watts meter-2');
   nc{'Wave_dissip'}.field = ncchar('Wave_dissip, scalar, series');
end
if(write_force)
   nc{'Wave_forcex'} = ncfloat('wave_time', 'eta_rho', 'xi_rho');
   nc{'Wave_forcex'}.long_name = ncchar('Wave_forcex');
   nc{'Wave_forcex'}.units = ncchar('Newtons meter-2');
   nc{'Wave_forcex'}.field = ncchar('Wave_forcex, scalar, series');

   nc{'Wave_forcey'} = ncfloat('wave_time', 'eta_rho', 'xi_rho');
   nc{'Wave_forcey'}.long_name = ncchar('Wave_forcey');
   nc{'Wave_forcey'}.units = ncchar('Newtons meter-2');
   nc{'Wave_forcey'}.field = ncchar('Wave_forcey, scalar, series');
end
if(write_hsig)
   nc{'Hwave'} = ncfloat('wave_time', 'eta_rho', 'xi_rho');
   nc{'Hwave'}.long_name = ncchar('Hwave');
   nc{'Hwave'}.units = ncchar('meter');
   nc{'Hwave'}.field = ncchar('Hwave, scalar, series');
end
if(write_setup),
   nc{'Wave_setup'} = ncfloat('wave_time', 'eta_rho', 'xi_rho');
   nc{'Wave_setup'}.long_name = ncchar('Wave_setup');
   nc{'Wave_setup'}.units = ncchar('meter');
   nc{'Wave_setup'}.field = ncchar('Wave_setup, scalar, series');
end
if(write_rtp)
   nc{'Pwave_top'} = ncfloat('wave_time', 'eta_rho', 'xi_rho');
   nc{'Pwave_top'}.long_name = ncchar('Pwave_top');
   nc{'Pwave_top'}.units = ncchar('second');
   nc{'Pwave_top'}.field = ncchar('Pwave_top, scalar, series');
end
if(write_tmbot)
   nc{'Pwave_bot'} = ncfloat('wave_time', 'eta_rho', 'xi_rho');
   nc{'Pwave_bot'}.long_name = ncchar('Pwave_bot');
   nc{'Pwave_bot'}.units = ncchar('second');
   nc{'Pwave_bot'}.field = ncchar('Pwave_bot, scalar, series');
end
if(write_ubot)
   nc{'Uwave_rms'} = ncfloat('wave_time', 'eta_rho', 'xi_rho');
   nc{'Uwave_rms'}.long_name = ncchar('Uwave_rms');
   nc{'Uwave_rms'}.units = ncchar('meter second-1');
   nc{'Uwave_rms'}.field = ncchar('Uwave_rms, scalar, series');
end
if(write_vel)
   nc{'velx'} = ncfloat('wave_time', 'eta_rho', 'xi_rho');
   nc{'velx'}.long_name = ncchar('velx');
   nc{'velx'}.units = ncchar('meter second-1');
   nc{'velx'}.field = ncchar('velx, scalar, series');

   nc{'vely'} = ncfloat('wave_time', 'eta_rho', 'xi_rho');
   nc{'vely'}.long_name = ncchar('vely');
   nc{'vely'}.units = ncchar('meter second-1');
   nc{'vely'}.field = ncchar('vely, scalar, series');
end
if(write_wdir)
   nc{'Dwave'} = ncfloat('wave_time', 'eta_rho', 'xi_rho');
   nc{'Dwave'}.long_name = ncchar('Dwave');
   nc{'Dwave'}.units = ncchar('degrees');
   nc{'Dwave'}.field = ncchar('Dwave, scalar, series');
end
if(write_wlen)
   nc{'Lwave'} = ncfloat('wave_time', 'eta_rho', 'xi_rho');
   nc{'Lwave'}.long_name = ncchar('Lwave');
   nc{'Lwave'}.units = ncchar('meter');
   nc{'Lwave'}.field = ncchar('Lwave, scalar, series');
end
if(write_break)
   nc{'Wave_break'} = ncfloat('wave_time', 'eta_rho', 'xi_rho');
   nc{'Wave_break'}.long_name = ncchar('Wave_break');
   nc{'Wave_break'}.units = ncchar('meter');
   nc{'Wave_break'}.field = ncchar('Wave_break, scalar, series');
end
if(write_xp)
   nc{'xp'} = ncfloat('eta_rho', 'xi_rho');
   nc{'xp'}.long_name = ncchar('xp');
   nc{'xp'}.units = ncchar('meter');
   nc{'xp'}.FillValue_ = ncfloat(100000.);
   nc{'xp'}.missing_value = ncfloat(100000.);
   nc{'xp'}.field = ncchar('xp, scalar, series');
end
if(write_yp)
   nc{'yp'} = ncfloat('eta_rho', 'xi_rho');
   nc{'yp'}.long_name = ncchar('yp');
   nc{'yp'}.units = ncchar('meter');
   nc{'yp'}.missing_value = ncfloat(100000.);
   nc{'yp'}.FillValue_ = ncfloat(100000.);
   nc{'yp'}.field = ncchar('yp, scalar, series');
end
if(write_wind)
   nc{'Windx'} = ncfloat('wave_time', 'eta_rho', 'xi_rho');
   nc{'Windx'}.long_name = ncchar('Windx');
   nc{'Windx'}.units = ncchar('meter second-1');
   nc{'Windx'}.field = ncchar('Windx, scalar, series');

   nc{'Windy'} = ncfloat('wave_time', 'eta_rho', 'xi_rho');
   nc{'Windy'}.long_name = ncchar('Windy');
   nc{'Windy'}.units = ncchar('meter second-1');
   nc{'Windy'}.field = ncchar('Windy, scalar, series');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now write the data from the arrays to the netcdf file

if(write_hsig)
% already read in above
   disp('   hsig -> Hwave')
%   F=load('hsig.mat');
   F=load(hsig_file);
   vname=fieldnames(F);
    for K=1:length(vname)
        disp(sprintf('%g/%g Hwave saved',K,length(vname)))
        
       Hsig=getfield(F,char(vname(K)));
       Hsig(find(Hsig == -99)) = 0.0;   % flag
       Hsig(find(isnan(Hsig)))=0.0;     % land
       nc{'Hwave'}(K,:,:)=squeeze(Hsig);
       clear Hsig;
   end
   clear F vname K
end

if(write_depth)
   disp('   depth -> swan_depth')
%   F=load('depth.mat');
   F=load(depth_file);
   vname=fieldnames(F);
   Depth=[];
   for K=1:length(vname)
        disp(sprintf('%g/%g swan_depth saved',K,length(vname)))
        
       Depth=getfield(F,char(vname(K)));
       Depth(find(Depth == -99)) = 0.0;     % flag
       Depth(find(isnan(Depth)))=0.0;       % land    
       nc{'Hwave'}(K,:,:)=squeeze(Depth);
       clear Depth;
   end
   clear F vname K
end

if(write_dissip)
   disp('   dissip -> Wave_dissip')
%   F=load('dissip.mat');
   F=load(dissip_file);
   vname=fieldnames(F);
   Dissip=[];
   for K=1:length(vname)
        disp(sprintf('%g/%g Dissip saved',K,length(vname)))
        
       Dissip=getfield(F,char(vname(K)));
       Dissip(find(Dissip == -9)) = 0.0;    % flag
       Dissp(find(isnan(Dissip))) = 0.0;    % land    
       nc{'Wave_dissip'}(K,:,:) = squeeze(Dissip);
       clear Dissip;
   end
   clear F vname K 
end

if(write_force)
   disp('   force => Wave_forcex, Wave_forcey')
%   F=load('force.mat');
   F=load(force_file);
   vname=fieldnames(F);
   for K=1:2:length(vname)-1
        disp(sprintf('%g/%g Force saved',K,length(vname)))
        
       Forcex=getfield(F,char(vname(K)));
       Forcey=getfield(F,char(vname(K+1)));
       Forcex(find(Forcex == 0)) = 0.0;     % flag
       Forcey(find(Forcey == 0)) = 0.0;     % flag
       Forcex(find(isnan(Forcex))) = 0.0;   % land
       Forcey(find(isnan(Forcey))) = 0.0;   % land
       nc{'Wave_forcex'}(K,:,:) = squeeze(Forcex);
       nc{'Wave_forcey'}(K,:,:) = squeeze(Forcey);
       clear Force*;
   end
   clear F vname K 
end
if(write_setup)
   disp('   setup -> Wave_setup')
%   F=load('setup.mat');
   F=load(setup_file);
   vname=fieldnames(F);
   Setup=[];
   for K=1:length(vname)
        disp(sprintf('%g/%g ''Setup'' saved',K,length(vname)))
        
       Setup=getfield(F,char(vname(K)));
       Setup(find(Setup == -9)) = 0.0;      % flag
       Setup(find(isnan(Setup))) = 0.0;     % land
       nc{'Wave_setup'}(K,:,:) = squeeze(Setup);
       clear Setup;
   end
   clear F vname K 
end
if(write_rtp)
   disp('   rtp - > Pwave_top')
%   F=load('rtp.mat');
   F=load(rtp_file);
   vname=fieldnames(F);
   for K=1:length(vname)
        disp(sprintf('%g/%g ''Pwave_top'' saved',K,length(vname)))
        
       Rtp=getfield(F,char(vname(K)));
       Rtp(find(Rtp == -9)) = 10.0;     % flag
       Rtp(find(isnan(Rtp))) = 0.0;     % land
       nc{'Pwave_top'}(K,:,:) = squeeze(Rtp);
       clear Rtp;
   end
   clear F vname K 
end
if(write_tmbot)
   disp('   tmbot - > Pwave_bot')
%   F=load('tmbot.mat');
   F=load(tmbot_file);
   vname=fieldnames(F);
   for K=1:length(vname)
        disp(sprintf('%g/%g ''Pwave_bottom'' saved',K,length(vname)))
        
       Tmbot=getfield(F,char(vname(K)));
       Tmbot(find(Tmbot == -9)) = 10.0;     % flag
       Tmbot(find(isnan(Tmbot))) = 0.0;     % land
       nc{'Pwave_bot'}(K,:,:) = squeeze(Tmbot);
       clear Tmbot;
   end
   clear F vname K 
end
if(write_ubot)
   disp('   ubot -> Uwave_rms')
%   F=load('ubot.mat');
   F=load(ubot_file);
   vname=fieldnames(F);
   for K=1:length(vname)
        disp(sprintf('%g/%g ''Uwave_rms'' saved',K,length(vname)))
        
       Ubot=getfield(F,char(vname(K)));
       Ubot(find(Ubot == -10)) = 0.0001;    % flag
       Ubot(find(isnan(Ubot))) = 0.0;       % land
       nc{'Uwave_rms'}(K,:,:) = squeeze(Ubot);
       clear Ubot;
   end
   clear F vname K 
end

if(write_vel)
   disp('   vel -> velx, vely')
%   F=load('vel.mat');
   F=load(vel_file);
   vname=fieldnames(F);
   for K=1:2:length(vname)-1
        disp(sprintf('%g/%g ''Velocity'' saved',K,length(vname)))
        
       Velx=getfield(F,char(vname(K)));
       Vely=getfield(F,char(vname(K+1)));
       Velx(find(Velx == 0)) = 0.0;         % flag
       Vely(find(Vely == 0)) = 0.0;         % flag
       Velx(find(isnan(Velx))) = 0.0;       % land
       Vely(find(isnan(Vely))) = 0.0;       % land
       nc{'velx'}(K,:,:) = squeeze(Velx);
       nc{'vely'}(K,:,:) = squeeze(Vely);
       clear Vel*;
   end
   clear F vname K 
end
if(write_wdir)
   disp('   wdir -> Dwave')
%   F=load('wdir.mat');
   F=load(wdir_file);
   vname=fieldnames(F);
   for K=1:length(vname)
        disp(sprintf('%g/%g ''Dwave'' saved',K,length(vname)))
        
       Wdir=getfield(F,char(vname(K)));
       Wdir(find(Wdir == -999)) = 0.0;  % flag
       Wdir(find(isnan(Wdir))) = 0.0;   % land
       nc{'Dwave'}(K,:,:) = squeeze(Wdir);
       clear Wdir;
   end
   clear F vname K 
end

if(write_wlen)
   disp('   wlen -> Lwave')
%   F=load('wlen.mat');
   F=load(wlen_file);
   vname=fieldnames(F);
   for K=1:length(vname)
        disp(sprintf('%g/%g ''Lwave'' saved',K,length(vname)))
        
       Wlen=getfield(F,char(vname(K)));
       Wlen(find(Wlen == -9)) = 10.0;   % flag
       Wlen(find(isnan(Wlen))) = 0.0;   % land
       nc{'Lwave'}(K,:,:) = squeeze(Wlen);
       clear Wlen;
   end
   clear F vname K 
end

if(write_break)
   disp('   qb -> Wave break')
%   F=load('wlen.mat');
   F=load(break_file);
   vname=fieldnames(F);
   for K=1:length(vname)
        disp(sprintf('%g/%g ''Qb'' saved',K,length(vname)))
        
       Qb=getfield(F,char(vname(K)));
       Qb(find(Qb == -9)) = 10.0;   % flag
       Qb(find(isnan(Qb))) = 0.0;   % land
       nc{'Wave_break'}(K,:,:) = squeeze(Qb);
       clear Qb;
   end
   clear F vname K 
end

if(write_xp)
   disp('   xp')
   F=load('xp.mat');
   vname=fieldnames(F);
   for K=1:1 %length(vname)
        disp(sprintf('%g/%g ''Xp'' saved',K,length(vname)))
        
        Xp=getfield(F,char(vname(K)));
        nc{'xp'}(:,:) = squeeze(Xp);
       clear Xp;
   end
   clear F vname K 
end

if(write_yp)
   disp('   yp')
   F=load('yp.mat');
   vname=fieldnames(F);
   for K=1:1 %length(vname)
        disp(sprintf('%g/%g ''Yp'' saved',K,length(vname)))
        
        Yp=getfield(F,char(vname(K)));
        nc{'yp'}(:,:) = squeeze(Yp);
        clear Yp;
   end
   clear F vname K 
end

if(write_wind)
   disp('   wind -> Windx, Windy')
%   F=load('wind.mat');
   F=load(wind_file);
   vname=fieldnames(F);
   count=1;
   for K=1:2:length(vname)-1
       disp(sprintf('%g/%g ''Wind'' saved',K,length(vname)))
       Wdirx=getfield(F,char(vname(K)));
       Wdiry=getfield(F,char(vname(K+1)));
       Wdirx(find(Wdirx == -9)) = 0.0;      % flag
       Wdiry(find(Wdiry == -9)) = 0.0;      % flag
       Wdirx(find(isnan(Wdirx))) = 0.0;     % land
       Wdiry(find(isnan(Wdiry))) = 0.0;     % land
       nc{'Windx'}(count,:,:) = squeeze(Wdirx);
       nc{'Windy'}(count,:,:) = squeeze(Wdiry);
       clear Wdir*;
       count=count+1;
   end
   clear F vname K 
end


%jdshift=datenum(1968,5,23,0,0,0);
mjd_days=jdmat0+[0:(nt-1)]*(dt/(24*3600));%-jdshift;  %time in datenum
% write time
nc{'wave_time'}(:) = mjd_days;   % Modified Julian Days = (Julian Day - 2440000.) = days since 1968-05-23 00:00 UTC

close(nc)

%Plot the values, if this was requested - in the plots, defaults are replaced with NaN's.
if exist('plot_key')
if plot_key >= 1;

   if nt == 1
      ncload(fname)
   else
      tidx =plot_key;
      nc=netcdf(fname);
      depth=squeeze(nc{'depth'}(:,:));
      Wave_dissip=squeeze(nc{'Wave_dissip'}(tidx,:,:));
      Wave_forcex=squeeze(nc{'Wave_forcex'}(tidx,:,:));
      Wave_forcey=squeeze(nc{'Wave_forcey'}(tidx,:,:));
      Hwave=squeeze(nc{'Hwave'}(tidx,:,:));
      Wave_setup=squeeze(nc{'Wave_setup'}(tidx,:,:));
      Pwave_top=squeeze(nc{'Pwave_top'}(tidx,:,:));
      Pwave_bot=squeeze(nc{'Pwave_bot'}(tidx,:,:));
      Uwave_rms=squeeze(nc{'Uwave_rms'}(tidx,:,:));
      velx=squeeze(nc{'velx'}(tidx,:,:));
      vely=squeeze(nc{'vely'}(tidx,:,:));
      Dwave=squeeze(nc{'Dwave'}(tidx,:,:));
      Lwave=squeeze(nc{'Lwave'}(tidx,:,:));
      Windx=squeeze(nc{'Windx'}(tidx,:,:));
      Windy=squeeze(nc{'Windy'}(tidx,:,:));
      xp=nc{'xp'}(:);
      yp=nc{'yp'}(:);
   end

   % load (x,y) grid points
   F=load('xp.mat'); vname=fieldnames(F);
   xp=getfield(F,char(vname(1)));
   clear F vname
   
   F=load('yp.mat'); vname=fieldnames(F);
   up=getfield(F,char(vname(1)));
   clear F vname
   
   % plotting 
   if write_hsig
   figure;
   %   hsig(find(hsig == 0)) = NaN;
   pcolor(xp,yp,squeeze(Hwave(:,:))); shading interp; colorbar;
   title('Hs Output of Swan');
   end
   
   if write_ubot
   figure;
   %    Uwave_rms(find(Uwave_rms == 0.001)) = NaN;
   pcolor(xp,yp,squeeze(Uwave_rms(:,:))); shading interp; colorbar;
   title('Ub Output of Swan');
   end
   
   if write_tmbot
   figure;
   Pwave_bot(find(Pwave_bot == 9999)) = NaN;
   pcolor(xp,yp,squeeze(Pwave_bot(:,:))); shading interp;
   title('Pwave_bot Output of Swan');
   caxis([2 15]);colorbar
   end
   
   if write_wdir
   figure;
   %    Dwave(find(Dwave == 0)) = NaN;
   pcolor(xp,yp,squeeze(Dwave(:,:))); shading interp; colorbar;
   title('Dwave Output of Swan');
   end

   if write_setup 
   figure;
   %    Dwave(find(Dwave == 0)) = NaN;
   pcolor(xp,yp,squeeze(Wave_setup(:,:))); shading interp; colorbar;
   title('setup Output of Swan');
   end

   if write_wlen 
   figure;
   %    Dwave(find(Dwave == 0)) = NaN;
   pcolor(xp,yp,squeeze(Lwave(:,:))); shading interp;
   title('Lwave Output of Swan');
   caxis([0 75]);colorbar
   end
   
   if write_dissip
   figure
   %    Dwave(find(Dwave == 0)) = NaN;
   pcolor(xp,yp,squeeze(Wave_dissip(:,:))); shading interp; colorbar;
   title('Wave dissipation')
   end

   if write_force 
   figure
   subplot(121)
   %    Dwave(find(Dwave == 0)) = NaN;
   pcolor(xp,yp,squeeze(Wave_forcex(:,:))); shading interp; colorbar;
   title('forcex Output of Swan');

   subplot(122)
   %    Dwave(find(Dwave == 0)) = NaN;
   pcolor(xp,yp,squeeze(Wave_forcey(:,:))); shading interp; colorbar;
   title('forcey Output of Swan');
   end

   if write_force 
   figure
   subplot(121)
%    vel(find(vel == 0)) = NaN;
   pcolor(xp,yp,squeeze(velx(:,:))); shading interp; colorbar;
   title('velx Output of Swan');
   
   subplot(122)
   %    vel(find(vel == 0)) = NaN;
   pcolor(xp,yp,squeeze(vely(:,:))); shading interp; colorbar;
   title('vely Output of Swan');
   end

   if write_wind
   figure
   subplot(121)
%    vel(find(vel == 0)) = NaN;
   pcolor(xp,yp,squeeze(Windx(:,:))); shading interp; colorbar;
   title('wind x Output of Swan');
   
   subplot(122)
   %    vel(find(vel == 0)) = NaN;
   pcolor(xp,yp,squeeze(Windy(:,:))); shading interp; colorbar;
   title('wind y Output of Swan');
   end
end
end
