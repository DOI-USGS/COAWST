function swanmat2roms(fname);
% swanmat2roms - Convert SWAN output matlab files to ROMS input netCDF files
% Usage: swanmat2roms(fname);
% 
% Inputs: fname = ROMS forcing file name (.e.g. 'swan_frc.nc')
%
% Function to convert SWAN output data files of interest to a ROMS
%          netCDF forcing file ('fname').
%
% jcwarner 11Sept2013 and based on several earlier versions
%
disp('  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('                                                         ')
disp(' PLEASE EDIT swanmat2roms AND LIST NAMES of *.mat FILES .')
disp('                                                         ')
disp('  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

% User needs to specify the SWAN Output data file names.
% Place a '1' for the files to get, or a '0' to not include them.
 get_hsig=1;         hsig_file   = 'hsig.mat';
 get_dissip=1;       dissip_file = 'dissip.mat';
 get_rtp=1;          rtp_file    = 'rtp.mat';
 get_tmbot=1;        tmbot_file  = 'tmbot.mat';
 get_ubot=1;         ubot_file   = 'ubot.mat';
 get_wdir=1;         wdir_file   = 'wdir.mat';
 get_wlen=1;         wlen_file   = 'wlen.mat';
 get_break=1;        break_file  = 'qb.mat';
 get_xp=1;           xp_file     = 'xp.mat';
 get_yp=1;           yp_file     = 'yp.mat';
 get_dissip_break=0; dissip_break_file= 'dissip_break.mat';
 get_dissip_wcap=0;  dissip_wcap_file= 'dissip_wcap.mat';
 get_dissip_fric=0;  dissip_fric_file= 'dissip_fric.mat';

% See Kumar et al 2012 OM 47, 65-95.
% here are the files you need for each cpp:
% WEC_MELLOR           hsig, wdir, wlen
% WEC_VF               hsig, wdir, wlen
%
% need to select method to compute dissip_break and dissip_wcap
% WDISS_THORGUZA       hsig, rtp
% WDISS_CHURTHOR       hsig, rtp
% if you do not select either of those 2, then you need to provide:
% dissip_break and dissip_wcap
%
% ROLLER_SVENDSEN      hsig, wdir, wlen, break, dissip_break
% ROLLER_MONO          hsig, wdir, wlen, break, dissip_break
% ROLLER_RENIERS       hsig, wdir, wlen, dissip_break
%
% BOTTOM_STREAMING     hsig, wdir, wlen, dissip_fric
% SURFACE_STREAMING    hsig, wdir, wlen
%
% WAVE_MIXING          dissip_break and dissip_wcap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load hsig file to get grid size and number of time steps.
F=load(hsig_file);
vname=fieldnames(F);
Hsig=getfield(F,char(vname(1)));
[ysize,xsize]=size(Hsig)
numsteps=length(vname);
%get times
for ii=1:length(vname)
  date=vname{ii};
  yr(ii)=str2num(date(6:9));
  mo(ii)=str2num(date(10:11));
  dy(ii)=str2num(date(12:13));
  hr(ii)=str2num(date(15:16));
  mm(ii)=str2num(date(17:18));
  ss(ii)=str2num(date(19:20));
end
clear F Hsig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create netcdf file.
nc_swan=netcdf.create(fname,'clobber');
if isempty(nc_swan), return, end

%% Global attributes:
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_swan,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by swanmat2roms on ' datestr(now)]);
netcdf.putAtt(nc_swan,netcdf.getConstant('NC_GLOBAL'),'type', 'roms forcing file of swan data');

%% Dimensions:
disp(' ## Defining Dimensions...')
xrhodimID = netcdf.defDim(nc_swan,'xrho',xsize);
erhodimID = netcdf.defDim(nc_swan,'erho',ysize);
wtdimID = netcdf.defDim(nc_swan,'wave_time',numsteps);

%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')

%define time array
wtID = netcdf.defVar(nc_swan,'wave_time','double',wtdimID);
netcdf.putAtt(nc_swan,wtID,'long_name','wave field time');
netcdf.putAtt(nc_swan,wtID,'units','days');
netcdf.putAtt(nc_swan,wtID,'field','wave_time, scalar, series');

if(get_hsig)
  ID = netcdf.defVar(nc_swan,'Hwave','double',[xrhodimID erhodimID wtdimID]);
  netcdf.putAtt(nc_swan,ID,'long_name','Hwave');
  netcdf.putAtt(nc_swan,ID,'units','meter');
  netcdf.putAtt(nc_swan,ID,'coordinates','xp yp');
  netcdf.putAtt(nc_swan,ID,'field','Hwave, scalar, series');
end
if(get_dissip)
  ID = netcdf.defVar(nc_swan,'Wave_dissip','double',[xrhodimID erhodimID wtdimID]);
  netcdf.putAtt(nc_swan,ID,'long_name','Wave_dissip');
  netcdf.putAtt(nc_swan,ID,'units','Watts meter-2');
  netcdf.putAtt(nc_swan,ID,'coordinates','xp yp');
  netcdf.putAtt(nc_swan,ID,'field','Wave_dissip, scalar, series');
end
if(get_rtp)
  ID = netcdf.defVar(nc_swan,'Pwave_top','double',[xrhodimID erhodimID wtdimID]);
  netcdf.putAtt(nc_swan,ID,'long_name','Pwave_top');
  netcdf.putAtt(nc_swan,ID,'units','second');
  netcdf.putAtt(nc_swan,ID,'coordinates','xp yp');
  netcdf.putAtt(nc_swan,ID,'field','Pwave_top, scalar, series');
end
if(get_tmbot)
  ID = netcdf.defVar(nc_swan,'Pwave_bot','double',[xrhodimID erhodimID wtdimID]);
  netcdf.putAtt(nc_swan,ID,'long_name','Pwave_bot');
  netcdf.putAtt(nc_swan,ID,'units','second');
  netcdf.putAtt(nc_swan,ID,'coordinates','xp yp');
  netcdf.putAtt(nc_swan,ID,'field','Pwave_bot, scalar, series');
end
if(get_ubot)
  ID = netcdf.defVar(nc_swan,'Uwave_rms','double',[xrhodimID erhodimID wtdimID]);
  netcdf.putAtt(nc_swan,ID,'long_name','Uwave_rms');
  netcdf.putAtt(nc_swan,ID,'units','meter second-1');
  netcdf.putAtt(nc_swan,ID,'coordinates','xp yp');
  netcdf.putAtt(nc_swan,ID,'field','Uwave_rms, scalar, series');
end
if(get_wdir)
  ID = netcdf.defVar(nc_swan,'Dwave','double',[xrhodimID erhodimID wtdimID]);
  netcdf.putAtt(nc_swan,ID,'long_name','Dwave');
  netcdf.putAtt(nc_swan,ID,'units','degrees');
  netcdf.putAtt(nc_swan,ID,'coordinates','xp yp');
  netcdf.putAtt(nc_swan,ID,'field','Dwave, scalar, series');
end
if(get_wlen)
  ID = netcdf.defVar(nc_swan,'Lwave','double',[xrhodimID erhodimID wtdimID]);
  netcdf.putAtt(nc_swan,ID,'long_name','Lwave');
  netcdf.putAtt(nc_swan,ID,'units','meter');
  netcdf.putAtt(nc_swan,ID,'coordinates','xp yp');
  netcdf.putAtt(nc_swan,ID,'field','Lwave, scalar, series');
end
if(get_break)
  ID = netcdf.defVar(nc_swan,'Wave_break','double',[xrhodimID erhodimID wtdimID]);
  netcdf.putAtt(nc_swan,ID,'long_name','Wave_break');
  netcdf.putAtt(nc_swan,ID,'units','decimal percent');
  netcdf.putAtt(nc_swan,ID,'coordinates','xp yp');
  netcdf.putAtt(nc_swan,ID,'field','Wave_break, scalar, series');
end
if(get_xp)
  ID = netcdf.defVar(nc_swan,'xp','double',[xrhodimID erhodimID]);
  netcdf.putAtt(nc_swan,ID,'long_name','xp');
  netcdf.putAtt(nc_swan,ID,'units','meter');
  netcdf.putAtt(nc_swan,ID,'coordinates','xp yp');
  netcdf.putAtt(nc_swan,ID,'field','xp, scalar, series');
% netcdf.putAtt(nc_swan,ID,'FillValue_',ncfloat(100000.));
% netcdf.putAtt(nc_swan,ID,'missing_value',ncfloat(100000.));
end
if(get_yp)
  ID = netcdf.defVar(nc_swan,'yp','double',[xrhodimID erhodimID]);
  netcdf.putAtt(nc_swan,ID,'long_name','yp');
  netcdf.putAtt(nc_swan,ID,'units','meter');
  netcdf.putAtt(nc_swan,ID,'coordinates','xp yp');
  netcdf.putAtt(nc_swan,ID,'field','yp, scalar, series');
% netcdf.putAtt(nc_swan,ID,'FillValue_',ncfloat(100000.));
% netcdf.putAtt(nc_swan,ID,'missing_value',ncfloat(100000.));
end
netcdf.close(nc_swan)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now write the data from the arrays to the netcdf file
if(get_hsig)
   disp('   hsig -> Hwave')
   F=load(hsig_file);
   vname=fieldnames(F);
   for K=1:length(vname)
       disp(sprintf('%g/%g Hwave saved',K,length(vname)))
       Hsig=getfield(F,char(vname(K)));
       Hsig(find(Hsig == -99)) = 0.0;   % flag
       Hsig(find(isnan(Hsig)))=0.0;     % land
       ncwrite(fname,'Hwave',squeeze(Hsig).',[1 1 K])
       clear Hsig;
   end
   clear F vname K
end

if(get_dissip)
   disp('   dissip -> Wave_dissip')
   F=load(dissip_file);
   vname=fieldnames(F);
   for K=1:length(vname)
       disp(sprintf('%g/%g Dissip saved',K,length(vname)))
       Dissip=getfield(F,char(vname(K)));
       Dissip(find(Dissip == -9)) = 0.0;    % flag
       Dissip(find(isnan(Dissip))) = 0.0;   % land
       ncwrite(fname,'Wave_dissip',squeeze(Dissip).',[1 1 K])
       clear Dissip;
   end
   clear F vname K 
end

if(get_rtp)
   disp('   rtp - > Pwave_top')
   F=load(rtp_file);
   vname=fieldnames(F);
   for K=1:length(vname)
       disp(sprintf('%g/%g ''Pwave_top'' saved',K,length(vname)))
       Rtp=getfield(F,char(vname(K)));
       Rtp(find(Rtp == -9)) = 10.0;     % flag
       Rtp(find(isnan(Rtp))) = 0.0;     % land
       ncwrite(fname,'Pwave_top',squeeze(Rtp).',[1 1 K])
       clear Rtp;
   end
   clear F vname K 
end

if(get_tmbot)
   disp('   tmbot - > Pwave_bot')
   F=load(tmbot_file);
   vname=fieldnames(F);
   for K=1:length(vname)
       disp(sprintf('%g/%g ''Pwave_bottom'' saved',K,length(vname)))
       Tmbot=getfield(F,char(vname(K)));
       Tmbot(find(Tmbot == -9)) = 10.0;     % flag
       Tmbot(find(isnan(Tmbot))) = 0.0;     % land
       ncwrite(fname,'Pwave_bot',squeeze(Tmbot).',[1 1 K])
       clear Tmbot;
   end
   clear F vname K 
end

if(get_ubot)
   disp('   ubot -> Uwave_rms')
   F=load(ubot_file);
   vname=fieldnames(F);
   for K=1:length(vname)
       disp(sprintf('%g/%g ''Uwave_rms'' saved',K,length(vname)))
       Ubot=getfield(F,char(vname(K)));
       Ubot(find(Ubot == -10)) = 0.0001;    % flag
       Ubot(find(isnan(Ubot))) = 0.0;       % land
       ncwrite(fname,'Uwave_rms',squeeze(Ubot).',[1 1 K])
       clear Ubot;
   end
   clear F vname K 
end

if(get_wdir)
   disp('   wdir -> Dwave')
   F=load(wdir_file);
   vname=fieldnames(F);
   for K=1:length(vname)
       disp(sprintf('%g/%g ''Dwave'' saved',K,length(vname)))
       Wdir=getfield(F,char(vname(K)));
       Wdir(find(Wdir == -999)) = 0.0;  % flag
       Wdir(find(isnan(Wdir))) = 0.0;   % land
       ncwrite(fname,'Dwave',squeeze(Wdir).',[1 1 K])
       clear Wdir;
   end
   clear F vname K 
end

if(get_wlen)
   disp('   wlen -> Lwave')
   F=load(wlen_file);
   vname=fieldnames(F);
   for K=1:length(vname)
       disp(sprintf('%g/%g ''Lwave'' saved',K,length(vname)))
       Wlen=getfield(F,char(vname(K)));
       Wlen(find(Wlen == -9)) = 10.0;   % flag
       Wlen(find(isnan(Wlen))) = 0.0;   % land
       ncwrite(fname,'Lwave',squeeze(Wlen).',[1 1 K])
       clear Wlen;
   end
   clear F vname K 
end

if(get_break)
   disp('   qb -> Wave break')
   F=load(break_file);
   vname=fieldnames(F);
   for K=1:length(vname)
       disp(sprintf('%g/%g ''Qb'' saved',K,length(vname)))
       Qb=getfield(F,char(vname(K)));
       Qb(find(Qb == -9)) = 10.0;   % flag
       Qb(find(isnan(Qb))) = 0.0;   % land
       ncwrite(fname,'Wave_break',squeeze(Qb).',[1 1 K])
       clear Qb;
   end
   clear F vname K 
end

if(get_xp)
   disp('   xp')
   F=load(xp_file');
   vname=fieldnames(F);
   disp(sprintf('%g/%g ''Xp'' saved',length(vname)))
   Xp=getfield(F,char(vname(1)));
   ncwrite(fname,'xp',squeeze(Xp).')
   clear Xp
   clear F vname
end

if(get_yp)
   disp('   yp')
   F=load(yp_file');
   vname=fieldnames(F);
   disp(sprintf('%g/%g ''Yp'' saved',length(vname)))
   Yp=getfield(F,char(vname(1)));
   ncwrite(fname,'yp',squeeze(Yp).')
   clear Yp
   clear F vname
end

jdshift=datenum(1968,5,23,0,0,0);
times=datenum(yr,mo,dy,hr,mm,ss);
% Modified Julian Days = (Julian Day - 2440000.) = days since 1968-05-23 00:00 UTC
mjd_days=times-jdshift;
ncwrite(fname,'wave_time',mjd_days);
   

