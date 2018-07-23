% THIS FILE CONTAINS THE DEFINITIONS OF THE PARAMETERS NEEDED TO CREATE THE 
% INPUT FILES FOR THE INWAVE MODEL for the InWave_shoreface test case:
%
% InWave_shoreface_grd.nc
% InWave_shoreface_ini.nc
% InWave_shoreface_bry.nc
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                DEFINE WHICH FILES TO GENERATE                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_InWave_grd=0;
make_InWave_ini=1;
make_InWave_bry=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%        GENERAL PARAMETERS: THESE NEED TO BE DEFINED ALWAYS       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lp = 220;                % number of total rho points in the xi direction
Mp = 8;                  % number of total rho points in the eta direction
TA = 5.0;                % representative absolute wave period (sec)
theta=270;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                GRID AND BATHYMETRY DEFINITION                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath='..\..\..\Projects\Lip\';

if (make_InWave_grd)

%  Need to define: 1) grid_name
%                  2,3,4,5) x, y, dx, dy
%                  6) depth
%                  7) angle
%                  8) mask
%                  9) f
%                  10) spherical

%   1) grid_name
%   I made a grid file with mfiles/mtools/create_roms_xygrid.m 
%   That grid is for roms and we can use same grid for InWave.
%   So we really dont need to do this step.
  netcdf_load([filepath,'lip_roms_grid.nc'])

% 2,3,4,5) x, y, dx, dy - Grid dimensions

% 6) depth - Bathymetry characteristics

% 7) angle - set grid angle

% 8) mask - set masking
  
% 9) f - set coriolis f

% 10) spherical - set if use spherical (F=no; T=yes)
% spherical='F';
  
% 11 - These are other fun things that might be needed for the open boundary.
  % compute wave number
  L0=(9.81*TA.^2)./(2*pi);
  L=zeros(size(L0));
  
  for nn=1:100
      L=L0.*tanh((h.*2*pi)./L);
  end
  
  k=(2*pi)./L;
  
  % compute wave celerity
  
  C=L./TA;
  
  % compute wave group celerity
  
  Cg=(C.*0.5).*(1+((k.*2).*h)./sinh((k.*2).*h));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 INITIAL CONDITION DEFINITION                    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (make_InWave_ini)

  Nbins= 1;                        % number of directional bins considered in the simulation
  Bindirs_centers = [270];         % center angles of the directional bins.
  
  ini_file=strcat(filepath,'InWave_lip_ini.nc');  % name of the initial file

  Ac=ones(Lp  ,Mp  ,Nbins).*0;
  Cx=ones(Lp-1,Mp  ,Nbins).*0;
  Cy=ones(Lp  ,Mp-1,Nbins).*0;
  Ct=ones(Lp  ,Mp  ,Nbins+1).*0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 BOUNDARY CONDITION DEFINITION                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (make_InWave_bry)

  bry_file=strcat(filepath,'InWave_lip_bry.nc');  % name of the boundary file

  % Duration of the simulation and time increment for the boundaries
  
 %dt= 1.0;            % time increment for the boundaries (seconds)
  dt= 0.2;            % time increment for the boundaries (seconds)
  drtn= 3600*6;       % total duration of the simulation (seconds)

  time=[0:dt:drtn];
  
  % Specify by 1 the open boundary: N E S W
  obc=[0 0 0 1];

  Nbins_bnd = 3;              % number of directional bins with energy at the boundaries
  dir_bnd = [269:1:271];    % center angle (degrees of the bin containing energy at the boudaries
  
  bin_error=0;
  
% SPECIFY CHARCATERISTICS OF Boundary waves

% call the jonswap m file
% For the Lip 2e case I used:
  df=0.001;
  f   = [0.001:df:1.0];
  Ohm = 2*pi./f;
  Hs  = 1.4;
  Tp  = 5;
  [ S, Amp, Phase ] = jonswap_lip( Ohm, Hs, Tp);
  figure;   plot(f,S)
%
% Create inwave envelope time series
% 
  fmin=f(1);
  fmax=f(end);
  nfreq=length(f);
  Insteps=length(time);
  CompFn=zeros(1,Insteps);
  for j=1:nfreq
    phase(j)=rand(1);
    amp(j)=sqrt(2*S(j)*df);
  end
  for i=1:Insteps
    for j=1:nfreq
      cff=(i-1)*dt;
      CompFn(i)=CompFn(i)+amp(j)*cos(2*pi*f(j)*cff+phase(j)*2*pi);
    end
  end
  figure
  hold on
  plot(time,CompFn,'r-')
  x_hil=hilbert(CompFn);
  amp_hil=sqrt(real(x_hil).^2+imag(x_hil).^2);
  plot(time,amp_hil,'g','LineWidth',2)
  %  avg for 5 secs
  if (0)
    cff=amp_hil;
    for mm=3:length(amp_hil)-3
       amp_hil(mm)=sum(cff(mm-2:mm+2))/5;
    end
    amp_hil(1:2)=amp_hil(3);
    amp_hil(end-2:end)=amp_hil(end-3);
    plot(time,amp_hil,'b','LineWidth',2)
  end
  
 % West
   Ac_west=zeros(Mp,Nbins_bnd,length(time));
   for i=1:Mp
     for j=1:Nbins_bnd
       cff=0.5*9.81*1025.0*(amp_hil.^2);
       Ac_west(i,j,:)=cff*TA/(2*pi);
     end
   end
  
 % if obc(1)==1
    Ac_north=zeros(Lp,Nbins_bnd,length(time));
 % end

 % if obc(3)==1
    Ac_south=zeros(Lp,Nbins_bnd,length(time));
%  end
    
%  if obc(2)==1
    Ac_east=zeros(Mp,Nbins_bnd,length(time));
%  end

  if obc(4)==1
    TA_west=TA;
  end
  if obc(2)==1
    TA_east=TA;
  end
  if obc(3)==1
    TA_south=TA;
  end
  if obc(1)==1
    TA_north=TA;
  end

end



