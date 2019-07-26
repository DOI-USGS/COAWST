% THIS FILE CONTAINS THE DEFINITIONS OF THE PARAMETERS NEEDED TO CREATE THE 
% INPUT FILES FOR THE INWAVE MODEL:
%
% InWave_grd.nc
% InWave_ini.nc
% InWave_bnd.nc
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                DEFINE WHICH FILES TO GENERATE                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_InWave_grd=1;
make_InWave_ini=1;
make_InWave_bry=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                GRID AND BATHYMETRY DEFINITION                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath='Projects\Inlet_test\InWave\';

if (make_InWave_grd)
%  Need to define: 1) grid_name
%                  2,3,4,5) x, y, dx, dy
%                  6) depth
%                  7) angle
%                  8) mask
%                  9) f
%                  10) spherical

% 1) grid_name
    grd_file=strcat(filepath,'InWave_inlet_test_grd.nc');  % name of the grid file

% 2,3,4,5) x, y, dx, dy - Grid dimensions
    ncellsx=77;  dx=20;
    ncellsy=72;  dy=20;
    xx=[-10:dx:dx*(ncellsx-1)-10];
    yy=[-10:dy:dy*(ncellsy-1)-10];
  % 
    x=repmat(xx',1,length(yy));
    y=repmat(yy,length(xx),1);
    [Lp,Mp]=size(x);

% 6) depth - Bathymetry characteristics
    depth=zeros(size(x))+4;
    zz=[4+0.016*dy.*([37:72]-36)];
    depth(:,37:72)=repmat(zz,length(xx),1);

% 7) angle - set grid angle
     roms_angle=zeros(size(depth));

% 8) mask - set masking
    mask_rho=ones(size(depth));
    mask_rho(:,1)=0;
    mask_rho(1,1:36)=0;
    mask_rho(end,1:36)=0;
    mask_rho(1:33,36)=0;
    mask_rho(45:end,36)=0;

% 9) f - set coriolis f
    f=zeros(size(depth));

% 10) spherical - set if use spherical (F=no; T=yes)
    spherical='F';

else
  grd_file=strcat(filepath,'inlet_test_grid.nc');  % name of the grid file
  ang=ncread(grd_file,'angle')*180/pi;
  depth=ncread(grd_file,'h');
  [Lp,Mp]=size(depth);
  TA= 8.3;                % representative absolute wave period (sec)
  theta=270;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 INITIAL CONDITION DEFINITION                    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (make_InWave_ini)

% Nbins= 36;                     % number of computational directional bins
% Bindirs_centers = [5:10:355];  % center angles of the directional bins, 
%                                  size Nbins. Directions coming from.

  Nbins= 19;                     % number of computational directional bins
  Bindirs_centers = [-90:10:90]; % center angles of the directional bins,
%                                  size Nbins. Directions coming from.
  TA= 8.3;                       % representative absolute wave period (sec)

  ini_file=strcat(filepath,'InWave_inlet_test_ini.nc');  % name of the initial file

  Ac=ones(Lp  ,Mp  ,Nbins).*0;
  Cx=ones(Lp-1,Mp  ,Nbins).*0;
  Cy=ones(Lp  ,Mp-1,Nbins).*0;
  Ct=ones(Lp  ,Mp  ,Nbins+1).*0;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 BOUNDARY CONDITION DEFINITION                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (make_InWave_bry)
  %%% This Inlet_test case does not use a bry file, instead it reads a 2dspec file. %%%

  bry_file=strcat(filepath,'InWave_inlet_test_bry.nc');  % name of the boundary file

  % Specify by 1 the open boundary: N E S W
  obc=[1 1 0 1];

 % Specify number of directions at the boundaries (we have to specify at
 % least 2)
   Nbins_bnd= 3;          % number of directional bins with energy at the boundaries
   dir_bnd= [355 360 365];  % center angle (degrees of the bin containing energy at the boudaries

   if sum(ismember(dir_bnd,Bindirs_c)) ~= Nbins_bnd; 
     bin_error=1;
   else
     bin_error=0;
   end      


% compute wave number
  if (0)
    L0=(9.81*TA.^2)./(2*pi);
    L=zeros(size(L0));
  
    for nn=1:100
        L=L0.*tanh((depth.*2*pi)./L);
    end
  
    k=(2*pi)./L;
  
    % compute wave celerity
  
    C=L./TA;
  
    % compute wave group celerity
  
    Cg=(C.*0.5).*(1+((k.*2).*depth)./sinh((k.*2).*depth));

      % SPECIFY WAVE GROUPS
    end


    Tm=TA;

    close all
    Ac_south=zeros(LP,Nbins_bnd,length(time));
    Ac_west=zeros(MP,Nbins_bnd,length(time));
    Ac_north=zeros(LP,Nbins_bnd,length(time));
    Ac_east=zeros(MP,Nbins_bnd,length(time));

    E=1/8*1025*9.81*(2.*eta_tot).^2;
    Ac1=E./(2*pi/Tm);
    Ac1(isnan(Ac1)==1)=0.0;

for mm=1:Nbins_bnd
  if obc(4)==1
      for i=1:MP
          Ac_west(i,mm,:)=Ac1(:);
      end
  end
  
  if obc(1)==1
      for i=1:LP
        Ac_north(i,mm,:)=Ac1(:);
      end
  end

  if obc(3)==1
      for i=1:LP
          Ac_south(i,mm,:)=Ac1(:);
      end
  end
    
  if obc(2)==1
      for i=1:MP
        Ac_east(i,mm,:)=Ac1(:);
      end
  end
end


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



