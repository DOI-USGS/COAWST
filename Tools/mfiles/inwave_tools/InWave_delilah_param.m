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

filepath='E:\data\models\InWave\readswan\Projects\Delilah\coawst\';

if (make_InWave_grd)
%  Need to define: 1) grid_name
%                  2,3,4,5) x, y, dx, dy
%                  6) depth
%                  7) angle
%                  8) mask
%                  9) f
%                  10) spherical

% 1) grid_name
%   grd2 used 12Oct bathy, and mirrored north and south
%   grd_file=strcat(filepath,'InWave_delilah_grd2.nc');  % name of the grid file
%   grd3 used 07Oct bathy, and mirrored north and south
%   grd_file=strcat(filepath,'InWave_delilah_grd3.nc');  % name of the grid file
%   grd_file=strcat(filepath,'InWave_delilah_grd3_per.nc');  % name of the grid file
    grd_file=strcat(filepath,'InWave_delilah_grd4_per.nc');  % name of the grid file

% 2,3,4,5) x, y, dx, dy - Grid dimensions
    x=[50:2:950];
    y=[600:5:1500];
    x=repmat(x,length(y),1)';
    y=repmat(y,size(x,1),1);
    dx=gradient(x);
    dy=gradient(y')';

% 6) depth - Bathymetry characteristics: 2 sources
    % this is the bathy I interpolated from Fig 2 in 1990_DELILAH_Nearshore_Experiment_Summary_Report.pdf
    z1=load([filepath,'jcw_8m_bathy_400.txt']);
    %
    % this is high res data just for the Delilah array area.
    % the X and Y positions stay the same, so i just list them here.
    clear X Y Z
    X=[50.0 58.8 67.5  76.3 85.0 93.8 102.5 111.3 120.0 128.8 137.5 146.3 155.0 163.8 ...
       172.5 181.3 190.0 198.8 207.5 216.3 225.0 233.8 242.5 251.3 260.0 268.8 277.5 ...
       286.3 295.0 303.8 312.5 321.3 330.0 338.8 347.5 356.3 365.0 373.8 382.5 391.3 400.0];
    Y=[725 748 772 796 820 844 868 892 916 940 964 988 1012 1036 1060 1084 1108 1132 ...
       1156 1180 1204 1228 1252 1276 1300];
    %now get the Z data from a specific file
    fid=fopen([filepath,'FRF_11Oct1990.grid']);   %grd4
    %fid=fopen([filepath,'FRF_07Oct1990.grid']);    %grd3
    line=fgetl(fid);
    for mm=1:100000
      line=fgetl(fid);
      if (line==-1); break; end
      junk=str2num(line);
      Z(mm,:)=-junk(2:end);
    end
    fclose(fid);
    Z(1,end)=Z(1,end-1);  %this point looked bad
    Z(:,1)=Z(:,2);  %repeat the bottom row, it looked strange
    % here we repeat the first and last rows for X from 0 to 400.
%    Y=[600 Y 1500];
    Y=[600 Y 1400 1500];
    X=repmat(X,length(Y),1)';
    Y=repmat(Y,size(X,1),1);
%   Z=[Z(:,1) Z Z(:,end)];
    Z=[Z(:,1) Z Z(:,1) Z(:,1)];
    %combine both data sets
    XX=[X(:);z1(:,1)];
    YY=[Y(:);z1(:,2)];
    ZZ=[Z(:);z1(:,3)];
    F = TriScatteredInterp(XX,YY,ZZ);
    z=F(x,y);
    z(isnan(z))=-8;
    depth=z;
    for ii=2:size(z,1)-1
      for jj=2:size(z,2)-1
        depth(ii,jj)=0.2*(z(ii-1,jj)+z(ii,jj)+z(ii+1,jj)+z(ii,jj-1)+z(ii,jj+1));
      end
    end
    figure
    pcolorjw(x,y,depth)
    colorbar
    colormap('jet')

% 7) angle - set grid angle
    roms_angle=zeros(size(z))+18.2*pi/180.;

% 8) mask - set masking
    mask_rho=ones(size(z));

% 9) f - set coriolis f
    f=zeros(size(z));

% 10) spherical - set if use spherical (F=no; T=yes)
    spherical='F';


%%%%%%%%  here we make a taller grid to be used for swan
  if (0) 
    grd_file=strcat(filepath,'InWave_delilah_grd4_tall.nc');  % name of the grid file

% 2,3,4,5) x, y, dx, dy - Grid dimensions
    xs=[repmat(x(:,1),1,24), x];
    ys=[repmat([0:25:575],451,1), y];
    zs=[repmat(depth(:,1),1,24)*0, depth];
    for mm=1:24
     zs(:,mm)=zs(:,25);
    end
    dxs=gradient(x);
    dys=gradient(y')';

% 7) angle - set grid angle
    roms_angles=zeros(size(zs))+18.2*pi/180.;

% 8) mask - set masking
    mask_rhos=ones(size(zs));

% 9) f - set coriolis f
    fs=zeros(size(zs));

% 10) spherical - set if use spherical (F=no; T=yes)
    spherical='F';
    create_InWave_grd(xs,ys,dxs,dys,zs,roms_angles,mask_rhos,fs,spherical,grd_file)
    roms2swan(grd_file);
    !move swan_coord.grd delilah_swan_grd4_tall.grd
    !move swan_bathy.bot delilah_swan_grd4_tall.bot
  end

else
  grd_file=strcat(filepath,'delilah_grid.nc');  % name of the grid file
  ang=ncread(grd_file,'angle')*180/pi;
  depth=ncread(grd_file,'h');
  [Lp,Mp]=size(depth);
  TA= 8.3;                % representative absolute wave period (sec)
  theta=90;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 INITIAL CONDITION DEFINITION                    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (make_InWave_ini)

  Bindirs_centers = [10:5:170];  % center angles of the directional bins,
  Nbins= length(Bindirs_centers);% number of computational directional bins
%                                  size Nbins. Directions coming from.
  TA= 8.3;                       % representative absolute wave period (sec)

  ini_file=strcat(filepath,'InWave_delilah_ini.nc');  % name of the initial file

  [Lp,Mp]=size(x);
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



