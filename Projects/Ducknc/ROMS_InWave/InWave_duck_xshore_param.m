% THIS FILE CONTAINS THE DEFINITIONS OF THE PARAMETERS NEEDED TO CREATE THE 
% INPUT FILES FOR THE INWAVE MODEL:
%
% InWave_grd.nc
% InWave_ini.nc
% InWave_bry.nc
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                DEFINE WHICH FILES TO GENERATE                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_InWave_grd=1;
make_InWave_ini=1;
make_InWave_bry=0;

%cd 'd:/Projects/Duck_swash/model/run_files/';
cd D:\models\COAWST_updates\COAWST_v3.9\Ducknc_xshore_test\Projects\Ducknc\ROMS_InWave

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                GRID AND BATHYMETRY DEFINITION                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (make_InWave_grd)
%  Need to define: 1) grid_name
%                  2,3,4,5) x, y, dx, dy
%                  6) depth
%                  7) angle
%                  8) mask
%                  9) f
%                  10) spherical

  load bath.mat

% 1) grid_name
%   used get_duck_data from
%   url='https://chldata.erdc.dren.mil/thredds/dodsC/frf/geomorphology/DEMs/surveyDEM/data/FRF_geomorphology_DEMs_surveyDEM_20231109.nc'
    grd_file='duck_xshore_grd.nc';                   % name of the grid file decrease dy to make it 3 cells wide, increas 2m distance

% 2,3,4,5) x, y, dx, dy - Grid dimensions
    xx=xFRF;
    clear x1 x2
    % make x more refined at coast
    x1=[50:2:300];
    dx1=gradient(x1');
    x2(1)=x1(end)+dx1(end)*1.1;
    dx2(1)=x2(1)-x1(end);
    for mm=2:44
      x2(mm)=x2(mm-1)+dx2(mm-1)*1.1;
      dx2(mm)=min(x2(mm)-x2(mm-1),20);
    end
    x=[x1 x2]';
    y=[0:20:80]';
    x=repmat(x,1,length(y));
    y=repmat(y,1,size(x,1))';
    dx=gradient(x'); dx=dx.';
    dy=gradient(y);

% 6) depth - Bathymetry characteristics: 2 sources
    figure
    plot(xFRF,elev(:,42),'k-')
    %
    depth=zeros(size(x));
    for yy=1:size(x,2)
      depth(:,yy)=interp1(xFRF,-elev(:,42),x(:,1));
    end
    hold on
    plot(x(:,3),-depth(:,3),'r*-')
    legend('xFRF data','New Xsection')

% 7) angle - set grid angle
    roms_angle=zeros(size(depth));

% 8) mask - set masking
    mask_rho=ones(size(depth));
    mask_rho(:,1)=0;
    mask_rho(:,2)=0;
    mask_rho(:,end-1)=0;
    mask_rho(:,end)=0;

% 9) f - set coriolis f
    f=zeros(size(depth));

% 10) spherical - set if use spherical (F=no; T=yes)
    spherical='F';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 INITIAL CONDITION DEFINITION                    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (make_InWave_ini)

  %Bindirs_centers = [60:10:120];  % center angles of the directional bins,
  Bindirs_centers = [90];  % center angles of the directional bins,
  Nbins= length(Bindirs_centers); % number of computational directional bins
                                  % size Nbins. Directions coming from.
  TA= 12.0;                       % representative absolute wave period (sec)

  ini_file='duck_xshore_InWave_ini.nc';  % name of the initial file

  [Lp,Mp]=size(depth);
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



