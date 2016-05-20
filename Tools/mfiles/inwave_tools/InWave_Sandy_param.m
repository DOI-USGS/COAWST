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

make_InWave_grd=0;   % this needs to be created, or loaded if it exists.
make_InWave_ini=1;
make_InWave_bnd=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                GRID AND BATHYMETRY DEFINITION                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%filepath='Projects\Sandy_InWave\';
filepath='C:\work\models\COAWST_tests\inwave_test2\Projects\Sandy_InWave\';
if (make_InWave_grd)

else
% grd_file=strcat(filepath,'FI_breach1.nc');  % name of the grid file
  grd_file=strcat(filepath,'FI_breach3.nc');  % name of the grid file
  ang=ncread(grd_file,'angle')*180/pi;
  depth=ncread(grd_file,'h');
%?depth=ones(size(depth)).*mean(depth(:));
  [LP,MP]=size(depth);
  TA=14.6198;                % representative absolute wave period (sec)
  theta=270;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 INITIAL CONDITION DEFINITION                    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (make_InWave_ini || make_InWave_bnd )  

%  Nbins= 21;                 % number of directional bins considered in the simulation, Nbins
  Bindirs_c = [90:10:270];   % center angles of the directional bins, Nbins
  Nbins=length(Bindirs_c);
  Bindirs = [90-5:10:270+5]; %   Nbins+1
  pd=ones(size(Bindirs)).*5./(Nbins); % Nbins+1
  
end  

if (make_InWave_ini)  
    
  ini_file=strcat(filepath,'InWave_Sandy_ini_grd3.nc');  % name of the initial file

  Ac=ones(LP,MP,Nbins).*0;
  Cx=ones(LP-1,MP,Nbins).*0;
  Cy=ones(LP,MP-1,Nbins).*0;
  Ct=ones(LP,MP,Nbins+1).*0;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 BOUNDARY CONDITION DEFINITION                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (make_InWave_bnd)

  bnd_file=strcat(filepath,'InWave_Sandy_bnd_grd3.nc');  % name of the boundary file

  % Duration of the simulation and time increment for the boundaries

  % call to read swan and create the eta time series
  CREATE_ENVELOPE_SWAN
%
%  load('realizations.mat');
  dt= time(2)-time(1);                % time increment for the boundaries (seconds)
  drtn= time(end)-time(1);         % total duration of the simulation (seconds)
%  time=[0:dt:drtn];


  % Specify by 1 the open boundary: N E S W
  obc=[0 0 1 0];

 % Specify number of directions at the boundaries (we have to specify at
 % least 2)
  %Nbins_bnd= 19;          % number of directional bins with energy at the boundaries
 % dir_bnd= [90:10:270];  % center angle (degrees of the bin containing energy at the boudaries
%  dumd=find(dir_bnd==theta);

% if sum(ismember(dir_bnd,Bindirs_c)) ~= Nbins_bnd;
%   bin_error=1;
% else
%   bin_error=0;
% end


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
    Ac_south=zeros(LP,Nbins,length(time));
    Ac_west=zeros(MP,Nbins,length(time));
    Ac_north=zeros(LP,Nbins,length(time));
    Ac_east=zeros(MP,Nbins,length(time));

    E=1/8*1025*9.81*(2.*eta_tot).^2;
    Ac1=E./(2*pi/Tm);
    Ac1(isnan(Ac1)==1)=0.0;
% hsig = sqrt2 hrms

% now redistribute the AC1 into directional bins
dfreq=diff([0;FA1(:,1)]).';
zz=(dfreq*SSS);
figure
plot(DIR1(1,:),zz,'+')
[DIR2,I] = sort(DIR1(1,:),2);
zz2=zz(I);
hold on
plot(DIR2,zz2,'r+')
zz3=interp1(DIR2,zz2,Bindirs_c);
plot(Bindirs_c,zz3,'go')
wghts=zz3./sum(zz3);


for mm=1:Nbins
  if obc(4)==1
      for i=1:MP
          Ac_west(i,mm,:)=Ac1(:).*wghts(mm);
      end
  end
  
  if obc(1)==1
      for i=1:LP
        Ac_north(i,mm,:)=Ac1(:).*wghts(mm);
      end
  end

  if obc(3)==1
      for i=1:LP
          Ac_south(i,mm,:)=Ac1(:).*wghts(mm);
      end
  end
    
  if obc(2)==1
      for i=1:MP
        Ac_east(i,mm,:)=Ac1(:).*wghts(mm);
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



