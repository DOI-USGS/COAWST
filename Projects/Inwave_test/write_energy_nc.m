% Program create_Energy(ncname,grdname,obc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program create_Energy(ncname,grdname,obc);
%
%   This function creates the Energy boudnary condition for the InWave
%   model
%
%   Input:
%
%   ncname      Netcdf energy file name (character string).
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

cd 'F:\MAITANE\PROJECTS\INWAVE\matlab'

ncname='energy.nc';


% S E N W
obc=[0 0 0 1];

disp(' ')
disp([' Creating the file : ',ncname])
disp(' ')

Lp=20;   % xi_rho
Mp=20;   % eta_rho
%Nangle=101;  % Number of directional bins
Nangle=11;  % Number of directional bins

L=Lp-1;
M=Mp-1;


%
%  Create the boundary file
%

type = 'BOUNDARY file for the InWave model' ; 
nc = netcdf(ncname,'clobber');
result = redef(nc);

%
%  Read first energy file
%

file=strcat('F:\MAITANE\PROJECTS\INWAVE\InWave_Boundary\dir1.dat');
energy=load(file);
time=[0:1:length(energy)-1]';

file=strcat('F:\MAITANE\PROJECTS\INWAVE\InWave_Boundary\Directions.dat');
dir1=load(file);
dir1=dir1-90;
dum=find(dir1<=0);
dir1(dum)=360+dir1(dum);
dir=dir1(64:74)

    
%
%  Create dimensions
%

nc('xi_u') = L;
nc('xi_rho') = Lp;
nc('eta_v') = M;
nc('eta_rho') = Mp;
nc('ocean_time') = length(time);
nc('energy_angle') = Nangle;

%
%  Create variables and attributes
%

nc{'ocean_time'} = ncdouble('ocean_time') ;
nc{'ocean_time'}.long_name = ncchar('time for energy envelope');
nc{'ocean_time'}.long_name = 'time for energy envelope';
nc{'ocean_time'}.units = ncchar('seconds');
nc{'ocean_time'}.units = 'seconds';

nc{'energy_angle'} = ncdouble('energy_angle') ;
nc{'energy_angle'}.long_name = ncchar('direction respect to the north of the bin');
nc{'energy_angle'}.long_name = 'direction respect to the north of the bin';
nc{'energy_angle'}.units = ncchar('degrees');
nc{'energy_angle'}.units = 'degrees';

%

if obc(1)==1

%   Southern boundary

  nc{'energy_south'} = ncdouble('energy_angle','ocean_time','xi_rho') ;
  nc{'energy_south'}.long_name = ncchar('southern boundary wave energy envelope');
  nc{'energy_south'}.long_name = 'southern boundary wave energy envelope';
  nc{'energy_south'}.units = ncchar('Joules');
  nc{'energy_south'}.units = 'Joules';
  
end

if obc(2)==1

%   Eastern boundary

  nc{'energy_east'} = ncdouble('energy_angle','ocean_time','eta_rho') ;
  nc{'energy_east'}.long_name = ncchar('eastern boundary wave energy envelope');
  nc{'energy_east'}.long_name = 'eastern boundary wave energy envelope';
  nc{'energy_east'}.units = ncchar('Joules');
  nc{'energy_east'}.units = 'Joules';

end

if obc(3)==1
%
%   Northern boundary
%
  nc{'energy_north'} = ncdouble('energy_angle','ocean_time','xi_rho') ;
  nc{'energy_north'}.long_name = ncchar('northern boundary wave energy envelope');
  nc{'energy_north'}.long_name = 'northern boundary wave energy envelope';
  nc{'energy_north'}.units = ncchar('Joules');
  nc{'energy_north'}.units = 'Joules';
%
end
%
if obc(4)==1
%
%   Western boundary
%
  nc{'energy_west'} = ncdouble('energy_angle','ocean_time','eta_rho') ;
  nc{'energy_west'}.long_name = ncchar('western boundary wave energy envelope');
  nc{'energy_west'}.long_name = 'western boundary wave energy envelope';
  nc{'energy_west'}.units = ncchar('Joules');
  nc{'energy_west'}.units = 'Joules';

end
%
%
% Create global attributes
%

nc.date = ncchar(date);
nc.date = date;
nc.clim_file = ncchar(ncname);
nc.clim_file = ncname;
nc.type = ncchar(type);
nc.type = type;

%
% Leave define mode
%

result = endef(nc);

%
% Write variables
%

nc{'ocean_time'}(:) =  time; 
nc{'energy_angle'}(:) =  dir; 
    
%for j=1:Nangle
    
for j=64:74
    
    j1=j-63;
    
    file=strcat('F:\MAITANE\PROJECTS\INWAVE\InWave_Boundary\dir',num2str(j),'.dat');
    energy=load(file);

    if obc(1)==1
      energy_s=zeros(1,size(time),Lp);
      for i=1:Lp
          energy_s(1,:,i)=energy(:,1);
      end
      nc{'energy_south'}(j1,i,:) =  energy_s(1,:,:); 
    end 
    
    if obc(2)==1
      energy_e=zeros(1,size(time),Mp);
      for i=1:Mp
          energy_e(1,:,i)=energy(:,1);
      end
      nc{'energy_east'}(j1,:,:) =  energy_e(1,:,:); 
    end 
    
    if obc(3)==1
      energy_n=zeros(1,size(time),Lp);
      for i=1:Lp
          energy_n(1,:,i)=energy(:,1);
      end
      nc{'energy_north'}(j1,i,:) =  energy_n(1,:,:); 
    end 
    
    if obc(4)==1
      energy_w=zeros(1,size(time),Mp);
      for i=1:Mp
          energy_w(1,:,i)=energy(:,1);
      end
      nc{'energy_west'}(j1,:,:) =  energy_w(1,:,:); 
    end
    
end

close(nc)

return


