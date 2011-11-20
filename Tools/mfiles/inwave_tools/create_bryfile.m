function create_bryfile(bryname,grdname,time,zeta,ubar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function create_bryfile(bryname,grdname,title,obc...
%                          theta_s,theta_b,hc,N,...
%                          time,cycle,clobber);
%
%   This function create the header of a Netcdf climatology 
%   file.
%
%   Input:
%
%   bryname      Netcdf climatology file name (character string).
%   grdname      Netcdf grid file name (character string).
%   obc          open boundaries flag (1=open , [S E N W]).
%   theta_s      S-coordinate surface control parameter.(Real)
%   theta_b      S-coordinate bottom control parameter.(Real)
%   hc           Width (m) of surface or bottom boundary layer 
%                where higher vertical resolution is required 
%                during stretching.(Real)
%   N            Number of vertical levels.(Integer)
%   time         time.(vector)
%   cycle        Length (days) for cycling the climatology.(Real)
%   clobber      Switch to allow or not writing over an existing
%                file.(character string)
% 
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2001-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp([' Creating the file : ',bryname])
disp(' ')
%
%  Read the grid file and check the topography
%
nc = netcdf(grdname, 'nowrite');
h=nc{'h'}(:);
maskr=nc{'mask_rho'}(:);
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
status=close(nc);
L=Lp-1;
M=Mp-1;

%
%  Create the boundary file
%
type = 'BOUNDARY file' ; 
history = 'ROMS' ;
nc = netcdf(bryname,'clobber');
result = redef(nc);
%
%  Create dimensions
%
nc('xi_u') = L;
nc('xi_rho') = Lp;
nc('eta_v') = M;
nc('eta_rho') = Mp;
nc('ocean_time') = length(time);
nc('bry_time') = length(time);
nc('zeta_time') = length(time);
nc('v2d_time') = length(time);
nc('one') = 1;

%
nc{'bry_time'} = ncdouble('bry_time') ;
nc{'bry_time'}.long_name = ncchar('time for boundary');
nc{'bry_time'}.long_name = 'time for boundary';
nc{'bry_time'}.units = ncchar('seconds');
nc{'bry_time'}.units = 'seconds';
nc{'bry_time'}.cycle_length = 1000000;%

nc{'ocean_time'} = ncdouble('ocean_time') ;
nc{'ocean_time'}.long_name = ncchar('time for boundary');
nc{'ocean_time'}.long_name = 'time for boundary';
nc{'ocean_time'}.units = ncchar('seconds');
nc{'ocean_time'}.units = 'seconds';
nc{'ocean_time'}.cycle_length = 1000000;%

nc{'zeta_time'} = ncdouble('zeta_time') ;
nc{'zeta_time'}.long_name = ncchar('time for boundary');
nc{'zeta_time'}.long_name = 'time for boundary';
nc{'zeta_time'}.units = ncchar('seconds');
nc{'zeta_time'}.units = 'seconds';
nc{'zeta_time'}.cycle_length = 1000000;%

nc{'v2d_time'} = ncdouble('v2d_time') ;
nc{'v2d_time'}.long_name = ncchar('time for boundary');
nc{'v2d_time'}.long_name = 'time for boundary';
nc{'v2d_time'}.units = ncchar('seconds');
nc{'v2d_time'}.units = 'seconds';
nc{'v2d_time'}.cycle_length = 1000000;%
%
nc{'zeta_east'} = ncdouble('zeta_time','eta_rho') ;
nc{'zeta_east'}.long_name = ncchar('eastern boundary sea surface height');
nc{'zeta_east'}.long_name = 'eastern boundary sea surface height';
nc{'zeta_east'}.units = ncchar('meter');
nc{'zeta_east'}.units = 'meter';

nc{'ubar_east'} = ncdouble('v2d_time','eta_rho') ;
nc{'ubar_east'}.long_name = ncchar('eastern boundary barotropic current');
nc{'ubar_east'}.long_name = 'eastern boundary barotropic current';
nc{'ubar_east'}.units = ncchar('meter');
nc{'ubar_east'}.units = 'meter';

nc{'vbar_east'} = ncdouble('v2d_time','eta_v') ;
nc{'vbar_east'}.long_name = ncchar('eastern boundary barotropic current');
nc{'vbar_east'}.long_name = 'eastern boundary barotropic current';
nc{'vbar_east'}.units = ncchar('meter');
nc{'vbar_east'}.units = 'meter';


%
% Create global attributes
%
nc.date = ncchar(date);
nc.date = date;
nc.clim_file = ncchar(bryname);
nc.clim_file = bryname;
nc.grd_file = ncchar(grdname);
nc.grd_file = grdname;
nc.type = ncchar(type);
nc.type = type;
nc.history = ncchar(history);
nc.history = history;
%
% Leave define mode
%
result = endef(nc);
 
nc{'bry_time'}(:) =  time(:); 
nc{'ocean_time'}(:) =  time(:);
nc{'zeta_time'}(:) =  time(:); 
nc{'v2d_time'}(:) =  time(:); 

for ii=1:Mp
    zeta1(:,ii)=zeta(1,:)';
end

for ii=1:Mp
    ubar1(:,ii)=ubar(1,:)';
end

nc{'zeta_east'}(:) =  zeta1(:,:)-mean(zeta); 
nc{'ubar_east'}(:) =  ubar1(:,:)-mean(ubar); 
nc{'vbar_east'}(:) =  0; 


close(nc)
return


