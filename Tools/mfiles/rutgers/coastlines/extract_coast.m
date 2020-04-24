%
%  EXTRACT_COAST:  Driver to extract coastline data
%
%  This is a user modifiable script that can be used to extract coastline
%  data from the GSHHS database at the specified coordinates.
%

% svn $Id: extract_coast.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%


 job='seagrid';            % Prepare coastlines for SeaGrid
%job='ploting';            % Prepare coastlines for NCAR ploting programs

%database='full';          % Full resolution database
%database='high';          % High resolution database
 database='intermediate';  % Intermediate resolution database
%database='low';           % Low resolution database
%database='crude';         % crude resolution database

switch job
  case 'seagrid'
    Oname='uswest_coast.mat';
  case 'ploting'
    Oname='uswest_inter.cst';
end

 GSHHS_DIR='~/ocean/GSHHS/Version_1.2/';
%GSHHS_DIR='~/ocean/GSHHS/Version_1.5/';
%GSHHS_DIR='~/ocean/GSHHS/Version_2.3.6/';

switch database
  case 'full'
    Cname=strcat(GSHHS_DIR,'gshhs_f.b');
    name='gshhs_f.b';
  case 'high'
    Cname=strcat(GSHHS_DIR,'gshhs_h.b');
    name='gshhs_h.b';
  case 'intermediate'
    Cname=strcat(GSHHS_DIR,'gshhs_i.b');
    name='gshhs_i.b';
  case 'low'
    Cname=strcat(GSHHS_DIR,'gshhs_l.b');
    name='gshhs_l.b';
  case 'crude'
    Cname=strcat(GSHHS_DIR,'gshhs_c.b');
    name='gshhs_c.b';
end

Llon=-134.0;              % Left   corner longitude     % US west Coast
Rlon=-118.0;              % Right  corner longitude
Blat=33.0;                % Bottom corner latitude
Tlat=49.0;                % Top    corner latitude

spval=999.0;              % Special value

%--------------------------------------------------------------------------
%  Extract coastlines from GSHHS database.
%--------------------------------------------------------------------------

disp(['Reading GSHHS database: ',name]);
[Coast]=r_gshhs(Llon,Rlon,Blat,Tlat,Cname);

disp('Processing read coastline data');
switch job
  case 'seagrid'
    [C]=x_gshhs(Llon,Rlon,Blat,Tlat,Coast,'patch');
  case 'ploting'
    [C]=x_gshhs(Llon,Rlon,Blat,Tlat,Coast,'on');
end

%--------------------------------------------------------------------------
%  Save extrated coastlines.
%--------------------------------------------------------------------------

lon=C.lon;
lat=C.lat;

switch job
  case 'seagrid'
    save(Oname,'lon','lat');
  case 'ploting'
    x=lon;
    y=lat;
    ind=find(isnan(x));
    if (~isempty(ind))
      if (length(ind) == length(C.type))
        x(ind)=C.type;
        y(ind)=spval;

% Cliping of out-of-range values failed. Try original values.
      
      elseif (length(ind) == length(Coast.type)) 
	x(ind)=Coast.type;
        y(ind)=spval;
      else      
	x(ind)=1;
        y(ind)=spval;
      end
    end
    fid=fopen(Oname,'w');
    if (fid ~= -1)
      for i=1:length(x)
        fprintf(fid,'%11.6f  %11.6f\n',y(i),x(i));
      end
      fclose(fid);
    end
end
