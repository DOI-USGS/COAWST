function [Fout]=mercator2roms(Vname,S,Finp,lon,lat,mask,depth);

%
% MERCATOR2ROMS:  Interpolates data from Mercator to ROMS grid.
%
% [Fout]=mercator2roms(Vname,S,Finp,lon,lat,mask,depth)
%
% This function interpolates requested variable form Mercator to
% ROMS grid.  All the fields are interpolated on the RHO-points
% grid including velocity fields. This is done in the event the
% velocity vectors must be rotated into a ROMS curvilinear grid.
% Such processing is done elsewhere.
%
% On Input:
%
%    Vname       ROMS variable to process (string)
%    S           ROMS grid data (structure array):
%
%                  S.lon_rho          Longitude at RHO-points
%                  S.lat_rho          Latitude  at RHO-points
%                  S.z_r              3D depths at RHO-points
%                  S.mask_rho         land/sea mask on RHO-points, if any
%
%    Finp        Mercator data (2D or 3D array)
%    lon         Mercator data longitude (2D array)
%    lat         Mercator data latitude (2D array)
%    mask        Mercator data Land/Sea mask (2D array)
%    depth       Mercator data depth (1D array, only needed for 3D fields)
%
% On Output:
%
%    Fout        Interpolated data on ROMS grid (2D or 3D array)
%

% svn $Id: mercator2roms.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Determine if processing a 2D or 3D field.

switch Vname
  case ('zeta')
    is2d=1;
  otherwise
    is2d=0;
end,

if (~is2d & nargin < 7),
  error([ 'MERCATOR2ROMS - missing Mercator depth argument: depth']);
end,

% Check input ROMS grid structure for required fields.

if (isfield(S,'lon_rho')),
  rlon=S.lon_rho;
else,
  error([ 'MERCATOR2ROMS - cannot find longitude field: lon_rho, ', ...
          'in structure array S']);
end,

if (isfield(S,'lat_rho')),
  rlat=S.lat_rho;
else,
  error([ 'MERCATOR2ROMS - cannot find latitude field: lat_rho, ', ...
          'in structure array S']);
end,

if (isfield(S,'mask_rho')),
  rmask=S.mask_rho;
else,
  rmask=ones(size(rlon));
end,

if (~is2d),
  if (isfield(S,'z_r')),
    z_r=S.z_r;
  else,
    error([ 'MERCATOR2ROMS - cannot find depth field: z_r, ', ...
            'in structure array S']);
  end,
end,

%----------------------------------------------------------------------------
%  Trim Mercator data form large input grid to make interpolations more
%  tractable.  This possible since the Mercator grid has a North-South
%  orientation.
%----------------------------------------------------------------------------

ii=min(find(lon(:,1) > (min(rlon(:))-1))): ...
   max(find(lon(:,1) < (max(rlon(:))+1)));

jj=min(find(lat(1,:) > (min(rlat(:))-1))): ...
   max(find(lat(1,:) < (max(rlat(:))+1)));

x=lon(ii,jj);
y=lat(ii,jj);

%----------------------------------------------------------------------------
%  Interpolate 2D field.
%----------------------------------------------------------------------------

if (is2d),

  M=mask(ii,jj);           % Sample Mercator Land/Sea mask
  wet=find(M > 0);         % Mercator wet points
  
  F=Finp(ii,jj);
  
  Fout=rmask.*griddata(x(wet),y(wet),F(wet),rlon,rlat,'linear');

%  Replace NaNs with neareast neighbor value.

  indnan=find(rmask == 1 & isnan(Fout) == 1);
  indwet=find(rmask == 1 & isnan(Fout) == 0);

  icnan=length(indnan);
%  disp(['  ', sprintf('%4s',Vname), ', icNaN = ', num2str(icnan)]);

  for m=1:length(indnan),
    d=sqrt((rlon(indnan(m))-rlon(indwet)).*(rlon(indnan(m))-rlon(indwet)) + ...
           (rlat(indnan(m))-rlat(indwet)).*(rlat(indnan(m))-rlat(indwet)));
    near=find(d == min(d));
    if ~isempty(near),
      Fout(indnan(m))=Fout(indwet(near));
    end,
  end,

  indwet=find(rmask == 1);
  Fmin=min(min(Fout(indwet)));
  Fmax=max(max(Fout(indwet)));

end,

%----------------------------------------------------------------------------
%  Interpolate 3D field.
%----------------------------------------------------------------------------

if (~is2d),

  icnan=0;

%  Determine how many level of Mercator data to process. Sometimes the
%  bottom level(s) have zero values.

  Nlev=size(Finp,3);
  
  for k=size(Finp,3):-1:20,
    if (min(min(Finp(:,:,k))) == 0 & max(max(Finp(:,:,k))) == 0),
      Nlev=k-1;
    end,
  end,
  
%  Horizontal interpolation into Mercator Z-levels.

  for k=1:Nlev,

    M=mask(ii,jj,k);         % Sample Mercator Land/Sea mask
    wet=find(M > 0);         % Mercator wet points

    F=squeeze(Finp(ii,jj,k));
    Fwork=rmask.*griddata(x(wet),y(wet),F(wet),rlon,rlat,'linear');
  
%  Replace NaNs with neareast neighbor value.

    indnan=find(rmask == 1 & isnan(Fwork) == 1);
    indwet=find(rmask == 1 & isnan(Fwork) == 0);

    icnan=icnan+length(indnan);
%   disp(['  ', sprintf('%4s',Vname), ', Level = ', sprintf('%2i',k), ...
%         ' icNaN = ', num2str(icnan)]);
    
    for m=1:length(indnan),
      d=sqrt((rlon(indnan(m))-rlon(indwet)).*(rlon(indnan(m))-rlon(indwet)) + ...
             (rlat(indnan(m))-rlat(indwet)).*(rlat(indnan(m))-rlat(indwet)));
      near=find(d == min(d));
      if ~isempty(near),
        Fwork(indnan(m))=Fwork(indwet(near));
      end,
    end,

    Flev(:,:,k)=Fwork;
  
  end,

%  Vertical interpolation into ROMS grid. Notice that Mercator depths
%  are positive.

  [Im,Jm,N]=size(z_r);

  Fout=zeros(size(z_r));
  depth=-abs(depth(1:Nlev));
  depth(1)=0.0;                % reset surface level to zero (0.49 orinally)
  
  for j=1:Jm,
    for i=1:Im,
      if (rmask(i,j) == 1),
	Fout(i,j,:)=interp1(depth,squeeze(Flev(i,j,:)),squeeze(z_r(i,j,:)));
      end,
    end,
  end,

  indwet=find((repmat(rmask,[1,1,N])) == 1);
  Fmin=min(min(min(Fout(indwet))));
  Fmax=max(max(max(Fout(indwet))));
  
end,

disp([ 'Interpolated  ', sprintf('%4s',Vname), ...
       '  number of NaNs replaced with nearest value = ',num2str(icnan)]);
disp([ '                    Min=', sprintf('%12.5e',Fmin), ...
       '  Max=', sprintf('%12.5e',Fmax)]);

return
