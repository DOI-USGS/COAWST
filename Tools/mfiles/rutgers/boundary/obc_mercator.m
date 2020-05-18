function [Fout]=obc_mercator(Vname,S,Finp,lon,lat,mask,depth);

%
% OBC_MERCATOR:  Interpolates OBC data from Mercator to ROMS grid.
%
% [Fout]=obc_mercator(Vname,S,Finp,lon,lat,mask,depth)
%
% This function interpolates requested variable open boundary
% conditions form Mercator to ROMS grid boundary edges. All
% fields are interpolated on the RHO-points grid including
% velocity fields. This is done in the event the velocity
% vectors must be rotated into a ROMS curvilinear grid. Such
% processing is done elsewhere.
%
% On Input:
%
%    Vname       ROMS variable to process (string)
%    S           ROMS grid data (structure array):
%
%                  S.boundary(1)      Process western  boundary
%                  S.boundary(2)      Process eastern  boundary
%                  S.boundary(3)      Process southern boundary
%                  S.boundary(4)      Process northern boundary
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
%    Fout        Interpolated data on ROMS grid (structure array)
%
%                  Fout.west          Western  boundary conditions
%                  Fout.east          Eastern  boundary conditions
%                  Fout.south         Southern boundary conditions
%                  Fout.north         Northern boundary conditions
%

% svn $Id: obc_mercator.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Set cheap nearest neighbor search.
  
CHEAP=0;

% Determine if processing a 2D or 3D field.

switch Vname
  case ('zeta')
    is2d=1;
  otherwise
    is2d=0;
end,

if (~is2d & nargin < 7),
  error([ 'OBC_MERCATOR - missing Mercator depth argument: depth']);
end,

% Check input ROMS grid structure for required fields.

if (isfield(S,'lon_rho')),
  rlon=S.lon_rho;
else,
  error([ 'OBC_MERCATOR - cannot find longitude field: lon_rho, ', ...
          'in structure array S']);
end,

if (isfield(S,'lat_rho')),
  rlat=S.lat_rho;
else,
  error([ 'OBC_MERCATOR - cannot find latitude field: lat_rho, ', ...
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
    error([ 'OBC_MERCATOR - cannot find depth field: z_r, ', ...
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
  
  Fwrk=rmask.*griddata(x(wet),y(wet),F(wet),rlon,rlat,'linear');

%  Replace NaNs with neareast neighbor value.

  icnan=0;
  Bmin(1:4)=NaN;
  Bmax(1:4)=NaN;
  
  for ib=1:4,
  
    if (S.boundary(ib)),

      switch ib
        case {1}                            % western edge
          Fobc=Fwrk(1,:);
          blon=rlon(1,:);
          blat=rlat(1,:);
          bmask=rmask(1,:);
        case {2}                            % eastern edge
          Fobc=Fwrk(end,:);
          blon=rlon(end,:);
          blat=rlat(end,:);
          bmask=rmask(end,:);
        case {3}                            % southern edge
          Fobc=Fwrk(:,1);
          blon=rlon(:,1);
          blat=rlat(:,1);
          bmask=rmask(:,1);
        case {4}                            % northern edge
          Fobc=Fwrk(:,end);
          blon=rlon(:,end);
          blat=rlat(:,end);
          bmask=rmask(:,end);
      end,

      if (CHEAP),
        indnan=find(bmask > 0 & isnan(Fobc) == 1);
        indwet=find(bmask > 0 & isnan(Fobc) == 0);
        icnan=icnan+length(indnan);
        for m=1:length(indnan),
          d=sqrt((blon(indnan(m))-blon(indwet)).*(blon(indnan(m))-blon(indwet)) + ...
                 (blat(indnan(m))-blat(indwet)).*(blat(indnan(m))-blat(indwet)));
          near=find(d == min(d));
          if ~isempty(near),
            Fobc(indnan(m))=Fobc(indwet(near));
          end,
        end,
      else,
        indnan=find(bmask > 0 & isnan(Fobc) == 1);
        indwet=find(rmask > 0 & isnan(Fwrk) == 0);
        icnan=icnan+length(indnan);
        for m=1:length(indnan),
          d=sqrt((blon(indnan(m))-rlon(indwet)).*(blon(indnan(m))-rlon(indwet)) + ...
                 (blat(indnan(m))-rlat(indwet)).*(blat(indnan(m))-rlat(indwet)));
          near=find(d == min(d));
          if ~isempty(near),
            Fobc(indnan(m))=Fwrk(indwet(near));
          end,
        end,
      end,
      
      switch ib
        case {1}
          Fout.west=Fobc;
        case {2}
          Fout.east=Fobc;
        case {3}
          Fout.south=Fobc;
        case {4}
          Fout.north=Fobc;
      end,

      indwet=find(bmask > 0);
      Bmin(ib)=min(Fobc(indwet));
      Bmax(ib)=max(Fobc(indwet));
      
    end,
  end,
end,

%----------------------------------------------------------------------------
%  Interpolate 3D field.
%----------------------------------------------------------------------------

if (~is2d),

  icnan=0;
  Bmin(1:4)=NaN;
  Bmax(1:4)=NaN;

%  Determine how many level of Mercator data to process. Sometimes the
%  bottom level(s) have zero values.

  Nlev=size(Finp,3);
  
  for k=size(Finp,3):-1:2,
    if (min(min(Finp(:,:,k))) == 0 & max(max(Finp(:,:,k))) == 0),
      Nlev=k-1;
    end,
  end,
  
%  Horizontal interpolation into Mercator Z-levels.

  for k=1:Nlev,

    M=mask(ii,jj,k);         % Sample Mercator Land/Sea mask
    wet=find(M > 0);         % Mercator wet points

    F=squeeze(Finp(ii,jj,k));
    Fwrk=rmask.*griddata(x(wet),y(wet),F(wet),rlon,rlat,'linear');
  
%  Replace NaNs with neareast neighbor value.

    for ib=1:4,

      if (S.boundary(ib)),

        switch ib
          case {1}                                   % western edge
            if (Vname == 'u' | Vname == 'v'),
              Fobc=Fwrk(1:2,:);
              blon=rlon(1:2,:);
              blat=rlat(1:2,:);
              bmask=rmask(1:2,:);
            else,
              Fobc=Fwrk(1,:);
              blon=rlon(1,:);
              blat=rlat(1,:);
              bmask=rmask(1,:);
            end,
          case {2}                                   % eastern edge
            if (Vname == 'u' | Vname == 'v'),
              Fobc=Fwrk(end-1:end,:);
              blon=rlon(end-1:end,:);
              blat=rlat(end-1:end,:);
              bmask=rmask(end-1:end,:);
            else,
              Fobc=Fwrk(end,:);
              blon=rlon(end,:);
              blat=rlat(end,:);
              bmask=rmask(end,:);
            end,
          case {3}                                   % southern edge
            if (Vname == 'u' | Vname == 'v'),
              Fobc=Fwrk(:,1:2);
              blon=rlon(:,1:2);
              blat=rlat(:,1:2);
              bmask=rmask(:,1:2);
            else,
              Fobc=Fwrk(:,1);
              blon=rlon(:,1);
              blat=rlat(:,1);
              bmask=rmask(:,1);
            end,
          case {4}                                   % northern edge
            if (Vname == 'u' | Vname == 'v'),
              Fobc=Fwrk(:,end-1:end);
              blon=rlon(:,end-1:end);
              blat=rlat(:,end-1:end);
              bmask=rmask(:,end-1:end);
            else,
              Fobc=Fwrk(:,end);
              blon=rlon(:,end);
              blat=rlat(:,end);
              bmask=rmask(:,end);
            end,
        end,

        if (CHEAP),
          indnan=find(bmask > 0 & isnan(Fobc) == 1);
          indwet=find(bmask > 0 & isnan(Fobc) == 0);
          icnan=icnan+length(indnan);
          for m=1:length(indnan),
            d=sqrt((blon(indnan(m))-blon(indwet)).*(blon(indnan(m))-blon(indwet))+ ...
                   (blat(indnan(m))-blat(indwet)).*(blat(indnan(m))-blat(indwet)));
            near=find(d == min(d));
            if ~isempty(near),
              Fobc(indnan(m))=Fobc(indwet(near));
            end,
          end,
        else,
          indnan=find(bmask > 0 & isnan(Fobc) == 1);
          indwet=find(rmask > 0 & isnan(Fwrk) == 0);
          icnan=icnan+length(indnan);
          for m=1:length(indnan),
            d=sqrt((blon(indnan(m))-rlon(indwet)).*(blon(indnan(m))-rlon(indwet))+ ...
                   (blat(indnan(m))-rlat(indwet)).*(blat(indnan(m))-rlat(indwet)));
            near=find(d == min(d));
            if ~isempty(near),
              Fobc(indnan(m))=Fwrk(indwet(near));
            end,
          end,
        end,

        switch ib
          case {1}
            if (Vname == 'u' | Vname == 'v'),
              Flev.west(:,:,k)=Fobc;
            else,
              Flev.west(:,k)=Fobc;
            end,
          case {2}
            if (Vname == 'u' | Vname == 'v'),
              Flev.east(:,:,k)=Fobc;
            else,
              Flev.east(:,k)=Fobc;
            end,
          case {3}
            if (Vname == 'u' | Vname == 'v'),
              Flev.south(:,:,k)=Fobc;
            else,
              Flev.south(:,k)=Fobc;
            end,
          case {4}
            if (Vname == 'u' | Vname == 'v'),
              Flev.north(:,:,k)=Fobc;
            else,
              Flev.north(:,k)=Fobc;
            end,
        end,

      end,
    end,
  end,
  

%  Vertical interpolation into ROMS grid. Notice that Mercator depths
%  are positive. Notice that "u" and "v" velocities are interpolated
%  at two RHO-points next to the boundary to facilitate rotation and
%  average to their respective C-grid locations.

  [Im,Jm,N]=size(z_r);
  depth=-abs(depth(1:Nlev));
  depth(1)=0.0;                % reset surface level to zero (0.49 orinally)
  
  for ib=1:4,

    if (S.boundary(ib)),

      switch ib
        case {1}
          if (Vname == 'u' | Vname == 'v'),
            bmask=rmask(1:2,:);
            z=squeeze(z_r(1:2,:,:));
            Fout.west=zeros(size(z));
            for j=1:Jm,
              for i=1:2,
                if (bmask(i,j) > 0),
                  Fout.west(i,j,:)=interp1(depth,squeeze(Flev.west(i,j,:)),squeeze(z(i,j,:)));
                end,
              end,
            end,
            indwet=find((repmat(bmask,[1,1,N])) > 0);
            Bmin(ib)=min(min(Fout.west(indwet)));
            Bmax(ib)=max(max(Fout.west(indwet)));
          else,
            bmask=rmask(1,:);
            z=squeeze(z_r(1,:,:));
            Fout.west=zeros(size(z));
            for j=1:Jm,
              if (bmask(j) > 0),
                Fout.west(j,:)=interp1(depth,squeeze(Flev.west(j,:)),squeeze(z(j,:)));
              end,
            end,
            indwet=find((repmat(bmask,[1,N])) > 0);
            Bmin(ib)=min(min(Fout.west(indwet)));
            Bmax(ib)=max(max(Fout.west(indwet)));
          end,
        case {2}
          if (Vname == 'u' | Vname == 'v'),
            bmask=rmask(end-1:end,:);
            z=squeeze(z_r(end-1:end,:,:));
            Fout.east=zeros(size(z));
            for j=1:Jm,
              for i=1:2,
                if (bmask(i,j) > 0),
                  Fout.east(i,j,:)=interp1(depth,squeeze(Flev.east(i,j,:)),squeeze(z(i,j,:)));
                end,
              end,
            end,
            indwet=find((repmat(bmask,[1,1,N])) > 0);
            Bmin(ib)=min(min(Fout.east(indwet)));
            Bmax(ib)=max(max(Fout.east(indwet)));
          else,
            bmask=rmask(end,:);
            z=squeeze(z_r(end,:,:));
            Fout.east=zeros(size(z));
            for j=1:Jm,
              if (bmask(j) > 0),
                Fout.east(j,:)=interp1(depth,squeeze(Flev.east(j,:)),squeeze(z(j,:)));
              end,
            end,
            indwet=find((repmat(bmask,[1,N])) > 0);
            Bmin(ib)=min(min(Fout.east(indwet)));
            Bmax(ib)=max(max(Fout.east(indwet)));
          end,
        case {3}
          if (Vname == 'u' | Vname == 'v'),
            bmask=rmask(:,1:2);
            z=squeeze(z_r(:,1:2,:));
            Fout.south=zeros(size(z));
            for j=1:2,
              for i=1:Im,
                if (bmask(i,j) > 0),
                  Fout.south(i,j,:)=interp1(depth,squeeze(Flev.south(i,j,:)),squeeze(z(i,j,:)));
                end,
              end,
            end,
            indwet=find((repmat(bmask,[1,1,N])) > 0);
            Bmin(ib)=min(min(Fout.south(indwet)));
            Bmax(ib)=max(max(Fout.south(indwet)));
          else,
            bmask=rmask(:,1);
            z=squeeze(z_r(:,1,:));
            Fout.south=zeros(size(z));
            for i=1:Im,
              if (bmask(i) > 0),
                Fout.south(i,:)=interp1(depth,squeeze(Flev.south(i,:)),squeeze(z(i,:)));
              end,
            end,
            indwet=find((repmat(bmask,[1,N])) > 0);
            Bmin(ib)=min(min(Fout.south(indwet)));
            Bmax(ib)=max(max(Fout.south(indwet)));
          end,
        case {4}
          if (Vname == 'u' | Vname == 'v'),
            bmask=rmask(:,end-1:end);
            z=squeeze(z_r(:,end-1:end,:));
            Fout.north=zeros(size(z));
            for j=1:2,
              for i=1:Im,
                if (bmask(i,j) > 0),
                  Fout.north(i,j,:)=interp1(depth,squeeze(Flev.north(i,j,:)),squeeze(z(i,j,:)));
                end,
              end,
            end,
            indwet=find((repmat(bmask,[1,1,N])) > 0);
            Bmin(ib)=min(min(Fout.north(indwet)));
            Bmax(ib)=max(max(Fout.north(indwet)));
          else,
            bmask=rmask(:,end);
            z=squeeze(z_r(:,end,:));
            Fout.north=zeros(size(z));
            for i=1:Im,
              if (bmask(i) > 0),
                Fout.north(i,:)=interp1(depth,squeeze(Flev.north(i,:)),squeeze(z(i,:)));
              end,
            end,
            indwet=find((repmat(bmask,[1,N])) > 0);
            Bmin(ib)=min(min(Fout.north(indwet)));
            Bmax(ib)=max(max(Fout.north(indwet)));
          end,
      end,
    
    end,
  end,
end,

disp([ 'Interpolated  ', sprintf('%4s',Vname), ...
       '  number of NaNs replaced with nearest value = ',num2str(icnan)]);
disp([ '                    Min=', sprintf('%12.5e',min(Bmin)), ...
       '  Max=', sprintf('%12.5e',max(Bmax))]);

return
