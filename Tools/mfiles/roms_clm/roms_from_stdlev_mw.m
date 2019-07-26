function roms = roms_from_stdlev_mw(lon,lat,zlev,data,grd,CgridPos,do_interp2)
% Interpolate a 3D gridded (e.g. climatology or synthetic t(z),s(z) analysis
% to a 3D ROMS grid
%
% roms = roms_from_stdlev(lon,lat,zlev,data,grd,CgridPos,do_interp2)
%
% Inputs: 
%         standard level data:
%            data is dimensioned (time,zlev,x,y) 
%            lon is (x,y)
%            lat is (x,y) [like from meshgrid]
%            zlev(z) vector of standard levels; zlev(1) being the deepest 
%               if isempty(zlev) the interpolation is assumed to be for 
%               a 2D (e.g. surface or depth-averaged) dataset, and 
%               interpolation is to the horizontal grid only
%            grd = roms grid structure
%            CgridLocation = target grid 'rho' (default), 'u' or 'v'
%
% Outputs: 
%         roms = data interpolated to the ROMS grid, preserving the time
%            dimension of the input data
%
% John Wilkin Dec 2000

if nargin < 6
  CgridPos = 'rho';
end

switch CgridPos 
  case 'rho'    
    romslon = grd.lon_rho;
    romslat = grd.lat_rho;
    romsmask = grd.mask_rho;

  case 'u'    
    romslon = grd.lon_u;
    romslat = grd.lat_u;
    romsmask = grd.mask_u;
    
  case 'v'    
    romslon = grd.lon_v;
    romslat = grd.lat_v;
    romsmask = grd.mask_v;

  otherwise    
    error('Input CgridPos must be ''rho'', ''u'' or ''v'' ')    
end

land = find(romsmask==0);

% preallocate arrays for data -> roms_stdlev([time],stdlev,lat,lon)
% preallocate arrays for data -> roms([time],scoord,lat,lon)

roms = [];

if isempty(zlev)

  if ndims(data) == 3 
    Nt = size(data,1);
    roms_stdlev = zeros([Nt size(romslon)]);
  
  elseif ndims(data) == 2;
    Nt = 1;
    roms_stdlev = zeros([size(romslon)]);
  
  else
    error('Check the size of ''data'' ')
  end

else

  if ndims(data) == 4
    Nt = size(data,1);
    Nz = size(data,2);
    roms_stdlev = zeros([Nt Nz size(romslon)]);
    Ns = size(grd.z_r,1);
    roms = zeros([Nt Ns size(romslon)]);
  
  elseif ndims(data) == 3;
    Nt = 1;%number of time steps
    Nz = size(data,1);%number of z levels
    roms_stdlev = zeros([Nz size(romslon)]);
    Ns = size(grd.z_r,1);
    roms = zeros([Ns size(romslon)]);    
  
  elseif ndims(data) == 2;
    Nt = 1;
    Nz = 1;
    roms_stdlev = zeros([size(romslon)]);
  
  else
    error('Check the size of ''data'' ')
  end
  
  % force zlev to be negative values to be consistent with ROMS
  % convention (i.e. trap the case of z given as standard depths (>0) 
  zlev = -abs(zlev(:));
  
  if Nz>1    
    % check that the data is ordered from deep to shallow
    if any(diff(zlev)<0)
      % depths are arranged from shallowest to deepest, so flip
      warning([ ' reversing zlev to arrange as deep to shallow']);
      zlev = flipud(zlev);
      data = flipdim(data,ndims(data)-2);  % check ndims to get right z index 
    end
  end
      
end

% INTERPOLATE TO HORIZONTAL GRID ----------------------------------------------
% -----------------------------------------------------------------------------


if nargin < 7
  do_interp2 = 1;
end

if do_interp2

  disp('Interpolating to the ROMS horizontal grid at standard levels') % -----
  switch ndims(roms_stdlev)
    
    case 4
      for l = 1:Nt
        disp([' Doing time ' int2str(l) ' of ' int2str(Nt)])
        for k = 1:Nz
          datawrk = squeeze(data(l,k,:,:));
          tmp = interp2(lon,lat,datawrk,romslon,romslat,'spline');
          if ~isempty(land),tmp(land)=0;end
          roms_stdlev(l,k,:,:) = tmp;
        end
      end
      
    case 3
      for k = 1:max([Nt Nz])
        datawrk = squeeze(data(k,:,:));
        try
          tmp = interp2(lon,lat,datawrk,romslon,romslat,'spline');
          if ~isempty(land),tmp(land)=0;end
          roms_stdlev(k,:,:) = tmp;
        catch
          tmp = griddata(lon,lat,datawrk,romslon,romslat);
          if ~isempty(land),tmp(land)=0;end
          roms_stdlev(k,:,:) = tmp;
        end
      end
      
    case 2
      datawrk = data;
      try
        tmp = interp2(lon,lat,datawrk,romslon,romslat,'spline');
        if ~isempty(land),tmp(land)=0;end
        roms_stdlev(:,:) = tmp;
      catch
        tmp = griddata(lon,lat,datawrk,romslon,romslat);
        if ~isempty(land),tmp(land)=0;end
        roms_stdlev(:,:) = tmp;
      end
      
  end
  
else
  
  disp('Input is assumed to be on the ROMS horizontal rho grid') % ------------
  switch CgridPos
    
    case 'rho'
      roms_stdlev = data;
      
    case 'u'
      disp(' but is being averaged to the u-points grid ')
      [LL,MM] = size(romslon);
      switch ndims(roms_stdlev)
        case 4
          roms_stdlev = 0.5*(data(:,:,1:LL,:)+data(:,:,2:LL+1,:));
        case 3
          roms_stdlev = 0.5*(squeeze(data(:,1:LL,:))+squeeze(data(:,2:LL+1,:)));
        case 2
          roms_stdlev = 0.5*(data(1:LL,:)+data(2:LL+1,:));
      end
      
    case 'v'
      disp(' but is being averaged to the v-points grid ')
      [LL,MM] = size(romslon);
      switch ndims(roms_stdlev)
        case 4
          roms_stdlev = 0.5*(data(:,:,:,1:MM)+data(:,:,:,2:MM+1));
        case 3
          roms_stdlev = 0.5*(data(:,:,1:MM)+data(:,:,2:MM+1));
        case 2
          roms_stdlev = 0.5*(data(:,1:MM)+data(:,2:MM+1));
      end
      
  end
  
end

% INTERPOLATE TO VERTICAL GRID ------------------------------------------------
% -----------------------------------------------------------------------------

if isempty(roms)
  
  disp('No vertical interpolation is required')
  roms = roms_stdlev;
  
else
  
  disp('Interpolating to ROMS s-coordinates') % -------------------------------
  
  % if necessary, average the roms z_r to velocity points
  z_ = grd.z_r;
  switch CgridPos
    case 'u'
      L = size(z_,2);
      z_ = 0.5*(z_(:,1:L-1,:)+z_(:,2:L,:));
    case 'v'
      M = size(z_,3);
      z_ = 0.5*(z_(:,:,1:M-1)+z_(:,:,2:M));
  end
  
  switch ndims(roms)
    
    case 4
      for l=1:Nt
        disp([' Doing time ' int2str(l) ' of ' int2str(Nt)])
        
        % interpolate a y-z plane each time
        Nx = size(roms,3);
        for i=1:Nx % x index
          
          if ~rem(i,20)
            disp(['  Doing i = ' int2str(i) ' of ' int2str(Nx)])
          end
          
          z = squeeze(z_(:,i,:));
          M = size(z,2);
          x = repmat(1:M,[Ns 1]);
          
          % There may be ROMS z values outside the stdlev z range, so pad
          % above and below before interp2 (just like in roms_zslice)
          % (Hmmm ... there may still be a catch if there are some very
          % deep depths with NaNs in the data)
          [xa,za] = meshgrid(1:M,[-10000; -abs(zlev); 10]);
          
          data = squeeze(roms_stdlev(l,:,i,:));
          data = [data(1,:); data; data(Nz,:)];
          roms(l,:,i,:) = interp2(xa,za,data,x,z,'spline');
          
        end
      end
      
    case 3
    clear lat lon romslat romslon romsmask   
      % interpolate a y-z plane each time
      Nx = size(roms,2);
      for i=1:Nx % x index
        
        if ~rem(i,20)
          disp(['  Doing i = ' int2str(i) ' of ' int2str(Nx)])
        end
        
        z = squeeze(z_(:,i,:));
        M = size(z,2);
        x = repmat(1:M,[Ns 1]);
        
        % There may be ROMS z values outside the stdlev z range, so pad
        % above and below before interp2 (just like in roms_zslice)
        % (Hmmm ... there may still be a catch if there are some very
        % deep depths with NaNs in the data)
        [xa,za] = meshgrid(1:M,[-10000; -abs(zlev); 10]);
        
        data = squeeze(roms_stdlev(:,i,:));
        data = [data(1,:); data; data(Nz,:)];
        data((isnan(data)==1))=0;
        roms(:,i,:) = interp2(xa,za,data,x,z,'spline');
        
      end
  end
  
end


