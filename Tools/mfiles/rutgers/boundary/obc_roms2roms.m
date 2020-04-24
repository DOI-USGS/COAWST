function B = obc_roms2roms(ncfile,D,R,VarList,Tindex,boundary,varargin)

%
% OBC_ROMS2ROMS: Interpolates ROMS LBC variable to specified ROMS grid
%
% B = obc_roms2roms(ncfile,D,R,VarList,Tindex,boundary, ...
%                   Hmethod,offset,RemoveNaN);
%
% This function interpolates lateral boundary conditions (LBC) variables
% between two ROMS application grids. The receiver grid must be inside
% of the donor grid.
%
% This function is intended for down-scaling or nesting applications.
% The horizontal/vertical coordinates for the donor and the receiver
% grids are specified with array structures 'D' and 'R', which are
% builded elsewhere using script 'get_roms_grid.m' for efficiency
% and functionality. It uses 'TriScatteredInterp' for interpolating
% lateral boundary variables.
% 
% On Input:
%
%    ncfile        Donor NetCDF file/URL name (string) containing
%                    variable to process
%
%    D             Donor grid structure containing all horizontal
%                    and vertical variables (struct array)
%
%    R             Receiver grid structure containing all horizontal
%                    and vertical variables (struct array)
%
%    VarList       List of variable names to process (cell array)
%
%    Tindex        Time record index to process (scalar).
%
%    boundary      Lateral boundary condition switches of the grid
%                    edges to process (struct array)
%
%                    boundary.west        Western  edge
%                    boundary.east        Eastern  edge
%                    boundary.south       Southern edge
%                    boundary.north       Northern edge
%
%    Hmethod       Horizontal Interpolation method in 'TriScatteredInterp'
%                    (string):
%
%                    'natural'     natural neighbor interpolation
%                    'linear'      linear interpolation (default)
%                    'nearest'     nearest-neighbor interpolation
%
%    offset        Number of extra points to use when sampling the
%                    donor grid so it is large enough to contain
%                    the receiver grid  (default 5)
%
%    RemoveNaN     Switch to remove NaN values from interpolated 
%                    variable with a second interpolation step
%                    using the nearest-neighbor method
%                    (default false)
%
% On Output:
%
%    B             Interpolated requested 2D/3D variable boundary
%                    conditions (struct array)
%
 
% svn $Id: obc_roms2roms.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%  Set optional arguments.

Hmethod = 'linear';
Vmethod = 'linear';
offset = 5;
RemoveNaN = false;

switch numel(varargin)
  case 1
    Hmethod = varargin{1};
  case 2
    Hmethod = varargin{1};
    offset = varargin{2};
  case 3
    Hmethod = varargin{1};
    offset = varargin{2};
    RemoveNaN = varargin{3};
end

%  Check if 'vector_rotation' field is available in the donor and
%  receiver Grid structures.Determine. Vector variables require special
%  processing.

if (~isfield(D,'vector_rotation'))
  if (isfield(D,'angle'))
    if (min(D.angle(:)) == 0 && max(D.angle(:)) == 0)
      D.vector_rotation = false;
    else
      D.vector_rotation = true;
    end
  else
    D.vector_rotation = false;
  end
end

if (~isfield(R,'vector_rotation'))
  if (isfield(R,'angle'))
    if (min(R.angle(:)) == 0 && max(R.angle(:)) == 0)
      R.vector_rotation = false;
    else
      R.vector_rotation = true;
    end
  else
    R.vector_rotation = false;
  end
end

Ucomponent = false;
Vcomponent = false;

%==========================================================================
%  Process every field in the input variable list.
%==========================================================================

icount = 0;

for var = VarList

  Vname  = char(var);
  icount = icount+1;

%  Initialize.

  got.Mname = false;
  got.Xname = false;
  got.Yname = false;
  got.Zname = false;

  isr3d = false;
  isw3d = false;
  isvec = false;

  Ncount = 0;

%  Get information about variable to process.

  Info = nc_vinfo(ncfile,Vname);

  nvdims = length(Info.Dimensions);

%  Check variable dimensions and determine horizontal/vertical
%  coordinates and Land/Sea mask arrays.

  if (nvdims > 0)
    for n=1:nvdims
      dimnam = char(Info.Dimensions(n).Name);
      switch dimnam
        case 's_rho'
          isr3d = true;
        case 's_w'
          isw3d = true;
        case {'xi_rho','eta_rho'}
          Mname = 'mask_rho';
          got.Mname = true;
          if (~(got.Xname || got.Yname))
            if (D.spherical)
              Xname = 'lon_rho';
              Yname = 'lat_rho';
            else
              Xname = 'x_rho';
              Yname = 'y_rho';
            end
            Zname = 'z_r';
            got.Xname = true;
            got.Yname = true;
            got.Zname = true;
          end
        case {'xi_psi','eta_psi'}
          Mname = 'mask_psi';
          got.Mname = true;
          if (~(got.Xname || got.Yname))
            if (D.spherical)
              Xname = 'lon_psi';
              Yname = 'lat_psi';
            else
              Xname = 'x_psi';
              Yname = 'y_psi';
            end
            got.Xname = true;
            got.Yname = true;       
          end
        case {'xi_u','eta_u'}
          Mname = 'mask_u';
          got.Mname = true;
          Ucomponent = true;
          if (~(got.Xname || got.Yname))
            if (D.spherical)
              Xname = 'lon_u';
              Yname = 'lat_u';
            else
              Xname = 'x_u';
              Yname = 'y_u';
            end 
            Zname = 'z_u';
            got.Xname = true;
            got.Yname = true;        
            got.Zname = true;
          end
          isvec = true;
        case {'xi_v','eta_v'}
          Mname = 'mask_v';
          got.Mname = true;
          Vcomponent = true;
          if (~(got.Xname || got.Yname))
            if (D.spherical)
              Xname = 'lon_v';
              Yname = 'lat_v';
            else
              Xname = 'x_v';
              Yname = 'y_v';
            end
            Zname = 'z_v';
            got.Xname = true;
            got.Yname = true;        
            got.Zname = true;
          end
          isvec = true;
      end
    end
    if (isw3d)
      Zname = 'z_w';
    end  
  end

  is3d = isr3d || isw3d;

%--------------------------------------------------------------------------
%  Get horizontal and vertical coordinates from donor and receiver grids.
%--------------------------------------------------------------------------

%  Donor grid.

  if (isfield(D,Xname))
    if (~isempty(D.(Xname)))
      XD = D.(Xname);
    else
      error([' OBC_ROMS2ROMS - field '', Xname, ''',                    ...
             ' is empty in donor grid structure: D']);
    end
  else
    error([' OBC_ROMS2ROMS - unable to find field '', Xname, ''',       ...
           ' in donor grid structure: D']);
  end

  if (isfield(D,Yname))
    if (~isempty(D.(Yname)))
      YD = D.(Yname);
    else
      error([' OBC_ROMS2ROMS - field '', Yname, ''',                    ...
             ' is empty in donor grid structure: D']);
    end
  else
    error([' OBC_ROMS2ROMS - unable to find field '', Yname, ''',       ...
           ' in donor grid structure: D']);
  end

  if (is3d)
    if (isfield(D,Zname))
      if (~isempty(D.(Zname)))
        ZD = D.(Zname);
      else
        error([' OBC_ROMS2ROMS - field '', Zname, ''',                  ...
               ' is empty in donor grid structure: D']);
      end
    else
      error([' OBC_ROMS2ROMS - unable to find field '', Zname, ''',     ...
             ' in donor grid structure: D']);
    end
  end

  if (isfield(D,Mname))
    if (~isempty(D.(Mname)))
      Dmask = D.(Mname);
    else
      error([' OBC_ROMS2ROMS - field '', Mname, ''',                    ...
             ' is empty in donor grid structure: D']);
    end
  else
    error([' OBC_ROMS2ROMS - unable to find field '', Mname, ''',       ...
           ' in donor grid structure: D']);
  end

%  Receiver grid.

  if (isvec && R.vector_rotation)
    if (D.spherical)
      Xname = 'lon_rho';           % If requested, interpolate
      Yname = 'lat_rho';           % U- and V-points variables
    else                           % to RHO-points instead to
      Xname = 'x_rho';             % facilitate the curvilinear
      Yname = 'y_rho';             % rotation to true East and
    end                            % North.  This rotation is
    Mname = 'mask_rho';            % done somewhere else.
    if (is3d)
      Zname = 'z_r';
    end
  end

  if (isfield(R,Xname))
    if (~isempty(R.(Xname)))  
      XR = R.(Xname);
    else
      error([' OBC_ROMS2ROMS - field '', Xname, ''',                    ...
             ' is empty in receiver grid structure: R']);
    end
  else
    error([' OBC_ROMS2ROMS - unable to find field '', Xname, ''',       ...
           ' in receiver grid structure: R']);
  end

  if (isfield(R,Yname))
    if (~isempty(R.(Yname)))
      YR = R.(Yname);
    else
      error([' OBC_ROMS2ROMS - field '', Yname, ''',                    ...
             ' is empty in receiver grid structure: R']);
    end
  else
    error([' OBC_ROMS2ROMS - unable to find field '', Yname, ''',       ...
           ' in receiver grid structure: R']);
  end

  if (is3d)
    if (isfield(R,Zname))
      if (~isempty(R.(Zname)))
        ZR = R.(Zname);
      else
        error([' OBC_ROMS2ROMS - field '', Zname, ''',                  ...
               ' is empty in receiver grid structure: R']);
      end
    else
      error([' OBC_ROMS2ROMS - unable to find field '', Zname, ''',     ...
             ' in receiver grid structure: R']);
    end
  end
  
  if (isfield(R,Mname))
    if (~isempty(R.(Mname)))
      Rmask = R.(Mname);
    else
      error([' OBC_ROMS2ROMS - field '', Mname, ''',                    ...
             ' is empty in receiver grid structure: R']);
    end
  else
    error([' OBC_ROMS2ROMS - unable to find field '', Mname, ''',       ...
           ' in receiver grid structure: R']);
  end

%--------------------------------------------------------------------------
%  Read in requested variable from donor NetCDF file.
%--------------------------------------------------------------------------

  ReplaceValue = NaN;
  PreserveType = false;

  VD = nc_read(ncfile,Vname,Tindex,ReplaceValue,PreserveType);

%--------------------------------------------------------------------------
%  Set donor grid sampling indices to accelerate the interpolation.
%  The horizontally sampled donor grid is large enough to contain
%  the receiver grid. The parameter 'offset' is used to add extra
%  points when computing the sampling indices (Istr:Iend,Jstr:Jend).
%  That is, the sampled grid is 'offset' points larger in all sides.
%  This is done to resolve well the interpolation near the boundaries
%  of the receiver grid.
%--------------------------------------------------------------------------

  [Istr,Iend,Jstr,Jend] = sample_grid(XD,YD,XR,YR,offset);
  
%--------------------------------------------------------------------------
%  Set receiver grid lateral boundary conditions locations.  If donor or
%  receiver grids are curvilinear, extract two 
%--------------------------------------------------------------------------

  if (isvec && (D.vector_rotation || R.vector_rotation))
    RX.west  = XR(1:2,:);
    RY.west  = YR(1:2,:);
    RM.west  = Rmask(1:2,:);

    RX.east  = XR(end-1:end,:);
    RY.east  = YR(end-1:end,:);
    RM.east  = Rmask(end-1:end,:);

    RX.south = XR(:,1:2);
    RY.south = YR(:,1:2);
    RM.south = Rmask(:,1:2);

    RX.north = XR(:,end-1:end);
    RY.north = YR(:,end-1:end);
    RM.north = Rmask(:,end-1:end);

    if (is3d)
      RZ.west  = ZR(1:2,:,:);
      RZ.east  = ZR(end-1:end,:,:);
      RZ.south = ZR(:,1:2,:);
      RZ.north = ZR(:,end-1:end,:);
    end
  else
    RX.west  = XR(1,:);
    RY.west  = YR(1,:);
    RM.west  = Rmask(1,:);

    RX.east  = XR(end,:);
    RY.east  = YR(end,:);
    RM.east  = Rmask(end,:);

    RX.south = XR(:,1);
    RY.south = YR(:,1);
    RM.south = Rmask(:,1);

    RX.north = XR(:,end);
    RY.north = YR(:,end);
    RM.north = Rmask(:,end);

    if (is3d)
      RZ.west  = squeeze(ZR(1,:,:));
      RZ.east  = squeeze(ZR(end,:,:));
      RZ.south = squeeze(ZR(:,1,:));
      RZ.north = squeeze(ZR(:,end,:));
    end
  end
  
%--------------------------------------------------------------------------
%  Interpolate lateral boundary conditions to receiver grid: Build
%  interpolation data structure, I.
%--------------------------------------------------------------------------
   
  I.VarList  = VarList;
  I.Vname    = Vname;
  I.nvdims   = nvdims-1;
  I.boundary = boundary;

%  Determine if processing 2D or 3D ROMS state variables.

  switch (nvdims-1)
 
    case 2

      [ImD,JmD]=size(XD);
      [ImR,JmR]=size(XR);

      disp(' ');
      disp(['Interpolating 2D OBC variable(s): ', Vname,                ...
            ' (', num2str(ImR), 'x', num2str(JmR),') from donor ',      ...
            '(', num2str(ImD), 'x', num2str(JmD),') ...']);
      disp(' ');

      I.VD    = VD(Istr:1:Iend,Jstr:1:Jend);
      
      I.Dmask = Dmask(Istr:1:Iend,Jstr:1:Jend);
      I.XD    = XD(Istr:1:Iend,Jstr:1:Jend);
      I.YD    = YD(Istr:1:Iend,Jstr:1:Jend);
      I.ZD    = [];

      I.Rmask = RM;
      I.XR    = RX;
      I.YR    = RY;
      I.ZR    = [];

      I.Zsur  = [];
      I.Zbot  = [];
      
    case 3

      [ImD,JmD,KmD]=size(ZD);
      [ImR,JmR,KmR]=size(ZR);

      disp(' ');
      disp(['Interpolating 3D OBC variable(s): ', Vname,                ...
            ' (', num2str(ImR), 'x', num2str(JmR), 'x',                 ...
                  num2str(KmR), ') from donor ',                        ...
             '(', num2str(ImD), 'x', num2str(JmD), 'x',                 ...
                  num2str(KmD), ') ...']);
      disp(' ');
  
      I.VD    = VD(Istr:1:Iend,Jstr:1:Jend,:);
      
      I.Dmask = Dmask(Istr:1:Iend,Jstr:1:Jend);
      I.XD    = XD(Istr:1:Iend,Jstr:1:Jend);
      I.YD    = YD(Istr:1:Iend,Jstr:1:Jend);
      I.ZD    = ZD(Istr:1:Iend,Jstr:1:Jend,:);

      I.Rmask = RM;
      I.XR    = RX;
      I.YR    = RY;
      I.ZR    = RZ;

      I.Zsur  = max(R.z_w(:))+eps;
      I.Zbot  = min(R.z_w(:))-eps;
  end

%  Interpolate lateral boundaries for requested state variable.

  S = interp_boundary(I,Hmethod,Vmethod,RemoveNaN);

%  Concatenate structures.

  if (icount == 1)
    B = S;
  else
    B = cell2struct(cat(1, struct2cell(B), struct2cell(S)),             ...
                    cat(1, fieldnames(B),  fieldnames(S)), 1);
  end

%  Process next state variable in the list.  Clear interpolarion
%  structure from memory.

  clear I

end

%--------------------------------------------------------------------------
%  Rotate vector components.
%--------------------------------------------------------------------------

%  Vector components that require rotation were interpolated at
%  RHO-points.  They are interpolated over two points adjacent to
%  the boundary edge to allow averaging to the appropriate C-grid
%  location.

if ((Ucomponent && Vcomponent) &&                                       ...
    (D.vector_rotation || R.vector_rotation))
  B = rotate_vectors(B,D,R,VarList,boundary);
end

return

function B = rotate_vectors(B,D,R,VarList,boundary)

% This function rotates vector components to receiver application grid.
% In order to allow rotation, the vector components were interpolated
% at RHO-points. There is data available over the two points adjacent
% (columns/rows) to the receiver application grid boundary.  After the
% rotations are carried, the vector components are averaged to their
% respective C-grid locations.
%
% This function assumes that the 'parent_angle' field is in the
% Receiver Grid structure, R, and was interpolated from the Donor
% Grid elsewhere.
%
% On Input:
%
%    B             Lateral boundary conditions structure containing
%                    interpolated variables (struct array)
%
%    D             Donor grid structure containing all horizontal
%                    and vertical variables (struct array)
%
%    R             Receiver grid structure containing all horizontal
%                    and vertical variables (struct array)
%
%    VarList       List of variables names to process (cell array)
%
%    boundary      Lateral boundary condition switches of the grid
%                    edges to process (struct array)
%
%                    boundary.west        Western  edge
%                    boundary.east        Eastern  edge
%                    boundary.south       Southern edge
%                    boundary.north       Northern edge
%
% On Output:
%
%    Bnew          Update lateral boundary conditions structure
%                    (struct array)
%

%  Set receiver grid dimensions at RHO-points.

Lr = R.Lm + 2;
Mr = R.Mm + 2;
N  = R.N;

%--------------------------------------------------------------------------
%  Rotate 2D momentum vector components.
%--------------------------------------------------------------------------

if (max(strcmp(VarList,'ubar')) && max(strcmp(VarList,'vbar')))

%  Donor Grid Orientation: rotate from (XI,ETA) coordinates to
%                           TRUE East/North.

  if (D.vector_rotation)
    angle.west  = R.parent_angle(1:2,:);
    angle.east  = R.parent_angle(end-1:end,:);
    angle.south = R.parent_angle(:,1:2);
    angle.north = R.parent_angle(:,end-1:end);

    for var = {'west','east','south','north'}
      edge = char(var);
      ufield = strcat('ubar','_',edge);
      vfield = strcat('vbar','_',edge);
      if (boundary.(edge))
        if (isfield(B,ufield) && isfield(B,vfield))
          ubar = B.(ufield);
          vbar = B.(vfield);

          B.(ufield) = ubar .* cos(angle.(edge)) -                      ...
                       vbar .* sin(angle.(edge));
          B.(vfield) = vbar .* cos(angle.(edge)) +                      ...
                       ubar .* sin(angle.(edge));
        else
          error([' Unable to find fields: ',ufield,' or ',vfield,       ...
                 ' in structure array ''B''']);
        end
      end
    end
  end

%  Receiver Grid Orientation: rotate from TRUE East/North to
%                           (XI,ETA) coordinates.

  if (R.vector_rotation)
    angle.west  = R.angle(1:2,:);
    angle.east  = R.angle(end-1:end,:);
    angle.south = R.angle(:,1:2);
    angle.north = R.angle(:,end-1:end);

    for var = {'west','east','south','north'}
      edge = char(var);
      ufield = strcat('ubar','_',edge);
      vfield = strcat('vbar','_',edge);
      if (boundary.(edge))
        if (isfield(B,ufield) && isfield(B,vfield))
          ubar = B.(ufield);
          vbar = B.(vfield);

          B.(ufield) = ubar .* cos(angle.(edge)) +                      ...
                       vbar .* sin(angle.(edge));
          B.(vfield) = vbar .* cos(angle.(edge)) -                      ...
                       ubar .* sin(angle.(edge));
        else
          error([' Unable to find fields: ',ufield,' or ',vfield,       ...
                 ' in structure array ''B''']);
        end
      end
    end
  end

%  Average vector components to U-point and V-points locations and
%  apply Lan/Sea mask.

  if (D.vector_rotation || R.vector_rotation)
    umask.west  = R.mask_u(1,:);
    umask.east  = R.mask_u(end,:);
    umask.south = R.mask_u(:,1);
    umask.north = R.mask_u(:,end);

    vmask.west  = R.mask_v(1,:);
    vmask.east  = R.mask_v(end,:);
    vmask.south = R.mask_v(:,1);
    vmask.north = R.mask_v(:,end);
    
    for var = {'west','east','south','north'}
      edge = char(var);
      ufield = strcat('ubar','_',edge);
      vfield = strcat('vbar','_',edge);
      if (boundary.(edge))
        ubar = B.(ufield);
        vbar = B.(vfield);
        switch edge
          case {'west','east'}
            B.(ufield) = 0.5 .* (ubar(1,:) +                            ...
                                 ubar(2,:));
            B.(vfield) = 0.5 .* (vbar(1,1:Mr-1) +                       ...
                                 vbar(1,2:Mr  ));

            B.(ufield) = B.(ufield) .* umask.(edge);
            B.(vfield) = B.(vfield) .* vmask.(edge);
          case {'south','north'}
            B.(ufield) = 0.5.*(ubar(1:Lr-1,1)+                          ...
                               ubar(2:Lr  ,1));
            B.(vfield) = 0.5.*(vbar(:,1)+                               ...
                               vbar(:,2));

            B.(ufield) = B.(ufield) .* umask.(edge);
            B.(vfield) = B.(vfield) .* vmask.(edge);
        end
      end
    end
  end
end

%--------------------------------------------------------------------------
%  Rotate 3D momentum vector components.
%--------------------------------------------------------------------------

if (max(strcmp(VarList,'u')) && max(strcmp(VarList,'v')))


%  Donor Grid Orientation: rotate from (XI,ETA) coordinates to
%                           TRUE East/North.

  if (D.vector_rotation)
    angle.west  = repmat(R.parent_angle(1:2,:), [1,1,N]);
    angle.east  = repmat(R.parent_angle(end-1:end,:), [1,1,N]);
    angle.south = repmat(R.parent_angle(:,1:2), [1,1,N]);
    angle.north = repmat(R.parent_angle(:,end-1:end), [1,1,N]);

    for var = {'west','east','south','north'}
      edge = char(var);
      ufield = strcat('u','_',edge);
      vfield = strcat('v','_',edge);
      if (boundary.(edge))
        if (isfield(B,ufield) && isfield(B,vfield))
          u = B.(ufield);
          v = B.(vfield);

          B.(ufield) = u .* cos(angle.(edge)) -                         ...
                       v .* sin(angle.(edge));
          B.(vfield) = v .* cos(angle.(edge)) +                         ...
                       u .* sin(angle.(edge));
        else
          error([' Unable to find fields: ',ufield,' or ',vfield,       ...
                 ' in structure array ''B''']);
        end
      end
    end
  end

%  Receiver Grid Orientation: rotate from TRUE East/North to
%                             (XI,ETA) coordinates.

  if (R.vector_rotation)
    angle.west  = repmat(R.angle(1:2,:), [1,1,N]);
    angle.east  = repmat(R.angle(end-1:end,:), [1,1,N]);
    angle.south = repmat(R.angle(:,1:2), [1,1,N]);
    angle.north = repmat(R.angle(:,end-1:end), [1,1,N]);

    for var = {'west','east','south','north'}
      edge = char(var);
      ufield = strcat('u','_',edge);
      vfield = strcat('v','_',edge);
      if (boundary.(edge))
        if (isfield(B,ufield) && isfield(B,vfield))
          u = B.(ufield);
          v = B.(vfield);

          B.(ufield) = u .* cos(angle.(edge)) +                         ...
                       v .* sin(angle.(edge));
          B.(vfield) = v .* cos(angle.(edge)) -                         ...
                       u .* sin(angle.(edge));
        else
          error([' Unable to find fields: ',ufield,' or ',vfield,       ...
                 ' in structure array ''B''']);
        end
      end
    end
  end

%  Average vector components to U-point and V-points locations and
%  apply Land/Sea mask.

  if (D.vector_rotation || R.vector_rotation)
    umask.west  = squeeze(repmat(R.mask_u(1,:), [1,1,N]));
    umask.east  = squeeze(repmat(R.mask_u(end,:), [1,1,N]));
    umask.south = squeeze(repmat(R.mask_u(:,1), [1,1,N]));
    umask.north = squeeze(repmat(R.mask_u(:,end), [1,1,N]));

    vmask.west  = squeeze(repmat(R.mask_v(1,:), [1,1,N]));
    vmask.east  = squeeze(repmat(R.mask_v(end,:), [1,1,N]));
    vmask.south = squeeze(repmat(R.mask_v(:,1), [1,1,N]));
    vmask.north = squeeze(repmat(R.mask_v(:,end), [1,1,N]));
    
    for var = {'west','east','south','north'}
      edge = char(var);
      ufield = strcat('u','_',edge);
      vfield = strcat('v','_',edge);
      if (boundary.(edge))
        u = B.(ufield);
        v = B.(vfield);
        switch edge
          case {'west','east'}
            B.(ufield) = squeeze(0.5 .* (u(1,:,:) +                     ...
                                         u(2,:,:)));
            B.(vfield) = squeeze(0.5 .* (v(1,1:Mr-1,:) +                ...
                                         v(1,2:Mr  ,:)));

            B.(ufield) = B.(ufield) .* umask.(edge);
            B.(vfield) = B.(vfield) .* vmask.(edge);
          case {'south','north'}
            B.(ufield) = squeeze(0.5.*(u(1:Lr-1,1,:)+                   ...
                                       u(2:Lr  ,1,:)));
            B.(vfield) = squeeze(0.5.*(v(:,1,:)+                        ...
                                       v(:,2,:)));

            B.(ufield) = B.(ufield) .* umask.(edge);
            B.(vfield) = B.(vfield) .* vmask.(edge);
        end
      end
    end
  end
end

return
