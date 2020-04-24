function B = interp_boundary(I,varargin)

%
% INTERP_BOUNDARY:  Interpolates lateral boundary conditions for ROMS
%
% B = interp_boundary(I,Hmethod,Vmethod,RemoveNaN)
%
% This function interpolates lateral boundary conditions for a ROMS
% generic 2D or 3D state variable from a Donor to Receiver.
%  
% The horizontal interpolation is done with 'scatteredInterpolant'. If the
% 'RemoveNaN' switch is activated (true), a second interpolation pass is
% carried out with 'scatteredInterpolant' using the nearest-neighbor method
% to remove unbounded (NaN) values.
%  
% If 3D (vertical) interpolation, the Donor Grid data is interpolated first
% to the Receiver Grid horizontal locations using 'scatteredInterpolant' at
% each of the Donor Grid vertical levels.  Then,  'interp1'  is used to
% interpolate to Receiver Grid vertical locations.
%  
% Notice that it is possible for the Donor and Receiver Grids to have or
% not a vertical terrain-following coordinates distribution. The strategy
% is to duplicate the shallowest and deepest levels of data to correspond
% to bogus depths:  ZD = I.Zsur for shallowest and ZD = I.Zbot for deepest
% (depths are negative in ROMS).  This is done to facilitate an easy
% extrapolation when the Receiver Grid depths are outside the Donor Grid
% range.
%
% On Input:
%
%    I             Interpolation data (struct array):
%
%                    I.VarList   List of all state boundary variables
%                                  (Optional, cell array)
%                    I.Vname     Field variable name (string)
%                    I.nvdims    Number of variable dimensions (scalar)
%
%                    I.boundary  Lateral boundary condition switches of the
%                                   grid edges to process (struct array)
%
%                                I.boundary.west        Western  edge
%                                I.boundary.east        Eastern  edge
%                                I.boundary.south       Southern edge
%                                I.boundary.north       Northern edge
%
%                    I.VD        Donor Grid variable data (2D/3D array)
%
%                    I.Dmask     Donor Grid land/sea masking (2D array)
%                    I.XD        Donor Grid X-locations (2D array)
%                    I.YD        Donor Grid Y-locations (2D array)
%                    I.ZD        Donor Grid Z-locations (3D array)
%
%                    I.Rmask     Receiver Grid land/sea masking (2D array) 
%                    I.XR        Receiver Grid X-locations (2D array)
%                    I.YR        Receiver Grid Y-locations (2D array)
%                    I.ZR        Receiver Grid Z-locations (3D array
%
%                    I.Zsur      shallowest depth for extracpolation (scalar)
%                    I.Zbot      deepest depth for extracpolation (scalar)
%                     
%    Hmethod       Horizontal interpolation method for 'scatteredInterpolant'
%                    (string):
%
%                    'natural'   natural neighbor interpolation
%                    'linear'    linear interpolation (default)
%                    'nearest'   nearest-neighbor interpolation
%
%    Vmethod       Vertical interpolation method for 'interp1' (string):
%
%                    'nearest'   nearest neighbor interpolation
%                    'linear'    bilinear interpolation (default)
%                    'spline'    spline interpolation
%                    'cubic'     bicubic interpolation as long as the
%                                  data is uniformly spaced, otherwise
%                                  the same as 'spline'
%
%    RemoveNaN     Switch to remove NaN values from the interpolated 
%                    variable with a second interpolation step
%                    using the nearest-neighbor method
%                    (default false)
%
% On Output:
%
%    B             Interpolated lateral boundary conditions (struct array)
%

% svn $Id: interp_boundary.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license           Hernan G. Arango      %
%    See License_ROMS.txt                           John Wilkin           %
%=========================================================================%  

%  Set optional arguments.
  
Hmethod = 'linear';
Vmethod = 'linear';
RemoveNaN = false;

switch numel(varargin)
  case 1
    Hmethod = varargin{1};
  case 2
    Hmethod = varargin{1};
    Vmethod = varargin{2};
  case 3
    Hmethod = varargin{1};
    Vmethod = varargin{2};
    RemoveNaN = varargin{3};
end

%  Check interpolatiom data structure.  

hvar_list = {'Vname', 'nvdims', 'boundary',                             ...
             'VD', 'Dmask', 'XD', 'YD',                                 ...
                   'Rmask', 'XR', 'YR'};

zvar_list = {'ZD', 'ZR', 'Zsur', 'Zbot'};

for var = hvar_list
  field = char(var);
  if (~isfield(I, field))
    error(['INTERP_BOUNDARY: unable to find field: "',field,'" in ',    ...
           'input structure,  I']);
  end
end

if (I.nvdims > 2)
  for var = zvar_list
    field = char(var);
    if (~isfield(I, field))
      error(['INTERP_BOUNDARY: unable to find field: "',field,'" in ',  ...
             'input structure,  I']);
    end
  end
end

%  Set report format for all state boundary variables to allow easy
%  read of reported information

if (isfield(I,'VarList'))
  lstr = 0;
  for var = I.VarList
    lstr = max(lstr, length(char(var)));
  end
  frmt = strcat('%',num2str(lstr+6),'s');
else
  frmt = '%10s';
end

%--------------------------------------------------------------------------
%  Interpolate field variable from Donor to Receiver Grid.
%--------------------------------------------------------------------------

switch (I.nvdims)
 
  case 2                                % 2D field variable

    Ncount = 0;

    x = I.XD(:);                        % scatteredInterpolant wants 1-D
    y = I.YD(:);                        % vectors as inputs
    v = I.VD(:);

    Dind = find(I.Dmask < 0.5);    
    if (~isempty(Dind))
      x(Dind) = [];                     % remove land points, if any
      y(Dind) = [];
      v(Dind) = [];
    end
    Dmin = min(v);
    Dmax = max(v);
   
    F = scatteredInterpolant(x,y,v,Hmethod);
    if (RemoveNaN)
      FN = scatteredInterpolant(x,y,v,'nearest');
    end
    
    for side = {'west','east','south','north'}

      edge = char(side);
      field = strcat(char(I.Vname),'_',edge);

      if (I.boundary.(edge))
        B.(field) = F(I.XR.(edge), I.YR.(edge));
        Rmin = min(B.(field)(:));
        Rmax = max(B.(field)(:));
    
        Rind = find(I.Rmask.(edge) < 0.5);
        if (~isempty(Rind))
          B.(field)(Rind) = 0;
        end

%  If applicable, remove interpolated variable NaNs values with a
%  nearest neighbor interpolant.

        ind = find(isnan(B.(field)));

        if (~isempty(ind))
          if (RemoveNaN)
            B.(field)(ind) = FN(I.XR.(edge)(ind), I.YR.(edge)(ind));
            Rmin = min(Rmin, min(B.(field)(ind)));
            Rmax = max(Rmax, max(B.(field)(ind)));

            ind = find(isnan(B.(field)));
            if (~isempty(ind))
              Ncount = length(ind);
            end       
          else
            Ncount = length(ind);
          end
        end

        disp(['   ',sprintf(frmt,field),':  ',                          ...
              '   Donor Min = ', sprintf('%12.5e',Dmin), '   ',         ...
              '   Donor Max = ', sprintf('%12.5e',Dmax)]);
        disp(['   ',sprintf(frmt,field),':  ',                          ...
              'Receiver Min = ', sprintf('%12.5e',Rmin), '   ',         ...
              'Receiver Max = ', sprintf('%12.5e',Rmax), '   ',         ...
              'Nan count = ',  num2str(Ncount)]);
      end
    end
 
  case 3                                % 3D field variable 

    Kcount = 0;
  
    KmD = size(I.ZD,3);

%  First, perform horizontal interpolation using 'scatteredInterpolant'
%  at each Donor grid level.
   
    x = I.XD(:);                       % scatteredInterpolant wants 1-D
    y = I.YD(:);                       % vectors as inputs

    Dind = find(I.Dmask < 0.5);
    if ~isempty(Dind)
      x(Dind) = [];                    % remove land points, if any
      y(Dind) = [];
    end  

%  Initialize 'interpolation objects' for null data (ones).
   
    F  = scatteredInterpolant(x,y,ones(size(x)),Hmethod);
    FN = scatteredInterpolant(x,y,ones(size(x)),'nearest');
    
%  Horizontal interpolation at Donor Grid vertical levels

    Dmin = Inf;
    Dmax = -Inf;

    for k = 1:KmD
     
      vk = I.VD(:,:,k);                % Donor Grid data for level k
      vk = vk(:);

      zk = I.ZD(:,:,k);                % Donor Grid depths for level k
      zk = zk(:);
     
      if (~isempty(Dind))
        vk(Dind) = [];
        zk(Dind) = [];
      end
      Dmin = min(Dmin,min(vk));
      Dmax = max(Dmax,max(vk));
     
      F.Values  = vk;                  % place boundary data to interpolate
      FN.Values = vk;                  % level k in scatteredInterpolant
                                       % object

      for side = {'north','south','east','west'}

        edge = char(side);             % interpolate boundary data at
                                       % Receiver Grid horizontal            
        if (I.boundary.(edge))         % locations but Donor levels
          XR = I.XR.(edge);
          YR = I.YR.(edge);
          Vk = F(XR, YR);
                                       
          ind = find(isnan(Vk));
          if (~isempty(ind))
            if (RemoveNaN)
              Vk(ind) = FN(XR(ind), YR(ind));

              ind = find(isnan(Vk));
              if (~isempty(ind))
                Kcount = Kcount+length(ind);
              end       
            else
              Kcount = Kcount+length(ind);
            end
          end
          V.(edge)(:,k) = Vk(:);
        end
      end

      F.Values  = zk;                  % place depths to interpolate for
      FN.Values = zk;                  % level k in scatteredInterpolant
                                       % object

      for side = {'north','south','east','west'}

        edge = char(side);             % interpolate Donor Grid depths
                                       % at Receiver Grid horizontal          
        if (I.boundary.(edge))         % locations but Donor levels
          XR = I.XR.(edge);
          YR = I.YR.(edge);
          Zk = F(XR, YR);
                                       
          ind = find(isnan(Zk));
          if (~isempty(ind))
            if RemoveNaN
              Zk(ind) = FN(XR(ind), YR(ind));
            end
          end
          Z.(edge)(:,k) = Zk(:);
        end
      end

    end

%  Lastly, perform vertical interpolation.

    for side = {'north','south','east','west'}
      edge = char(side);
      field = strcat(I.Vname,'_',edge);

      if (I.boundary.(edge))

        Ncount = 0;
        dsize = size(I.ZR.(edge));
        KmR = dsize(end);
	
        VD = V.(edge);
        ZD = Z.(edge);
        ZR = reshape(I.ZR.(edge), [numel(I.XR.(edge)) KmR]);

%  Duplicate shallowest and deepest levels of data to correspond
%  to bogus ZD = 0 and ZR = -I.hmaxR. This forces simple vertical 
%  extrapolation if the Receiver Grid has depths outside the Donor
%  range. This will generally be the case if the receiver has more
%  levels or a stretching that emphasizes boundary layers.

        VD = VD(:,[1 1:end end]);
        ZD = [I.Zbot*ones([length(ZD) 1]),ZD,I.Zsur*ones([length(ZD) 1])];

%  Vertically interpolate.

        VR = NaN([size(ZR,1) KmR]);
        for i=1:size(ZR,1)
          VR(i,:) = interp1(ZD(i,:), VD(i,:), ZR(i,:), Vmethod);
        end
        B.(field) = squeeze(reshape(VR, dsize));

        Rmin = min(VR(:));
        Rmax = max(VR(:));

        Rind = find(repmat(I.Rmask.(edge),[1,1,KmR]) < 0.5);
        if (~isempty(Rind))
          B.(field)(Rind) = 0;
        end
	  
        ind = find(isnan(VR));
        if (~isempty(ind))
          Ncount = length(ind);
        end

        disp(['   ',sprintf(frmt,field),':  ',                          ...
              '   Donor Min = ', sprintf('%12.5e',Dmin), '   ',         ...
              '   Donor Max = ', sprintf('%12.5e',Dmax)]);
        disp(['   ',sprintf(frmt,field),':  ',                          ...
              'Receiver Min = ', sprintf('%12.5e',Rmin), '   ',         ...
              'Receiver Max = ', sprintf('%12.5e',Rmax), '   ',         ...
              'Nan count = ',  num2str(Ncount)]);
      
      end
    end

end

return
