function V = interp_field(I,varargin)

%
% INTERP_FIELD:  Interpolates 2D or 3D variable from Donor to Receiver Grid
%
% V = interp_field(I,Hmethod,Vmethod,RemoveNaN)
%
% This function interpolates a generic 2D or 3D field variable from a Donor
% to Receiver Grid.
%  
% The horizontal interpolation is done with 'TriScatteredInterp'. If the
% 'RemoveNaN' switch is activated (true), a second interpolation pass is
% carried out with 'TriScatteredInterp' using the nearest-neighbor method
% to remove unbounded (NaN) values.
%  
% If 3D interpolation, the Donor Grid data is interpolated first to the
% Receiver Grid horizontal locations using 'TriScatteredInterp' at each
% of the Donor Grid vertical levels.  Then,  'interp1'  is used to
% interpolate to Receiver Grid vertical locations.
%  
% Notice that it is possible for the Donor and Receiver Grids to have or
% not a vertical terrain-following coordinates distribution. The strategy
% is to duplicate the shallowest and deepest levels of data to correspond
% to bogus depths:  ZD = I.Zsur for shallowest and ZD = I.Zbot for deepest
% (depths are negative in ROMS).  This is done to facilitate an easy
% extrapolation when the Reciever Grid depths are outside the Donor Grid
% range.
%
% On Input:
%
%    I             Interpolation data (struct array):
%
%                    I.Vname     field variable name
%                    I.nvdims    number of variable dimensions
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
%    Hmethod       Horizontal interpolation method for 'TriScatteredInterp'
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
%    V             Interpolated 2D/3D field
%

% svn $Id: interp_field.m 996 2020-01-10 04:28:56Z arango $
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

hvar_list = {'Vname', 'nvdims', 'VD', 'Dmask', 'XD', 'YD',              ...
                                      'Rmask', 'XR', 'YR'};

zvar_list = {'ZD', 'ZR', 'Zsur', 'Zbot'};

for var = hvar_list,
  field = char(var);
  if (~isfield(I, field)),
    error(['INTERP_FIELD: unable to find field: "',field,'" in ',        ...
           'input structure,  I']);
  end
end

if (I.nvdims > 2),
  for var = zvar_list,
    field = char(var);
    if (~isfield(I, field)),
      error(['INTERP_FIELD: unable to find field: "',field,'" in ',      ...
             'input structure,  I']);
    end
  end
end

%--------------------------------------------------------------------------
%  Interpolate field variable from Donor to Receiver Grid.
%--------------------------------------------------------------------------

switch (I.nvdims),
 
 case 2                                % 2D field variable

   Ncount = 0;

   x = I.XD(:);                        % TriScatteredInterp wants 1-D
   y = I.YD(:);                        % vectors as inputs
   v = I.VD(:);

   Dind = find(I.Dmask < 0.5);    
   if (~isempty(Dind)),
     x(Dind) = [];                     % remove land points, if any
     y(Dind) = [];
     v(Dind) = [];
   end

   Dind = isnan(v);                    % remove NaN's
   if (any(Dind)),
     x(Dind) = [];
     y(Dind) = [];
     v(Dind) = [];
   end

   Dmin = min(v);
   Dmax = max(v);

   F = TriScatteredInterp(x,y,v,Hmethod);
   V = F(I.XR,I.YR);
   Rmin = min(V(:));
   Rmax = max(V(:));
    
   Rind = find(I.Rmask < 0.5);
   if (~isempty(Rind)),
     V(Rind) = 0;
   end,

%  If applicable, remove interpolated variable NaNs values with a
%  nearest neighbor interpolant.

   ind = find(isnan(V));

   if (~isempty(ind)),
     if (RemoveNaN),
       FN = TriScatteredInterp(x,y,v,'nearest');
       V(ind) = FN(I.XR(ind),I.YR(ind));
       Rmin = min(Rmin, min(V(ind)));
       Rmax = max(Rmax, max(V(ind)));

       ind = find(isnan(V));
       if (~isempty(ind)),
         Ncount = length(ind);
       end       
     else
       Ncount = length(ind);
     end
   end
   
   disp(['      Donor Min = ', sprintf('%12.5e',Dmin), '  ',            ...
         '    Donor Max = ', sprintf('%12.5e',Dmax)]);
   disp(['   Receiver Min = ', sprintf('%12.5e',Rmin), '  ',            ...
           ' Receiver Max = ', sprintf('%12.5e',Rmax), '  ',            ...
           ' NaN count = ',  num2str(Ncount)]);        
 
 case 3                                % 3D field variable

   Kcount = 0;
   Ncount = 0;
  
   [ImD,JmD,KmD]=size(I.ZD);
   [ImR,JmR,KmR]=size(I.ZR);

%  First, perform horizontal interpolation using 'TriScatteredInterp'
%  at each Donor grid level.
   
   x = I.XD(:);                        % TriScatteredInterp wants 1-D
   y = I.YD(:);                        % vectors as inputs

   Dind = find(I.Dmask < 0.5);
   if ~isempty(Dind)
     x(Dind) = [];                     % remove land points, if any
     y(Dind) = [];
   end  

%  Initialize 'TriScatteredInterp objects' for null data (ones).
   
   F  = TriScatteredInterp(x,y,ones(size(x)),Hmethod);
   FN = TriScatteredInterp(x,y,ones(size(x)),'nearest');

%  Horizontal interpolation at Donor Grid vertical level

   Dmin = Inf;
   Dmax = -Inf;

   V = NaN([size(I.XR) KmD]);          % variable data
   Z = NaN([size(I.XR) KmD]);          % Donor Grid depths at Receiver Grid

   for k = 1:KmD
     
     vk = I.VD(:,:,k);                 % Donor Grid data for level k
     vk = vk(:);

     zk = I.ZD(:,:,k);                 % Donor Grid depths for level k
     zk = zk(:);
     
     if (~isempty(Dind)),
       vk(Dind) = [];
       zk(Dind) = [];
     end
     Dmin = min(Dmin,min(vk));
     Dmax = max(Dmax,max(vk));
     
     F.V = vk;                         % place data to interpolate for
                                       % level k in TriScatteredInterp
				       % object

     Vk = F(I.XR,I.YR);                % interpolate on receiver 
                                       % horizontal grid
     
     ind = find(isnan(Vk));
     if (~isempty(ind)),
       if (RemoveNaN),
         FN.V = vk;
         Vk(ind) = FN(I.XR(ind),I.YR(ind));

         ind = find(isnan(Vk));
         if (~isempty(ind)),
           Kcount = Kcount+length(ind);
         end       
       else
         Kcount = Kcount+length(ind);
       end
     end
     V(:,:,k) = Vk;                    % load into full array


     F.V = zk;                         % place depths to interpolate for
                                       % level k in TriScatteredInterp
				       % object

     Zk = F(I.XR,I.YR);                % interpolate on receiver 
                                       % horizontal grid
     
     ind = find(isnan(Zk));
     if (~isempty(ind)),
       FN.V = zk;
       Zk(ind) = FN(I.XR(ind),I.YR(ind));
     end
     Z(:,:,k) = Zk;                    % load depths into full array
     
   end

%  Reshape interpolation data to facilitate vectorization

   VD = reshape(V   , [numel(I.XR) KmD]);
   ZD = reshape(Z   , [numel(I.XR) KmD]);
   ZR = reshape(I.ZR, [numel(I.XR) KmR]);
 
%  Duplicate shallowest and deepest levels of data to correspond
%  to bogus ZD = 0 and ZR = -I.hmaxR. This forces simple vertical 
%  extrapolation if the Receiver Grid has depths outside the Donor
%  range. This will generally be the case if the receiver has more
%  levels or a stretching that emphasizes boundary layers.

   VD = VD(:,[1 1:end end]);
   ZD = [I.Zbot*ones([length(ZD) 1]), ZD, I.Zsur*ones([length(ZD) 1])];

% Lastly, perform vertical interpolation.

   V = NaN([numel(I.XR) KmR]);

   for i=1:size(ZR,1)
     V(i,:) = interp1(ZD(i,:),VD(i,:),ZR(i,:),Vmethod);
   end
   V = reshape(V, size(I.ZR));
   
   Rmin = min(V(:));
   Rmax = max(V(:));
    
   Rind = find(repmat(I.Rmask,[1,1,KmR]) < 0.5);
   if (~isempty(Rind)),
     V(Rind) = 0;
   end

   ind = find(isnan(V));
   if (~isempty(ind)),
       Ncount = length(ind);
   end   

   disp(['      Donor Min = ', sprintf('%12.5e',Dmin), '  ',            ...
         '    Donor Max = ', sprintf('%12.5e',Dmax)]);
   disp(['   Receiver Min = ', sprintf('%12.5e',Rmin), '  ',            ...
           ' Receiver Max = ', sprintf('%12.5e',Rmax), '  ',            ...
           ' NaN count = ',  num2str(Ncount)]);        

end

return
