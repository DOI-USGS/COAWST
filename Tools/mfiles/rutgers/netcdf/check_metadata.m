function S = check_metadata(Sinp)

%
% CHECK_METADATA:  Checks metadata structure for a ROMS NetCDF file
%
% Sout = check_metadata(Sinp)
%
% This function checks ROMS metadata structure for consistency and
% fills unassigned fields. This structure will be used elsewhere to
% create NetCDF files in a compact way (see "nc_create.m") or to
% append new variables to existing NetCDF files (see "nc_append.m").
%
% This structure has similar schema as that returned by:
%
%    S = nc_inq(ncdife)
% or
%    S = ncinfo(ncfile)          native Matlab function
%
% On Input:
%
%    Sinp       ROMS metadata structure (struct array)
%
% On Ouput:
%
%    S          Updated ROMS metadata structure (struct array)
%

% svn $Id: check_metadata.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMVariables.txt                   Hernan G. Arango      %
%=========================================================================%

if (~isstruct(Sinp)),
  error('CHECK_METADATA: input argument ''Sinp'' is not a structure');
end
  
% Initialize output structure.

S = Sinp;

%--------------------------------------------------------------------------
% Check ROMS metadata structure.
%--------------------------------------------------------------------------

nvars = length(S.Variables);

for n=1:nvars,

% Fill variable dimension(s) length and variable size.
  
  nvdims = length(S.Variables(n).Dimensions);
  vsize  = [];
  
  for i=1:nvdims,
    dname   = char(S.Variables(n).Dimensions(i).Name);
    dindex  = strcmp({S.Dimensions.Name}, dname);

    RenameDim = false;
    if (any(dindex)),
      dsize = S.Dimensions(dindex).Length;
    else
      RenameDim = true;
    end

% If applicable, rename time dimension. In 'roms_metadata' all the state
% variables have 'ocean_time' as default.  However, it posible to have
% other names like 'time', 'clm_time', etc when processing climatology
% or other type of files that have ROMS state variables names.  Here the
% file dimensions are scanned for the substring 'time' and the structure
% is corrected with the appropriate available dimension.

    if (RenameDim),
      foundit = false;
      if (strcmp(dname, 'ocean_time'))
	    dindex = strfind({S.Dimensions.Name}, 'time');
        dindex = ~cellfun(@isempty, dindex);
        if (any(dindex)),
          dname = S.Dimensions(dindex).Name;
          dsize = S.Dimensions(dindex).Length;
          S.Variables(n).Dimensions(i).Name = dname;
          foundit = true;
        end
      end
      if (~foundit),
        error(['CHECK_METADATA: dimension "',dname,'" is not ',         ...
               'available for variable "', char(S.Variables(n).Name),'"']);
      end
    end     
    
    Dcel{i} = dname;                                % horizontal cell array

    if (S.Variables(n).Dimensions(i).Unlimited);
      S.Variables(n).Dimensions(i).Length = 0;
      vsize = [vsize 0];
    else      
      S.Variables(n).Dimensions(i).Length = dsize;
      vsize = [vsize dsize];
    end
  end
  if (nvdims > 0),
    S.Variables(n).Size = vsize;
  end

% If applicable, check the 'coordinates' attribute. If the 'coordinates'
% attribute string "Astr", saved also as cell array "Acel", is not equal
% to the variable dimensions cell array "Dcel", replace its value with
% the correct attribute string "Avalue".

   iatt = strcmp({S.Variables(n).Attributes.Name}, 'coordinates');
   if (any(iatt)),
     Astr = S.Variables(n).Attributes(iatt).Value;

     Acel = textscan(Astr, '%s');                   % vertical cell array

     if (~isequal(Acel{1}, Dcel')),                 % notice transpose
       Avalue = sprintf('%s ',Dcel{:});
       S.Variables(n).Attributes(iatt).Value = Avalue;
     end
   end
   clear iatt

% If applicable, check the 'time' attribute. If the 'time' attribute
% string "Astr" is not equal to time record dimensions string Dcel{end},
% replace its value with the correct attribute string.

   iatt = strcmp({S.Variables(n).Attributes.Name}, 'time');
   if (any(iatt)),
     Astr = S.Variables(n).Attributes(iatt).Value;

     if (~isequal(Astr, char(Dcel{end}))),
       S.Variables(n).Attributes(iatt).Value = char(Dcel{end});
     end
   end
   
end  
  
return