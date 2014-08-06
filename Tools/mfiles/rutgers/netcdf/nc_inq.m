function I = nc_inq(ncfile, Lprint)

%
% NC_INQ:  Inquire about the contents of a NetCDF file
%
% I = nc_inq(ncfile, Lprint)
%
% This gets and prints the contents of a NetCDF file.  It displays the
% dimensions variables.
%
% On Input:
%
%    ncfile      NetCDF file name or URL file name (string)
%    Lprint      Switch to print information (optional)
%
% On Ouput:
%
%    I           NetCDF information (struct array):
%
%                 I.Filename    NetCDF file name (string)
%
%                 I.Attributes  NetCDF global attributes (struct array):
%
%                                 I.Attributes(:).Name
%                                 I.Attributes(:).Value
%
%                 I.Dimensions  NetCDF file dimensions (struct array):
%
%                                 I.Dimensions(:).Name
%                                 I.Dimensions(:).Length
%                                 I.Dimensions(:).Unlimited
%
%                 I.Variables   NetCDF file variables (struct array):
%
%                                 I.Variables(:).Name
%
%                                 I.Variables(:).Dimensions(:).Name
%                                 I.Variables(:).Dimensions(:).Length
%                                 I.Variables(:).Dimensions(:).Unlimited
%
%                                 I.Variables(:).Size
%
%                                 I.Variables(:).Attributes(:).Name
%                                 I.Variables(:).Attributes(:).Value
%
%                                 I.Variables(:).Cgridtype.Name
%                                 I.Variables(:).Cgridtype.Value
%
%                                 I.Variables(:).Datatype
%
%                                 I.Variables(:).ncType
%

% svn $Id: nc_inq.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

if (nargin < 2),
  Lprint = false;
end

%--------------------------------------------------------------------------
% Build NetCDF infomation structures.
%--------------------------------------------------------------------------

A = nc_getatt(ncfile);            % NetCDF file global attributes
D = nc_dinfo (ncfile);            % NetCDF file dimensions
V = nc_vnames(ncfile);            % NetCDF file variables

if (isfield(D, 'dimid')),
  D = rmfield(D,'dimid');        % Remove dimension IDs
end

% Build output structure.

I.Filename   = ncfile;
I.Attributes = A;
I.Dimensions = D;
I.Variables  = V.Variables;

%--------------------------------------------------------------------------
% Report NetCDF information.
%--------------------------------------------------------------------------

% Report dimensions.

if (Lprint),
  disp(' ')
  disp('Available dimensions and values:');
  disp(' ')

  ndims = length(I.Dimensions);
  unlimited=0;

  for i=1:ndims,
    dimnam = I.Dimensions(i).Name;
    dimsiz = I.Dimensions(i).Length;
    unlimited = unlimited + I.Dimensions(i).Unlimited;

    if (i > 9),
      s = [' '  int2str(i) ') ' dimnam ' = ' int2str(dimsiz)];
    else
      s = ['  ' int2str(i) ') ' dimnam ' = ' int2str(dimsiz)];
    end,
    disp(s)
  end,

  if (unlimited == 0),
    disp(' ')
    disp('     None of the dimensions is unlimited.')
  else
    disp(' ')
    for i=1:ndims,
      if (I.Dimensions(i).Unlimited),
        dimnam = I.Dimensions(i).Name;
        s = ['     ' dimnam ' is unlimited in length.'];
        disp(s)
      end    
    end
  end

% Report available variables.

  nvars=length(I.Variables);

  disp(' ')
  disp('Available Variables:');
  disp(' ')

  for i=1:3:nvars

    stri = int2str(i);

    if (length(stri) == 1)
      stri=[ ' ' stri];
    end
    varnam = I.Variables(i).Name;
    s = [ '  ' stri ') ' varnam ];
    addit = 26 - length(s);
    for j=1:addit
      s = [ s ' '];
    end
   
    if (i < nvars)
      stri = int2str(i+2);
      if (length(stri) == 1)
        stri = [ ' ' stri];
      end
      varnam = I.Variables(i+1).Name;
      s = [ s '  ' stri ') ' varnam ];
      addit = 52-length(s);
      for j=1:addit
        s = [ s ' '];
      end
    end 
   
    if (i < nvars - 1)
      stri = int2str(i+3);
      if (length(stri) == 1)
        stri = [ ' ' stri];
      end
      varnam = I.Variables(i+2).Name;
      s = [ s '  ' stri ') ' varnam ];
    end 
    disp(s)
  end
end

return
