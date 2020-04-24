function S = nc_stats(ncfile, Vname, varargin)

%
% NC_STATS:  Computes the statistics of NetCDF variable
%
% S = nc_stats(ncfile, Vname, Recs, Scale)
%
% This function computes the statistics of requested NetCDF variable for
% the specified record(s).  If land/sea mask is available, it computes
% the statistic for the water points only.  I can be used to determine
% global colormap dynamical range.
%
% On Input:
%
%    ncfile     NetCDF file name or URL name (character string)
%
%    Vname      Field variable name (character string)
%
%    Recs       Time record(s) to process (integer vector, optional)
%                 If ommited or Inf it computes statistics over all
%                 available records
%
%    Scale      Field scale value to change units (scalar, optional)
%
% On Output:
%
%    S          Requested variable statistics (struct array):
%
%                 S.ncfile      NetCDF file name (string)
%                 S.varname     Variable (string)
%                 S.long_name   Variable description (string)
%                 S.units       Variable units
%                 S.Nrecs       Number of time records processed
%                 S.Npts        Number of values processed
%                 S.min         Overall Minimum value
%                 S.max         Overall Maximum value
%                 S.avg         Overall Average value
%
%

% svn $Id: nc_stats.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialze.

S.ncfile = ncfile;
S.varname = Vname;
S.long_name = blanks(1);
S.units = blanks(1);
S.Nrecs = NaN;
S.Npts = NaN;
S.min = NaN;
S.max = NaN;
S.avg = NaN;

% Inquire about the NetCDF file.

I = nc_inq(ncfile);

% Check about requested variable.

is2d = false;
is3d = false;

if (any(strcmp({I.Variables.Name}, Vname)))
  ivar = find(strcmp({I.Variables.Name}, Vname));
  nvdims = length(I.Variables(ivar).Dimensions);
  vdnames = {I.Variables(ivar).Dimensions.Name};

  if (any(strcmp(vdnames, 's_rho')) ||  any(strcmp(vdnames, 's_w')))
    is3d = true;
    Km = I.Variables(ivar).Dimensions(3).Length;
  else
    is2d = true;
  end

% Inquire about time records.  
  
  idims = contains({I.Variables(ivar).Dimensions.Name},'time');
  if (any(idims))
    Nrec = I.Variables(ivar).Dimensions(idims).Length;
  else
    Nrec = I.Variables(ivar).Dimensions(nvdims).Length;
  end

% Get long_name and unit attributes.

ind = strcmp({I.Variables(ivar).Attributes.Name}, 'long_name');
if (any(ind))
  S.long_name = I.Variables(ivar).Attributes(ind).Value;
end

ind = strcmp({I.Variables(ivar).Attributes.Name}, 'units');
if (any(ind))
  S.units = I.Variables(ivar).Attributes(ind).Value;
else
  S.units = 'unitless';
end

% Inquire about land/sea mask

  imsk = find(contains({I.Variables.Name}, 'mask'));
  if (~isempty(imsk))
    for n=1:length(imsk)
      mdnames = {I.Variables(imsk(n)).Dimensions.Name};
      if (strcmp(vdnames(1), mdnames(1)) &&                             ...
          strcmp(vdnames(2), mdnames(2)))
        got_mask=true;
        mask_var = I.Variables(imsk(n)).Name;
        break
      end
    end
  else
    got_mask=false;
  end
else

  disp(blanks(1));
  disp(['Variable:  ', Vname, ' not founded.']);
  disp(blanks(1));
  disp('Available Variables:');
  disp(blanks(1));

  nvars=length(I.Variables);

  for i=1:3:nvars

    s = [];
    
    stri = int2str(i);

    if (length(stri) == 1)
      stri=[ ' ' stri];
    end
    nvdims = length(I.Variables(i).Dimensions);
    if (nvdims > 2)
      ind = contains({I.Variables(i).Dimensions.Name},'time');
      if (any(ind))
        varnam = I.Variables(i).Name;
        s = [ '  ' stri ') ' varnam ];
        addit = 26 - length(s);
        for j=1:addit
          s = [ s ' '];
        end
      end
    end
    
    if (i < nvars)
      stri = int2str(i+2);
      if (length(stri) == 1)
        stri = [ ' ' stri];
      end
      nvdims = length(I.Variables(i+1).Dimensions);
      if (nvdims > 2)
        ind = contains({I.Variables(i+1).Dimensions.Name},'time');
        if (any(ind))
          varnam = I.Variables(i+1).Name;
          s = [ s '  ' stri ') ' varnam ];
          addit = 52-length(s);
          for j=1:addit
            s = [ s ' '];
          end
        end
      end 
    end
    
    if (i < nvars - 1)
      stri = int2str(i+3);
      if (length(stri) == 1)
        stri = [ ' ' stri];
      end
      nvdims = length(I.Variables(i+2).Dimensions);
      if (nvdims > 2)
        ind = contains({I.Variables(i+2).Dimensions.Name},'time');
        if (any(ind))
          varnam = I.Variables(i+2).Name;
          s = [ s '  ' stri ') ' varnam ];
        end
      end 
    end

    if (~isempty(s))
      disp(s)
    end
  end
  return
end

% Set optional arguments.

switch numel(varargin)
  case 0
    Recs  = [1:Nrec];
    Scale = 1.0d0; 
  case 1
    if (isinf(varargin{1}))
      Recs = [1:Nrec];
    else
      Recs = varargin{1};
    end
     Scale = 1.0d0;
  case 2
    if (isinf(varargin{1}))
      Recs = [1:Nrec];
    else
      Recs = varargin{1};
    end
    Scale = varargin{2};
end

S.Nrecs = length(Recs);

%--------------------------------------------------------------------------
% Compute statistics.
%--------------------------------------------------------------------------

% Set the land/sea mask.

if (got_mask)
  mask = nc_read(ncfile, mask_var);
  if (is3d)
    mask = repmat(mask, [1 1 Km]);
  end
  ind = find(mask > 0);  
end

% Process requested variable.

Fmin = NaN(1,length(Recs));
Fmax = NaN(1,length(Recs));
Favg = NaN(1,length(Recs));
Fstd = NaN(1,length(Recs));

Npts = 0;

for rec = Recs
  
  F = nc_read(ncfile, Vname, rec);
  F = F(:) .* Scale;
  
  if (got_mask)
    Npts = Npts + length(ind);
    
    Fmin(rec) = min(F(ind));
    Fmax(rec) = max(F(ind));
    Favg(rec) = sum(F(ind));
  else
    Npts = Npts + length(F);

    Fmin(rec) = min(F(ind));
    Fmax(rec) = max(F(ind));
    Favg(rec) = sum(F(ind));
  end

end
 
%  Fill structure.

S.Npts = Npts;
S.min = min(Fmin);
S.max = max(Fmax);
S.avg = sum(Favg)/Npts;

disp(blanks(1));
disp(['   NetCDF filename: ', S.ncfile]);
disp(['          Variable: ', S.varname]);
disp(['       Description: ', S.long_name]);
disp(['             Units: ', S.units]);
disp([' Number of Records: ', num2str(S.Nrecs)]);
disp(['  Number of Values: ', num2str(S.Npts)]);
disp(['   Overall Minimum: ', num2str(S.min)]);
disp(['   Overall Maximum: ', num2str(S.max)]);
disp(['   Overall Average: ', num2str(S.avg)]);
disp(blanks(1));

return
