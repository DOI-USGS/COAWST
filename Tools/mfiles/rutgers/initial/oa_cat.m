function oa_cat(outfile,annfile,monfile)

% OA_CAT:  Concatenates monthly and annual OA data
%
% oa_cat(outfile,annfile,monfile)
%
% This function reads Annual and Monthly objective analysis (OA) files
% of the Levitus climatology and appends bottom Annual levels to Monthly
% OA fields. The Monthly Levitus Climatology is available from the 
% surface to 1000m (1994 dataset) or 1500m (1998 dataset). The annual
% fields are used for the missing bottom levels.
%
% The strategy here is to copy the annual file as 'outfile' and replace
% the upper water column levels with the monthly data.
%
% On Input:
%
%    outfile     New monthly climatology NetCDF file name (string)
%
%    annfile     Annual  climatology NetCDF file name (string)
%
%    monfile     Monthly climatology NetCDF file name (string)
%
% Notice:
%
% The OA ocean depths (zout) are assumed to be negative.  Therefore, the
% order of levels are from bottom to surface.  For example, the depths of
% Levitus standard levels are:
%
% zout = [-5500, -5000, -4500, -4000, -3500, -3000, -2500, -2000,       ...
%         -1750, -1500, -1400, -1300, -1200, -1100, -1000,  -900,       ...  
%          -800,  -700,  -600,  -500,  -400,  -300,  -250,  -200,       ...
%          -150,  -125,  -100,   -75,   -50,   -30,   -20,   -10,       ...
%             0]  
%
% Make sure that the zout variable in the OA NetCDF follows similar order
% from bottom to surface.  It is always a good idea to OA a particular
% application to the relevant standard levels. Then, any ROMS vertical
% coordinate stretching can be interpolated from the standard levels.
%

% svn $Id: oa_cat.m 996 2020-01-10 04:28:56Z arango $
%=========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Inquire about the annual and Monthly files.
  
A = nc_inq(annfile);
M = nc_inq(monfile);

% Determine upper water column levels to replace.

TindA = strcmp({A.Variables.Name},'temp');
TindM = strcmp({M.Variables.Name},'temp');

if (any(TindA)),
  Asize  = [A.Variables(TindA).Dimensions.Length];
  KmA    = Asize(3);
  Adepth = nc_read(annfile,'zout');
else
  error(['OA_CAT: unable to find ''temp'' in:', annfile]);
end

if (any(TindM)),
  Msize = [M.Variables(TindM).Dimensions.Length];
  KmM   = Msize(3);
  Mdepth = nc_read(monfile,'zout');
else
  error(['OA_CAT: unable to find ''temp'' in:', monfile]);
end

for n=1:2
  if (Asize(n) ~= Msize(n)),
    error(['OA_CAT: Dimensions mismatch, ',                             ...
           A.Variables(IndA).Dimensions(n).Name, ' = ',                 ...
           numstr(Asize(n)), '  ', num2str(Msize(n))]);
  end
end

if (abs(Adepth(1)) < abs(Adepth(end))),
  error(['OA_CAT: not supported order OA depth levels in: ', annfile]);
end

if (abs(Mdepth(1)) < abs(Mdepth(end))),
  error(['OA_CAT: not supported order OA depth levels in: ', monfile]);
end
  
Nrecs = Asize(end);

%--------------------------------------------------------------------------
% Copy annual file to new monthly climatology file.
%--------------------------------------------------------------------------

[success,~,~] = copyfile(annfile, outfile);

if (~success),
  error(['OA_CAT: unable to create file: ', outfile]);
end

%--------------------------------------------------------------------------
% Replace upper levels for all 3D fields.
%--------------------------------------------------------------------------

Nvars = length({A.Variables.Name});

for n=1:Nvars,
  if (length(A.Variables(n).Size) > 3),
    vname = char(A.Variables(n).Name);

    for Rec=1:Nrecs
      Fann = nc_read(annfile, vname, Rec);
      Fmon = nc_read(annfile, vname, Rec);

      for level=0:KmM-1
	Lann = KmA-level;
	Lmon = KmM-level;
        Fann(:,:,Lann) = Fmon(:,:,Lmon);
      end
      
      status = nc_write(outfile, vname, Fann, Rec);
      if (status ~= 0),
        error(['OA_CAT: error while processing: ', vname]);
      end
    end
  
  end
end

% Write appropriate time.

time = nc_read(monfile,'time');
status = nc_write(outfile,'time',time);
if (status ~= 0),
  error('OA_CAT: error while writing time');
end

%--------------------------------------------------------------------------
% Update global attributes
%--------------------------------------------------------------------------
      
Natts = length({M.Attributes.Name});

for n=1:Natts,
  Aname  = char(M.Attributes(n).Name);
  Avalue = M.Attributes(n).Value;

  if (strcmp(Aname, 'out_file')),
    Avalue = outfile;
  end
  status = nc_attadd(outfile, Aname, Avalue);
    
  if (status ~= 0),
    error(['OA_CAT: error while writing global attribute: ', Aname]);
  end
end

return

