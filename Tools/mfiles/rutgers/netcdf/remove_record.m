function remove_record (InpName, OutName, Rec)

% REMOVE_RECORD:  Removes specified record in NetCDF file
%
% remove_record (InpName, OutName, Rec)
%
% Creates a new netcdf with requested time record removed.  It is
% used to remove repeated days in multi-files.
%
% On Input:
%
%    InpName    Input NetCDF filename to read (string)
%
%    OutName    Input NetCDF filename to create (string)
%
%    Rec        Data time record to remove (string)
%
%                 'first'      removes first record of data
%                 'last'       removes last  record of data
%

%--------------------------------------------------------------------------
% Get information about file.
%--------------------------------------------------------------------------

I = nc_inq(InpName);

%--------------------------------------------------------------------------
% Check unlimited or time record dimension.
%--------------------------------------------------------------------------

got_unlimited = false;

RecDim=[I.Dimensions.Unlimited];
if (any(RecDim))
  got_unlimited = true;
  Rname = I.Dimensions(RecDim).Name;
  Nrecs = I.Dimensions(RecDim).Length;
end

if (~got_unlimited)
  RecDim = contains({I.Dimensions.Name}, 'time');
  Rname  = I.Dimensions(RecDim).Name;
  Nrecs  = I.Dimensions(RecDim).Length;
end

%--------------------------------------------------------------------------
% Set new global attribute.
%--------------------------------------------------------------------------

% Get reference time.

ind1=strcmp({I.Variables.Name}, Rname);
if (any(ind1))
  ind2=strcmp({I.Variables(ind1).Attributes.Name}, 'units');
  if (any(ind2))
    Tunits=I.Variables(ind1).Attributes(ind2).Value;
  end
end

if (~isempty(Tunits))
  if (contains(Tunits, 'seconds since'))
    epoch_str=strsplit(Tunits, 'seconds since ');
    epoch=datenum(epoch_str(2));
    scale=1/86400;
  end
  if (contains(Tunits, 'days since'))
    epoch_str=strsplit(Tunits, 'days since ');
    epoch=datenum(epoch_str(2));
    scale=1;
  end
else
  epoch=0;
end

% Set 'modified' attribute.

if (ischar(Rec))
  switch lower(Rec)
   case 'last' 
     time=nc_read(InpName, Rname, Nrecs);
     time=time*scale;
     time_str=datestr(epoch+time);
   case 'first'
     time=nc_read(InpName, Rname, 1);
     time=time*scale;
     time_str=datestr(epoch+time);
  end
end

Avalue=['Removed record for ', time_str, ' using ', which(mfilename),   ...
	' on ', date_stamp];

%--------------------------------------------------------------------------
% Create new NetCDF file.
%--------------------------------------------------------------------------

% Modify information structure.

I.Dimensions(strcmp({I.Dimensions.Name}, Rname)).Length = Nrecs-1;

% Create NetCDF file.

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

nc_create(OutName, mode, I);

% Add 'modified' global attribute.

status = nc_attadd(OutName, 'modified', Avalue);

%--------------------------------------------------------------------------
% Write out NetCDF Variables.
%--------------------------------------------------------------------------

% Set records to process.

if (ischar(Rec))
  switch lower(Rec)
   case 'last'
     recs = 1:Nrecs-1;
   case 'first'
     recs = 2:Nrecs;
  end
end

% Write out variables.

nvars = length(I.Variables);

for n=1:nvars
  vname = char(I.Variables(n).Name);
  vsize = length(I.Variables(n).Size);

  if (isempty(I.Variables(n).Dimensions))

    Vinp = nc_read(InpName, vname);                 % scalar data
    if (~isempty(Vinp))
      status = ncwrite(OutName, vname, Vinp);
    end

  else

    if (any([I.Variables(n).Dimensions.Unlimited]) ||                   ...
        any(contains({I.Variables(n).Dimensions.Name}, 'time')))

      for m = 1:Nrecs-1                             % record data
        InpRec = recs(m);
        Vinp = nc_read(InpName, vname, InpRec);
        status = nc_write(OutName, vname, Vinp, m);
      end

    else

      Vinp = nc_read(InpName, vname);               % nonrecord data
      status = nc_write(OutName, vname, Vinp);

    end
  
  end

end

%--------------------------------------------------------------------------
% Check information about time records.
%--------------------------------------------------------------------------

check_records ({InpName, OutName}, Rname);

return
