function check_records (ncnames, vname)

% CHECK_RECORDS:  Checks multi-file records for monotonicity
%
% check_records(ncnames)
%
% Checks multi-files for increasely monotinic time coordinate. In
% order for the multi-files to work while running the adjoint model
% the time coordinate has to be sctrictly monotonoc increasing for
% the backward interpolation to be bounded.
%
% On Input:
%
%    ncnames    Data NetCDF filenames (cell array)
%
%    vname      Record time variable to check (string)
%

% Get reference time.

I=nc_inq(char(ncnames(1)));

ind1=strcmp({I.Variables.Name}, vname);
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

%  Print upper and lower available time data information.

Nfiles=length(ncnames);

disp(' ');
for n=1:Nfiles,
  time=nc_read(char(ncnames(n)), vname);
  time=time.*scale;
  time_str=datestr(epoch+time);
  disp(['File : ', char(ncnames(n))]);
  disp(['       Record(1)     = ', time_str(1,:),                       ...
        '  (', num2str(time(1)),')']);
  disp(['       Record(2)     = ', time_str(2,:),                       ...
        '  (', num2str(time(2)),')']);
  disp(['       Record(end-1) = ', time_str(end-1,:),                   ...
        '  (', num2str(time(end-1)),')' ]);
  disp(['       Record(end  ) = ', time_str(end,:),                     ...
        '  (', num2str(time(end)), ')']);
  disp(' ');
end

return

