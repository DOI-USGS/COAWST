function netcdf_load(ncfile);

% function netcdf_load('ncfile');
% ncloads the whole netcdf file
% jcwarner 31Jan2012
%
ncid = netcdf.open(ncfile,'NC_NOWRITE');

[ndims, nvars, ngatts, unlimdimid] = netcdf.inq(ncid);

for m=0:nvars-1
  [varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,m);
% eval([varname,' = netcdf.getVar(ncid,m);'])
  val = netcdf.getVar(ncid,m);
  assignin('base',varname,val);
end

netcdf.close(ncid)
