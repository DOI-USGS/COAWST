function nc_test(ncfile, Interface, Lplot)

%
% NC_TEST:  Test reading/writing NetCDF data with various iterfaces
%
% nc_test(ncfile, interface, Lplot)
%
% This function creates a NetCDF using data from the peaks(40) function.
% Several datatype variables are created original double precision data
% and packed into byte (int8), short (int16), integer (int32), and single
% precission.  It is used to test the reading/writing NetCDF interface
% to Matlab. It computes and reports the processing and truncation RMSE.
%
% On Input:
%
%    ncfile      NetCDF file name to create (string)
%
%    Interface   NetCDF interface for Matlab (string):
%
%                  'native'     Native Matlab interface
%                  'mexnc'      MEXNC interface
%                  'roms'       ROMS 'nc_write' and 'nc_read' interface      
%                  'snctools'   SNCTOOLS inteface
%
%    Lplot       Switch to plot the data (default false, Optional)
%

% svn $Id: nc_test.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Initialize.
  
if (nargin < 3),
  Lplot = false;
end

original = false;        % switch for SNCTOOLS: original 'nc_varput' does
                         % not work... This is a nasty bug!

%  Create NetCDF file.

mode = netcdf.getConstant('CLOBBER');
ncid = netcdf.create(ncfile,mode);

%  Clear and set persistent switch to process data in column-major order.

if (strncmpi(Interface, 'SNCTOOLS', 3)),
  if (ispref('SNCTOOLS','PRESERVE_FVD')),
    saved_preserve_fvd = getpref('SNCTOOLS','PRESERVE_FVD');
  else
    saved_preserve_fvd = false;             % default value in SNCTOOLS
  end
  setpref('SNCTOOLS','PRESERVE_FVD',true);
end

%  Create data set, double precision.

A = peaks(40);

[Im,Jm] = size(A);
Amin = min(min(A));
Amax = max(max(A));
fill_double = 1.0d+37;

%  Single precision data: fillvalue=1.0E+37

A_single = single(A);
fill_single = single(fill_double);
A_singleMin = min(A_single(:));
A_singleMax = max(A_single(:));

%  Integer type data, int32:  fillvalue=999

scale_int = double(0.0001);
fill_int = int32(9990000);
A_int = int32(A./scale_int);
A_intMin = min(A_int(:));
A_intMax = max(A_int(:));

%  Short type data, int16: -32768:32767; fillvalue=-32768

max_short = 65534;
offset_short = double(0.5*(Amin+Amax));
scale_short = double((Amax-Amin)/max_short);
fill_short = int16(-32768);
A_short = int16((A-offset_short)./scale_short);
A_shortMin = min(A_short(:));
A_shortMax = max(A_short(:));

%  Byte type data, int8: -128:127; fillvalue=-128

max_byte = 254;
offset_byte = double(0.5*(Amin+Amax));
scale_byte = double((Amax-Amin)/max_byte);
fill_byte = int8(-128);
A_byte = int8((A-offset_byte)./scale_byte);
A_byteMin = min(A_byte(:));
A_byteMax = max(A_byte(:));

disp(' ');
disp(['   double range: ', sprintf('%15.8f',Amin), '  ',        ...
                           sprintf('%15.8f',Amax)]);
disp(['   single range: ', sprintf('%15.5f',A_singleMin), '  ', ...
                           sprintf('%15.5f',A_singleMax)]);
disp(['  integer range: ', sprintf('%15i',A_intMin), '  ',      ...
                           sprintf('%15i',A_intMax)]);
disp(['    short range: ', sprintf('%15i',A_shortMin),'  ',     ...
                           sprintf('%15i',A_shortMax)]);
disp(['     byte range: ', sprintf('%15i',A_byteMin), '  ',     ...
                           sprintf('%15i',A_byteMax)]);
disp(' ');

%  Set fill values.

B = A;
B(5:7,4:12) = NaN;
ind = isnan(B);

A(ind)        = fill_double;
A_single(ind) = fill_single;
A_int(ind)    = fill_int;
A_short(ind)  = fill_short;
A_byte(ind)   = fill_byte;

%  Create several variables of different type containing truncated
%  original double precision data.

did.Im=netcdf.defDim(ncid,'x',Im);
did.Jm=netcdf.defDim(ncid,'y',Jm);

varid=netcdf.defVar(ncid,'A_byte',netcdf.getConstant('nc_byte'),        ...
                    [did.Im did.Jm]);
netcdf.putAtt(ncid,varid,'long_name'   ,'a byte variable, int8');
netcdf.putAtt(ncid,varid,'scale_factor',scale_byte);
netcdf.putAtt(ncid,varid,'add_offset'  ,offset_byte);
netcdf.putAtt(ncid,varid,'_FillValue'  ,fill_byte);

varid=netcdf.defVar(ncid,'A_short',netcdf.getConstant('nc_short'),      ...
                    [did.Im did.Jm]);
netcdf.putAtt(ncid,varid,'long_name'   ,'a short variable, int16');
netcdf.putAtt(ncid,varid,'scale_factor',scale_short);
netcdf.putAtt(ncid,varid,'add_offset'  ,offset_short);
netcdf.putAtt(ncid,varid,'_FillValue'  ,fill_short);

varid=netcdf.defVar(ncid,'A_int',netcdf.getConstant('nc_int'),          ...
                    [did.Im did.Jm]);
netcdf.putAtt(ncid,varid,'long_name'   ,'an integer variable, int32');
netcdf.putAtt(ncid,varid,'scale_factor',scale_int);
netcdf.putAtt(ncid,varid,'_FillValue'  ,fill_int);

varid=netcdf.defVar(ncid,'A_single',netcdf.getConstant('nc_float'),     ...
                    [did.Im did.Jm]);
netcdf.putAtt(ncid,varid,'long_name'   ,'a single precision variable');
netcdf.putAtt(ncid,varid,'_FillValue'  ,fill_single);

varid=netcdf.defVar(ncid,'A_double',netcdf.getConstant('nc_double'),    ...
                    [did.Im did.Jm]);
netcdf.putAtt(ncid,varid,'long_name'   ,'a double precision variable');
netcdf.putAtt(ncid,varid,'_FillValue'  ,fill_double);

%  Leave definition mode and close NetCDF file.

netcdf.endDef(ncid);
netcdf.close(ncid);

%--------------------------------------------------------------------------
%  Write out and read in data.
%--------------------------------------------------------------------------

if (strncmpi(Interface,'NATIVE',3)),
  ncwrite(ncfile, 'A_byte'  , A_byte  );
  ncwrite(ncfile, 'A_short' , A_short );
  ncwrite(ncfile, 'A_int'   , A_int   );
  ncwrite(ncfile, 'A_single', A_single);
  ncwrite(ncfile, 'A_double', A       );

  R.byte   = ncread(ncfile, 'A_byte'  );
  R.short  = ncread(ncfile, 'A_short' );
  R.int    = ncread(ncfile, 'A_int'   );
  R.single = ncread(ncfile, 'A_single');
  R.double = ncread(ncfile, 'A_double');
  
elseif (strncmpi(Interface,'MEXNC',3)),
  ncid = mexnc('ncopen', ncfile, 'nc_write');

  status = mexnc('put_var_schar' , ncid, 0, A_byte  );
  status = mexnc('put_var_short' , ncid, 1, A_short );
  status = mexnc('put_var_int'   , ncid, 2, A_int   );
  status = mexnc('put_var_float' , ncid, 3, A_single);
  status = mexnc('put_var_double', ncid, 4, A       );

  status = mexnc('ncclose',ncid);

  ncid = mexnc('ncopen', ncfile, 'nc_nowrite');

  R.byte   = mexnc('get_var_schar' , ncid, 0);
  my_ind   = find(abs(R.byte-fill_byte) < 4*eps(double(fill_byte)));
  R.byte   = double(R.byte).*scale_byte + offset_byte;
             R.byte(my_ind) = NaN;

  R.short  = mexnc('get_var_short' , ncid, 1);
  my_ind   = find(abs(R.short-fill_short) < 4*eps(double(fill_short)));
  R.short  = double(R.short).*scale_short + offset_short;
	     R.short(my_ind) = NaN;

  R.int    = mexnc('get_var_int'   , ncid, 2);
  my_ind   = find(abs(R.int-fill_int) < 4*eps(double(fill_int)));
  R.int    = double(R.int).*scale_int;
	     R.int(my_ind) = NaN;

  R.single = mexnc('get_var_float' , ncid, 3);
  my_ind   = find(abs(R.single-fill_single) < 4*eps(double(fill_single)));
	     R.single(my_ind) = NaN;

  R.double = mexnc('get_var_double', ncid, 4);
  my_ind   = find(abs(R.double-fill_double) < 4*eps(fill_double));
	     R.double(my_ind) = NaN;

  status = mexnc('ncclose',ncid);

  disp(' ');

elseif (strncmpi(Interface,'ROMS',4)),
  status = nc_write(ncfile, 'A_byte'  , A_byte  );
  status = nc_write(ncfile, 'A_short' , A_short );
  status = nc_write(ncfile, 'A_int'   , A_int   );
  status = nc_write(ncfile, 'A_single', A_single);
  status = nc_write(ncfile, 'A_double', A       );

  ReplaceValue = NaN;
  PreserveType = false;
  
  R.byte   = nc_read(ncfile, 'A_byte'  , [], ReplaceValue, PreserveType);
  R.short  = nc_read(ncfile, 'A_short' , [], ReplaceValue, PreserveType);
  R.int    = nc_read(ncfile, 'A_int'   , [], ReplaceValue, PreserveType);
  R.single = nc_read(ncfile, 'A_single', [], ReplaceValue, PreserveType);
  R.double = nc_read(ncfile, 'A_double', [], ReplaceValue, PreserveType);

  disp(' ');

elseif (strncmpi(Interface,'SNCTOOLS',3)),
  if (original),
    nc_varput(ncfile, 'A_byte'  , A_byte  );    % SNCTOOLS fails to do this
    nc_varput(ncfile, 'A_short' , A_short );    % correctly... The logic is
    nc_varput(ncfile, 'A_int'   , A_int   );    % unwisely biased towards
    nc_varput(ncfile, 'A_single', A_single);    % double precision. Even
    nc_varput(ncfile, 'A_double', A       );    % the reverse linear
  else                                          % transformation is wrong!
    nc_write(ncfile, 'A_byte'  , A_byte  );
    nc_write(ncfile, 'A_short' , A_short );
    nc_write(ncfile, 'A_int'   , A_int   );
    nc_write(ncfile, 'A_single', A_single);
    nc_write(ncfile, 'A_double', A       );
  end

  R.byte   = nc_varget(ncfile, 'A_byte'  );
  R.short  = nc_varget(ncfile, 'A_short' );
  R.int    = nc_varget(ncfile, 'A_int'   );
  R.single = nc_varget(ncfile, 'A_single');
  R.double = nc_varget(ncfile, 'A_double');
 
  disp(' ');
end

% Set persistent switch back to user or SNCTOOLS default value.

if (strncmpi(Interface, 'SNCTOOLS', 3)),
 setpref('SNCTOOLS','PRESERVE_FVD', saved_preserve_fvd);
end

%  Report.

ncdisp(ncfile);

info = ncinfo(ncfile);

%  Plot results.

if (Lplot),
  figure(1)
  B=A;  B(ind)=NaN;
  D=R.double;
  pcolor(B'); shading interp; colorbar;
  title('double precision data');
  B(ind)=[];
  D(ind)=[]; r_rms=sqrt(mean(B(:)-D(:)).^2);
  xlabel(['Processing RMSE = ', num2str(r_rms)]);
  
  figure(2)
  C=double(A_single);  C(ind)=NaN;
  D=R.single;
  pcolor(C'); shading interp; colorbar;
  title('single precision data');
  C(ind)=[]; p_rms=sqrt(mean(C(:)-B(:)).^2);
  D(ind)=[]; r_rms=sqrt(mean(C(:)-D(:)).^2);
  xlabel(['Truncation RMSE = ', num2str(p_rms), ',   ',                 ...
          'Processing RMSE = ', num2str(r_rms)]);
  
  figure(3)
  C=double(A_int).*scale_int; C(ind)=NaN;
  D=R.int;
  pcolor(C'); shading interp; colorbar;
  title('unpacked integer data, int32');
  C(ind)=[]; p_rms=sqrt(mean(C(:)-B(:)).^2);
  D(ind)=[]; r_rms=sqrt(mean(C(:)-D(:)).^2);
  xlabel(['Truncation RMSE = ', num2str(p_rms), ',   ',                 ...
          'Processing RMSE = ', num2str(r_rms)]);

  figure(4)
  C=double(A_short).*scale_short+offset_short; C(ind)=NaN;
  D=R.short;
  pcolor(C'); shading interp; colorbar;
  title('unpacked short data, int16');
  C(ind)=[]; p_rms=sqrt(mean(C(:)-B(:)).^2);
  D(ind)=[]; r_rms=sqrt(mean(C(:)-D(:)).^2);
  xlabel(['Truncation RMSE = ', num2str(p_rms), ',   ',                 ...
          'Processing RMSE = ', num2str(r_rms)]);

  figure(5)
  C=double(A_byte).*scale_byte+offset_byte; C(ind)=NaN;
  D=R.byte;
  pcolor(C'); shading interp; colorbar;
  title('unpacked byte data, int8');
  C(ind)=[]; p_rms=sqrt(mean(C(:)-B(:)).^2);
  D(ind)=[]; r_rms=sqrt(mean(C(:)-D(:)).^2);
  xlabel(['Truncation RMSE = ', num2str(p_rms), ',   ',                 ...
          'Processing RMSE = ', num2str(r_rms)]);
end
