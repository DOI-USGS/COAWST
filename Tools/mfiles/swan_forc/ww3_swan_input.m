% Download, convert and interpolate WW3 data to TPAR input for SWAN
% Written 4/23/09 by Brandy Armstrong


% *************USER INPUT****************************

% working directory, location of ww3 grb files downloaded from 
% ftp://polar.ncep.noaa.gov/pub/history/waves/
working_drive='D:';
% ***WARNING*** TPAR files are saved in the working directory and are named
% generically (numbered), so make sure your working directory does not already have
% TPAR files with the same names
working_dir='Temporary\ww3\';
eval(['cd ',working_drive,'\',working_dir,';']);

% dates of data requested
yearww3='2005';%input year of data yyyy 
mmww3='12';%input month of data mm
% for recent data, last 7 days, input exact day, time period covered is 
% approximately 7.5 days
% if older than the last 7 days keep ddww3='00'; for archived grb files
ddww3='00';% dd

% Enter the ww3 area needed to cover the grid
% akw, enp, nah, nph, nww3, wna
ww3_area='wna'; %western north atlantic

% location/name of netcdf grid
modelgrid='D:\BNA\data\Grids\USeast_grd12.nc';

% determine spec points and interpolate spec for TPAR files
% ww3_specpoints assumes masking of 0 for land and NaN for water
if str2num(ddww3)==0; %assumes these are grb files
    ww3gb_2TPAR(modelgrid,yearww3,mmww3,working_dir,working_drive,ww3_area,ddww3)
else
    ww3nc_2TPAR(modelgrid,yearww3,mmww3,working_dir,working_drive,ww3_area,ddww3)
end


