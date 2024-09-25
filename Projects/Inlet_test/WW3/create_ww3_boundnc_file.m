% create_ww3_boundnc_file
%
% create a netcdf file with boundary 2D spec for WW3
%
% jcw 20Sep2024
%

fn='ww3_bound_spec.nc'

%Dimensions:
time      = 2     %# of times
station   = 1     %# of stations
string40  = 40    %leave this alone
frequency = 30    %# of freqs
direction = 30    %# of dirs

%% create bndry file
nc_bndry=netcdf.create(fn,bitor(0,4096));   %JBZ update for NC4 files, equivalent to 'clobber' + 'NETCDF4'
if isempty(nc_bndry), return, end

%% Global attributes:
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_bndry,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by create_ww3_boundnc on ' datestr(now)]);
%% Dimensions:

disp(' ## Defining Dimensions...')
 
timedimID = netcdf.defDim(nc_bndry,'time',time);
stationdimID = netcdf.defDim(nc_bndry,'station',station);
string40dimID = netcdf.defDim(nc_bndry,'string40',string40);
frequencydimID = netcdf.defDim(nc_bndry,'frequency',frequency);
directiondimID = netcdf.defDim(nc_bndry,'direction',direction);
 
%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')
 
timeID = netcdf.defVar(nc_bndry,'time','double',[timedimID]);
netcdf.putAtt(nc_bndry,timeID,'long_name','julian day (UT)');
netcdf.putAtt(nc_bndry,timeID,'standard_name','time');
netcdf.putAtt(nc_bndry,timeID,'units','days since 1990-01-01 00:00:00');
netcdf.putAtt(nc_bndry,timeID,'conventions','Relative julian days with decimal part (as parts of the day)');
netcdf.putAtt(nc_bndry,timeID,'axis','T');
netcdf.putAtt(nc_bndry,timeID,'calendar','standard');

stationID = netcdf.defVar(nc_bndry,'station','nc_int',[stationdimID]);
netcdf.putAtt(nc_bndry,stationID,'long_name','station id');
%netcdf.putAtt(nc_bndry,stationID,'_FillValue',-2147483647);
netcdf.putAtt(nc_bndry,stationID,'axis','X');

%string40ID = netcdf.defVar(nc_bndry,'string40','int32',string40dimID);
%netcdf.putAtt(nc_bndry,string40ID,'long_name','station_name number of characters');
%netcdf.putAtt(nc_bndry,string40ID,'_FillValue',-2147483647);
%netcdf.putAtt(nc_bndry,string40ID,'axis','W');

station_nameID = netcdf.defVar(nc_bndry,'station_name','char',[stationdimID string40dimID]);
netcdf.putAtt(nc_bndry,station_nameID,'long_name','station name');
netcdf.putAtt(nc_bndry,station_nameID,'content','XW');
netcdf.putAtt(nc_bndry,station_nameID,'associates','station string40');

longitudeID = netcdf.defVar(nc_bndry,'longitude','nc_double',[stationdimID timedimID]);
netcdf.putAtt(nc_bndry,longitudeID,'long_name','longitude');
netcdf.putAtt(nc_bndry,longitudeID,'standard_name','longitude');
netcdf.putAtt(nc_bndry,longitudeID,'globwave_name','longitude');
netcdf.putAtt(nc_bndry,longitudeID,'units','degree_east');
netcdf.putAtt(nc_bndry,longitudeID,'scale_factor',1);
netcdf.putAtt(nc_bndry,longitudeID,'add_offset',0);
netcdf.putAtt(nc_bndry,longitudeID,'valid_min',-180);
netcdf.putAtt(nc_bndry,longitudeID,'valid_max',360);
netcdf.putAtt(nc_bndry,longitudeID,'_FillValue',9.969209968386869e+36);
netcdf.putAtt(nc_bndry,longitudeID,'content','TX');
netcdf.putAtt(nc_bndry,longitudeID,'associates','time station');

latitudeID = netcdf.defVar(nc_bndry,'latitude','nc_double',[stationdimID timedimID]);
netcdf.putAtt(nc_bndry,latitudeID,'long_name','latitude');
netcdf.putAtt(nc_bndry,latitudeID,'standard_name','latitude');
netcdf.putAtt(nc_bndry,latitudeID,'globwave_name','latitude');
netcdf.putAtt(nc_bndry,latitudeID,'units','degree_east');
netcdf.putAtt(nc_bndry,latitudeID,'scale_factor',1);
netcdf.putAtt(nc_bndry,latitudeID,'add_offset',0);
netcdf.putAtt(nc_bndry,latitudeID,'valid_min',-180);
netcdf.putAtt(nc_bndry,latitudeID,'valid_max',360);
netcdf.putAtt(nc_bndry,latitudeID,'_FillValue',9.969209968386869e+36);
netcdf.putAtt(nc_bndry,latitudeID,'content','TX');
netcdf.putAtt(nc_bndry,latitudeID,'associates','time station');

xID = netcdf.defVar(nc_bndry,'x','nc_double',[stationdimID timedimID]);
netcdf.putAtt(nc_bndry,xID,'long_name','x');
netcdf.putAtt(nc_bndry,xID,'standard_name','x');
netcdf.putAtt(nc_bndry,xID,'globwave_name','x');
netcdf.putAtt(nc_bndry,xID,'units','distance m');
netcdf.putAtt(nc_bndry,xID,'scale_factor',1);
netcdf.putAtt(nc_bndry,xID,'add_offset',0);
netcdf.putAtt(nc_bndry,xID,'valid_min',-180);
netcdf.putAtt(nc_bndry,xID,'valid_max',360);
netcdf.putAtt(nc_bndry,xID,'_FillValue',9.969209968386869e+36);
netcdf.putAtt(nc_bndry,xID,'content','TX');
netcdf.putAtt(nc_bndry,xID,'associates','time station');

yID = netcdf.defVar(nc_bndry,'y','nc_double',[stationdimID timedimID]);
netcdf.putAtt(nc_bndry,yID,'long_name','y');
netcdf.putAtt(nc_bndry,yID,'standard_name','y');
netcdf.putAtt(nc_bndry,yID,'globwave_name','y');
netcdf.putAtt(nc_bndry,yID,'units','distance m');
netcdf.putAtt(nc_bndry,yID,'scale_factor',1);
netcdf.putAtt(nc_bndry,yID,'add_offset',0);
netcdf.putAtt(nc_bndry,yID,'valid_min',-180);
netcdf.putAtt(nc_bndry,yID,'valid_max',360);
netcdf.putAtt(nc_bndry,yID,'_FillValue',9.969209968386869e+36);
netcdf.putAtt(nc_bndry,yID,'content','TX');
netcdf.putAtt(nc_bndry,yID,'associates','time station');

frequencyID = netcdf.defVar(nc_bndry,'frequency','nc_double',[frequencydimID]);
netcdf.putAtt(nc_bndry,frequencyID,'long_name','frequency of center band');
netcdf.putAtt(nc_bndry,frequencyID,'standard_name','sea_surface_wave_frequency');
netcdf.putAtt(nc_bndry,frequencyID,'globwave_name','frequency');
netcdf.putAtt(nc_bndry,frequencyID,'units','s-1');
netcdf.putAtt(nc_bndry,frequencyID,'scale_factor',1);
netcdf.putAtt(nc_bndry,frequencyID,'add_offset',0);
netcdf.putAtt(nc_bndry,frequencyID,'valid_min',0);
netcdf.putAtt(nc_bndry,frequencyID,'valid_max',10);
netcdf.putAtt(nc_bndry,frequencyID,'_FillValue',9.969209968386869e+36);
netcdf.putAtt(nc_bndry,frequencyID,'axis','Y');

frequency1ID = netcdf.defVar(nc_bndry,'frequency1','nc_double',[frequencydimID]);
netcdf.putAtt(nc_bndry,frequency1ID,'long_name','frequency of lower band');
netcdf.putAtt(nc_bndry,frequency1ID,'standard_name','frequency_of_lower_band');
netcdf.putAtt(nc_bndry,frequency1ID,'globwave_name','frequency_lower_band');
netcdf.putAtt(nc_bndry,frequency1ID,'units','s-1');
netcdf.putAtt(nc_bndry,frequency1ID,'scale_factor',1);
netcdf.putAtt(nc_bndry,frequency1ID,'add_offset',0);
netcdf.putAtt(nc_bndry,frequency1ID,'valid_min',0);
netcdf.putAtt(nc_bndry,frequency1ID,'valid_max',10);
netcdf.putAtt(nc_bndry,frequency1ID,'_FillValue',9.969209968386869e+36);
netcdf.putAtt(nc_bndry,frequency1ID,'axis','Y');
netcdf.putAtt(nc_bndry,frequency1ID,'associates','frequency');

frequency2ID = netcdf.defVar(nc_bndry,'frequency2','nc_double',[frequencydimID]);
netcdf.putAtt(nc_bndry,frequency2ID,'long_name','frequency of upper band');
netcdf.putAtt(nc_bndry,frequency2ID,'standard_name','frequency_of_upper_band');
netcdf.putAtt(nc_bndry,frequency2ID,'globwave_name','frequency_upper_band');
netcdf.putAtt(nc_bndry,frequency2ID,'units','s-1');
netcdf.putAtt(nc_bndry,frequency2ID,'scale_factor',1);
netcdf.putAtt(nc_bndry,frequency2ID,'add_offset',0);
netcdf.putAtt(nc_bndry,frequency2ID,'valid_min',0);
netcdf.putAtt(nc_bndry,frequency2ID,'valid_max',10);
netcdf.putAtt(nc_bndry,frequency2ID,'_FillValue',9.969209968386869e+36);
netcdf.putAtt(nc_bndry,frequency2ID,'axis','Y');
netcdf.putAtt(nc_bndry,frequency2ID,'associates','frequency');

directionID = netcdf.defVar(nc_bndry,'direction','nc_double',[directiondimID]);
netcdf.putAtt(nc_bndry,directionID,'long_name','sea surface wave to direction');
netcdf.putAtt(nc_bndry,directionID,'standard_name','sea_surface_wave_to_direction');
netcdf.putAtt(nc_bndry,directionID,'globwave_name','direction');
netcdf.putAtt(nc_bndry,directionID,'units','degree');
netcdf.putAtt(nc_bndry,directionID,'scale_factor',1);
netcdf.putAtt(nc_bndry,directionID,'add_offset',0);
netcdf.putAtt(nc_bndry,directionID,'valid_min',0);
netcdf.putAtt(nc_bndry,directionID,'valid_max',360);
netcdf.putAtt(nc_bndry,directionID,'_FillValue',9.969209968386869e+36);
netcdf.putAtt(nc_bndry,directionID,'axis','Z');

efthID = netcdf.defVar(nc_bndry,'efth','nc_double',[directiondimID frequencydimID stationdimID timedimID]);
netcdf.putAtt(nc_bndry,efthID,'long_name','sea surface wave directional variance spectral density');
netcdf.putAtt(nc_bndry,efthID,'standard_name','sea_surface_wave_directional_variance_spectral_density');
netcdf.putAtt(nc_bndry,efthID,'globwave_name','directional_variance_spectral_density');
netcdf.putAtt(nc_bndry,efthID,'units','m2 s rad-1');
netcdf.putAtt(nc_bndry,efthID,'scale_factor',1);
netcdf.putAtt(nc_bndry,efthID,'add_offset',0);
netcdf.putAtt(nc_bndry,efthID,'valid_min',0);
netcdf.putAtt(nc_bndry,efthID,'valid_max',1.000000020040877e+20);
netcdf.putAtt(nc_bndry,efthID,'_FillValue',9.969209968386869e+36);
netcdf.putAtt(nc_bndry,efthID,'content','TXYZ');
netcdf.putAtt(nc_bndry,efthID,'associates','time station frequency direction');

%close file
netcdf.close(nc_bndry)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now fill the data

times = [3652 3653] ;
station  = 1
station_name ='north'
longitude = 0;
latitude = 0;
% copy l ine from ww3_grid.inp
%$ Frequency increment factor and first frequency (Hz) ---------------- $
%$ number of frequencies (wavenumbers) and directions, relative offset
%$ of first direction in terms of the directional increment [-0.5,0.5].
%$ In versions 1.18 and 2.22 of the model this value was by definiton 0,
%$ it is added to mitigate the GSE for a first order scheme. Note that
%$ this factor is IGNORED in the print plots in ww3_outp.
%$
%   1.1  0.04118  30  30  6.
%
freq(1) = 0.04118;%  from ww3_grid line
for mm=2:frequency
  freq(mm)=freq(mm-1)*1.1;
end

freq1(1)=freq(1);
for mm=2:frequency
  freq1(mm)=freq(mm)*0.9545;
end

for mm=1:frequency-1
  freq2(mm)=freq(mm)*1.05;
end
freq2(frequency)=freq(end);

clear dir
dir(1)=6;             %  from ww3_grid line
dtheta=360/direction
for mm=2:direction
  dir(mm)=dir(mm-1)+dtheta;
end
dir

%         --  freqs---
%       d
%       i
%       r
%       s
%
scale=  0.001;
%scale=  0.01;
data=[
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     1    25   203   775  1795  3526  8828  9713  4278  2458  1797  1296   906   619   416   277   182   120    78    51    33    21    14     9     6     4 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 ];
efth=zeros([direction frequency station time]);
efth(:,:,1,1)=data.*scale;
efth(:,:,1,2)=data.*scale;
x=[7500 7500];
y=[14000 14000];

%now write the data out
ncwrite(fn,'time',times)
ncwrite(fn,'station',station)
ncwrite(fn,'station_name',station_name)
ncwrite(fn,'longitude',longitude)
ncwrite(fn,'latitude',latitude)
ncwrite(fn,'x',x)
ncwrite(fn,'y',y)
ncwrite(fn,'frequency',freq)
ncwrite(fn,'frequency1',freq1)
ncwrite(fn,'frequency2',freq2)
ncwrite(fn,'direction',dir)
ncwrite(fn,'efth',efth)

% thenedit ww3_bounc.inp and  run ww3_bounc
% this will create a nest.ww3


