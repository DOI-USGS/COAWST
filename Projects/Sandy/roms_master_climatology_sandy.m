% roms_master_climatology_coawst_mw  
% Edited to be specific to the Sandy test case.
%
% This routine :
%  - creates climatology, boundary, and initial condition files for ROMS: 
%    coawst_clm.nc ; coawst_bdy.nc ; coawst_ini.nc 
%    on a user-defined grid for a user-defined date.
%
% This is currently set up to use opendap calls to acquire data
%  from HYCOM + NCODA Global 1/12 Degree Analysis and interp to US_East grid.
%
% Before running this routine, user needs to setup "nctoolbox" within Matlab.
%  
% based on efforts by:
% written by Mingkui Li, May 2008
% Modified by Brandy Armstrong March 2009
% jcwarner April 20, 2009
% Ilgar Safak modified on June 27, 2012 such that now:
% - HYCOM url is a user-definition
% - "hc" is called from the structure "gn".(still needs to be tested with wet/dry).
% - updatinit_coawst_mw.m modified to get desired time (T1) as a variable;
%    ocean_time=T1-datenum(1858,11,17,0,0,0)

%%%%%%%%%%%%%%%%%%%%%   START OF USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%

% (1) Enter start date (T1) and number of days to get climatology data 
T1=datenum(2012,10,28,12,0,0); %start date
%number of days to create clm for
numdays=4;
dayFrequency=1;

% (2) Enter URL of the HYCOM catalog for the requested time, T1
%     see http://tds.hycom.org/thredds/catalog.html
url='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9';      % Jan 2011 - Present

% (3) Enter working directory (wdr)
wdr='c:\work\models\COAWST\Projects\Sandy2';

% (4) Enter path and name of the ROMS grid (modelgrid)
%modelgrid='Sandy_roms_grid.nc'
modelgrid='Sandy_roms_grid.nc'

% (5) Enter grid vertical coordinate parameters --These need to be consistent with the ROMS setup. 
theta_s=5;
theta_b=0.4;
Tcline=50;
N=16;

%%%%%%%%%%%%%%%%%%%%%   END OF USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%
eval(['cd ',wdr])
eval(['gridname=''',modelgrid,''';']);

disp('getting roms and hycom grid dimensions ...');
gn=roms_get_grid_mw(gridname,[theta_s theta_b Tcline N]);

tic

% Call to get HYCOM indices for the defined ROMS grid
disp('getting hycom indices')
get_ijrg(gn,url)

% Call to create the climatology (clm) file
disp('going to create clm file')
fn=updatclim_coawst_mw(T1,gn,'coawst_clm.nc',wdr,url)

% Call to create the boundary (bdy) file
disp('going to create bndry file')
updatbdry_coawst_mw(fn,gn,'coawst_bdy.nc',wdr)

% Call to create the initial (ini) file
disp('going to create init file')
updatinit_coawst_mw(fn,gn,'coawst_ini.nc',wdr,T1)

%% Call to create the long climatology (clm) file
if numdays>1
    disp('going to create more days of clm and bnd files')
    eval(['!copy coawst_clm.nc coawst_clm_',datestr(T1,'yyyymmdd'),'.nc'])
    eval(['!copy coawst_bdy.nc coawst_bdy_',datestr(T1,'yyyymmdd'),'.nc'])
    for it=dayFrequency:dayFrequency:numdays-1      %1st day already created, NEED to set number of days at top!
        fname=['coawst_clm_',datestr(T1+it,'yyyymmdd'),'.nc']
        fn=updatclim_coawst_mw(T1+it,gn,fname,wdr,url)
        fname=['coawst_bdy_',datestr(T1+it,'yyyymmdd'),'.nc'];
        updatbdry_coawst_mw(fn,gn,fname,wdr)
    end
    %% get an organized list of dated files
    Dclm=dirsort('coawst_clm_*.nc');
    Dbdy=dirsort('coawst_bdy_*.nc');
    %names for merged climatology/boundary files
    fout='merged_coawst_clm.nc';
    foutb='merged_coawst_bdy.nc';
    %create netcdf files to merge climatology into
    create_roms_netcdf_clm_mwUL(fout,gn,length(Dclm));% converted to BI functions
    create_roms_netcdf_bndry_mwUL(foutb,gn,length(Dbdy));% converted to BI functions
    %% fill merged climatology files with data from each clm file
    % each file must contain only ONE time step
    %get variable names
    vinfo=ncinfo(fout);
    for nf=1:length(Dclm)
        fin=Dclm(nf).name;
        for nv=1:length({vinfo.Variables.Name})
            if length({vinfo.Variables(nv).Dimensions.Name})==4;
                eval(['ncwrite(fout,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[1 1 1 nf]);']);
            elseif length({vinfo.Variables(nv).Dimensions.Name})==3;
                eval(['ncwrite(fout,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[1 1 nf]);']);
            elseif length({vinfo.Variables(nv).Dimensions.Name})==2;
                try
                    eval(['ncwrite(fout,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[1 nf]);']);
                catch
                    display([vinfo.Variables(nv).Name ' is a dimension and has already been written to the file.'])
                end
            elseif length({vinfo.Variables(nv).Dimensions.Name})==1;
                try
                    eval(['ncwrite(fout,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[nf]);']);
                catch
                    display([vinfo.Variables(nv).Name ' is a dimension and has already been written to the file.'])
                end
            end
        end
    end
    
    vinfo=ncinfo(foutb);
    for nf=1:length(Dbdy)
        for nv=1:length({vinfo.Variables.Name})
            fin=Dbdy(nf).name;
            if length({vinfo.Variables(nv).Dimensions.Name})==4;
                eval(['ncwrite(foutb,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[1 1 1 nf]);']);
            elseif length({vinfo.Variables(nv).Dimensions.Name})==3;
                eval(['ncwrite(foutb,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[1 1 nf]);']);
            elseif length({vinfo.Variables(nv).Dimensions.Name})==2;
                try
                    eval(['ncwrite(foutb,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[1 nf]);']);
                catch
                    display([vinfo.Variables(nv).Name ' is a dimension and has already been written to the file.'])
                end
                
            elseif length({vinfo.Variables(nv).Dimensions.Name})==1;
                try
                    eval(['ncwrite(foutb,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[nf]);']);
                catch
                    display([vinfo.Variables(nv).Name ' is a dimension and has already been written to the file.'])
                end
            end
        end
    end
end

toc

