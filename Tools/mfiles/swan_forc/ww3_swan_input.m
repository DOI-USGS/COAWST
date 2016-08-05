% script ww3_swan_input.m
%
% This script is the main driver to 
% download, convert and interpolate WW3 data to TPAR input for SWAN.
%
% Written 4/23/09 by Brandy Armstrong
% some mods, jcwarner Arpil 27, 2009
%
% READ THE INSTRUCTIONS IN THE COAWST MANUAL 
% FOR SWAN BC's SECTION 10.
%
% First, acquire the necessary grib files from
% ftp://polar.ncep.noaa.gov/pub/history/waves/
%

% ************* BEGIN USER INPUT   ****************************

%1) Enter WORKING DIRECTORY.
% This is the location of ww3 grb files downloaded and the
% location of output for the TPAR files to be created.
% ***WARNING***
% The TPAR files created are saved in the working directory and are named
% generically (numbered). Any existing TPAR files will be overwritten !!!!
%
working_dir='D:\COAWST_tests\sandy_for_workshop\Projects\Sandy'
eval(['cd ',working_dir,';']);

%2) Enter dates of data requested.
yearww3='2012';    %input year of data yyyy 
mmww3='10';       %input month of data mm
ddww3='00';        %keep this as '00'

%3) Enter the ww3 grid area
ww3_area='multi_1.at_10m';    %western north atlantic

%4) Enter path\name of SWAN netcdf grid. This is typically the same
% as the roms grid.
modelgrid='D:/COAWST_tests/sandy_for_workshop/Projects/Sandy/Sandy_roms_grid.nc';

%5) Enter the spacings of TPAR file locations around the perimeter
% of the grid. One TPAR file every 'specres' point.
% ww3_specpoints assumes masking of 0 for land and NaN for water
specres=10; % spec point resolution

% flag for simulations that change month (e.g., goes from October 10 to
% November 20)
long_run=0;  % DEFAULT
%long_run=1;
if long_run
    total_number_months=2; % total number of months 
                        % (not the length of the run, but rather the number
                        % of grib files)
end

% *************END OF USER INPUT****************************

if ~long_run
    % Call routine to compute TPAR files.
    ww3gb_2TPAR(modelgrid,yearww3,mmww3,ww3_area,ddww3,specres)

    % After creating the TPAR files, tell the user what info is needed to 
    % be added to INPUT file.
    % Write out boundary file lines for INPUT file
    bdry_com %script writes out file Bound_spec_command to working directory
    display('BOUNDSPEC command lines can be found in the file Bound_spec_command');
    display('You must copy the lines from that file into your SWAN INPUT file');
else
    eval(['!mkdir ',yearww3,mmww3])
    
    % Call routine to compute TPAR files.
    ww3gb_2TPAR(modelgrid,yearww3,mmww3,ww3_area,ddww3,specres)

    % After creating the TPAR files, tell the user what info is needed to 
    % be added to INPUT file.
    % Write out boundary file lines for INPUT file
    bdry_com %script writes out file Bound_spec_command to working directory
    display('BOUNDSPEC command lines can be found in the file Bound_spec_command');
    display('You must copy the lines from that file into your SWAN INPUT file');
    
    eval(['!mv TPAR*.txt ',yearww3,mmww3,'/.'])

    for im=1:total_number_months-1
        yymm=datestr(datenum(str2double(yearww3),str2double(mmww3)+im,1),'yyyymm');
        eval(['!mkdir ',yymm])

        % Call routine to compute TPAR files.
        ww3gb_2TPAR(modelgrid,yymm(1:4),yymm(5:6),ww3_area,ddww3,specres)

        bdry_com %script writes out file Bound_spec_command to working directory
        display('BOUNDSPEC command lines can be found in the file Bound_spec_command');
        display('You must copy the lines from that file into your SWAN INPUT file');
        
        eval(['!mv TPAR*.txt ',yymm,'/.'])
    end   
    
    D=dir([yearww3,mmww3,'/TPAR*.txt']);
    dn={D.name}';
    [y,i]=sort(dn);
    D=D(i);
    
    for itp=1:length(D)
        data=[];
        [pfid,message]=fopen([yearww3,mmww3,'/',D(itp).name]);
        tline = fgetl(pfid);
        while 1
            tline = fgetl(pfid);
            if ~ischar(tline), break, end
            data=[data;str2num(tline)];
        end
%         dum=fgets(pfid);
%         data=fscanf(pfid,'%f %f %f %f %f',[5 inf])';        
%         C=textscan(pfid,'','headerlines',1);
%         if ~isempty(find(isnan(C{3}),1))
%             C=textscan(pfid,'%f%f%f%f.%f','headerlines',1);
%         end
%         data=cell2mat(C);
        fclose(pfid);
        for im=1:total_number_months-1
            yymm=datestr(datenum(str2double(yearww3),str2double(mmww3)+im,1),'yyyymm');
            [pfid,message]=fopen([yymm,'/',D(itp).name]);
            tline = fgetl(pfid);
            while 1
                tline = fgetl(pfid);
                if ~ischar(tline), break, end
                data=[data;str2num(tline)];
            end
            fclose(pfid);
        end           
        [CA,IA,IC] = unique(data(:,1));
        data=data(IA,:);
        
        [pfid,message]=fopen(D(itp).name,'w');
        fprintf(fid,'TPAR \n');
        fprintf(pfid,'%8.4f %3.2f %3.2f %3.f. %2.f\n',data');
        fclose(pfid);
    end           
end

