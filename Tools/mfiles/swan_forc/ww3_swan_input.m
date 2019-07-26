% script ww3_swan_input.m
%
% This script is the main driver to 
% download, convert and interpolate WW3 data to TPAR input for SWAN.
%
% Written 4/23/09 by Brandy Armstrong
% some mods, jcwarner Arpil 27, 2009
% added partition file call, jcw, 04Feb2019
%
% First, acquire the necessary grib, ascii (gz), or nc files from
% ftp://polar.ncep.noaa.gov/pub/history/waves/
%


% ************* BEGIN USER INPUT   ****************************

% 1) Enter WORKING DIRECTORY.
% This is the location of ww3 files downloaded and the
% location of output for the forcing files to be created.
%
working_dir='F:\data\models\COAWST\Projects\Sandy\ww3'

% 2) Enter start dates of data requested.
yearww3='2012';    %input year of data yyyy 
mmww3='10';        %input month of data mm
ddww3='00';        %keep this as '00'
%
% flag for simulations that span several months. 
% 0 = one month, 1 = more than one month
% This long run flag is only coded for TPAR files for now.
long_run=0;
if long_run
    total_number_months=2; % total number of months 
                        % (not the length of the run, but rather the number
                        % of ww3 files)
end

% 3) Enter path\name of SWAN grid. This is set up to use the roms grid as the same for swan.
modelgrid='F:\data\models\COAWST\Projects\Sandy\Sandy_roms_grid.nc';

% 4) Enter the spacings of the forcing file locations around the perimeter
% of the grid. One forcings file spans between the 'specres' points.
specres=20; % spec point resolution

% 5)
% Here you decide to use gridded field output (to make TPAR) or spectral partition data (to make 2Dspec). Read this:
% ftp://polar.ncep.noaa.gov/pub/history/waves/multi_1/00README
%
% Pick one of these:
use_field_output=0;     % will create TPAR files.
use_partition_data=1;   % will create 2D spec files.
%
if (use_field_output)
  % Enter the ww3 grid area
  ww3_area='multi_1.at_10m';    %western north atlantic
end
if (use_partition_data)
  % Enter name of wave watch 3 partition file
  partfile='multi_1.partition.glo_30m.201210';
end

% *************END OF USER INPUT****************************

eval(['cd ',working_dir,';']);

% call to get the spectral points
[specpts]=ww3_specpoints(modelgrid,specres);

if ~long_run
    % Call routine to compute TPAR files.
    if (use_field_output)
      ww3gb_2TPAR(modelgrid,yearww3,mmww3,ww3_area,ddww3,specpts)
    end
    if (use_partition_data)
      ww3partition_2TPAR(partfile,specpts,yearww3,mmww3)
    end
    % After creating the TPAR files, tell the user what info is needed to 
    % be added to INPUT file.
    % Write out boundary file lines for INPUT file
    bdry_com %script writes out file Bound_spec_command to working directory
    display('BOUNDSPEC command lines can be found in the file Bound_spec_command');
    display('You must copy the lines from that file into your SWAN INPUT file');
else
    eval(['!mkdir ',yearww3,mmww3])
    
    % Call routine to compute TPAR files.
    ww3gb_2TPAR(modelgrid,yearww3,mmww3,ww3_area,ddww3,specpts)

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
        ww3gb_2TPAR(modelgrid,yymm(1:4),yymm(5:6),ww3_area,ddww3,specpts)

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

