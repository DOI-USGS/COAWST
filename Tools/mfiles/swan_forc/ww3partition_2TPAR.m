function ww3partition_2TPAR(partfile,specpts,yearww3,mmww3)
%
% reads a Wave Watch III partition file 
% to extract wave paramaters near a grid edge.
%
% 04Feb2019: jcwarner 
%
%  Based on the approach presented here:
%  Kumar, N., D.L. Cahl, S.C. Crosby, and G. Voulgaris, 2017: 
%   Bulk versus Spectral Wave Parameters: Implications on Stokes Drift Estimates,
%   Regional Wave Modeling, and HF Radars Applications. J. Phys. Oceanogr., 47, 
%   1413–1431, https://doi.org/10.1175/JPO-D-16-0203.1 

%
% First get all the grid points.
%
fid=fopen(partfile);
tidx=0;
clear lon lat
tic
while 1
  tline = fgetl(fid);
  if (strcmp(tline(2:9),[yearww3,mmww3,'01']))
    tidx=tidx+1;
    lat(tidx)=str2num(tline(18:24));
    lon(tidx)=str2num(tline(26:32));
  end
  if (strcmp(tline(2:9),[yearww3,mmww3,'02']))
    break;
  end
  if ~ischar(tline), break, end
end
fclose(fid);
disp(['There are ', num2str(tidx), ' number of points'])
toc
%
% Now pick which lon/lat pairs to use for each specpt
%
for ii=1:size(specpts,1)
  for mm=1:tidx
   [dist(mm),phaseangle] = sw_dist([specpts(ii,2),lat(mm)],[specpts(ii,1),lon(mm)]);
  end
  [B,I]=sort(dist);
  ww3_index(ii)=I(1);
end
figure
plot(lon,lat,'+')
hold on
plot(lon(ww3_index),lat(ww3_index),'r+')
title('WW3 points in blue, and the output pts in red')
%
% Now read thru the file, and when we get to a point
% that we want to save, then write those lines to 
% a temp file.
%
fname='temp_partfile.dat';
fid=fopen(partfile);
fidout=fopen(fname,'w');
tic
while 1
  tline = fgetl(fid);
  if ~ischar(tline), break, end
  if (strcmp(tline(2:5),yearww3))
    [C,IA,IB]=intersect(str2num(tline(18:24)),lat(ww3_index));
    if (~isempty(C))
      if intersect(str2num(tline(26:32)),lon(ww3_index(IB)))
        disp(['writing ', tline])
        tline(35:46)='   99999    ';    % replace 'grid_point' with a flag
        fprintf(fidout,'%s\n',tline);
        for mm=1:str2num(tline(48:49))+1
          tline = fgetl(fid);
          fprintf(fidout,'%s\n',tline);
        end
      end
    end
  end
end
fclose(fid);
fclose(fidout);
toc
%
% Now create the full spectra from the partition data 
% that we saved in the temp file.
% Loop for each point.
%
A     = importdata(fname);
for mm=1:length(ww3_index)
  Sfth=[];
  disp(['Working on position ',num2str(mm),' of ', num2str(length(ww3_index)), ...
        '  at lon = ',num2str(lon(ww3_index(mm))),' lat = ', num2str(lat(ww3_index(mm)))])
  Prt   = A(:,5);    % flag for start of partition data
  Lon   = A(:,4);
  Lat   = A(:,3);
% j will be the indices in A that have the lat/lon we want and are the main data lines with the time stamp.
  j     = find(Prt==99999 & Lon==lon(ww3_index(mm)) & Lat==lat(ww3_index(mm)) );
  toadd = ['00',num2str(mm)]; toadd=toadd(end-1:end);
  sname = ['WW3_Partitions','_',toadd,'.mat'];
  if (~isempty(j))
%
% Compute the combined spectra.
%
    extract_and_save_WW3_partitioned_data
%
% Write out the 2D spec files.
%
    Specname = ['WW3part','_',toadd,'.spc2d'];
    createswan2Dspec(YYYYMODD,HHMMSS,[lon(ww3_index(mm)) lat(ww3_index(mm))],Freq,Dir,Sfth,Specname)
  else
    disp(['no data for position ', num2str(mm),' at lat lon of ',num2str(lat(ww3_index(mm))),' ',num2str(lon(ww3_index(mm))),' ']);
  end
%
end

disp('here')



 