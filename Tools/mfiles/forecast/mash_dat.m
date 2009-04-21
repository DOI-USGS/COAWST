function mash_dat(modelgrid,windstart,windstart_gfs);
cd D:\Temporary\COAWST\SWAN\forcing
 %Grid to interpolate to
    ncg=netcdf(modelgrid);
    [eta xi]=size(ncg{'lon_rho'}(:));
    ncclose
    
%open wind data
fid_gfs = fopen('gfswind_update.dat');
fid_nam = fopen('namwind_update.dat');
fid_comb = fopen('comb_wind_update.dat','w'); %name of data file to write

%figure out if there is a time difference between start time for dat files
dt=str2num(windstart)-str2num(windstart_gfs);
if dt~=0
    if dt>0
        for i=1:dt*2*8
            [A_gfs,count_gfs]=fscanf(fid_gfs,'%f',xi*eta);
        end
    elseif dt<0
        for i=1:(-dt)*2*8
            [A_nam,count_nam]=fscanf(fid_nam,'%f',xi*eta);
        end
    end
else
end

%combine grids
for tidx=1:9999999
  tidx
  %read data
  [A_nam,count_nam]=fscanf(fid_nam,'%f',xi*eta);
  [A_gfs,count_gfs]=fscanf(fid_gfs,'%f',xi*eta);
  if count_nam==0; break; end
  %reshape data
  wind_temp_nam=reshape(A_nam,xi,eta);
  wind_temp_nam=wind_temp_nam';
  
  wind_temp_gfs=reshape(A_gfs,xi,eta);
  wind_temp_gfs=wind_temp_gfs';

  %interpolate wind outside of ETA/NAM using global data
  wind_temp_comb=smoothwinds_1(wind_temp_nam,wind_temp_gfs);

  %check
%   pcolorjw(wind_temp_nam);
%   pause(1);
  
  %write wind
  for index = 1:eta;
      fprintf(fid_comb,'%3.2f\n',wind_temp_comb(index,:));
  end
  %fprintf(fid_comb,'%f\n',wind_temp_nam');
end
%close files
fclose(fid_gfs);
fclose(fid_nam);
fclose(fid_comb);

