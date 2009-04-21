function updatbdry(fn,gn)

% written by Mingkui Li, May 2008
% further adaptions by jcw, April 18, 2009

%Use clm file to get the data
nc_clm=netcdf(fn);
%get number of time steps in clm file
ocean_time=nc_clm{'ocean_time'}(:);

%file name to be created
bndry_file='USE_bndry.nc';
%call to create this bndry file.
create_roms_netcdf_bndry(bndry_file,gn,ocean_time)
nc_bndry=netcdf(bndry_file,'w');

%now write the data from the arrays to the netcdf file
disp(' ## Filling Variables in netcdf file with data...')

nc_bndry{'zeta_time'}(:) = ocean_time;
nc_bndry{'v2d_time'}(:) = ocean_time;
nc_bndry{'v3d_time'}(:) = ocean_time;
nc_bndry{'salt_time'}(:) = ocean_time;
nc_bndry{'temp_time'}(:) = ocean_time;

%zeta
if (length(ocean_time)==1)
  zeta_south=nc_clm{'zeta'}(1,:);
else
  zeta_south=nc_clm{'zeta'}(:,1,:);
end
nc_bndry{'zeta_south'}(:) = zeta_south;
clear zeta_south

if (length(ocean_time)==1)
  zeta_east=nc_clm{'zeta'}(:,end);
else
  zeta_east=nc_clm{'zeta'}(:,:,end);
end
nc_bndry{'zeta_east'}(:) = zeta_east;
clear zeta_east

%ubar
if (length(ocean_time)==1)
  ubar_south=nc_clm{'ubar'}(1,:);
else
  ubar_south=nc_clm{'ubar'}(:,1,:);
end
nc_bndry{'ubar_south'}(:) = ubar_south;
clear ubar_south

if (length(ocean_time)==1)
  ubar_east=nc_clm{'ubar'}(:,end);
else
  ubar_east=nc_clm{'ubar'}(:,:,end);
end
nc_bndry{'ubar_east'}(:) = ubar_east;
clear ubar_east

%vbar
if (length(ocean_time)==1)
  vbar_south=nc_clm{'vbar'}(1,:);
else
  vbar_south=nc_clm{'vbar'}(:,1,:);
end
nc_bndry{'vbar_south'}(:) = vbar_south;
clear vbar_south

if (length(ocean_time)==1)
  vbar_east=nc_clm{'vbar'}(:,end);
else
  vbar_east=nc_clm{'vbar'}(:,:,end);
end
nc_bndry{'vbar_east'}(:) = vbar_east;
clear vbar_east

%u
if (length(ocean_time)==1)
  u_south=nc_clm{'u'}(:,1,:);
else
  u_south=nc_clm{'u'}(:,:,1,:);
end
nc_bndry{'u_south'}(:) = u_south;
clear u_south

if (length(ocean_time)==1)
  u_east=nc_clm{'u'}(:,:,end);
else
  u_east=nc_clm{'u'}(:,:,:,end);
end
nc_bndry{'u_east'}(:) = u_east;
clear u_east

%v
if (length(ocean_time)==1)
  v_south=nc_clm{'v'}(:,1,:);
else
  v_south=nc_clm{'v'}(:,:,1,:);
end
nc_bndry{'v_south'}(:) = v_south;
clear v_south

if (length(ocean_time)==1)
  v_east=nc_clm{'v'}(:,:,end);
else
  v_east=nc_clm{'v'}(:,:,:,end);
end
nc_bndry{'v_east'}(:) = v_east;
clear v_east

%temp
if (length(ocean_time)==1)
  temp_south=nc_clm{'temp'}(:,1,:);
else
  temp_south=nc_clm{'temp'}(:,:,1,:);
end
nc_bndry{'temp_south'}(:) = temp_south;
clear temp_south

if (length(ocean_time)==1)
  temp_east=nc_clm{'temp'}(:,:,end);
else
  temp_east=nc_clm{'temp'}(:,:,:,end);
end
nc_bndry{'temp_east'}(:) = temp_east;
clear temp_east

%salt
if (length(ocean_time)==1)
  salt_south=nc_clm{'salt'}(:,1,:);
else
  salt_south=nc_clm{'salt'}(:,:,1,:);
end
nc_bndry{'salt_south'}(:) = salt_south;
clear salt_south

if (length(ocean_time)==1)
  salt_east=nc_clm{'salt'}(:,:,end);
else
  salt_east=nc_clm{'salt'}(:,:,:,end);
end
nc_bndry{'salt_east'}(:) = salt_east;
clear salt_east

close(nc_clm)
close(nc_bndry)

