function updatinit_coawst_mw(fn,gn,ini,wdr,T1)
% script create_roms_init
%
% Create a netcdf file that contains initialization data for ROMS.
% Initializes temp, salt, u, v, ubar, vbar, and 
% all sediment parameters.
%
% In cppdefs.h you should have 
% #undef ana_initial
% #undef ana_sediment
% This m file is set to work from roms_master_climatology.
%
% jcw 5-25-2005
% jcw 3-7-07 add L and M for get grid
% 28July2016 Big update to be consistent with create_roms_init.
%

%1) Enter name of netcdf initial file to be created.
%   If it already exists it will be overwritten!!.
    init_file=[ini];

%2) If you want to read the initial data from another netcdf file, then 
%   enter read_int=1 and provide a name of history or rst file to read data from.
%   Otherwise set read_int=0.
    read_init=1;
%   Use clm file to get the data
    nc_clm=netcdf.open(fn,'NC_NOWRITE');

%3) Enter start time of initial file, in seconds.
%   This time needs to be consistent with model time (ie dstart and time_ref).
%   See *.in files for more detail. 
    if (read_init)
      ocean_time=T1-datenum(1858,11,17,0,0,0);
    else
      ocean_time=0;
    end

%4) Set values of theta_s, theta_b, Tcline, and N from your *.in file.
    theta_s   = gn.theta_s;
    theta_b   = gn.theta_b;
    Tcline    = gn.Tcline;
    N         = gn.N;
    Vtransform= gn.Vtransform;
    Vstretching=gn.Vstretching;

%5) Obtain grid information.
   h=gn.h;
   hmin=min(h(:));
   hc=min([hmin,Tcline]);
   [LP,MP]=size(gn.lon_rho);
   L=LP-1;
   Lm=L-1;
   M=MP-1;
   Mm=M-1;
   L  = Lm+1;
   M  = Mm+1;
   xi_psi  = L;
   xi_rho  = LP;
   xi_u    = L;
   xi_v    = LP;
   eta_psi = M;
   eta_rho = MP;
   eta_u   = MP;
   eta_v   = M;
   s       = gn.N;
   Cs_r=gn.Cs_r;
   Cs_w=gn.Cs_w;
   sc_r=gn.s_rho;
   sc_w=gn.s_w;
    
%6) Initialize zeta, salt, temp, u, v, ubar, vbar.
%   moved this to fill part below
%% 
NAT=2;            % Number of active tracers. Usually 2 (temp + salt).

% 
%7) Enter number of mud sediments (NCS) and number of sand sediments (NNS).
%   These values should be the same as in mod_param.F
    NCS = 0;   %number of cohesive sed classes
    NNS = 1;   %number of non-cohesive sed classes
%
% calc sed parameters. Do not alter.
%
   NST = NCS + NNS;     % total number of sed tracers.
   NT = NAT+NST;        % total number of tracers.

%8) Enter number of bed layers
%   This value should be the same as in mod_param.F
    Nbed = 1;

%9) Sediment class properties (in order, mud first then sand).
%  These values should coincide with your sediment.in file.
  mud_Srho=ones(1,NCS)*2650;        %kg m-3, NCS values
  mud_Sd50=[0.01 0.06]/1000;        %m,      NCS values
  mud_Wsed=[0.10 1.0]/1000;         %m s-1,  NCS values
  mud_tau_ce=[0.05 0.05];           %N m-2,  NCS values
  mud_Erate=[5 5]*1e-5;             %kg m-2 s-1, NCS values
  sand_Srho=ones(1,NNS)*2650;       %kg m-3, NNS values
  sand_Sd50=[1.0]/1000;             %m,      NNS values
  sand_Wsed=[1.0]/1000;             %m s-1,  NNS values
  sand_tau_ce=[0.07];               %N m-2,  NNS values
  sand_Erate=[1]*1e-4;              %kg m-2 s-1, NNS values
%
% make some combined arrays.  Do not alter.
%
  Srho=  [mud_Srho,sand_Srho];
  Sd50=  [mud_Sd50,sand_Sd50];
  Wsed=  [mud_Wsed,sand_Wsed];
  tau_ce=[mud_tau_ce,sand_tau_ce];
  Erate= [mud_Erate,sand_Erate];


%9) Provide initial sediment properties in water column.
%
  display('Initializing suspended sediments.')
%
% mud.
%
for idmud=1:NCS
  count=['0',num2str(idmud)];
  count=count(end-1:end);
  eval(['mud_',count,'(1:xi_rho,1:eta_rho,1:N,1:length(ocean_time)) = 0;'])   %mud conc in water column
end
%
% sand.
%
for isand=1:NNS
  count=['0',num2str(isand)];
  count=count(end-1:end);
  eval(['sand_',count,'(1:xi_rho,1:eta_rho,1:N,1:length(ocean_time)) = 0;'])  %sand conc in water column
end

%10) Provide initial sediment properties in bed.
%
% bed properties
  display('Initializing sediment bed.')
%
for k=1:Nbed
  for time=1:length(ocean_time)
    bed_thickness(1:xi_rho,1:eta_rho,k,time) = 0.1;
    bed_age(1:xi_rho,1:eta_rho,k,time)       = ocean_time(1);
    bed_porosity(1:xi_rho,1:eta_rho,k,time)  = 0.9;
    bed_biodiff(1:xi_rho,1:eta_rho,k,time)   = 0.0;
  end
end
%
% for mud
%
for idsed=1:NCS
  count=['0',num2str(idsed)];
  count=count(end-1:end);
  for k=1:Nbed
    for time=1:length(ocean_time)
      eval(['mudfrac_',count,'(1:xi_rho,1:eta_rho,k,time) = 1/NST;'])      %fraction of each sed class in each bed cell
      eval(['mudmass_',count,'(1:xi_rho,1:eta_rho,k,time) = squeeze(bed_thickness(1:xi_rho,1:eta_rho,k,time)).*Srho(idsed).*(1.0-squeeze(bed_porosity(1:xi_rho,1:eta_rho,k,time))).*squeeze(mudfrac_',count,'(1:xi_rho,1:eta_rho,k,time));'])          %mass of each sed class in each bed cell
    end
  end
end
%
% for sand
%
for isand=1:NNS
 count=['0',num2str(isand)];
 count=count(end-1:end);
  for k=1:Nbed
    for time=1:length(ocean_time)
      eval(['sandfrac_',count,'(1:xi_rho,1:eta_rho,k,time) = 1/NST;'])      %fraction of each sed class in each bed cell
      eval(['sandmass_',count,'(1:xi_rho,1:eta_rho,k,time) = squeeze(bed_thickness(1:xi_rho,1:eta_rho,k,time)).*Srho(isand).*(1.0-squeeze(bed_porosity(1:xi_rho,1:eta_rho,k,time))).*squeeze(sandfrac_',count,'(1:xi_rho,1:eta_rho,k,time));'])          %mass of each sed class in each bed cell
      eval(['bedload_Usand_',count,'(1:xi_u,1:eta_u,time) = 0;'])              %bed load
      eval(['bedload_Vsand_',count,'(1:xi_v,1:eta_v,time) = 0;'])              %bed load
    end
  end
end

%11)
%
% set some surface properties
  display('Initializing sediment surface properties.')
%
    cff1=1.0;
    cff2=1.0;
    cff3=1.0;
    cff4=1.0;
    for ised=1:NCS
      count=['0',num2str(ised)];
      count=count(end-1:end);
      eval(['cff1=cff1*mud_Sd50(ised)^squeeze(mudfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
      eval(['cff2=cff2*mud_Srho(ised)^squeeze(mudfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
      eval(['cff3=cff3*mud_Wsed(ised)^squeeze(mudfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
      eval(['cff4=cff4*mud_tau_ce(ised)^squeeze(mudfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
    end
    for ised=1:NNS
      count=['0',num2str(ised)];
      count=count(end-1:end);
      eval(['cff1=cff1*sand_Sd50(ised).^squeeze(sandfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
      eval(['cff2=cff2*sand_Srho(ised).^squeeze(sandfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
      eval(['cff3=cff3*sand_Wsed(ised).^squeeze(sandfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
      eval(['cff4=cff4*sand_tau_ce(ised).^squeeze(sandfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
    end
    grain_diameter(1:xi_rho,1:eta_rho,time)=cff1;
    grain_density(1:xi_rho,1:eta_rho,time)=cff2;
    settling_vel(1:xi_rho,1:eta_rho,time)=cff3;
    erosion_stress(1:xi_rho,1:eta_rho,time)=cff4;
    ripple_length(1:xi_rho,1:eta_rho,time)=0.10;
    ripple_height(1:xi_rho,1:eta_rho,time)=0.01;
    dmix_offset(1:xi_rho,1:eta_rho,time)=0.0;
    dmix_slope(1:xi_rho,1:eta_rho,time)=0.0;
    dmix_time(1:xi_rho,1:eta_rho,time)=0.0;

%12)
%
% set vegetation properties 
  NVEG=1; 
  plant_density=zeros(xi_rho,eta_rho,NVEG);
  plant_height=zeros(xi_rho,eta_rho,NVEG);
  plant_diameter=zeros(xi_rho,eta_rho,NVEG);
  plant_thickness=zeros(xi_rho,eta_rho,NVEG);
  marsh_mask=zeros(xi_rho,eta_rho,NVEG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now create the netcdf file first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

create_roms_netcdf_init_mw(init_file,gn,Nbed,NNS,NCS,NVEG)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now write the data from the arrays to the netcdf file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nc_init=netcdf.open(init_file,'NC_WRITE');

disp(' ## Filling Variables in netcdf file with data...')

% copy clm time to init time
  tempid = netcdf.inqVarID(nc_init,'ocean_time');  %get id
  netcdf.putVar(nc_init,tempid,ocean_time);        %set variable

%
% *** Init values, ONLY WORKS FOR VALUES ON SAME GRID ****************
% *** ADD Interpolation if grids are not the same ********************
vars2d={'zeta','ubar','vbar'};
vars3d={'u','v','temp','salt'};
%% 2D variables
for i=1:length(vars2d)
    eval(['tempclmid = netcdf.inqVarID(nc_clm,''',vars2d{i},''');']);%get id
    eval([vars2d{i},'=netcdf.getVar(nc_clm,tempclmid);']);%get data   
    eval(['tempid = netcdf.inqVarID(nc_init,''',vars2d{i},''');']);%get id
    eval(['netcdf.putVar(nc_init,tempid,',vars2d{i},');']);%set variable
end
%% 3D variables
for i=1:length(vars3d)
    eval(['tempclmid = netcdf.inqVarID(nc_clm,''',vars3d{i},''');']);%get id
    eval([vars3d{i},'=netcdf.getVar(nc_clm,tempclmid);']);%get data
    eval(['tempid = netcdf.inqVarID(nc_init,''',vars3d{i},''');']);%get id
    eval(['netcdf.putVar(nc_init,tempid,',vars3d{i},');']);%set variable
    clear temp3 tempv tempt
end

spherical=1;
morvars={'theta_s','theta_b','Tcline','Cs_r','Cs_w','sc_w','sc_r','hc', ...
         'Vtransform','Vstretching','spherical','bed_thickness',...
         'bed_age','bed_porosity','bed_biodiff','grain_diameter', ...
         'grain_density','settling_vel','erosion_stress',...
         'ripple_height','ripple_length','dmix_offset','dmix_slope','dmix_time', ...
         'ocean_time'};
for i=1:length(morvars)
    eval(['tempid = netcdf.inqVarID(nc_init,''',morvars{i},''');']);%get id
    eval(['netcdf.putVar(nc_init,tempid,',morvars{i},');']);%set variable
    eval(['clear ',morvars{i},';']);
end

for mm=1:NCS
    count=['00',num2str(mm)];
    count=count(end-1:end);
    
    eval(['tempid = netcdf.inqVarID(nc_init,''mud_',count,''');']);%sed conc in water column
    eval(['netcdf.putVar(nc_init,tempid,mud_',count,');']);
    
    eval(['tempid = netcdf.inqVarID(nc_init,''mudfrac_',count,''');']); %fraction of each sed class in each bed cell
    eval(['netcdf.putVar(nc_init,tempid,mudfrac_',count,');']);
    
    eval(['tempid = netcdf.inqVarID(nc_init,''mudmass_',count,''');']);%mass of each sed class in each bed cell
    eval(['netcdf.putVar(nc_init,tempid,mudmass_',count,');']);
end
for mm=1:NNS
    count=['00',num2str(mm)];
    count=count(end-1:end);
    
    eval(['tempid = netcdf.inqVarID(nc_init,''sand_',count,''');']);%sed conc in water column
    eval(['netcdf.putVar(nc_init,tempid,sand_',count,');']);
    
    eval(['tempid = netcdf.inqVarID(nc_init,''sandfrac_',count,''');']); %fraction of each sed class in each bed cell
    eval(['netcdf.putVar(nc_init,tempid,sandfrac_',count,');']);
    
    eval(['tempid = netcdf.inqVarID(nc_init,''sandmass_',count,''');']);%mass of each sed class in each bed cell
    eval(['netcdf.putVar(nc_init,tempid,sandmass_',count,');']);

    eval(['tempid = netcdf.inqVarID(nc_init,''bedload_Usand_',count,''');']);%mass of each sed class in each bed cell
    eval(['netcdf.putVar(nc_init,tempid,bedload_Usand_',count,');']);
    eval(['tempid = netcdf.inqVarID(nc_init,''bedload_Vsand_',count,''');']);%mass of each sed class in each bed cell
    eval(['netcdf.putVar(nc_init,tempid,bedload_Vsand_',count,');']);
end

%Veg
ncwrite(init_file,'plant_height',plant_height);
ncwrite(init_file,'plant_diameter',plant_diameter);
ncwrite(init_file,'plant_density',plant_density);
ncwrite(init_file,'plant_thickness',plant_thickness);
ncwrite(init_file,'marsh_mask',marsh_mask);

%close file
netcdf.close(nc_init)
netcdf.close(nc_clm)



