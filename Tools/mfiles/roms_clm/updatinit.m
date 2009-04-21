function updatinit(fn,gn)
% script create_roms_init
%
% Create a netcdf file that contains initialization data for ROMS.
% Initializes temp, salt, u, v, ubar, vbar, and 
% all sediment parameters.
%
% In cppdefs.h you should have 
% #undef ana_initial
% #undef ana_sediment
% This m file is set to initalize US_East_grid.nc.
%
% jcw 5-25-2005
% jcw 3-7-07 add L and M for get grid
%

%Use clm file to get the data
nc_clm=netcdf(fn);
%get number of time steps in clm file
t_clim=nc_clm{'ocean_time'}(:);

%1) Enter name of netcdf initial file to be created.
%   If it already exists it will be overwritten!!.
    init_file='USE_init.nc'
%   create_roms_netcdf_init(init_file,gn,t_clim, Nbed, NNS,NCS)
    
%2) Enter start time of initial file, in seconds and time step if file
% has more than one
%   This time needs to be consistent with model time (ie dstart and time_ref).
%   See *.in files for more detail. 
    tidx=1;
    init_time=nc_clm{'ocean_time'}(tidx);

%3) Enter number of vertical sigma levels in model.
%   This will be same value as entered in mod_param.F
    N=gn.N;

%4) Enter the values of theta_s, theta_b, and Tcline from your *.in file.
    theta_s = gn.theta_s;  
    theta_b = gn.theta_b;  
    Tcline =  gn.Tcline;

  [MP,LP]=size(gn.lon_rho);
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

%5)
   hmin=0;
   hc=min([hmin,Tcline]);
   if (theta_s~=0.0)
     cff1=1.0/sinh(theta_s);
     cff2=0.5/tanh(0.5*theta_s);
   end
   sc_w(1)=-1.0;
   Cs_w(1)=-1.0;
   cff=1.0/N;
   for k=1:N
     sc_w(k+1)=cff*(k-N);
     sc_r(k)=cff*((k-N)-0.5);
     if (theta_s~=0)
       Cs_w(k+1)=(1.0-theta_b)*cff1*sinh(theta_s*sc_w(k+1))+   ...
                      theta_b*(cff2*tanh(theta_s*(sc_w(k+1)+0.5))-0.5);
       Cs_r(k)  =(1.0-theta_b)*cff1*sinh(theta_s*sc_r(k))+   ...
                      theta_b*(cff2*tanh(theta_s*(sc_r(k)+0.5))-0.5);
     else
       Cs_w(k+1)=sc_w(k+1);
       Cs_r(k)=sc_r(k)
      end
    end
 
    
%6) Initialize zeta, salt, temp, u, v, ubar, vbar.
%
% Init values for zeta.
  display('Initializing zeta')
%

%interpolate values to new grid
zeta=nc_clm{'zeta'}(tidx,:,:);

display('Initializing u, v, ubar, and vbar')

ubar=nc_clm{'ubar'}(tidx,:,:);

vbar=nc_clm{'vbar'}(tidx,:,:);

u=nc_clm{'u'}(tidx,:,:,:);

v=nc_clm{'v'}(tidx,:,:,:);

display('Initializing temp and salt')

salt=nc_clm{'salt'}(tidx,:,:,:);

temp=nc_clm{'temp'}(tidx,:,:,:);

NAT=2;                            % Number of active tracers. Usually 2 (temp + salt). Same 

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
  eval(['mud_',count,'(1:length(init_time),1:N,1:eta_rho,1:xi_rho) = 0;'])               %mud conc in water column
end
%
% sand.
%
for isand=1:NNS
  count=['0',num2str(isand)];
  count=count(end-1:end);
  eval(['sand_',count,'(1:length(init_time),1:N,1:eta_rho,1:xi_rho) = 0;'])               %mud conc in water column
end

%10) Provide initial sediment properties in bed.
%
% bed properties
  display('Initializing sediment bed.')
%
for k=1:Nbed
  for time=1:length(init_time)
    bed_thickness(time,k,1:eta_rho,1:xi_rho) = 0.1;
    bed_age(time,k,1:eta_rho,1:xi_rho)       = init_time(1);
    bed_porosity(time,k,1:eta_rho,1:xi_rho)  = 0.9;
    bed_biodiff(time,k,1:eta_rho,1:xi_rho)   = 0.0;
  end
end
%
% for mud
%
for idsed=1:NCS
  count=['0',num2str(idsed)];
  count=count(end-1:end);
  for k=1:Nbed
    for time=1:length(init_time)
      eval(['mudfrac_',count,'(time,k,1:eta_rho,1:xi_rho) = 1/NST;'])      %fraction of each sed class in each bed cell
      eval(['mudmass_',count,'(time,k,1:eta_rho,1:xi_rho) = squeeze(bed_thickness(time,k,1:eta_rho,1:xi_rho)).*Srho(idsed).*(1.0-squeeze(bed_porosity(time,k,1:eta_rho,1:xi_rho))).*squeeze(mudfrac_',count,'(time,k,1:eta_rho,1:xi_rho));'])          %mass of each sed class in each bed cell
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
    for time=1:length(init_time)
      eval(['sandfrac_',count,'(time,k,1:eta_rho,1:xi_rho) = 1/NST;'])      %fraction of each sed class in each bed cell
      eval(['sandmass_',count,'(time,k,1:eta_rho,1:xi_rho) = squeeze(bed_thickness(time,k,1:eta_rho,1:xi_rho)).*Srho(isand).*(1.0-squeeze(bed_porosity(time,k,1:eta_rho,1:xi_rho))).*squeeze(sandfrac_',count,'(time,k,1:eta_rho,1:xi_rho));'])          %mass of each sed class in each bed cell
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
      eval(['cff1=cff1*mud_Sd50(ised)^squeeze(mudfrac_',count,'(1,1,1:eta_rho,1:xi_rho));'])
      eval(['cff2=cff2*mud_Srho(ised)^squeeze(mudfrac_',count,'(1,1,1:eta_rho,1:xi_rho));'])
      eval(['cff3=cff3*mud_Wsed(ised)^squeeze(mudfrac_',count,'(1,1,1:eta_rho,1:xi_rho));'])
      eval(['cff4=cff4*mud_tau_ce(ised)^squeeze(mudfrac_',count,'(1,1,1:eta_rho,1:xi_rho));'])
    end
    for ised=1:NNS
      count=['0',num2str(ised)];
      count=count(end-1:end);
      eval(['cff1=cff1*sand_Sd50(ised).^squeeze(sandfrac_',count,'(1,1,1:eta_rho,1:xi_rho));'])
      eval(['cff2=cff2*sand_Srho(ised).^squeeze(sandfrac_',count,'(1,1,1:eta_rho,1:xi_rho));'])
      eval(['cff3=cff3*sand_Wsed(ised).^squeeze(sandfrac_',count,'(1,1,1:eta_rho,1:xi_rho));'])
      eval(['cff4=cff4*sand_tau_ce(ised).^squeeze(sandfrac_',count,'(1,1,1:eta_rho,1:xi_rho));'])
    end
    grain_diameter(time,1:eta_rho,1:xi_rho)=cff1;
    grain_density(time,1:eta_rho,1:xi_rho)=cff2;
    settling_vel(time,1:eta_rho,1:xi_rho)=cff3;
    erosion_stress(time,1:eta_rho,1:xi_rho)=cff4;
    ripple_length(time,1:eta_rho,1:xi_rho)=0.10;
    ripple_height(time,1:eta_rho,1:xi_rho)=0.01;
    dmix_offset(time,1:eta_rho,1:xi_rho)=0.0;
    dmix_slope(time,1:eta_rho,1:xi_rho)=0.0;
    dmix_time(time,1:eta_rho,1:xi_rho)=0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now create the netcdf file first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

create_roms_netcdf_init(init_file,gn,t_clim,Nbed,NNS,NCS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now write the data from the arrays to the netcdf file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ## Filling Variables in netcdf file with data...')

nc_init=netcdf(init_file,'w');
nc_init{'theta_s'}(:) = theta_s;
nc_init{'theta_b'}(:) = theta_b;
nc_init{'Tcline'}(:)  = Tcline;
nc_init{'Cs_r'}(:) = Cs_r;
nc_init{'Cs_w'}(:) = Cs_w;
nc_init{'sc_w'}(:) = sc_w;
nc_init{'sc_r'}(:) = sc_r;
nc_init{'hc'}(:) = hc;

nc_init{'ocean_time'}(:) = init_time;
nc_init{'temp'}(:) = temp;
nc_init{'salt'}(:) = salt;
nc_init{'u'}(:)    = u;
nc_init{'ubar'}(:) = ubar;
nc_init{'v'}(:)    = v;
nc_init{'vbar'}(:) = vbar;
nc_init{'zeta'}(:) = zeta;

for mm=1:NCS
  count=['00',num2str(mm)];
  count=count(end-1:end);
  eval(['nc_init{''mud_',count,'''}(:)     = mud_',count,';'])           %sed conc in water column
  eval(['nc_init{''mudfrac_',count,'''}(:) = mudfrac_',count,';'])       %fraction of each sed class in each bed cell
  eval(['nc_init{''mudmass_',count,'''}(:) = mudmass_',count,';'])       %mass of each sed class in each bed cell
end
for mm=1:NNS
  count=['00',num2str(mm)];
  count=count(end-1:end);
  eval(['nc_init{''sand_',count,'''}(:)     = sand_',count,';'])           %sed conc in water column
  eval(['nc_init{''sandfrac_',count,'''}(:) = sandfrac_',count,';'])       %fraction of each sed class in each bed cell
  eval(['nc_init{''sandmass_',count,'''}(:) = sandmass_',count,';'])       %mass of each sed class in each bed cell
end

nc_init{'bed_thickness'}(:) = bed_thickness;
nc_init{'bed_age'}(:)       = bed_age;
nc_init{'bed_porosity'}(:)  = bed_porosity;
nc_init{'bed_biodiff'}(:)   = bed_biodiff;

nc_init{'grain_diameter'}(:) = grain_diameter;
nc_init{'grain_density'}(:)  = grain_density;
nc_init{'settling_vel'}(:)   = settling_vel;
nc_init{'erosion_stress'}(:) = erosion_stress;
nc_init{'ripple_height'}(:) = ripple_height;
nc_init{'ripple_length'}(:) = ripple_length;
nc_init{'dmix_offset'}(:) = dmix_offset;
nc_init{'dmix_slope'}(:) = dmix_slope;
nc_init{'dmix_time'}(:) = dmix_time;

%close file
close(nc_init)
close(nc_clm)



