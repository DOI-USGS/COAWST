% create_InWave_files.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates all the nc files needed to run the InWave model: 
% InWave_grd.nc
% InWave_ini.nc
% InWave_bnd.nc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% USER SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     The user has to define which application is going to run
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1) SET THE CASE TO = 1.

STEP=1;
SOLITARY=0;
BEACH=0;
SHOAL=0;
JET=0;
DIR_TRANSPORT=0;
MY_APP=1;


%2) PRovide the name of the LOAD THE CORRESPONDING APPLICATION CONFIGURATION FILE

if (STEP) 
  inwave_gen_file='InWave_step_param';
end
elseif (SOLITARY)
    InWave_solitary_param
elseif (BEACH)
    InWave_beach_param
elseif (SHOAL)
    InWave_shoal_param
elseif (JET)
    InWave_jet_param
elseif (DIR_TRANSPORT)
    InWave_dir_transport_param
elseif (MY_APP)
%    InWave_myapp_param
    InWave_test_willap_myapp_param
end

%%%%%%%%%%%%%%%%%%%%% END OF USER SECTION %%%%%%%%%%%%%%%%%%%%%%%%%

%2.5
eval([inwave_gen_file])

%3) CREATE InWave GRID FILE

create_inwave_grid(x,y,dx,dy,depth,roms_angle,mask_rho,f,grd_file)

%4) CREATE InWave INI FILE

create_inwave_ini(Lm,Mm,Nangle,Ac,Cx,Cy,Ct,TA,ini_file)

%5) CREATE InWave BND FILE

create_inwave_bnd(Lm, Mm, Nangle_bnd, dir_bnd, obc, ...
    Ac_west, Ac_east, Ac_south, Ac_north,...
    TA_west, TA_east, TA_south, TA_north,...
    time, bnd_file)

%6) END OF FILE GENERATION

disp(['INWAVE FILES CREATED'])

