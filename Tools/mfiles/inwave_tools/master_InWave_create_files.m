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

clear all
close all

%1) SET THE CASE TO = 1.

INWAVE_SHOREFACE=0;
INLET_TEST=0;
LIP=1;
MY_APP=0;

%2) Provide the name of the mfile containing configuration parameters of
% the run;

if (INWAVE_SHOREFACE)
    inwave_gen_file='InWave_shoreface_param';
elseif (INLET_TEST)
    inwave_gen_file='InWave_inlet_test_param';
elseif (LIP)
    inwave_gen_file='InWave_lip_param';
elseif (MY_APP)
    inwave_gen_file='InWave_myapp_param';
end

%%%%%%%%%%%%%%%%%%%%% END OF USER SECTION %%%%%%%%%%%%%%%%%%%%%%%%%

%3) Evaluate the configuration parameters file

eval([inwave_gen_file])

%4) CREATE InWave GRID FILE
if (make_InWave_grd)
  create_InWave_grd(x,y,dx,dy,depth,roms_angle,mask_rho,f,spherical,grd_file)
end

%5) CREATE InWave INI FILE
if (make_InWave_ini)
  create_InWave_ini(Lp,Mp,Nbins,Bindirs_centers,Ac,Cx,Cy,Ct,TA,ini_file)
end

%6) CREATE InWave BRY FILE
if (make_InWave_bry)

  if (bin_error==1)

    disp([' ERROR WHEN CREATING INWAVE BOUNDARY FILE:'])
    disp([' You need to change the direction of the bins containing the energy'])
    disp([' At least one of them does not coincide with the considered central bin angles'])
    disp([' Select one of the following:'])
    disp([Bindirs])

  else

    create_InWave_bry(Lp, Mp, Nbins_bnd, dir_bnd, obc, ...
                      Ac_north,Ac_east,Ac_south,Ac_west,TA,time, bry_file)
  end

  % END OF INWAVE FILE GENERATION

  disp(['INWAVE FILES CREATED'])

end

