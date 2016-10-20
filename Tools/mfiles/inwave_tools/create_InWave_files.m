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

STEP=0;
SOLITON=0;
BEACH=0;
SHOAL=0;
JET=0;
BICRO=0;
DIR_TRANSPORT=0;
DUCK85=0;
EDGE=0;
INWAVE_SHOREFACE=0;
INWAVE_SHORE=0;
INLET_TEST=0;
INWAVE_SANDY=1;
MY_APP=0;

%2) Provide the name of the mfile containing configuration parameters of
% the run;

if (STEP)
    inwave_gen_file='InWave_step_param';
elseif (SOLITON)
    inwave_gen_file='InWave_soliton_param';
elseif (BEACH)
    inwave_gen_file='InWave_beach_param';
elseif (SHOAL)
    inwave_gen_file='InWave_shoal_param';
elseif (BICRO)
    inwave_gen_file='InWave_bicro_param';
elseif (JET)
    inwave_gen_file='InWave_jet_param';
elseif (DIR_TRANSPORT)
    inwave_gen_file='InWave_dir_transport_param';
elseif (DUCK85)
    inwave_gen_file='InWave_duck85_param';
elseif (EDGE)
    inwave_gen_file='InWave_edge_param';
elseif (INWAVE_SHOREFACE)
    inwave_gen_file='InWave_shoreface_param';
elseif (INWAVE_SHORE)
    inwave_gen_file='InWave_shore_param';
elseif (INWAVE_SANDY)
    inwave_gen_file='InWave_Sandy_param_strip';
elseif (INLET_TEST)
    inwave_gen_file='InWave_inlet_test_param';
elseif (MY_APP)
    inwave_gen_file='InWave_myapp_param';
end

%%%%%%%%%%%%%%%%%%%%% END OF USER SECTION %%%%%%%%%%%%%%%%%%%%%%%%%

%2.5 Includes the configuration parameters

eval([inwave_gen_file])

%3) CREATE InWave GRID FILE
if (make_InWave_grd)
  create_InWave_grid(x,y,dx,dy,depth,roms_angle,mask_rho,f,grd_file)
end

%4) CREATE InWave INI FILE
if (make_InWave_ini)
  create_InWave_ini(LP,MP,Nbins,Bindirs,Bindirs_c,pd,Ac,Cx,Cy,Ct,TA,ini_file)
end

%5) CREATE InWave BND FILE
if (make_InWave_bnd)

  if (bin_error==1)

    disp([' ERROR WHEN CREATING INWAVE BOUNDARY FILE:'])
    disp([' You need to change the direction of the bins containing the energy'])
    disp([' At least one of them does not coincide with the considered central bin angles'])
    disp([' Select one of the following:'])
    disp([Bindirs])

  else

%     if obc(1)==1
%         create_inwave_bnd(Lm, Mm, Nbins_bnd, dir_bnd, obc, ...
%             Ac_north,TA,time, bnd_file)
%     elseif obc(2)==1
%         create_inwave_bnd(Lm, Mm, Nbins_bnd, dir_bnd, obc, ...
%             Ac_east,TA,time, bnd_file)
%     elseif obc(3)==1
%         create_inwave_bnd(Lm, Mm, Nbins_bnd, dir_bnd, obc, ...
%             Ac_south,TA,time, bnd_file)
%     elseif obc(4)==1
%         create_inwave_bnd(Lm, Mm, Nbins_bnd, dir_bnd, obc, ...
%             Ac_west,TA,time, bnd_file)
%     end

  create_InWave_bnd(LP, MP, Nbins, Bindirs_c, obc, ...
    Ac_north,Ac_east,Ac_south,Ac_west,TA,time, bnd_file)

    %6) END OF INWAVE FILE GENERATION

    if (DUCK85)
        for ii=1:3
            if ii==1
                zeta=zeta_list;
                ubar=ubar_list;
                bryname=strcat(bry_file(1:end-3),'_list.nc');
            elseif ii==2
                zeta=zeta_on;
                ubar=ubar_on;
                bryname=strcat(bry_file(1:end-3),'_on.nc');
            elseif ii==3
                zeta=zeta_don';
                ubar=ubar_don';
                bryname=strcat(bry_file(1:end-3),'_don.nc');
            end


            grdname=grd_file;
            create_bryfile(bryname,grdname,time,zeta,ubar)
        end
    end
  end

    %7) END OF INWAVE FILE GENERATION

    disp(['INWAVE FILES CREATED'])


end

