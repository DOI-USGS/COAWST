function create_roms_netcdf_init_mw(init_file,gn,Nbed,NNS,NCS,NVEG)

%create init file
nc_init=netcdf.create(init_file,'clobber');
 
%% Global attributes:

disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_init,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by updatclim on ' datestr(now)]);
netcdf.putAtt(nc_init,netcdf.getConstant('NC_GLOBAL'),'type', 'initial forcing file from http://hycom.coaps.fsu.edu:8080/thredds/dodsC/glb_analysis');


%% Dimensions:

disp(' ## Defining Dimensions...')
 
%get some grid info
  [LP,MP]=size(gn.lon_rho);
  L=LP-1;
  Lm=L-1;
  M=MP-1;
  Mm=M-1;
  L  = Lm+1;
  M  = Mm+1;
  xpsi  = L;
  xrho  = LP;
  xu    = L;
  xv    = LP;
  epsi = M;
  erho = MP;
  eu   = MP;
  ev   = M;
  N       = gn.N;
  
psidimID = netcdf.defDim(nc_init,'xpsi',L);
xrhodimID = netcdf.defDim(nc_init,'xrho',LP);
xudimID = netcdf.defDim(nc_init,'xu',L);
xvdimID = netcdf.defDim(nc_init,'xv',LP);

epsidimID = netcdf.defDim(nc_init,'epsi',M);
erhodimID = netcdf.defDim(nc_init,'erho',MP);
eudimID = netcdf.defDim(nc_init,'eu',MP);
evdimID = netcdf.defDim(nc_init,'ev',M);

s_rhodimID = netcdf.defDim(nc_init,'sc_r',N);
s_wdimID = netcdf.defDim(nc_init,'sc_w',N+1);
NbeddimID = netcdf.defDim(nc_init,'Nbed',Nbed);
timedimID = netcdf.defDim(nc_init,'time',1);
NvegdimID = netcdf.defDim(nc_init,'Nveg',NVEG);

%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')

sphericalID = netcdf.defVar(nc_init,'spherical','short',timedimID);
netcdf.putAtt(nc_init,sphericalID,'long_name','grid type logical switch');
netcdf.putAtt(nc_init,sphericalID,'flag_meanings','spherical Cartesian');
netcdf.putAtt(nc_init,sphericalID,'flag_values','1, 0');

VtransformID = netcdf.defVar(nc_init,'Vtransform','long',timedimID);
netcdf.putAtt(nc_init,VtransformID,'long_name','vertical terrain-following transformation equation');

VstretchingID = netcdf.defVar(nc_init,'Vstretching','long',timedimID);
netcdf.putAtt(nc_init,VstretchingID,'long_name','vertical terrain-following stretching function');
 
theta_bID = netcdf.defVar(nc_init,'theta_b','double',timedimID);
netcdf.putAtt(nc_init,theta_bID,'long_name','S-coordinate bottom control parameter');
netcdf.putAtt(nc_init,theta_bID,'units','1');

theta_sID = netcdf.defVar(nc_init,'theta_s','double',timedimID);
netcdf.putAtt(nc_init,theta_sID,'long_name','S-coordinate surface control parameter');
netcdf.putAtt(nc_init,theta_sID,'units','1');

tcline_ID = netcdf.defVar(nc_init,'Tcline','double',timedimID);
netcdf.putAtt(nc_init,tcline_ID,'long_name','S-coordinate surface/bottom layer width');
netcdf.putAtt(nc_init,tcline_ID,'units','meter');

hc_ID = netcdf.defVar(nc_init,'hc','double',timedimID);
netcdf.putAtt(nc_init,hc_ID,'long_name','S-coordinate parameter, critical depth');
netcdf.putAtt(nc_init,hc_ID,'units','meter');

Cs_rID = netcdf.defVar(nc_init,'Cs_r','double',s_rhodimID);
netcdf.putAtt(nc_init,Cs_rID,'long_name','S-coordinate stretching curves at RHO-points');
netcdf.putAtt(nc_init,Cs_rID,'units','1');
netcdf.putAtt(nc_init,Cs_rID,'valid_min',-1);
netcdf.putAtt(nc_init,Cs_rID,'valid_max',0);
netcdf.putAtt(nc_init,Cs_rID,'field','Cs_r, scalar');

Cs_wID = netcdf.defVar(nc_init,'Cs_w','double',s_wdimID);
netcdf.putAtt(nc_init,Cs_wID,'long_name','S-coordinate stretching curves at W-points');
netcdf.putAtt(nc_init,Cs_wID,'units','1');
netcdf.putAtt(nc_init,Cs_wID,'valid_min',-1);
netcdf.putAtt(nc_init,Cs_wID,'valid_max',0);
netcdf.putAtt(nc_init,Cs_wID,'field','Cs_w, scalar');

sc_rID = netcdf.defVar(nc_init,'sc_r','double',s_rhodimID);
netcdf.putAtt(nc_init,sc_rID,'long_name','S-coordinate at RHO-points');
netcdf.putAtt(nc_init,sc_rID,'units','1');
netcdf.putAtt(nc_init,sc_rID,'valid_min',-1);
netcdf.putAtt(nc_init,sc_rID,'valid_max',0);
netcdf.putAtt(nc_init,sc_rID,'field','sc_r, scalar');

sc_wID = netcdf.defVar(nc_init,'sc_w','double',s_wdimID);
netcdf.putAtt(nc_init,sc_wID,'long_name','S-coordinate at W-points');
netcdf.putAtt(nc_init,sc_wID,'units','1');
netcdf.putAtt(nc_init,sc_wID,'valid_min',-1);
netcdf.putAtt(nc_init,sc_wID,'valid_max',0);
netcdf.putAtt(nc_init,sc_wID,'field','sc_w, scalar');

ocean_timeID = netcdf.defVar(nc_init,'ocean_time','double',timedimID);
netcdf.putAtt(nc_init,ocean_timeID,'long_name','time since initialization');
netcdf.putAtt(nc_init,ocean_timeID,'units','days');
netcdf.putAtt(nc_init,ocean_timeID,'field','ocean_time, scalar, series');

saltID = netcdf.defVar(nc_init,'salt','float',[xrhodimID erhodimID s_rhodimID timedimID]);
netcdf.putAtt(nc_init,saltID,'long_name','salinity');
netcdf.putAtt(nc_init,saltID,'units','PSU');
netcdf.putAtt(nc_init,saltID,'field','salinity, scalar, series');

tempID = netcdf.defVar(nc_init,'temp','float',[xrhodimID erhodimID s_rhodimID timedimID]);
netcdf.putAtt(nc_init,tempID,'long_name','temperature');
netcdf.putAtt(nc_init,tempID,'units','C');
netcdf.putAtt(nc_init,tempID,'field','temperature, scalar, series');

uID = netcdf.defVar(nc_init,'u','float',[xudimID eudimID s_rhodimID timedimID]);
netcdf.putAtt(nc_init,uID,'long_name','u-momentum component');
netcdf.putAtt(nc_init,uID,'units','meter second-1');
netcdf.putAtt(nc_init,uID,'field','u-velocity, scalar, series');

ubarID = netcdf.defVar(nc_init,'ubar','float',[xudimID eudimID timedimID]);
netcdf.putAtt(nc_init,ubarID,'long_name','vertically integrated u-momentum component');
netcdf.putAtt(nc_init,ubarID,'units','meter second-1');
netcdf.putAtt(nc_init,ubarID,'field','ubar-velocity, scalar, series');

vID = netcdf.defVar(nc_init,'v','float',[xvdimID evdimID s_rhodimID timedimID]);
netcdf.putAtt(nc_init,vID,'long_name','v-momentum component');
netcdf.putAtt(nc_init,vID,'units','meter second-1');
netcdf.putAtt(nc_init,vID,'field','v-velocity, scalar, series');

vbarID = netcdf.defVar(nc_init,'vbar','float',[xvdimID evdimID timedimID]);
netcdf.putAtt(nc_init,vbarID,'long_name','vertically integrated v-momentum component');
netcdf.putAtt(nc_init,vbarID,'units','meter second-1');
netcdf.putAtt(nc_init,vbarID,'field','vbar-velocity, scalar, series');
 
zetaID = netcdf.defVar(nc_init,'zeta','float',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,zetaID,'long_name','free-surface');
netcdf.putAtt(nc_init,zetaID,'units','meter');
netcdf.putAtt(nc_init,zetaID,'field','free-surface, scalar, series');
 
for mm=1:NCS
    count=['00',num2str(mm)];
    count=count(end-1:end);

    eval(['mud_',count,'ID = netcdf.defVar(nc_init,''mud_',count,''',''double'',[xrhodimID erhodimID s_rhodimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,mud_',count,'ID,''long_name'',''suspended cohesive sediment, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,mud_',count,'ID,''units'',''kilogram meter-3'');'])
    eval(['netcdf.putAtt(nc_init,mud_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,mud_',count,'ID,''field'',''mud_',count,', scalar, series'');'])

    eval(['mudfrac_',count,'ID = netcdf.defVar(nc_init,''mudfrac_',count,''',''double'',[xrhodimID erhodimID NbeddimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,mudfrac_',count,'ID,''long_name'',''cohesive sediment fraction, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,mudfrac_',count,'ID,''units'',''nondimensional'');'])
    eval(['netcdf.putAtt(nc_init,mudfrac_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,mudfrac_',count,'ID,''field'',''mudfrac_',count,', scalar, series'');'])

     eval(['mudmass_',count,'ID = netcdf.defVar(nc_init,''mudmass_',count,''',''double'',[xrhodimID erhodimID NbeddimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,mudmass_',count,'ID,''long_name'',''cohesive sediment mass, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,mudmass_',count,'ID,''units'',''kilogram meter-2'');'])
    eval(['netcdf.putAtt(nc_init,mudmass_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,mudmass_',count,'ID,''field'',''mudmass_',count,', scalar, series'');'])

end
for mm=1:NNS
    count=['00',num2str(mm)];
    count=count(end-1:end);

    eval(['sand_',count,'ID = netcdf.defVar(nc_init,''sand_',count,''',''double'',[xrhodimID erhodimID s_rhodimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,sand_',count,'ID,''long_name'',''suspended noncohesive sediment, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,sand_',count,'ID,''units'',''kilogram meter-3'');'])
    eval(['netcdf.putAtt(nc_init,sand_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,sand_',count,'ID,''field'',''sand_',count,', scalar, series'');'])

    eval(['sandfrac_',count,'ID = netcdf.defVar(nc_init,''sandfrac_',count,''',''double'',[xrhodimID erhodimID NbeddimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,sandfrac_',count,'ID,''long_name'',''noncohesive sediment fraction, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,sandfrac_',count,'ID,''units'',''nondimensional'');'])
    eval(['netcdf.putAtt(nc_init,sandfrac_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,sandfrac_',count,'ID,''field'',''sandfrac_',count,', scalar, series'');'])

    eval(['sandmass_',count,'ID = netcdf.defVar(nc_init,''sandmass_',count,''',''double'',[xrhodimID erhodimID NbeddimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,sandmass_',count,'ID,''long_name'',''noncohesive sediment mass, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,sandmass_',count,'ID,''units'',''kilogram meter-2'');'])
    eval(['netcdf.putAtt(nc_init,sandmass_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,sandmass_',count,'ID,''field'',''sandmass_',count,', scalar, series'');'])

    eval(['bedload_Usand_',count,'ID = netcdf.defVar(nc_init,''bedload_Usand_',count,''',''double'',[xudimID eudimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,bedload_Usand_',count,'ID,''long_name'',''bed load flux of sand in U-direction, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,bedload_Usand_',count,'ID,''units'',''kilogram meter-1 s-1'');'])
    eval(['netcdf.putAtt(nc_init,bedload_Usand_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,bedload_Usand_',count,'ID,''field'',''bedload_Usand_',count,', scalar, series'');'])

    eval(['bedload_Vsand_',count,'ID = netcdf.defVar(nc_init,''bedload_Vsand_',count,''',''double'',[xvdimID evdimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,bedload_Vsand_',count,'ID,''long_name'',''bed load flux of sand in V-direction, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,bedload_Vsand_',count,'ID,''units'',''kilogram meter-1 s-1'');'])
    eval(['netcdf.putAtt(nc_init,bedload_Vsand_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,bedload_Vsand_',count,'ID,''field'',''bedload_Vsand_',count,', scalar, series'');'])
end

bed_thicknessID = netcdf.defVar(nc_init,'bed_thickness','double',[xrhodimID erhodimID NbeddimID timedimID]);
netcdf.putAtt(nc_init,bed_thicknessID,'long_name','sediment layer thickness');
netcdf.putAtt(nc_init,bed_thicknessID,'units','meter');
netcdf.putAtt(nc_init,bed_thicknessID,'time','ocean_time');
netcdf.putAtt(nc_init,bed_thicknessID,'field','bed thickness, scalar, series');

bed_ageID = netcdf.defVar(nc_init,'bed_age','double',[xrhodimID erhodimID NbeddimID timedimID]);
netcdf.putAtt(nc_init,bed_ageID,'long_name','sediment layer age');
netcdf.putAtt(nc_init,bed_ageID,'units','day');
netcdf.putAtt(nc_init,bed_ageID,'time','ocean_time');
netcdf.putAtt(nc_init,bed_ageID,'field','bed age, scalar, series');

bed_porosityID = netcdf.defVar(nc_init,'bed_porosity','double',[xrhodimID erhodimID NbeddimID timedimID]);
netcdf.putAtt(nc_init,bed_porosityID,'long_name','sediment layer porosity');
netcdf.putAtt(nc_init,bed_porosityID,'units','nondimensional');
netcdf.putAtt(nc_init,bed_porosityID,'time','ocean_time');
netcdf.putAtt(nc_init,bed_porosityID,'field','bed porosity, scalar, series');

bed_biodiffID = netcdf.defVar(nc_init,'bed_biodiff','double',[xrhodimID erhodimID NbeddimID timedimID]);
netcdf.putAtt(nc_init,bed_biodiffID,'long_name','biodiffusivity at bottom of each layer');
netcdf.putAtt(nc_init,bed_biodiffID,'units','meter2 second-1');
netcdf.putAtt(nc_init,bed_biodiffID,'time','ocean_time');
netcdf.putAtt(nc_init,bed_biodiffID,'field','bed biodiffusivity, scalar, series');

grain_diameterID = netcdf.defVar(nc_init,'grain_diameter','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,grain_diameterID,'long_name','sediment median grain diameter size');
netcdf.putAtt(nc_init,grain_diameterID,'units','meter');
netcdf.putAtt(nc_init,grain_diameterID,'time','ocean_time');
netcdf.putAtt(nc_init,grain_diameterID,'field','grain diameter, scalar, series');

grain_densityID = netcdf.defVar(nc_init,'grain_density','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,grain_densityID,'long_name','sediment median grain density');
netcdf.putAtt(nc_init,grain_densityID,'units','kilogram meter-3');
netcdf.putAtt(nc_init,grain_densityID,'time','ocean_time');
netcdf.putAtt(nc_init,grain_densityID,'field','grain density, scalar, series');

settling_velID = netcdf.defVar(nc_init,'settling_vel','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,settling_velID,'long_name','sediment median grain settling velocity');
netcdf.putAtt(nc_init,settling_velID,'units','meter second-1');
netcdf.putAtt(nc_init,settling_velID,'time','ocean_time');
netcdf.putAtt(nc_init,settling_velID,'field','settling vel, scalar, series');

erosion_stressID = netcdf.defVar(nc_init,'erosion_stress','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,erosion_stressID,'long_name','sediment median critical erosion stress');
netcdf.putAtt(nc_init,erosion_stressID,'units','meter2 second-2');
netcdf.putAtt(nc_init,erosion_stressID,'time','ocean_time');
netcdf.putAtt(nc_init,erosion_stressID,'field','erosion stress, scalar, series');

ripple_lengthID = netcdf.defVar(nc_init,'ripple_length','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,ripple_lengthID,'long_name','bottom ripple length');
netcdf.putAtt(nc_init,ripple_lengthID,'units','meter');
netcdf.putAtt(nc_init,ripple_lengthID,'time','ocean_time');
netcdf.putAtt(nc_init,ripple_lengthID,'field','ripple length, scalar, series');

ripple_heightID = netcdf.defVar(nc_init,'ripple_height','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,ripple_heightID,'long_name','bottom ripple height');
netcdf.putAtt(nc_init,ripple_heightID,'units','meter');
netcdf.putAtt(nc_init,ripple_heightID,'time','ocean_time');
netcdf.putAtt(nc_init,ripple_heightID,'field','ripple height, scalar, series');

dmix_offsetID = netcdf.defVar(nc_init,'dmix_offset','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,dmix_offsetID,'long_name','dmix erodibility profile offset');
netcdf.putAtt(nc_init,dmix_offsetID,'units','meter');
netcdf.putAtt(nc_init,dmix_offsetID,'time','ocean_time');
netcdf.putAtt(nc_init,dmix_offsetID,'field','dmix_offset, scalar, series');

dmix_slopeID = netcdf.defVar(nc_init,'dmix_slope','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,dmix_slopeID,'long_name','dmix erodibility profile slope');
netcdf.putAtt(nc_init,dmix_slopeID,'units','_');
netcdf.putAtt(nc_init,dmix_slopeID,'time','ocean_time');
netcdf.putAtt(nc_init,dmix_slopeID,'field','dmix_slope, scalar, series');

dmix_timeID = netcdf.defVar(nc_init,'dmix_time','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,dmix_timeID,'long_name','dmix erodibility profile time scale');
netcdf.putAtt(nc_init,dmix_timeID,'units','seconds');
netcdf.putAtt(nc_init,dmix_timeID,'time','ocean_time');
netcdf.putAtt(nc_init,dmix_timeID,'field','dmix_time, scalar, series');

vegID = netcdf.defVar(nc_init,'plant_height','double',[xrhodimID erhodimID NvegdimID timedimID]);
netcdf.putAtt(nc_init,vegID,'long_name','plant height');
netcdf.putAtt(nc_init,vegID,'units','meter');
netcdf.putAtt(nc_init,vegID,'time','ocean_time');
netcdf.putAtt(nc_init,vegID,'field','plant_height, scalar, series');

vegID = netcdf.defVar(nc_init,'plant_density','double',[xrhodimID erhodimID NvegdimID timedimID]);
netcdf.putAtt(nc_init,vegID,'long_name','plant density');
netcdf.putAtt(nc_init,vegID,'units','plant-meter2');
netcdf.putAtt(nc_init,vegID,'time','ocean_time');
netcdf.putAtt(nc_init,vegID,'field','plant_density, scalar, series');

vegID = netcdf.defVar(nc_init,'plant_diameter','double',[xrhodimID erhodimID NvegdimID timedimID]);
netcdf.putAtt(nc_init,vegID,'long_name','plant diameter');
netcdf.putAtt(nc_init,vegID,'units','meter');
netcdf.putAtt(nc_init,vegID,'time','ocean_time');
netcdf.putAtt(nc_init,vegID,'field','plant_diameter, scalar, series');

vegID = netcdf.defVar(nc_init,'plant_thickness','double',[xrhodimID erhodimID NvegdimID timedimID]);
netcdf.putAtt(nc_init,vegID,'long_name','plant thickness');
netcdf.putAtt(nc_init,vegID,'units','meter');
netcdf.putAtt(nc_init,vegID,'time','ocean_time');
netcdf.putAtt(nc_init,vegID,'field','plant_thickness, scalar, series');

vegID = netcdf.defVar(nc_init,'marsh_mask','double',[xrhodimID erhodimID NvegdimID timedimID]);
netcdf.putAtt(nc_init,vegID,'long_name','marsh mask');
netcdf.putAtt(nc_init,vegID,'units','nondimensional');
netcdf.putAtt(nc_init,vegID,'time','ocean_time');
netcdf.putAtt(nc_init,vegID,'field','marsh_mask, scalar, series');

netcdf.close(nc_init)



