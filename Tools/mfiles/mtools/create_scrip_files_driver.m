% script create_scrip_input_files
%
% This m file is part of a series of routines driven by "create_scrip_weights_master.m"
% Calls scrip for each combination of grids to create the weights files.
%
% jcwarner 14July2014
%

%
%%%%%%%%%%%%%%%%%%%%%   OCN2WAV    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for mw=1:Ngrids_swan
  for mo=1:Ngrids_roms
%
% create the source file
%
    gf=roms_grids{mo};
    mask=mask_rho_o{mo};
    source_file='src.nc';
    create_scrip_roms_infile(gf,mask,source_file);
%
% create the destination file
%
    bf=bath_file{mw};
    gf=swan_coord{mw};
    nx=Numx{mw};
    ny=Numy{mw};
    eval(['mask=O',num2str(mo),'src_W',num2str(mw),'dst;'])
    cart=cartesian{mw};
    dest_file='dst.nc';
    create_scrip_swan_infile(bf,gf,nx,ny,cart,mask,dest_file);
%
% create the scrip weights file
    interp_file=['ocn',num2str(mo),'_to_','wav',num2str(mw),'_weights.nc']
    % create scrip in file
    fid=fopen('scrip_in','w');
    fprintf(fid,'&remap_inputs\n');
    fprintf(fid,'num_maps = 1\n');
    fprintf(fid,['grid1_file = ','''',source_file,'''','\n']);
    fprintf(fid,['grid2_file = ','''',dest_file,'''','\n']);
    fprintf(fid,['interp_file1 = ','''',interp_file,'''','\n']);
    fprintf(fid,'map1_name = ''OCN2WAV distwgt Mapping''\n');
    fprintf(fid,'map_method = ''conservative''\n');
    fprintf(fid,'normalize_opt = ''fracarea''\n');
    fprintf(fid,'output_opt = ''scrip''\n');
    fprintf(fid,'restrict_type = ''latitude''\n');
    fprintf(fid,'num_srch_bins = 90 \n');
    fprintf(fid,'luse_grid1_area = .false.\n');
    fprintf(fid,'luse_grid2_area = .false.\n');
    fprintf(fid,'/\n');
    fclose(fid);
    % call scrip
    check_scrip_src_dst
    if (use_scrip_exe)
      eval(['!"',scrip_exe,'"'])
    else
     scrip
    end
  end
end
%
%%%%%%%%%%%%%%%%%%%%%   WAV2OCN    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for mo=1:Ngrids_roms
  for mw=1:Ngrids_swan
%
% create the src file
%
    bf=bath_file{mw};
    gf=swan_coord{mw};
    nx=Numx{mw};
    ny=Numy{mw};
    mask=mask_rho_w{mw};
    cart=cartesian{mw};
    source_file='src.nc';
    create_scrip_swan_infile(bf,gf,nx,ny,cart,mask,source_file);
%
% create the destination file
%
    gf=roms_grids{mo};
    eval(['mask=W',num2str(mw),'src_O',num2str(mo),'dst;'])
    dest_file='dst.nc';
    create_scrip_roms_infile(gf,mask,dest_file);
%
% create the scrip weights file
    interp_file=['wav',num2str(mw),'_to_','ocn',num2str(mo),'_weights.nc']
    % create scrip in file
    fid=fopen('scrip_in','w');
    fprintf(fid,'&remap_inputs\n');
    fprintf(fid,'num_maps = 1\n');
    fprintf(fid,['grid1_file = ','''',source_file,'''','\n']);
    fprintf(fid,['grid2_file = ','''',dest_file,'''','\n']);
    fprintf(fid,['interp_file1 = ','''',interp_file,'''','\n']);
    fprintf(fid,'map1_name = ''WAV2OCN distwgt Mapping''\n');
    fprintf(fid,'map_method = ''conservative''\n');
    fprintf(fid,'normalize_opt = ''fracarea''\n');
    fprintf(fid,'output_opt = ''scrip''\n');
    fprintf(fid,'restrict_type = ''latitude''\n');
    fprintf(fid,'num_srch_bins = 90 \n');
    fprintf(fid,'luse_grid1_area = .false.\n');
    fprintf(fid,'luse_grid2_area = .false.\n');
    fprintf(fid,'/\n');
    fclose(fid);
    % call scrip
    check_scrip_src_dst
    if (use_scrip_exe)
      eval(['!"',scrip_exe,'"'])
    else
     scrip
    end
  end
end
%
%%%%%%%%%%%%%%%%%%%%%   ATM2OCN    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for mo=1:Ngrids_roms
  for ma=1:Ngrids_wrf
%
% create the src file
%
    gf=wrf_grids{ma};
    mask=mask_rho_a{ma};
    source_file='src.nc';
    create_scrip_wrf_infile(gf,mask,source_file);
%
% create the destination file
%
    gf=roms_grids{mo};
    eval(['mask=A',num2str(ma),'src_O',num2str(mo),'dst;'])
    dest_file='dst.nc';
    create_scrip_roms_infile(gf,mask,dest_file);
%
% create the scrip weights file
    interp_file=['atm',num2str(ma),'_to_','ocn',num2str(mo),'_weights.nc']
    % create scrip in file
    fid=fopen('scrip_in','w');
    fprintf(fid,'&remap_inputs\n');
    fprintf(fid,'num_maps = 1\n');
    fprintf(fid,['grid1_file = ','''',source_file,'''','\n']);
    fprintf(fid,['grid2_file = ','''',dest_file,'''','\n']);
    fprintf(fid,['interp_file1 = ','''',interp_file,'''','\n']);
    fprintf(fid,'map1_name = ''ATM2OCN distwgt Mapping''\n');
    fprintf(fid,'map_method = ''conservative''\n');
    fprintf(fid,'normalize_opt = ''fracarea''\n');
    fprintf(fid,'output_opt = ''scrip''\n');
    fprintf(fid,'restrict_type = ''latitude''\n');
    fprintf(fid,'num_srch_bins = 90 \n');
    fprintf(fid,'luse_grid1_area = .false.\n');
    fprintf(fid,'luse_grid2_area = .false.\n');
    fprintf(fid,'/\n');
    fclose(fid);
    % call scrip
    check_scrip_src_dst
    if (use_scrip_exe)
      eval(['!"',scrip_exe,'"'])
    else
     scrip
    end
  end
end
%
%%%%%%%%%%%%%%%%%%%%%   OCN2ATM    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for ma=1:Ngrids_wrf
  for mo=1:Ngrids_roms
%
% create the source file
%
    gf=roms_grids{mo};
    mask=mask_rho_o{mo};
    source_file='src.nc';
    create_scrip_roms_infile(gf,mask,source_file);
%
% create the destination file
%
    gf=wrf_grids{ma};
    eval(['mask=O',num2str(mo),'src_A',num2str(ma),'dst;'])
    dest_file='dst.nc';
    create_scrip_wrf_infile(gf,mask,dest_file);
%
% create the scrip weights file
    interp_file=['ocn',num2str(mo),'_to_','atm',num2str(ma),'_weights.nc']
    % create scrip in file
    fid=fopen('scrip_in','w');
    fprintf(fid,'&remap_inputs\n');
    fprintf(fid,'num_maps = 1\n');
    fprintf(fid,['grid1_file = ','''',source_file,'''','\n']);
    fprintf(fid,['grid2_file = ','''',dest_file,'''','\n']);
    fprintf(fid,['interp_file1 = ','''',interp_file,'''','\n']);
    fprintf(fid,'map1_name = ''OCN2ATM distwgt Mapping''\n');
    fprintf(fid,'map_method = ''conservative''\n');
    fprintf(fid,'normalize_opt = ''fracarea''\n');
    fprintf(fid,'output_opt = ''scrip''\n');
    fprintf(fid,'restrict_type = ''latitude''\n');
    fprintf(fid,'num_srch_bins = 90 \n');
    fprintf(fid,'luse_grid1_area = .false.\n');
    fprintf(fid,'luse_grid2_area = .false.\n');
    fprintf(fid,'/\n');
    fclose(fid);
    % call scrip
    check_scrip_src_dst
    if (use_scrip_exe)
      eval(['!"',scrip_exe,'"'])
    else
     scrip
    end
  end
end
%
%%%%%%%%%%%%%%%%%%%%%   ATM2WAV    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for mw=1:Ngrids_swan
  for ma=1:Ngrids_wrf
%
% create the source file
%
    gf=wrf_grids{ma};
    mask=mask_rho_a{ma};
    source_file='src.nc';
    create_scrip_wrf_infile(gf,mask,source_file);
%
% create the destination file
%
    bf=bath_file{mw};
    gf=swan_coord{mw};
    nx=Numx{mw};
    ny=Numy{mw};
    eval(['mask=A',num2str(ma),'src_W',num2str(mw),'dst;'])
    cart=cartesian{mw};
    dest_file='dst.nc';
    create_scrip_swan_infile(bf,gf,nx,ny,cart,mask,dest_file);
%
% create the scrip weights file
    interp_file=['atm',num2str(ma),'_to_','wav',num2str(mw),'_weights.nc']
    % create scrip in file
    fid=fopen('scrip_in','w');
    fprintf(fid,'&remap_inputs\n');
    fprintf(fid,'num_maps = 1\n');
    fprintf(fid,['grid1_file = ','''',source_file,'''','\n']);
    fprintf(fid,['grid2_file = ','''',dest_file,'''','\n']);
    fprintf(fid,['interp_file1 = ','''',interp_file,'''','\n']);
    fprintf(fid,'map1_name = ''ATM2WAV distwgt Mapping''\n');
    fprintf(fid,'map_method = ''conservative''\n');
    fprintf(fid,'normalize_opt = ''fracarea''\n');
    fprintf(fid,'output_opt = ''scrip''\n');
    fprintf(fid,'restrict_type = ''latitude''\n');
    fprintf(fid,'num_srch_bins = 90 \n');
    fprintf(fid,'luse_grid1_area = .false.\n');
    fprintf(fid,'luse_grid2_area = .false.\n');
    fprintf(fid,'/\n');
    fclose(fid);
    % call scrip
    check_scrip_src_dst
    if (use_scrip_exe)
      eval(['!"',scrip_exe,'"'])
    else
     scrip
    end
  end
end
%
%%%%%%%%%%%%%%%%%%%%%   WAV2ATM    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for ma=1:Ngrids_wrf
  for mw=1:Ngrids_swan
%
% create the src file
%
    bf=bath_file{mw};
    gf=swan_coord{mw};
    nx=Numx{mw};
    ny=Numy{mw};
    mask=mask_rho_w{mw};
    cart=cartesian{mw};
    source_file='src.nc';
    create_scrip_swan_infile(bf,gf,nx,ny,cart,mask,source_file);
%
% create the destination file
%
    gf=wrf_grids{ma};
    eval(['mask=W',num2str(mw),'src_A',num2str(ma),'dst;'])
    dest_file='dst.nc';
    create_scrip_wrf_infile(gf,mask,dest_file);
%
% create the scrip weights file
    interp_file=['wav',num2str(mw),'_to_','atm',num2str(ma),'_weights.nc']
    % create scrip in file
    fid=fopen('scrip_in','w');
    fprintf(fid,'&remap_inputs\n');
    fprintf(fid,'num_maps = 1\n');
    fprintf(fid,['grid1_file = ','''',source_file,'''','\n']);
    fprintf(fid,['grid2_file = ','''',dest_file,'''','\n']);
    fprintf(fid,['interp_file1 = ','''',interp_file,'''','\n']);
    fprintf(fid,'map1_name = ''WAV2ATM distwgt Mapping''\n');
    fprintf(fid,'map_method = ''conservative''\n');
    fprintf(fid,'normalize_opt = ''fracarea''\n');
    fprintf(fid,'output_opt = ''scrip''\n');
    fprintf(fid,'restrict_type = ''latitude''\n');
    fprintf(fid,'num_srch_bins = 90 \n');
    fprintf(fid,'luse_grid1_area = .false.\n');
    fprintf(fid,'luse_grid2_area = .false.\n');
    fprintf(fid,'/\n');
    fclose(fid);
    % call scrip
    check_scrip_src_dst
    if (use_scrip_exe)
      eval(['!"',scrip_exe,'"'])
    else
     scrip
    end
  end
end
%

