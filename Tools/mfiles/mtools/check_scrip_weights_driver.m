% script create_scrip_input_files
%
% This m file is part of a series of routines driven by "create_scrip_weights_master.m"
% Checks the weigths in all the files to make sure that all 
% non-masked cells get data from a src grid.
%
% jcwarner 15Sept2014
%

%
%%%%%%%%%%%%%%%%%%%%%   OCN2WAV    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if (0)
  for mw=1:Ngrids_swan
    for mo=1:Ngrids_roms
      interp_file=['ocn',num2str(mo),'_to_','wav',num2str(mw),'_weights.nc']
      src_lon=lon_rho_o{mo};
      src_lat=lat_rho_o{mo};
      src_mask=mask_rho_o{mo};
      dst_lon=lon_rho_w{mw};
      dst_lat=lat_rho_w{mw};
      dst_mask=mask_rho_w{mw};
      add_new_weights(interp_file,src_lon,src_lat,src_mask,dst_lon,dst_lat,dst_mask)
    end
  end
end
%
%%%%%%%%%%%%%%%%%%%%%   WAV2OCN    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if (0)
  for mo=1:Ngrids_roms
    for mw=1:Ngrids_swan
      interp_file=['wav',num2str(mw),'_to_','ocn',num2str(mo),'_weights.nc']
      src_lon=lon_rho_w{mw};
      src_lat=lat_rho_w{mw};
      src_mask=mask_rho_w{mw};
      dst_lon=lon_rho_o{mo};
      dst_lat=lat_rho_o{mo};
      dst_mask=mask_rho_o{mo};
      add_new_weights(interp_file,src_lon,src_lat,src_mask,dst_lon,dst_lat,dst_mask)
    end
  end
end
%
%%%%%%%%%%%%%%%%%%%%%   ATM2OCN    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if (1)
  for mo=1:Ngrids_roms
    for ma=1:Ngrids_wrf
      interp_file=['atm',num2str(ma),'_to_','ocn',num2str(mo),'_weights.nc']
      src_lon=lon_rho_a{ma};
      src_lat=lat_rho_a{ma};
      src_mask=mask_rho_a{ma};
      dst_lon=lon_rho_o{mo};
      dst_lat=lat_rho_o{mo};
      dst_mask=mask_rho_o{mo};
      add_new_weights(interp_file,src_lon,src_lat,src_mask,dst_lon,dst_lat,dst_mask)
    end
  end
end
%
%%%%%%%%%%%%%%%%%%%%%   OCN2ATM    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if (1)
  for ma=1:Ngrids_wrf
    for mo=1:Ngrids_roms
      interp_file=['ocn',num2str(mo),'_to_','atm',num2str(ma),'_weights.nc']
      src_lon=lon_rho_o{mo};
      src_lat=lat_rho_o{mo};
      src_mask=mask_rho_o{mo};
      dst_lon=lon_rho_a{ma};
      dst_lat=lat_rho_a{ma};
      dst_mask=mask_rho_a{ma};
      add_new_weights(interp_file,src_lon,src_lat,src_mask,dst_lon,dst_lat,dst_mask)
    end
  end
end
%
%%%%%%%%%%%%%%%%%%%%%   ATM2WAV    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if (0)
  for mw=1:Ngrids_swan
    for ma=1:Ngrids_wrf
      interp_file=['atm',num2str(ma),'_to_','wav',num2str(mw),'_weights.nc']
      src_lon=lon_rho_a{ma};
      src_lat=lat_rho_a{ma};
      src_mask=mask_rho_a{ma};
      dst_lon=lon_rho_w{mw};
      dst_lat=lat_rho_w{mw};
      dst_mask=mask_rho_w{mw};
      add_new_weights(interp_file,src_lon,src_lat,src_mask,dst_lon,dst_lat,dst_mask)
    end
  end
end
%
%%%%%%%%%%%%%%%%%%%%%   WAV2ATM    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if (0)
  for ma=1:Ngrids_wrf
    for mw=1:Ngrids_swan
      interp_file=['wav',num2str(mw),'_to_','atm',num2str(ma),'_weights.nc']
      src_lon=lon_rho_w{mw};
      src_lat=lat_rho_w{mw};
      src_mask=mask_rho_w{mw};
      dst_lon=lon_rho_a{ma};
      dst_lat=lat_rho_a{ma};
      dst_mask=mask_rho_a{ma};
      add_new_weights(interp_file,src_lon,src_lat,src_mask,dst_lon,dst_lat,dst_mask)
    end
  end
end
%



