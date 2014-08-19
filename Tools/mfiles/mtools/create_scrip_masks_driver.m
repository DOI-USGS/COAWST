% script create_scrip_masks_driver
%
% This m file is part of a series of routines driven by "create_scrip_weights_master.m"
%
% This routine computes the source and destination masks for each model.
%
% jcwarner14July2014
%

eval(['cd ',wdir,' '])

%%%%%%%%%%%%%%%%%%%%%   SWAN Model   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For SWAN, load some grid information into a structure.
for mw=1:Ngrids_swan
  bf=bath_file{mw};
  bath{mw}=textread(bf)';
  zz=ones(size(bath{mw}));
  zz(bath{mw}==9999)=0;
  mask_rho_w{mw}=zz;
%
  nx=Numx{mw};
  ny=Numy{mw};
  grd=textread(swan_coord{mw});
  gridsize=length(grd)/2;
  zz=grd(1:gridsize);
  lon_rho_w{mw}=reshape(zz,nx,ny);
  zz=grd(gridsize+1:end);
  lat_rho_w{mw}=reshape(zz,nx,ny);
end

% For SWAN, find out start/end indices of where the children grids are.
% Istrc,Iendc and Jstrc,Jendc are the bounding indices in the parent grid
% where the child is located.
 for mw=1:Ngrids_swan-1
  nx=Numx{mw};
  ny=Numy{mw};
  zz=[1:1:nx*ny];
  zz=reshape(zz,nx,ny);
  zt=griddata(lon_rho_w{mw},lat_rho_w{mw},zz,lon_rho_w{mw+1}(1,1),lat_rho_w{mw+1}(1,1),'nearest');
  [Istr_w{mw},Jstr_w{mw}]=ind2ij(zz,zt);
  zt=griddata(lon_rho_w{mw},lat_rho_w{mw},zz,lon_rho_w{mw+1}(end,end),lat_rho_w{mw+1}(end,end),'nearest');
  [Iend_w{mw},Jend_w{mw}]=ind2ij(zz,zt);
end

%%%%%%%%%%%%%%%%%%%%%   ROMS Model   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For ROMS, load some grid information into a structure.
for mo=1:Ngrids_roms
  rg=roms_grids{mo};
  spherical=ncread(rg,'spherical');
  if ((spherical=='F')||(spherical=='f'))||(spherical==0)
    lon_rho_o{mo}=ncread(rg,'x_rho');
    lat_rho_o{mo}=ncread(rg,'y_rho');
  else
    lon_rho_o{mo}=ncread(rg,'lon_rho');
    lat_rho_o{mo}=ncread(rg,'lat_rho');
  end
  mask_rho_o{mo}=ncread(rg,'mask_rho');
end

% For ROMS, find out start/end indices of where the children grids are.
% Istr,Iend and Jstr,Jend are the bounding indices in the parent grid
% where the child is located.
for mo=1:Ngrids_roms-1
  [nx,ny]=size(mask_rho_o{mo});
  zz=[1:1:nx*ny];
  zz=reshape(zz,nx,ny);
  zt=griddata(lon_rho_o{mo},lat_rho_o{mo},zz,lon_rho_o{mo+1}(1,1),lat_rho_o{mo+1}(1,1),'nearest');
  [Istr_o{mo},Jstr_o{mo}]=ind2ij(zz,zt);
  zt=griddata(lon_rho_o{mo},lat_rho_o{mo},zz,lon_rho_o{mo+1}(end,end),lat_rho_o{mo+1}(end,end),'nearest');
  [Iend_o{mo},Jend_o{mo}]=ind2ij(zz,zt);
end

%%%%%%%%%%%%%%%%%%%%%   WRF Model   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For WRF, load some grid information into a structure.
for ma=1:Ngrids_wrf
  wg=wrf_grids{ma};
  lon_rho_a{ma}=double(ncread(wg,'XLONG'));
  lat_rho_a{ma}=double(ncread(wg,'XLAT'));
  mask_rho_a{ma}=double(1-ncread(wg,'LANDMASK'));
end

% For WRF, find out start/end indices of where the children grids are.
% Istr,Iend and Jstr,Jend are the bounding indices in the parent grid
% where the child is located.
for ma=1:Ngrids_wrf-1
  [nx,ny]=size(mask_rho_a{ma});
  zz=[1:1:nx*ny];
  zz=reshape(zz,nx,ny);
  zt=griddata(lon_rho_a{ma},lat_rho_a{ma},zz,lon_rho_a{ma+1}(1,1),lat_rho_a{ma+1}(1,1),'nearest');
  [Istr_a{ma},Jstr_a{ma}]=ind2ij(zz,zt);
  zt=griddata(lon_rho_a{ma},lat_rho_a{ma},zz,lon_rho_a{ma+1}(end-1,end-1),lat_rho_a{ma+1}(end,end),'nearest');
  [Iend_a{ma},Jend_a{ma}]=ind2ij(zz,zt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Now create the masks between each model %%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%   Ocn2Wav    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%: Determine destination masks for the swan grids uisng ocn as the source.
for mw=1:Ngrids_swan
  tot_mask{mw}=0;
  for mo=1:Ngrids_roms
      src_mask=mask_rho_o{mo};
      if Ngrids_roms>mo  %mask out child portion of this ocn grid
        src_mask(Istr_o{mo}:Iend_o{mo},Jstr_o{mo}:Jend_o{mo})=0;
      end
%
%  Compute destination masking on destination grid. Only use locations on 
%  dst grid where the src grid masking = 1.
%
      dst_mask=griddata(lon_rho_o{mo},lat_rho_o{mo},src_mask,lon_rho_w{mw},lat_rho_w{mw},'nearest');
%
% Set all points in the dst mask that are outside the src grid to = 0.
%
    X=[lon_rho_o{mo}(1,1) lon_rho_o{mo}(1,end) lon_rho_o{mo}(end,end) lon_rho_o{mo}(end,1)];
    Y=[lat_rho_o{mo}(1,1) lat_rho_o{mo}(1,end) lat_rho_o{mo}(end,end) lat_rho_o{mo}(end,1)];
    ind=inpolygon(lon_rho_w{mw}(:),lat_rho_w{mw}(:),X,Y);
    dst_mask(ind==0)=0;
    tot_mask{mw}=tot_mask{mw}+dst_mask;
%
    eval(['O',num2str(mo),'src_W',num2str(mw),'dst=dst_mask;'])
  end
%   Add land/sea mask to total mask
    tot_mask{mw}=tot_mask{mw}+(1-mask_rho_w{mw});
    if(sum(1-tot_mask{mw}(:))~=0)
      disp('had to adjust totmask for O1src_W1dst')
      zz=find(tot_mask{mw}~=1);
      eval(['O1src_W',num2str(mw),'dst(zz)=1;'])
      for mo=2:Ngrids_roms
        eval(['O',num2str(mo),'src_W',num2str(mw),'dst(zz)=0;'])
      end
      tot_mask{mw}(zz)=1;
    end
end

%plot out SWAN masks
for mw=1:Ngrids_swan
  figure
  for mo=1:Ngrids_roms
    subplot(1,Ngrids_roms+2,mo)
    eval(['dst=O',num2str(mo),'src_W',num2str(mw),'dst;'])
    pcolorjw(lon_rho_w{mw},lat_rho_w{mw},dst)
    title(['mask from OCN grid ',num2str(mo)])
    colorbar
  end
  subplot(1,Ngrids_roms+2,Ngrids_roms+1)
    pcolorjw(lon_rho_w{mw},lat_rho_w{mw},1-mask_rho_w{mw})
    title(['Landmask in SWAN grid ',num2str(mw)])
    colorbar
  subplot(1,Ngrids_roms+2,Ngrids_roms+2)
    pcolorjw(lon_rho_w{mw},lat_rho_w{mw},tot_mask{mw})
    title(['TOTAL masking for SWAN grid ',num2str(mw)])
    colorbar
  set(gcf,'name',['SWAN grid ',num2str(mw),' destination masks'])
  set(gcf,'renderer','painters')
end

%%%%%%%%%%%%%%%%%%%%%   Wav2Ocn    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determine destination masks for the ROMS grids.
%
for mo=1:Ngrids_roms
  tot_mask{mo}=0;
  for mw=1:Ngrids_swan
% Get the src mask.
      src_mask=mask_rho_w{mw};
      if Ngrids_swan>mw  %mask out child portion of this wave grid
        src_mask(Istr_w{mw}:Iend_w{mw},Jstr_w{mw}:Jend_w{mw})=0;
      end
%
%  Compute destination masking on destination grid. Only use locations on 
%  dst grid where the src grid masking = 1.
%
      dst_mask=griddata(lon_rho_w{mw},lat_rho_w{mw},src_mask,lon_rho_o{mo},lat_rho_o{mo},'nearest');
%
% Set all points in the dst mask that are outside the src grid to = 0.
%
    X=[lon_rho_w{mw}(1,1) lon_rho_w{mw}(1,end) lon_rho_w{mw}(end,end) lon_rho_w{mw}(end,1)];
    Y=[lat_rho_w{mw}(1,1) lat_rho_w{mw}(1,end) lat_rho_w{mw}(end,end) lat_rho_w{mw}(end,1)];
    ind=inpolygon(lon_rho_o{mo}(:),lat_rho_o{mo}(:),X,Y);
    dst_mask(ind==0)=0;
    tot_mask{mo}=tot_mask{mo}+dst_mask;
    eval(['W',num2str(mw),'src_O',num2str(mo),'dst=dst_mask;'])
  end
% Add land/sea mask to total mask
    tot_mask{mo}=tot_mask{mo}+(1-mask_rho_o{mo});
    if(sum(1-tot_mask{mo}(:))~=0)
      disp(['had to adjust totmask for W1src_O',num2str(mo),'dst'])
      zz=find(tot_mask{mo}~=1);
      eval(['W1src_O',num2str(mo),'dst(zz)=1;'])
      for mw=2:Ngrids_swan
        eval(['W',num2str(mw),'src_O',num2str(mo),'dst(zz)=0;'])
      end
      tot_mask{mo}(zz)=1;
    end
end

%plot out ROMS masks
for mo=1:Ngrids_roms
  figure
  for mw=1:Ngrids_swan
    subplot(1,Ngrids_swan+2,mw)
    eval(['dst=W',num2str(mw),'src_O',num2str(mo),'dst;'])
    pcolorjw(lon_rho_o{mo},lat_rho_o{mo},dst)
    title(['mask from SWAN grid ',num2str(mw)])
    colorbar
  end
  subplot(1,Ngrids_swan+2,Ngrids_swan+1)
    pcolorjw(lon_rho_o{mo},lat_rho_o{mo},1-mask_rho_o{mo})
    title(['Landmask in OCN grid ',num2str(mo)])
    colorbar
  subplot(1,Ngrids_swan+2,Ngrids_swan+2)
    pcolorjw(lon_rho_o{mo},lat_rho_o{mo},tot_mask{mo})
    title(['TOTAL masking for OCN grid ',num2str(mo)])
    colorbar
    set(gcf,'name',['ROMS grid ',num2str(mo),' destination masks'])
    set(gcf,'renderer','painters')
end

%%%%%%%%%%%%%%%%%%%%%   Atm2Ocn    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determine destination masks for the ROMS grids using atm as the src.
%
for mo=1:Ngrids_roms
  tot_mask{mo}=0;
  for ma=1:Ngrids_wrf
% Get the src mask.
      src_mask=mask_rho_a{ma};
      if Ngrids_wrf>ma  %mask out child portion of this atm grid
        src_mask(Istr_a{ma}:Iend_a{ma},Jstr_a{ma}:Jend_a{ma})=0;
      end
%
%  Compute destination masking on destination grid. Only use locations on 
%  dst grid where the src grid masking = 1.
%
    dst_mask=griddata(lon_rho_a{ma},lat_rho_a{ma},src_mask,lon_rho_o{mo},lat_rho_o{mo},'nearest');
%
% Set all points in the dst mask that are outside the src grid to = 0.
%
    X=[lon_rho_a{ma}(1,1) lon_rho_a{ma}(1,end) lon_rho_a{ma}(end,end) lon_rho_a{ma}(end,1)];
    Y=[lat_rho_a{ma}(1,1) lat_rho_a{ma}(1,end) lat_rho_a{ma}(end,end) lat_rho_a{ma}(end,1)];
    ind=inpolygon(lon_rho_o{mo}(:),lat_rho_o{mo}(:),X,Y);
    dst_mask(ind==0)=0;
    tot_mask{mo}=tot_mask{mo}+dst_mask;
%
    eval(['A',num2str(ma),'src_O',num2str(mo),'dst=dst_mask;'])
  end
%  Add land/sea mask to total mask
    tot_mask{mo}=tot_mask{mo}+(1-mask_rho_o{mo});
    if(sum(1-tot_mask{mo}(:))~=0)
      disp('had to adjust totmask for A1src_O1dst')
      zz=find(tot_mask{mo}~=1);
      eval(['A1src_O',num2str(mo),'dst(zz)=1;'])
      for ma=2:Ngrids_wrf
        eval(['A',num2str(ma),'src_O',num2str(mo),'dst(zz)=0;'])
      end
      tot_mask{mo}(zz)=1;
    end
end

%plot out ROMS masks
for mo=1:Ngrids_roms
  figure
  for ma=1:Ngrids_wrf
    subplot(1,Ngrids_wrf+2,ma)
    eval(['dst=A',num2str(ma),'src_O',num2str(mo),'dst;'])
    pcolorjw(lon_rho_o{mo},lat_rho_o{mo},dst)
    title(['mask from WRF grid ',num2str(ma)])
    colorbar
  end
  subplot(1,Ngrids_wrf+2,Ngrids_wrf+1)
    pcolorjw(lon_rho_o{mo},lat_rho_o{mo},1-mask_rho_o{mo})
    title(['Landmask in ROMS grid ',num2str(mo)])
    colorbar
  subplot(1,Ngrids_wrf+2,Ngrids_wrf+2)
    pcolorjw(lon_rho_o{mo},lat_rho_o{mo},tot_mask{mo})
    title(['TOTAL masking for ROMS grid ',num2str(mo)])
    colorbar
  set(gcf,'name',['ROMS grid ',num2str(mo),' destination masks'])
  set(gcf,'renderer','painters')
end


%%%%%%%%%%%%%%%%%%%%%   Ocn2Atm    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determine destination masks for the WRF grids using roms as the src.
%
for ma=1:Ngrids_wrf
  tot_mask{ma}=0;
  for mo=1:Ngrids_roms

% Get the src mask.
      src_mask=mask_rho_o{mo};
      if Ngrids_roms>mo  %mask out child portion of this ocean grid
        src_mask(Istr_o{mo}:Iend_o{mo},Jstr_o{mo}:Jend_o{mo})=0;
      end
%
%  Compute destination masking on destination grid. Only use locations on 
%  dst grid where the src grid masking = 1.
%
      dst_mask=griddata(lon_rho_o{mo},lat_rho_o{mo},src_mask,lon_rho_a{ma},lat_rho_a{ma},'nearest');
%
% Set all points in the dst mask that are outside the src grid to = 0.
%
    X=[lon_rho_o{mo}(1,1) lon_rho_o{mo}(1,end) lon_rho_o{mo}(end,end) lon_rho_o{mo}(end,1)];
    Y=[lat_rho_o{mo}(1,1) lat_rho_o{mo}(1,end) lat_rho_o{mo}(end,end) lat_rho_o{mo}(end,1)];
    ind=inpolygon(lon_rho_a{ma}(:),lat_rho_a{ma}(:),X,Y);
    dst_mask(ind==0)=0;
    tot_mask{ma}=tot_mask{ma}+dst_mask;
    eval(['O',num2str(mo),'src_A',num2str(ma),'dst=dst_mask;'])
%
%   modify the wrfinput_d** cplmask  - Hold off on this for now. 
%   if (mo==1)
%     eval(['cplmask{ma}=ncread(''wrfinput_d0',num2str(ma),''',''CPLMASK'');'])
      cplmask{ma}=1-dst_mask;
%     eval(['ncwrite(''wrfinput_d0',num2str(ma),''',''CPLMASK'',cplmask{ma});'])
%   end
  end
% Add land/sea mask to total mask
%   tot_mask{ma}=tot_mask{ma}+(1-mask_rho_a{ma});
    tot_mask{ma}=tot_mask{ma}+cplmask{ma};
    if(sum(1-tot_mask{ma}(:))~=0)
      disp(['had to adjust totmask for O1src_A',num2str(ma),'dst'])
      zz=find(tot_mask{ma}~=1);
      eval(['O1src_A',num2str(ma),'dst(zz)=1;'])
      for mo=2:Ngrids_roms
        eval(['O',num2str(mo),'src_A',num2str(ma),'dst(zz)=0;'])
      end
      tot_mask{ma}(zz)=1;
    end
end

%plot out WRF masks
for ma=1:Ngrids_wrf
  figure
  for mo=1:Ngrids_roms
    subplot(1,Ngrids_roms+2,mo)
    eval(['dst=O',num2str(mo),'src_A',num2str(ma),'dst;'])
    pcolorjw(lon_rho_a{ma},lat_rho_a{ma},dst)
    title(['mask from ROMS grid ',num2str(mo)])
    colorbar
  end
  subplot(1,Ngrids_roms+2,Ngrids_roms+1)
    pcolorjw(lon_rho_a{ma},lat_rho_a{ma},cplmask{ma})
    title(['CPLMASK in ATM grid ',num2str(ma)])
    colorbar
  subplot(1,Ngrids_roms+2,Ngrids_roms+2)
    pcolorjw(lon_rho_a{ma},lat_rho_a{ma},tot_mask{ma})
    title(['TOTAL masking for WRF grid ',num2str(ma)])
    colorbar
  set(gcf,'name',['WRF grid ',num2str(ma),' destination masks'])
  set(gcf,'renderer','painters')
end


%%%%%%%%%%%%%%%%%%%%%   Atm2Wav    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determine destination masks for the SWAN grids using atm as the source.
%
for mw=1:Ngrids_swan
  tot_mask{mw}=0;
  for ma=1:Ngrids_wrf
      src_mask=mask_rho_a{ma};
      if Ngrids_wrf>ma  %mask out child portion of this atm grid
        src_mask(Istr_a{ma}:Iend_a{ma},Jstr_a{ma}:Jend_a{ma})=0;
      end
%
%  Compute destination masking on destination grid. Only use locations on 
%  dst grid where the src grid masking = 1.
%
      dst_mask=griddata(lon_rho_a{ma},lat_rho_a{ma},src_mask,lon_rho_w{mw},lat_rho_w{mw},'nearest');
%
% Set all points in the dst mask that are outside the src grid to = 0.
%
    X=[lon_rho_a{ma}(1,1) lon_rho_a{ma}(1,end) lon_rho_a{ma}(end,end) lon_rho_a{ma}(end,1)];
    Y=[lat_rho_a{ma}(1,1) lat_rho_a{ma}(1,end) lat_rho_a{ma}(end,end) lat_rho_a{ma}(end,1)];
    ind=inpolygon(lon_rho_w{mw}(:),lat_rho_w{mw}(:),X,Y);
    dst_mask(ind==0)=0;
    tot_mask{mw}=tot_mask{mw}+dst_mask;
%
    eval(['A',num2str(ma),'src_W',num2str(mw),'dst=dst_mask;'])
  end
%   Add land/sea mask to total mask
    tot_mask{mw}=tot_mask{mw}+(1-mask_rho_w{mw});
    if(sum(1-tot_mask{mw}(:))~=0)
      disp('had to adjust totmask for A1src_W1dst')
      zz=find(tot_mask{mw}~=1);
      eval(['A1src_W',num2str(mw),'dst(zz)=1;'])
      for ma=2:Ngrids_wrf
        eval(['A',num2str(ma),'src_W',num2str(mw),'dst(zz)=0;'])
      end
      tot_mask{mw}(zz)=1;
    end
end

%plot out SWAN masks
for mw=1:Ngrids_swan
  figure
  for ma=1:Ngrids_wrf
    subplot(1,Ngrids_wrf+2,ma)
    eval(['dst=A',num2str(ma),'src_W',num2str(mw),'dst;'])
    pcolorjw(lon_rho_w{mw},lat_rho_w{mw},dst)
    title(['mask from WRF grid ',num2str(ma)])
    colorbar
  end
  subplot(1,Ngrids_wrf+2,Ngrids_wrf+1)
    pcolorjw(lon_rho_w{mw},lat_rho_w{mw},1-mask_rho_w{mw})
    title(['Landmask in SWAN grid ',num2str(mw)])
    colorbar
  subplot(1,Ngrids_wrf+2,Ngrids_wrf+2)
    pcolorjw(lon_rho_w{mw},lat_rho_w{mw},tot_mask{mw})
    title(['TOTAL masking for SWAN grid ',num2str(mw)])
    colorbar
  set(gcf,'name',['SWAN grid ',num2str(mw),' destination masks'])
  set(gcf,'renderer','painters')
end


%%%%%%%%%%%%%%%%%%%%%   Wav2Atm    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%: Determine destination masks for the atm grids.
for ma=1:Ngrids_wrf
  tot_mask{ma}=0;
  for mw=1:Ngrids_swan
% Get the src mask.
      src_mask=mask_rho_w{mw};
      if Ngrids_swan>mw  %mask out child portion of this wave grid
        src_mask(Istr_w{mw}:Iend_w{mw},Jstr_w{mw}:Jend_w{mw})=0;
      end
%
%  Compute destination masking on destination grid. Only use locations on 
%  dst grid where the src grid masking = 1.
%
      dst_mask=griddata(lon_rho_w{mw},lat_rho_w{mw},src_mask,lon_rho_a{ma},lat_rho_a{ma},'nearest');
%
% Set all points in the dst mask that are outside the src grid to = 0.
%
    X=[lon_rho_w{mw}(1,1) lon_rho_w{mw}(1,end) lon_rho_w{mw}(end,end) lon_rho_w{mw}(end,1)];
    Y=[lat_rho_w{mw}(1,1) lat_rho_w{mw}(1,end) lat_rho_w{mw}(end,end) lat_rho_w{mw}(end,1)];
    ind=inpolygon(lon_rho_a{ma}(:),lat_rho_a{ma}(:),X,Y);
    dst_mask(ind==0)=0;
    tot_mask{ma}=tot_mask{ma}+dst_mask;
    eval(['W',num2str(mw),'src_A',num2str(ma),'dst=dst_mask;'])
%
%   modify the wrfinput_d** cplmask
%   if (mw==1)
%     eval(['cplmask{ma}=ncread(''wrfinput_d0',num2str(ma),''',''CPLMASK'');'])
      cplmask{ma}=1-dst_mask;
%     eval(['ncwrite(''wrfinput_d0',num2str(ma),''',''CPLMASK'',cplmask{ma});'])
%   end
  end
% Add land/sea mask to total mask
%   tot_mask{ma}=tot_mask{ma}+(1-mask_rho_a{ma});
    tot_mask{ma}=tot_mask{ma}+cplmask{ma};
    if(sum(1-tot_mask{ma}(:))~=0)
      disp(['had to adjust totmask for W1src_A',num2str(ma),'dst'])
      zz=find(tot_mask{ma}~=1);
      eval(['W1src_A',num2str(ma),'dst(zz)=1;'])
      for mw=2:Ngrids_swan
        eval(['W',num2str(mw),'src_A',num2str(ma),'dst(zz)=0;'])
      end
      tot_mask{ma}(zz)=1;
    end
end

%plot out WRF masks
for ma=1:Ngrids_wrf
  figure
  for mw=1:Ngrids_swan
    subplot(1,Ngrids_swan+2,mw)
    eval(['dst=W',num2str(mw),'src_A',num2str(ma),'dst;'])
    pcolorjw(lon_rho_a{ma},lat_rho_a{ma},dst)
    title(['mask from SWAN grid ',num2str(mw)])
    colorbar
  end
  subplot(1,Ngrids_swan+2,Ngrids_swan+1)
    pcolorjw(lon_rho_a{ma},lat_rho_a{ma},cplmask{ma})
    title(['CPLMASK in ATM grid ',num2str(ma)])
    colorbar
  subplot(1,Ngrids_swan+2,Ngrids_swan+2)
    pcolorjw(lon_rho_a{ma},lat_rho_a{ma},tot_mask{ma})
    title(['TOTAL masking for WRF grid ',num2str(ma)])
    colorbar
  set(gcf,'name',['WRF grid ',num2str(ma),' destination masks'])
  set(gcf,'renderer','painters')
end


