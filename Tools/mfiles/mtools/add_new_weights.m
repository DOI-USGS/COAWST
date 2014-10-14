function add_new_weights(interp_file,src_lon,src_lat,src_mask,dst_lon,dst_lat,dst_mask)
%
% This m file is part of a series of routines driven by "create_scrip_weights_master.m"
% Checks the weigths in all the files to make sure that all 
% non-masked cells get data from a src grid.
% This is called from check_scrip_weights_driver.m
%
% jcwarner 15Sept2014
%

%   Get some info from the interp file
    dst_address=ncread(interp_file,'dst_address');
    src_address=ncread(interp_file,'src_address');
    src_grid_dims=ncread(interp_file,'src_grid_dims');

%   Check if all the dst non-masked cells have weights that sum to =1.
%   zz are all the indices of the dst grid that need weights.
%   dst_address are the dst grid cells that already have weights.
%   need_wts are the non-masked cells that dont have weights yet.
    zz=find(dst_mask==1);
%   limit these points to be inside the src grid.
    zz2=[];
    X=[src_lon(1,1) src_lon(1,end) src_lon(end,end) src_lon(end,1)];
    Y=[src_lat(1,1) src_lat(1,end) src_lat(end,end) src_lat(end,1)];
    for mm=1:length(zz)
      ind=inpolygon(dst_lon(zz(mm)),dst_lat(zz(mm)),X,Y);
      if (ind); zz2=[zz2 zz(mm)]; end
    end
%
    add_dst_address=setdiff(zz2,dst_address);
    if (~isempty(add_dst_address))
%
%     plot the points we are getting weights for.
      figure
      pcolorjw(dst_lon,dst_lat,dst_mask)
      hold on
      for mm=1:length(add_dst_address)
        [i,j]=ind2ij(dst_mask,add_dst_address(mm));
        plot(dst_lon(i,j),dst_lat(i,j),'go')
      end
%
%     reuse zz. Now compute zz to be the src grid points that we can use.
      zz=find(src_mask==1);
%     set up an interpolation to find these grid indices. only use non-masked
%     points.
      x=src_lon(zz); y=src_lat(zz);
      src_grid_address=[1:1:src_grid_dims(1)*src_grid_dims(2)];
      v=double(src_grid_address(zz)');
      F=TriScatteredInterp(x,y,v,'nearest');
      for mm=1:length(add_dst_address)
        add_src_address(mm)=F(dst_lon(add_dst_address(mm)),dst_lat(add_dst_address(mm)));
        add_remap_matrix(mm)=1;
      end
%
%     plot where the weights are coming from.
      for mm=1:length(add_dst_address)
        [i,j]=ind2ij(src_mask,add_src_address(mm));
        plot(src_lon(i,j),src_lat(i,j),'m+')
      end
      title(['Green circles are the points in ',interp_file,' that had weights added.'])
%
%     put new src_address, dst_address, and src_weights into the netcdf file.
      nc=netcdf.open(interp_file,'NC_WRITE');
%
      disp(' ## Defining Dimensions...')
      netcdf.reDef(nc)
      add_wts = netcdf.defDim(nc,'add_wts',length(add_dst_address));
%
      disp(' ## Defining Variables')
      v1 = netcdf.defVar(nc,'add_src_address','int',add_wts);
      netcdf.putAtt(nc,v1,'long_name','src_address_for_masked_src_cells')
      netcdf.putAtt(nc,v1,'units','---')
%
      v2 = netcdf.defVar(nc,'add_dst_address','int',add_wts);
      netcdf.putAtt(nc,v2,'long_name','dst_address_for_masked_src_cells')
      netcdf.putAtt(nc,v2,'units','---')
%
      v3 = netcdf.defVar(nc,'add_remap_matrix','double',add_wts);
      netcdf.putAtt(nc,v1,'long_name','remap_matrix_for_masked_src_cells')
      netcdf.putAtt(nc,v1,'units','---')
      netcdf.endDef(nc)
      netcdf.close(nc);
%
      disp(' ## Adding new weights')
      ncwrite(interp_file,'add_src_address',single(add_src_address));
      ncwrite(interp_file,'add_dst_address',add_dst_address);
      ncwrite(interp_file,'add_remap_matrix',add_remap_matrix);
    end
end
 