function [specpts]=ww3_specpoints(modelgrid,spec_res);
%*************************************************************************
%this function uses the model grid to determine specpoints along the 
%boundary of the grid which are then used in the loaddap_ww3_2TPAR function
%
%modelgrid=name/location of grid used in the model
%spec_res=number of grid points between spec points/resolution of spec
%points
%
%the output for this functions is saved in the working directory as
%specpts.mat
%
%created by Brandy Armstrong 4/8/09 for use in COAWST system v2
%modified jcw 05Feb2012 to use native matlab netcdf
%*************************************************************************

%load in grid info
gx=ncread(modelgrid,'lon_rho');gx2=gx;%gx2(gx<0) = gx2(gx<0)+360;
gy=ncread(modelgrid,'lat_rho');
[xi eta]=size(gy);
mask=ncread(modelgrid,'mask_rho');
mask(isnan(mask))=1;
if size(mask)==0;
    mask=ones(xi,eta);
end
gx3=gx2.*mask;
gy3=gy.*mask;
rr=1;
offset=0.001;
for ct=1:spec_res:eta
    if gx3(1,ct)~=0 %determine if the point is masked/land
        if ct+spec_res<eta
            specpts(rr,1)=gx2(1,ct)+offset;
            specpts(rr,2)=gy(1,ct)+offset;
            specpts(rr,3)=gx2(1,ct+spec_res)+offset;
            specpts(rr,4)=gy(1,ct+spec_res)+offset;
        else
            specpts(rr,1)=gx2(1,ct)+offset;
            specpts(rr,2)=gy(1,ct)+offset;
            specpts(rr,3)=gx2(1,end)+offset;
            specpts(rr,4)=gy(1,end)+offset;
        end
        rr=rr+1;
    end

end
for ct=1:spec_res:eta
    if gx3(end,ct)~=0 %determine if the point is masked/land
        if ct+spec_res<eta
            specpts(rr,1)=gx2(end,ct)-offset;
            specpts(rr,2)=gy(end,ct)-offset;
            specpts(rr,3)=gx2(end,ct+spec_res)-offset;
            specpts(rr,4)=gy(end,ct+spec_res)-offset;
        else
            specpts(rr,1)=gx2(end,ct)-offset;
            specpts(rr,2)=gy(end,ct)-offset;
            specpts(rr,3)=gx2(end,end)-offset;
            specpts(rr,4)=gy(end,end)-offset;
        end
        rr=rr+1;
    end

end
for ct=1:spec_res:xi
    if gx3(ct,1)~=0 %determine if the point is masked/land
        if ct+spec_res<xi
            specpts(rr,1)=gx2(ct,1);
            specpts(rr,2)=gy(ct,1)+offset;
            specpts(rr,3)=gx2(ct+spec_res,1);
            specpts(rr,4)=gy(ct+spec_res,1)+offset;
        else
            specpts(rr,1)=gx2(ct,1);
            specpts(rr,2)=gy(ct,1)+offset;
            specpts(rr,3)=gx2(end,1)-offset;
            specpts(rr,4)=gy(end,1)+offset;
        end
        rr=rr+1;
    end

end
for ct=1:spec_res:xi
    if gx3(ct,end)~=0 %determine if the point is masked/land
        if ct+spec_res<xi
            specpts(rr,1)=gx2(ct,end)-offset;
            specpts(rr,2)=gy(ct,end);
            specpts(rr,3)=gx2(ct+spec_res,end)-offset;
            specpts(rr,4)=gy(ct+spec_res,end);
        else
            specpts(rr,1)=gx2(ct,end)-offset;
            specpts(rr,2)=gy(ct,end);
            specpts(rr,3)=gx2(end,end)-offset;
            specpts(rr,4)=gy(end,end);
        end
        rr=rr+1;
    end
end

%flip specpoints so that it is consistent with what ncgeodataset is look for.
%specpts=specpts.';
%save file as mat files
save specpts.mat specpts
