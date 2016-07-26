function out=maplev(a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file fills the land points or invalid data points using 5-point
% laplacian filter. For water points(valid points), their original
% values are unchanged
%
% im - row of the data
% jm - column of the data
% a - data needed to process (land points are set to be NaN)
% out- results
%
% R.He original on 12/30/99
% R.He revised  on 05/06/01
% R.He revised  on 10/12/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[im,jm]=size(a);

ii=find(~isnan(a));
av=sum(a(ii))/length(ii); %% calculate the mean value of the valid data points

jj=find(isnan(a));
if (length(ii)<2)
 a(jj)=av;           %% find the invalid point and fill them with the mean av
else
%jcw 9/12/2015 - dont use avg to replace, use closest point
 [X,Y]=meshgrid([1:jm],[1:im]);
 a(jj)=griddata(X(ii),Y(ii),a(ii),X(jj),Y(jj),'nearest');
end
b=a;                %% define a working arrey

%% Five point laplacian filter smoothing.
lpp=100;            %% do 100 times smoothing, you may need another number
                    %% More loop, a more smooth field.

for k=1:lpp
 i=[2:im-1]; j=[2:jm-1];
 cc(i,j)=b(i,j)+0.5/4*(b(i+1,j)+b(i,j-1)+b(i-1,j)+b(i,j+1)-4*b(i,j));

 % set the boundary equal ot the next interior points
 cc(1,:) =cc(2,:);
 cc(im,:)=cc(im-1,:);
 cc(:,1) =cc(:,2);
 cc(:,jm)=cc(:,jm-1);

 b(jj)=cc(jj);
end
 a(jj)=cc(jj);    %% only change the invalid data points, keep valid points
                  %% unchange

out=a;

