function [xr,yr]=rot(x,y,theta,ya);
%function [xr,yr]=rot(x,y,theta,[ya]);
%or  [xr]=rot(x,theta,[ya]); (where x is complex)
%      rotates a vector counterclockwise theta degrees
% OR   rotates the coord system clockwise theta degrees
%    rot(1,0,90) returns (0,1)
%   if 4 arguments then it uses xa and ya as orientation

if (nargin==2)|(length(x)~=length(y))
   if nargin==3; ya=theta; narg=4; end
   theta=y;
   y=imag(x);
   x=real(x);
   comp=1;
   narg=2;
else; 
   comp=0; 
   narg=nargin;
end

if narg==4
theta=-atan2(ya,theta)*180/pi;
end
if imag(theta)~=0
theta=-atan2(imag(theta),real(theta))*180/pi;
end
costheta=cos(theta/180*pi);
sintheta=sin(theta/180*pi);

if length(theta)==1
    xr=x*costheta-y*sintheta;
    yr=x*sintheta+y*costheta;
else
    xr=x.*costheta-y.*sintheta;
    yr=x.*sintheta+y.*costheta;
end

if comp==1;
    xr=xr+i*yr;
    yr=[];
end
