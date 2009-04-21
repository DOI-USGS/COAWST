function pos = roms_cgridpos(s,grd)
% $Id: roms_cgridpos.m 358 2008-04-07 14:15:03Z zhang $
% function pos = roms_cgridpos(s,grd)
%
% From the size of a variable determine whether it is defined on the
% u, v, rho or psi location on the ROMS Arakawa-C grid
% 
% John Wilkin

if ndims(s)~=2
  s = size(s);
else
  ss = size(s);
  if ss(1)~= 1
    s = size(s);
  end
end

su = size(grd.lon_u);
sv = size(grd.lon_v);
sr = size(grd.lon_rho);
sp = size(grd.lon_psi);
if sr(1)==sr(2)
  warning([ 'The grid is square so the result of ' ...
	which(mfilename) ' may be unreliable']);
end

% only the right-most 2 dimensions matter
s = s([end-1 end]);
if all(~(s-sv))
  pos = 'v';
elseif all(~(s-su))
  pos = 'u';
elseif all(~(s-sr))
  pos = 'rho';
elseif all(~(s-sp))
  pos = 'psi';
else
  error([ 'Could not determine C-grid location'])
end

