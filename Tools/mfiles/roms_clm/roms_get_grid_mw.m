function grd = roms_get_grid_mw(grd_file,scoord,tindex,calc_zuv)
% $Id$
% grd = roms_get_grid(grd_file)
% grd = roms_get_grid(grd_file,scoord);
% grd = roms_get_grid(grd_file,outfile);
% grd = roms_get_grid(grd_file,outfile,zeta_input,calc_zuv);
%
% Gets the lon,lat,mask,depth [and z coordinates] from netcdf roms netcdf 
% grd_file or output file
% 
% Input:
%     grd_file: The roms netcdf grid file name
%           or, an existing grd structure to which the vertical coordinates 
%               are to be added or updated
%
% Optional inputs:
%     scoord:   ROMS his/rst/avg file from which the s-coord params can be
%               determined
%            or 4-element vector [theta_s theta_b Tcline N]
%
%     zeta_in:  How to obtain zeta information to use
%               when including free surface height in calculating z-coordinates
%            0 implies assume zeta=0
%            integer implies use this time index into his/rst/avg
%               file to read zeta
%            2-d array of zeta values
%
%     calc_zuv: If present, this argument (any value) activates computing
%               the depths z_u and z_v on the u and v points of the 
%               ROMS C-grid
%            
% Output is a structure containing all the grid information
%
% John Wilkin
% Updated (Sept 2002) to correct scoordinate formulation and optionally
% include zeta in the calculation

if isstruct(grd_file)
  % if the first input is already a grd structure
  % the intention is to add vertical coordinates below
  grd = grd_file;
  
else
  % get the horizontal grid information from a ROMS grid file (or
  % history, average etc file)
  grd.grd_file = grd_file;
  
  varlist = ...
    { 'mask_rho','mask_psi','mask_u','mask_v',...
    'h','angle','pm','pn','f'};
  % open the grid file
  for v = varlist
    vname = char(v);
    try
      tmp=ncread(grd_file,vname);
      grd = setfield(grd,vname,tmp);
    catch
      warning(['Variable not found: ' vname])
      if strcmp(vname,'angle')
        grd = setfield(grd,vname,zeros(size(grd.h)));
      end
    end
  end
  
  varlist = {'x_rho','y_rho','x_u','y_u','x_v','y_v','x_psi','y_psi'};
    for v = varlist
      vname = char(v);
      try
        tmp=ncread(grd_file,vname);
        grd = setfield(grd,vname,tmp);
      catch
        warning(['Variable not found: ' vname])
      end
    end
  
  varlist = ...
    { 'lon_rho','lat_rho','lon_psi','lat_psi',...
    'lon_v','lat_v','lon_u','lat_u'};
  for v = varlist
    vname = char(v);
    try
      tmp=ncread(grd_file,vname);
      grd = setfield(grd,vname,tmp);
    catch
      warning([vname ' not found. Substituting x/y coords instead'])
      if strcmp(vname(1:3),'lon')
        usevname = strrep(vname,'lon','x');
      else
        usevname = strrep(vname,'lat','y');
      end
      grd = setfield(grd,vname,getfield(grd,usevname));
    end
  end
    
  if isfield(grd,'mask_rho')
    grd.mask_rho_nan = grd.mask_rho;
    land = find(grd.mask_rho_nan==0);
    grd.mask_rho_nan(land) = NaN;
  else
    % there is no mask information in the file so create unit masks in case
    % code tries to use them
    grd.mask_rho = ones(size(grd.h));
    grd.mask_rho_nan = grd.mask_rho;
    grd.mask_u = ones(size(grd.h(:,2:end)));
    grd.mask_v = ones(size(grd.h(2:end,:)));
    grd.mask_psi = ones(size(grd.h(2:end,2:end)));
  end
  
  % If the grid file includes coastline data, such as a file being used
  % with the Rutgers version of editmask.m, load this too
  try
     grd.lon_coast=ncread(grd_file,'lon_coast');
  catch
  end
  try
      grd.lat_coast=ncread(grd_file,'lat_coast');
  catch
  end
  
end

if nargin > 1  
  
  % get z_r and z_w for the given s-coordinate parameters
  
  if ~ischar(scoord)
   
    theta_s = scoord(1);
    theta_b = scoord(2);
    Tcline  = scoord(3);
    N       = scoord(4);
    h = grd.h;
    
    % code lifted from hernan's scoord3.m
    c1=1.0;
    c2=2.0;
    p5=0.5;
    Np=N+1;
    ds=1.0/N;
    hmin=min(min(h));
    hmax=max(max(h));
    hc=min(hmin,Tcline);
    [Mp Lp]=size(h);
    % rho points
    Nlev=N;
    lev=1:N;
    sc=-c1+(lev-p5).*ds;
    Ptheta=sinh(theta_s.*sc)./sinh(theta_s);
    Rtheta=tanh(theta_s.*(sc+p5))./(c2*tanh(p5*theta_s))-p5;
    Cs=(c1-theta_b).*Ptheta+theta_b.*Rtheta;
    sc_r = sc(:);
    Cs_r = Cs(:);    
    % w points
    Nlev=Np;
    lev=0:N;
    sc=-c1+lev.*ds;
    Ptheta=sinh(theta_s.*sc)./sinh(theta_s);
    Rtheta=tanh(theta_s.*(sc+p5))./(c2*tanh(p5*theta_s))-p5;
    Cs=(c1-theta_b).*Ptheta+theta_b.*Rtheta;
    sc_w = sc(:);
    Cs_w = Cs(:);

  else
  
    % input 'scoord' is the name of a his/avg/rst file name
    % so attempt to get s-coord params from the file
    theta_s = ncread(grd_file,'theta_s');
    theta_b = ncread(grd_file,'theta_b');
    Tcline = ncread(grd_file,'Tcline');
    sc_r = ncread(grd_file,'s_rho');
    sc_w = ncread(grd_file,'s_w'); 
    Cs_w = ncread(grd_file,'Cs_w');
    Cs_r = ncread(grd_file,'Cs_r');
    N = length(sc_r);
    Np = N+1;
    hc = ncread(grd_file,'hc');
    if isempty(hc)
      hc = min(h(:));
    end
  
  end

  % For compatibility with new syntax in CF compliant files
  s_w = sc_w;
  s_rho = sc_r;
  
  % zeta  
  zeta = zeros(size(grd.h)); % default
  if nargin > 2 % option to include zeta in z calculation
    
    if tindex == 0 % if tindex==0 zeta defaults to zero
      % do nothing
      
    else % if tindex==0 zeta defaults to zero
      if length(tindex)==1
        % tindex is a single index to zeta in a roms output file
        if ~ischar(scoord)
          error([ 'Can''t process zeta from file in the case that ' ...
            ' scoord parameters are input as a vector'])'
        end
        zeta = nc.data('zeta',[tindex-1 0 0],[1 -1 -1]);
        if isempty(zeta)
          warning([ 'zeta not found in ' scoord '. Assuming zeta=0.'])
          zeta = zeros(size(grd.h));
        end
      else
        % tindex should be a 2-d field of zeta values
        if any(size(tindex)-size(grd.h))
          % sizes of zeta and h don't match
          error('input tindex as zeta does not match grid dimensions')
        else
          zeta = tindex;
        end
      end
    end
    
  end
  grd.zeta = zeta;
  
  % rho-points
  h = grd.h;
  scmCshc = (sc_r-Cs_r)*hc;
  z_r = repmat(scmCshc,[1 length(h(:))]) + Cs_r*h(:)';
  if any(zeta(:)~=0)
    z_r = z_r + scmCshc*[zeta(:)./h(:)]' + (1+Cs_r)*zeta(:)';
  end
  grd.z_r = reshape(z_r,[N size(h)]);
  
  % w-points
  scmCshc_w = (sc_w-Cs_w)*hc;
  z_w = repmat(scmCshc_w,[1 length(h(:))]) + Cs_w*h(:)';
  if any(zeta(:)~=0)
    z_w = z_w + scmCshc_w*[zeta(:)./h(:)]' + (1+Cs_w)*zeta(:)';
  end
  grd.z_w = reshape(z_w,[Np size(h)]);
  clear z_r z_w
  
  if nargin > 3
    
    % u-points (cell centres in vertical)
    hu = 0.5*(h(:,1:end-1)+h(:,2:end));
    zu = 0.5*(zeta(:,1:end-1)+zeta(:,2:end));
    z_u = repmat(scmCshc,[1 length(hu(:))]) + Cs_r*hu(:)';
    if any(zu(:)~=0)
      z_u = z_u + scmCshc*[zu(:)./hu(:)]' + (1+Cs_r)*zu(:)';
    end
    grd.z_u = reshape(z_u,[N size(hu)]);
    clear z_u;
    
    % u-points (cell edges in vertical)
    z_uw = repmat(scmCshc_w,[1 length(hu(:))]) + Cs_w*hu(:)';
    if any(zu(:)~=0)
      z_uw = z_uw + scmCshc_w*[zu(:)./hu(:)]' + (1+Cs_w)*zu(:)';
    end
    grd.z_uw = reshape(z_uw,[Np size(hu)]);
    clear z_uw;
    
    % v-points (cell centres in vertical)
    hv = 0.5*(h(1:end-1,:)+h(2:end,:));
    zv = 0.5*(zeta(1:end-1,:)+zeta(2:end,:));
    z_v = repmat(scmCshc,[1 length(hv(:))]) + Cs_r*hv(:)';
    if any(zeta(:)~=0)
      z_v = z_v + scmCshc*[zv(:)./hv(:)]' + (1+Cs_r)*zv(:)';
    end
    grd.z_v = reshape(z_v,[N size(hv)]);
    
    % v-points (cell edges in vertical)
    z_vw = repmat(scmCshc_w,[1 length(hv(:))]) + Cs_w*hv(:)';
    if any(zv(:)~=0)
      z_vw = z_vw + scmCshc_w*[zv(:)./hv(:)]' + (1+Cs_w)*zv(:)';
    end
    grd.z_vw = reshape(z_vw,[Np size(hv)]);
    clear z_vw;
    
  end 
  grd.theta_s = theta_s;
  grd.theta_b = theta_b;
  grd.Tcline = Tcline;
  grd.N = N;
  grd.hc = hc;
  grd.sc_w = sc_w;
  grd.Cs_w = Cs_w;
  grd.sc_r = sc_r;
  grd.Cs_r = Cs_r;
  grd.s_w = s_w;
  grd.s_rho = s_rho;
  
end

