function  [X,Y,vals,labI]=mp_tmerc(optn,varargin)
% MP_TMERC   Transverse Mercator projection
%           This function should not be used directly; instead it is
%           is accessed by various high-level functions named M_*.


% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
% Mathematical formulas for the projections and their inverses are taken from
%
%      Snyder, John P., Map Projections used by the US Geological Survey, 
%      Geol. Surv. Bull. 1532, 2nd Edition, USGPO, Washington D.C., 1983.
%
% It was handy to include both the TM (a cylindrical projection) and the
% so-called "pseudo-cylindrical" sinusoidal projection here, and since we have
% the sinusoidal, I started including the other "global" projections like the hammer.
%
% Transverse Mercator - cylindrical conformal
% Sinusoidal - cylindrical equal-area
% Gall-Peters - an area-conserving rectangular projection.
% Hammer-Aitoff - another equal-area projection (for global pics)

% 7/2/98 - Added defaults for Sinusoidal projection
% 15/4/98 - Gall-Peters projection
% 8/8/98 - Hammer-Aitoff
% 17/8/98 - Mollweide projection.
%  7/6/99 - fixed tendency to re-define .ulongs if .clong set by user

global MAP_PROJECTION MAP_VAR_LIST

name={'Transverse Mercator','Sinusoidal','Gall-Peters','Hammer-Aitoff','Mollweide'};

pi180=pi/180;

switch optn,

  case 'name',

     X=name;

  case {'usage','set'}

     X=char({['     ''' varargin{1} ''''],...
              '     <,''lon<gitude>'',[min max]>',...
              '     <,''lat<itude>'',[min max]>',...
              '     <,''clo<ngitude>'',value>',...
              '     <,''rec<tbox>'', ( ''on'' | ''off'' )>'});

  case 'get',

     X=char([' Projection: ' MAP_PROJECTION.name '  (function: ' MAP_PROJECTION.routine ')'],...
            [' longitudes: ' num2str(MAP_VAR_LIST.ulongs) ' (centered at ' num2str(MAP_VAR_LIST.clong) ')'],...
            [' latitudes: ' num2str(MAP_VAR_LIST.ulats) ],...
            [' Rectangular border: ' MAP_VAR_LIST.rectbox ]); 

  case 'initialize',

    MAP_VAR_LIST=[];
    MAP_PROJECTION.name=varargin{1};
    switch MAP_PROJECTION.name,
      case name(1),
        MAP_VAR_LIST.ulongs=[-125 -122];
        MAP_VAR_LIST.ulats=[47 51];
      case name(2),
        MAP_VAR_LIST.ulongs=[-90 30];
        MAP_VAR_LIST.ulats=[-65 65];
      case name(3),
        MAP_VAR_LIST.ulongs=[30 390];
        MAP_VAR_LIST.ulats=[-65 65];
      case {name{4},name{5}}
        MAP_VAR_LIST.ulongs=[-300 60];
        MAP_VAR_LIST.ulats=[-90 90];
    end;
    MAP_VAR_LIST.clong=NaN;
    MAP_VAR_LIST.rectbox='off';
    k=2;longs_def=0;
    while k<length(varargin),   
      switch varargin{k}(1:3),
         case 'lon',
           MAP_VAR_LIST.ulongs=varargin{k+1};longs_def=1;
         case 'clo',
           MAP_VAR_LIST.clong=varargin{k+1};
         case 'lat',
           MAP_VAR_LIST.ulats=varargin{k+1};
         case 'rec',
           MAP_VAR_LIST.rectbox=varargin{k+1};
         otherwise
           disp(['Unknown option: ' varargin{k}]);
         end;
       k=k+2;
     end;
    if isnan(MAP_VAR_LIST.clong),  MAP_VAR_LIST.clong=mean(MAP_VAR_LIST.ulongs); 
    elseif ~longs_def, MAP_VAR_LIST.ulongs=MAP_VAR_LIST.clong+[-180 180];   end;
    
    MAP_VAR_LIST.clat=mean(MAP_VAR_LIST.ulats);
   
    % Get X/Y and (if rectboxs are desired) recompute lat/long limits.

    mu_util('xylimits');
    if strcmp(MAP_VAR_LIST.rectbox,'on'), mu_util('lllimits'); end;

  case 'll2xy',

    long=varargin{1};
    lat=varargin{2};

    % Clip out-of-range values (lat/long boxes)
    
    if  strcmp(MAP_VAR_LIST.rectbox,'off') & ~strcmp(varargin{4},'off'),
        [long,lat]=mu_util('clip',varargin{4},long,MAP_VAR_LIST.longs(1),long<MAP_VAR_LIST.longs(1),lat);
        [long,lat]=mu_util('clip',varargin{4},long,MAP_VAR_LIST.longs(2),long>MAP_VAR_LIST.longs(2),lat);
        [lat,long]=mu_util('clip',varargin{4},lat,MAP_VAR_LIST.lats(1),lat<MAP_VAR_LIST.lats(1),long);
        [lat,long]=mu_util('clip',varargin{4},lat,MAP_VAR_LIST.lats(2),lat>MAP_VAR_LIST.lats(2),long);
    end;

    switch MAP_PROJECTION.name,
      case name(1),
        B=cos(lat*pi180).*sin((long-MAP_VAR_LIST.clong)*pi180);
        X=atanh(B);
        Y=atan2(tan(lat*pi180),cos((long-MAP_VAR_LIST.clong)*pi180)) - MAP_VAR_LIST.clat*pi180;
      case name(2),
        Y=lat*pi180;
        X=pi180*(long-MAP_VAR_LIST.clong).*cos(Y)+MAP_VAR_LIST.clong*pi180;
      case name(3),
        X=pi180*(long-MAP_VAR_LIST.clong)*cos(45*pi180);
        Y=sin(lat*pi180)/cos(45*pi180);
      case name(4),
        z=sqrt((1+cos(lat*pi180).*cos((long-MAP_VAR_LIST.clong)*(pi180/2)))/2);
        X=2*cos(lat*pi180).*sin((long-MAP_VAR_LIST.clong)*(pi180/2))./z;
        Y=sin(lat*pi180)./z;
      case name(5),
        
        % Have to use interative scheme to get intermediate variable "theta".
        theta=(asin(lat/90)+lat*pi180)/2;
        dt=-(2*theta+sin(2*theta)-pi*sin(lat*pi180))./(1+cos(theta))/2;
        k=1;
        while any(abs(dt(:))>.001) & k<15,
          theta=theta+dt;
          dt=-(2*theta+sin(2*theta)-pi*sin(lat*pi180))./(1+cos(theta))/2;
          k=k+1;
 %% fprintf('%f %f\n',max(theta(:))/pi180,max(abs(dt(:))));
        end;
        if k==15, warning('Iterative coordinate conversion is not converging!'); end;
        theta=theta+dt;
        
        X=((long-MAP_VAR_LIST.clong).*cos(theta)+MAP_VAR_LIST.clong)/90;
        Y=sin(theta);
    end;

    % Clip out-of-range values (rectboxes)

    if strcmp(MAP_VAR_LIST.rectbox,'on')  & ~strcmp(varargin{4},'off'),
        [X,Y]=mu_util('clip',varargin{4},X,MAP_VAR_LIST.xlims(1),X<MAP_VAR_LIST.xlims(1),Y);
        [X,Y]=mu_util('clip',varargin{4},X,MAP_VAR_LIST.xlims(2),X>MAP_VAR_LIST.xlims(2),Y);
        [Y,X]=mu_util('clip',varargin{4},Y,MAP_VAR_LIST.ylims(1),Y<MAP_VAR_LIST.ylims(1),X);
        [Y,X]=mu_util('clip',varargin{4},Y,MAP_VAR_LIST.ylims(2),Y>MAP_VAR_LIST.ylims(2),X);
    end;

  case 'xy2ll',
    
    switch MAP_PROJECTION.name,
      case name(1),
        D=varargin{2}+MAP_VAR_LIST.clat*pi180;
        X=MAP_VAR_LIST.clong+atan2(sinh(varargin{1}),cos(D))/pi180;
        Y=asin(sin(D)./cosh(varargin{1}))/pi180;
      case name(2),
        Y=varargin{2}/pi180;
        X=MAP_VAR_LIST.clong+(varargin{1}-MAP_VAR_LIST.clong*pi180)./cos(varargin{2})/pi180;
      case name(3),
        X=varagin{1}/cos(45*pi180)/pi180 + MAP_VAR_LIST.clong;
        Y=asin(varargin{2}*cos(45*pi180))/pi180; 
      case name(4),
        z=sqrt(1-(varargin{1}/4).^2-(varargin{2}/2).^2);
        X=MAP_VAR_LIST.clong+2*atan2(z.*varargin{1},2*(2*z.^2-1))/pi180;
        Y=asin(varargin{2}.*z)/pi180;
      case name(5),
        theta=asin(varargin{2});
        Y=asin((2*theta+sin(2*theta))/pi)/pi180;
        X=(varargin{1}*90-MAP_VAR_LIST.clong)./cos(theta)+MAP_VAR_LIST.clong;
    end;
        
  case 'xgrid',
   
    [X,Y,vals,labI]=mu_util('xgrid',MAP_VAR_LIST.longs,MAP_VAR_LIST.lats,varargin{1},31,varargin{2});

  case 'ygrid',
   
    [X,Y,vals,labI]=mu_util('ygrid',MAP_VAR_LIST.lats,MAP_VAR_LIST.longs,varargin{1},31,varargin{2});

  case 'box',

    [X,Y]=mu_util('box',31);

end;


