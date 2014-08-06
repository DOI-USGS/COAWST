function [z,s,C]=scoord(h, x, y, Vtransform, Vstretching, ...
                        theta_s, theta_b, hc,             ...
                        N, kgrid, column, index, plt, Zzoom);
%
% SCOORD:  Compute and plot ROMS vertical stretched coordinates
%
% [z,s,C]=scoord(h, x, y, Vtransform, Vstretching, theta_s, theta_b, ...
%                hc, N, kgrid, column, index, plt, Zzoom)
%
% Given a batymetry (h) and terrain-following stretching parameters,
% this function computes the depths of RHO- or W-points for a vertical
% grid section along columns (ETA-axis) or rows (XI-axis). Check the
% following link for details:
%
%    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
%
% On Input:
%
%    h             Bottom depth 2D array, h(1:Lp,1:Mp), m, positive
%    x             X-coordinate 2D array, x(1:Lp,1:Mp), m or degrees_east
%    y             Y-coordinate 2D array, y(1:Lp,1:Mp), m or degrees_north
%    Vtransform    Vertical transformation equation:
%                    Vtransform = 1,   original transformation
%
%                      z(x,y,s,t)=Zo(x,y,s)+zeta(x,y,t)*[1+Zo(x,y,s)/h(x,y)]
%
%                      Zo(x,y,s)=hc*s+[h(x,y)-hc]*C(s)
%
%                    Vtransform = 2,   new transformation
%
%                      z(x,y,s,t)=zeta(x,y,t)+[zeta(x,y,t)+h(x,y)]*Zo(x,y,s)
%
%                       Zo(x,y,s)=[hc*s(k)+h(x,y)*C(k)]/[hc+h(x,y)]
%    Vstretching   Vertical stretching function:
%                    Vstretching = 1,  original (Song and Haidvogel, 1994)
%                    Vstretching = 2,  A. Shchepetkin (UCLA-ROMS, 2005)
%                    Vstretching = 3,  R. Geyer BBL refinement
%                    Vstretching = 4,  A. Shchepetkin (UCLA-ROMS, 2010)
%    theta_s       S-coordinate surface control parameter (scalar)
%    theta_b       S-coordinate bottom control parameter (scalar)
%    hc            Width (m) of surface or bottom boundary layer in which
%                    higher vertical resolution is required during
%                    stretching (scalar)
%    N             Number of vertical levels (scalar)
%    kgrid         Depth grid type logical switch:
%                    kgrid = 0,        depths of RHO-points
%                    kgrid = 1,        depths of W-points
%    column        Grid direction logical switch:
%                    column = 1,       column section
%                    column = 0,       row section
%    index         Column or row to compute (scalar)
%                    if column = 1,    then   1 <= index <= Lp
%                    if column = 0,    then   1 <= index <= Mp
%    plt           Switch to plot scoordinate (scalar):
%                    plt = 0,          do not plot
%                    plt = 1,          plot
%                    plt = 2,          plot 2 pannels with zoom
%    Zzoom         If plt=2, maximum depth of the zoom in upper pannel
%
% On Output:
%
%    z             Depths (m) of RHO- or W-points (matrix)
%    s             S-coordinate independent variable, [-1 <= s <= 0] at
%                    vertical RHO- or W-points (vector)
%    C             Nondimensional, monotonic, vertical stretching function,
%                    C(s), 1D array, [-1 <= C(s) <= 0]
%

% svn $Id: scoord.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

z=[];
s=[];
C=[];

%----------------------------------------------------------------------------
%  Set several parameters.
%----------------------------------------------------------------------------

if (hc > min(min(h)) & Vtransform == 1),
  disp(' ');
  disp([setstr(7),'*** Error:  SCOORD - critical depth exceeds minimum' ...
       ' bathymetry value.',setstr(7)]);
  disp([setstr(7),'                     Vtranform = ',num2str(Vtransform),setstr(7)]);
  disp([setstr(7),'                     hc        = ',num2str(hc),setstr(7)]);
  disp([setstr(7),'                     hmax      = ',num2str(min(min(h))), ...
        setstr(7)]);
  disp(' ');
  return
end,

if (Vtransform < 1 | Vtransform > 2),
  disp(' ');
  disp([setstr(7),'*** Error:  SCOORD - Illegal parameter Vtransform = ' ...
        num2str(Vtransfrom), setstr(7)]);
  return
end,

if (Vstretching < 1 | Vstretching > 4),
  disp(' ');
  disp([setstr(7),'*** Error:  SCOORD - Illegal parameter Vstretching = ' ...
        num2str(Vstretching), setstr(7)]);
  return
end,

[Lp Mp]=size(h);
hmin=min(min(h));
hmax=max(max(h));
havg=0.5*(hmax+hmin);

%----------------------------------------------------------------------------
% Test input to see if it's in an acceptable form.
%----------------------------------------------------------------------------

if (nargin < 12),
  disp(' ');
  disp([setstr(7),'*** Error:  SCOORD - too few arguments.',setstr(7)]);
  disp([setstr(7),'                     number of supplied arguments: ',...
       num2str(nargin),setstr(7)]);
  disp([setstr(7),'                     number of required arguments: 8',...
       setstr(7)]);
  disp(' ');
  return
end,

if (column),

  if (index < 1 | index > Lp),
    disp(' ');
    disp([setstr(7),'*** Error:  SCOORD - illegal column index.',setstr(7)]);
    disp([setstr(7),'                     valid range:  1 <= index <= ',...
         num2str(Lp),setstr(7)]);
    disp(' ');
    return
  end,

else,

  if (index < 1 | index > Mp),
    disp(' ');
    disp([setstr(7),'*** Error:  SCOORD - illegal row index.',setstr(7)]);
    disp([setstr(7),'                     valid range:  1 <= index <= ',...
         num2str(Mp),setstr(7)]);
    disp(' ');
    return
  end,

end,

%----------------------------------------------------------------------------
% Compute vertical stretching function, C(k):
%----------------------------------------------------------------------------

report=0;

[s,C]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid, report);

if (kgrid == 1),
  Nlev=N+1;
else,
  Nlev=N;
end,

if (Vtransform == 1),

  for k=Nlev:-1:1,
    zhc(k)=hc*s(k);
    z1 (k)=zhc(k)+(hmin-hc)*C(k);
    z2 (k)=zhc(k)+(havg-hc)*C(k);
    z3 (k)=zhc(k)+(hmax-hc)*C(k);
  end,

elseif (Vtransform == 2),

  for k=Nlev:-1:1,
    if (hc > hmax),
      zhc(k)=hmax*(hc*s(k)+hmax*C(k))/(hc+hmax);
    else
      zhc(k)=0.5*min(hc,hmax)*(s(k)+C(k));
    end,
    z1 (k)=hmin*(hc*s(k)+hmin*C(k))/(hc+hmin);
    z2 (k)=havg*(hc*s(k)+havg*C(k))/(hc+havg);
    z3 (k)=hmax*(hc*s(k)+hmax*C(k))/(hc+hmax);
  end,

end,

report=1;

if (report),
  disp(' ');
  if (Vtransform == 1),
    disp(['Vtransform  = ',num2str(Vtransform), '   original ROMS']);
  elseif (Vtransform == 2),
    disp(['Vtransform  = ',num2str(Vtransform), '   ROMS-UCLA']);
  end,
  if (Vstretching == 1),
    disp(['Vstretching = ',num2str(Vstretching), '   Song and Haidvogel (1994)']);
  elseif (Vstretching == 2),
    disp(['Vstretching = ',num2str(Vstretching), '   Shchepetkin (2005)']);
  elseif (Vstretching == 3),
    disp(['Vstretching = ',num2str(Vstretching), '   Geyer (2009), BBL']);
  elseif (Vstretching == 4),
    disp(['Vstretching = ',num2str(Vstretching), '   Shchepetkin (2010)']);
  end,
  if (kgrid == 1)
    disp(['   kgrid    = ',num2str(kgrid), '   at vertical W-points']);
  else,
    disp(['   kgrid    = ',num2str(kgrid), '   at vertical RHO-points']);
  end,
  disp(['   theta_s  = ',num2str(theta_s)]);
  disp(['   theta_b  = ',num2str(theta_b)]);
  disp(['   hc       = ',num2str(hc)]);

  disp(' ');
  disp(' S-coordinate curves: ')
  disp(' ');
  disp([' level     S-coord    Cs-Curve   Z  at hmin      ', ...
        ' at hc    half way     at hmax']);
  disp(' ');

  if (kgrid == 1),
    for k=Nlev:-1:1,
      disp(['   ', ...
            sprintf('%3i',k-1      ), ' ', ...
            sprintf('%12.7f',s(k)  ), ...
            sprintf('%12.7f',C(k)  ), ...
            sprintf('%12.3f',z1(k) ), ...
            sprintf('%12.3f',zhc(k)), ...
            sprintf('%12.3f',z2(k) ), ...
            sprintf('%12.3f',z3(k) )]);
    end,
  else
    for k=Nlev:-1:1,
      disp(['   ', ...
            sprintf('%3i',k        ), ' ', ...
            sprintf('%12.7f',s(k)  ), ...
            sprintf('%12.7f',C(k)  ), ...
            sprintf('%12.3f',z1(k) ), ...
            sprintf('%12.3f',zhc(k)), ...
            sprintf('%12.3f',z2(k) ), ...
            sprintf('%12.3f',z3(k) )]);
    end,
  end,
  disp(' ');

end,

%============================================================================
% Compute depths at requested grid section.  Assume zero free-surface.
%============================================================================

zeta=zeros(size(h));

%----------------------------------------------------------------------------
% Column section: section along ETA-axis.
%----------------------------------------------------------------------------

if (column),

  if (Vtransform == 1),

    z=zeros(Mp,Nlev);
    for k=1:Nlev,
      z0=hc.*(s(k)-C(k))+h(index,:)*C(k);
      z(:,k)=z0+zeta(index,:).*(1.0+z0/h(index,:));
    end,

  elseif (Vtransform == 2),

    z=zeros(Mp,Nlev);
    for k=1:Nlev,
      z0=(hc.*s(k)+C(k).*h(index,:))./(h(index,:)+hc);
      z(:,k)=zeta(index,:)+(zeta(index,:)+h(index,:)).*z0;
    end,

  end,

%----------------------------------------------------------------------------
% Row section: section along XI-axis.
%----------------------------------------------------------------------------

else,

  if (Vtransform == 1),

    z=zeros(Lp,Nlev);
    for k=1:Nlev,
      z0=hc.*(s(k)-C(k))+h(:,index)*C(k);
      z(:,k)=z0+zeta(:,index).*(1.0+z0./h(:,index));
    end,

  elseif (Vtransform == 2),

    z=zeros(Lp,Nlev);
    for k=1:Nlev,
      z0=(hc.*s(k)+C(k).*h(:,index))./(h(:,index)+hc);
      z(:,k)=zeta(:,index)+(zeta(:,index)+h(:,index)).*z0;
    end,

  end,

end,

%============================================================================
% Plot grid section.
%============================================================================

if nargin < 13,
  plt = 1;
end

if (plt > 0),

  figure;

  if (column),

    set(gcf,'Units','Normalized',...
        'Position',[0.2 0.1 0.6 0.8],...
        'PaperOrientation', 'landscape', ...
        'PaperUnits','Normalized',...
        'PaperPosition',[0.2 0.1 0.6 0.8]);

    eta=y(index,:);
    eta2=[eta(1) eta eta(Mp)];
    hs=-h(index,:);
    zmin=min(hs);
    hs=[zmin hs zmin];

    if (plt == 2),
      h1=subplot(2,1,1);
      p1=get(h1,'pos');
      p1(2)=p1(2)+0.1;
      p1(4)=p1(4)-0.1;
      set (h1,'pos',p1);
    end,

    hold off;
    fill(eta2,hs,[0.6 0.7 0.6]);
    hold on;
    han1=plot(eta',z);
    set(han1,'color', [0.5 0.5 0.5]);
    if (plt == 2),
      set(gca,'xlim',[-Inf Inf],'ylim',[-abs(Zzoom) 0]);
      ylabel('depth  (m)');
    else
      set(gca,'xlim',[-Inf Inf],'ylim',[zmin 0]);
    end,

    if (kgrid == 0),
      title(['Grid (\rho-points) Section at  \xi = ',num2str(index)]);
    else

      title(['Grid (W-points) Section at  \xi = ',num2str(index)]);
    end,

    if (plt == 2),
      h2=subplot(2,1,2);
      p2=get(h2,'pos');
      p2(4)=p2(4)+0.2;
      set (h2,'pos',p2);

      hold off;
      fill(eta2,hs,[0.6 0.7 0.6]);
      hold on;
      han2=plot(eta',z);
      set(han2,'color', [0.5 0.5 0.5]);
      set(gca,'xlim',[-Inf Inf],'ylim',[zmin 0]);
    end,

    xlabel({['Vcoord = ', num2str(Vtransform),',',num2str(Vstretching), ...
             '   \theta_s = ' num2str(theta_s), ...
             '   \theta_b = ' num2str(theta_b), ...
             '   hc  = ' num2str(hc), ...
             '   N = ' num2str(N)],['y-axis']});
    ylabel('depth  (m)');

  else,

    set(gcf,'Units','Normalized',...
        'Position',[0.2 0.1 0.6 0.8],...
        'PaperOrientation', 'landscape', ...
        'PaperUnits','Normalized',...
        'PaperPosition',[0.2 0.1 0.6 0.8]);

    xi=x(:,index)';
    xi2=[xi(1) xi xi(Lp)];
    hs=-h(:,index)';
    zmin=min(hs);
    hs=[zmin hs zmin];

    if (plt == 2),
      h1=subplot(2,1,1);
      p1=get(h1,'pos');
      p1(2)=p1(2)+0.1;
      p1(4)=p1(4)-0.1;
      set (h1,'pos',p1);
    end,

    hold off;
    fill(xi2,hs,[0.6 0.7 0.6]);
    hold on;
    han1=plot(xi,z);
    set(han1,'color', [0.5 0.5 0.5]);
    if (plt == 2),
      set(gca,'xlim',[-Inf Inf],'ylim',[-abs(Zzoom) 0]);
      ylabel('depth  (m)');
    else
      set(gca,'xlim',[-Inf Inf],'ylim',[zmin 0]);
    end,

    if (kgrid == 0),
      title(['Grid (\rho-points) Section at  \eta = ',num2str(index)]);
    else
      title(['Grid (W-points) Section at  \eta = ',num2str(index)]);
    end,

    if (plt == 2),
      h2=subplot(2,1,2);
      p2=get(h2,'pos');
      p2(4)=p2(4)+0.2;
      set (h2,'pos',p2);

      hold off;
      fill(xi2,hs,[0.6 0.7 0.6]);
      hold on;
      han2=plot(xi,z);
      set(han2,'color', [0.5 0.5 0.5]);
      set(gca,'xlim',[-Inf Inf],'ylim',[zmin 0]);
    end,

    xlabel({['Vcoord = ', num2str(Vtransform),',',num2str(Vstretching), ...
             '   \theta_s = ' num2str(theta_s), ...
             '   \theta_b = ' num2str(theta_b), ...
             '   hc  = ' num2str(hc), ...
             '   N = ' num2str(N)],['x-axis']});
    ylabel('depth  (m)');

  end,

end,

return


