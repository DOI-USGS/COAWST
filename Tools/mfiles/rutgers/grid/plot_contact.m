function plot_contact(G, S)

%
% PLOT_CONTACT:  Plots various ROMS Nested Grids Contact Points Figures
%
% h = plot_contact(G, S)
%
% This function plot contact points for each contact region.
%
% On Input:
%
%    G          Information grids structure (1 x Ngrids struct array)
%
%                 G(ng) = get_roms_grid ( char(Gnames(ng)) )
%
%    S          Nested Grids Structure (struct array)
%
%                 [S, G] = contact (Gnames, Cname, ...)
%

% svn $Id: plot_contact.m 738 2014-10-14 21:49:14Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Plot perimeters and boundary edged conectivity.
%--------------------------------------------------------------------------
  
for cr=1:S.Ncontact,
  dg = S.contact(cr).donor_grid;
  rg = S.contact(cr).receiver_grid;

  if (S.spherical),
    XminD = min(min(G(dg).lon_rho)); XmaxD = max(max(G(dg).lon_rho));
    YminD = min(min(G(dg).lat_rho)); YmaxD = max(max(G(dg).lat_rho));
    XminR = min(min(G(rg).lon_rho)); XmaxR = max(max(G(rg).lon_rho));
    YminR = min(min(G(rg).lat_rho)); YmaxR = max(max(G(rg).lat_rho));
    Xmin  = min(XminD, XminR);       Xmax  = max(XmaxD, XmaxR);
    Ymin  = min(YminD, YminR);       Ymax  = max(YmaxD, YmaxR);
  else
    XminD = min(min(G(dg).x_rho));   XmaxD = max(max(G(dg).x_rho));
    YminD = min(min(G(dg).y_rho));   YmaxD = max(max(G(dg).y_rho));
    XminR = min(min(G(rg).x_rho));   XmaxR = max(max(G(rg).x_rho));
    YminR = min(min(G(rg).y_rho));   YmaxR = max(max(G(rg).y_rho));
    Xmin  = min(XminD, XminR);       Xmax  = max(XmaxD, XmaxR);
    Ymin  = min(YminD, YminR);       Ymax  = max(YmaxD, YmaxR);
  end
  
  figure;
  ph1 = plot(S.grid(dg).perimeter.X_psi,                                ...
             S.grid(dg).perimeter.Y_psi, 'k:',                          ...
             S.grid(rg).perimeter.X_psi,                                ...
             S.grid(rg).perimeter.Y_psi, 'k:');
  axis([Xmin Xmax Ymin Ymax]);

  hold on;

  if (dg < rg),
    ph2 = plot(S.grid(dg).perimeter.X_psi,                              ...
               S.grid(dg).perimeter.Y_psi, 'r+',                        ...
               S.grid(rg).perimeter.X_psi,                              ...
               S.grid(rg).perimeter.Y_psi, 'bo');

    ph3 = plot(S.contact(cr).interior.Xdg,                              ...
               S.contact(cr).interior.Ydg,  'r^');
%   set(ph3,'MarkerSize',8);

    title(['PSI-Points Perimeters:',blanks(4),                          ...
           'Contact Region = ',num2str(cr),blanks(4),                   ...
           'Donor Grid = ',num2str(dg),blanks(2),'(red)',blanks(4),     ...
           'Receiver Grid = ',num2str(rg),blanks(2),'(blue)']);
    xlabel('triangles: RHO-points,  Other: PSI-points'); 
    set(gca,'fontsize',14,'fontweight','bold');
  else
    ph2 = plot(S.grid(dg).perimeter.X_psi,                              ...
               S.grid(dg).perimeter.Y_psi, 'b+',                        ...
               S.grid(rg).perimeter.X_psi,                              ...
               S.grid(rg).perimeter.Y_psi, 'ro');

    ph3 = plot(S.contact(cr).interior.Xdg,                              ...
               S.contact(cr).interior.Ydg,  'b^');
%   set(ph3,'MarkerSize',8);

    title(['PSI-Points Perimeters:',blanks(4),                          ...
           'Contact Region = ', num2str(cr), blanks(4),                 ...
           'Donor Grid = ',num2str(dg), blanks(2), '(blue)', blanks(4), ...
           'Receiver Grid = ', num2str(rg), blanks(2), '(red)']);
    xlabel('triangles: RHO-points,  Other: PSI-points'); 
    set(gca,'fontsize',14,'fontweight','bold');
  end      
      
  if (S.contact(cr).corners.okay),
    ph3 = plot(S.contact(cr).corners.Xdg,                               ...
               S.contact(cr).corners.Ydg, 'ks');
    set(ph3, 'MarkerSize', 12);
  end

  for ib=1:4,
    if (S.contact(cr).boundary(ib).okay),
      ph4 = plot(S.contact(cr).boundary(ib).Xdg,                         ...
                 S.contact(cr).boundary(ib).Ydg, 'ks');
      set(ph4, 'MarkerSize', 12);
    end
  end
  
  hold off;
end

%--------------------------------------------------------------------------
%  Plot contact points (Spherical or Cartesian Coordinates).
%--------------------------------------------------------------------------

for cr=1:S.Ncontact,
  dg = S.contact(cr).donor_grid;
  rg = S.contact(cr).receiver_grid;

  figure;

  if (S.spherical),
    XminD = min(min(G(dg).lon_rho)); XmaxD = max(max(G(dg).lon_rho));
    YminD = min(min(G(dg).lat_rho)); YmaxD = max(max(G(dg).lat_rho));
    XminR = min(min(G(rg).lon_rho)); XmaxR = max(max(G(rg).lon_rho));
    YminR = min(min(G(rg).lat_rho)); YmaxR = max(max(G(rg).lat_rho));
    Xmin  = min(XminD, XminR);       Xmax  = max(XmaxD, XmaxR);
    Ymin  = min(YminD, YminR);       Ymax  = max(YmaxD, YmaxR);
      
    ph1 = pcolor(G(dg).lon_rho, G(dg).lat_rho, G(dg).mask_rho);
    shading faceted; colormap(gray);
    axis([Xmin Xmax Ymin Ymax]);
    hold on

    ph2 = pcolor(G(rg).lon_rho, G(rg).lat_rho, G(rg).mask_rho);
    shading faceted; colormap(gray);
      
    set(gcf,'renderer','OpenGL');
    set(ph1,'EdgeColor',[0.6 0.6 0.6]);
    set(ph2,'EdgeColor',[0.6 0.6 0.6]);
    alpha(ph1, 0.2);
    alpha(ph2, 0.2);
  else
    XminD = min(min(G(dg).x_rho)); XmaxD = max(max(G(dg).x_rho));
    YminD = min(min(G(dg).y_rho)); YmaxD = max(max(G(dg).y_rho));
    XminR = min(min(G(rg).x_rho)); XmaxR = max(max(G(rg).x_rho));
    YminR = min(min(G(rg).y_rho)); YmaxR = max(max(G(rg).y_rho));
    Xmin  = min(XminD, XminR);     Xmax  = max(XmaxD, XmaxR);
    Ymin  = min(YminD, YminR);     Ymax  = max(YmaxD, YmaxR);
      
    ph1 = pcolorjw(G(dg).x_rho, G(dg).y_rho, G(dg).mask_rho);
    shading faceted; colormap(gray);
    axis([Xmin Xmax Ymin Ymax]);
    hold on;

    ph2 = pcolorjw(G(rg).x_rho, G(rg).y_rho, G(rg).mask_rho);
    shading faceted; colormap(gray);
      
    set(gcf,'renderer','OpenGL');
    set(ph1,'EdgeColor',[0.6 0.6 0.6]);
    set(ph2,'EdgeColor',[0.6 0.6 0.6]);
    alpha(ph1, 0.2);
    alpha(ph2, 0.2);
  end  

  xlabel(['PSI-points = circles,    RHO-points = squares,    ',         ...
          'U-points = diamods,    V-points = triangles']);
  
  if (dg < rg),
    ph3 = plot(S.grid(dg).perimeter.X_psi,                              ...
               S.grid(dg).perimeter.Y_psi, 'o');
    set(ph3,'MarkerEdgeColor','none','MarkerFaceColor','r');

    ph4 = plot(S.grid(rg).perimeter.X_psi,                              ...
               S.grid(rg).perimeter.Y_psi, 'o');
    set(ph4,'MarkerEdgeColor','none','MarkerFaceColor','b');

    ph5 = plot(S.contact(cr).point.Xrg_rho,                             ...
              S.contact(cr).point.Yrg_rho, 's');
    set(ph5,'MarkerEdgeColor','none','MarkerFaceColor','k');

    ph6 = plot(S.contact(cr).point.Xrg_u,                               ...
               S.contact(cr).point.Yrg_u, 'd');
    set(ph6,'MarkerEdgeColor','none','MarkerFaceColor','g',             ...
            'MarkerSize',8);

    ph7 = plot(S.contact(cr).point.Xrg_v,                               ...
              S.contact(cr).point.Yrg_v, 'v');
    set(ph7,'MarkerEdgeColor','none','MarkerFaceColor','m',             ...
            'MarkerSize',8);

    title(['Contact Points:  ',                                         ...
           'Contact Region = ',num2str(cr),',  ',                       ...
           'Donor Grid = ',num2str(dg),',  (red),  ',                   ...
           'Receiver Grid = ',num2str(rg),'  (blue),  '                 ...
           'RHO-mask']);
    set(gca,'fontsize',12,'fontweight','bold');
  
  else

    ph3 = plot(S.grid(dg).perimeter.X_psi,                              ...
               S.grid(dg).perimeter.Y_psi, 'o');
    set(ph3,'MarkerEdgeColor','none','MarkerFaceColor','b');

    ph4 = plot(S.grid(rg).perimeter.X_psi,                              ...
               S.grid(rg).perimeter.Y_psi, 'o');
    set(ph4,'MarkerEdgeColor','none','MarkerFaceColor','r');

    ph5 = plot(S.contact(cr).point.Xrg_rho,                             ...
               S.contact(cr).point.Yrg_rho, 's');
    set(ph5,'MarkerEdgeColor','none','MarkerFaceColor','k');

    ph6 = plot(S.contact(cr).point.Xrg_u,                               ...
               S.contact(cr).point.Yrg_u, 'd');
    set(ph6,'MarkerEdgeColor','none','MarkerFaceColor','g',             ...
            'MarkerSize',8);

    ph7 = plot(S.contact(cr).point.Xrg_v,                               ...
               S.contact(cr).point.Yrg_v, 'v');
    set(ph7,'MarkerEdgeColor','none','MarkerFaceColor','m',             ...
            'MarkerSize',8);

    title(['Contact Points:  ',                                         ...
           'Contact Region = ',num2str(cr),',  ',                       ...
           'Donor Grid = ',num2str(dg),',  (blue),  ',                  ...
           'Receiver Grid = ',num2str(rg),'  (red),  '                  ...
           'RHO-mask']);
    set(gca,'fontsize',12,'fontweight','bold');
  end
  hold off;
end

%--------------------------------------------------------------------------
%  Plot contact points in Fractional Coordinates.
%--------------------------------------------------------------------------

for cr=1:S.Ncontact,
  dg = S.contact(cr).donor_grid;
  rg = S.contact(cr).receiver_grid;

  if (S.contact(cr).refinement),
    if (S.grid(rg).refine_factor > 0),
      delta = 1.0/S.grid(rg).refine_factor;

      Io = min(S.contact(cr).corners.Idg);    % (Io,Jo) are the coarse grid
      Jo = min(S.contact(cr).corners.Jdg);    % coordinates (PSI-points)
    end                                       % used to extract refinement
  end                                         % grid
  
  figure;

  if (S.contact(cr).refinement),
    if (dg < rg),
      ph1 = pcolorjw(S.grid(dg).I_rho+0.5, S.grid(dg).J_rho+0.5,        ...
                     G(dg).mask_rho);
      shading faceted; colormap(gray);
      hold on;
      
%     Need to convert finer grid contact points (ROMS indices
%     Ir = -3:1:L+2, Jr = -3:1:M+2) to coarse grid fractional
%     coordinates.

      Irho = Io + (S.contact(cr).point.Irg_rho - 0.5).*delta;
      Ipsi = Io + (S.grid(rg).I_psi            - 1.0).*delta;
      Iu   = Io + (S.contact(cr).point.Irg_u   - 1.0).*delta;
      Iv   = Io + (S.contact(cr).point.Irg_v   - 0.5).*delta;
                                                
      Jrho = Jo + (S.contact(cr).point.Jrg_rho - 0.5).*delta;
      Jpsi = Jo + (S.grid(rg).J_psi            - 1.0).*delta;
      Ju   = Jo + (S.contact(cr).point.Jrg_u   - 0.5).*delta;
      Jv   = Jo + (S.contact(cr).point.Jrg_v   - 1.0).*delta;

      ph2 = plot(Ipsi, Jpsi, 'k-', Ipsi', Jpsi', 'k-');
    
      ph2 = plot(Ipsi(:,1)  , Jpsi(:,1)  , 'o',                         ...
                 Ipsi(end,:), Jpsi(end,:), 'o',                         ...
                 Ipsi(:,end), Jpsi(:,end), 'o',                         ...
                 Ipsi(1,:)  , Jpsi(1,:)  , 'o');
      set(ph2,'MarkerEdgeColor','none','MarkerFaceColor','b');
    else
      ph1 = pcolorjw(S.grid(rg).I_rho+0.5, S.grid(rg).J_rho+0.5,        ...
                     G(rg).mask_rho);
      shading faceted; colormap(gray);
      hold on;

      Irho = S.contact(cr).point.Irg_rho + 0.5;    % Need to convert
      Ipsi = S.grid(rg).I_psi;
      Iu   = S.contact(cr).point.Irg_u;            % coarse grid
      Iv   = S.contact(cr).point.Irg_v   + 0.5;    % contact points
                                                   % inside of finer
      Jrho = S.contact(cr).point.Jrg_rho + 0.5;    % grid (ROMS indices)
      Jpsi = S.grid(rg).J_psi;
      Ju   = S.contact(cr).point.Jrg_u   + 0.5;    % to coarse grid
      Jv   = S.contact(cr).point.Jrg_v;            % fractional
                                                   % coordinates
      ph2 = plot(Ipsi, Jpsi, 'k-', Ipsi', Jpsi', 'k-');
    
      ph2 = plot(Ipsi(:,1)  , Jpsi(:,1)  , 'o',                         ...
                 Ipsi(end,:), Jpsi(end,:), 'o',                         ...
                 Ipsi(:,end), Jpsi(:,end), 'o',                         ...
                 Ipsi(1,:)  , Jpsi(1,:)  , 'o');
      set(ph2,'MarkerEdgeColor','none','MarkerFaceColor','b');
    end
  else
    ph1 = pcolorjw(S.grid(dg).I_rho+0.5, S.grid(dg).J_rho+0.5,          ...
                   G(dg).mask_rho);
    shading faceted; colormap(gray);

    V = axis;
    V(1) = V(1) - 0.5;                        % insure drawing of all
    V(2) = V(2) + 0.5;                        % axis and meaningfull
    V(3) = V(3) - 0.5;                        % labels
    V(4) = V(4) + 0.5;
    axis (V);
    hold on;

    Irho = S.contact(cr).point.Idg_rho + 0.5;       % Convert to donor
    Jrho = S.contact(cr).point.Jdg_rho + 0.5;       % grid fractional
    Iu   = S.contact(cr).point.Idg_u;               % coordiates
    Ju   = S.contact(cr).point.Jdg_u   + 0.5;
    Iv   = S.contact(cr).point.Idg_v   + 0.5;
    Jv   = S.contact(cr).point.Jdg_v;

    ph2 = plot(S.grid(dg).I_psi(1,:),                                   ...
               S.grid(dg).J_psi(1,:), 'o',                              ...
               S.grid(dg).I_psi(:,1),                                   ...
               S.grid(dg).J_psi(:,1), 'o',                              ...
               S.grid(dg).I_psi(end,:),                                 ...
               S.grid(dg).J_psi(end,:), 'o',                            ...
               S.grid(dg).I_psi(:,end),                                 ...
               S.grid(dg).J_psi(:,end), 'o');
    if (dg < rg),
      set(ph2,'MarkerEdgeColor','none','MarkerFaceColor','r');
    else
      set(ph2,'MarkerEdgeColor','none','MarkerFaceColor','b');
    end      
  end
    
  set(gcf,'renderer','OpenGL');
  set(ph1,'EdgeColor',[0.6 0.6 0.6]);
  alpha(ph1, 0.2);

  xlabel(['PSI-points = circles,    RHO-points = squares,    ',         ...
          'U-points = diamods,    V-points = triangles']);
  
  ph3 = plot(Irho, Jrho, 's');
  set(ph3,'MarkerEdgeColor','none','MarkerFaceColor','k');

  ph4 = plot(Iu, Ju, 'd');
  set(ph4,'MarkerEdgeColor','none','MarkerFaceColor','g',               ...
          'MarkerSize',8);

  ph5 = plot(Iv, Jv, 'v');
  set(ph5,'MarkerEdgeColor','none','MarkerFaceColor','m',               ...
          'MarkerSize',8);

  title(['Contact Points:  ',                                           ...
         'Contact Region = ',num2str(cr),',  ',                         ...
         'Donor Grid = ',num2str(dg),',  ',                             ...
         'Receiver Grid = ',num2str(rg),'  ',                           ...
         'RHO-mask']);
  set(gca,'fontsize',12,'fontweight','bold');
  hold off;

% Plot contact points adjacent to refinement grid boundaries.
  
  if (S.contact(cr).refinement && dg < rg),

% RHO-contact points.

    figure;

    ph1 = pcolorjw(S.grid(dg).I_rho+0.5, S.grid(dg).J_rho+0.5,          ...
                   G(dg).mask_rho);
    shading faceted; colormap(gray);
    hold on;

    ph2 = plot(S.grid(dg).I_psi(:), S.grid(dg).J_psi(:), 'o');
    set(ph2,'MarkerEdgeColor','none','MarkerFaceColor','r');
    
    Ipsi = Io + (S.grid(rg).I_psi - 1.0).*delta;
    Jpsi = Jo + (S.grid(rg).J_psi - 1.0).*delta;
    ph3 = plot(Ipsi(:,1)  , Jpsi(:,1)  , 'o',                           ...
               Ipsi(end,:), Jpsi(end,:), 'o',                           ...
               Ipsi(:,end), Jpsi(:,end), 'o',                           ...
               Ipsi(1,:)  , Jpsi(1,:)  , 'o');
    set(ph3,'MarkerEdgeColor','none','MarkerFaceColor','b');
    
    X = Io + (S.refined(cr).Irg_rho - 0.5).*delta;
    Y = Jo + (S.refined(cr).Jrg_rho - 0.5).*delta;

    ph4 = plot(X, Y, 'k-', X', Y', 'k-');
    
    set(gcf,'renderer','OpenGL');
    set(ph1,'EdgeColor',[0.6 0.6 0.6]);
    alpha(ph1, 0.2);
    
    Npoints = length(S.contact(cr).point.Irg_rho);
    
    for n=1:Npoints,
      Xr = Io + (S.contact(cr).point.Irg_rho(n) - 0.5).*delta;
      Yr = Jo + (S.contact(cr).point.Jrg_rho(n) - 0.5).*delta;
      plot(Xr, Yr, '');
      strdonor = [num2str(S.contact(cr).point.Idg_rho(n)), ',',         ...
                  num2str(S.contact(cr).point.Jdg_rho(n))];
      text(Xr, Yr, strdonor, 'HorizontalAlignment','center');
    end

    title(['RHO-Contact Points:  ',                                     ...
           'Contact Region = ',num2str(cr),',  ',                       ...
           'Donor Grid = ',num2str(dg),',  ',                           ...
           'Receiver Grid = ',num2str(rg),',  ',                        ...
           'RHO-mask']);
    xlabel('Circles: PSI-points');
    set(gca,'fontsize',12,'fontweight','bold');
    hold off;
    
% U-contact points.
     
    figure;

    ph1 = pcolorjw(S.grid(dg).I_u, S.grid(dg).J_u+0.5,                  ...
                   G(dg).mask_u); 
    shading faceted; colormap(gray);
    hold on;

    ph2 = plot(S.grid(dg).I_psi(:)-0.5, S.grid(dg).J_psi(:), 'o');
    set(ph2,'MarkerEdgeColor','none','MarkerFaceColor','r');

    Ipsi = Io + (S.grid(rg).I_psi - 1.0).*delta;
    Jpsi = Jo + (S.grid(rg).J_psi - 1.0).*delta;

    ph3 = plot(Ipsi(:,1)  , Jpsi(:,1)  , 'o',                           ...
               Ipsi(end,:), Jpsi(end,:), 'o',                           ...
               Ipsi(:,end), Jpsi(:,end), 'o',                           ...
               Ipsi(1,:)  , Jpsi(1,:)  , 'o');
    set(ph3,'MarkerEdgeColor','none','MarkerFaceColor','b');
    
    X = Io + (S.refined(cr).Irg_u - 1.0).*delta;
    Y = Jo + (S.refined(cr).Jrg_u - 0.5).*delta;

    ph4 = plot(X, Y, 'k-', X', Y', 'k-');
    
    set(gcf,'renderer','OpenGL');
    set(ph1,'EdgeColor',[0.6 0.6 0.6]);
    alpha(ph1, 0.2);
    
    Npoints = length(S.contact(cr).point.Irg_u);
    
    for n=1:Npoints,
      Xu = Io + (S.contact(cr).point.Irg_u(n) - 1.0).*delta;
      Yu = Jo + (S.contact(cr).point.Jrg_u(n) - 0.5).*delta;
      plot(Xu, Yu, '');
      strdonor = [num2str(S.contact(cr).point.Idg_u(n)), ',',           ...
                  num2str(S.contact(cr).point.Jdg_u(n))];
      text(Xu, Yu, strdonor, 'HorizontalAlignment','center');
    end

    title(['U-Contact Points:  ',                                       ...
           'Contact Region = ',num2str(cr),',  ',                       ...
           'Donor Grid = ',num2str(dg),',  ',                           ...
           'Receiver Grid = ',num2str(rg),',  ',                        ...
           'U-mask']);
    xlabel('Circles: PSI-points');
    set(gca,'fontsize',12,'fontweight','bold');
    hold off
    
% V-contact points.
     
    figure;

    ph1 = pcolorjw(S.grid(dg).I_v+0.5, S.grid(dg).J_v,                  ...
                   G(dg).mask_v );
    shading faceted; colormap(gray);
    hold on;

    ph2 = plot(S.grid(dg).I_psi(:), S.grid(dg).J_psi(:)-0.5, 'o');
    set(ph2,'MarkerEdgeColor','none','MarkerFaceColor','r');
    
    Ipsi = Io + (S.grid(rg).I_psi - 1.0).*delta;
    Jpsi = Jo + (S.grid(rg).J_psi - 1.0).*delta;

    ph3 = plot(Ipsi(:,1)  , Jpsi(:,1)  , 'o',                           ...
               Ipsi(end,:), Jpsi(end,:), 'o',                           ...
               Ipsi(:,end), Jpsi(:,end), 'o',                           ...
               Ipsi(1,:)  , Jpsi(1,:)  , 'o');
    set(ph3,'MarkerEdgeColor','none','MarkerFaceColor','b');
    
    X = Io + (S.refined(cr).Irg_v - 0.5).*delta;
    Y = Jo + (S.refined(cr).Jrg_v - 1.0).*delta;

    ph4 = plot(X, Y, 'k-', X', Y', 'k-');
    
    set(gcf,'renderer','OpenGL');
    set(ph1,'EdgeColor',[0.6 0.6 0.6]);
    alpha(ph1, 0.2);
    
    Npoints = length(S.contact(cr).point.Irg_v);
    
    for n=1:Npoints,
      Xv = Io + (S.contact(cr).point.Irg_v(n) - 0.5).*delta;
      Yv = Jo + (S.contact(cr).point.Jrg_v(n) - 1.0).*delta;
      plot(Xv, Yv, '');
      strdonor = [num2str(S.contact(cr).point.Idg_v(n)), ',',           ...
                  num2str(S.contact(cr).point.Jdg_v(n))];
      text(Xv, Yv, strdonor, 'HorizontalAlignment','center');
    end

    title(['V-Contact Points:  ',                                       ...
           'Contact Region = ',num2str(cr),',  ',                       ...
           'Donor Grid = ',num2str(dg),',  ',                           ...
           'Receiver Grid = ',num2str(rg),',  ',                        ...
           'V-mask']);
    xlabel('Circles: PSI-points');
    set(gca,'fontsize',12,'fontweight','bold');
  end

end

%--------------------------------------------------------------------------
%  Plot contact points associated with lateral boundary conditions.
%--------------------------------------------------------------------------

for cr=1:S.Ncontact,
  dg = S.contact(cr).donor_grid;
  rg = S.contact(cr).receiver_grid;

% Contact points on RHO-boundary.

  ONr = S.contact(cr).point.boundary_rho;

  if (any(ONr)),
    figure;
    if (S.spherical)
      ph1 = pcolorjw(G(dg).lon_rho, G(dg).lat_rho, G(dg).mask_rho);
      dx  = G(dg).lon_rho(2) - G(dg).lon_rho(1);
      dy  = G(dg).lat_rho(2) - G(dg).lat_rho(1);
    else
      ph1 = pcolorjw(G(dg).x_rho, G(dg).y_rho, G(dg).mask_rho);
      dx  = G(dg).x_rho(2) - G(dg).x_rho(1);
      dy  = G(dg).y_rho(2) - G(dg).y_rho(1);
    end
    shading faceted; colormap(gray);
    hold on
    
    V = axis;
    if (~S.contact(cr).refinement),
      V(1) = V(1) - 0.5*dx;                   % insure drawing of all
      V(2) = V(2) + 0.5*dx;                   % axis and meaningfull
      V(3) = V(3) - 0.5*dy;                   % labels
      V(4) = V(4) + 0.5*dy;
      axis (V);    
    end
    
    ph2 = plot(S.contact(cr).point.Xrg_rho,                             ...
               S.contact(cr).point.Yrg_rho, 'ro',                       ...
               S.grid(rg).perimeter.X_psi,                              ...
               S.grid(rg).perimeter.Y_psi, 'k:');

    ph3 = plot(S.grid(rg).perimeter.X_psi,                              ...
               S.grid(rg).perimeter.Y_psi, 's');
    set(ph3,'MarkerEdgeColor','none','MarkerFaceColor','k');

    ph4 = plot(S.contact(cr).point.Xrg_rho(ONr),                        ...
               S.contact(cr).point.Yrg_rho(ONr), 'o');
    set(ph4,'MarkerEdgeColor','none','MarkerFaceColor','r');

    set(gcf,'renderer','OpenGL');
    set(ph1,'EdgeColor',[0.6 0.6 0.6]);
    alpha(ph1, 0.2);
    
    title(['RHO-Contact Points on boundary:  ',                         ...
           'Contact Region = ',num2str(cr),',  ',                       ...
           'Donor Grid = ',num2str(dg),',  ',                           ...
           'Receiver Grid = ',num2str(rg),',  ',                        ...
           'RHO-mask']);
    xlabel('Squares: PSI-points, Circles: RHO-points');
  end

% Contact points on U-boundary.

  ONu = S.contact(cr).point.boundary_u;

  if (any(ONu)),
    figure;
    if (S.spherical)
      ph1 = pcolorjw(G(dg).lon_u, G(dg).lat_u, G(dg).mask_u);
      dx  = G(dg).lon_u(2) - G(dg).lon_u(1);
      dy  = G(dg).lat_u(2) - G(dg).lat_u(1);
    else
      ph1 = pcolorjw(G(dg).x_u, G(dg).y_u, G(dg).mask_u);
      dx  = G(dg).x_u(2) - G(dg).x_u(1);
      dy  = G(dg).y_u(2) - G(dg).y_u(1);
    end
    shading faceted; colormap(gray);
    hold on

    V = axis;
    if (~S.contact(cr).refinement),
      V(1) = V(1) - 0.5*dx;                   % insure drawing of all
      V(2) = V(2) + 0.5*dx;                   % axis and meaningfull
      V(3) = V(3) - 0.5*dy;                   % labels
      V(4) = V(4) + 0.5*dy;
      axis (V);    
    end

    ph2 = plot(S.contact(cr).point.Xrg_u,                               ...
               S.contact(cr).point.Yrg_u, 'ro',                         ...
               S.grid(rg).perimeter.X_psi,                              ...
               S.grid(rg).perimeter.Y_psi, 'k:');

    ph3 = plot(S.grid(rg).perimeter.X_psi,                              ...
               S.grid(rg).perimeter.Y_psi, 's');
    set(ph3,'MarkerEdgeColor','none','MarkerFaceColor','k');

    ph4 = plot(S.contact(cr).point.Xrg_u(ONu),                          ...
               S.contact(cr).point.Yrg_u(ONu), 'o');
    set(ph4,'MarkerEdgeColor','none','MarkerFaceColor','r');

    set(gcf,'renderer','OpenGL');
    set(ph1,'EdgeColor',[0.6 0.6 0.6]);
    alpha(ph1, 0.2);
    
    title(['U-Contact Points on boundary:  ',                           ...
           'Contact Region = ',num2str(cr),',  ',                       ...
           'Donor Grid = ',num2str(dg),',  ',                           ...
           'Receiver Grid = ',num2str(rg),',  ',                        ...
           'U-mask']);
    xlabel('Squares: PSI-points, Circles: U-points');
  end

% Contact points on V-boundary.

  ONv = S.contact(cr).point.boundary_v;

  if (any(ONv)),
    figure;
    if (S.spherical)
      ph1 = pcolorjw(G(dg).lon_v, G(dg).lat_v, G(dg).mask_v);
      dx  = G(dg).lon_v(2) - G(dg).lon_v(1);
      dy  = G(dg).lat_v(2) - G(dg).lat_v(1);
    else
      ph1 = pcolorjw(G(dg).x_v, G(dg).y_v, G(dg).mask_v);
      dx  = G(dg).x_v(2) - G(dg).x_v(1);
      dy  = G(dg).y_v(2) - G(dg).y_v(1);
    end
    shading faceted; colormap(gray);
    hold on

    V = axis;
    if (~S.contact(cr).refinement),
      V(1) = V(1) - 0.5*dx;                   % insure drawing of all
      V(2) = V(2) + 0.5*dx;                   % axis and meaningfull
      V(3) = V(3) - 0.5*dy;                   % labels
      V(4) = V(4) + 0.5*dy;
      axis (V);    
    end

    ph2 = plot(S.contact(cr).point.Xrg_v,                               ...
               S.contact(cr).point.Yrg_v, 'ro',                         ...
               S.grid(rg).perimeter.X_psi,                              ...
               S.grid(rg).perimeter.Y_psi, 'k:');

    ph3 = plot(S.grid(rg).perimeter.X_psi,                              ...
               S.grid(rg).perimeter.Y_psi, 's');
    set(ph3,'MarkerEdgeColor','none','MarkerFaceColor','k');

    ph4 = plot(S.contact(cr).point.Xrg_v(ONv),                          ...
               S.contact(cr).point.Yrg_v(ONv), 'o');
    set(ph4,'MarkerEdgeColor','none','MarkerFaceColor','r');

    set(gcf,'renderer','OpenGL');
    set(ph1,'EdgeColor',[0.6 0.6 0.6]);
    alpha(ph1, 0.2);
    
    title(['V-Contact Points on boundary:  ',                           ...
           'Contact Region = ',num2str(cr),',  ',                       ...
           'Donor Grid = ',num2str(dg),',  ',                           ...
           'Receiver Grid = ',num2str(rg),',  ',                        ...
           'V-mask']);
    xlabel('Squares: PSI-points, Circles: V-points');
  end
end

return
