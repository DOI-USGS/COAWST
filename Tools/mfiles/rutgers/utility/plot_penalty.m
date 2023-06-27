function [J]=plot_penalty(ncname, varargin)

%
% PLOT_PENALTY:  Plots ROMS 4DVAR penalty (cost) function
%
% [J]=plot_penalty(ncname, oname)
%
% This script plots 4DVAR penalty function or cost function.
% It uses the global attribute "Algorithm" to determine the
% type of penalty/cost functional to process.  If the optional
% file is provided, it also draws J for comparison.
%
% On Input:
%
%    ncname        4D-Var DAV/MOD output NetCDF filename (string)
%
%    oname         Optional NetCDF filename for comparisons (string)
%
% On Output:
%
%    J             Penalty/Cost function structure (struct)
%
  
% svn $Id: plot_penalty.m 1156 2023-02-18 01:44:37Z arango $
%=========================================================================%
%  Copyright (c) 2002-2023 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

switch numel(varargin)
  case 0
    got_overlay = false;
  case 1
    got_overlay = true;
    oname = varargin{1};
end
  
% Inquire NetCDF file.

I = nc_inq(ncname);

index = strcmp({I.Attributes.Name},'Algorithm');
if (any(index))
  Algorithm = I.Attributes(index).Value;
else
  disp(['Global Attribute ''Algorithm'' not found in: ',ncname']);
  return
end

index  = strcmp({I.Dimensions.Name},'Nouter');
Nouter = I.Dimensions(index).Length;

index  = strcmp({I.Dimensions.Name},'Ninner');
Ninner = I.Dimensions(index).Length;

% Initialize.

J.ncname    = ncname;
J.dir       = pwd;
J.algorithm = Algorithm;
J.Nouter    = Nouter;
J.Ninner    = Ninner;

%--------------------------------------------------------------------------
%  Read in cost function data.
%--------------------------------------------------------------------------

switch Algorithm
 
 case {'I4DVAR', 'IS4DVAR', 'SPLIT_I4DVAR'}

    J.Nobs = length(nonzeros(nc_read(ncname, 'obs_scale')));
    J.o    = nc_read(ncname, 'TLcost_function');
    J.b    = nc_read(ncname, 'back_function');
    J.NL   = nc_read(ncname, 'NLcost_function');
    J.act  = J.b + J.o;

    Niter  = length(J.act);
    
 case {'RBL4DVAR', 'SPLIT_RBL4DVAR', 'W4DPSAS'}

    half=0.5;   % by definition of penalty functional
							     
    J.Nobs = length(nonzeros(nc_read(ncname, 'obs_scale')));
    Jobs   = nc_read(ncname, 'Jobs', [], NaN);
    Jb     = nc_read(ncname, 'Jb', [], NaN);
    Jmod   = nc_read(ncname, 'Jmod', [], NaN);
    Jact   = nc_read(ncname, 'Jact', [], NaN);
    Jopt   = nc_read(ncname, 'Jopt', [], NaN);
    J.NL   = nc_read(ncname, 'NL_fDataPenalty');
    if (Nouter == 1)
      Jobs(1)=[]; Jb(1)=[]; Jmod(1)=[]; Jact(1)=[]; Jopt(1)=[];
    else
      Jobs(1,:)=[]; Jb(1,:)=[]; Jmod(1,:)=[]; Jact(1,:)=[]; Jopt(1,:)=[];
    end
    J.o   = Jobs.*half;
    J.b   = Jb.*half;
    J.mod = Jmod.*half;
    J.act = Jact.*half;
    J.opt = Jopt.*half;

    Niter = J.Ninner*J.Nouter;

 case {'R4DVAR', 'SPLIT_R4DVAR', 'W4DVAR'}
  
    half=0.5;   % by definition of penalty functional

    J.Nobs = length(nonzeros(nc_read(ncname, 'obs_scale')));
    Jobs   = nc_read(ncname, 'Jobs', [], NaN);
    Jb     = nc_read(ncname, 'Jb', [], NaN);
    Jmod   = nc_read(ncname, 'Jmod', [], NaN);
    Jact   = nc_read(ncname, 'Jact', [], NaN);
    Jopt   = nc_read(ncname, 'Jopt', [], NaN);
    J.RP   = nc_read(ncname, 'RP_fDataPenalty');
    if (Nouter == 1)
      Jobs(1)=[]; Jb(1)=[]; Jmod(1)=[]; Jact(1)=[]; Jopt(1)=[];
    else
      Jobs(1,:)=[]; Jb(1,:)=[]; Jmod(1,:)=[]; Jact(1,:)=[]; Jopt(1,:)=[];
    end
    
    J.o   = Jobs.*half;
    J.b   = Jb.*half;
    J.act = Jact.*half;
    
    Niter = J.Ninner*J.Nouter;
end

%--------------------------------------------------------------------------
% Read in overlay cost function data, if any. It is usually use for
% comparisons.
%--------------------------------------------------------------------------

if (got_overlay)
 I2 = nc_inq(oname);

 index = strcmp({I2.Attributes.Name},'Algorithm');
 if (any(index))
   Algorithm2 = I2.Attributes(index).Value;
 else
  disp(['Global Attribute ''Algorithm'' not found in: ',oname']);
  return
 end

  index   = strcmp({I2.Dimensions.Name},'Nouter');
  Nouter2 = I2.Dimensions(index).Length;

  index  = strcmp({I2.Dimensions.Name},'Ninner');
  Ninner2 = I2.Dimensions(index).Length;

  switch Algorithm2
 
    case {'I4DVAR', 'IS4DVAR', 'SPLIT_I4DVAR'}

      J2.o   = nc_read(oname, 'TLcost_function');
      J2.b   = nc_read(oname, 'back_function');
      J2.NL  = nc_read(oname, 'NLcost_function');
      J2.act = J2.b + J2.o;

      Niter2 = length(J2.act);
    
    case {'RBL4DVAR', 'SPLIT_RBL4DVAR', 'W4DPSAS'}

      half=0.5;   % by definition of penalty functional
							     
      Jobs   = nc_read(oname, 'Jobs', [], NaN);
      Jb     = nc_read(oname, 'Jb', [], NaN);
      Jmod   = nc_read(oname, 'Jmod', [], NaN);
      Jact   = nc_read(oname, 'Jact', [], NaN);
      Jopt   = nc_read(oname, 'Jopt', [], NaN);
      J2.NL  = nc_read(oname, 'NL_fDataPenalty');
      if (Nouter == 1)
        Jobs(1)=[]; Jb(1)=[]; Jmod(1)=[]; Jact(1)=[]; Jopt(1)=[];
      else
        Jobs(1,:)=[]; Jb(1,:)=[]; Jmod(1,:)=[]; Jact(1,:)=[]; Jopt(1,:)=[];
      end
      J2.o   = Jobs.*half;
      J2.b   = Jb.*half;
      J2.mod = Jmod.*half;
      J2.act = Jact.*half;
      J2.opt = Jopt.*half;

      Niter2 = Ninner2*Nouter2;

    case {'R4DVAR', 'SPLIT_R4DVAR', 'W4DVAR'}
  
      half=0.5;   % by definition of penalty functional

      Jobs   = nc_read(oname, 'Jobs', [], NaN);
      Jb     = nc_read(oname, 'Jb', [], NaN);
      Jmod   = nc_read(oname, 'Jmod', [], NaN);
      Jact   = nc_read(oname, 'Jact', [], NaN);
      Jopt   = nc_read(oname, 'Jopt', [], NaN);
      J2.RP  = nc_read(oname, 'RP_fDataPenalty');
      if (Nouter == 1)
        Jobs(1)=[]; Jb(1)=[]; Jmod(1)=[]; Jact(1)=[]; Jopt(1)=[];
      else
        Jobs(1,:)=[]; Jb(1,:)=[]; Jmod(1,:)=[]; Jact(1,:)=[]; Jopt(1,:)=[];
      end
    
      J2.o   = Jobs.*half;
      J2.b   = Jb.*half;
      J2.act = Jact.*half;
    
      Niter2 = Ninner2*Nouter2;
  end

  iter2 = [1:Niter2];
  J.act2 = J2.act;
  J.o2 = J2.o;
  J.b2 = J2.b;
  
end

%--------------------------------------------------------------------------
% Plot penalty/cost function.
%--------------------------------------------------------------------------

Jmin  = J.Nobs/2;
iter  = [1:Niter];

figure;

ha=plot(iter, log10(J.act(:)), 'c-', 'LineWidth', 2);
hold on;
if (got_overlay)
  Ocolor=[0.2353 0.7019 0.4431];    % sea green
  ha2=plot(iter2, log10(J.act2(:)), 'k^', 'MarkerSize', 8);
      set(ha2, 'markerfacecolor', Ocolor)
end

ho=plot(iter, log10(J.o(:)),   'rs', 'MarkerSize', 4);
   set(ho, 'markerfacecolor', get(ho, 'color'))
hb=plot(iter, log10(J.b(:)),   'b-', 'LineWidth', 2);

line([1 Niter],[log10(Jmin) log10(Jmin)],                               ...
     'LineStyle','--','Color',[0 0 0],'LineWidth',2);

if (isfield(J,'NL'))
  hnl=plot(Niter,log10(J.NL(1,end)+J.b(end)),'kd','MarkerSize',8);
      set(hnl, 'markerfacecolor', 'm')
end
if (isfield(J,'RP'))
  hrp=plot(Niter,log10(J.RP(1,end)+J.b(end)),'kd','MarkerSize',8);
      set(hrp, 'markerfacecolor', 'm')
end


A=axis;
axis([1 Niter A(3) A(4)]);

grid on;

ylabel('log_{10}(J)')
xlabel('Iteration number')

if (got_overlay)
  if (strcmp(Algorithm, Algorithm2))
    label1='J';
    label2='J_{compare}';
  else
    label1=strcat('J_{',untexlabel(Algorithm),'}');
    label2=strcat('J_{',untexlabel(Algorithm2),'}');
  end
else
  label1='J';
end

if (isfield(J,'NL'))
  if (got_overlay)
    legend(label1,label2,'J_o','J_b','J_{min}','J_{NL}','Location','Southeast')
  else
    legend(label1,'J_o','J_b','J_{min}','J_{NL}','Location','Southeast')
  end
elseif (isfield(J,'RP'))
  if (got_overlay)
    legend(label1,label2,'J_o','J_b','J_{min}','J_{RP}','Location','Southeast')
  else
    legend(label1,'J_o','J_b','J_{min}','J_{RP}','Location','Southeast')
  end
else
  if (got_overlay)
    legend(label1,label2,'J_o','J_b','J_{min}','Location','Southeast')
  else
    legend(label1,'J_o','J_b','J_{min}','Location','Southeast')
  end
end

if (got_overlay)
  title([untexlabel(Algorithm),'/',untexlabel(Algorithm2),                 ...
	' Cost Functions, ',                                               ...
        ' Ninner=',num2str(Ninner),', Nouter=',num2str(Nouter)]);
else
  title([untexlabel(Algorithm), ' Cost Functions, ',                       ...
        ' Ninner=',num2str(Ninner),', Nouter=',num2str(Nouter)]);
end

return
