function [x_full_grid,y_full_grid]=create_extra_rho_grid(x_rho,y_rho)
  %create an extra set of bounding rho points.
  zsize=size(x_rho);
  x_full_grid=ones(zsize+2);
  x_full_grid(2:end-1,2:end-1)=x_rho;
  x_full_grid(2:end-1,1)=x_rho(:,1)-(x_rho(:,2)-x_rho(:,1));
  x_full_grid(2:end-1,end)=x_rho(:,end)+(x_rho(:,end)-x_rho(:,end-1));
  x_full_grid(1,2:end-1)=x_rho(1,:)-(x_rho(2,:)-x_rho(1,:));
  x_full_grid(end,2:end-1)=x_rho(end,:)+(x_rho(end,:)-x_rho(end-1,:));
  x_full_grid(1,1)=x_full_grid(1,2)-(x_full_grid(1,3)-x_full_grid(1,2));
  x_full_grid(1,end)=x_full_grid(1,end-1)+(x_full_grid(1,end-1)-x_full_grid(1,end-2));;
  x_full_grid(end,1)=x_full_grid(end,2)-(x_full_grid(end,3)-x_full_grid(end,2));
  x_full_grid(end,end)=x_full_grid(end,end-1)+(x_full_grid(end,end-1)-x_full_grid(end,end-2));;
%
  y_full_grid=ones(zsize+2);
  y_full_grid(2:end-1,2:end-1)=y_rho;
  y_full_grid(2:end-1,1)=y_rho(:,1)-(y_rho(:,2)-y_rho(:,1));
  y_full_grid(2:end-1,end)=y_rho(:,end)+(y_rho(:,end)-y_rho(:,end-1));
  y_full_grid(1,2:end-1)=y_rho(1,:)-(y_rho(2,:)-y_rho(1,:));
  y_full_grid(end,2:end-1)=y_rho(end,:)+(y_rho(end,:)-y_rho(end-1,:));
  y_full_grid(1,1)=y_full_grid(1,2)-(y_full_grid(1,3)-y_full_grid(1,2));
  y_full_grid(1,end)=y_full_grid(1,end-1)+(y_full_grid(1,end-1)-y_full_grid(1,end-2));;
  y_full_grid(end,1)=y_full_grid(end,2)-(y_full_grid(end,3)-y_full_grid(end,2));
  y_full_grid(end,end)=y_full_grid(end,end-1)+(y_full_grid(end,end-1)-y_full_grid(end,end-2));;


