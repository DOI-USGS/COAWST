function create_inwave_grid(x,y,dx,dy,depth,roms_angle,mask_rho,f,grd_file)

  rho.x=x;
  rho.y=y;  
  rho.dx=dx;
  rho.dy=dy;
  rho.depth=depth;
  rho.angle = roms_angle;
  rho.mask = mask_rho
  rho.f = f;
  projection='mercator';
  spherical='T';
  save temp_jcw33.mat
  eval(['mat2roms_mw(''temp_jcw33.mat'',''',grd_file,''');'])
  !del temp_jcw33.mat
  disp(['Created roms grid -->   ',grd_file])
  
end