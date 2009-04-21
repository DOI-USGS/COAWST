function SWANInstat(timestart,timeend,gridname,r,specpts,windstart,windend,modelgrid)
ncg=netcdf(modelgrid);
[eta xi]=size(ncg{'lat_rho'}(:));
ncclose

cd D:\Temporary\COAWST\SWAN\runstat\
ofile=['INPUT'];
fid=fopen(ofile,'w');
fprintf(fid,'PROJECT ''Ctest runstat'' ''''\n');
fprintf(fid,'''''\n');
fprintf(fid,'''''\n');
fprintf(fid,'''''\n');
fprintf(fid,'MODE NONSTATIONARY TWODIMENSIONAL\n');
fprintf(fid,'SET DEPMIN 0.10 INRHOG 1 NAUTICAL\n');
fprintf(fid,'COORDINATES SPHERICAL\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'CGRID CURVILINEAR ')
fprintf(fid,'%3g ',xi-1);
fprintf(fid,'%3g ',eta-1);
fprintf(fid,'EXC 9.999000e+003 &\n');
fprintf(fid,'CIRCLE 60 0.04 1.0 20\n');
fprintf(fid,'\n');
fprintf(fid, 'READGRID COORDINATES 1 ''../forcings/grid_coord_')
fprintf(fid,gridname)
fprintf(fid,'.grd'' 4 0 0 FREE\n');
fprintf(fid,'\n');
fprintf(fid,'INPGRID BOTTOM CURVILINEAR 0 0 ')
fprintf(fid,'%3g ',xi-1);
fprintf(fid,'%3g ',eta-1);
fprintf(fid,'EXC 9.999000e+003\n');
fprintf(fid, 'READINP BOTTOM 1 ''../forcings/roms_bathy_')
fprintf(fid,gridname)
fprintf(fid,'.bot'' 4 0 FREE\n');
fprintf(fid,'\n');

fprintf(fid,'INPGRID WIND CURVILINEAR 0 0 ')
fprintf(fid,'%3g ',xi-1); 
fprintf(fid,'%3g ',eta-1); 
fprintf(fid,'EXC 9.999000e+003 &\n');
 
fprintf(fid,'NONSTATIONARY ')
fprintf(fid,'%8.6f',str2num(windstart))% 20030901.013000 
fprintf(fid,' 3 HR ')
fprintf(fid,'%8.6f',str2num(windend))
fprintf(fid,'\n');

 fprintf(fid,'READINP WIND 1 ''../forcings/comb_wind_update.dat')
 fprintf(fid,''' 4 0 FREE\n' );
fprintf(fid,'\n');

fprintf(fid,'& Boundary files  ****************************************\n');
fprintf(fid,'BOUND SHAPESPEC JONSWAP PEAK DSPR DEGREES\n');
fprintf(fid,'\n');
for bd=2:9
   fprintf(fid,'BOUNDSPEC SEGMENT XY '); 
   fprintf(fid,'%3.4f %3.4f %3.4f %3.4f ',specpts(bd,:));
   fprintf(fid,'VARIABLE FILE 0 ''../forcings/TPAR')
   fprintf(fid,'%1g',bd);
   fprintf(fid,'.txt''\n');
end
for bd=10:r
   fprintf(fid,'BOUNDSPEC SEGMENT XY '); 
   fprintf(fid,'%3.4f %3.4f %3.4f %3.4f ',specpts(bd,:));
   fprintf(fid,'VARIABLE FILE 0 ''../forcings/TPAR')
   fprintf(fid,'%2g',bd);
   fprintf(fid,'.txt''\n');
end
fprintf(fid,'\n');
fprintf(fid,'& Restart name **********************************\n');
fprintf(fid,'& INITIAL HOTSTART ''ctest.hot''\n');
fprintf(fid,'\n');
fprintf(fid,'& PHYSICS  ****************************************\n');
fprintf(fid,'&BREAKING CONSTANT 1.0 0.73\n');
fprintf(fid,'FRICTION MADSEN 0.05\n');
fprintf(fid,'& OFF QUAD\n');%ADD & here for wind
fprintf(fid,'GEN3\n');
fprintf(fid,'PROP BSBT\n');
fprintf(fid,'\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''depth.mat'' LAY 4 DEPTH 1.   OUTPUT ')
fprintf(fid,'%8.6f ',timestart);
fprintf(fid,'1 HR\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''dissip.mat'' LAY 4 DISSIP 1.   OUTPUT ')
fprintf(fid,'%8.6f ',timestart);
fprintf(fid,'1 HR\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''force.mat'' LAY 4 FORCE 1.   OUTPUT ')
fprintf(fid,'%8.6f ',timestart);
fprintf(fid,'1 HR\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''hsig.mat'' LAY 4 HSIGN 1.   OUTPUT ')
fprintf(fid,'%8.6f ',timestart);
fprintf(fid,'1 HR\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''rtp.mat'' LAY 4 RTP 1.   OUTPUT ')
fprintf(fid,'%8.6f ',timestart);
fprintf(fid,'1 HR\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''setup.mat'' LAY 4 SETUP 1.   OUTPUT ')
fprintf(fid,'%8.6f ',timestart);
fprintf(fid,'1 HR\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''tmbot.mat'' LAY 4 TMBOT 1.   OUTPUT ')
fprintf(fid,'%8.6f ',timestart);
fprintf(fid,'1 HR\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''ubot.mat'' LAY 4 UBOT 1.   OUTPUT ')
fprintf(fid,'%8.6f ',timestart);
fprintf(fid,'1 HR\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''vel.mat'' LAY 4 VEL 1.   OUTPUT ')
fprintf(fid,'%8.6f ',timestart);
fprintf(fid,'1 HR\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''wdir.mat'' LAY 4 DIR 1.   OUTPUT ')
fprintf(fid,'%8.6f ',timestart);
fprintf(fid,'1 HR\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''wlen.mat'' LAY 4 WLEN 1.   OUTPUT ')
fprintf(fid,'%8.6f ',timestart);
fprintf(fid,'1 HR\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''wind.mat'' LAY 4 WIND 1.   OUTPUT ')
fprintf(fid,'%8.6f ',timestart);
fprintf(fid,'1 HR\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''wind.mat'' LAY 4 WIND 1.   OUTPUT ')
fprintf(fid,'%8.6f ',timestart);
fprintf(fid,'1 HR\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''qb.mat'' LAY 4 QB 1.   OUTPUT ')
fprintf(fid,'%8.6f ',timestart);
fprintf(fid,'1 HR\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''xp.mat'' LAY 4 XP 1.\n');
fprintf(fid,'BLOCK ''COMPGRID'' NOHEADER ''yp.mat'' LAY 4 YP 1.\n');
fprintf(fid,'\n');
fprintf(fid,'COMPUTE STAT ');
fprintf(fid,'%8.6f ',timestart);
fprintf(fid,'\n');
fprintf(fid,'HOTFILE ''../forcings/ctest.hot''\n');
fprintf(fid,'\n');
fprintf(fid,'STOP\n');
fprintf(fid,'\n');
fclose(fid);
cd ../../;