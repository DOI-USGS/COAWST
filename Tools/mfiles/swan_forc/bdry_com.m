% bdry_com
% This script writes out the boundary statements for TPAR files created
% using ww3_swan_input.m, ww3_specpoints.m, ww3gb_2TPAR.m and ww3nc_2TPAR.m

load specpts.mat;
lth=length(specpts);
ofile='Bound_spec_command';
fid=fopen(ofile,'w');

fprintf(fid,'& **   COPY THESE LINES TO SWAN INPUT FILE   ********\n');
fprintf(fid,'& Boundary files  ***********************************\n');
fprintf(fid,'BOUND SHAPESPEC JONSWAP MEAN DSPR DEGREES\n');
for bd=1:lth
   fprintf(fid,'BOUNDSPEC SEGMENT XY '); 
   fprintf(fid,'%3.4f %3.4f %3.4f %3.4f ',specpts(bd,:));
   fprintf(fid,'VARIABLE FILE 0 ''../forcings/TPAR');
   pg=num2str(length(num2str(bd)));
   eval(['fprintf(fid,''%',pg,'g'',bd);']);
   fprintf(fid,'.txt''\n');
end
fprintf(fid,'\n');
fclose(fid);
