% bdry_com
% This script writes out the boundary statements for TPAR files created
% using ww3_swan_input.m, ww3_specpoints.m, ww3gb_2TPAR.m and ww3nc_2TPAR.m

load specpts.mat;
specsize=size(specpts,1);
ofile='Bound_spec_command';
%
if (use_field_output)
  fid=fopen(ofile,'w');
  fprintf(fid,'& **   COPY THESE LINES TO SWAN INPUT FILE   ********\n');
  fprintf(fid,'& TPAR Boundary files  ******************************\n');
  fprintf(fid,'BOUND SHAPESPEC JONSWAP MEAN DSPR DEGREES\n');
  for bd=1:specsize
     fprintf(fid,'BOUNDSPEC SEGMENT XY '); 
     fprintf(fid,'%3.4f %3.4f %3.4f %3.4f ',specpts(bd,1:4));
     fprintf(fid,'VARIABLE FILE 0 ''../forcings/TPAR');
     pg=num2str(length(num2str(bd)));
     eval(['fprintf(fid,''%',pg,'g'',bd);']);
     fprintf(fid,'.txt''\n');
  end
  fprintf(fid,'\n');
  fclose(fid);
end
%
if (use_partition_data)
  fid=fopen(ofile,'w');
  fprintf(fid,'& **   COPY THESE LINES TO SWAN INPUT FILE   ********\n');
  fprintf(fid,'& 2D Spec Boundary files  ****************************\n');
  for bd=1:specsize
     fprintf(fid,'BOUNDSPEC SEGMENT IJ '); 
     fprintf(fid,'%5i %5i %5i %5i ',specpts(bd,5:8));
     fprintf(fid,'VARIABLE FILE 0 ''../forcings/Spc2d');
     pg=num2str(length(num2str(bd)));
     eval(['fprintf(fid,''%',pg,'g'',bd);']);
     fprintf(fid,'.spc2d''\n');
  end
  fprintf(fid,'\n');
  fclose(fid);
end
