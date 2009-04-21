
ofile=['COAWST_generator.php'];
fid=fopen(ofile,'w');
fprintf(fid,'function mymap(mymap)\n');
fprintf(fid,'{var handler=( { directory: [')
for i=1:3*24
fprintf(fid,'"USEAST_COAWST_')
fprintf(fid,num2str(i))
fprintf(fid,'.gif",'); 
end
fprintf(fid,'] } ) ;\n');
fprintf(fid,'return handler(mymap);}\n');