function [F,data] = mpi_dump(fname)

% MPI_DUMP: Reads ROMS ascii tiled file for parallel debugging
%
% [F,data] = mpi_dump(fname)
%
% On Input:
%
%    fname      ASCII filename(s), (string, cell string)
%                 For example, fname = 'fort.100'
%                              fname = {'fort.100', 'fort.101', ...}
%
% On Output:
%
%    F          File data (sctruc or struct vector)
%
%                 F.name            F(n).name
%                 F.Istr            F(n).Istr
%                 F.Iend            F(n).Iend
%                 F.Jstr            F(n).Jstr
%                 F.Iend            F(n).Jend
%                 F.value           F(n).value
%  
%    data       Concatenated data, if multi-files.
%

data = [];

if iscell(fname)
  nfiles = length(fname);
  Iadd = 0;
  Jadd = 0;
  for n = 1:nfiles
    fid = fopen(fname{n}, 'r');
    [fdat count] = fscanf(fid, '%g', inf);
    fclose(fid); 

    Istr = fdat(1);
    Iend = fdat(2);
    Jstr = fdat(3);
    Jend = fdat(4);
    fdat = fdat(5:count);
    
    F(n).name = fname{n};
    F(n).Istr = Istr;
    F(n).Iend = Iend;
    F(n).Jstr = Jstr;
    F(n).Jend = Jend;

    Imax = Iend-Istr+1;
    Jmax = Jend-Jstr+1;
    F(n).value = reshape(fdat,Imax,Jmax);

    if (Istr == 0)
      Iadd=1;
    end
    if (Jstr == 0)
      Jadd=1;
    end
  end

  Imax = max(F(:).Iend)+Iadd;
  Jmax = max(F(:).Jend)+Jadd;
  
  Aglobal = nan(Imax,Jmax);
  for n = 1:nfiles
    Is = F(n).Istr+Iadd;
    Ie = F(n).Iend+Iadd;
    Js = F(n).Jstr+Jadd;
    Je = F(n).Jend+Jadd;
    data(Is:Ie,Js:Je) = F(n).value;
  end

  x=(1:Imax)-Iadd;
  y=(1:Jmax)-Jadd;  
  X = repmat(x', [1 Jmax]);
  Y = repmat(y,  [Imax 1]);

  figure;
  pcolorjw(X, Y, data);
  colorbar;
  shading flat;
  colormap(flipud(mpl_Paired(256)));

else

  fid = fopen(fname, 'r');
    [fdat count] = fscanf(fid, '%g', inf);
    fclose(fid);

    F.name = fname;
    F.Istr = fdat(1);
    F.Iend = fdat(2);
    F.Jstr = fdat(3);
    F.Jend = fdat(4);

    fdat = fdat(5:count);
    Imax = F.Iend-F.Istr+1;
    Jmax = F.Jend-F.Jstr+1;

    F.value = reshape(fdat,Imax,Jmax);
end

return
