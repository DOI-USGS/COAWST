function mapout = cmap_odv(varargin)
% OceanDataView (ODV) colormaps
%
% map = cmap_odv([string],[number])
%
% All inputs optional in any order
%
%  string = name of OVD palette 
%          'list' or 'pal' will dir the available palette names
%          default is 'Odv_437';
%  number = rescale map to have this many color entries
%
% John Wilkin - March 2019

pal = 'Odv_437';
n = [];

List = dir(fullfile(cd, '**', 'odv_palettes', '.'));
if (length(List) > 0)
  odv_pal = List(1).folder;
else
  odv_pal = '/home/arango/ocean/repository/git/roms_matlab/colormaps/odv_palettes';
end

for k=1:nargin
  opt = varargin{k};
  if ischar(opt)
    if any(strcmp(opt,{'list','pal'}))
      unix(['dir ', odv_pal]);
    else
      pal = opt;
    end
  else
    n = opt;
  end
end
palfile = fullfile(odv_pal, [pal '.pal']);
nrgb = load(palfile);

nrgb(:,1)=[];
nrgb(nrgb<0)=0;
nrgb(nrgb>1)=1;

nrgb = nrgb(33:145,:);

map = nrgb;
if ~isempty(n)
  t = linspace(1,length(map),n);
  map = interp1(1:length(map),map,t);
end

if nargout > 0
  mapout = map;
end

