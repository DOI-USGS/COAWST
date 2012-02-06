function m_proj(proj,varargin)
% M_PROJ  Initializes map projections info, putting the result into a structure
%
%         M_PROJ('set') tells you the current state
%         M_PROJ('get') gives you a list of all possibilities
%         M_PROJ('get','proj name') gives info about a projection in the 
%                                   'get' list.
%         M_PROJ('proj name','property',value,...) initializes a projection.
%
%         see also M_GRID, M_LL2XY, M_XY2LL.

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.


global MAP_PROJECTION MAP_VAR_LIST

% Get all the projections
projections=m_getproj;

if nargin==0, proj='usage'; end;

proj=lower(proj);

switch proj,

  case 'get',              % Print out their names
    if nargin==1,
      disp(' ');
      disp('Available projections are:'); 
      for k=1:length(projections),
        disp(['     ' projections(k).name]);
      end;
    else
      k=m_match(varargin{1},projections(:).name);
      eval(['X=' projections(k).routine '(''set'',projections(k).name);']);
      disp(X);
    end;

  case 'set',              % Get the values of all set parameters
    if nargin==1,
      if isempty(MAP_PROJECTION),
         disp('No map projection initialized');
         m_proj('usage');
      else
         k=m_match(MAP_PROJECTION.name,projections(:).name);
         eval(['X=' projections(k).routine '(''get'');']);
         disp('Current mapping parameters -');
         disp(X);
      end;
    else
      k=m_match(varargin{1},projections(:).name);
      eval(['X=' projections(k).routine '(''get'');']);
      disp(X);
    end;

  case 'usage',
    disp(' ');
    disp('Possible calling options are:');
    disp('  ''usage''                    - this list');
    disp('  ''set''                      - list of projections');
    disp('  ''set'',''projection''         - list of properties for projection');
    disp('  ''get''                      - get current mapping parameters (if defined)');
    disp('  ''projection'' <,properties> - initialize projection\n');

 otherwise                % If a valid name, give the usage.
    k=m_match(proj,projections(:).name);
    MAP_PROJECTION=projections(k);
    eval([ projections(k).routine '(''initialize'',projections(k).name,varargin{:});']);

end;

%---------------------------------------------------------
function projections=m_getproj;
% M_GETPROJ Gets a list of the different projection routines
%           and returns a structure containing both their
%           names and the formal name of the projection.
%           (used by M_PROJ).

% Rich Pawlowicz (rich@ocgy.ubc.ca) 9/May/1997
%
% 9/May/97 - fixed paths for Macs (thanks to Dave Enfield)
%
% 7/05/98 - VMS pathnames (thanks to Helmut Eggers)

% Get all the projections

lpath=which('m_proj');
fslashes=findstr(lpath,'/');
bslashes=findstr(lpath,'\');
colons=findstr(lpath,':');
closparantheses=findstr(lpath,']');
if ~isempty(fslashes),
  lpath=[ lpath(1:max(fslashes)) 'private/'];
elseif ~isempty(bslashes),
  lpath=[ lpath(1:max(bslashes)) 'private\'];
elseif ~isempty(closparantheses),       % for VMS computers only, others don't use ']' in filenames
  lpath=[ lpath(1:max(closparantheses)-1) '.private]'];
else,
  lpath=[ lpath(1:max(colons)) 'private:'];
end;

w=dir([lpath 'mp_*.m']);

l=1;
projections=[];
for k=1:length(w),
 funname=w(k).name(1:(findstr(w(k).name,'.'))-1);
 projections(l).routine=funname;
 eval(['names= ' projections(l).routine '(''name'');']);
 for m=1:length(names);
   projections(l).routine=funname;
   projections(l).name=names{m};
   l=l+1;
 end;
end;


%----------------------------------------------------------
function match=m_match(arg,varargin);
% M_MATCH Tries to match input string with possible options

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997

match=strmatch(lower(arg),cellstr(lower(char(varargin))));

if length(match)>1,
  error(['Projection ''' arg ''' not a unique specification']);
elseif isempty(match),
  error(['Projection ''' arg ''' not recognized']);
end;
