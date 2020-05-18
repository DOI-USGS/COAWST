function newstr = untexlabel(str)

% newstr = untexlabel(str)
% 
% Produces a TeX format from a character string. It Builds a string
% which used with a TeX interpreter and shows exactly the string as
% it is.
%
%    If str is a one row string, UNTEXLABEL returns a one row string.
%
%    If str is a multiple row string, UNTEXLABEL is applied to each row,
%    and the vertical concatenation is returned (see strvcat).
%
%    If str is a cell array of strings, UNTEXLABEL returns a cell array of
%    the same size, each containing the correspondant result.
%
%    Example:
%
%    >> untexlabel('c:\matlab\temp\')
%    ans =
%    c:\\matlab\\temp\\
%
%    See also TEXLABEL (matlab function)
%
% Adapted from Giuseppe Ridino 08-Jul-2004 script.
%

% svn $Id: untexlabel.m 996 2020-01-10 04:28:56Z arango $

% Initialize output

newstr = '';
if ~isempty(str),
  if isa(str,'char'),
    newstr = untexlabel_local(str);
  elseif iscellstr(str),
    elements = prod(size(str));
    newstr = cell(size(str));
    for index=1:elements,
      newstr{index} = untexlabel_local(str{index});
    end
  else
    error('Argument must be a char array or a cell array of strings.')
  end
end

%==========================================================================
function newstr = untexlabel_local(str)
%==========================================================================

% This is for a multiple string line

newstr = '';
for index = 1:size(str,1),
  newstr = strvcat(newstr,untexstring(str(index,:)));
end

%==========================================================================
function newstr = untexstring(str)
%==========================================================================

% This is for a single string line

newstr = '';

if ~isempty(str),
  index1 = find(str=='^');        % get '^' index
  index2 = find(str=='_');        % get '_' index
  index3 = find(str=='\');        % get '\' index

  index_end   = [sort([index1,index2,index3]-1) length(str)];  % merge
  index_begin = [1,index_end(1:end-1)+1];

% Build new string

   for counter = 1:length(index_end)-1,
     tok = str(index_begin(counter):index_end(counter));
     newstr = strcat(newstr,tok,'\');
   end

% Add end of str

   counter = length(index_end);
   tok = str(index_begin(counter):index_end(counter));
   newstr = strcat(newstr,tok);
end
