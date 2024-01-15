function out=replace_nan(var, value)
%
% jcw 24Jan2024
%

% for some reason matlab did not know isnan
%
out=var;
out(isnan(out))=value;

