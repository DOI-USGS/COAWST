function S=fstats(field)

%
% FSTAST: Computes the field statistics
%
% S=fstast(field)
%
% It computes the statistics of requested field.
%
% On Input:
%
%    field       Field to process (vector or array)
%
% On Output:
%
%    S           Field statistics (struct)
%

% svn $Id: fstats.m 996 2020-01-10 04:28:56Z arango $
%===========================================================================%
%  Copyright (c) 2002-2020 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

S.min    = min(field(:));
S.max    = max(field(:));
S.mean   = mean(field(:));
S.median = median(field(:));
S.std    = std(field(:));
S.var    = var(field(:));

return
