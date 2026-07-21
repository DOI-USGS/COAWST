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

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                   %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.md                            Hernan G. Arango        %
%===========================================================================%

ind = ~isnan(field);

S.checksum = bitcount(field(ind));
S.min      = min(field(ind));
S.max      = max(field(ind));
S.mean     = mean(field(ind));
S.median   = median(field(ind));
S.std      = std(field(ind));
S.var      = var(field(ind));

return
