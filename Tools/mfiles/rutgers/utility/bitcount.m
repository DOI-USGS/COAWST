function B = bitcount(A)
  
% BITCOUNT: Count the number of set bits in each column of the input
%
% B = bitcount(A)
%
% Count the number of set bits in each column of the input array,
% typecast as a bit vector.
%
% On Input:
%
%    A          MxNx... input array.
%
% On Output:
%
%    B          1xNx... output array of bit counts.
%
% For example, to get the global sum value of an array use:
%
%    checksum = bitcount( A(:) )
  
persistent lutable

% Generate the lookup table.

if isempty(lutable)
  lutable = uint8(sum(dec2bin(0:255) - '0', 2));
end

% Convert to an index into the lookup table.

sz = size(A);
sz(1) = 1;
A = reshape(typecast(A(:), 'uint8'), [], prod(sz));

% Look up the number of set bits for each byte.

try
  A = intlut(A, lutable);
catch
  A = uint16(A) + uint16(1);
  A = lutable(A);
end

% Sum the number of set bits per column.

if (size(A, 1) < 32)
  A = sum(A, 1, 'native');
else
  A = sum(A, 1);
end

B = reshape(A, sz);
