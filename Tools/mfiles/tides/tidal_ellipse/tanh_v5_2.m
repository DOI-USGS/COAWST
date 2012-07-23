function y = tanh(xin)
%TANH   Hyperbolic tangent.
%   TANH(X) is the hyperbolic tangent of the elements of X.

%   C. Moler, 6-17-92
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.4 $  $Date: 1997/11/21 23:28:39 $

%   21 Oct 92 (jmh) -- fixed for matrix inputs
%        2 Sep 93 (jmh) -- fixed for null matrix input

% Complex argument
if ~isreal(xin) %  2 Sep 93 (jmh)
    tr = tanh(real(xin));
    ti = i*tan(imag(xin));
    y = (tr + ti)./(1 + tr.*ti);

% Real argument
else

   % Reference: W. J. Cody and W. Waite, "Software Manual
   % for the Elementary Functions", 1980, chapter 13.
   
   xbig = 18.715;                       % 0.5*log(4/eps)
   xmid = 0.5493061443340549;           % log(3)/2
   xsmall = 1.49e-8;                    % sqrt(eps)
   x = abs(xin);
   y = zeros(size(x));

   % tanh(x) == 1 to working precision
   k = find(xbig <= x);
   if ~isempty(k)           % 21 Oct 92 (jmh)
      y(k) = ones(size(find(k)));       % added FIND -- 20 Jul 92 (jmh)
   end

   % 1/2 < tanh(x) < 1
   k = find((xmid < x) & (x < xbig));
   if ~isempty(k)           % 21 Oct 92 (jmh)
      y(k) = 1 - 2./(exp(2*x(k))+1);
   end

   % x < tanh(x) <= 1/2
   k = find((xsmall < x) & (x <= xmid));
   if ~isempty(k)           % 21 Oct 92 (jmh)
      p = [-0.16134119023996228e4 -0.99225929672236083e2 -0.96437492777225470];
      q = [0.48402357071988689e4 0.22337720718962313e4  0.11274474380534949e3];
      xx = x(k).^2;
      y(k) = x(k) + x(k).*(xx.*(((p(3)*xx + p(2)).*xx + p(1))) ./ ...
                          (((xx + q(3)).*xx + q(2)).*xx + q(1)));
   end

   % tanh(x) == x to working precision
   k = find((x <= xsmall) | isnan(x));
   if ~isempty(k)       % 21 Oct 92 (jmh)
      y(k) = x(k);
   end

   y = sign(xin).*y;
end
