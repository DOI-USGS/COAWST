function k = waven(T,h)
% WAVEN  Wavenumber k for gravity waves, iterating with Newtons method
%
% k = waven(T,h)
%
% Input: period T [s] and water depth h [m]
%
% Solve:  0 = kh*tanh(kh)-w**2*h/g using Newton's method

% Chris Sherwood, USGS
% March 17, 1999
% (after Press et al., 1986)
%
%...MAX 20 ITERATIONS, SHOULD BE CORRECT TO .0001

jmax=20;
xacc=0.0001;
g = 9.80665;
w=2*pi/T;
w2h = w*w*h/g;
%...BOUNDS ARE 0 AND SLIGHTLY MORE THAN W2
x1=0.;
x2=w2h+1.;
%...BEST INITIAL GUESS IS W2
kh=w2h;
for j=1:jmax,
%       ...COMPUTE F(X)
        tx=tanh(kh);
        f=kh*tx-w2h;
%       ...COMPUTE F'(X)
        df=kh-kh*tx*tx + tx;
        dx=f/df;
        kh=kh-dx;
        if((x1-kh)*(kh-x2)<0),error('WAVEN out of brackets'),end;
        if(abs(dx)<xacc),break,end;
end
if(j==jmax),error('waven exceeded maximum iterations'),end
k=kh/h;






