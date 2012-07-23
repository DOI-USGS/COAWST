% Demonstrate how to use ap2ep and ep2ap
%
Au=rand(4,3,2);         % so 4x3x2 multi-dimensional matrices are used for the 
Av=rand(4,3,2);         % demonstration.
Phi_v=rand(4,3,2)*360;  % phase lags inputs are expected to be in degrees.
Phi_u=rand(4,3,2)*360;

figure(1)
clf
[SEMA ECC INC PHA w]=ap2ep(Au, Phi_u, Av, Phi_v, [2 3 1]);
figure(2)
clf
[rAu rPhi_u rAv rPhi_v rw]=ep2ap(SEMA, ECC, INC, PHA, [2 3 1]);

%check if ep2ap has recovered Au, Phi_u, Av, Phi_v
max(max(max(abs(rAu-Au))))               %  = 9.9920e-16
max(max(max(abs(rAv-Av))))               %  = 6.6613e-16
max(max(max(abs(rPhi_u-Phi_u))))         %  = 4.4764e-13
max(max(max(abs(rPhi_v-Phi_v))))         %  = 1.1369e-13
max(max(max(max(abs(w-rw)))))            %  = 1.3710e-15
% for the random realization I had, the differences are listed on the right 
% hand of the above column. What are yours?

% The above example function calls have already plotted an ellipse for you. 
% To plot an ellipse separately, you may do
%
figure(3)
clf
plot(real(squeeze(w(2,3,1,:))), imag(squeeze(w(2,3,1,:))));

%here squeeze is needed because w is a multiple dimensional array.


%Zhigang Xu
%Nov. 12, 2000
