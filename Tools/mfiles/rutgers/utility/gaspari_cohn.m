function F = gaspari_cohn(c)

% GASPARI_COHN: Gaspari-Cohn and Gausian Correlations
%
% F = gaspari_cohn(c)
%
% It computes and plots the Gaspari-Cohn (1999) correlation function
% based on the support radius, c. It also plots the Gaussian correlation
% with length-scale, L.  For isotropic covatiances L = c * sqrt(0.3).
%
% On Input:
%
%    c         Gaspari-Cohn support radius (scalar, km)
%
% On Output:
%
%    F         Isotropic correlation functions (struct)
%

% Set distance (km).

d = -400:1:400;

% Compute Gaussian correlation function.

L = c * sqrt(0.3);

G = exp(- 0.5 * (d ./ L).^2);

% Compute Gaspari-Cohn correlation function.

GC = zeros(size(d));

ind1 = find(0 <= abs(d) & abs(d) <= c);
ind2 = find(c <= abs(d) & abs(d) <= 2*c);
ind3 = find(2*c <= abs(d));

doc  = d ./ c;
adoc = abs(d) ./ c;

GC(ind1) = - (1.0/4.0) * adoc(ind1).^5 +  ...
             (1.0/2.0) *  doc(ind1).^4 +  ...
             (5.0/8.0) * adoc(ind1).^3 -  ...
             (5.0/3.0) *  doc(ind1).^2 +  ...
             1.0;

GC(ind2) =   (1.0/12.0) * adoc(ind2).^5 - ...
             (1.0/2.0)  *  doc(ind2).^4 + ...
             (5.0/8.0)  * adoc(ind2).^3 + ...
             (5.0/3.0)  *  doc(ind2).^2 - ...
             5.0        * adoc(ind2)    + ...
             4.0                        - ...
             (2.0/3.0)  *  1./adoc(ind2);

GC(ind3) = 0.0;

% Gaspari-Cohn covolution (gc) and box car (bc) functions.

ind     = find (-c <= d & d <= c);
bc      = zeros(size(d));
bc(ind) = 1;

% Plot correlations.

figure;

ymin = min(GC(:));
ymax = max(GC(:));
ymax = max(ymax, 1.1);

h=plot(d,G,'b-',                      ...
       d,GC,'r-',                     ...
       [L L],[ymin ymax], 'b--',      ...
       [L/2 L/2],[ymin ymax], 'b--',  ...
       [c c],[ymin ymax], 'r--');


axis([-Inf Inf 0 ymax]);
grid on;

Llabel  = ['L = ', num2str(L) ' km'];
L2label = ['L/2 = ', num2str(L/2) ' km'];
clabel  = ['c = ', num2str(c) ' km'];
legend('Gaussian, L', 'Gaspar-Cohn, c',                   ...
       Llabel, L2label, clabel,                           ...
       'Location', 'northwest');
xlabel('Distance (km)');
ylabel('Correlation Scale');

hold on;
text(L, 0.1, 'L', 'Color','b', 'FontWeight', 'bold',      ...
     'HorizontalAlignment','center','FontSize',10);
text(L/2, 0.1, 'L/2', 'Color','b', 'FontWeight', 'bold',  ...
     'HorizontalAlignment','center','FontSize',10);
text(c, 0.1, 'c', 'Color','r', 'FontWeight', 'bold',      ...
     'HorizontalAlignment','center','FontSize',10);
text(-c, 0.1, '-c', 'Color','r', 'FontWeight', 'bold',    ...
     'HorizontalAlignment','center','FontSize',10);
hold off


F.c = c;
F.L = L;
F.d = d;
F.G = G;
F.GC = GC;

return
