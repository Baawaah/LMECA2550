% P Chatelain
% Look-up for the behavior of the NACA16-509


%%% Initialization in the main matlab script:
%%% Code snipplet to load the airfoil data in main file:
% global alphas_naca16_509_m06 cls_naca16_509_m06 cds_naca16_509_m06
% load 'naca16-509-m06_clcd.mat'





function [cl,cd] = naca16_509_m06(a)
global alphas_naca16_509_m06 cls_naca16_509_m06 cds_naca16_509_m06

a = 180/pi * a;

cl = interp1(alphas_naca16_509_m06,cls_naca16_509_m06,a,'spline',NaN);

oi = find(isnan(cl));

% Polynomial fits

%p1 = polyfit([alphas_naca16_509_m06(1) -80 -90 -100], [-0.246163682864450 -0.02 0 -0.02] , 3);
p1 = [ 0.000001993645465   0.000338284275449   0.012246220243973  -0.184575265468565];
%p2 = polyfit([alphas_naca16_509_m06(end) 12 70 90 110], [0.966112531969310 1.1 0.1 0 0.1] , 4);
p2 = [-0.000000079103690   0.000023256824678  -0.002153261849550   0.055200854989657   0.709111947532640];


cl(oi) = (a(oi)<0) * (-0.246163682864450) + ... 
(a(oi)>0) * (0.966112531969310);


cl(oi) = (a(oi)<0) .* (polyval(p1,a(oi))) + ... 
(a(oi)>0) .* (polyval(p2,a(oi)));


cd = interp1(alphas_naca16_509_m06,cds_naca16_509_m06,a,'spline',NaN);

coef1 =  -0.067317398298573;
coef2 =  0.023120640006175;



oi = find(isnan(cd));
cd(oi) = (a(oi)<0) .* (cds_naca16_509_m06(1)-0.5 * coef1 * (a(oi)-alphas_naca16_509_m06(1)).^2) + ... 
(a(oi)>0) .* (cds_naca16_509_m06(end) + 0.5 * coef2 * (a(oi)-alphas_naca16_509_m06(end)).^2);

cdmax = 2;
cd = min(cd,cdmax*ones(size(cd)));
