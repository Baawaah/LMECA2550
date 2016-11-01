%
% LMECA2550 Aircraft Propulsion
% by Thanh-Son TRAN - 8116-12-00
%
%% ==========================================
%Equation Box
% Thrust Coeff
% kt = T/(ro*n*n*(D^4));
% Torque Coeff
% kq = Q/(ro*n*n*(D^5));
% Power Coeff
% kp = P/(ro*n*n*n*(D^5)); for P = Q*omega
% Advance Ratio
% J = V/(n*D);
% Reynold's Number
% ro*n*(D^2)/mu;
% Tip Mach Number;
% M_tip = V/a ou (n*(D/2))/a
%% ==========================================
% Equation à résoudre:
% S = (B*chord)/(2*pi*r)
% dT = pi*S*(r^3)*ro*(omega^2)*((1-a_prime)^2)*lambda1*(1/((cos(psi))^2))*dr
% dQ = pi*S*(r^4)*ro*(omega^2)*((1-a_prime)^2)*lambda2*(1/((cos(psi))^2))*dr
% Pour cela on a besoin de résoudre le système suivant:
%
% a_new       = (1+a      )*(0.5)*S*lambda1/(1-cos(2*psi))
% a_prime_new = (1-a_prime)*(0.5)*S*lambda2/(  sin(2*psi))
% Relaxing for omega = 0.3
% ak = (1-omega)*a(k-1)+omega*a_new
% 
% lambda1 = (Cl*cos(psi)) - (Cd*sin(psi));
% lambda2 = (Cl*sin(psi)) - (Cd*cos(psi));
%
% pour psi connu par :
%  tan(psi) = (V/(omega*r)) * (1+a)/(1-a_prime)

% Calculation of "a" and "a_prime" with error gestion
% r is said to be known

%Flight parameter
global pdata alphas_naca16_509_m06 cls_naca16_509_m06 cds_naca16_509_m06
load 'naca16-509-m06_clcd.mat'

%% Init the data
pdata.M          = 0.5;
pdata.z          = 20000;
pdata.B          = 4;
pdata.chord      = 0.25;
pdata.RPM        = 3000;
pdata.omega      = 0.3;
pdata.n          = pdata.RPM/60.0;
pdata.big_omega  = pdata.n*2.0*pi;
pdata.Beta0      = 0.26;
pdata.R          = 1.7;
 [p0,ro,T,gamma] = stdatm(pdata.z);
pdata.V          = pdata.M*sqrt(gamma*(p0/ro));
pdata.ro         = ro;
pdata.a(1)       = 0.4;
pdata.a_prime(1) = 0.4;
% %% =============
%Sucre syntaxique
r               = 1.7; %parameter
z               = pdata.z;
omega           = pdata.omega; 
%B               = pdata.B; % Number of Blade
chord           = pdata.chord; % Chord length
big_omega       = pdata.big_omega;
Beta0           = pdata.Beta0;
R               = pdata.R;
V               = pdata.V;% [m/s]

%% Gauss-Legendre
A = [0.1,1.7]; % Interval
N = 10;
h = 1/N;
X = [0.7746 0 -0.7746];
W = [0.5556 0.8889 0.5556];
L = A(2)-A(1);
B = linspace(A(1),A(2),N);

M= 0.5;
n = linspace(50,100,20);
V = M*sqrt(gamma*(p0/ro));
for j = 1 : length(n);
   
 acc = 0;
 for i = 1 : length(B)-1
  alpha = (B(i+1)-B(i))/2;
  beta  = (B(i+1)+B(i))/2; 
  big_omega = n(j)*2.0*pi;
  T(j) =  alpha*( W(1)*fTprime( alpha*X(1) + beta,V,big_omega) + W(2)*fTprime( alpha*X(2) + beta,V,big_omega) + W(3)*fTprime( alpha*X(3) + beta,V,big_omega)) + acc;
 end
 J(j) = V/( n(j)*(pdata.R*2) );
 kt(j) = T(j)/(pdata.ro*n(j)*n(j)*((2*pdata.R)^4)); 
end

plot(J,kt);