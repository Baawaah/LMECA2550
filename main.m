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


clear all
%Flight parameter
global pdata alphas_naca16_509_m06 cls_naca16_509_m06 cds_naca16_509_m06
load 'naca16-509-m06_clcd.mat'

%% Init the data
pdata.M          = 0.5;
pdata.z          = 20000;
pdata.B          = 4;
pdata.chord      = 0.25;
pdata.omega      = 0.3;
%pdata.Beta0      = 0.26;
pdata.R          = 1.7;
 [p0,ro,T,gamma] = stdatm(pdata.z);
pdata.V          = pdata.M*sqrt(gamma*(p0/ro));
pdata.ro         = ro;
pdata.a(1)       = 0.4;
pdata.a_prime(1) = 0.4;

%% Gauss-Legendre
A = [0.1,1.7]; % Interval
N = 20;
h = 1/N;
X = [0.7746 0 -0.7746];
W = [0.5556 0.8889 0.5556];
L = A(2)-A(1);
B = linspace(A(1),A(2),N);

nN= 30;
M = 0.5;
n = linspace(10,150,nN);
V = M*sqrt(gamma*(p0/ro));
set = [15 25 35 45 55 62]*(2*pi/360);
T  = zeros(length(set),nN);
Q  = zeros(length(set),nN);
kt = zeros(length(set),nN);
kq = zeros(length(set),nN);
kp = zeros(length(set),nN); 
np = zeros(length(set),nN);
J  = zeros(length(set),nN);
for k = 1 : length(set);
  pdata.Beta0 = set(k);
for j = 1 : length(n);
   
  T(k,j) = 0;
  Q(k,j) = 0;
 for i = 1 : length(B)-1
  alpha = (B(i+1)-B(i))/2;
  beta  = (B(i+1)+B(i))/2; 
  big_omega = n(j)*2.0*pi;
  [T1,Q1] = fTQprime( alpha*X(1) + beta,V,big_omega);
  [T2,Q2] = fTQprime( alpha*X(2) + beta,V,big_omega);
  [T3,Q3] = fTQprime( alpha*X(3) + beta,V,big_omega);
  T(k,j) =  alpha.*( W(1)*T1 + W(2)*T2 + W(3)*T3) + T(k,j);
  Q(k,j) =  alpha.*( W(1)*Q1 + W(2)*Q2 + W(3)*Q3) + Q(k,j);
 end
 J(k,j) = V/( n(j)*(pdata.R*2) );
 kt(k,j) = T(k,j)/(pdata.ro*n(j)^2*((2*pdata.R)^4)); 
 kq(k,j) = Q(k,j)/(pdata.ro*n(j)^2*((2*pdata.R)^5));
 kp(k,j) = (Q(k,j)*big_omega)/(ro*n(j)^3*((2*pdata.R)^5));
 np(k,j) = J(k,j)*kt(k,j)/kq(k,j);
end
end
% Plotting

for k = 1 : length(set);
    subplot(2,2,1);
    plot(J(k,:)',kt(k,:));
    title('kt');
    axis([0 inf 0 inf])
    hold on;
    subplot(2,2,2);
    plot(J(k,:),kq(k,:));
    title('kq');
    axis([0 inf 0 inf])
    hold on;
    subplot(2,2,3);
    plot(J(k,:),kp(k,:));
    title('kp');
    axis([0 inf 0 inf])
    hold on;
    subplot(2,2,4);
    plot(J(k,:),np(k,:));
    title('np');
    axis([0 inf 0 inf])
    hold on;
end

%% 