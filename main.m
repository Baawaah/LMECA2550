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

%% Init the data
clear all
%Flight parameter
global pdata alphas_naca16_509_m06 cls_naca16_509_m06 cds_naca16_509_m06
load 'naca16-509-m06_clcd.mat'

pdata.M          = 0.5;
pdata.z          = 20000;
pdata.B          = 4;
pdata.chord      = 0.25;
pdata.omega      = 0.2;
%pdata.Beta0      = 0.26;
pdata.R          = 1.7;
 [p0,ro,Temp,gamma] = stdatm(pdata.z);
pdata.V          = pdata.M*sqrt(gamma*(p0/ro));
pdata.ro         = ro;
pdata.a(1)       = 0.4;
pdata.a_prime(1) = 0.4;

%% BEM Theory resolving via Gauss-Legendre
A = [0.1,1.7]; % Interval
N = 3;
h = 1/N;
X = [0.7746 0 -0.7746];
W = [0.5556 0.8889 0.5556];
L = A(2)-A(1);
B = linspace(A(1),A(2),N);

M = 0.5;
V = M*sqrt(gamma*(p0/ro));
set = [15 25 35 45 55 62]*(2*pi/360);

nNmax = 130;
T  = zeros(length(set),nNmax);
Q  = zeros(length(set),nNmax);
kt = zeros(length(set),nNmax);
kq = zeros(length(set),nNmax);
kp = zeros(length(set),nNmax); 
np = zeros(length(set),nNmax);
J  = zeros(length(set),nNmax);

nNrange = [100 200 200 300 350 400];
nrange = [ 45 350; 28 350; 19.3 350; 13.5 350; 9.4 350; 7.4 350]; %
for k = 1 : length(set);
  pdata.a = 0.2;
  pdata.a_prime = 0.2;
  pdata.Beta0 = set(k);
  nN = nNrange(k);
  n = linspace(nrange(k,1),nrange(k,2),nN);
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
 kp(k,j) = (Q(k,j)*big_omega)/(pdata.ro*n(j)^3*((2*pdata.R)^5));
 np(k,j) = J(k,j)*kt(k,j)/kp(k,j);

 if J(k,j) < 0.15 
    kt(k,j) = NaN; 
    kq(k,j) = NaN;
    kp(k,j) = NaN;
    np(k,j) = NaN;
 end
end
end

%% Plotting the result
figure; 
for k = 1 : length(set);    
    %subplot(2,2,1);
    plot(J(k,:),kt(k,:));
    %title('kt');
    axis([0 5 0 inf])
    legend('15°','25°','35°','45°','55°', '62°');
    hold on;
end
figure;
for k = 1 : length(set);
    %subplot(2,2,2);
    plot(J(k,:),kq(k,:));
    %title('kq');
    axis([0 5 0 inf])
    hold on;
end    
figure;
for k = 1 : length(set);    
    %subplot(2,2,3);
    plot(J(k,:),kp(k,:));
    %title('kp');
    axis([0 5 0 inf])
    hold on;
end    
figure;
for k = 1 : length(set);    
    %subplot(2,2,4);
    plot(J(k,:),np(k,:));
    title('np');
    axis([0 5 0 inf])
    hold on;
end
%% Maximum speed
P  = 1.118550*10^6; %Watt 
nmotor = 49.9985;
nvmax  = nmotor * 0.477;
kpvmax = P/(pdata.ro*nvmax*nvmax*nvmax*((pdata.R*2)^5)); 
Jvmax = [0.15 0.55 1.03 1.64 2.43 3.71] %Founded via graph
Vvmax = Jvmax *(nvmax*pdata.R*2);
Mvmax = Vvmax/sqrt(gamma*(p0/ro));

A = [0.1,1.7]; % Rayon de la palme
N = 10;
X = [0.7746 0 -0.7746];
W = [0.5556 0.8889 0.5556];
B   = linspace(A(1),A(2),N);
nNmax = 100;
M45 = linspace(0.01,1,nNmax);
%T2  = zeros(1,nNmax);
%Q2  = zeros(1,nNmax);
pdata.a = 0.1;
pdata.a_prime = 0.1;
pdata.Beta0 = 35*(2*pi/360);
pdata.omega = 0.1;
Tvmax = zeros(1,nNmax);
Qvmax = zeros(1,nNmax);
for j = 1 : length(M45);
  V = M45(j)*sqrt(gamma*(p0/ro));  
  Tvmax(j) = 0;
  Qvmax(j) = 0;
 for i = 1 : length(B)-1
  % Calcule du Thrust par rapport au Mach   
  alpha = (B(i+1)-B(i))/2;
  beta  = (B(i+1)+B(i))/2; 
  big_omega = nvmax*2.0*pi;
  [T1,Q1] = fTQprime( alpha*X(1) + beta,V,big_omega);
  [T2,Q2] = fTQprime( alpha*X(2) + beta,V,big_omega);
  [T3,Q3] = fTQprime( alpha*X(3) + beta,V,big_omega);
  Tvmax(j) =  alpha.*( W(1)*T1 + W(2)*T2 + W(3)*T3) + Tvmax(j);
  Qvmax(j) =  alpha.*( W(1)*Q1 + W(2)*Q2 + W(3)*Q3) + Qvmax(j);
 end
end
Betavmax = zeros(1,length(B));
for i = 1 : length(B)
% Calcule de la distribution de beta0 en fonction du rayon.
  Betavmax(i) = atan( 2*pi*0.75*pdata.R*tan(pdata.Beta0)/(2*pi*B(i)));
end

figure;
plot(M45,Tvmax);
axis([0 1 0 inf])
figure;
plot(B,Betavmax);

%%
