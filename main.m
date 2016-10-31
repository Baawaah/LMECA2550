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
global alphas_naca16_509_m06 cls_naca16_509_m06 cds_naca16_509_m06
load 'naca16-509-m06_clcd.mat'

V = 50; % [m/s]

r = 1.5; %parameter

omega = 0.3 ; 
B     = 4   ; % Number of Blade
chord = 0.25; % Chord length

a_old       = 0.5;
a_prime_old = 0.5;
psi_old     =   2;
error_a       = 1; % error between a's
error_a_prime = 1;
error_psi     = 1;

a = a_old;
a_prime = a_prime_old;
psi = psi_old;
while abs(error_psi) >= 0.001  
    % Iteration of a and a'
    while (abs(error_a) >= 0.001) || (abs(error_a_prime) >= 0.001)
        % Getting the Cl and Cd
        [Cl,Cd] = naca16_509_m06(a_new);
        
        % Computation of S
        S = (B*chord)/(2*pi*r);
        
        % Computation of lambda1 and lambda2 
        lambda1 = (Cl*cos(psi)) - (Cd*sin(psi));
        lambda2 = (Cl*sin(psi)) + (Cd*cos(psi));
        
        % Computation of a and a'
        a_new_young       = (1+a_new      )*(0.5)*S*lambda1/(1-cos(2*psi));
        a_prime_new_young = (1-a_prime_new)*(0.5)*S*lambda2/(  sin(2*psi));
        
        % Relaxing the value
        a_new       = (1-omega)* a_old       + omega*a_new_young;
        a_prime_new = (1-omega)* a_prime_old + omega*a_prime_new_young;
        
        % Error computation
        error_a       = a_new - a_old;
        error_a_prime = a_prime_new - a_prime_old;
        a_old = a_new;
        a_prime_old = a_prime_new;
    end
    % Computation of psi
    psi = atan((V/(omega*r)) * (1+a)/(1-a_prime));
    
    % Error computation
    error_psi = psi - psi_old;
    psi_old = psi;
end
disp('       a        a_p       err_a    err_a_p   psi   err_psi');
disp([a_new,a_prime_new,error_a,error_a_prime,psi,error_psi]);

