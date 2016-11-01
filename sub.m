
function [a,a_prime,psi,lambda1,lambda2] = sub(r,V,a_old,a_prime_old)
global pdata;


p_ref0 = 2*pi*0.75*pdata.R*tan(pdata.Beta0);
Beta = atan(p_ref0/(2*pi*r));

error_a       = 1; 
error_a_prime = 1;

a_new = a_old;
a_prime_new = a_prime_old;

% Iteration of a and a'
while (abs(error_a) >= 0.00001) || (abs(error_a_prime) >= 0.00001)
    psi = atan( (V/( pdata.big_omega*r)) * ((1+a_new)/(1-a_prime_new)) );
    % Getting the Cl and Cd
    [Cl,Cd] = naca16_509_m06(Beta-psi);
    % Computation of S
    S = (pdata.B*pdata.chord)/(2*pi*r);
    % Computation of lambda1 and lambda2
    lambda1 = (Cl*cos(psi)) - (Cd*sin(psi));
    lambda2 = (Cl*sin(psi)) + (Cd*cos(psi));    
    % Computation of a and a'
    a_new_young       = (1 + a_new      )*(0.5)*S*lambda1/(1-cos(2*psi));
    a_prime_new_young = (1 - a_prime_new)*(0.5)*S*lambda2/(  sin(2*psi));
    % Relaxing the value
    a_new       = (1-pdata.omega)* a_old       + pdata.omega*a_new_young;
    a_prime_new = (1-pdata.omega)* a_prime_old + pdata.omega*a_prime_new_young;
    % Error computation
    error_a       = a_new - a_old;
    error_a_prime = a_prime_new - a_prime_old;
    a_old = a_new;
    a_prime_old = a_prime_new;
end
    a = a_new;
    a_prime = a_prime_new;
end

