function [ dQ ] = fQprime(r,V,big_omega)
global pdata;
  [a,a_prime,psi,lambda1,lambda2,S]= sub(r,V,big_omega,0.2,0.2);
  dQ = pi*S*(r^4)*ro*(big_omega^2)*((1-a_prime)^2)*lambda2*(1/((cos(psi))^2))*dr;
end

