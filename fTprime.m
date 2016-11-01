function [ dT ] = fTprime(r,V,big_omega)
global pdata;
  [a,a_prime,psi,lambda1,lambda2,S]= sub(r,V,big_omega,0.2,0.2);
  dT = pi*S*(r*r*r)*pdata.ro*(pdata.big_omega^2)*((1-a_prime)^2)*lambda1*(1/((cos(psi))^2));
end

