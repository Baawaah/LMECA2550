% Philippe Chatelain
% MECA 2550
% Standard atmosphere

function [p,rho,T,gamma,R,cp] = stdatm(z)

T0 = 288.15;% K
g = 9.81; %	m/s2
lambda = -6.5E-3; %	K/m
gamma =  1.4; %
R = 287.0529; % J/kg/K
cp = 1005;% J/kg/K
gammad = -9.8E-3;% K/m
rho0 = 1.225; % kg/m3
p0 = 101325; % Pa


T = (1 + lambda*z/T0) * T0;

p = ((T/T0)^(-g/(R*lambda))) * p0;

rho = ((T/T0)^(-g/(R*lambda) -1 )) * rho0;

p0 / (rho0*T0);

p - rho*R*T;

p0 - rho0*R*T0;
