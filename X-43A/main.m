% X-43A
%% Flight condition
z = 30000  ; %[m]
M0 = 6.84 ; 
T_0 =  227  ; %[K]
p_0 = 1171  ; %[Pa]
[p,rho,T,gamma,R,cp] = stdatm(z);
C = sqrt(gamma*R*T);

%%
% By the geometry of the problem, we can deduce theta and an design "beta"
% where we will add the pitch.

%%

options = optimset('Display','off');

%X = [ 2.5]; 
X = linspace(0,5,300);
Pi_c = zeros(length(X),1);
mdot =  zeros(length(X),1);
for i = 1 : length(X);
theta_1 = abs(pi/180*( -2.9196 - X(i)));
theta_2 = abs(pi/180*( -8.986 - X(i)));
beta_1_0  = abs(pi/180*(  -8.21 - X(i)));
beta_2_0  = abs(pi/180*( -16.45 - X(i)));
    
M1 = M0;


theta =theta_1;
betasolve = @(beta) 2*cot(beta)*(M1^2*sin(beta)^2 -1)/(M1^2*(gamma+cos(2*beta))+2 ) - tan(theta);
%betasolve = @(beta) cot(theta)*cot(beta)+1 - (gamma+1)/2 * ( sin(beta)^2 - 1 / M1^2)^(-1);
%betasolve = @(beta) tan(beta)*((gamma+1)*M1^2 / (2*( M1^2 *sin(beta_1)-1 ) ) -1) - cot(theta);

beta_1 = fsolve(betasolve,beta_1_0,options);
%pb_pa = @(Ma) (1 + 2*(gamma)/(gamma+1) * ( Ma^2 * sin(beta_1)^2 -1 ));
pb_pa = @(Ma) (2*gamma* Ma^2 *sin(beta_1)^2 - (gamma-1))/(gamma+1);
%ptb_pta = @(Ma,Mb) (1 + 2*(gamma)/(gamma+1) * ( Ma^2 * sin(beta_1)^2 -1 )) * ((1+ (gamma-1)/2 * Mb^2 )^(gamma/(gamma-1)))/((1+(gamma-1)/2 * Ma^2 )^(gamma/(gamma-1)));
ptb_pta = @(Ma,Mb) ( (gamma+1)/(2*gamma*Ma^2 * sin(beta_1)^2 - (gamma -1)) )^(1/(gamma-1)) * ( ((gamma+1)*Ma^2 * sin(beta_1)^2)/(2+(gamma-1)*Ma^2 * sin(beta_1)^2) )^(gamma/(gamma-1));
%ptb_pta = @(Ma,Mb) ( (gamma+1)*Ma^2 * sin(beta_1)^2/((gamma-1)*Ma^2 * sin(beta_1)^2+2 ) )^(gamma/(gamma-1)) * ( ((gamma+1)*Ma^2 * sin(beta_1)^2)/(2+(gamma-1)*Ma^2 * sin(beta_1)^2) )^(gamma/(gamma-1));

rhob_roha = @(Ma) ( (gamma+1)* Ma^2 *sin(beta_1)^2 ) / ( (gamma-1)* Ma^2 * sin(beta_1)^2 + 2  ) ;  
M2 = 1/(sin(beta_1 - theta_1 )) * sqrt( ( 2+ (gamma-1)*(M1^2)*(sin(beta_1)^2)) / ( 2*gamma*(M1^2)*(sin(beta_1)^2) - (gamma-1)) );

pt2_pt1 = ptb_pta(M1,M2);
p2_p1 = ((1+ (gamma-1)/2 * (M1^2)*(sin(beta_1)^2) )^(gamma/(gamma-1)))/((1+ (gamma-1)/2 * (M2^2)*(sin(beta_1)^2) )^(gamma/(gamma-1))) * pt2_pt1;
%Pip2_p1 = pb_pa(M1); 
rho2_rho1 = rhob_roha(M1);

Mbeta = M2^2;
theta =theta_2;
betasolve = @(beta) 2*cot(beta)*(Mbeta*sin(beta)^2-1)/(Mbeta*(gamma+cos(2*beta))+2 ) - tan(theta);
%betasolve = @(beta) cot(theta)*cot(beta)+1 - (gamma+1)/2 * ( sin(beta)^2 - 1 / M2^2)^(-1);
beta_2 = fsolve(betasolve,beta_2_0,options);
pb_pa = @(Ma) (1 + 2*(gamma)/(gamma+1) * ( Ma^2 * sin(beta_2)^2 -1 ));
%pb_pa = @(Ma) (2*gamma* Ma^2 *sin(beta_2)^2 - (gamma-1))/(gamma+1);
%ptb_pta = @(Ma,Mb) (1 + 2*(gamma)/(gamma+1) * ( Ma^2 * sin(beta_2)^2 -1 )) * ((1+ (gamma-1)/2 * Mb^2 )^(gamma/(gamma-1)))/((1+(gamma-1)/2 * Ma^2 )^(gamma/(gamma-1)));
ptb_pta = @(Ma,Mb) ( (gamma+1)/(2*gamma*Ma^2 * sin(beta_2)^2 - (gamma-1)) )^(1/(gamma-1)) * ( ((gamma+1)*Ma^2 * sin(beta_2)^2)/(2+(gamma-1)*Ma^2 * sin(beta_2)^2) )^(gamma/(gamma-1));

rhob_roha = @(Ma) ( (gamma+1)* Ma^2 *sin(beta_2)^2 ) / ( (gamma-1)* Ma^2 * sin(beta_2)^2 + 2  ) ; 
M3 = 1/(sin(beta_2 - theta_2 )) * sqrt( ( 2+ (gamma-1)*(M2^2)*(sin(beta_2)^2)) / ( 2*gamma*(M2^2)*(sin(beta_2)^2) - (gamma-1)) );

pt3_pt2 = ptb_pta(M2,M3);
p3_p2 = ((1+ (gamma-1)/2 * (M2^2)*(sin(beta_2)^2) )^(gamma/(gamma-1)))/((1+ (gamma-1)/2 * (M3^2)*(sin(beta_2)^2) )^(gamma/(gamma-1))) * pt3_pt2;
%p3_p2 = pb_pa(M2);
rho3_rho2 = rhob_roha(M2);

Pi_c(i) = pt2_pt1 * pt3_pt2 ; 
%tau_c = ( 1+ (gamma-1)/2 * M3^2 )/( 1 + (gamma-1)/2 * M1^2 ) * p3_p2/rho3_rho2 * p2_p1/rho2_rho1;
tau_c = ( 1+ (gamma-1)/2 * M3^2 )/( 1 + (gamma-1)/2 * M1^2 ) * (2*gamma*(M2^2)*(sin(beta_1)^2)-(gamma-1))*((gamma-1)*(M2^2)*(sin(beta_2)^2)+2)/((gamma+1)^2 * (M2^2)*(sin(beta_2)^2)) * (2*gamma*(M1^2)*(sin(beta_1)^2)-(gamma-1))*((gamma-1)*(M1^2)*(sin(beta_1)^2)+2)/((gamma+1)^2 * (M1^2)*(sin(beta_1)^2));

% isentropic efficiency
eta_c = (Pi_c(i).^((gamma-1)/gamma) - 1) / (tau_c-1);
% polytropic efficiency
poly_eff = @(e_p) (( (Pi_c(i))^((gamma-1)/gamma) - 1 ) / ( (Pi_c(i))^((gamma-1)/(gamma*e_p)) - 1 )) - eta_c;
e_c = fsolve(poly_eff,eta_c,options);
% mass flow rate
% mdot = A * rho_fin * M*C
if     beta_1 < 0.1014
    A = 0;
elseif beta_1 > 0.1014 && beta_1 < 0.1433
    A = (1.76*tan(beta_1) - 0.174) * 0.45;
else
    A = 0.08*0.45;
end
mdot(i) = A*rho3_rho2*rho2_rho1*rho*M3*C;
 
% %%
% R2D = 180/pi;
% disp('       M1         M2       M3      Pic     beta1    beta2        eta_c     e_c');
% disp([M0,M2,M3,Pi_c(i),beta_1*R2D,beta_2*R2D,eta_c,e_c]);
end
figure;
plot(X,Pi_c);
xlabel('Angle [degree]');
ylabel('Pi_c');
title('Evolution of Pi_c following the pitch angle');

figure;
plot(X,mdot);
xlabel('Angle [degree]');
ylabel('mdot');
title('Evolution of mdot following the pitch angle');