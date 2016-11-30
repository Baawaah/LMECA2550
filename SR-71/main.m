%%
% Thanh-Son TRAN
%
%%
% Flight Condition
z     = 20000      ; %[m]
p0    = 5475       ; %[Pa]   
T0    = 216        ; %[K]
rho   = 0.088      ; %[kg/m^3]
gamma = 1.4        ;
T_t4  = 1400       ; %[K]
T_t7  = 2300       ; %[K]
LHV   = 43.19*10^6 ; %[K]
c0 = sqrt(gamma*p0/rho);
[~,~,~,~,~,cp] = stdatm(z); 
%% Without losses
M0 = 0.5;
M = linspace(2,3.5,100);
F_m0dot    = zeros(length(M),1);
F_m0dotRAM = zeros(length(M),1);
S          = zeros(length(M),1);
Sram       = zeros(length(M),1);

Pi_c= 20;
for i = 1 : length(M)
        tau_L = T_t4/T0; 
        tau_L_AB = T_t7/T0; 
        tau_r =   1 + (gamma-1)/2 * M(i)^2;
        Pi_r  = ( 1 + (gamma-1)/2 * M(i)^2 )^(gamma/(gamma-1));
        tau_c = Pi_c^((gamma-1)/gamma);
        tau_t = 1 - tau_r/tau_L * (tau_c - 1);
        
        f    = cp*T0/LHV * (tau_r*tau_c) * ( (tau_L/(tau_r*tau_c))*tau_t*((tau_L_AB/(tau_L*tau_t))-1)+(tau_L/(tau_r*tau_c))-1);
        fram = cp*T0/LHV * tau_r*( tau_L_AB/tau_r -1 );
        
        F_m0dot(i)    = c0 *(sqrt( 2/(gamma-1)*(tau_r*tau_c*tau_t -1)*(tau_L/(tau_r*tau_c))*(tau_L_AB/(tau_L*tau_t)) )- M(i) );
        F_m0dotRAM(i) = c0 *(sqrt( 2/(gamma-1)*(tau_r-1)*(tau_L_AB/tau_r) )- M(i) );

        S(i)    = f/F_m0dot(i);
        Sram(i) = fram/F_m0dotRAM(i);
end
figure;
hold on;
plot(M,F_m0dot);
plot(M,F_m0dotRAM);
hold off;

figure;
hold on;
plot(M,S);
plot(M,Sram);
hold off;
%% With losses
%% Without losses
M0 = 0.5;
M = linspace(2,3.5,100);
F_m0dot    = zeros(length(M),1);
F_m0dotRAM = zeros(length(M),1);
S          = zeros(length(M),1);
Sram       = zeros(length(M),1);

Pi_c= 20;
e_c = 0.9;
e_t = 0.9;
for i = 1 : length(M)
        tau_L = T_t4/T0; 
        tau_L_AB = T_t7/T0; 
        tau_r =   1 + (gamma-1)/2 * M(i)^2;
        Pi_r  = ( 1 + (gamma-1)/2 * M(i)^2 )^(gamma/(gamma-1));
        tau_c = Pi_c^((gamma-1)/(e_c*gamma));
        tau_t = 1 - tau_r/tau_L * (tau_c - 1);
        
        f    = cp*T0/LHV * (tau_r*tau_c) * ( (tau_L/(tau_r*tau_c))*tau_t*((tau_L_AB/(tau_L*tau_t))-1)+(tau_L/(tau_r*tau_c))-1);
        fram = cp*T0/LHV * tau_r*( tau_L_AB/tau_r -1 );
        
        F_m0dot(i)    = c0 *(sqrt( 2/(gamma-1)*(tau_r*tau_c*tau_t -1)*(tau_L/(tau_r*tau_c))*(tau_L_AB/(tau_L*tau_t)) )- M(i) );
        F_m0dotRAM(i) = c0 *(sqrt( 2/(gamma-1)*(tau_r-1)*(tau_L_AB/tau_r) )- M(i) );

        S(i)    = f/F_m0dot(i);
        Sram(i) = fram/F_m0dotRAM(i);
end
figure;
hold on;
plot(M,F_m0dot);
plot(M,F_m0dotRAM);
hold off;

figure;
hold on;
plot(M,S);
plot(M,Sram);
hold off;

%% Multicase
% M0 = 0.5;
% M = linspace(0.01,5,100);
% F_m0dot = zeros(length(M),1);
% 
% Pi_c = [20 14 5 2 1];
% for i = 1 : length(M)
%     for n = 1 : length(Pi_c)
%         tau_L = T_t4/T0; 
%         tau_L_AB = T_t7/T0; 
%         tau_r =   1 + (gamma-1)/2 * M(i)^2;
%         Pi_r  = ( 1 + (gamma-1)/2 * M(i)^2 )^(gamma/(gamma-1));
%         tau_c = Pi_c(n)^((gamma-1)/gamma);
%         tau_t = 1 - tau_r/tau_L * (tau_c - 1);
%         F_m0dot(i,n)  = c0 *(sqrt( 2/(gamma-1)*(tau_r*tau_c*tau_t -1)*(tau_L/(tau_r*tau_c))*(tau_L_AB/(tau_L*tau_t)) )- M(i) ); 
%     end
% end
% figure;
% hold on;
% for i = 1 : length(Pi_c)
%     plot(M,F_m0dot(:,i));
%     axis([0 inf 0 inf]);
% end
% hold off;



