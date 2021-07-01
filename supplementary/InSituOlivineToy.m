% In-situ Olivine melting toy model
% component end-members
% bulk composition
clear all; close all
CMg_r   = 0.3;  CFe_r   = 0.7;  % reference bulk compositions
Tm_Mg   = 2000; Tm_Fe   = 1800; % reference melting points
Phi0    = 0;
T0      = 1700;   P0    = 0;
clap    = 6*1e-8;

perT = 0; perCs = 1; perCf = 1; PhDg = 0.5;

T = linspace(T0,2100,100);
T = (T-Tm_Fe)./(Tm_Mg-Tm_Fe);
Phi = zeros(1,length(T));
figure(100)
plot(Phi+CMg_r,T,'--')
hold on
% for i = 1:1:length(T)
[phi,cs,cf]  =  equilibrium(T,CMg_r,P0,perT,perCs,perCf,clap,PhDg);
Phi = phi;
c_out = phi.*cf +(1-phi).*cs;
plot(cs,T,'or')
hold on
plot(cf,T,'ok')
hold on
plot(phi,T)
hold on
ylabel('T')
xlabel('c_i [%]')
% [phi,cs,cf]  =  equilibrium(T,CMg_r,P0,perT,perCs,perCf,clap,PhDg);
% Phi = phi;
% plot(cs,T,'or')
% hold on
% plot(cf,T,'ok')
% hold on
% plot(phi,T)
% % end
figure()
plot(CMg_r,T)
hold on
plot(c_out,T)
function [f,cs,cf]  =  equilibrium(T,C,P,perT,perCs,perCf,clap,PhDg)
% T ranging from 0-1, T=0 melting point of the lowest end-member at P0,
% T = 1 = highest end-member melting point at P0
% P = P_l + P_s+ P_cmp
% clap = clapeyron slope, assumed
T   = T - P*clap;

cs1 = max(0,min(1,          perCs .*erfc((2+PhDg).*(T-perT)./(1-perT))));
cs2 = max(0,min(1, perCs+(1-perCs).*erfc((2+PhDg).*(T     )./   perT) ));

cs = zeros(size(T));
cs(T>=perT) = cs1(T>=perT);
cs(T< perT) = cs2(T< perT);

cf1 = max(0,min(1,          perCf .*erf(1.5.*(1-T)         ./(1-perT))./erf(1.5)));
cf2 = max(0,min(1, perCf+(1-perCf).*erf(1.5.*(1-T-(1-perT))./(  perT))./erf(1.5)));


cf = zeros(size(T));
cf(T>=perT) = cf1(T>=perT);
cf(T< perT) = cf2(T< perT);

f = max(1e-16,min(1-1e-16, (C-cs) ./ (cf-cs) )); % f = phi
end
