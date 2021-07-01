function [f,cs,cf,C]  =  equilibrium(T,C,P,perT,perCs,perCf,clap,PhDg)
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
f(T>=1) = 1; % force melt fraction to go to full melt , otherwise numerically unstable at T>1 
% C(T>=1)
end