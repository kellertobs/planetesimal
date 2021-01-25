function [Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha] = ...
    Heat_generation(Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,nx,nz,dx,dz,...
    vx_out,vz_out,Eta_out,Eta_mid,Rho_vz,Alpha_mid,T_mid,gz)

for j = 1:1:nx
for i = 1:1:nz
    Epsxz(i,j)  = 1/2 * ((vx_out(i+1,j) - vx_out(i,j))/dz...
                       + (vz_out(i,j+1) - vz_out(i,j))/dx);   %shear strain on normal grid
    Sigxz(i,j)  =       2*Eta_out(i,j)*Epsxz(i,j);          %deviatoric stress on normal grid
end 
end
for j = 2:1:nx %calculate normal stress/strain components on the middle of grid/pressure nodes
for i = 2:1:nz
    Epsxx(i,j)  = (vx_out(i,j) - vx_out(i,j-1))/dx; %Normal strain rate
    Sigxx(i,j)  = 2*Eta_mid(i,j)*Epsxx(i,j); %deviatoric stress
end
end
    
%compute shear heating and adiabatic heating
for j = 2:1:nx 
for i = 2:1:nz
    Hxz         = (Epsxz(i-1,j-1)*Sigxz(i-1,j-1)    + Epsxz(i-1,j)*Sigxz(i-1,j) +...
                   Epsxz(i-1,j)*Sigxz(i-1,j)        + Epsxz(i,j)*Sigxz(i,j))/4; %average of xz products
    Hs(i,j)     = 2*Hxz + 2*Epsxx(i,j)*Sigxx(i,j); %shear heating
    % now compute adiabatic heating
    Ha(i,j)     = (vz_out(i,j) + vz_out(i-1,j))/2 * (Rho_vz(i,j) + Rho_vz(i-1,j))/2 ...
        *T_mid(i,j)*gz*Alpha_mid(i,j);
end
end