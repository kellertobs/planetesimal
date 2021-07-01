function [Rho_vx,Rho_vz,Rho_mid] = Thermal_expansion(Rho0,...
    T_mid,T_vx,T_vz,Alpha_vx,Alpha_vz,Alpha_mid,T_top)

% Tvx = zeros(nz1,nx1); Tvz = zeros(nz1,nx1);
% Tvx(:,1:nx1-1) = (T_mid(:,1:end-1) + T_mid(:,2:end)) / 2; % midpoints between each node in a row
% Tvz(1:nz1-1,:) = (T_mid(1:end-1,:) + T_mid(2:end,:)) / 2; % midpoints between each node in a collumn

Rho_mid = Rho0.*(1 - Alpha_mid.*(T_mid  -T_top));
Rho_vx  = Rho0.*(1 - Alpha_vx.*(T_vx     -T_top));
Rho_vz  = Rho0.*(1 - Alpha_vz.*(T_vz     -T_top));

