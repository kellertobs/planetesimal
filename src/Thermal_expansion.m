function [Rho_vx,Rho_vz,Rho_mid] = Thermal_expansion(Material,Rho_mantle,Rho_air,...
    T_mid,Alpha_vx,Alpha_vz,Alpha_mid,T_air,nx1,nz1)

Rho0 = Material;
Rho0(Rho0==1) = Rho_mantle;
Rho0(Rho0==2) = Rho_air;
Tvx = zeros(nz1,nx1); Tvz = zeros(nz1,nx1);
Tvx(:,1:nx1-1) = (T_mid(:,1:end-1) + T_mid(:,2:end)) / 2; % midpoints between each node in a row
Tvz(1:nz1-1,:) = (T_mid(1:end-1,:) + T_mid(2:end,:)) / 2; % midpoints between each node in a collumn

Rho_mid = Rho0.*(1 - Alpha_mid.*(T_mid  -T_air));
Rho_vx  = Rho0.*(1 - Alpha_vx.*(Tvx     -T_air));
Rho_vz  = Rho0.*(1 - Alpha_vz.*(Tvz     -T_air));

