function [dT,Tm] = dTsubgrid(dT,dt,Tm,T_mid,Rho0,T_top,dsubgridt,Nx,Nz,...
    nx1,nz1,xp,zp,xm,zm,dx,dz,Alpham,Cpm,Kappam,marknum,Rho_mid,Cp_mid)
TSUM        =zeros(Nz,Nx);
WTSUM    =zeros(Nz,Nx);
for m = 1:1:marknum
    % Define i,j indexes for the upper left node
    j       = fix((xm(m)-xp(1))/dx)+1;
    i       = fix((zm(m)-zp(1))/dz)+1;
    if(j<2)
        j   = 2;
    elseif(j>nx1)
        j   = nx1;
    end
    if(i<2)
        i   = 2;
    elseif(i>nz1)
        i   = nz1;
    end
    
    % Compute distances
    dxm1    = xm(m) - xp(j);
    dzm1    = zm(m) - zp(i);
    % Compute weights
    wtij    = (1-dxm1/dx) * (1-dzm1/dz);
    wti1j   = (1-dxm1/dx) * (dzm1/dz);    
    wtij1   = (dxm1/dx)   * (1-dzm1/dz);
    wti1j1  = (dxm1/dx)   * (dzm1/dz);
    Rhom    = Rho0.*(1 - Alpham(m).*(Tm(m) -T_top));
    
    dTm0    = Tm(m)-(T_mid(i,j)*wtij+T_mid(i+1,j)*wti1j+...
            T_mid(i,j+1)*wtij1+T_mid(i+1,j+1)*wti1j1);
    tdiff   = Rhom*Cpm(m)/Kappam(m)/((2/dx^2)+(2/dz^2));
    dTm1    = dTm0*exp(-dsubgridt*dt/tdiff);
    ddtkm   = dTm1-dTm0;
    Tm(m)   = Tm(m)+ddtkm;
    
    %i,j
    TSUM(i,j)       = ddtkm*Rho_mid(i,j)*Cp_mid(i,j)*wtij;
    WTSUM(i,j)      = Rho_mid(i,j)*Cp_mid(i,j)*wtij;
    %i+1,j
    TSUM(i+1,j)     = ddtkm*Rho_mid(i+1,j)*Cp_mid(i+1,j)*wti1j;
    WTSUM(i+1,j)    = Rho_mid(i+1,j)*Cp_mid(i+1,j)*wti1j;
    %i,j+1
    TSUM(i,j+1)     = ddtkm*Rho_mid(i,j+1)*Cp_mid(i,j+1)*wtij1;
    WTSUM(i,j+1)    = Rho_mid(i,j+1)*Cp_mid(i,j+1)*wtij1;
    %i+1,j+1
    TSUM(i+1,j+1)   = ddtkm*Rho_mid(i+1,j+1)*Cp_mid(i+1,j+1)*wtij1;
    WTSUM(i+1,j+1)  = Rho_mid(i+1,j+1)*Cp_mid(i+1,j+1)*wtij1;    
end
dT_subgrid  = TSUM./WTSUM;
dT          = dT - dT_subgrid;

for m = 1:1:marknum
    % Define i,j indexes for the upper left node
    j       = fix((xm(m)-xp(1))/dx)+1;
    i       = fix((zm(m)-zp(1))/dz)+1;
    if(j<2)
        j   = 2;
    elseif(j>nx1)
        j   = nx1;
    end
    if(i<2)
        i   = 2;
    elseif(i>nz1)
        i   = nz1;
    end
    
    % Compute distances
    dxm1    = xm(m) - xp(j);
    dzm1    = zm(m) - zp(i);
    % Compute weights
    wtij   = (1-dxm1/dx) * (1-dzm1/dz);
    wti1j  = (1-dxm1/dx) * (dzm1/dz);    
    wtij1  = (dxm1/dx)   * (1-dzm1/dz);
    wti1j1 = (dxm1/dx)   * (dzm1/dz);
    
    Tm(m) = Tm(m) + (dT(i,j)*wtij+dT(i+1,j)*wti1j+...
            dT(i,j+1)*wtij1+dT(i+1,j+1)*wti1j1);
end
%interpolate dT into markers

