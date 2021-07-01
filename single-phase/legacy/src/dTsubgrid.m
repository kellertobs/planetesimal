function dT = dTsubgrid(nx,nz,nx1,nz1,xp,zp,dx,dz,xm,zm,...
    dsubgridt,T_mid,dT,Tm,dt,marknum,...
    RhoCpm,Kappam)
TSUM        =zeros(nz1,nx1);
RHOCPSUM    =zeros(nz1,nx1);

for m = 1:1:marknum
    % Define i,j indexes for the upper left node
    j       = fix((xm(m)-xp(1))/dx)+1;
    i       = fix((zm(m)-zp(1))/dz)+1;
    if(j<1)
        j   = 1;
    elseif(j>nx)
        j   = nx;
    end
    if(i<1)
        i   = 1;
    elseif(i>nz)
        i   = nz;
    end
    % Compute distances
    dxm1    = xm(m) - xp(j);
    dzm1    = zm(m) - zp(i);
    % Compute weights
    wtmij   = (1-dxm1/dx) * (1-dzm1/dz);
    wtmi1j  = (1-dxm1/dx) * (dzm1/dz);    
    wtmij1  = (dxm1/dx)   * (1-dzm1/dz);
    wtmi1j1 = (dxm1/dx)   * (dzm1/dz);
    % Compute marker-node T difference
    dTm0    = Tm(m) - (T_mid(i,j)  *wtmij   + T_mid(i+1,j)  *wtmi1j+...
                       T_mid(i,j+1)*wtmij1  + T_mid(i+1,j+1)*wtmi1j1);
    % Relax temperature difference
    dTm1    = dTm0*exp(-dsubgridt*Kappam(m)*dt/RhoCpm(m)*(2/dx^2+2/dz^2));
    % Correct marker temperature
    ddTm    =dTm1-dTm0;
    Tm(m)   =Tm(m)+ddTm;
    % Update subgrid diffusion on nodes
    % i,j Node
    TSUM(i,j)           = TSUM(i,j)         + ddTm*RhoCpm(m)*wtmij;
    RHOCPSUM(i,j)       = RHOCPSUM(i,j)     + RhoCpm(m)*wtmij;
    % i+1,j Node
    TSUM(i+1,j)         = TSUM(i+1,j)       + ddTm*RhoCpm(m)*wtmi1j;
    RHOCPSUM(i+1,j)     = RHOCPSUM(i+1,j)   + RhoCpm(m)*wtmi1j;
    % i,j+1 Node
    TSUM(i,j+1)         = TSUM(i,j+1)       + ddTm*RhoCpm(m)*wtmij1;
    RHOCPSUM(i,j+1)     = RHOCPSUM(i,j+1)   + RhoCpm(m)*wtmij1;
    % i+1,j+1 Node
    TSUM(i+1,j+1)       = TSUM(i+1,j+1)     + ddTm*RhoCpm(m)*wtmi1j1;
    RHOCPSUM(i+1,j+1)   = RHOCPSUM(i+1,j+1) + RhoCpm(m)*wtmi1j1;
end
% Compute DTsubgrid
dTsubgrid   = zeros(nz1,nx1);
% P-nodes
for j=1:1:nx1
for i=1:1:nz1
    if(RHOCPSUM(i,j)>0)
        dTsubgrid(i,j)  = TSUM(i,j)/RHOCPSUM(i,j);
    end
end
end
% Correct DT
dT=dT-dTsubgrid;

    
