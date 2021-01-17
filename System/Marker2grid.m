function [Eta_out,Eta_mid,Rho_vz,Rho_vx,Rho_mid,k_vz,k_vx,RhoCp_mid,T_mid,Alpha_mid,Hr]...
    = Marker2grid(marknum,nx,nx1,nz,nz1,dx,dz,xm,zm,...
    x,z,xvx,zvx,xvz,zvz,xp,zp,...
    Etam,Rhom,Kappam,Tm,RhoCpm,Alpham,Hrm)

%% create arrays
%ordinary nodes
ETASUM          = zeros(nz,nx); % viscosity         
WTSUM           = zeros(nz,nx); % weighted average
% vz grid
RHOSUMvz        = zeros(nz,nx1); % Density
KSUMvz          = zeros(nz,nx1); % Thermal conductivity
WTSUMvz         = zeros(nz,nx1); % weighted average
%vx grid
RHOSUMvx        = zeros(nz1,nx); % density
KSUMvx          = zeros(nz1,nx); % Thermal conductivity
WTSUMvx         = zeros(nz1,nx); % weighted average
%Pgrid/middle of blocks grid
TSUMp           = zeros(nz1,nx1); % Temperature
RHOSUMp         = zeros(nz1,nx1); % Density
RHOCPSUMp       = zeros(nz1,nx1); % Heat capacity
ETASUMp         = zeros(nz1,nx1); % Viscosity
ASUMp           = zeros(nz1,nx1); % Thermal expansion
HRSUMp          = zeros(nz1,nx1); % Radiogenic heating
WTSUMp          = zeros(nz1,nx1); % weighted average


% loop through each marker
for m = 1:1:marknum    
    %% ordinary grid   
    % Define i,j indexes for the upper left node
    j           = fix((xm(m)-x(1))/dx)+1;
    i           = fix((zm(m)-z(1))/dz)+1;
    % make sure markers are within the grid
    if(j<1)
        j=1;
    elseif(j>nx-1)
        j=nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>nz-1)
        i=nz-1;
    end
    
    %distance between nodes and markers
    dxm1        = abs(xm(m)-x(j));
    dzm1        = abs(zm(m)-z(i));
    %compute weights of distances from real nodes
    wtij        = (1-dxm1/dx)   *(1-dzm1/dz);
    wtij1       = (dxm1/dx)     *(1-dzm1/dz);
    wti1j       = (1-dxm1/dx)   *(dzm1/dz);
    wti1j1      = dxm1*dzm1/dx/dz;
    % update properties
    % [i,j]
    ETASUM(i,j)     = ETASUM(i,j)       + Etam(m)*wtij;
    WTSUM(i,j)      = WTSUM(i,j)        + wtij;
    % [i,j+1]
    ETASUM(i,j+1)   = ETASUM(i,j+1)     + Etam(m)*wtij1;
    WTSUM(i,j+1)    = WTSUM(i,j+1)      + wtij1;
    % [i+1,j]
    ETASUM(i+1,j)   = ETASUM(i+1,j)     + Etam(m)*wti1j;
    WTSUM(i+1,j)    = WTSUM(i+1,j)      + wti1j;
    % [i+1,j+1]
    ETASUM(i+1,j+1) = ETASUM(i+1,j+1)   + Etam(m)*wti1j1;
    WTSUM(i+1,j+1)  = WTSUM(i+1,j+1)    + wti1j1;
    
    %% staggered vx grid
    % interpolate density onto the vx grid (for the x stokes equation)
    j           = fix((xm(m)-xvx(1))/dx)+1;
    i           = fix((zm(m)-zvx(1))/dz)+1;
    %make sure markers are within grid
    if(j<1)
        j=1;
    elseif(j>nx-1)
        j=nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>nz)
        i=nz;
    end
    %compute distances
    dxm1        = abs(xm(m)-xvx(j));
    dzm1        = abs(zm(m)-zvx(i));
    %compute weights of distances from real nodes
    wtij        = (1-dxm1/dx)   *(1-dzm1/dz);
    wtij1       = (dxm1/dx)     *(1-dzm1/dz);
    wti1j       = (1-dxm1/dx)   *(dzm1/dz);
    wti1j1      = dxm1*dzm1/dx/dz;
    % update properties
    % [i,j]
    RHOSUMvx(i,j)     = RHOSUMvx(i,j)       + Rhom(m)*wtij;
    KSUMvx(i,j)       = KSUMvx(i,j)         + Kappam(m)*wtij;
    WTSUMvx(i,j)      = WTSUMvx(i,j)        + wtij;
    % [i,j+1]
    RHOSUMvx(i,j+1)   = RHOSUMvx(i,j+1)     + Rhom(m)*wtij1;
    KSUMvx(i,j+1)     = KSUMvx(i,j+1)       + Kappam(m)*wtij1;
    WTSUMvx(i,j+1)    = WTSUMvx(i,j+1)      + wtij1;
    % [i+1,j]
    RHOSUMvx(i+1,j)   = RHOSUMvx(i+1,j)     + Rhom(m)*wti1j;
    KSUMvx(i+1,j)     = KSUMvx(i+1,j)       + Kappam(m)*wti1j;
    WTSUMvx(i+1,j)    = WTSUMvx(i+1,j)      + wti1j;
    % [i+1,j+1]
    RHOSUMvx(i+1,j+1) = RHOSUMvx(i+1,j+1)   + Rhom(m)*wti1j1;
    KSUMvx(i+1,j+1)   = KSUMvx(i+1,j+1)     + Kappam(m)*wti1j1;
    WTSUMvx(i+1,j+1)  = WTSUMvx(i+1,j+1)    + wti1j1;
    
    %% staggered vz grid
    % interpolate density onto the vz grid (for the x stokes equation)
    j           = fix((xm(m)-xvz(1))/dx)+1;
    i           = fix((zm(m)-zvz(1))/dz)+1;
    %make sure markers are within grid
    if(j<1)
        j=1;
    elseif(j>nx)
        j=nx;
    end
    if(i<1)
        i=1;
    elseif(i>nz-1)
        i=nz-1;
    end
    %compute distances
    dxm1        = abs(xm(m)-xvz(j));
    dzm1        = abs(zm(m)-zvz(i));
    %compute weights of distances from real nodes
    wtij        = (1-dxm1/dx)   *(1-dzm1/dz);
    wtij1       = (dxm1/dx)     *(1-dzm1/dz);
    wti1j       = (1-dxm1/dx)   *(dzm1/dz);
    wti1j1      = dxm1*dzm1/dx/dz;
    % update properties
    % [i,j]
    RHOSUMvz(i,j)     = RHOSUMvz(i,j)       + Rhom(m)*wtij;
    KSUMvz(i,j)       = KSUMvz(i,j)         + Kappam(m)*wtij;
    WTSUMvz(i,j)      = WTSUMvz(i,j)        + wtij;
    % [i,j+1]
    RHOSUMvz(i,j+1)   = RHOSUMvz(i,j+1)     + Rhom(m)*wtij1;
    KSUMvz(i,j+1)     = KSUMvz(i,j+1)       + Kappam(m)*wtij1;
    WTSUMvz(i,j+1)    = WTSUMvz(i,j+1)      + wtij1;
    % [i+1,j]
    RHOSUMvz(i+1,j)   = RHOSUMvz(i+1,j)     + Rhom(m)*wti1j;
    KSUMvz(i+1,j)     = KSUMvz(i+1,j)       + Kappam(m)*wti1j;
    WTSUMvz(i+1,j)    = WTSUMvz(i+1,j)      + wti1j;
    % [i+1,j+1]
    RHOSUMvz(i+1,j+1) = RHOSUMvz(i+1,j+1)   + Rhom(m)*wti1j1;
    KSUMvz(i+1,j+1)   = KSUMvz(i+1,j+1)     + Kappam(m)*wti1j1;
    WTSUMvz(i+1,j+1)  = WTSUMvz(i+1,j+1)    + wti1j1;
    
    %% staggered P/mid grid
    % interpolate density onto the P grid 
    j           = fix((xm(m)-xp(1))/dx)+1;
    i           = fix((zm(m)-zp(1))/dz)+1;
    %make sure markers are within grid
    if(j<1)
        j=1;
    elseif(j>nx)
        j=nx;
    end
    if(i<1)
        i=1;
    elseif(i>nz)
        i=nz;
    end
    %compute distances
    dxm1        = abs(xm(m)-xp(j));
    dzm1        = abs(zm(m)-zp(i));
    %compute weights of distances from real nodes
    wtij        = (1-dxm1/dx)   *(1-dzm1/dz);
    wtij1       = (dxm1/dx)     *(1-dzm1/dz);
    wti1j       = (1-dxm1/dx)   *(dzm1/dz);
    wti1j1      = dxm1*dzm1/dx/dz;
    % update properties
    % [i,j]
    ETASUMp(i,j)            = ETASUMp(i,j)          + Etam(m)*wtij;
    RHOSUMp(i,j)            = RHOSUMp(i,j)          + Rhom(m)*wtij;
    RHOCPSUMp(i,j)          = RHOCPSUMp(i,j)        + RhoCpm(m)*wtij;
    TSUMp(i,j)              = TSUMp(i,j)            + Tm(m)*wtij;
    ASUMp(i,j)              = ASUMp(i,j)            + Alpham(m)*wtij;
    HRSUMp(i,j)             = HRSUMp(i,j)           + Hrm(m)*wtij;
    WTSUMp(i,j)             = WTSUMp(i,j)           + wtij;
    % [i,j+1]
    ETASUMp(i,j+1)          = ETASUMp(i,j+1)        + Etam(m)*wtij1;
    RHOSUMp(i,j+1)          = RHOSUMp(i,j+1)        + Rhom(m)*wtij1;
    RHOCPSUMp(i,j+1)        = RHOCPSUMp(i,j+1)      + RhoCpm(m)*wtij1;
    TSUMp(i,j+1)            = TSUMp(i,j+1)          + Tm(m)*wtij1;
    ASUMp(i,j+1)            = ASUMp(i,j+1)          + Alpham(m)*wtij1;
    HRSUMp(i,j+1)           = HRSUMp(i,j+1)         + Hrm(m)*wtij1;
    WTSUMp(i,j+1)           = WTSUMp(i,j+1)         + wtij1;
    % [i+1,j]
    ETASUMp(i+1,j)          = ETASUMp(i+1,j)        + Etam(m)*wti1j;
    RHOSUMp(i+1,j)          = RHOSUMp(i+1,j)        + Rhom(m)*wti1j;
    RHOCPSUMp(i+1,j)        = RHOCPSUMp(i+1,j)      + RhoCpm(m)*wti1j;
    TSUMp(i+1,j)            = TSUMp(i+1,j)          + Tm(m)*wti1j;
    ASUMp(i+1,j)            = ASUMp(i+1,j)          + Alpham(m)*wti1j;
    HRSUMp(i+1,j)           = HRSUMp(i+1,j)         + Hrm(m)*wti1j;
    WTSUMp(i+1,j)           = WTSUMp(i+1,j)         + wti1j;
    % [i+1,j+1]
    ETASUMp(i+1,j+1)        = ETASUMp(i+1,j+1)      + Etam(m)*wti1j1;
    RHOSUMp(i+1,j+1)        = RHOSUMp(i+1,j+1)      + Rhom(m)*wti1j1;
    RHOCPSUMp(i+1,j+1)      = RHOCPSUMp(i+1,j+1)    + RhoCpm(m)*wti1j1;
    TSUMp(i+1,j+1)          = TSUMp(i+1,j+1)        + Tm(m)*wti1j1;
    ASUMp(i+1,j+1)          = ASUMp(i+1,j+1)        + Alpham(m)*wti1j1;
    HRSUMp(i+1,j+1)         = HRSUMp(i+1,j+1)       + Hrm(m)*wti1j1;
    WTSUMp(i+1,j+1)         = WTSUMp(i+1,j+1)       + wti1j1;
end

%% interpolate markers into arrays

% ordinary nodes
Eta_out     = ETASUM./WTSUM;        % Viscosity

% vx grid
Rho_vx      = RHOSUMvx./WTSUMvx;    % Density 
k_vx        = KSUMvx./WTSUMvx;      % Thermal conductivity

% vz grid
Rho_vz      = RHOSUMvz./WTSUMvz;    % Density
k_vz        = KSUMvz./WTSUMvz;      % Thermal conductivity

% P/mid grid
Eta_mid     = ETASUMp./WTSUMp;      % Viscosity
Rho_mid     = RHOSUMp./WTSUMp;      % Density
RhoCp_mid   = RHOCPSUMp./WTSUMp;    % Heat capacity
T_mid       = TSUMp./WTSUMp;        % Tempaerature
Alpha_mid   = ASUMp./WTSUMp;        % Thermal expansion
Hr          = HRSUMp./WTSUMp;       % radiogenic heating

