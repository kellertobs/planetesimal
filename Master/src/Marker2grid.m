function [Eta_out,Eta_mid,k_vz,k_vx,Cp_mid,T_mid,Alpha_mid,Hr,Material,Alpha_vx,Alpha_vz,T_vx,T_vz]...
    = Marker2grid(marknum,nx,nx1,Nx,nz,nz1,Nz,dx,dz,xm,zm,...
    x,z,xvx,zvx,xvz,zvz,xp,zp,Etam,Kappam,Tm,Cpm,Alpham,Hrm,Mtype)

%% create arrays
%ordinary nodes
ETASUM          = zeros(Nz,Nx); % viscosity         
WTSUM           = zeros(Nz,Nx); % weighted average
% vz grid
KSUMvz          = zeros(Nz,Nx); % Thermal conductivity
ASUMvz          = zeros(Nz,Nx);% Thermal expansion
WTSUMvz         = zeros(Nz,Nx); % weighted average
TSUMvz          = zeros(Nz,Nx);
%vx grid
KSUMvx          = zeros(Nz,Nx); % Thermal conductivity
ASUMvx          = zeros(Nz,Nx);% Thermal expansion
WTSUMvx         = zeros(Nz,Nx); % weighted average
TSUMvx          = zeros(Nz,Nx); 
%Pgrid/middle of blocks grid
TSUMp           = zeros(Nz,Nx); % Temperature
RHOCPSUMp       = zeros(Nz,Nx); % Heat capacity
ETASUMp         = zeros(Nz,Nx); % Viscosity
ASUMp           = zeros(Nz,Nx); % Thermal expansion
HRSUMp          = zeros(Nz,Nx); % Radiogenic heating
WTSUMp          = zeros(Nz,Nx); % weighted average
MARKSUMp        = zeros(Nz,Nx);


% loop through each marker
for m = 1:1:marknum    
    %% edge point grid   
    % Define i,j indexes for the upper left node
    j           = fix((xm(m)-x(1))/dx)+1;
    i           = fix((zm(m)-z(1))/dz)+1;
    % make sure markers are within the grid
    if(j<1)
        j=1;
    elseif(j>nx1)
        j=nx1;
    end
    if(i<1)
        i=1;
    elseif(i>nz1)
        i=nz1;
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
    elseif(j>nx)
        j=nx;
    end
    if(i<1)
        i=1;
    elseif(i>nz1)
        i=nz1
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
    KSUMvx(i,j)       = KSUMvx(i,j)         + Kappam(m)*wtij;
    ASUMvx(i,j)       = ASUMvx(i,j)         + Alpham(m)*wtij;
    WTSUMvx(i,j)      = WTSUMvx(i,j)        + wtij;
    TSUMvx(i,j)       = TSUMvx(i,j)         + Tm(m)*wtij;
    % [i,j+1]
    KSUMvx(i,j+1)     = KSUMvx(i,j+1)       + Kappam(m)*wtij1;
    ASUMvx(i,j+1)     = ASUMvx(i,j+1)       + Alpham(m)*wtij1;
    WTSUMvx(i,j+1)    = WTSUMvx(i,j+1)      + wtij1;
    TSUMvx(i,j+1)     = TSUMvx(i,j+1)       + Tm(m)*wtij1;
    % [i+1,j]
    KSUMvx(i+1,j)     = KSUMvx(i+1,j)       + Kappam(m)*wti1j;
    ASUMvx(i+1,j)     = ASUMvx(i+1,j)       + Alpham(m)*wti1j;
    WTSUMvx(i+1,j)    = WTSUMvx(i+1,j)      + wti1j;
    TSUMvx(i+1,j)     = TSUMvx(i+1,j)       + Tm(m)*wti1j;
    % [i+1,j+1]
    KSUMvx(i+1,j+1)   = KSUMvx(i+1,j+1)     + Kappam(m)*wti1j1;
    ASUMvx(i+1,j+1)   = ASUMvx(i+1,j+1)     + Alpham(m)*wti1j1;
    WTSUMvx(i+1,j+1)  = WTSUMvx(i+1,j+1)    + wti1j1;
    TSUMvx(i+1,j+1)   = TSUMvx(i+1,j+1)     + Tm(m)*wti1j1;
    
    %% staggered vz grid
    % interpolate density onto the vz grid (for the x stokes equation)
    j           = fix((xm(m)-xvz(1))/dx)+1;
    i           = fix((zm(m)-zvz(1))/dz)+1;
    %make sure markers are within grid
    if(j<1)
        j=1;
    elseif(j>nx1)
        j=nx1;
    end
    if(i<1)
        i=1;
    elseif(i>nz)
        i=nz;
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
    KSUMvz(i,j)       = KSUMvz(i,j)         + Kappam(m)*wtij;
    ASUMvz(i,j)       = ASUMvz(i,j)         + Alpham(m)*wtij;
    WTSUMvz(i,j)      = WTSUMvz(i,j)        + wtij;
    TSUMvz(i,j)       = TSUMvz(i,j)         + Tm(m)*wtij;
    % [i,j+1]
    KSUMvz(i,j+1)     = KSUMvz(i,j+1)       + Kappam(m)*wtij1;
    ASUMvz(i,j+1)     = ASUMvz(i,j+1)       + Alpham(m)*wtij1;
    WTSUMvz(i,j+1)    = WTSUMvz(i,j+1)      + wtij1;
    TSUMvz(i,j+1)     = TSUMvz(i,j+1)       + Tm(m)*wtij1;
    % [i+1,j]
    KSUMvz(i+1,j)     = KSUMvz(i+1,j)       + Kappam(m)*wti1j;
    ASUMvz(i+1,j)     = ASUMvz(i+1,j)       + Alpham(m)*wti1j;
    WTSUMvz(i+1,j)    = WTSUMvz(i+1,j)      + wti1j;
    TSUMvz(i+1,j)     = TSUMvz(i+1,j)       + Tm(m)*wti1j;
    % [i+1,j+1]
    KSUMvz(i+1,j+1)   = KSUMvz(i+1,j+1)     + Kappam(m)*wti1j1;
    ASUMvz(i+1,j+1)   = ASUMvz(i+1,j+1)     + Alpham(m)*wti1j1;
    WTSUMvz(i+1,j+1)  = WTSUMvz(i+1,j+1)    + wti1j1;
    TSUMvz(i+1,j+1)   = TSUMvz(i+1,j+1)     + Tm(m)*wti1j1;
    
    %% staggered P/mid grid
    % interpolate density onto the P grid 
    j           = fix((xm(m)-xp(1))/dx)+1;
    i           = fix((zm(m)-zp(1))/dz)+1;
    %make sure markers are within grid
    if(j<2)
        j=2;
    end
    if(j>nx1)
        j=nx1;
    end
    if(i<2)
        i=2;
    end
    if(i>nz1)
        i=nz1;
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
    RHOCPSUMp(i,j)          = RHOCPSUMp(i,j)        + Cpm(m)*wtij;
    TSUMp(i,j)              = TSUMp(i,j)            + Tm(m)*wtij;
    ASUMp(i,j)              = ASUMp(i,j)            + Alpham(m)*wtij;
    HRSUMp(i,j)             = HRSUMp(i,j)           + Hrm(m)*wtij;
    MARKSUMp(i,j)           = MARKSUMp(i,j)         + Mtype(m)*wtij;
    WTSUMp(i,j)             = WTSUMp(i,j)           + wtij;
    % [i,j+1]
    ETASUMp(i,j+1)          = ETASUMp(i,j+1)        + Etam(m)*wtij1;
    RHOCPSUMp(i,j+1)        = RHOCPSUMp(i,j+1)      + Cpm(m)*wtij1;
    TSUMp(i,j+1)            = TSUMp(i,j+1)          + Tm(m)*wtij1;
    ASUMp(i,j+1)            = ASUMp(i,j+1)          + Alpham(m)*wtij1;
    HRSUMp(i,j+1)           = HRSUMp(i,j+1)         + Hrm(m)*wtij1;
    MARKSUMp(i,j+1)         = MARKSUMp(i,j+1)       + Mtype(m)*wtij1;
    WTSUMp(i,j+1)           = WTSUMp(i,j+1)         + wtij1;
    % [i+1,j]
    ETASUMp(i+1,j)          = ETASUMp(i+1,j)        + Etam(m)*wti1j;
    RHOCPSUMp(i+1,j)        = RHOCPSUMp(i+1,j)      + Cpm(m)*wti1j;
    TSUMp(i+1,j)            = TSUMp(i+1,j)          + Tm(m)*wti1j;
    ASUMp(i+1,j)            = ASUMp(i+1,j)          + Alpham(m)*wti1j;
    HRSUMp(i+1,j)           = HRSUMp(i+1,j)         + Hrm(m)*wti1j;
    MARKSUMp(i+1,j)         = MARKSUMp(i+1,j)       + Mtype(m)*wti1j;
    WTSUMp(i+1,j)           = WTSUMp(i+1,j)         + wti1j;
    % [i+1,j+1]
    ETASUMp(i+1,j+1)        = ETASUMp(i+1,j+1)      + Etam(m)*wti1j1;
    RHOCPSUMp(i+1,j+1)      = RHOCPSUMp(i+1,j+1)    + Cpm(m)*wti1j1;
    TSUMp(i+1,j+1)          = TSUMp(i+1,j+1)        + Tm(m)*wti1j1;
    ASUMp(i+1,j+1)          = ASUMp(i+1,j+1)        + Alpham(m)*wti1j1;
    HRSUMp(i+1,j+1)         = HRSUMp(i+1,j+1)       + Hrm(m)*wti1j1;
    MARKSUMp(i+1,j+1)       = MARKSUMp(i+1,j+1)     + Mtype(m)*wti1j1;
    WTSUMp(i+1,j+1)         = WTSUMp(i+1,j+1)       + wti1j1;
end

%% interpolate markers into arrays

% ordinary nodes
Eta_out     = ETASUM./WTSUM;        % Viscosity

% vx grid
k_vx        = KSUMvx./WTSUMvx;      % Thermal conductivity
Alpha_vx    = ASUMvx./WTSUMvx;      % thermal expansion staggered
T_vx        = TSUMvx./WTSUMvx;

% vz grid
k_vz        = KSUMvz./WTSUMvz;      % Thermal conductivity
Alpha_vz    = ASUMvz./WTSUMvz;      % thermal expansion staggered
T_vz        = TSUMvz./WTSUMvz;

% P/mid grid
Eta_mid     = ETASUMp./WTSUMp;      % Viscosity
Cp_mid   = RHOCPSUMp./WTSUMp;    % Heat capacity
T_mid       = TSUMp./WTSUMp;        % Tempaerature
Alpha_mid   = ASUMp./WTSUMp;        % Thermal expansion
Hr          = HRSUMp./WTSUMp;       % radiogenic heating
Material    = round(MARKSUMp./WTSUMp);     % material type

