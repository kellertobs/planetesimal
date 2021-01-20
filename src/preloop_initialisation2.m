%numerical setup file
%preferably not user editable

%% setup numerical grid
[x,z,xvx,zvx,xvz,zvz,xp,zp] =   vectorsetup(D+Air,L,dx,dz); % get vectors
[X,Z]                       =   meshgrid(x,z);              % get meshgrid

%% setup marker grid
marknum         = nxm_all*nzm_all;  % Total number of markers
xm              = zeros(1,marknum); % Horizontal coordinates, m
zm              = zeros(1,marknum); % Vertical coordinates, m
Rhom            = zeros(1,marknum); % density at M
Etam            = zeros(1,marknum); % Viscosity at M
Alpham          = zeros(1,marknum); % Thermal expansion at M
RhoCpm          = zeros(1,marknum); % Heat capacity at M
Kappam          = zeros(1,marknum); % Thermal conductivity at M
Hrm             = zeros(1,marknum); % Rafiogenic heating at M
Tm              = zeros(1,marknum); % Temperature at marker
m               = 1;                % initialise marker counter

%loop to set markers
for jm=1:1:nxm_all
for im=1:1:nzm_all
    % Define marker coordinates
    xm(m)   = dxm/2+(jm-1)*dxm+(rand-0.5)*dxm; % randomly distributed
    zm(m)   = dzm/2+(im-1)*dzm+(rand-0.5)*dzm;
    % Input marker properties
    if(zm(m)>Air)
        Tm(m)       = T_top + abs(zm(m)-Air)/D*(T_bot-T_top);
        Rhom(m)     = Rho_mantle;    % Density 
        Etam(m)     = Eta_mantle;    % Viscosity
        Alpham(m)   = Alpha_mantle;  % Thermal expansion
        RhoCpm(m)   = RhoCp_mantle;  % Heat capacity
        Kappam(m)   = Kappa_mantle;  % Thermal conductivity
        Hrm(m)      = Hr_mantle;     % Radiogenic heating
    else
        Tm(m)       = T_air;
        Rhom(m)     = Rho_air;       % Density 
        Etam(m)     = Eta_air;       % Viscosity
        Alpham(m)   = Alpha_air;     % Thermal expansion
        RhoCpm(m)   = RhoCp_air;     % Heat capacity
        Kappam(m)   = Kappa_air;     % Thermal conductivity
        Hrm(m)      = Hr_air;        % Radiogenic heating
    end
    % Update marker counter
    m=m+1;
end
end

%% initialise arrays
Epsxz           = zeros(nz,nx);     % strain rate on the ordinary grid
Sigxz           = zeros(nz,nx);     % deviatoric stress on the ordinary grid
Epsxx           = Epsxz;            % strain rate in 
Sigxx           = Sigxz;            % deviatoric stress the middle of grid/pressure nodes
Hs              = Sigxx;            % shear heating, on the pressure nodes
Ha              = Hs;               % adiabatic heating, on pressure nodes
Hr              = Hs;               % radiogenic heating


%% create indexing system
Number  = numsetup(nz1,nx1);% ordinary grid
indP    = Number.*3;        % A matrix indexing of P
indvx   = indP-2;           % A matrix indexing of vx grid
indvz   = indP-1;           % A matrix indexing of vz grid