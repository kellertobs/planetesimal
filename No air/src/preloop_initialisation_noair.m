%numerical setup file
%preferably not user editable

%% setup numerical grid
% setup vectors
[x,z,xvx,zvx,xvz,zvz,xp,zp] = vectorsetup(D,L,dx,dz);

% setup meshes

%% setup marker grid
marknum         = nxm_all*nzm_all;  % Total number of markers
xm              = zeros(1,marknum); % Horizontal coordinates, m
zm              = zeros(1,marknum); % Vertical coordinates, m
Etam            = zeros(1,marknum); % Viscosity at M
Alpham          = zeros(1,marknum); % Thermal expansion at M
RhoCpm          = zeros(1,marknum); % Heat capacity at M
Kappam          = zeros(1,marknum); % Thermal conductivity at M
Hrm             = zeros(1,marknum); % Rafiogenic heating at M
Tm              = zeros(1,marknum); % Temperature at marker
Mtype           = zeros(1,marknum);
m               = 1;                % initialise marker counter

%loop to set markers
for jm=1:1:nxm_all
for im=1:1:nzm_all
    % Define marker coordinates
    xm(m)   = dxm/2+(jm-1)*dxm+(rand-0.5)*dxm; % randomly distributed
    zm(m)   = dzm/2+(im-1)*dzm+(rand-0.5)*dzm;
    % Input marker properties
    Tm(m)       = T_top + abs(zm(m))/D*(T_bot-T_top); 
    Etam(m)     = Eta_mantle;    % Viscosity
    Alpham(m)   = Alpha_mantle;  % Thermal expansion
    RhoCpm(m)   = RhoCp_mantle;  % Heat capacity
    Kappam(m)   = Kappa_mantle;  % Thermal conductivity
    Hrm(m)      = Hr_mantle;     % Radiogenic heating
    if zm(m)>D*0.9
        Mtype(m) = 1;
    elseif zm(m)>D*0.8
        Mtype(m) = 2;
    elseif zm(m)>D*0.7
        Mtype(m) = 1;
    elseif zm(m)>D*0.6
        Mtype(m) = 2;
    elseif zm(m)>D*0.5
        Mtype(m) = 1;
    elseif zm(m)>D*0.3
        Mtype(m) = 2;
    elseif zm(m)>D*0.1
        Mtype(m) = 1;
    else
        Mtype(m) = 2;
    end
        

    % Update marker counter
    m=m+1;
end
end

%% initialise arrays
Epsxz           = zeros(Nz,Nx);     % strain rate on the ordinary grid
Sigxz           = zeros(Nz,Nx);     % deviatoric stress on the ordinary grid
Epsxx           = Epsxz;            % strain rate in 
Sigxx           = Sigxz;            % deviatoric stress the middle of grid/pressure nodes
Hs              = Sigxx;            % shear heating, on the pressure nodes
Ha              = Hs;               % adiabatic heating, on pressure nodes
Hr              = Hs;               % radiogenic heating


%% create indexing system
Number  = numsetup(Nz,Nx);  % ordinary grid
