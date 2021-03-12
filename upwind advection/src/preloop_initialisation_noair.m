%numerical setup file
%preferably not user editable

%% setup numerical grid
% setup vectors
[x,z,xvx,zvx,xvz,zvz,xp,zp] = vectorsetup(D,L,dx,dz);

% setup meshes
[x2d,z2d] = meshgrid(x,z); %edgepoint grid
[xvx2d,zvx2d] = meshgrid(xvx,zvx); % staggered vx grid
[xvz2d,zvz2d] = meshgrid(xvz,zvz); % staggered vz grid
[xp2d,zp2d] = meshgrid(xp,zp); % centre clock grod

%% setup material grids
% temperature
switch Ambtype
    case 'constant'
        T_mid = zeros(Nz,Nx)+T_top; 
        T_mid(Nz:nz-1,:) = T_bot;
        T_mid(Nz,:) = T_bot+300;
    case 'linear'
        T_mid = T_top + abs(zp2d)./D.*(T_bot-T_top);
    case 'gaussian'
        radius = L/7;
%         T_amp = T_bot-T_top+300;
%         T_amp = 600;
        T_mid = zeros(Nz,Nx)+T_top; 
        T_mid = T_mid + (T_bot-T_top+300).*exp(- (xp2d-L/2).^2./radius.^2 - (zp2d-L/2).^2./radius.^2 );
        T_mid(Nz,:) = T_bot;
end

    Material    = ones(Nz,Nx); % artificial colours for visual representation of deformation
    Material(1:nz/10,:) = 2; Material(0.2*nz:0.3*nz,:) = 2; Material(0.4*nz:0.5*nz,:) = 2;
    Material(0.6*nz:0.7*nz,:) = 2; Material(0.8*nz:0.9*nz,:) = 2;
    Alpha_mid   = zeros(Nz,Nx) + Alpha_mantle;
    Eta_mid     = zeros(Nz,Nx) + Eta_mantle; Eta_out = Eta_mid; % viscosity
    Kappa_mid   = zeros(Nz,Nx) + Kappa_mantle; % thermal conductivity
    Cp_mid      = zeros(Nz,Nx) + Cp_mantle; % heat capacity
    Hr_mid      = zeros(Nz,Nx) + Hr_mantle; % radiogenic heat production

    %initialise staggered grids
    T_vx    = zeros(Nz,Nx); T_vz    = zeros(Nz,Nx);
    k_vx    = zeros(Nz,Nx); k_vz = zeros(Nz,Nx);
    Alpha_vx= zeros(Nz,Nx); Alpha_vz = zeros(Nz,Nx);


%% initialise arrays
Epsxz           = zeros(Nz,Nx);     % strain rate on the ordinary grid
Sigxz           = zeros(Nz,Nx);     % deviatoric stress on the ordinary grid
Epsxx           = Epsxz;            % strain rate in 
Sigxx           = Sigxz;            % deviatoric stress the middle of grid/pressure nodes
Hs              = Sigxx;            % shear heating, on the pressure nodes
Ha              = Hs;               % adiabatic heating, on pressure nodes
Hr              = Hs;               % radiogenic heating
adv_T = 0;     lapl_T = 0;   dMdt = 0;

%% create indexing system
Number  = numsetup(Nz,Nx);  % ordinary grid
