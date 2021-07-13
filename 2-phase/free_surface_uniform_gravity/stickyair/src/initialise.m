% planetesimal: initialise model run

% print initialisation header
fprintf(1,'  ---  initialise model run \n\n');


%% setup numerical grid
% set dimensions of staggered/ghosted 2D grid
NUM.nxC  =  NUM.nx+1;                           % number of corner nodes in x-direction
NUM.nzC  =  NUM.nz+1;                           % number of corner nodes in z-direction
NUM.nxP  =  NUM.nx+2;                           % number of centre nodes in x-direction
NUM.nzP  =  NUM.nz+2;                           % number of centre nodes in z-direction
NUM.nxW  =  NUM.nx+2;                           % number of z-face nodes in x-direction
NUM.nzW  =  NUM.nz+1;                           % number of z-face nodes in z-direction
NUM.nxU  =  NUM.nx+1;                           % number of x-face nodes in x-direction
NUM.nzU  =  NUM.nz+2;                           % number of x-face nodes in z-direction
NUM.NC   =  NUM.nxC*NUM.nzC;                    % total number of corner nodes
NUM.NP   =  NUM.nxP*NUM.nzP;                    % total number of corner nodes
NUM.NW   =  NUM.nxW*NUM.nzW;                    % total number of z-face nodes
NUM.NU   =  NUM.nxU*NUM.nzU;                    % total number of x-face nodes
NUM.NDOF =  NUM.NP+NUM.NW+NUM.NU;               % total number of all degrees of freedum

% set coordinate vectors
NUM.xC   =  0:NUM.dx:NUM.L;                     % Horizontal coordinates of corner nodes [m]
% NUM.zC   =  0:NUM.dz:NUM.D;                     % Vertical   coordinates of corner nodes [m]
NUM.zC   =  -NUM.A:NUM.dz:NUM.D;                     % Vertical   coordinates of corner nodes [m]
NUM.xP   =  -NUM.dx/2:NUM.dx:NUM.L+NUM.dx/2;	% Horizontal coordinates of centre nodes [m]
% NUM.zP   =  -NUM.dz/2:NUM.dz:NUM.D+NUM.dz/2;	% Vertical   coordinates of centre nodes [m]
NUM.zP   =  -NUM.dz/2-NUM.A:NUM.dz:NUM.D+NUM.dz/2;	% Vertical   coordinates of centre nodes [m]
NUM.xW   =  NUM.xP;                             % Horizontal coordinates of z-face nodes [m]
NUM.zW   =  NUM.zC;                             % Vertical   coordinates of z-face nodes [m]
NUM.xU   =  NUM.xC;                             % Horizontal coordinates of x-face nodes [m]
NUM.zU   =  NUM.zP;                             % Vertical   coordinates of x-face nodes [m]

% set 2D coordinate grids
[NUM.XC,NUM.ZC] = meshgrid(NUM.xC,NUM.zC);      % corner nodes grid
[NUM.XP,NUM.ZP] = meshgrid(NUM.xP,NUM.zP);      % centre nodes grid
[NUM.XW,NUM.ZW] = meshgrid(NUM.xW,NUM.zW);      % z-face nodes grid
[NUM.XU,NUM.ZU] = meshgrid(NUM.xU,NUM.zU);      % x-face nodes grid


%% setup mapping arrays
NUM.MapW  =  reshape(1:NUM.NW,NUM.nzW,NUM.nxW);
NUM.MapU  =  reshape(1:NUM.NU,NUM.nzU,NUM.nxU) + NUM.NW;
NUM.MapP  =  reshape(1:NUM.NP,NUM.nzP,NUM.nxP); % + NUM.NW + NUM.NU;


%% setup initial condition for thermo-chemical solution arrays
% set temperature initial condition
pert = -NUM.dx/2.*cos(NUM.XP*2*pi/NUM.D);
switch SOL.Ttype
    case 'constant'     % constant temperature
        SOL.T  = zeros(NUM.nzP,NUM.nxP) + SOL.T1;
        SOL.T(1:11,:) = SOL.Ta;  
    case 'linear'       % linear temperature gradient with depth
        SOL.T  = SOL.T0 + abs(NUM.ZP+pert)./NUM.D.*(SOL.T1-SOL.T0);
    case 'gaussian'     % constant temperature with gaussian plume
        SOL.T = zeros(NUM.nzP,NUM.nxP) + SOL.T0;
        SOL.T = SOL.T + SOL.dT.*exp(-(NUM.XP-SOL.xT).^2./SOL.rT.^2 - (NUM.ZP-SOL.zT).^2./SOL.rT.^2 );
    case 'hot bottom'
        SOL.T = zeros(NUM.nzP,NUM.nxP) + SOL.T0;
        SOL.T = SOL.T + 1./(1+exp(-(NUM.ZP-SOL.zT+pert)./(NUM.D/50))) .* (SOL.T1-SOL.T0);
end

% convert from potential to natural temperature
% SOL.T = SOL.T.*exp(PHY.aT0*(PHY.gz.*(NUM.ZP+pert) + PHY.gx.*NUM.XP)./PHY.Cp0);



% zero gradient boundaries
SOL.T([1 end],:) = SOL.T([2 end-1],:);
SOL.T(:,[1 end]) = SOL.T(:,[2 end-1]);


switch SOL.PhiType
    case 'gaussian'
        % set liquid fraction initial condition (uniform background + gaussian)
        SOL.phi = zeros(NUM.nzP,NUM.nxP) + SOL.phi0;
        SOL.phi = SOL.phi + SOL.dphi.*exp(-(NUM.XP-SOL.xT).^2./SOL.rT.^2 - (NUM.ZP-SOL.zT).^2./SOL.rT.^2 );
        SOL.phi(1:11,:) = SOL.philim;
    case 'constant'
        SOL.phi = zeros(NUM.nzP,NUM.nxP) + SOL.phi0 + SOL.dphi;
end


%% setup velocity-pressure solution arrays

SOL.W.s    = zeros(NUM.nzW,NUM.nxW); SOL.W.l    = zeros(NUM.nzW,NUM.nxW);   % z-velocity on z-face nodes
SOL.U.s    = zeros(NUM.nzU,NUM.nxU); SOL.U.l    = zeros(NUM.nzU,NUM.nxU);   % x-velocity on x-face nodes
SOL.P.s    = zeros(NUM.nzP,NUM.nxP); SOL.P.l    = zeros(NUM.nzP,NUM.nxP);   % Dynamic pressure on centre nodes
SOL.Pt.s   = zeros(NUM.nzP,NUM.nxP); SOL.Pt.l   = zeros(NUM.nzP,NUM.nxP);   % pressure on centre nodes

% project velocities to centre nodes
SOL.UP.s = zeros(NUM.nzP,NUM.nxP);  SOL.UP.l = zeros(NUM.nzP,NUM.nxP);
SOL.WP.s = zeros(NUM.nzP,NUM.nxP);  SOL.WP.l = zeros(NUM.nzP,NUM.nxP);

SOL.UP.s(:,2:end-1) = SOL.U.s(:,1:end-1)+SOL.U.s(:,2:end)./2;
SOL.WP.s(2:end-1,:) = SOL.W.s(1:end-1,:)+SOL.W.s(2:end,:)./2;
SOL.UP.l(:,2:end-1) = SOL.U.l(:,1:end-1)+SOL.U.l(:,2:end)./2;
SOL.WP.l(2:end-1,:) = SOL.W.l(1:end-1,:)+SOL.W.l(2:end,:)./2;

% set segregation velocities and compaction pressures
SOL.Useg    = SOL.U.l;   SOL.Wseg   = SOL.W.l;  SOL.Pcmp = SOL.P.s;
SOL.UP.seg  = SOL.UP.s;  SOL.WP.seg = SOL.WP.s;

%% setup material property arrays
MAT.Rho.s	= zeros(NUM.nzP,NUM.nxP) + PHY.Rho0.s;  	  % solid density
MAT.Rho.s	= zeros(NUM.nzP,NUM.nxP) + PHY.Rho0.l;        % liquid density
MAT.Eta.s	= zeros(NUM.nzP,NUM.nxP) + PHY.Eta0.s;   	  % viscosity
MAT.Eta.l	= zeros(NUM.nzP,NUM.nxP) + PHY.Eta0.l;   	  % liquid viscosity
MAT.aT	= zeros(NUM.nzP,NUM.nxP) + PHY.aT0;               % thermal expansivity
MAT.kT	= zeros(NUM.nzP,NUM.nxP) + PHY.kT0;               % thermal conductivity
MAT.Cp	= zeros(NUM.nzP,NUM.nxP) + PHY.Cp0;   	          % heat capacity
MAT.k   = zeros(NUM.nzP,NUM.nxP) + PHY.k0 .* SOL.phi0^3;  % permeability

% set sticky air properties
MAT.Rho.s(1:11,:)   = PHY.Rho0.a;
MAT.Eta.s(1:11,:)   = PHY.Eta0.a;
MAT.aT(1:11,:)      = PHY.aT0a;
MAT.kT(1:11,:)      = PHY.kT0a;

%% setup deformation property arrays
DEF.ups = zeros(NUM.nzP,NUM.nxP);               % velocity divergence on centre nodes
DEF.exx	= zeros(NUM.nzP,NUM.nxP);               % x-normal strain rate on centre nodes
DEF.ezz = zeros(NUM.nzP,NUM.nxP);               % z-normal strain rate on centre nodes
DEF.exz = zeros(NUM.nzC,NUM.nxC);               % xz-shear strain rate on corner nodes
DEF.eII = zeros(NUM.nzP,NUM.nxP);               % strain rate magnitude on centre nodes
DEF.txx = zeros(NUM.nzP,NUM.nxP);               % x-normal stress on centre nodes
DEF.tzz = zeros(NUM.nzP,NUM.nxP);               % z-normal stress on centre nodes
DEF.txz = zeros(NUM.nzC,NUM.nxC);               % xz-shear stress on corner nodes
DEF.tII = zeros(NUM.nzP,NUM.nxP);               % stress magnitude on centre nodes


%% setup heating rates
SOL.dTdt = zeros(NUM.nz ,NUM.nx );              % temperature rate of change
SOL.Hs   = zeros(NUM.nzP,NUM.nxP);          	% shear heating rate
SOL.Ha   = zeros(NUM.nzP,NUM.nxP);          	% adiabatic heating rate
MAT.Hr	 = zeros(NUM.nzP,NUM.nxP) + PHY.Hr0;	% radiogenic heating


%% update nonlinear material properties
up2date;

SOL.Lc = max(max((MAT.k.*MAT.Eta.s./MAT.Eta.l./SOL.phi).^(1/2)));
%% output initial condition
RUN.frame = 0;     % initialise output frame count
NUM.time  = 0;      % initialise time count
NUM.step  = 0;      % initialise time step count
output;
