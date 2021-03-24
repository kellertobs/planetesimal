% Planetesimal user startup
% valid for 2D rectangular grid
clear; 
close all
% profile off
% profile on
% profile off
RunID       = '200x200 imping plume fromm';   % run identifier

%% setup model domain (user editable)
L           = 500*1e3;                          % length of x domain
D           = 500*1e3;                          % length of z domain
nx          = 200;                               % number of real x block nodes
nz          = 200;                               % number of real z block nodes
nx1         = nx+1;                             % number of x gridpoints
nz1         = nz+1;                             % number of z gridpoints
Nx          = nx+2;                             % number of x edgepoints (+ ghost)
Nz          = nz+2;                             % number of z edgepoints (+ ghost)
dx          = L/nx;                             % spacing of x coordinates
dz          = D/nz;                             % spacing of z coordinates
rplume      = 100000;                           % radius of plume (if modelling round plume)   
%% setup physical parameters
gz          = 10;                               % vertical/z gravitational constant 
gx          = 0;                                % horizontal/x gravitational constant
Ambtype    = 'constant'                        % constant ambient background temperature
% Ambtype    = 'linear'                          % linear temperaure profile between top and bottom
% Ambtype     = 'gaussian'



%% setup material parameters
%        Mantle      |1
Rho_mantle  =   3300;   % density
Eta_mantle  =   1e20;   % viscosity  
Alpha_mantle=   3e-5;   % thermal expansion   
Cp_mantle   =   1000;   % Volumetric heat capacity   
Kappa_mantle=   10;     % Thermal conductivity
Hr_mantle   =   2e-8;   % Radiogenic heat production

% stable temperature profile, constant gradient between base and top of
% crust
switch Ambtype
%           bot      ||          top         ||
case 'constant'
T_bot       =   1200;   T_top       =   1200; 
case 'linear'
T_bot       =   1500;   T_top       =   1200;  
case 'gaussian'
T_bot       =   1500;   T_top       =   1200;
end

Rho0        = Rho_mantle; %baseline density = minimum density

%% setup grid parameters // markers now obsolete
% markers per grid block = nxm x nzm
% nxm             = 6;                % number of x markers within each block
% nzm             = 6;                % number of z markers within each block
% dxm             = dx/nxm;           % average spacing between each x marker
% dzm             = dz/nzm;           % average spacing between each z marker
% nxm_all         = nxm*nx;           % total number of x markers in x vector
% nzm_all         = nzm*nz;           % total number of z markers in z vector

%% choose boundary conditions
% Boundary conditions: free slip=-1; No Slip=1
bcleft          = -1;
bcright         = -1;
bctop           = -1;
bcbottom        = -1;
% thermal boundary conditions, for constant boundary conditions
% T_top           = T_air;
% T_bot           = T_plume;

Ra = Rho0*Alpha_mantle*300*(L^3)*gz/Eta_mantle/(Kappa_mantle/Rho0/Cp_mantle)
%% Loop settings
% % =================================================================% %
% choose temperature solver regime %
% Tsolver         = 'implicit';
Tsolver         = 'explicit'; %explicit recommended, implicit currently isn't compatible
% % =================================================================% %
% choose advection regime %
AdvRegime       = 'fromm'
% AdvRegime       = 'first upwind'
% AdvRegime       = 'second upwind'
% AdvRegime       = 'third upwind'
% AdvRegime       = 'flxdiv'
% % =================================================================% %

nt              = 10000;      % number of loop iterations
max_time        = 7e14;
vpratio         = 1/3;      % Weight of averaged velocity for moving markers
dt              = 1e10;     % initial time-stepping (variable within code)
CFL             = 1/4;        % Courant number to limit advection time step
theta           = 0.5;      % 0 = backwards Euler, 0.5 = Crank-Nicholson, 1 = Forward Euler
% dtkoef          = 1.2;      % timestep increment
dTmax           = 50;       % maximum temperature increase
dsubgridt       = 0;        % subgrid for temperature advection, 1=with, 0=without
Restol          = 1e-3;     %residual tolerance

% % =======
% test a maximum time
% max_time = D/2/2e-9;

if ~exist(['../out/',RunID],'dir'); mkdir(['../out/',RunID]); end
addpath('../usr/cbrewer');
cm = cbrewer('seq', 'YlOrRd',30); % colour map
addpath('../src')
run('preloop_initialisation_noair');
run('main_loop');
% profile report