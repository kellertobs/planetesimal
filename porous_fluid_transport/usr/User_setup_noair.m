% Planetesimal user startup
% valid for 2D rectangular grid
clear; 
close all
% profile off
% profile on
% profile off
RunID       = 'Test';   % run identifier

%% setup model domain (user editable)
L           = 500*1e3;                          % length of x domain
D           = 500*1e3;                          % length of z domain
nx          = 100;                               % number of real x block nodes
nz          = 100;                               % number of real z block nodes
nx1         = nx+1;                             % number of x gridpoints
nz1         = nz+1;                             % number of z gridpoints
nx2         = nx+2;                             % number of x edgepoints (+ ghost)
nz2         = nz+2;                             % number of z edgepoints (+ ghost)
dx          = L/nx;                             % spacing of x coordinates
dz          = D/nz;                             % spacing of z coordinates

%% setup physical parameters
gz          = 10;                               % vertical/z gravitational constant 
gx          = 0;                                % horizontal/x gravitational constant

% initial ambent conditions
rplume      = 100000;                           % radius of plume (if modelling round plume)  
% Ambtype     = 'constant'                        % constant ambient background temperature
% Ambtype    = 'linear'                      % linear temperaure profile between top and bottom
% Ambtype    = 'gaussian'                    % Gaussian central plume
Ambtype     = 'hot bottom'                      % a hot lower boundary, skips the initial T diffusion stage

% bottom boundary conditions
BotBoundary     = 'uniform'                     % uniform, constant bottom boundary
% BotBoundary     = 'hell portal'                 % the centre on the bottom boundary hotter than the other

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
case 'hot bottom'
T_bot       =   1500;   T_top       =   1200;
end


%% setup material parameters
%        Mantle      ||
Rho_mantle   = 3300;   % reference mantle density
Eta_mantle   = 1e20;   % reference mantle viscosity  
Alpha_mantle = 3e-5;   % thermal expansion   
Cp_mantle    = 1000;   % Volumetric heat capacity   
Kappa_mantle = 10;     % Thermal conductivity
Hr_mantle    = 2e-8;   % Radiogenic heat production


%% choose boundary conditions
% Boundary conditions: free slip=-1; No Slip=1
bcleft          = -1;
bcright         = -1;
bctop           = -1;
bcbottom        = -1;

Ra = Rho_mantle*Alpha_mantle*(T_bot-T_top)*(L^3)*gz/Eta_mantle/(Kappa_mantle/Rho_mantle/Cp_mantle)
% Ra = 1
if Ra<1e3
    disp('Rayleigh number too low, no free convection')
    return
end


%% Solver settings
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

nt              = 10000;    % number of time steps
max_time        = 7e14;     % model stopping time [s]
vpratio         = 1/3;      % Weight of averaged velocity for moving markers
dt              = 1e10;     % initial time-stepping (variable within code)
CFL             = 1/4;      % Courant number to limit advection time step
theta           = 0.5;      % 0 = backwards Euler, 0.5 = Crank-Nicholson, 1 = Forward Euler
Restol          = 1e-3;     % residual tolerance
cstab           = 1e-7;     % stabilising coefficient for P-diagonal

% % =======
% test a maximum time
% max_time = D/2/2e-9;

% create output directory
if ~exist(['../out/',RunID],'dir'); mkdir(['../out/',RunID]); end

% use color brewer to create colormaps
addpath('../usr/cbrewer');
cm1 = cbrewer('seq','YlOrRd',30); % colour map
cm2 = cbrewer('div','RdBu',30); % colour map

% add path to source directory
addpath('../src')

% run model
run('preloop_initialisation_noair');
run('main_loop');
% profile report