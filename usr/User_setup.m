% Planetesimal user startup
% valid for 2D rectangular grid
clear; 
close all

%% setup model domain (user editable)
L           = 500*1e3;                          % length of x domain
D           = 500*1e3;                          % length of z domain
Air         = 100*1e3;                          % length of sticky air, z domain
nx          = 51;                               % number of x nodes
nz          = 61;                               % number of z nodes
nx1         = nx+1;                             % x vector plus ghost nodes
nz1         = nz+1;                             % z vector plus ghost nodes
dx          = L/(nx-1);                         % spacing of x coordinates
dz          = (D+Air)/(nz-1);                   % spacing of z coordinates
rplume      = 100000;                           % radius of plume (if modelling round plume)   
%% setup physical parameters
gz          = 10;                               % vertical/z gravitational constant 
gx          = 0;                                % horizontal/x gravitational constant
% model_type  = 'round plume';                          % model type, 'layer' or 'round plume'
model_type  = 'layer';  

%% setup material parameters
%        Mantle      ||         Plume        ||         Air
Rho_mantle  =   3300;   Rho_plume   =   3200;   Rho_air     =   1;      % density
Eta_mantle  =   1e21;   Eta_plume   =   1e20;   Eta_air     =   1e18;   % viscosity  
T_mantle    =   1500;   T_plume     =   1800;   T_air       =   273;    % Temperture
Alpha_mantle=   3e-5;   Alpha_plume =   2e-5;   Alpha_air   =   0;      % thermal expansion   
RhoCp_mantle=   3.3e6;  RhoCp_plume =   3.2e6;  RhoCp_air   =   3.3e6;  % Volumetric heat capacity   
Kappa_mantle=   3;      Kappa_plume =   2;      Kappa_air   =   0.1;   % Thermal conductivity
Hr_mantle   =   2e-8;   Hr_plume    =   3e-8;   Hr_air      =   0;      % Radiogenic heat production

%% setup grid parameters
% markers per grid block = nxm x nzm
nxm             = 5;                % number of x markers within each block
nzm             = 5;                % number of z markers within each block
dxm             = dx/nxm;           % spacing between each x marker
dzm             = dz/nzm;           % spacing between each z marker
nxm_all         = nxm*(nx-1);       % total number of x markers in x vector
nzm_all         = nzm*(nz-1);       % total number of z markers in z vector

%% choose boundary conditions
% Boundary conditions: free slip=-1; No Slip=1
bcleft          = -1;
bcright         = -1;
bctop           = -1;
bcbottom        = -1;
% thermal boundary conditions, for constant boundary conditions
T_top           = T_air;
T_bot           = T_mantle;

%% Loop settings
nt              = 200;      % number of loop iterations
vpratio         = 1/3;      % Weight of averaged velocity for moving markers
dt              = 1e10;     % initial time-stepping (variable within code)
dtkoef          = 1.2;      % timestep increment
dxzmax          = 0.5;      % maximum advected movement
dTmax           = 20;       % maximum temperature increase
dsubgridt       = 1;        % subgrid for temperature advection, 1=with, 0=without