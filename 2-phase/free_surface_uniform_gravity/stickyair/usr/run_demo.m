%% 'Sticky Air' method of free surface
%
% This current build is only stable in the darcy's flow limit. 
% the fluid migration velocities is sufficiently fast so that the
% time-stepping does not significantly deform the free surface. Some minor
% bugs include high shear heating at the edges of the planetesimal domain,
% largely fixed via the level-set zero heating at space.
% currently the free-slip boundary condition is preferred over the normal
% stress free boundary which shows asymmetry.

% Does not work for the Stoke's flow limit, The free surface deforms to a
% square-shape. 
% Thermally driven solid convection cannot be shown as the velocity fields
% of the space layer overwhelm the thermal expansion induced momentum.

% planetesimal: user control script
clear; close all;


%% set model run options
RUN.ID       =  'thick_planetesimal_demo_NSF';          % run identifier
RUN.plot     =  1;               % switch on to plot live output
RUN.save     =  0;               % switch on to save output files
RUN.nop      =  1;               % output every 'nop' grid steps of transport
RUN.nup      =  1;               % update every 'nup' grid steps of transport
RUN.selfgrav = 1;

%% set model timing
NUM.yr      =  3600*24*365.25;  % seconds per year
NUM.maxstep =  1e4;             % maximum number of time steps
NUM.tend    =  1e8*NUM.yr;      % model stopping time [s]

% [do not modify]
NUM.dt      =  1e3*NUM.yr;      % (initial) time step [s]


%% set model domain
% NUM.D       =  60*1e3;        % length of z domain
% NUM.L       =  60*1e3;        % length of x domain
% NUM.A       =  NUM.D*0.1;     % length of free surface
NUM.D       =  150*1e3;        % length of z domain
NUM.L       =  150*1e3;        % length of x domain
% NUM.A       =  NUM.D*0.1;     % length of free surface

NUM.nz      =  100;             % number of real z block nodes
NUM.nx      =  100;          	% number of real x block nodes

% [do not modify]
NUM.dx      =  NUM.L/NUM.nx;    % spacing of x coordinates
NUM.dz      =  NUM.D/NUM.nz;    % spacing of z coordinates


%% set physicsal parameters
%        solid       ||       liquid         ||     Sticky Air
PHY.Rho0.s  =  3300;    PHY.Rho0.l  =  2600;    PHY.Rho0.a  =   0;       % reference density [kg/m3]
PHY.Eta0.s  =  1e20;    PHY.Eta0.l  =  1e1;     PHY.Eta0.a  =   1e17;    % reference viscosity [Pas]
PHY.aT0     =  3e-5;                            PHY.aT0a    =   0;       % thermal expansivity [1/K]
PHY.kT0     =  1;                               PHY.kT0a    =   300;    % Thermal conductivity [W/m/K]
PHY.Cp0     =  1000;            % Volumetric heat capacity [J/kg/K]
PHY.Hr0     =  1e-8;            % Radiogenic heat productivity [W/m3]
PHY.gz      =  10;              % z-gravity 
PHY.gx      =  0;               % x-gravity
PHY.k0      =  1e-7;            % background permeability

%% set initial condition
SOL.T0      =  298;           	% reference/top potential temperature [C]
SOL.T1      =  298;           	% bottom potential temperature (if different from top) [C]
SOL.Ta      =  298;             % Air temperature
SOL.dT      =  200;           	% temperature perturbation amplitude [C]
SOL.rT      =  NUM.L/6;         % radius of hot plume [m]
SOL.zT      =  0;         % z-position of hot plume [m]
SOL.xT      =  0;         % x-position of hot plume [m]

% SOL.phi0    =  0.01;            % background liquid fraction [vol]
SOL.phi0    =  0.05;            % background liquid fraction [vol]
SOL.dphi    =  0;           	% liquid fraction perturbation amplitude [vol]
SOL.philim  =  1e-4;            % limit liquid fraction for numerical stability

SOL.Ttype   = 'constant';       % constant ambient background temperature
% SOL.Ttype   = 'linear';         % linear temperaure profile between top and bottom
% SOL.Ttype   = 'gaussian';       % Gaussian central plume
% SOL.Ttype   = 'hot bottom';     % hot deep layer, skips the initial T diffusion stage

SOL.PhiType = 'constant';
% SOL.PhiType = 'wet bottom'
% SOL.PhiType = 'gaussian';
% SOL.PhiType = 'dry';


%% set boundary conditions
% Temperature boundary conditions
% SOL.BCTempTop     = 'isothermal';    	% 'isothermal' or 'insulating' bottom boundaries
% SOL.BCTempBot     = 'isothermal';    	% 'isothermal' or 'insulating' bottom boundaries
% SOL.BCTempSides   = 'insulating';    	% 'isothermal' or 'insulating' bottom boundaries

% Velocity boundary conditions: free slip = -1; no slip = 1
SOL.BCleft  = -1;               % left side boundary
SOL.BCright = -1;               % right side boundary
SOL.BCtop   = -1;               % top boundary
SOL.BCbot   = -1;               % bottom boundary


%% set solver options
% advection scheme
NUM.AdvnScheme  = 'fromm';
% NUM.AdvnScheme  = 'first upwind'
% NUM.AdvnScheme  = 'second upwind';
% NUM.AdvnScheme  = 'third upwind'
% NUM.AdvnScheme  = 'flxdiv'

NUM.CFL         = 0.2;   	% Courant number to limit physical time step
NUM.theta     	= 0.5;      % 0 = backwards Euler, 0.5 = Crank-Nicholson, 1 = Forward Euler
NUM.restol    	= 1e-3;     % residual tolerance for nonlinear iterations
NUM.cstab     	= 1e-11;     % stabilising coefficient for P-diagonal


%% start model run [do not modify]

% check thermal Rayleigh number
SOL.Ra = PHY.Rho0.l*PHY.aT0.*(SOL.T1-SOL.T0)*(NUM.L^3)*PHY.gz/PHY.Eta0.l/(PHY.kT0/PHY.Rho0.l/PHY.Cp0);
% check characteristic fluid velocity
SOL.W0 = (PHY.k0*(SOL.phi0+SOL.dphi)^3)*(PHY.Rho0.s-PHY.Rho0.l)*PHY.gz/PHY.Eta0.l/(SOL.phi0+SOL.dphi);
% check compaction length
SOL.Lc = ((PHY.k0*(SOL.phi0+SOL.dphi)^3)*PHY.Eta0.s/(SOL.phi0+SOL.dphi)/PHY.Eta0.l)^(1/2);
if SOL.Ra < 1e3
    disp('WARNING: Rayleigh number too low, no free convection')
end

if NUM.D < SOL.Lc
   disp('WARNING: compaction length too low, no fluid migration') 
end

% create output directory
if ~exist(['../out/',RUN.ID],'dir'); mkdir(['../out/',RUN.ID]); end

% add path to source directory
addpath('../src')
addpath('../src/cbrewer/')

% use color brewer to create colormaps
cm1 =        cbrewer('seq','YlOrRd',30) ; % sequential colour map
cm2 = flipud(cbrewer('div','RdBu'  ,30)); % divergent colour map

% run model
run('main');
