% planetesimal: user control script
clear; close all;


%% set model run options
RUN.ID      =  'demo 3';          % run identifier
RUN.plot    =  1;               % switch on to plot live output
RUN.save    =  0;               % switch on to save output files
RUN.nop     =  1;               % output every 'nop' grid steps of transport
RUN.nup     =  1;               % update every 'nup' grid steps of transport


%% set model timing
NUM.yr      =  3600*24*365.25;  % seconds per year
NUM.maxstep =  1e4;             % maximum number of time steps
NUM.tend    =  1e8*NUM.yr;      % model stopping time [s]

% [do not modify]
NUM.dt      =  1e3*NUM.yr;      % (initial) time step [s]


%% set model domain
NUM.D       =  60*1e3;        % length of z domain
NUM.L       =  60*1e3;        % length of x domain
% NUM.D       =  100*1e3;        % length of z domain
% NUM.L       =  100*1e3;        % length of x domain
NUM.nz      =  100;             % number of real z block nodes
NUM.nx      =  100;          	% number of real x block nodes

% [do not modify]
NUM.dx      =  NUM.L/NUM.nx;    % spacing of x coordinates
NUM.dz      =  NUM.D/NUM.nz;    % spacing of z coordinates


%% set physicsal parameters
%        solid       ||       liquid
PHY.Rho0.s  =  3300;    PHY.Rho0.l  =  2600;        % reference density [kg/m3]
PHY.Eta0.s  =  1e20;    PHY.Eta0.l  =  1e2;         % reference viscosity [Pas]
PHY.aT0     =  3e-5;            % thermal expansivity [1/K]
PHY.kT0     =  10;              % Thermal conductivity [W/m/K]
PHY.Cp0     =  1000;            % Volumetric heat capacity [J/kg/K]
PHY.Hr0     =  1e-6;            % Radiogenic heat productivity [W/m3]
PHY.gz      =  10;              % z-gravity 
PHY.gx      =  0;               % x-gravity
PHY.k0      =  1e-7;            % background permeability

%% set initial condition
SOL.T0      =  100;           	% reference/top potential temperature [C]
SOL.T1      =  2000;           	% bottom potential temperature (if different from top) [C]
SOL.dT      =  100;           	% temperature perturbation amplitude [C]
SOL.rT      =  NUM.L/5;         % radius of hot plume [m]
SOL.zT      =  NUM.D/2;         % z-position of hot plume [m]
SOL.xT      =  NUM.L/2;         % x-position of hot plume [m]

% SOL.phi0    =  0.01;            % background liquid fraction [vol]
SOL.phi0    =  0.01;            % background liquid fraction [vol]
SOL.dphi    =  0.01;           	% liquid fraction perturbation amplitude [vol]
SOL.philim  =  1e-4;            % limit liquid fraction for numerical stability

SOL.Ttype   = 'constant';       % constant ambient background temperature
% SOL.Ttype   = 'linear';         % linear temperaure profile between top and bottom
% SOL.Ttype   = 'gaussian';       % Gaussian central plume
% SOL.Ttype   = 'hot bottom';     % hot deep layer, skips the initial T diffusion stage

% SOL.PhiType = 'constant';
SOL.PhiType = 'gaussian';


%% set boundary conditions
% Temperature boundary conditions
SOL.BCTempTop     = 'isothermal';    	% 'isothermal' or 'insulating' bottom boundaries
SOL.BCTempBot     = 'isothermal';    	% 'isothermal' or 'insulating' bottom boundaries
SOL.BCTempSides   = 'insulating';    	% 'isothermal' or 'insulating' bottom boundaries

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

NUM.CFL         = 0.5;   	% Courant number to limit physical time step
NUM.theta     	= 0.5;      % 0 = backwards Euler, 0.5 = Crank-Nicholson, 1 = Forward Euler
NUM.restol    	= 1e-3;     % residual tolerance for nonlinear iterations
NUM.cstab     	= 1e-11;     % stabilising coefficient for P-diagonal


%% start model run [do not modify]

% check thermal Rayleigh number
SOL.Ra = PHY.Rho0.l*PHY.aT0*(SOL.T1-SOL.T0)*(NUM.L^3)*PHY.gz/PHY.Eta0.l/(PHY.kT0/PHY.Rho0.l/PHY.Cp0);
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
