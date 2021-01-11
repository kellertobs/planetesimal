% Planetesimal user startup
% valid for 2D rectangular grid
clear; close all

%% setup model domain (user editable)
L           = 500*1e3;                          % length of x domain
D           = 500*1e3;                          % length of z domain
nx          = 51;                               % number of x nodes
nz          = 51;                               % number of z nodes
dx          = L/(nx-1);                         % spacing of x coordinates
dz          = D/(nz-1);                         % spacing of z coordinates

%% setup physical parameters
gz          = 10;                               % vertical/z gravitational constant 
gx          = 0;                                % horizontal/x gravitational constant

%% setup material parameters
%               1mantle 2plume  3air 
rhom        =   [3300   3200    1   ];
etam        =   [1e21   1e20    1e18];

