function [x,z,xvx,zvx,xvz,zvz,xp,zp] = vectorsetup(zsize,xsize,dx,dz)
%setup x and y vectors for vx, vz & P
x   = 0:dx:xsize;                   % Horizontal coordinates of basic grid points, m
z   = 0:dz:zsize;                   % Vertical coordinates of basic grid points, m
xvx = 0:dx:xsize;                   % Horizontal coordinates of vx grid points, m
zvx = -dz/2:dz:zsize+dz/2;          % Vertical coordinates of vx grid points, m
xvz = -dx/2:dx:xsize+dx/2;          % Horizontal coordinates of vy grid points, m
zvz = 0:dz:zsize;                   % Vertical coordinates of vy grid points, m
xp  = -dx/2:dx:xsize+dx/2;          % Horizontal coordinates of P grid points, m
zp  = -dz/2:dz:zsize+dz/2;          % Vertical coordinates of P grid points, m
end
