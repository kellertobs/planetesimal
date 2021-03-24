function [gx,gz] = Poisson_gravity(xp,zp,dx,dz,nx1,nz1,nx,nz,Number,xsize,zsize,...
    G,K,Rho_mid,gx,gz)

PHI = zeros(nz1,nx1);
ind     = Number;
% setup matrix of coefficients
A       = sparse(nx1*nz1,nx1*nz1);
RHS     = zeros(nx1*nz1);

for j = 1:1:nx1
for i = 1:1:nz1
    Rxz     = ((xp(j)-xsize/2)^2+(zp(i)-zsize/2)^2)^0.5;    % radial position of current node
    
    % set boundary conditions
    if (Rxz>xsize/2 || j == 1 || j == nx1 || i == 1 || j == nz1)
        %PHI = 0
        A(ind(i,j),ind(i,j))    = 1;
        RHS(ind(i,j))           = 0;
    else
    
    % internal points
    % d^2PHI/dx^2 - d^2PHI/dz^2 = 4KPiGRho
    % fill A matrix coeffs
    A(ind(i,j),ind(i,j-1))      = 1/dx^2;       % left of current node
    A(ind(i,j),ind(i,j+1))      = 1/dx^2;       % right of current node
    A(ind(i,j),ind(i-1,j))      = 1/dz^2;       % above current node
    A(ind(i,j),ind(i+1,j))      = 1/dz^2;       % below current node
    A(ind(i,j),ind(i,j))        = - 2/dx^2 ...
                                  - 2/dz^2;     % current node
    
    % Fill RHS vector
    RHS(ind(i,j))               = 4*K*pi*G*Rho_mid(i,j);
    end
end
end

% solve matrices to output PHI
c = A\RHS;

for j = 1:1:nx1
for i = 1:1:nz1
    PHI(i,j)    = c(ind(i,j));
end
end 
% solve gx and gz vectors
for j = 1:1:nx
for i = 1:1:nz
    gx(i,j)     = (PHI(i,j+1) - PHI(i,j))/dx;
    gz(i,j)     = (PHI(i+1,j) - PHI(i,j))/dz;                       
end
end
        
