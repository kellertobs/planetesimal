function [P_out,vx_out,vz_out] = stokes_continuity(nx,nz,nx1,nz1,...
    indvx,indvz,indP,dx,dz,Pscale,...
    Eta_out,Eta_mid,Rho_vx,Rho_vz,gx,gz,dt,bctop,bcbottom,bcleft,bcright)

N = nx1*nz1*3; %total number of unknowns to solve (vx,vz,P for each node)
%We solve implicitly in the form A*c = RHS, where c is the solution
%implemented a sticky air solver on the vz grid
A = sparse(N,N);
RHS = zeros(N,1);

for j = 1:1:nx1
for i = 1:1:nz1
    %% x-Stokes eq. ETA*(d2Vx/dx^2+d2Vx/dy^2)-dP/dx=0
    % solve x equation
    
    % boundary conditions of vx
    if(i==1 || i==nz1 || j==1 || j==nx || j==nx1)
            A(indvx(i,j),indvx(i,j))    = 1; % A matrix coefficient
            RHS(indvx(i,j))             = 0; % RHS
        %Top boundary
        if (i==1 && j>1 && j<nx)
            A(indvx(i,j),indvx(i+1,j))  = bctop; %only solve for the bottom of the top boundary
        end
        %Bottom boundary    
        if (i==nz1 && j>1 && j<nx)
            A(indvx(i,j),indvx(i-1,j))  = bcbottom; % above the bottom boundary
        end
    % now solve internal points on the real grid
    else
    % A matrix coefficients
    Eta1    = Eta_out(i-1,j);
    Eta2    = Eta_out(i,j);
    EtaP1   = Eta_mid(i,j);
    EtaP2   = Eta_mid(i,j+1);
    drhodx  = (Rho_vx(i,j+1) - Rho_vx(i,j-1))/2/dx;
    drhodz  = (Rho_vx(i+1,j) - Rho_vx(i-1,j))/2/dz;
    
    
    A(indvx(i,j),indvx(i,j-1))      = 2*EtaP1/dx^2;                 % vx left of current node
    A(indvx(i,j),indvx(i-1,j))      = Eta1/dz^2;                    % vx  above current node
    A(indvx(i,j),indvx(i,j))        = -2*(EtaP1+EtaP2)/dx^2-...
                                        (Eta1+Eta2)/dz^2 -...
                                        drhodx*gx*dt;               % vx current node
    A(indvx(i,j),indvx(i+1,j))      = Eta2/dz^2;                    % vx below current node
    A(indvx(i,j),indvx(i,j+1))      = 2*EtaP2/dx^2;                 % vx right of current node
    
    A(indvx(i,j),indvz(i,j))        = -Eta2/dx/dz   -drhodz*gx*dt;  % vz bottomleft
    A(indvx(i,j),indvz(i,j+1))      = Eta2/dx/dz    -drhodz*gx*dt;  % vz bottomright      
    A(indvx(i,j),indvz(i-1,j))      = Eta1/dx/dz    -drhodz*gx*dt;  % vz topleft
    A(indvx(i,j),indvz(i-1,j+1))    = -Eta1/dx/dz   -drhodz*gx*dt;  % vz topright     
    A(indvx(i,j),indP(i,j))         = Pscale/dx;                    % P1; current node
    A(indvx(i,j),indP(i,j+1))       = -Pscale/dx;                   % P2; right of current node
    % RHS
    RHS(indvx(i,j))                 = Rho_vx(i,j)*gx;               % x direction gravity
    end
    
    %% z-Stokes eq. ETA*(d2Vy/dx^2+d2Vy/dy^2)-dP/dy=-RHO*gy
    % solve z equation
    % boundary conditions of vz
    if(j==1 || j==nx1 || i==1 || i==nz || i==nz1)
            A(indvz(i,j),indvz(i,j))    = 1;        % A matrix coefficient
            RHS(indvz(i,j))             = 0;        % RHS
        %left boundary
        if (j==1 && i>1 && i<nz)
            A(indvz(i,j),indvz(i,j+1))  = bcleft;   %solve for right of the leftmost bodes
        end
        %right boundary    
        if (j==nx1 && i>1 && i<nz)
            A(indvz(i,j),indvz(i,j-1))  = bcright;  % above the bottom boundary
        end
    % solve internal points    
    else
        
        Eta1    = Eta_out(i,j-1);
        Eta2    = Eta_out(i,j);
        EtaP1   = Eta_mid(i,j);
        EtaP2   = Eta_mid(i+1,j);
        drhodx  = (Rho_vz(i,j+1) - Rho_vz(i,j-1))/2/dx;
        drhodz  = (Rho_vz(i+1,j) - Rho_vz(i-1,j))/2/dz;
        
    % A matrix coefficients
    A(indvz(i,j),indvz(i,j-1))      = Eta1/dx^2;                    % vx1 left of current node
    A(indvz(i,j),indvz(i-1,j))      = 2*EtaP1/dz^2;                 % vx2 above current node
    A(indvz(i,j),indvz(i,j))        = -2*(EtaP1+EtaP2)/dz^2-...
                                        (Eta1+Eta2)/dx^2-...
                                        drhodz*gz*dt;               % vx3 current node
    A(indvz(i,j),indvz(i+1,j))      = 2*EtaP2/(dz^2);               % vx4 below current node
    A(indvz(i,j),indvz(i,j+1))      = Eta2/(dx^2);                  % vx5 right of current node
    A(indvz(i,j),indvx(i,j))        = -Eta2/dx/dz   -drhodx*gz*dt/4;% topright
    A(indvz(i,j),indvx(i+1,j))      = Eta2/dx/dz    -drhodx*gz*dt/4;% bottomright
    A(indvz(i,j),indvx(i,j-1))      = Eta1/dx/dz    -drhodx*gz*dt/4;% topleft
    A(indvz(i,j),indvx(i+1,j-1))    = -Eta1/dx/dz   -drhodx*gz*dt/4;% bottomleft
    
    A(indvz(i,j),indP(i,j))         = Pscale/dz;                    %P1; current node
    A(indvz(i,j),indP(i+1,j))       = -Pscale/dz;                   %P2; bottom of current node
    % RHS
    RHS(indvz(i,j))                 = -gz*Rho_vz(i,j);
    end
    
    %% P-Stokes eq. dVx/dx+dVy/dy=0
    % boundary points
    % P equation External points
    if(i==1 || j==1 || i==nz1 || j==nx1 ||...
            (i==2 && j==2))
        % Boundary Condition
        % 1*P=0
            A(indP(i,j),indP(i,j))  = 1;        % Left part
            RHS(indP(i,j))          = 0;        % Right part
         % Real BC
        if(i==2 && j==2)
            A(indP(i,j),indP(i,j))  = 1*Pscale; % Left part
            RHS(indP(i,j))          = 1e+9;     % Right part   
        end
    % now solve internal points    
    else
    % A matrix coefficients    
    A(indP(i,j),indvx(i,j-1))       = -1/dx;    % left of current node
    A(indP(i,j),indvx(i,j))         = 1/dx;     % current node
    A(indP(i,j),indvz(i-1,j))       = -1/dz;    % above current node
    A(indP(i,j),indvz(i,j))         = 1/dz;     % below current node
    
    % RHS
    RHS(indP(i,j)) = 0;
    end       
end 
end

%% Solve stokes matrix and convert output vector to matrices
c = A\RHS; %get solution vector
%extrapolate into individual matrices
for j = 1:1:nx1
for i = 1:1:nz1
    
    P_out(i,j)  = c(indP(i,j))*Pscale;  %output pressure
    vx_out(i,j) = c(indvx(i,j));        %output vx
    vz_out(i,j) = c(indvz(i,j));        %output vz
    
end
end

