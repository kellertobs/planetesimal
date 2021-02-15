function [P_out,vx_out,vz_out,vx_mid,vz_mid] = stokes_continuity_noair(nx,nz,Nx,Nz,...
    nx1,nz1,dx,dz,Pscale,...
    Eta_out,Eta_mid,Rho_vx,Rho_vz,gx,gz,bctop,bcbottom,bcleft,bcright)

NP      = Nx*Nz; % total number of P nodes to solve + ghost nodes
NU      = Nx*Nz; % total number of vx nodes to solve + ghost nodes
NW      = Nx*Nz; % total number of vz nodes to solve + ghost nodes
N_all   = NP+NU+NW;

%indexing of unknowns
indvx   = reshape(1:2:(NU+NW),Nz,Nx);
indvz   = reshape(2:2:(NW+NU),Nz,Nx);
indP    = reshape(1:NP,Nz,Nx) + NU + NW;

% setup A matrix and RHS vector
II  = [];
JJ  = [];
AA  = [];
RR  = [];

for j = 1:1:Nx
for i = 1:1:Nz
    %% x-Stokes eq. ETA*(d2Vx/dx^2+d2Vx/dy^2)-dP/dx=0
    % solve x equation
    
    % boundary conditions of vx
    if(i==1 || i==Nz || j==1 || j==nx1 || j==Nx)
            A(indvx(i,j),indvx(i,j))    = 1; % A matrix coefficient
            RHS(indvx(i,j))             = 0; % RHS
        %Top boundary
        if (i==1 && j>1 && j<nx1)
            A(indvx(i,j),indvx(i+1,j))  = bctop; %only solve for the bottom of the top boundary
        end
        %Bottom boundary    
        if (i==Nz && j>1 && j<nx1)
            A(indvx(i,j),indvx(i-1,j))  = bcbottom; % above the bottom boundary
        end
    % now solve internal points on the real grid
    else
    % A matrix coefficients
    Eta1    = Eta_out(i-1,j);
    Eta2    = Eta_out(i,j);
    EtaP1   = Eta_mid(i,j);
    EtaP2   = Eta_mid(i,j+1);
    
    II = [II,indvx(i,j)]; JJ = [JJ,indvx(i,j-1)];  AA = [AA,2*EtaP1/dx^2];
    II = [II,indvx(i,j)]; JJ = [JJ,indvx(i-1,j)];  AA = [AA,  Eta1 /dz^2];
    II = [II,indvx(i,j)]; JJ = [JJ,indvx(i+1,j)];  AA = [AA,  Eta2 /dz^2];
    RR = [RR,0];
    
%     A(indvx(i,j),indvx(i,j-1))      = 2*EtaP1/dx^2;                 % vx left of current node
%     A(indvx(i,j),indvx(i-1,j))      = Eta1/dz^2;                    % vx  above current node
    A(indvx(i,j),indvx(i,j))        = -2*(EtaP1+EtaP2)/dx^2 -...
                                        (Eta1+Eta2)/dz^2;           % vx current node
%     A(indvx(i,j),indvx(i+1,j))      = Eta2/dz^2;                    % vx below current node
    A(indvx(i,j),indvx(i,j+1))      = 2*EtaP2/dx^2;                 % vx right of current node
    
    A(indvx(i,j),indvz(i,j))        = -Eta2/dx/dz;                  % vz bottomleft
    A(indvx(i,j),indvz(i,j+1))      = Eta2/dx/dz;                   % vz bottomright      
    A(indvx(i,j),indvz(i-1,j))      = Eta1/dx/dz;                   % vz topleft
    A(indvx(i,j),indvz(i-1,j+1))    = -Eta1/dx/dz;                  % vz topright     
    A(indvx(i,j),indP(i,j))         = Pscale/dx;                    % P1; current node
    A(indvx(i,j),indP(i,j+1))       = -Pscale/dx;                   % P2; right of current node
    % RHS
    %RHS(indvx(i,j))                 = Rho_vx(i,j)*gx;               % x direction gravity
%     RHS(indvx(i,j))                 = 0;               % x direction gravity
    
    end
    
    %% z-Stokes eq. ETA*(d2Vy/dx^2+d2Vy/dy^2)-dP/dy=-RHO*gy
    % solve z equation
    % boundary conditions of vz
    if(j==1 || j==Nx || i==1 || i==nz1 || i==Nz)
            A(indvz(i,j),indvz(i,j))    = 1;        % A matrix coefficient
            RHS(indvz(i,j))             = 0;        % RHS
        %left boundary
        if (j==1 && i>1 && i<nz1)
            A(indvz(i,j),indvz(i,j+1))  = bcleft;   %solve for right of the leftmost bodes
        end
        %right boundary    
        if (j==Nx && i>1 && i<nz1)
            A(indvz(i,j),indvz(i,j-1))  = bcright;  % above the bottom boundary
        end
    % solve internal points    
    else
        
        Eta1    = Eta_out(i,j-1);
        Eta2    = Eta_out(i,j);
        EtaP1   = Eta_mid(i,j);
        EtaP2   = Eta_mid(i+1,j);
        
    % A matrix coefficients
    A(indvz(i,j),indvz(i,j-1))      = Eta1/dx^2;                    % vx1 left of current node
    A(indvz(i,j),indvz(i-1,j))      = 2*EtaP1/dz^2;                 % vx2 above current node
    A(indvz(i,j),indvz(i,j))        = -2*(EtaP1+EtaP2)/dz^2-...
                                        (Eta1+Eta2)/dx^2;           % vx3 current node
    A(indvz(i,j),indvz(i+1,j))      = 2*EtaP2/(dz^2);               % vx4 below current node
    A(indvz(i,j),indvz(i,j+1))      = Eta2/(dx^2);                  % vx5 right of current node
    A(indvz(i,j),indvx(i,j))        = -Eta2/dx/dz;                  % topright
    A(indvz(i,j),indvx(i+1,j))      = Eta2/dx/dz;                   % bottomright
    A(indvz(i,j),indvx(i,j-1))      = Eta1/dx/dz;                   % topleft
    A(indvz(i,j),indvx(i+1,j-1))    = -Eta1/dx/dz;                  % bottomleft
    
    A(indvz(i,j),indP(i,j))         = Pscale/dz;                    %P1; current node
    A(indvz(i,j),indP(i+1,j))       = -Pscale/dz;                   %P2; bottom of current node
    % RHS
    RHS(indvz(i,j))                 = -gz*Rho_vz(i,j);
    end
    
    %% P-Stokes eq. dVx/dx+dVy/dy=0
    % boundary points
    % P equation External points
    if(i==1 || j==1 || i==Nz || j==Nx ||...
            (i==2 && j==2))
        % Boundary Condition
        % 1*P=0
            A(indP(i,j),indP(i,j))  = 1;        % Left part
            RHS(indP(i,j))          = 0;        % Right part
         % Real BC
        if(i==2 && j==2)
            A(indP(i,j),indP(i,j))  = 1*Pscale; % Left part
            RHS(indP(i,j))          = 0;     % Right part   
        end
    % now solve internal points    
    else
    % A matrix coefficients    
    A(indP(i,j),indvx(i,j-1))       = -1/dx;    % left of current node
    A(indP(i,j),indvx(i,j))         = 1/dx;     % current node
    A(indP(i,j),indvz(i-1,j))       = -1/dz;    % above current node
    A(indP(i,j),indvz(i,j))         = 1/dz;     % below current node
    
    A(indP(i,j),indP(i,j))          = 0.001*dx*dz/Eta_mid(i,j);
    
    % RHS
    RHS(indP(i,j)) = 0;
    end       
end 
end

% Assemble coefficient matrix and right-hand side vector
A       = sparse(II,JJ,AA,N_all,N_all);
RHS     = sparse(II,ones(size(II)),RR,N_all,1);


%% Scale system of equations (diagonal preconditioning)
X           =  sqrt(abs(diag(A)));
X           =  diag(sparse(1./X));

A           =  X*A*X;
RHS         =  X*RHS;

%% Solve stokes matrix and convert output vector to matrices
c = X*(A\RHS); %get solution vector
%extrapolate into individual matrices
for j = 1:1:Nx
for i = 1:1:Nz
    
    P_out(i,j)  = c(indP(i,j))*Pscale;  %output pressure
    vx_out(i,j) = c(indvx(i,j));        %output vx
    vz_out(i,j) = c(indvz(i,j));        %output vz
    
end
end

    % averaging velocities on centre nodes
    vx_mid = zeros(Nz,Nx);
    vz_mid = zeros(Nz,Nx);    
    
    for j = 2:1:nx %solve for ordinary nodes only
    for i = 2:1:nz
        vx_mid(i,j) = 1/((1/vx_out(i,j)+1/vx_out(i,j-1))/2);% vx; (current+left)/2
        vz_mid(i,j) = 1/((1/vz_out(i,j)+1/vz_out(i-1,j))/2);% vz; (current+above)/2
    end
    end
    %applying free-slip boundary conditions
    %Top
    vx_mid(1,2:nx1)    = -bctop*vx_mid(2,2:nx1);
    vz_mid(1,:)         = -vz_mid(2,:);
    %bottom
    vx_mid(Nz,2:nx1)   = -bcbottom*vx_mid(nz1,2:nx1);
    vz_mid(Nz,:)        = -vz_mid(nz1,:);
    %left
    vx_mid(:,1)         = -vx_mid(:,2);
    vz_mid(2:nz,1)      = -bcleft*vz_mid(2:nz,2);
    %right
    vx_mid(:,Nx)       =-vx_mid(:,nx1);
    vz_mid(2:nz,Nx)  =-bcright*vz_mid(2:nz,nx1); % Free slip

