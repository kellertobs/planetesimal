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
% A = sparse(N_all,N_all);
% RHS = zeros(N_all,1);
II  = [];
JJ  = [];
AA  = [];
IR  = [];
RR  = [];

for j = 1:1:Nx
for i = 1:1:Nz
    %% x-Stokes eq. ETA*(d2Vx/dx^2+d2Vx/dy^2)-dP/dx=0
    % solve x equation
    
    % boundary conditions of vx
    if(i==1 || i==Nz || j==1 || j==nx1 || j==Nx)
        II = [II,indvx(i,j)]; JJ = [JJ,indvx(i,j)];   AA = [AA, 1];
        IR = [IR,indvx(i,j)]; RR = [RR, 0];
%             A(indvx(i,j),indvx(i,j))    = 1; % A matrix coefficient
%             RHS(indvx(i,j))             = 0; % RHS
        %Top boundary
        if (i==1 && j>1 && j<nx1)
            II = [II,indvx(i,j)]; JJ = [JJ,indvx(i+1,j)];   AA = [AA, bctop];
%             A(indvx(i,j),indvx(i+1,j))  = bctop; %only solve for the bottom of the top boundary
        end
        %Bottom boundary    
        if (i==Nz && j>1 && j<nx1)
            II = [II,indvx(i,j)]; JJ = [JJ,indvx(i-1,j)];   AA = [AA, bcbottom];
%             A(indvx(i,j),indvx(i-1,j))  = bcbottom; % above the bottom boundary
        end
    % now solve internal points on the real grid
    else
    % A matrix coefficients
    Eta1    = Eta_out(i-1,j);
    Eta2    = Eta_out(i,j);
    EtaP1   = Eta_mid(i,j);
    EtaP2   = Eta_mid(i,j+1);
    
    %vx
    II = [II,indvx(i,j)]; JJ = [JJ,indvx(i,j-1)];   AA = [AA, 2*EtaP1/dx^2]; % vx left of current node
    II = [II,indvx(i,j)]; JJ = [JJ,indvx(i-1,j)];   AA = [AA,    Eta1/dz^2]; % vx  above current node
    II = [II,indvx(i,j)]; JJ = [JJ,indvx(i,j)];     AA = [AA,-2*(EtaP1+EtaP2)/dx^2 - (Eta1+Eta2)/dz^2]; % vx current node
    II = [II,indvx(i,j)]; JJ = [JJ,indvx(i+1,j)];   AA = [AA,   Eta2 /dz^2]; % vx below current node
    II = [II,indvx(i,j)]; JJ = [JJ,indvx(i,j+1)];   AA = [AA, 2*EtaP2/dx^2]; % vx right of current node
    %vz
    II = [II,indvx(i,j)]; JJ = [JJ,indvz(i,j)];     AA = [AA,   -Eta2/dx/dz]; % vz bottomleft
    II = [II,indvx(i,j)]; JJ = [JJ,indvz(i,j+1)];   AA = [AA,    Eta2/dx/dz]; % vz bottomright
    II = [II,indvx(i,j)]; JJ = [JJ,indvz(i-1,j)];   AA = [AA,    Eta1/dx/dz]; % vz topleft
    II = [II,indvx(i,j)]; JJ = [JJ,indvz(i-1,j+1)]; AA = [AA,   -Eta1/dx/dz]; % vz topright
    
    II = [II,indvx(i,j)]; JJ = [JJ,indP(i,j)];      AA = [AA,    Pscale/dx]; % P1; current node
    II = [II,indvx(i,j)]; JJ = [JJ,indP(i,j+1)];    AA = [AA,   -Pscale/dx]; % P2; right of current node
    %RHS
    IR = [IR,indvx(i,j)]; RR = [RR, 0];
    
                                    
                    
%     A(indvx(i,j),indP(i,j))         = Pscale/dx;                    
%     A(indvx(i,j),indP(i,j+1))       = -Pscale/dx;                  
    % RHS
    %RHS(indvx(i,j))                 = Rho_vx(i,j)*gx;               % x direction gravity
%     RHS(indvx(i,j))                 = 0;               % x direction gravity
    
    end
    
    %% z-Stokes eq. ETA*(d2Vy/dx^2+d2Vy/dy^2)-dP/dy=-RHO*gy
    % solve z equation
    % boundary conditions of vz
    if(j==1 || j==Nx || i==1 || i==nz1 || i==Nz)
        II = [II,indvz(i,j)]; JJ = [JJ,indvz(i,j)];   AA = [AA, 1];
        IR = [IR,indvz(i,j)]; RR = [RR, 0];
%             A(indvz(i,j),indvz(i,j))    = 1;        % A matrix coefficient
%             RHS(indvz(i,j))             = 0;        % RHS
        %left boundary
        if (j==1 && i>1 && i<nz1)
            II = [II,indvz(i,j)]; JJ = [JJ,indvz(i,j+1)];   AA = [AA, bcleft];
%             A(indvz(i,j),indvz(i,j+1))  = bcleft;   %solve for right of the leftmost bodes
        end
        %right boundary    
        if (j==Nx && i>1 && i<nz1)
            II = [II,indvz(i,j)]; JJ = [JJ,indvz(i,j-1)];   AA = [AA, bcright];
%             A(indvz(i,j),indvz(i,j-1))  = bcright;  % above the bottom boundary
        end
    % solve internal points    
    else
        
        Eta1    = Eta_out(i,j-1);
        Eta2    = Eta_out(i,j);
        EtaP1   = Eta_mid(i,j);
        EtaP2   = Eta_mid(i+1,j);

    %vz
    II = [II,indvz(i,j)]; JJ = [JJ,indvz(i,j-1)];  AA = [AA,    Eta1/dx^2]; % vx1 left of current node
    II = [II,indvz(i,j)]; JJ = [JJ,indvz(i-1,j)];  AA = [AA,  2*EtaP1/dz^2]; % vx2 above current node
    II = [II,indvz(i,j)]; JJ = [JJ,indvz(i,j)];    AA = [AA,-2*(EtaP1+EtaP2)/dz^2-(Eta1+Eta2)/dx^2]; % vx3 current node
    II = [II,indvz(i,j)]; JJ = [JJ,indvz(i+1,j)];  AA = [AA,  2*EtaP2/dz^2]; % vx4 below current node
    II = [II,indvz(i,j)]; JJ = [JJ,indvz(i,j+1)];  AA = [AA,    Eta2/dx^2]; % vx5 right of current node
    %vx
    II = [II,indvz(i,j)]; JJ = [JJ,indvx(i,j)];     AA = [AA,   -Eta2/dx/dz]; % topright
    II = [II,indvz(i,j)]; JJ = [JJ,indvx(i+1,j)];   AA = [AA,    Eta2/dx/dz]; % bottomright
    II = [II,indvz(i,j)]; JJ = [JJ,indvx(i,j-1)];   AA = [AA,    Eta1/dx/dz]; % topleft
    II = [II,indvz(i,j)]; JJ = [JJ,indvx(i+1,j-1)]; AA = [AA,   -Eta1/dx/dz]; % bottomleft
    II = [II,indvz(i,j)]; JJ = [JJ,indP(i,j)];      AA = [AA,    Pscale/dx];  % P1; current node
    II = [II,indvz(i,j)]; JJ = [JJ,indP(i+1,j)];    AA = [AA,    -Pscale/dx]; % P2; bottom of current node  
        
    IR = [IR,indvz(i,j)]; RR = [RR, -gz*Rho_vz(i,j)];

    end
    
    %% P-Stokes eq. dVx/dx+dVy/dy=0
    % boundary points
    % P equation External points
    if(i==1 || j==1 || i==Nz || j==Nx)
        % Boundary Condition
        % 1*P=0
        II = [II,indP(i,j)]; JJ = [JJ,indP(i,j)];   AA = [AA, 1];
        IR = [IR,indP(i,j)]; RR = [RR, 0];
         % Real BC
    elseif (i==2 && j==2)
            II = [II,indP(i,j)]; JJ = [JJ,indP(i,j)];   AA = [AA, 1*Pscale];
            IR = [IR,indP(i,j)]; RR = [RR, 0];

    % now solve internal points    
    else
    
    II = [II,indP(i,j)]; JJ = [JJ,indvx(i,j-1)];  AA = [AA, -1/dx]; % left of current node
    II = [II,indP(i,j)]; JJ = [JJ,indvx(i,j)];    AA = [AA,  1/dx]; % current node
    II = [II,indP(i,j)]; JJ = [JJ,indvz(i-1,j)];  AA = [AA, -1/dz]; % above current node
    II = [II,indP(i,j)]; JJ = [JJ,indvz(i,j)];    AA = [AA,  1/dz]; % below current node
    II = [II,indP(i,j)]; JJ = [JJ,indP(i,j)];     AA = [AA, 0.001*dx*dz/Eta_mid(i,j)];
    % RHS
    IR = [IR,indP(i,j)]; RR = [RR, 0];
    
    end       
end 
end

% Assemble coefficient matrix and right-hand side vector
A       = sparse(II,JJ,AA,N_all,N_all);
RHS     = sparse(IR,ones(size(IR)),RR,N_all,1);


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

