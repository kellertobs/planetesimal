function [P_out,vx_out,vz_out,vx_mid,vz_mid] = stokes_continuity_noair(nx,nz,Nx,Nz,...
    nx1,nz1,dx,dz,Pscale,...
    Eta_out,Eta_mid,Rho_vx,Rho_vz,gx,gz,bctop,bcbottom,bcleft,bcright,cc)

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

%% solve internal points of x stokes equation
% top boundary (i==1 && j>1 && j<nx1)
ii = indvx(1,2:nx); jj1 = ii; jj2 = indvx(2,2:nx);
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];   AA = [AA, Asum+1];
II = [II, ii(:)']; JJ = [JJ, jj2(:)'];   AA = [AA, Asum+bctop];
% RHS
IR = [IR, ii(:)']; RR = [RR, Asum];

% bottom boundary (i==Nz && j>1 && j<nx1)
ii = indvx(Nz,2:nx); jj1 = ii; jj2 = indvx(nz1,2:nx);
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];   AA = [AA, Asum(:)'+1];
II = [II, ii(:)']; JJ = [JJ, jj2(:)'];   AA = [AA, Asum(:)'+bcbottom];
% RHS
IR = [IR, ii(:)']; RR = [RR, Asum(:)'];

%left & right (+ghost right)
ii = indvx(1:Nz,1); jj = ii;
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'+1];
IR = [IR, ii(:)']; RR = [RR, Asum(:)'];
ii = indvx(1:Nz,nx1:Nx); jj = ii;
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'+1];
IR = [IR, ii(:)']; RR = [RR, Asum(:)'];

% internal points
ii = indvx(2:nz1,2:nx);
Eta1 = Eta_out(1:nz,2:nx); Eta2 = Eta_out(2:nz1,2:nx); EtaP1 = Eta_mid(2:nz1,2:nx); EtaP2 = Eta_mid(2:nz1,3:nx1);
% vx nodes
% %           Left         ||          right          ||         top        ||          bottom
jj1 = indvx(2:nz1,1:nx-1); jj2 = indvx(2:nz1,3:nx+1); jj3 = indvx(1:nz,2:nx); jj4 = indvx(3:nz1+1,2:nx);
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];   AA = [AA, 2.*EtaP1(:)'./dx^2]; % vx left of current node
II = [II, ii(:)']; JJ = [JJ, jj3(:)'];   AA = [AA,    Eta1(:)'./dz^2]; % vx  above current node
Asum = -2.*(EtaP1+EtaP2)./dx^2 - (Eta1+Eta2)./dz^2;
II = [II, ii(:)']; JJ = [JJ, ii(:)'];    AA = [AA, Asum(:)']; % vx current node
II = [II, ii(:)']; JJ = [JJ, jj4(:)'];   AA = [AA,    Eta2(:)'./dz^2]; % vx below current node
II = [II, ii(:)']; JJ = [JJ, jj2(:)'];   AA = [AA, 2.*EtaP2(:)'./dx^2]; % vx right of current node

% vz nodes
% %       bottom Left    ||        bottom right     ||       top left       ||     top right
jj1 = indvz(2:nz1,2:nx); jj2 = indvz(2:nz1,3:nx+1); jj3 = indvz(1:nz,2:nx); jj4 = indvz(1:nz,3:nx1);
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];   AA = [AA,   -Eta2(:)'/dx/dz]; % vz bottomleft
II = [II, ii(:)']; JJ = [JJ, jj2(:)'];   AA = [AA,    Eta2(:)'/dx/dz]; % vz bottomright
II = [II, ii(:)']; JJ = [JJ, jj3(:)'];   AA = [AA,    Eta1(:)'/dx/dz]; % vz topleft
II = [II, ii(:)']; JJ = [JJ, jj4(:)'];   AA = [AA,   -Eta1(:)'/dx/dz]; % vz topright

% P nodes
Asum = zeros(size(ii));
jj1 = indP(2:nz1,2:nx); jj2 = indP(2:nz1,3:nx1);
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];   AA = [AA,    Asum(:)'+Pscale/dx]; % P1; current node
II = [II, ii(:)']; JJ = [JJ, jj2(:)'];   AA = [AA,    Asum(:)'-Pscale/dx]; % P2; right of current node

% RHS vector
Rsum = zeros(size(ii)) - gx*Rho_vx(2:nz1,2:nx);
IR = [IR, ii(:)']; RR = [RR, Rsum(:)'];

%% solve internal points of z stokes equation
% left boundary (j==1 && i>1 && i<nz1)
ii = indvz(2:nz,1); jj1 = ii; jj2 = indvz(2:nz,2);
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];   AA = [AA, Asum(:)'+1];
II = [II, ii(:)']; JJ = [JJ, jj2(:)'];   AA = [AA, Asum(:)'+bcleft];
% RHS
IR = [IR, ii(:)']; RR = [RR, Asum(:)'];

% right boundary (j==Nx && i>1 && i<nz1)
ii = indvz(2:nz,Nx); jj1 = ii; jj2 = indvz(2:nz,nx1);
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];   AA = [AA, Asum(:)'+1];
II = [II, ii(:)']; JJ = [JJ, jj2(:)'];   AA = [AA, Asum(:)'+bcright];
% RHS
IR = [IR, ii(:)']; RR = [RR, Asum(:)'];

%top & bottom (+ghost bottom)
ii = [indvz(1,1:Nx), indvz(nz1,1:Nx), indvz(Nz,1:Nx)]; jj = ii;
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'+1];
IR = [IR, ii(:)']; RR = [RR, Asum(:)'];

ii = indvz(2:nz,2:nx1);
% %       Left        ||          right          ||           top          ||          bottom
jj1 = indvz(2:nz,1:nx); jj2 = indvz(2:nz,3:nx1+1); jj3 = indvz(1:nz-1,2:nx1); jj4 = indvz(3:nz1,2:nx1);

Eta1 = Eta_out(2:nz,1:nx); Eta2 = Eta_out(2:nz,2:nx1); EtaP1 = Eta_mid(2:nz,2:nx1); EtaP2 = Eta_mid(3:nz1,2:nx1);

% vz
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];  AA = [AA,    Eta1(:)'/dx^2]; % vx1 left of current node
II = [II, ii(:)']; JJ = [JJ, jj3(:)'];  AA = [AA,  2*EtaP1(:)'/dz^2]; % vx2 above current node
Asum = -2*(EtaP1+EtaP2)/dz^2-(Eta1+Eta2)/dx^2;
II = [II, ii(:)']; JJ = [JJ, ii(:)'];   AA = [AA,    Asum(:)']; % vx3 current node
II = [II, ii(:)']; JJ = [JJ, jj4(:)'];  AA = [AA,  2*EtaP2(:)'/dz^2]; % vx4 below current node
II = [II, ii(:)']; JJ = [JJ, jj2(:)'];  AA = [AA,    Eta2(:)'/dx^2]; % vx5 right of current node

% %       bottom Left    ||        bottom right     ||       top left       ||     top right
jj1 = indvx(2:nz,2:nx1); jj2 = indvx(3:nz1,2:nx1); jj3 = indvx(2:nz,1:nx); jj4 = indvx(3:nz1,1:nx);
% vx
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];  AA = [AA,   -Eta2(:)'/dx/dz]; % topright
II = [II, ii(:)']; JJ = [JJ, jj2(:)'];  AA = [AA,    Eta2(:)'/dx/dz]; % bottomright
II = [II, ii(:)']; JJ = [JJ, jj3(:)'];  AA = [AA,    Eta1(:)'/dx/dz]; % topleft
II = [II, ii(:)']; JJ = [JJ, jj4(:)'];  AA = [AA,   -Eta1(:)'/dx/dz]; % bottomleft

% P nodes
Asum = zeros(size(ii));
jj1 = indP(2:nz,2:nx1); jj2 = indP(3:nz1,2:nx1);
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];   AA = [AA,   Asum(:)'+Pscale/dx]; % P1; current node
II = [II, ii(:)']; JJ = [JJ, jj2(:)'];   AA = [AA,   Asum(:)'-Pscale/dx]; % P2; right of current node

% RHS
Rsum = zeros(size(ii)) - gz*Rho_vz(2:nz,2:nx1);
IR = [IR, ii(:)']; RR = [RR, Rsum(:)'];

%% solve internal points of P continuity equation
% boundary points
ii = [indP(1,:), indP(Nz,:)]; %top & bottom
jj = ii;
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'+1];
IR = [IR, ii(:)']; RR = [RR, Asum(:)'];

ii = [indP(2:nz1,1), indP(2:nz1, Nx)]; % left & right
jj = ii;
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'+1];
IR = [IR, ii(:)']; RR = [RR, Asum(:)'];

%internal points
ii = indP(2:nz1,2:nx1);
% %       vx 1         ||          vx 2          ||           vz 1         ||          vz 2
jj1 = indvx(2:nz1,1:nx); jj2 = indvx(2:nz1,2:nx1); jj3 = indvz(1:nz,2:nx1); jj4 = indvz(2:nz1,2:nx1);
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];    AA = [AA, Asum(:)'-1/dx]; % left of current node
II = [II, ii(:)']; JJ = [JJ, jj2(:)'];    AA = [AA, Asum(:)'+1/dx]; % current node
II = [II, ii(:)']; JJ = [JJ, jj3(:)'];    AA = [AA, Asum(:)'-1/dz]; % above current node
II = [II, ii(:)']; JJ = [JJ, jj4(:)'];    AA = [AA, Asum(:)'+1/dz]; % below current node
Asum = cc*dx*dz./Eta_mid(2:nz1,2:nx1);
II = [II, ii(:)']; JJ = [JJ, ii(:)'];     AA = [AA, Asum(:)'];

% RHS
Rsum = zeros(size(ii));
IR = [IR, ii(:)'];
RR = [RR, Rsum(:)'];

% % % real boundary condition P(2,2) = Rho*gz*dz/2
RHOb = 3300;
II = [II,indP(2,2)]; JJ = [JJ,indP(2,2)];   AA = [AA, 1*Pscale];
IR = [IR,indP(2,2)]; RR = [RR, RHOb*dz/2*gz];

%% Assemble coefficient matrix and right-hand side vector
A       = sparse(II,JJ,AA,N_all,N_all);
RHS     = sparse(IR,ones(size(IR)),RR,N_all,1);


%% Scale system of equations (diagonal preconditioning)
X           =  sqrt(abs(diag(A)));
X           =  diag(sparse(1./X));

A           =  X*A*X;
RHS         =  X*RHS;
%
% tic
% %% Solve stokes matrix and convert output vector to matrices
c = X*(A\RHS); %get solution vector
% c = A\RHS;
% ccp = toc;
% 
% figure(6)
% subplot(1,2,1)
% semilogx(cc,ccp,'*')
% ylabel('time elapsed')
% xlabel('log(cc)')
% hold on

%extrapolate into individual matrices
P_out = reshape(c(indP(:)),Nz,Nx).*Pscale;%output pressure
vx_out = reshape(c(indvx(:)),Nz,Nx);
vz_out = reshape(c(indvz(:)),Nz,Nx);


% figure(5)
% imagesc(1:1:nx,1:1:nz,P_out(2:end-1,2:end-1)./1e6)
% colorbar
% title('Pressure MPa')

% averaging velocities on centre nodes
vx_mid = zeros(Nz,Nx);
vz_mid = zeros(Nz,Nx);

vx_mid(:,2:nx1) = vx_out(:,2:nx1)+vx_out(:,1:nx)./2;
vz_mid(2:nz1,:) = vz_out(2:nz1,:)+vz_out(1:nz,:)./2;


%applying free-slip boundary conditions
%Top
vx_mid(1,2:nx1)    = -bctop*vx_mid(2,2:nx1);
vz_mid(1,:)         = -vz_mid(2,:);
%bottom
vx_mid(Nz,2:nx1)   = -bcbottom*vx_mid(nz1,2:nx1);
vz_mid(nz1,:)        = -vz_mid(nz,:);
%left
vx_mid(:,1)         = -vx_mid(:,2);
vz_mid(2:nz,1)      = -bcleft*vz_mid(2:nz,2);
%right
vx_mid(:,nx1)       =-vx_mid(:,nx);
vz_mid(2:nz,Nx)  =-bcright*vz_mid(2:nz,nx1); % Free slip

% legacy loop to fill A matrix
% for j = 1:1:Nx
% for i = 1:1:Nz
%     %% x-Stokes eq. ETA*(d2Vx/dx^2+d2Vx/dy^2)-dP/dx=0
%     % solve x equation
%
%     % boundary conditions of vx
%     if(i==1 || i==Nz || j==1 || j==nx1 || j==Nx)
%         II = [II,indvx(i,j)]; JJ = [JJ,indvx(i,j)];   AA = [AA, 1];
%         IR = [IR,indvx(i,j)]; RR = [RR, 0];
% %             A(indvx(i,j),indvx(i,j))    = 1; % A matrix coefficient
% %             RHS(indvx(i,j))             = 0; % RHS
%         %Top boundary
%         if (i==1 && j>1 && j<nx1)
%             II = [II,indvx(i,j)]; JJ = [JJ,indvx(i+1,j)];   AA = [AA, bctop];
% %             A(indvx(i,j),indvx(i+1,j))  = bctop; %only solve for the bottom of the top boundary
%         end
%         %Bottom boundary
%         if (i==Nz && j>1 && j<nx1)
%             II = [II,indvx(i,j)]; JJ = [JJ,indvx(i-1,j)];   AA = [AA, bcbottom];
% %             A(indvx(i,j),indvx(i-1,j))  = bcbottom; % above the bottom boundary
%         end
%     % now solve internal points on the real grid
%     else
%     % A matrix coefficients
%     Eta1    = Eta_out(i-1,j);
%     Eta2    = Eta_out(i,j);
%     EtaP1   = Eta_mid(i,j);
%     EtaP2   = Eta_mid(i,j+1);
%
%     %vx
%     II = [II,indvx(i,j)]; JJ = [JJ,indvx(i,j-1)];   AA = [AA, 2*EtaP1/dx^2]; % vx left of current node
%     II = [II,indvx(i,j)]; JJ = [JJ,indvx(i-1,j)];   AA = [AA,    Eta1/dz^2]; % vx  above current node
%     II = [II,indvx(i,j)]; JJ = [JJ,indvx(i,j)];     AA = [AA,-2*(EtaP1+EtaP2)/dx^2 - (Eta1+Eta2)/dz^2]; % vx current node
%     II = [II,indvx(i,j)]; JJ = [JJ,indvx(i+1,j)];   AA = [AA,   Eta2 /dz^2]; % vx below current node
%     II = [II,indvx(i,j)]; JJ = [JJ,indvx(i,j+1)];   AA = [AA, 2*EtaP2/dx^2]; % vx right of current node
%     %vz
%     II = [II,indvx(i,j)]; JJ = [JJ,indvz(i,j)];     AA = [AA,   -Eta2/dx/dz]; % vz bottomleft
%     II = [II,indvx(i,j)]; JJ = [JJ,indvz(i,j+1)];   AA = [AA,    Eta2/dx/dz]; % vz bottomright
%     II = [II,indvx(i,j)]; JJ = [JJ,indvz(i-1,j)];   AA = [AA,    Eta1/dx/dz]; % vz topleft
%     II = [II,indvx(i,j)]; JJ = [JJ,indvz(i-1,j+1)]; AA = [AA,   -Eta1/dx/dz]; % vz topright
%
%     II = [II,indvx(i,j)]; JJ = [JJ,indP(i,j)];      AA = [AA,    Pscale/dx]; % P1; current node
%     II = [II,indvx(i,j)]; JJ = [JJ,indP(i,j+1)];    AA = [AA,   -Pscale/dx]; % P2; right of current node
%     %RHS
%     IR = [IR,indvx(i,j)]; RR = [RR, 0];
%
%
%
% %     A(indvx(i,j),indP(i,j))         = Pscale/dx;
% %     A(indvx(i,j),indP(i,j+1))       = -Pscale/dx;
%     % RHS
%     %RHS(indvx(i,j))                 = Rho_vx(i,j)*gx;               % x direction gravity
% %     RHS(indvx(i,j))                 = 0;               % x direction gravity
%
%     end
%
%     %% z-Stokes eq. ETA*(d2Vy/dx^2+d2Vy/dy^2)-dP/dy=-RHO*gy
%     % solve z equation
%     % boundary conditions of vz
%     if(j==1 || j==Nx || i==1 || i==nz1 || i==Nz)
%         II = [II,indvz(i,j)]; JJ = [JJ,indvz(i,j)];   AA = [AA, 1];
%         IR = [IR,indvz(i,j)]; RR = [RR, 0];
% %             A(indvz(i,j),indvz(i,j))    = 1;        % A matrix coefficient
% %             RHS(indvz(i,j))             = 0;        % RHS
%         %left boundary
%         if (j==1 && i>1 && i<nz1)
%             II = [II,indvz(i,j)]; JJ = [JJ,indvz(i,j+1)];   AA = [AA, bcleft];
% %             A(indvz(i,j),indvz(i,j+1))  = bcleft;   %solve for right of the leftmost bodes
%         end
%         %right boundary
%         if (j==Nx && i>1 && i<nz1)
%             II = [II,indvz(i,j)]; JJ = [JJ,indvz(i,j-1)];   AA = [AA, bcright];
% %             A(indvz(i,j),indvz(i,j-1))  = bcright;  % above the bottom boundary
%         end
%     % solve internal points
%     else
%
%         Eta1    = Eta_out(i,j-1);
%         Eta2    = Eta_out(i,j);
%         EtaP1   = Eta_mid(i,j);
%         EtaP2   = Eta_mid(i+1,j);
%
%     %vz
%     II = [II,indvz(i,j)]; JJ = [JJ,indvz(i,j-1)];  AA = [AA,    Eta1/dx^2]; % vx1 left of current node
%     II = [II,indvz(i,j)]; JJ = [JJ,indvz(i-1,j)];  AA = [AA,  2*EtaP1/dz^2]; % vx2 above current node
%     II = [II,indvz(i,j)]; JJ = [JJ,indvz(i,j)];    AA = [AA,-2*(EtaP1+EtaP2)/dz^2-(Eta1+Eta2)/dx^2]; % vx3 current node
%     II = [II,indvz(i,j)]; JJ = [JJ,indvz(i+1,j)];  AA = [AA,  2*EtaP2/dz^2]; % vx4 below current node
%     II = [II,indvz(i,j)]; JJ = [JJ,indvz(i,j+1)];  AA = [AA,    Eta2/dx^2]; % vx5 right of current node
%     %vx
%     II = [II,indvz(i,j)]; JJ = [JJ,indvx(i,j)];     AA = [AA,   -Eta2/dx/dz]; % topright
%     II = [II,indvz(i,j)]; JJ = [JJ,indvx(i+1,j)];   AA = [AA,    Eta2/dx/dz]; % bottomright
%     II = [II,indvz(i,j)]; JJ = [JJ,indvx(i,j-1)];   AA = [AA,    Eta1/dx/dz]; % topleft
%     II = [II,indvz(i,j)]; JJ = [JJ,indvx(i+1,j-1)]; AA = [AA,   -Eta1/dx/dz]; % bottomleft
%     II = [II,indvz(i,j)]; JJ = [JJ,indP(i,j)];      AA = [AA,    Pscale/dx];  % P1; current node
%     II = [II,indvz(i,j)]; JJ = [JJ,indP(i+1,j)];    AA = [AA,    -Pscale/dx]; % P2; bottom of current node
%
%     IR = [IR,indvz(i,j)]; RR = [RR, -gz*Rho_vz(i,j)];
%
%     end
%
%     %% P-Stokes eq. dVx/dx+dVy/dy=0
%     % boundary points
%     % P equation External points
%     if(i==1 || j==1 || i==Nz || j==Nx)
%         % Boundary Condition
%         % 1*P=0
%         II = [II,indP(i,j)]; JJ = [JJ,indP(i,j)];   AA = [AA, 1];
%         IR = [IR,indP(i,j)]; RR = [RR, 0];
%          % Real BC
%     elseif (i==2 && j==2)
%             II = [II,indP(i,j)]; JJ = [JJ,indP(i,j)];   AA = [AA, 1*Pscale];
%             IR = [IR,indP(i,j)]; RR = [RR, 0];
%
%     % now solve internal points
%     else
%
%     II = [II,indP(i,j)]; JJ = [JJ,indvx(i,j-1)];  AA = [AA, -1/dx]; % left of current node
%     II = [II,indP(i,j)]; JJ = [JJ,indvx(i,j)];    AA = [AA,  1/dx]; % current node
%     II = [II,indP(i,j)]; JJ = [JJ,indvz(i-1,j)];  AA = [AA, -1/dz]; % above current node
%     II = [II,indP(i,j)]; JJ = [JJ,indvz(i,j)];    AA = [AA,  1/dz]; % below current node
%     II = [II,indP(i,j)]; JJ = [JJ,indP(i,j)];     AA = [AA, 0.001*dx*dz/Eta_mid(i,j)];
%     % RHS
%     IR = [IR,indP(i,j)]; RR = [RR, 0];
%
%     end
% end
% end
