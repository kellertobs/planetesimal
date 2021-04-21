function [P_out,vx_out,vz_out,vx_mid,vz_mid] = stokes_continuity_noair(nx,nz,nx2,nz2,nx1,nz1,dx,dz,...
          Eta_out,Eta_mid,Rho_vx,Rho_vz,gx,gz,zp2d,bctop,bcbottom,bcleft,bcright,cstab)

NP      = nx2*nz2; % total number of P nodes to solve + ghost nodes
NU      = nx2*nz2; % total number of vx nodes to solve + ghost nodes
NW      = nx2*nz2; % total number of vz nodes to solve + ghost nodes
N_all   = NP+NU+NW;

%indexing of unknowns
indvx   = reshape(1:2:(NU+NW),nz2,nx2);
indvz   = reshape(2:2:(NW+NU),nz2,nx2);
indP    = reshape(1:NP,nz2,nx2) + NU + NW;

% setup A matrix and RHS vector
II  = [];
JJ  = [];
AA  = [];
IR  = [];
RR  = [];

Pscale = geomean(Eta_mid(:))/(dx*dz); % pressure scaling coefficient
RhoRef = mean(Rho_vz(:));


%% assemble coefficients of x momentum equation for vx
% top boundary (i==1 && j>1 && j<nx1)
ii = indvx(1,2:nx); jj1 = ii; jj2 = indvx(2,2:nx);
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];   AA = [AA, Asum+1];
II = [II, ii(:)']; JJ = [JJ, jj2(:)'];   AA = [AA, Asum+bctop];
% RHS
IR = [IR, ii(:)']; RR = [RR, Asum];

% bottom boundary (i==nz2 && j>1 && j<nx1)
ii = indvx(nz2,2:nx); jj1 = ii; jj2 = indvx(nz1,2:nx);
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];   AA = [AA, Asum(:)'+1];
II = [II, ii(:)']; JJ = [JJ, jj2(:)'];   AA = [AA, Asum(:)'+bcbottom];
% RHS
IR = [IR, ii(:)']; RR = [RR, Asum(:)'];

%left & right (+ghost right)
ii = indvx(1:nz2,1); jj = ii;
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'+1];
IR = [IR, ii(:)']; RR = [RR, Asum(:)'];
ii = indvx(1:nz2,nx1:nx2); jj = ii;
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
Rsum = zeros(size(ii)) - gx*(Rho_vx(2:nz1,2:nx)-RhoRef);
IR = [IR, ii(:)']; RR = [RR, Rsum(:)'];


%% assemble coefficients of z momentum equation for vz
% left boundary (j==1 && i>1 && i<nz1)
ii = indvz(2:nz,1); jj1 = ii; jj2 = indvz(2:nz,2);
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];   AA = [AA, Asum(:)'+1];
II = [II, ii(:)']; JJ = [JJ, jj2(:)'];   AA = [AA, Asum(:)'+bcleft];
% RHS
IR = [IR, ii(:)']; RR = [RR, Asum(:)'];

% right boundary (j==nx2 && i>1 && i<nz1)
ii = indvz(2:nz,nx2); jj1 = ii; jj2 = indvz(2:nz,nx1);
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];   AA = [AA, Asum(:)'+1];
II = [II, ii(:)']; JJ = [JJ, jj2(:)'];   AA = [AA, Asum(:)'+bcright];
% RHS
IR = [IR, ii(:)']; RR = [RR, Asum(:)'];

%top & bottom (+ghost bottom)
ii = [indvz(1,1:nx2), indvz(nz1,1:nx2), indvz(nz2,1:nx2)]; jj = ii;
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
Rsum = zeros(size(ii)) - gz*(Rho_vz(2:nz,2:nx1)-RhoRef);
IR = [IR, ii(:)']; RR = [RR, Rsum(:)'];


%% assemble coefficients of continuity equation for P
% boundary points
ii = [indP(1,:), indP(nz2,:)]; %top & bottom
jj = ii;
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'+1];
IR = [IR, ii(:)']; RR = [RR, Asum(:)'];

ii = [indP(2:nz1,1), indP(2:nz1, nx2)]; % left & right
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
Asum = cstab*dx*dz./Eta_mid(2:nz1,2:nx1);
II = [II, ii(:)']; JJ = [JJ, ii(:)'];     AA = [AA, Asum(:)'];

% RHS
Rsum = zeros(size(ii));
IR = [IR, ii(:)'];
RR = [RR, Rsum(:)'];

% real boundary condition P(2,2) = Rho*gz*dz/2
% RHOb = 3300;
% II = [II,indP(2,2)]; JJ = [JJ,indP(2,2)];   AA = [AA, 1*Pscale];
% IR = [IR,indP(2,2)]; RR = [RR, RHOb*dz/2*gz];


%% Assemble global coefficient matrix and right-hand side vector
A = sparse(II,JJ,AA,N_all,N_all);
R = sparse(IR,ones(size(IR)),RR,N_all,1);


%% Scale system of equations (diagonal preconditioning)
S  =  sqrt(abs(diag(A)));
S  =  diag(sparse(1./S));

A  =  S*A*S;
R  =  S*R;


%% Solve linear system of equations for vx, vz, P
X = S*(A\R); % get solution vector

% map solution vector to arrays
P_out  = reshape(X(indP(:)),nz2,nx2).*Pscale + RhoRef.*gz.*zp2d;
P_out  = P_out - mean(P_out(2,:)) + RhoRef*dz/2*gz;
vx_out = reshape(X(indvx(:)),nz2,nx2);
vz_out = reshape(X(indvz(:)),nz2,nx2);

% project velocities to centre nodes
vx_mid = zeros(nz2,nx2);
vz_mid = zeros(nz2,nx2);

vx_mid(:,2:nx1) = vx_out(:,2:nx1)+vx_out(:,1:nx)./2;
vz_mid(2:nz1,:) = vz_out(2:nz1,:)+vz_out(1:nz,:)./2;

% apply boundary conditions
%Top
vx_mid(1,2:nx1)   = -bctop*vx_mid(2,2:nx1);
vz_mid(1, :   )   = -      vz_mid(2, :   );
%bottom
vx_mid(nz2,2:nx1) = -bcbottom*vx_mid(nz1,2:nx1);
vz_mid(nz1, :   ) = -         vz_mid(nz , :   );
%left
vx_mid( :  ,1)    = -       vx_mid( :  ,2);
vz_mid(2:nz,1)    = -bcleft*vz_mid(2:nz,2);
%right
vx_mid( :  ,nx1)  = -        vx_mid( :  ,nx );
vz_mid(2:nz,nx2)  = -bcright*vz_mid(2:nz,nx1);

% get velocity divergence
Divv = diff(vx_out(2:end-1,1:end-1),1,2)./dx + diff(vz_out(1:end-1,2:end-1),1,1)./dz;

