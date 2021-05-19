% planetesimal: solve fluid mechanics equations

% print solver header
fprintf(1,'  ---  solve fluid mechanics equations \n');


%% prepare auxiliary variables

% get characteristic scales
Pscale = geomean(MAT.Eta.s(:))/(NUM.dx*NUM.dz); % pressure scaling coefficient
RhoRef = mean(MAT.Rho.t(:));                    % mean density for lithostatic pressure

% get mapping arrays
indU = NUM.MapU;
indW = NUM.MapW;
indP = NUM.MapP;

P_cmp = 0.*SOL.P.l(2:end-1,2:end-1);
w_seg = 0.*SOL.W.s(:,2:end-1);
u_seg = 0.*SOL.U.s(2:end-1,:);

%% Momentum equation; Row 1

% initialise global lists for vectorised assembly
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R
% % ================================================== % %
% Diagonal block [1,1]; UW_s Stress divergence
% % ================================================== % %
% Z Stokes equation
% boundary conditions
% left boundary
ii = indW(:,1); jj1 = ii; jj2 = indW(:,2);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa+SOL.BCleft];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% bottom boundary
ii = indW(:,end); jj1 = ii; jj2 = indW(:,end-1);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+SOL.BCright];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% top boundary
ii = indW(1,:).'; jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa+1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% bottom boundary
ii = indW(end,:).'; jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa+1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% internal points
ii    = indW(2:end-1,2:end-1);
EtaC1 = MAT.EtaC.s(2:end-1,1:end-1); EtaC2 = MAT.EtaC.s(2:end-1,2:end  );
EtaP1 = MAT.Eta.s (2:end-2,2:end-1); EtaP2 = MAT.Eta.s (3:end-1,2:end-1);

% coefficients multiplying z-velocities W
%             top          ||         bottom          ||           left            ||          right
jj1 = indW(1:end-2,2:end-1); jj2 = indW(3:end,2:end-1); jj3 = indW(2:end-1,1:end-2); jj4 = indW(2:end-1,3:end);

aa = -(EtaP1+EtaP2)./NUM.dz^2 - (EtaC1+EtaC2)./NUM.dx^2;
II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; aa(:)             ];      % W on stencil centre
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; EtaP1(:)./NUM.dz^2];      % W one above
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; EtaP2(:)./NUM.dz^2];      % W one below
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; EtaC1(:)./NUM.dx^2];      % W one to the left
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; EtaC2(:)./NUM.dx^2];      % W one to the right

% coefficients multiplying x-velocities U
%         top left         ||        bottom left          ||       top right       ||       bottom right
jj1 = indU(2:end-2,1:end-1); jj2 = indU(3:end-1,1:end-1); jj3 = indU(2:end-2,2:end); jj4 = indU(3:end-1,2:end);

II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; EtaC1(:)/NUM.dx/NUM.dz];  % W one to the top and left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA;-EtaC1(:)/NUM.dx/NUM.dz];  % W one to the bottom and left
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA;-EtaC2(:)/NUM.dx/NUM.dz];  % W one to the top and right
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; EtaC2(:)/NUM.dx/NUM.dz];  % W one to the bottom and right

% assemble RHS vector
rr = zeros(size(ii)) - PHY.gz*( (MAT.Rho.t(2:end-2,2:end-1)+MAT.Rho.t(3:end-1,2:end-1))/2 - RhoRef);
%     + diff(P_cmp,1,1)/NUM.dz;
IR = [IR; ii(:)];  RR = [RR; rr(:)];

% X Stokes equation
% top boundary
ii = indU(1,:).'; jj1 = ii; jj2 = indU(2,:).';
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa+SOL.BCtop];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% bottom boundary
ii = indU(end,:).'; jj1 = ii; jj2 = indU(end-1,:).';
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+SOL.BCbot];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% left side boundary
ii = indU(:,1); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa+1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% right side boundary
ii = indU(:,end); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa+1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% internal points
ii    = indU(2:end-1,2:end-1);
EtaC1 = MAT.EtaC.s(1:end-1,2:end-1); EtaC2 = MAT.EtaC.s(2:end  ,2:end-1);
EtaP1 = MAT.Eta.s (2:end-1,2:end-2); EtaP2 = MAT.Eta.s (2:end-1,3:end-1);

% coefficients multiplying x-velocities U
%            left          ||          right          ||           top             ||          bottom
jj1 = indU(2:end-1,1:end-2); jj2 = indU(2:end-1,3:end); jj3 = indU(1:end-2,2:end-1); jj4 = indU(3:end,2:end-1);

aa = -(EtaP1+EtaP2)./NUM.dx^2 - (EtaC1+EtaC2)./NUM.dz^2;
II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; aa(:)             ];      % U on stencil centre
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; EtaP1(:)./NUM.dx^2];      % U one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; EtaP2(:)./NUM.dx^2];      % U one to the right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; EtaC1(:)./NUM.dz^2];      % U one above
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; EtaC2(:)./NUM.dz^2];      % U one below

% coefficients multiplying z-velocities W
%         top left         ||        top right          ||       bottom left       ||       bottom right
jj1 = indW(1:end-1,2:end-2); jj2 = indW(1:end-1,3:end-1); jj3 = indW(2:end,2:end-2); jj4 = indW(2:end,3:end-1);

II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; EtaC1(:)/NUM.dx/NUM.dz];  % W one to the top and left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA;-EtaC1(:)/NUM.dx/NUM.dz];  % W one to the top and right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA;-EtaC2(:)/NUM.dx/NUM.dz];  % W one to the bottom and left
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; EtaC2(:)/NUM.dx/NUM.dz];  % W one to the bottom and right

% assemble RHS vector
rr = zeros(size(ii)) - PHY.gx*( (MAT.Rho.t(2:end-1,2:end-2)+MAT.Rho.t(2:end-1,3:end-1))/2 - RhoRef);
%     + diff(P_cmp,1,1)/NUM.dz;
IR = [IR; ii(:)];  RR = [RR; rr(:)];

% assemble global coefficient matrix of stress divergence
KV = sparse(II,JJ,AA,NUM.NW+NUM.NU,NUM.NW+NUM.NU);
RV = sparse(IR,ones(size(IR)),RR,NUM.NW+NUM.NU,1);

%% Block [2,2] Assemble permeability coefficients (W_seg, U_seg)
% reset global lists for vectorised assembly
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R

% % ================================================== % %
% coefficients multiplying z-velocities W
% % ================================================== % %
% boundary points
ii  = [indW(1,:).'; indW(end  ,:).']; % top & bottom
jj1 = ii;
jj2 = [indW(2,:).'; indW(end-1,:).'];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

ii  = [indW(:,1).'; indW(:,end  ).']; % left & right
jj1 = ii;
jj2 = [indW(:,2).'; indW(:,end-1).'];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii = indW(2:end-1,2:end-1);
aa = MAT.EtaW.l(2:end-1,2:end-1);
II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; aa(:)];      % W on stencil centre

% assemble RHS vector
rr = zeros(size(ii)) - PHY.gz*( (MAT.Rho.l(2:end-2,2:end-1)+MAT.Rho.l(3:end-1,2:end-1))/2 - RhoRef);
%     + diff(P_cmp,1,1)/NUM.dz;
IR = [IR; ii(:)];  RR = [RR; rr(:)];

% % ================================================== % %
% coefficients multiplying x-velocities U
% % ================================================== % %
% boundary points
ii  = [indU(1,:).'; indU(end  ,:).']; % top & bottom
jj1 = ii;
jj2 = [indU(2,:).'; indU(end-1,:).'];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

ii  = [indU(:,1).'; indU(:,end  ).']; % left & right
jj1 = ii;
jj2 = [indU(:,2).'; indU(:,end-1).'];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];
% internal points
ii = indU(2:end-1,2:end-1);
aa = MAT.EtaU.l(2:end-1,2:end-1);
II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; aa(:)];      % W on stencil centre

% assemble RHS vector
rr = zeros(size(ii)) - PHY.gx*( (MAT.Rho.l(2:end-1,2:end-2)+MAT.Rho.l(2:end-1,3:end-1))/2 - RhoRef);
IR = [IR; ii(:)];  RR = [RR; rr(:)];

% assemble global coefficient matrix of stress divergence
KD = sparse(II,JJ,AA,NUM.NW+NUM.NU,NUM.NW+NUM.NU);
RD = sparse(IR,ones(size(IR)),RR,NUM.NW+NUM.NU,1);

%% Block [3,3] Assemble pressure scaling coefficients (P_l)
% reset global lists for vectorised assembly
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R

% % ================================================== % %
% coefficients multiplying pressure P
% % ================================================== % %
% boundary points
ii  = [indP(1,:).'; indP(end  ,:).']; % top & bottom
jj1 = ii;
jj2 = [indP(2,:).'; indP(end-1,:).'];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

ii  = [indP(:,1).'; indP(:,end  ).']; % left & right
jj1 = ii;
jj2 = [indP(:,2).'; indP(:,end-1).'];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

%internal points
ii = indP(2:end-1,2:end-1);
aa = NUM.cstab*NUM.dx*NUM.dz./MAT.Eta.s(2:end-1,2:end-1);
II = [II; ii(:)]; JJ = [JJ; ii(:)];    AA = [AA; aa(:)];  % P on stencil centre

% RHS
rr = zeros(size(ii));
IR = [IR; ii(:)];
RR = [RR; rr(:)];

% assemble global coefficient matrix of stress divergence
KP = sparse(II,JJ,AA,NUM.NP,NUM.NP);
RP = sparse(IR,ones(size(IR)),RR,NUM.NP,1);

%% Block [4,4] Assemble compaction pressure coefficients (P_cmp)
% reset global lists for vectorised assembly
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R

% % ================================================== % %
% coefficients multiplying pressure P
% % ================================================== % %
% boundary points
ii  = [indP(1,:).'; indP(end  ,:).']; % top & bottom
jj1 = ii;
jj2 = [indP(2,:).'; indP(end-1,:).'];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

ii  = [indP(:,1).'; indP(:,end  ).']; % left & right
jj1 = ii;
jj2 = [indP(:,2).'; indP(:,end-1).'];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

%internal points
ii = indP(2:end-1,2:end-1);
aa = - max(SOL.philim,SOL.phi(2:end-1,2:end-1))./MAT.Eta.s(2:end-1,2:end-1);
II = [II; ii(:)]; JJ = [JJ; ii(:)];    AA = [AA; aa(:)];  % P on stencil centre

% RHS
rr = zeros(size(ii));
IR = [IR; ii(:)];
RR = [RR; rr(:)];

% assemble global coefficient matrix of stress divergence
KC = sparse(II,JJ,AA,NUM.NP,NUM.NP);
RC = sparse(IR,ones(size(IR)),RR,NUM.NP,1);
%% Pressure divergences on the momentum and Darcy equations
% % ================================================== % %
% Block [1,3]; Liquid Pressure divergence (P_l)
% Block [2,3]; Liquid Pressure divergence (P_l)
% Block [1,4]; Compaction Pressure divergence (P_cmp)
% coefficients multiplying pressures P_l / P_cmp
% % ================================================== % %
% reset global lists for vectorised assembly
II  = [];   JJ  = [];   AA  = [];

% Z-Stokes equation
ii    = indW(2:end-1,2:end-1);
%         top              ||          bottom
jj1 = indP(2:end-2,2:end-1); jj2 = indP(3:end-1,2:end-1);

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1/NUM.dz];     % P one to the top
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1/NUM.dz];     % P one to the bottom

% X-Stokes equation
ii    = indU(2:end-1,2:end-1);
%         left             ||           right
jj1 = indP(2:end-1,2:end-2); jj2 = indP(2:end-1,3:end-1);

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1/NUM.dx];     % P one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1/NUM.dx];     % P one to the right

GG = sparse(II,JJ,AA,NUM.NW+NUM.NU,NUM.NP);

%% velocity divergences on the continuity and compaction equations
% % ================================================== % %
% Block [3,1]; solid/reference velocity divergence (WU_s)
% Block [4,1]; solid/reference Pressure divergence (WU_s)
% Block [3,2]; segregation velocity divergence (WU_seg)
% coefficients multiplying velocities WU_s / WU_seg
% % ================================================== % %
% reset global lists for vectorised assembly
II  = [];   JJ  = [];   AA  = [];

ii    = indP(2:end-1,2:end-1);
%          left U          ||           right U       ||           top W           ||          bottom W
jj1 = indU(2:end-1,1:end-1); jj2 = indU(2:end-1,2:end); jj3 = indW(1:end-1,2:end-1); jj4 = indW(2:end,2:end-1);

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1/NUM.dx];  % U one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1/NUM.dx];  % U one to the right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; aa(:)+1/NUM.dz];  % W one above
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; aa(:)-1/NUM.dz];  % W one below

DD = sparse(II,JJ,AA,NUM.NP,NUM.NW+NUM.NU);

%% Assemble combined matrices
V0 = sparse(NUM.NW+NUM.NU,  NUM.NW+NUM.NU);
G0 = sparse(NUM.NW+NUM.NU,  NUM.NP);
P0 = sparse(NUM.NP       ,  NUM.NP);

A = [KV   V0   GG   GG;
     V0   KD   GG   G0;
     DD   DD   KP   P0;
     DD   G0.' P0   KC ];
 
R = [RV; RD; RP; RC];

%% Scale system of equations (diagonal preconditioning)
S  =  sqrt(abs(diag(A)));
S  =  diag(sparse(1./S));

A  =  S*A*S;
R  =  S*R;

%% Solve linear system of equations for vx, vz, P
X = S*(A\R); % get solution vector

%% Read out solution
% map solution vector to 2D arrays
SOL.P.l  = reshape((X(indP(:))-mean(X(indP(:)))),NUM.nzP,NUM.nxP);                 % dynamic fluid pressure
SOL.Pt.l = SOL.P.l + RhoRef.*PHY.gz.*NUM.ZP + RhoRef.*PHY.gx.*NUM.XP;   % total fluid pressure
SOL.Pt.l = SOL.Pt.l - mean(SOL.Pt.l(2,:)) + RhoRef*NUM.dz/2*PHY.gz;     % adjust total fluid pressure to set surface pressure
SOL.W.s  = reshape(X(indW(:)),NUM.nzW,NUM.nxW);                         % solid z-velocity
SOL.U.s  = reshape(X(indU(:)),NUM.nzU,NUM.nxU);                         % solid x-velocity

% temporary plots for testing purpouses
    figure(100); clf;
    
    subplot(2,3,1)
    imagesc(NUM.xU,NUM.zU,SOL.U.s); hold on;
    colormap(subplot(2,3,1),cm2)
    axis ij equal tight;
    colorbar
    title('matrix x-velocity [ms^-^1]')
    
    subplot(2,3,2)
    imagesc(NUM.xW,NUM.zW,-SOL.W.s); hold on;
    colormap(subplot(2,3,2),cm2)
    axis ij equal tight;
    colorbar
    title('matrix z-velocity [ms^-^1]')
    
    subplot(2,3,3)
    imagesc(NUM.xP,NUM.zP,SOL.P.l); hold on;
    colormap(subplot(2,3,3),cm2)
    axis ij equal tight;
    colorbar
    title('liquid pressure [Pa]')
    
    subplot(2,3,4)
    imagesc(NUM.xU,NUM.zU,u_seg); hold on;
    colormap(subplot(2,3,4),cm2)
    axis ij equal tight;
    colorbar
    title('segr. x-velocity [ms^-^1]')
    
    subplot(2,3,5)
    imagesc(NUM.xW,NUM.zW,-w_seg); hold on;
    colormap(subplot(2,3,5),cm2)
    axis ij equal tight;
    colorbar
    title('segr. z-velocity [ms^-^1]')
    
    subplot(2,3,6)
    imagesc(NUM.xP,NUM.zP,P_cmp); hold on;
    colormap(subplot(2,3,6),cm2)
    axis ij equal tight;
    colorbar
    title('compct. pressure [Pa]')
    
    drawnow;
%%
SOL.P.s(2:end-1,2:end-1) = P_cmp./(1-SOL.phi(2:end-1,2:end-1)) + SOL.P.l(2:end-1,2:end-1);                                 % dynamic solid pressure
SOL.Pt.s = SOL.P.s + RhoRef.*PHY.gz.*NUM.ZP + RhoRef.*PHY.gx.*NUM.XP;   % total solid pressure
SOL.Pt.s = SOL.Pt.s - mean(SOL.Pt.s(2,:)) + RhoRef*NUM.dz/2*PHY.gz;     % adjust total solid pressure to set surface pressure
SOL.W.l(:,2:end-1) = w_seg./SOL.phiW(:,2:end-1) + SOL.W.s(:,2:end-1);
SOL.U.l(2:end-1,:) = u_seg./SOL.phiU(2:end-1,:) + SOL.U.s(2:end-1,:);

% PP  = reshape(X(indP(:)),NUM.nzP,NUM.nxP).*Pscale;               % dynamic pressure
% PPt = PP + RhoRef.*PHY.gz.*NUM.ZP + RhoRef.*PHY.gx.*NUM.XP;   % total pressure
% PPt = PPt - mean(PPt(2,:)) + RhoRef*NUM.dz/2*PHY.gz;       % adjust total pressure to set surface pressure
% WW  = reshape(X(indW(:)),NUM.nzW,NUM.nxW);                       % z-velocity
% UU  = reshape(X(indU(:)),NUM.nzU,NUM.nxU);                       % x-velocity

% project velocities to centre nodes
SOL.UP.s(:,2:end-1) = SOL.U.s(:,1:end-1)+SOL.U.s(:,2:end)./2;
SOL.WP.s(2:end-1,:) = SOL.W.s(1:end-1,:)+SOL.W.s(2:end,:)./2;

SOL.UP.l(:,2:end-1) = SOL.U.l(:,1:end-1)+SOL.U.l(:,2:end)./2;
SOL.WP.l(2:end-1,:) = SOL.W.l(1:end-1,:)+SOL.W.l(2:end,:)./2;

toc_solve = toc;  % stop clock system solution
fprintf(1,'       solution time %1.4f s \n\n',toc_solve);

