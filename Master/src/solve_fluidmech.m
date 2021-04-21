% planetesimal: solve fluid mechanics equations

% print solver header
fprintf(1,'  ---  solve fluid mechanics equations \n');
    

%% prepare auxiliary variables

% get characteristic scales
Pscale = geomean(MAT.Eta(:))/(NUM.dx*NUM.dz); % pressure scaling coefficient
RhoRef = mean(MAT.Rho(:));                    % mean density for lithostatic pressure

% get mapping arrays
indU = NUM.MapU;
indW = NUM.MapW;
indP = NUM.MapP;

% initialise global lists for vectorised assembly
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R

tic;  % start clock system assembly


%% assemble coefficients of z momentum equation for z-velocity W

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
EtaC1 = MAT.EtaC(2:end-1,1:end-1); EtaC2 = MAT.EtaC(2:end-1,2:end  ); 
EtaP1 = MAT.Eta (2:end-2,2:end-1); EtaP2 = MAT.Eta (3:end-1,2:end-1);

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

% coefficients multiplying pressures P
%         top              ||          bottom
jj1 = indP(2:end-2,2:end-1); jj2 = indP(3:end-1,2:end-1);

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+Pscale/NUM.dz];     % P one to the top
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-Pscale/NUM.dz];     % P one to the bottom

% RHS vector
rr = zeros(size(ii)) - PHY.gz*( (MAT.Rho(2:end-2,2:end-1)+MAT.Rho(3:end-1,2:end-1))/2 - RhoRef);
IR = [IR; ii(:)];  RR = [RR; rr(:)];


%% assemble coefficients of x momentum equation for x-velocity U

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
EtaC1 = MAT.EtaC(1:end-1,2:end-1); EtaC2 = MAT.EtaC(2:end  ,2:end-1); 
EtaP1 = MAT.Eta (2:end-1,2:end-2); EtaP2 = MAT.Eta (2:end-1,3:end-1);

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

% coefficients multiplying pressures P
%         left             ||           right
jj1 = indP(2:end-1,2:end-2); jj2 = indP(2:end-1,3:end-1);

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+Pscale/NUM.dx];     % P one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-Pscale/NUM.dx];     % P one to the right

% RHS vector
rr = zeros(size(ii)) - PHY.gx*( (MAT.Rho(2:end-1,2:end-2)+MAT.Rho(2:end-1,3:end-1))/2 - RhoRef);
IR = [IR; ii(:)];  RR = [RR; rr(:)];


%% assemble coefficients of continuity equation for P
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

% coefficients multiplying pressure P
aa = NUM.cstab*NUM.dx*NUM.dz./MAT.Eta(2:end-1,2:end-1);
II = [II; ii(:)]; JJ = [JJ; ii(:)];    AA = [AA; aa(:)];  % P on stencil centre

% coefficients multiplying velocities U, W
%          left U          ||           right U       ||           top W           ||          bottom W
jj1 = indU(2:end-1,1:end-1); jj2 = indU(2:end-1,2:end); jj3 = indW(1:end-1,2:end-1); jj4 = indW(2:end,2:end-1);

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)-1/NUM.dx];  % U one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/NUM.dx];  % U one to the right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; aa(:)-1/NUM.dz];  % W one above
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; aa(:)+1/NUM.dz];  % W one below

% RHS
rr = zeros(size(ii));
IR = [IR; ii(:)];
RR = [RR; rr(:)];


%% assemble global coefficient matrix and right-hand side vector
A = sparse(II,JJ,AA,NUM.NDOF,NUM.NDOF);
R = sparse(IR,ones(size(IR)),RR,NUM.NDOF,1);


%% Scale system of equations (diagonal preconditioning)
S  =  sqrt(abs(diag(A)));
S  =  diag(sparse(1./S));

A  =  S*A*S;
R  =  S*R;

toc_assmb = toc;  % stop clock system assembly
fprintf(1,'       assembly time %1.4f s \n',toc_assmb);

tic;  %start clock system solution


%% Solve linear system of equations for vx, vz, P
X = S*(A\R); % get solution vector


%% Read out solution
% map solution vector to 2D arrays
SOL.P  = reshape(X(indP(:)),NUM.nzP,NUM.nxP).*Pscale;               % dynamic pressure
SOL.Pt = SOL.P + RhoRef.*PHY.gz.*NUM.ZP + RhoRef.*PHY.gx.*NUM.XP;   % total pressure
SOL.Pt = SOL.Pt - mean(SOL.Pt(2,:)) + RhoRef*NUM.dz/2*PHY.gz;       % adjust total pressure to set surface pressure
SOL.W  = reshape(X(indW(:)),NUM.nzW,NUM.nxW);                       % z-velocity
SOL.U  = reshape(X(indU(:)),NUM.nzU,NUM.nxU);                       % x-velocity

% project velocities to centre nodes
SOL.UP(:,2:end-1) = SOL.U(:,1:end-1)+SOL.U(:,2:end)./2;
SOL.WP(2:end-1,:) = SOL.W(1:end-1,:)+SOL.W(2:end,:)./2;

toc_solve = toc;  % stop clock system solution
fprintf(1,'       solution time %1.4f s \n\n',toc_solve);

