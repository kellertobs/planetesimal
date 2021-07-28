% planetesimal: solve self-gravity Poisson equation

% print solver header
fprintf(1,'  ---  solve self-gravity equation \n');
tic;


%% prepare auxiliary variables

% set gravitational constant
Grav = 6.67408e-11;   % [m3/(kg s2)]
% Grav = 10;
% get mapping arrays

ind = NUM.MapP;

% initialise global lists for vectorised assembly
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R

tic;  % start clock system assembly

%% assemble coefficients of poisson equation for gravitational potential G
%internal points
ii = ind(2:end-1,2:end-1);
% coefficients multiplying pressure G
%          left U          ||           right U       ||           top W           ||          bottom W
jj1 = ind(2:end-1,1:end-2); jj2 = ind(2:end-1,3:end); jj3 = ind(1:end-2,2:end-1); jj4 = ind(3:end,2:end-1);

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; ii(:)];    AA = [AA; aa(:)-2/NUM.dx.^2-2/NUM.dz.^2];  % G on stencil centre
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1/NUM.dx.^2];  % G one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/NUM.dx.^2];  % G one to the right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; aa(:)+1/NUM.dz.^2];  % G one above
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; aa(:)+1/NUM.dz.^2];  % G one below

% RHS
rr = 4*pi*Grav*MAT.Rho.t(ii);
IR = [IR; ii(:)];
RR = [RR; rr(:)];


%% assemble global coefficient matrix and right-hand side vector
A = sparse(II,JJ,AA,NUM.NP,NUM.NP);
R = sparse(IR,ones(size(IR)),RR,NUM.NP,1);


%% set equi-potential distance at D/2
ii = ind(NUM.PHI<0); ii = ii(:);
A(ii,:) = 0;
A(sub2ind(size(A),ii,ii)) = 1;
R(ii,1) = 0;


%% Solve linear system of equations for vx, vz, P
X = A\R; % get solution vector

% map solution to 2D array
SOL.G  = reshape(X(ind(:)),NUM.nzP,NUM.nxP);

% get gravity components gx, gz from gradient of potential G
PHY.gz = -2*diff(SOL.G,1,1)./NUM.dz;
PHY.gx = -2*diff(SOL.G,1,2)./NUM.dx;

% project velocities to centre nodes
PHY.gxP(:,2:end-1) = PHY.gx(:,1:end-1)+PHY.gx(:,2:end)./2;
PHY.gzP(2:end-1,:) = PHY.gz(1:end-1,:)+PHY.gz(2:end,:)./2;
PHY.gP = (PHY.gxP.^2 + PHY.gzP.^2).^0.5;
PHY.gP([1 end],:) = 0;
PHY.gP(:,[1 end]) = 0;

toc_solve = toc;  % stop clock system solution
fprintf(1,'       solution time %1.4f s \n\n',toc_solve);

figure()
subplot(1,3,1)
imagesc(NUM.xU,NUM.zU,PHY.gx); colorbar; axis ij image; title('gx')

subplot(1,3,2)
imagesc(NUM.xW,NUM.zW,PHY.gz); colorbar; axis ij image; title('gz')

subplot(1,3,3)
imagesc(NUM.xW,NUM.zW,SOL.G); colorbar; axis ij image; title('g')
