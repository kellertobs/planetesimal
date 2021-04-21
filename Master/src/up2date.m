% planetesimal: update material properties and stress/strain-rate fields

% print update header
fprintf(1,'  ---  update materials & deformation \n');
tic;  % start clock on update


%% update thermal parameters (currently held constant)


%% update T-dependent density
MAT.Rho = PHY.Rho0.*(1 - MAT.aT.*(SOL.T-SOL.T0));                   % density on centre nodes


%% update physical time step
V      = abs([SOL.U(:);SOL.W(:)]);
dtadvn = min( NUM.dx/2/max(V) , NUM.dz/2/max(V) );  % stable time step for T advection

kappa  = MAT.kT(:)./MAT.Rho(:)./MAT.Cp(:);
dtdiff = (max(NUM.dx,NUM.dz)/2)^2 / max(kappa);     % stable time step for T diffusion

NUM.dt = NUM.CFL * min(dtdiff,dtadvn);              % fraction of minimum stable time step


%% update T-dependent viscosity (to be implemented)
MAT.Eta  = MAT.Eta;                                                 % viscosity on centre nodes
MAT.EtaC = (MAT.Eta(1:end-1,1:end-1) ...
         +  MAT.Eta(2:end  ,1:end-1) ...
         +  MAT.Eta(1:end-1,2:end  ) ...
         +  MAT.Eta(2:end  ,2:end  ))/4;                            % interpolate to corner nodes


%% update strain-rate components
% get volumetric strain-rate (velocity divergence)
DEF.ups(2:end-1,2:end-1) = diff(SOL.U(2:end-1,:),1,2)./NUM.dx ...
                         + diff(SOL.W(:,2:end-1),1,1)./NUM.dz;      % velocity divergence
DEF.ups([1 end],:)       = DEF.ups([2 end-1],:);                    % apply boundary conditions
DEF.ups(:,[1 end])       = DEF.ups(:,[2 end-1]);

% get deviatoric strain rates
DEF.exx(:,2:end-1) = diff(SOL.U,1,2)./NUM.dx ...
                   - DEF.ups(:,2:end-1)./3;                         % x-normal strain rate
DEF.exx([1 end],:) = DEF.exx([2 end-1],:);                          % apply boundary conditions
DEF.exx(:,[1 end]) = DEF.exx(:,[2 end-1]);

DEF.ezz(2:end-1,:) = diff(SOL.W,1,1)./NUM.dz ...
                   - DEF.ups(2:end-1,:)./3;                         % z-normal strain rate
DEF.ezz([1 end],:) =  DEF.ezz([2 end-1],:);                         % apply boundary conditions
DEF.ezz(:,[1 end]) =  DEF.ezz(:,[2 end-1]);

DEF.exz            = (diff(SOL.U,1,1)./NUM.dz ...
                   +  diff(SOL.W,1,2)./NUM.dx)/2;                   % shear strain rate

% update strain-rate magnitude
DEF.eII(2:end-1,2:end-1) = (  (DEF.exx(2:end-1,2:end-1).^2 + DEF.ezz(2:end-1,2:end-1).^2 ...
                         + 2.*(DEF.exz(1:end-1,1:end-1).^2 + DEF.exz(2:end,1:end-1).^2 ...
                         +     DEF.exz(1:end-1,2:end  ).^2 + DEF.exz(2:end,2:end  ).^2)/4 )/2).^0.5 + 1e-16;
DEF.eII([1 end],:) = DEF.eII([2 end-1],:);                                 % apply boundaries
DEF.eII(:,[1 end]) = DEF.eII(:,[2 end-1]);


%% update stress components
DEF.txx = MAT.Eta  .* DEF.exx;                                      % x-normal stress
DEF.tzz = MAT.Eta  .* DEF.ezz;                                      % z-normal stress
DEF.txz = MAT.EtaC .* DEF.exz;                                      % xz-shear stress

% update strain-rate magnitude
DEF.tII(2:end-1,2:end-1) = (  (DEF.txx(2:end-1,2:end-1).^2 + DEF.tzz(2:end-1,2:end-1).^2 ...
                         + 2.*(DEF.txz(1:end-1,1:end-1).^2 + DEF.txz(2:end,1:end-1).^2 ...
                         +     DEF.txz(1:end-1,2:end  ).^2 + DEF.txz(2:end,2:end  ).^2)/4 )/2).^0.5 + 1e-16;
DEF.tII([1 end],:) = DEF.tII([2 end-1],:);                                 % apply boundaries
DEF.tII(:,[1 end]) = DEF.tII(:,[2 end-1]);


%% update heat source fields
% update shear heating
SOL.Hs = 2.*DEF.eII.*DEF.tII;

% update adiabatic heating
SOL.Ha = (SOL.WP.*PHY.gz + SOL.UP.*PHY.gx) .* MAT.Rho.*MAT.aT.*SOL.T;


toc_update = toc;
fprintf(1,'       update time %1.4f s \n\n',toc_update);