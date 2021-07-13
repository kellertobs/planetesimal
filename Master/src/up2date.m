% planetesimal: update material properties and stress/strain-rate fields

% print update header
fprintf(1,'  ---  update materials & deformation \n');
tic;  % start clock on update


%% update thermal parameters (currently held constant)


%% update permeability
MAT.k   = PHY.k0*max(SOL.philim,SOL.phi).^3;
% interpolate to staggered vz nodes
MAT.kW  = (MAT.k(1:end-1,:)+MAT.k(2:end,:))/2;
% interpolate to staggered vx nodes
MAT.kU  = (MAT.k(:,1:end-1)+MAT.k(:,2:end))/2;


%% update T-dependent density
MAT.Rho.s = PHY.Rho0.s.*(1 - MAT.aT.*(SOL.T-SOL.T0));    % solid density on centre nodes
MAT.Rho.l = PHY.Rho0.l.*(1 - MAT.aT.*(SOL.T-SOL.T0));    % liquid
MAT.Rho.t = SOL.phi.*MAT.Rho.l + (1-SOL.phi).*MAT.Rho.s; % total density
% interpolate to staggered vz nodes
MAT.RhoW.s = (MAT.Rho.s(1:end-1,:)+MAT.Rho.s(2:end,:))/2;
MAT.RhoW.l = (MAT.Rho.l(1:end-1,:)+MAT.Rho.l(2:end,:))/2;
% interpolate to staggered vx nodes
MAT.RhoU.s = (MAT.Rho.s(:,1:end-1)+MAT.Rho.s(:,2:end))/2;
MAT.RhoU.l = (MAT.Rho.l(:,1:end-1)+MAT.Rho.l(:,2:end))/2;


%% update physical time step
Vs      = abs([SOL.U.s(:);SOL.W.s(:)]);
dtadvns = min( NUM.dx/2/max(Vs) , NUM.dz/2/max(Vs) );
Vl      = abs([SOL.U.l(:);SOL.W.l(:)]);
dtadvnl = min( NUM.dx/2/max(Vl) , NUM.dz/2/max(Vl) );  % stable time step for T advection
dtadvn  = min(dtadvnl,dtadvns);

kappas  = MAT.kT(:)./MAT.Rho.s(:)./MAT.Cp(:);
dtdiffs = (max(NUM.dx,NUM.dz)/2)^2 / max(kappas);     % stable time step for T diffusion
kappal  = MAT.kT(:)./MAT.Rho.l(:)./MAT.Cp(:);
dtdiffl = (max(NUM.dx,NUM.dz)/2)^2 / max(kappal);
dtdiff  = min(dtdiffl,dtdiffs);

NUM.dt = NUM.CFL * min(dtdiff,dtadvn);              % fraction of minimum stable time step


%% update T-dependent viscosity (to be implemented)
MAT.Eta.s  = MAT.Eta.s;
MAT.Eta.l  = MAT.Eta.l;
% viscosity on corner nodes
MAT.EtaC.s = (MAT.Eta.s(1:end-1,1:end-1) ...
    +  MAT.Eta.s(2:end  ,1:end-1) ...
    +  MAT.Eta.s(1:end-1,2:end  ) ...
    +  MAT.Eta.s(2:end  ,2:end  ))/4;
MAT.EtaC.l = (MAT.Eta.l(1:end-1,1:end-1) ...
    +  MAT.Eta.l(2:end  ,1:end-1) ...
    +  MAT.Eta.l(1:end-1,2:end  ) ...
    +  MAT.Eta.l(2:end  ,2:end  ))/4;
% interpolate to staggered vz nodes
MAT.EtaW.s = (MAT.Eta.s(1:end-1,:)+MAT.Eta.s(2:end,:))/2;
MAT.EtaW.l = (MAT.Eta.l(1:end-1,:)+MAT.Eta.l(2:end,:))/2;
% interpolate to staggered vx nodes
MAT.EtaU.s = (MAT.Eta.s(:,1:end-1)+MAT.Eta.s(:,2:end))/2;
MAT.EtaU.l = (MAT.Eta.l(:,1:end-1)+MAT.Eta.l(:,2:end))/2;

%% update liquid fraction
SOL.phi     = SOL.phi;
SOL.phiC = (SOL.phi(1:end-1,1:end-1) ...
    +  SOL.phi(2:end  ,1:end-1) ...
    +  SOL.phi(1:end-1,2:end  ) ...
    +  SOL.phi(2:end  ,2:end  ))/4;
% interpolate to staggered vz nodes
SOL.phiW = (SOL.phi(1:end-1,:)+SOL.phi(2:end,:))/2;
% interpolate to staggered vx nodes
SOL.phiU = (SOL.phi(:,1:end-1)+SOL.phi(:,2:end))/2;


%% update strain-rate components
% % get volumetric strain-rate (velocity divergence)
DEF.ups(2:end-1,2:end-1) = diff(SOL.U.s(2:end-1,:),1,2)./NUM.dx ...
    + diff(SOL.W.s(:,2:end-1),1,1)./NUM.dz;      % velocity divergence
DEF.ups([1 end],:)       = DEF.ups([2 end-1],:);                    % apply boundary conditions
DEF.ups(:,[1 end])       = DEF.ups(:,[2 end-1]);

% % get deviatoric strain rates
DEF.exx(:,2:end-1) = diff(SOL.U.s,1,2)./NUM.dx ...
    - DEF.ups(:,2:end-1)./3;                         % x-normal strain rate
DEF.exx([1 end],:) = DEF.exx([2 end-1],:);                          % apply boundary conditions
DEF.exx(:,[1 end]) = DEF.exx(:,[2 end-1]);
%
DEF.ezz(2:end-1,:) = diff(SOL.W.s,1,1)./NUM.dz ...
    - DEF.ups(2:end-1,:)./3;                         % z-normal strain rate
DEF.ezz([1 end],:) = DEF.ezz([2 end-1],:);                         % apply boundary conditions
DEF.ezz(:,[1 end]) = DEF.ezz(:,[2 end-1]);
%
DEF.exz            = (diff(SOL.U.s,1,1)./NUM.dz ...
    +  diff(SOL.W.s,1,2)./NUM.dx)/2;                   % shear strain rate
%
% update strain-rate magnitude
DEF.eII(2:end-1,2:end-1) = (  (DEF.exx(2:end-1,2:end-1).^2 + DEF.ezz(2:end-1,2:end-1).^2 ...
    + 2.*(DEF.exz(1:end-1,1:end-1).^2 + DEF.exz(2:end,1:end-1).^2 ...
    +     DEF.exz(1:end-1,2:end  ).^2 + DEF.exz(2:end,2:end  ).^2)/4 )/2).^0.5 + 1e-16;
DEF.eII([1 end],:) = DEF.eII([2 end-1],:);                                 % apply boundaries
DEF.eII(:,[1 end]) = DEF.eII(:,[2 end-1]);


%% update stress components
DEF.txx = MAT.Eta.s  .* DEF.exx;                                      % x-normal stress
DEF.tzz = MAT.Eta.s  .* DEF.ezz;                                      % z-normal stress
DEF.txz = MAT.EtaC.s .* DEF.exz;                                      % xz-shear stress

% update strain-rate magnitude


DEF.tII(2:end-1,2:end-1) = (  (DEF.txx(2:end-1,2:end-1).^2 + DEF.tzz(2:end-1,2:end-1).^2 ...
    + 2.*(DEF.txz(1:end-1,1:end-1).^2 + DEF.txz(2:end,1:end-1).^2 ...
    +     DEF.txz(1:end-1,2:end  ).^2 + DEF.txz(2:end,2:end  ).^2)/4 )/2).^0.5 + 1e-16;
DEF.tII([1 end],:) = DEF.tII([2 end-1],:);                                 % apply boundaries
DEF.tII(:,[1 end]) = DEF.tII(:,[2 end-1]);

%% update heat source fields
% update shear heating
SOL.Hs = 2.*DEF.eII.*DEF.tII...
    + (SOL.UP.seg.^2).*MAT.Eta.l./MAT.k + (SOL.WP.seg.^2).*MAT.Eta.l./MAT.k ...
    + (SOL.Pcmp.^2).*SOL.phi./MAT.Eta.s;

% update adiabatic heating
SOL.Ha = (1-SOL.phi(2:end-1,2:end-1)).*SOL.T(2:end-1,2:end-1)    .*MAT.aT(2:end-1,2:end-1)...
    .*   (1-SOL.phi(2:end-1,2:end-1)).*MAT.Rho.s(2:end-1,2:end-1)...
    .*   (PHY.gz.*SOL.WP.s(2:end-1,2:end-1) + PHY.gx.*SOL.UP.s(2:end-1,2:end-1))...
    +     SOL.phi(2:end-1,2:end-1)   .*SOL.T(2:end-1,2:end-1).*MAT.aT(2:end-1,2:end-1)...
    .*    SOL.phi(2:end-1,2:end-1)   .*MAT.Rho.l(2:end-1,2:end-1)...
    .*   (PHY.gz.*SOL.WP.l(2:end-1,2:end-1) + PHY.gx.*SOL.UP.l(2:end-1,2:end-1));


toc_update = toc;
fprintf(1,'       update time %1.4f s \n\n',toc_update);