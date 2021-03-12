%main model loop, including functions below
time            = 0;     % initialise time loop
lapl_T0 = 0; adv_T0 = 0; dMdt0 = 0;
%% main loop
for ti = 1:nt
    
    Alpha_vx(:,2:nx1)   = (Alpha_mid(:,2:nx1)+Alpha_mid(:,1:nx))./2;
    Alpha_vz(2:nz1,:)   = (Alpha_mid(2:nz1,:)+Alpha_mid(1:nz,:))./2;
    k_vx(:,2:nx1)       = (Kappa_mid(:,2:nx1)+Kappa_mid(:,1:nx))./2;
    k_vz(2:nz1,:)       = (Kappa_mid(2:nz1,:)+Kappa_mid(1:nz,:))./2;
    
    % interpolate temperature onto staggered grid
    T_vx(:,2:nx1) = (T_mid(:,2:nx1)+T_mid(:,1:nx))./2;
    T_vz(2:nz1,:) = (T_mid(2:nz1,:)+T_mid(1:nz,:))./2;
    
    % Apply thermal boundary condition to interpolated nodes
    % upper boundary, constant temperature
    T_mid(1,:)   = 2*T_top   - T_mid(3,:);
    T_vx(1,:)    = 2*T_top   - T_vx(2,:);
    T_vz(1,:)   = 2*T_top   - T_vz(2,:);
    
    % lower boundary, constant temperature
    T_mid(end,:)= 2*(T_bot+300)   - T_mid(nz,:);
    T_vx(end,:)   = 2*(T_bot+300)   - T_vx(nz1,:);
    T_vz(end,:) = 2*(T_bot+300)   - T_vz(nz,:);
    %       T_mid(end,:) = T_bot;
    %       T_vx(end,:) = T_bot;
    %       T_vz(end,:) = T_bot;
    
    % left boundary, insulating boundary
    T_mid(2:nz1,2)  = T_mid(2:nz1,3);
    T_vx(:,1)       = T_vx(:,2);
    T_vz(1:nz1,1)   = T_vz(1:nz1,2);
    % right boundary, insulating boundary
    T_mid(2:nz1,nx1)= T_mid(2:nz1,nx);
    T_vx(:,nx1)     = T_vx(:,nx);
    T_vz(1:nz1,Nx)  = T_vz(1:nz1,nx1);
    %     plume beneath lithosphere
%     midx = fix(L/2);
%     indexP = find(abs(xp(1,:)) <= midx+ (L*0.25) & abs(xp(1,:)) >= midx - (L*0.25));
%     %     indexP = 1:Nx; indexvx = 1:Nx; indexvz = 1:Nx;
%     T_mid(end-1:end,indexP) = T_bot+300;
%     indexvx = find(abs(xvx(1,:)) <= midx+ (L*0.25) & abs(xvx(1,:)) >= midx - (L*0.25));
%     T_vx(end-1:end,indexvx) = T_bot+300;
%     indexvz = find(abs(xvz(1,:)) <= midx+ (L*0.25) & abs(xvz(1,:)) >= midx - (L*0.25));
%     T_vz(end-1:end,indexvz) = T_bot+300;
    
    
    [Rho_vx,Rho_vz,Rho_mid] = Thermal_expansion(Rho0,...
        T_mid,T_vx,T_vz,Alpha_vx,Alpha_vz,Alpha_mid,T_top);
    
    % ignore nonlinearity for the time being
    % % =============================================================
    % % Solve momentum and continuity equations
    % % =============================================================
    Pscale  = Eta_mantle/(dx+dz)*2; %pressure scaling coefficient, minimum viscosity
    [P_out,vx_out,vz_out,vx_mid,vz_mid] = stokes_continuity_noair(nx,nz,Nx,Nz,...
        nx1,nz1,dx,dz,Pscale,...
        Eta_out,Eta_mid,Rho_vx,Rho_vz,gx,gz,bctop,bcbottom,bcleft,bcright);
    
    % Define timestep dt_m for marker displacement
    dt1 = CFL * min( dx/2/max(abs(vx_out(:))) , dz/2/max(abs(vz_out(:))) ); %maximum timestep for advection
    dt2 = CFL *dx/2*dz/2/max(max(abs(Kappa_mantle./Rho_mid(2:nz1,2:nx1)./Cp_mid(2:nz1,2:nx1)))); %maximum timestep for T diffusion
    dt = min(dt2,dt1);
    
    %update temperature by diffusion and advection
    run('Update_Temperature')
%     [Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,T_out,lapl_T0,adv_T0,dMdt0,Material] =...
%         Update_Temperature(Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,nx1,nz1,dx,dz,...
%         vx_out,vz_out,vx_mid,vz_mid,Eta_out,Eta_mid,Rho_vz,Alpha_mid,T_mid,gz,...
%         Nx,Nz,k_vx,k_vz,Rho_mid,Cp_mid,Hr,T_top,dt,Tsolver,Material,ti,lapl_T0,adv_T0,dMdt0,AdvRegime);
    
    time = dt+time;
    
    if ~mod(ti,20) || ti==1
        run('plotfigures2')
        %     profile report
        %         saveas(figure(1),['../out/', RunID, '/Ti', num2str(ti), '.jpg'])
        drawnow
        ti
    end
    
    T_mid = T_out;
    
end

% profile report