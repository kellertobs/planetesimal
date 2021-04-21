%main model loop, including functions below
time            = 0;     % initialise time loop
ti = 0;
lapl_T0 = 0; adv_T0 = 0; dMdt0 = 0;
Alpha_vx(:,2:nx1)   = (Alpha_mid(:,2:nx1)+Alpha_mid(:,1:nx))./2;
Alpha_vz(2:nz1,:)   = (Alpha_mid(2:nz1,:)+Alpha_mid(1:nz,:))./2;
k_vx(:,2:nx1)       = (Kappa_mid(:,2:nx1)+Kappa_mid(:,1:nx))./2;
k_vz(2:nz1,:)       = (Kappa_mid(2:nz1,:)+Kappa_mid(1:nz,:))./2;
%% main loop
while time < max_time && ti<nt
    % for ti = 0:nt
    
    % interpolate temperature onto staggered grid
    T_vx(:,2:nx1) = (T_mid(:,2:nx1)+T_mid(:,1:nx))./2;
    T_vz(2:nz1,:) = (T_mid(2:nz1,:)+T_mid(1:nz,:))./2;
    
    % Apply thermal boundary condition to interpolated nodes
    % upper boundary, constant temperature
    T_vx(1,:)    = 2*T_top   - T_vx(2,:);
    T_vz(1,:)   = 2*T_top   - T_vz(2,:);
    
    % lower boundary, constant temperature
    T_vx(end,:)   = 2*(T_bot+300)   - T_vx(nz1,:);
    T_vz(end,:) = 2*(T_bot+300)   - T_vz(nz,:);
    %       T_vx(end,:) = T_bot;
    %       T_vz(end,:) = T_bot;
    
    % left boundary, insulating boundary
    T_vx(:,1)       = T_vx(:,2);
    T_vz(1:nz1,1)   = T_vz(1:nz1,2);
    % right boundary, insulating boundary
    T_vx(:,nx1)     = T_vx(:,nx);
    T_vz(1:nz1,Nx)  = T_vz(1:nz1,nx1);
%         plume beneath lithosphere
        midx = fix(L/2);
        indexP = find(abs(xp(1,:)) <= midx+ (L*0.25) & abs(xp(1,:)) >= midx - (L*0.25));
        %     indexP = 1:Nx; indexvx = 1:Nx; indexvz = 1:Nx;
        T_mid(end-1:end,indexP) = T_bot+300;
        indexvx = find(abs(xvx(1,:)) <= midx+ (L*0.25) & abs(xvx(1,:)) >= midx - (L*0.25));
        T_vx(end-1:end,indexvx) = T_bot+300;
        indexvz = find(abs(xvz(1,:)) <= midx+ (L*0.25) & abs(xvz(1,:)) >= midx - (L*0.25));
        T_vz(end-1:end,indexvz) = T_bot+300;
    [Rho_vx,Rho_vz,Rho_mid] = Thermal_expansion(Rho0,...
        T_mid,T_vx,T_vz,Alpha_vx,Alpha_vz,Alpha_mid,T_top);
    % ignore nonlinearity for the time being
    % % =============================================================
    % % Solve momentum and continuity equations
    % % =============================================================
    Pscale  = Eta_mantle/(dx+dz)*2; %pressure scaling coefficient, minimum viscosity
    cc = 1e-1;
    if ti ==0        
        [P_out,vx_out,vz_out,vx_mid,vz_mid] = stokes_continuity_noair(nx,nz,Nx,Nz,...
            nx1,nz1,dx,dz,Pscale,...
            Eta_out,Eta_mid,Rho_vx,Rho_vz,gx,gz,bctop,bcbottom,bcleft,bcright,cc);
        
    elseif (~mod(ti,round(1/CFL*2)))
        [P_out,vx_out,vz_out,vx_mid,vz_mid] = stokes_continuity_noair(nx,nz,Nx,Nz,...
            nx1,nz1,dx,dz,Pscale,...
            Eta_out,Eta_mid,Rho_vx,Rho_vz,gx,gz,bctop,bcbottom,bcleft,bcright,cc);
        
    end
    
%     P_a = zp2d.*Rho0.*gz;
%     P_out = full(P_out);
%     P_res = P_out-P_a;
%     P_res = P_res(2:end-1,2:end-1);
%     maxdP = max(abs(P_res(:)));
%     figure(7)
%     subplot(1,3,1)
%     imagesc(xp(2:end-1),zp(2:end-1),P_out(2:end-1,2:end-1))
%     colorbar
%     axis ij image
%     subplot(1,3,2)
%     imagesc(xp(2:end-1),zp(2:end-1),P_res)
%     colorbar
%     axis ij image
%     subplot(1,3,3)
%     imagesc(xp(2:end-1),zp(2:end-1),P_a(2:end-1,2:end-1))
%     colorbar
%     axis ij image
% %     
%     figure(6)
%     subplot(1,2,2)
%     loglog(cc,maxdP,'*')
%     xlabel('log(cc)')
%     ylabel('||res||')
%     hold on
    
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
    
    
    if ~mod(ti,round(nx/10/CFL*2)) || ti==0 % after 5 gridsteps advection
%     if ~mod(ti,20) || ti==0
        run('plotfigures2')
        %     profile report
                saveas(figure(1),['../out/', RunID, '/Ti', num2str(ti), '.jpg'])
        drawnow
        ti
    end
    ti = ti+1;
    T_mid = T_out;
    time = dt+time;
end

% profile report