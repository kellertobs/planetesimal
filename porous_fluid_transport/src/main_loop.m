%main model loop, including functions below
time            = 0;     % initialise time loop
ti = 0;
lapl_T0 = 0; adv_T0 = 0; dMdt0 = 0;
Alpha_vx(:,2:nx1)   = (Alpha_mid(:,2:nx1)+Alpha_mid(:,1:nx))./2;
Alpha_vz(2:nz1,:)   = (Alpha_mid(2:nz1,:)+Alpha_mid(1:nz,:))./2;
k_vx(:,2:nx1)       = (Kappa_mid(:,2:nx1)+Kappa_mid(:,1:nx))./2;
k_vz(2:nz1,:)       = (Kappa_mid(2:nz1,:)+Kappa_mid(1:nz,:))./2;


%% main time stepping loop
while time < max_time && ti < nt
    
    Temperature_boundary;
    
    % update T-dependent density
    [Rho_vx,Rho_vz,Rho_mid] = Thermal_expansion(Rho_mantle,...
    T_mid,T_vx,T_vz,Alpha_vx,Alpha_vz,Alpha_mid,T_top);

    % ignore nonlinearity for the time being
    % % =============================================================
    % % Solve momentum and continuity equations
    % % =============================================================
    
    % only solve stokes every one gridstep advected
    if ti == 0 || ~mod(ti,round(2/CFL))
        [P_out,vx_out,vz_out,vx_mid,vz_mid] = stokes_continuity_noair(nx,nz,nx2,nz2,nx1,nz1,dx,dz,...
            Eta_out,Eta_mid,Rho_vx,Rho_vz,gx,gz,zp2d,bctop,bcbottom,bcleft,bcright,cstab);
    end
    
    % Update timestep dt for advection-diffusion problem
    dt1 = min( dx/2/max(abs(vx_out(:))) , dz/2/max(abs(vz_out(:))) ); %maximum timestep for advection
    dt2 = dx/2*dz/2/max(max(abs(Kappa_mantle./Rho_mid(2:nz1,2:nx1)./Cp_mid(2:nz1,2:nx1)))); %maximum timestep for T diffusion
    dt  = CFL * min(dt2,dt1);
    
    %update temperature by diffusion and advection
    Update_Temperature;
    T_mid = T_out;
    
    if ~mod(ti,round(nx/10/CFL*2)) || ti==0 % after 5 gridsteps advection
%     if ~mod(ti,20) || ti==0
        plotfigures2;
        %     profile report
%                 saveas(figure(1),['../out/', RunID, '/Ti', num2str(ti), '.jpg'])
        drawnow
        ti
    end
    
    % increment time
    ti   = ti+1;
    time = dt+time;
end

% profile report