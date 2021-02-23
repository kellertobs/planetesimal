%main model loop, including functions below
time            = 0;     % initialise time loop

%% main loop
for ti = 1:nt

    % % =============================================================
    % % update grid
    % % =============================================================
    % interpolate markes onto eulerian grid
%     [Eta_out,Eta_mid,k_vz,k_vx,Cp_mid,T_mid,Alpha_mid,Hr,Material,Alpha_vx,Alpha_vz,T_vx,T_vz]...
%     = Marker2grid(marknum,nx,nx1,Nx,nz,nz1,Nz,dx,dz,xm,zm,...
%     x,z,xvx,zvx,xvz,zvz,xp,zp,Etam,Kappam,Tm,Cpm,Alpham,Hrm,Mtype);

    Alpha_vx(:,2:nx1) = (Alpha_mid(:,2:nx1)+Alpha_mid(:,1:nx))./2;
    Alpha_vz(2:nz1,:) = (Alpha_mid(2:nz1,:)+Alpha_mid(1:nz,:))./2;
    k_vx(:,2:nx1) = (Kappa_mid(:,2:nx1)+Kappa_mid(:,1:nx))./2;
    k_vz(2:nz1,:) = (Kappa_mid(2:nz1,:)+Kappa_mid(1:nz,:))./2;
    
    % interpolate temperature onto staggered grid
    T_vx(:,2:nx1) = (T_mid(:,2:nx1)+T_mid(:,1:nx))./2;
    T_vz(2:nz1,:) = (T_mid(2:nz1,:)+T_mid(1:nz,:))./2;   
    
    % Apply thermal boundary condition to interpolated nodes
    % upper boundary, constant temperature 
    T_mid(2,2:nx)   = 2*T_top   - T_mid(3,2:nx);
    T_vx(1,2:nx)    = 2*T_top   - T_vx(2,2:nx);
    T_vz(1,2:nx1)   = 2*T_top   - T_vz(2,2:nx1);
    % lower boundary, constant temperature
    T_mid(nz1,2:nx1)= 2*(T_bot)   - T_mid(nz,2:nx1);
    T_vx(Nz,2:nx)   = 2*(T_bot)   - T_vx(nz1,2:nx);
    T_vz(nz1,2:nx1) = 2*(T_bot)   - T_vz(nz,2:nx1);    
    
%     T_mid(nz1,2:nx1)= 2*(T_bot+300)   - T_mid(nz,2:nx1);
%     T_vx(nz1:Nz,2:nx)   = 2*(T_bot+300)   - T_vx(nz1,2:nx);
%     T_vz(nz1:Nz,2:nx1) = 2*(T_bot+300)   - T_vz(nz,2:nx1);  
    
    % left boundary, insulating boundary
    T_mid(2:nz1,2)  = T_mid(2:nz1,3);
    T_vx(:,1)       = T_vx(:,2);
    T_vz(1:nz1,1)   = T_vz(1:nz1,2);
    % right boundary, insulating boundary
    T_mid(2:nz1,nx1)= T_mid(2:nz1,nx);
    T_vx(:,nx1)     = T_vx(:,nx);
    T_vz(1:nz1,Nx)  = T_vz(1:nz1,nx1);
%     plume beneath lithosphere
    midx = fix(L/2); 
    
    indexP = find(abs(xp(1,:)) <= midx+ (L*0.23) & abs(xp(1,:)) >= midx - (L*0.23));
    T_mid(nz1,indexP) = T_bot+300;
    indexvx = find(abs(xvx(1,:)) <= midx+ (L*0.23) & abs(xvx(1,:)) >= midx - (L*0.23));
    T_vx(nz:Nz,indexvx) = T_bot+300;
    indexvz = find(abs(xvz(1,:)) <= midx+ (L*0.23) & abs(xvz(1,:)) >= midx - (L*0.23));
    T_vz(nz:Nz,indexvz) = T_bot+300;
   

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
        dt1 = CFL * min( dx/2/max(abs(vx_out(:))) , dz/2/max(abs(vz_out(:))) );
%         dt2 = (dx^2)/2/max(abs(Kappa_mantle*RhoCp_mantle));
%         dt = min(dt1,dt2);
%         dt      = dt*dtkoef;
%         maxvx   = max(max(abs(vx_out)));
%         maxvz   = max(max(abs(vz_out)));
%         if(dt*maxvx>dxzmax*dx)
%             dt  = dxzmax*dx/maxvx;
%         end
%         if(dt*maxvz>dxzmax*dz)
%             dt  = dxzmax*dz/maxvz;
%         end
%         
    dt2 = CFL *dx*dz/max(max(abs(Kappa_mantle./Rho_mid(2:nz1,2:nx1)./Cp_mid(2:nz1,2:nx1))));
    dt = min(dt2,dt1);

    % % =============================================================
    % % Update eulerian temperature, correct dt for maximum temperature
    % % =============================================================
    [Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,T_diff,dT] =...
    Update_Temperature(Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,nx1,nz1,dx,dz,...
    vx_out,vz_out,Eta_out,Eta_mid,Rho_vz,Alpha_mid,T_mid,gz,...
    Nx,Nz,k_vx,k_vz,Rho_mid,Cp_mid,Hr,T_top,dt);
    
    T_diff = full(T_diff);
    
%     Output = 'Fromm'
    [dTdt,dMdt] = Advection_solver(vx_out,vz_out,vx_mid,vz_mid,T_diff,Material,dz,dx,nx,nz,nx1,nz1,'fromm');
    dTdt = full(dTdt); dMdt = full(dMdt);
    if ti == 1
        T_mid(2:end-1,2:end-1)      = T_diff(2:end-1,2:end-1)   - dTdt.*dt;
        Material(2:end-1,2:end-1)   = round(Material(2:end-1,2:end-1) - dMdt.*dt);
    else
        T_mid(2:end-1,2:end-1)      = T_diff(2:end-1,2:end-1)   - (dTdt+dTdt0)./2.*dt;
        Material(2:end-1,2:end-1)   = round(Material(2:end-1,2:end-1) - (dMdt+dMdt0)./2.*dt);
    end
    dTdt0 = dTdt; dMdt0 = dMdt;
    
    time = dt+time; 
    
    if ~mod(ti,5) || ti==1
        
%         subplot(2,1,1)
%     pcolor(xp(2:nx1),zp(2:nz1),T_diff(2:nx1,2:nz1))
%     colormap(subplot(2,1,1),flipud(cm))
%     axis ij image;
%     shading flat; 
%     colorbar
%     title('colormap of Temperature')
%     hold on
%     contour(xp(2:nx1),zp(2:nz1),T_diff(2:nx1,2:nz1),[T_top:(T_bot+300-T_top)/10:T_bot+300],'k')
%     
%     subplot(2,1,2)
%     pcolor(xp(2:nx1),zp(2:nz1),dTdt*dt)
%     colormap(subplot(2,1,2),jet)
%     axis ij image;
%     shading flat; 
%     colorbar
%     title('colormap of Temperature')
%     hold on
%     quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')
    run('plotfigures2')
    pause(0.1)
    ti
    end
    
    
    
end

