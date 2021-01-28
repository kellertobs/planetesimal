%main model loop, including functions below
time            = 0;     % initialise time loop

%% main loop
for ti = 1:nt
    % % =============================================================
    % % update grid
    % % =============================================================
    % interpolate markes onto eulerian grid
    [Eta_out,Eta_mid,k_vz,k_vx,RhoCp_mid,T_mid,Alpha_mid,Hr,Material,Alpha_vx,Alpha_vz,T_vx,T_vz]...
    = Marker2grid(marknum,nx,nx1,Nx,nz,nz1,Nz,dx,dz,xm,zm,...
    x,z,xvx,zvx,xvz,zvz,xp,zp,Etam,Kappam,Tm,RhoCpm,Alpham,Hrm,Mtype);
    
%     subplot(2,2,1)
%     pcolor(xvx,zvx,k_vx(1:Nz,1:nx1))
%     axis ij image;
%     shading interp; 
%     colorbar
%     title('map of k_v_x')
%     
%     subplot(2,2,2)
%     pcolor(xvx,zvx,Alpha_vx(1:Nz,1:nx1))
%     axis ij image;
%     shading interp; 
%     colorbar
%     title('map of Alpha_v_x')
%     
%     subplot(2,2,3)
%     pcolor(xvz,zvz,k_vz(1:nz1,1:Nx))
%     axis ij image;
%     shading interp; 
%     colorbar
%     title('map of k_v_z')
%     
%     subplot(2,2,4)
%     pcolor(xvz,zvz,Alpha_vz(1:nz1,1:Nx))
%     axis ij image;
%     shading interp; 
%     colorbar
%     title('map of Alpha_v_z')

    % Apply thermal boundary condition to interpolated nodes
    % upper boundary, constant temperature 
    T_mid(1,2:nx)   = 2*T_top   - T_mid(2,2:nx);
    T_vx(1,2:nx)   = 2*T_top   - T_vx(2,2:nx);
    T_vz(1,2:nx)   = 2*T_top   - T_vz(2,2:nx);
    % lower boundary, constant temperature
    T_mid(Nz,2:nx)  = 2*T_bot   - T_mid(nz,2:nx);
    T_vx(Nz,2:nx)   = 2*T_bot   - T_vx(nz,2:nx);
    T_vz(Nz-1,2:nx) = 2*T_bot   - T_vz(nz,2:nx);
    % left boundary, insulating boundary
    T_mid(:,1)      = T_mid(:,2);
    T_vx(:,1)       = T_vx(:,2);
    T_vz(:,1)       = T_vz(:,2);
    % right boundary, insulating boundary
    T_mid(:,Nx)     = T_mid(:,nx);
    T_vx(:,Nx-1)    = T_vx(:,nx);
    T_vz(:,Nx)      = T_vz(:,nx);
    % plume beneath lithosphere
    midx = fix(L/2);
    indexP = find(abs(xp(1,:)) <= midx+ (L*0.2) & abs(xp(1,:)) >= midx - (L*0.2));
    T_mid(nz:Nz,indexP) = T_bot+300;
    indexvx = find(abs(xvx(1,:)) <= midx+ (L*0.2) & abs(xvx(1,:)) >= midx - (L*0.2));
    T_vx(nz:Nz,indexvx) = T_bot+300;
    indexvz = find(abs(xvz(1,:)) <= midx+ (L*0.2) & abs(xvz(1,:)) >= midx - (L*0.2));
    T_vz(nz:Nz,indexvz) = T_bot+300;
    
    [Rho_vx,Rho_vz,Rho_mid] = Thermal_expansion(Rho0,...
    T_mid,T_vx,T_vz,Alpha_vx,Alpha_vz,Alpha_mid,T_top);    

%     Resnorm = 1e3;
%     while Resnorm > Restol
        % % =============================================================
        % % Solve momentum and continuity equations
        % % =============================================================
        Pscale  = Eta_mantle/(dx+dz)*2; %pressure scaling coefficient, minimum viscosity
        [P_out,vx_out,vz_out,vx_mid,vz_mid] = stokes_continuity_noair(nx,nz,Nx,Nz,...
        nx1,nz1,dx,dz,Pscale,...
        Eta_out,Eta_mid,Rho_vx,Rho_vz,gx,gz,bctop,bcbottom,bcleft,bcright);
    
        % Define timestep dt_m for marker displacement
        dt      = dt*dtkoef;
        maxvx   = max(max(abs(vx_out)));
        maxvz   = max(max(abs(vz_out)));
        if(dt*maxvx>dxzmax*dx)
            dt  = dxzmax*dx/maxvx;
        end
        if(dt*maxvz>dxzmax*dz)
            dt  = dxzmax*dz/maxvz;
        end
        
    % % =============================================================
    % % Update eulerian temperature, correct dt for maximum temperature
    % % =============================================================
    [Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,T_diff,dT,dt] =...
    Update_Temperature(Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,nx1,nz1,dx,dz,...
    vx_out,vz_out,Eta_out,Eta_mid,Rho_vz,Alpha_mid,T_mid,gz,...
    Nx,Nz,k_vx,k_vz,RhoCp_mid,Hr,Number,T_top,dt,dTmax);
%     end
    
    dT1 = T_diff-T_mid;
    % apply dT subgrid correction    
    if(dsubgridt>0)
    dT = dTsubgrid(nx1,nz1,Nx,Nz,xp,zp,dx,dz,xm,zm,...
    dsubgridt,T_mid,dT,Tm,dt,marknum,...
    RhoCpm,Kappam);
    end
    
    % Update markers
    
    % update temperature markers
    Tm = Update_T_markers(Tm,T_diff,dT,marknum,nx1,nz1,dx,dz,xm,zm,xp,zp,ti);
    
    % advect markers by velocity
    [xm,zm] = marker_advection(marknum,xm,zm,nx1,nz1,dx,dz,xvx,zvx,xvz,zvz,xp,zp,...
    vx_out,vz_out,vx_mid,vz_mid,vpratio,dt);
    
    % output plots
    figure(1);colormap('Jet');clf
    subplot(2,2,1)
    pcolor(xp(2:nx1),zp(2:nz1),log10(Material(2:nz1,2:nx1)));
    shading flat;
    axis ij image;
    colorbar
    title('colormap of material')
    hold on
    quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')

    subplot(2,2,2)
    pcolor(xp(2:nx1),zp(2:nz1),Rho_mid(2:nz1,2:nx1));
    shading flat;
    axis ij image;
    colorbar
    title('colormap of RHO')
    hold on
    quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')
    
    subplot(2,2,3)
    pcolor(xp(2:nx1),zp(2:nz1),P_out(2:nz1,2:nx1));
    shading flat;
    axis ij image;
    colorbar
    title('colormap of pressure')
    hold on
    quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')
    
    subplot(2,2,4)
    pcolor(xp(2:nx1),zp(2:nz1),T_diff(2:nx1,2:nz1));
    axis ij image;
    shading interp; 
    colorbar
    title('colormap of Temperature')
    hold on
    contour(xp(2:nx1),zp(2:nz1),T_diff(2:nx1,2:nz1),[100:100:T_bot+300],'k')
    
    figure(2);colormap('Jet');clf
    subplot(2,1,1)
    pcolor(xp(2:end-1),zp(2:end-1),Hs(2:end-1,2:end-1));axis ij image;shading interp, colorbar
    title('Shear heat distribution')

    subplot(2,1,2)
    pcolor(xp(2:end-1),zp(2:end-1),Ha(2:end-1,2:end-1));axis ij image;shading interp, colorbar
    title('Adiabatic heat distribution')
    

    
end

