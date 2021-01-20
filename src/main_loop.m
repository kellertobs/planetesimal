%main model loop, including functions below
time            = 0;     % initialise time loop
%% main loop
for ti = 1:nt
    % interpolate marker to real grid
    [Eta_out,Eta_mid,Rho_vz,Rho_vx,Rho_mid,k_vz,k_vx,RhoCp_mid,T_mid,Alpha_mid,Hr]...
    = Marker2grid(marknum,nx,nx1,nz,nz1,dx,dz,xm,zm,...
    x,z,xvx,zvx,xvz,zvz,xp,zp,...
    Etam,Rhom,Kappam,Tm,RhoCpm,Alpham,Hrm);
    
    % Apply thermal boundary condition to interpolated nodes
    % upper boundary, constant temperature 
    T_mid(1,2:nx)   = 2*T_top   - T_mid(2,2:nx);
    % lower boundary, constant temperature
    T_mid(nz1,2:nx) = 2*T_bot   - T_mid(nz,2:nx);
    % left boundary, insulating boundary
    T_mid(:,1)      = T_mid(:,2);
    % right boundary, insulating boundary
    T_mid(:,nx1)    = T_mid(:,nx);

    % Imping plume beneath lithosphere
    midx = fix(L/2);
    indexx = find(abs(xp(1,:)) <= midx+ (L*0.2) & abs(xp(1,:)) >= midx - (L*0.2));
    T_mid(nz1,indexx) = T_bot+300;
    
    
    % Solve momentum and continuity equations
    Pscale  = Eta_mantle/(dx+dz)*2; %pressure scaling coefficient, minimum viscosity
    
    [P_out,vx_out,vz_out] = stokes_continuity(nx,nz,nx1,nz1,...
    indvx,indvz,indP,dx,dz,Pscale,...
    Eta_out,Eta_mid,Rho_vx,Rho_vz,gx,gz,dt,bctop,bcbottom,bcleft,bcright);
    % averaging velocities on centre nodes
    vx_mid = zeros(nz1,nx1);
    vz_mid = zeros(nz1,nx1);
    for j = 2:1:nx %solve for ordinary nodes only
    for i = 2:1:nz
        vx_mid(i,j) = 1/((1/vx_out(i,j)+1/vx_out(i,j-1))/2);% vx; (current+left)/2
        vz_mid(i,j) = 1/((1/vz_out(i,j)+1/vz_out(i-1,j))/2);% vz; (current+above)/2
    end
    end
    %applying free-slip boundary conditions
    %Top
    vx_mid(1,2:nx-1)    = -bctop*vx_mid(2,2:nx-1);
    vz_mid(1,:)         = -vz_mid(2,:);
    %bottom
    vx_mid(nz1,2:nx-1)  = -bcbottom*vx_mid(nz,2:nx-1);
    vz_mid(nz1,:)       = -vz_mid(nz,:);
    %left
    vx_mid(:,1)         = -vx_mid(:,2);
    vz_mid(2:nz-1,1)    = -bcleft*vz_mid(2:nz-1,2);
    %right
    vx_mid(:,nx1)       =-vx_mid(:,nx);
    vz_mid(2:nz-1,nx1)  =-bcright*vz_mid(2:nz-1,nx); % Free slip

    
    % Define timestep \inctm for marker displacement
    dt      = dt*dtkoef;
    maxvx   = max(max(abs(vx_out)));
    maxvz   = max(max(abs(vz_out)));
    if(dt*maxvx>dxzmax*dx)
        dt  = dxzmax*dx/maxvx;
    end
    if(dt*maxvz>dxzmax*dz)
        dt  = dxzmax*dz/maxvz;
    end
    
    % compute shear and adiabatic heating
    [Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha] = ...
    Heat_generation(Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,nx,nz,dx,dz,...
    vx_out,vz_out,Eta_out,Eta_mid,Rho_vz,Alpha_mid,T_mid,gz);
    
    % temperature diffusion
    [T_diff,dT,dt] = Thermal_diffusion(nx1,nz1,dx,dz,T_mid,k_vx,k_vz,RhoCp_mid,Number,Hs,Ha,Hr,T_top,T_bot,dt,dTmax);
%     dT1 = T_diff-T_mid; %calculate dT, Temperature difference
    dT1 = T_diff-T_mid;
    % apply dT subgrid correction    
    if(dsubgridt>0)
    dT1 = dTsubgrid(nx,nz,nx1,nz1,xp,zp,dx,dz,xm,zm,...
    dsubgridt,T_mid,dT1,Tm,dt,marknum,...
    RhoCpm,Kappam);
    end
    % update temperature markers
    Tm = Update_T_markers(Tm,T_diff,dT1,marknum,nx,nz,dx,dz,xm,zm,xp,zp,ti);
    
    % advect markers by velocity
    [xm,zm] = marker_advection(marknum,xm,zm,nx,nz,dx,dz,xvx,zvx,xvz,zvz,xp,zp,...
    vx_out,vz_out,vx_mid,vz_mid,vpratio,dt);
    
    % output plots
    figure(1);colormap('Jet');clf
    subplot(2,2,1)
    pcolor(X,Z,log10(Eta_out));
    shading flat;
    axis ij image;
    colorbar
    title('colormap of log_1_0Viscosity')
    hold on
    quiver(xp(3:5:nx1),zp(3:5:nz1),vx_mid(3:5:nz,3:5:nx1),vz_mid(3:5:nz1,3:5:nx1),'k')

    subplot(2,2,2)
    pcolor(xp(2:end-1),zp(2:end-1),Hs(2:end,2:end));axis ij image;shading interp, colorbar
    title('Shear heat distribution')

    subplot(2,2,3)
    pcolor(xp(2:end-1),zp(2:end-1),Ha(2:end,2:end));axis ij image;shading interp, colorbar
    title('Adiabatic heat distribution')
    
    subplot(2,2,4)
    pcolor(xp(2:end-1),zp(2:end-1),T_diff(2:end-1,2:end-1));axis ij image;shading interp, colorbar
    hold on
    contour(xp(2:end-1),zp(2:end-1),T_mid(2:end-1,2:end-1),[100:100:1800],'k')
    
end

