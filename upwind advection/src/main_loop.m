%main model loop, including functions below
time            = 0;     % initialise time loop

%% main loop
for ti = 1:nt
    % % =============================================================
    % % update grid
    % % =============================================================
    % interpolate markes onto eulerian grid
    [Eta_out,Eta_mid,k_vz,k_vx,Cp_mid,T_mid,Alpha_mid,Hr,Material,Alpha_vx,Alpha_vz,T_vx,T_vz]...
    = Marker2grid(marknum,nx,nx1,Nx,nz,nz1,Nz,dx,dz,xm,zm,...
    x,z,xvx,zvx,xvz,zvz,xp,zp,Etam,Kappam,Tm,Cpm,Alpham,Hrm,Mtype);
    

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
%     midx = fix(L/2);
%     
%     indexP = find(abs(xp(1,:)) <= midx+ (L*0.23) & abs(xp(1,:)) >= midx - (L*0.23));
%     T_mid(nz1,indexP) = T_bot+300;
%     indexvx = find(abs(xvx(1,:)) <= midx+ (L*0.23) & abs(xvx(1,:)) >= midx - (L*0.23));
%     T_vx(nz:Nz,indexvx) = T_bot+300;
%     indexvz = find(abs(xvz(1,:)) <= midx+ (L*0.23) & abs(xvz(1,:)) >= midx - (L*0.23));
%     T_vz(nz:Nz,indexvz) = T_bot+300;
   


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
        dt = CFL * min( dx/2/max(abs(vx_out(:))) , dz/2/max(abs(vz_out(:))) );
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
    dt1 = CFL *dx*dz/max(max(abs(Kappa_mantle./Rho_mid(2:nz1,2:nx1)./Cp_mid(2:nz1,2:nx1))));
    dt = min(dt,dt1);

    % % =============================================================
    % % Update eulerian temperature, correct dt for maximum temperature
    % % =============================================================
    [Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,T_diff,dT] =...
    Update_Temperature(Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,nx1,nz1,dx,dz,...
    vx_out,vz_out,Eta_out,Eta_mid,Rho_vz,Alpha_mid,T_mid,gz,...
    Nx,Nz,k_vx,k_vz,Rho_mid,Cp_mid,Hr,T_top,dt);

    

    % temporary temperature marker update, the Gerya temperature correction
    % seems to be slightly broken
    for m = 1:1:marknum
    % define i,j indeces for upper-left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((zm(m)-zp(1))/dz)+1;
    if(j<2)
        j   =2;
    elseif(j>nx1)
        j   =nx1;
    end
    if(i<2)
        i   =2;
    elseif(i>nz1)
        i   =nz1;
    end
    
    % compute distances
    dxm1    = xm(m) - xp(j);
    dzm1    = zm(m) - zp(i);
    
    % compute weights
    wtij    = (1 - dxm1/dx) * (1 - dzm1/dz);    % topleft node
    wtij1   =      dxm1/dx  * (1 - dzm1/dz);    % topright node
    wti1j   = (1 - dxm1/dx) *      dzm1/dz ;    % bottomleft node
    wti1j1  =      dxm1/dx  *      dzm1/dz;     % bottomright node
    
    Tm(m)   =   T_diff(i,j)*wtij    + T_diff(i+1,j)*wti1j +...
                    T_diff(i,j+1)*wtij1 + T_diff(i+1,j+1)*wti1j1; 
    end

%     dT1 = T_diff-T_mid;
%     % apply dT subgrid correction    
%     if(dsubgridt>0)
%     dT = dTsubgrid(nx1,nz1,Nx,Nz,xp,zp,dx,dz,xm,zm,...
%     dsubgridt,T_mid,dT,Tm,dt,marknum,...
%     RhoCpm,Kappam);
%     end
%     
%     % Update markers
%     
%     % update temperature markers
%     Tm = Update_T_markers(Tm,T_diff,dT,marknum,nx1,nz1,dx,dz,xm,zm,xp,zp,ti);
    
    % advect markers by velocity
    [xm,zm] = marker_advection(marknum,xm,zm,nx1,nz1,dx,dz,xvx,zvx,xvz,zvz,xp,zp,...
    vx_out,vz_out,vx_mid,vz_mid,vpratio,dt);

    time = dt+time; 
    run('plotfigures2')

%     figure(3)
%     plot(log10(time),log10(dt),'--o')
%     xlabel('cumulative time')
%     ylabel('dt')
%     hold on
%     
%     % output plots
%     if ~mod(ti,10)
%     
%     figure(1);
%     subplot(2,2,1)
%     pcolor(xp(2:nx1),zp(2:nz1),log10(Material(2:nz1,2:nx1)));
%     colormap(subplot(2,2,1),'Jet')
%     shading flat;
%     axis ij image;
%     title('colormap of material')
%     hold on
%     quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')
% 
%     subplot(2,2,2);
%     pcolor(xp(2:nx1),zp(2:nz1),Rho_mid(2:nz1,2:nx1));
%     colormap(subplot(2,2,2),cm)
%     shading flat;
%     axis ij image;
%     colorbar
%     title('colormap of RHO')
%     hold on
%     quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')
%     
%     subplot(2,2,3)
%     pcolor(xp(2:nx1),zp(2:nz1),P_out(2:nz1,2:nx1))
%     colormap(subplot(2,2,3),cm)
%     shading flat;
%     axis ij image;
%     colorbar
%     title('colormap of pressure')
%     hold on
%     quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')
%     
%     subplot(2,2,4)
%     pcolor(xp(2:nx1),zp(2:nz1),T_mid(2:nx1,2:nz1))
%     colormap(subplot(2,2,4),flipud(cm))
%     axis ij image;
%     shading interp; 
%     colorbar
%     title('colormap of Temperature')
%     hold on
%     contour(xp(2:nx1),zp(2:nz1),T_diff(2:nx1,2:nz1),[T_top:(T_top-T_bot+300)/10:T_bot+300],'k')
%     
%     figure(2);clf;
%     subplot(2,2,1);
%     pcolor(xp(2:end-1),zp(2:end-1),Hs(2:end-1,2:end-1)); colormap('jet')
%     axis ij image;shading interp, colorbar;
%     title('Shear heat distribution')
% 
%     subplot(2,2,2);
%     pcolor(xp(2:end-1),zp(2:end-1),Ha(2:end-1,2:end-1)); colormap('jet')
%     axis ij image;shading interp, colorbar
%     title('Adiabatic heat distribution')
%     
%     subplot(2,2,3)
%     pcolor(xp(2:end-1),zp(2:end-1),vx_mid(2:end-1,2:end-1));colormap('jet')
%     axis ij image;shading interp, colorbar
%     title('colourmap of v_x (ms^-^1)')
%     
%     subplot(2,2,4);
%     pcolor(xp(2:end-1),zp(2:end-1),vz_mid(2:end-1,2:end-1));colormap('jet');
%     axis ij image;shading interp, colorbar  
%     title('colourmap of v_z (ms^-^1)')
%     
%     
%     
%     ti
%     end
    
    
end

