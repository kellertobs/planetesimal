%plot figures
    figure(3)
    plot(log10(time),log10(dt),'--o')
    xlabel('cumulative time')
    ylabel('dt')
    hold on
    
    if ti==1
        figure(1);
    subplot(2,2,1)
    pcolor(xp(2:nx1),zp(2:nz1),log10(Material(2:nz1,2:nx1)));
    colormap(subplot(2,2,1),'Jet')
    shading flat;
    axis ij image;
    title('colormap of material')
    hold on
    quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')

    subplot(2,2,2);
    pcolor(xp(2:nx1),zp(2:nz1),Rho_mid(2:nz1,2:nx1));
    colormap(subplot(2,2,2),cm)
    shading flat;
    axis ij image;
    colorbar
    title('colormap of RHO')
    hold on
    quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')
    
    subplot(2,2,3)
    pcolor(xp(2:nx1),zp(2:nz1),P_out(2:nz1,2:nx1))
    colormap(subplot(2,2,3),cm)
    shading flat;
    axis ij image;
    colorbar
    title('colormap of pressure')
    hold on
    quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')
    
    subplot(2,2,4)
    pcolor(xp(2:nx1),zp(2:nz1),T_mid(2:nx1,2:nz1))
    colormap(subplot(2,2,4),flipud(cm))
    axis ij image;
    shading interp; 
    colorbar
    title('colormap of Temperature')
    hold on
    contour(xp(2:nx1),zp(2:nz1),T_diff(2:nx1,2:nz1),[T_top:(T_bot+300-T_top)/10:T_bot+300],'k')
    
    figure(2);clf;
    subplot(2,2,1);
    pcolor(xp(2:end-1),zp(2:end-1),Hs(2:end-1,2:end-1)); colormap('jet')
    axis ij image;shading interp, colorbar;
    title('Shear heat distribution')

    subplot(2,2,2);
    pcolor(xp(2:end-1),zp(2:end-1),Ha(2:end-1,2:end-1)); colormap('jet')
    axis ij image;shading interp, colorbar
    title('Adiabatic heat distribution')
    
    subplot(2,2,3)
    pcolor(xp(2:end-1),zp(2:end-1),vx_mid(2:end-1,2:end-1));colormap('jet')
    axis ij image;shading interp, colorbar
    title('colourmap of v_x (ms^-^1)')
    
    subplot(2,2,4);
    pcolor(xp(2:end-1),zp(2:end-1),vz_mid(2:end-1,2:end-1));colormap('jet');
    axis ij image;shading interp, colorbar  
    title('colourmap of v_z (ms^-^1)')
    
    saveas(figure(1),['../out/', RunID, '/MaterialPropertiesAt0Ma.jpg'])
    saveas(figure(2),['../out/', RunID, '/VelocitiesAt0Ma.jpg'])
    saveas(figure(3),['../out/', RunID, '/CumulativeTimeAt0Ma.jpg'])
    end
    
    % output plots
    if ~mod(ti,5)
    
    figure(1);
    subplot(2,2,1)
    pcolor(xp(2:nx1),zp(2:nz1),log10(Material(2:nz1,2:nx1)));
    colormap(subplot(2,2,1),'Jet')
    shading flat;
    axis ij image;
    title('colormap of material')
    hold on
    quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')

    subplot(2,2,2);
    pcolor(xp(2:nx1),zp(2:nz1),Rho_mid(2:nz1,2:nx1));
    colormap(subplot(2,2,2),cm)
    shading flat;
    axis ij image;
    colorbar
    title('colormap of RHO')
    hold on
    quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')
    
    subplot(2,2,3)
    pcolor(xp(2:nx1),zp(2:nz1),P_out(2:nz1,2:nx1))
    colormap(subplot(2,2,3),cm)
    shading flat;
    axis ij image;
    colorbar
    title('colormap of pressure')
    hold on
    quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')
    
    subplot(2,2,4)
    pcolor(xp(2:nx1),zp(2:nz1),T_mid(2:nx1,2:nz1))
    colormap(subplot(2,2,4),flipud(cm))
    axis ij image;
    shading interp; 
    colorbar
    title('colormap of Temperature')
    hold on
    contour(xp(2:nx1),zp(2:nz1),T_diff(2:nx1,2:nz1),[T_top:(T_bot+300-T_top)/10:T_bot+300],'k')

    
    figure(2);clf;
    subplot(2,2,1);
    pcolor(xp(2:end-1),zp(2:end-1),Hs(2:end-1,2:end-1)); colormap('jet')
    axis ij image;shading interp, colorbar;
    title('Shear heat distribution')

    subplot(2,2,2);
    pcolor(xp(2:end-1),zp(2:end-1),Ha(2:end-1,2:end-1)); colormap('jet')
    axis ij image;shading interp, colorbar
    title('Adiabatic heat distribution')
    
    subplot(2,2,3)
    pcolor(xp(2:end-1),zp(2:end-1),vx_mid(2:end-1,2:end-1));colormap('jet')
    axis ij image;shading interp, colorbar
    title('colourmap of v_x (ms^-^1)')
    
    subplot(2,2,4);
    pcolor(xp(2:end-1),zp(2:end-1),vz_mid(2:end-1,2:end-1));colormap('jet');
    axis ij image;shading interp, colorbar  
    title('colourmap of v_z (ms^-^1)')
    
    saveas(figure(1),['../out/', RunID, '/MaterialPropertiesAt', num2str(time/31536000/1e6),'Ma.jpg'])
    saveas(figure(2),['../out/', RunID, '/VelocitiesAt', num2str(time/31536000/1e6),'Ma.jpg'])
    saveas(figure(3),['../out/', RunID, '/CumulativeTimeAt', num2str(time/31536000/1e6),'Ma.jpg'])
    ti
    end