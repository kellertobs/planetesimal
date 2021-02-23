%plot figures
    
    if ti==1
    figure(1);
    subplot(3,3,1)
    pcolor(xp(2:nx1),zp(2:nz1),Material(2:nz1,2:nx1));
    colormap(subplot(3,3,1),'gray')
    shading flat;
    axis ij image;
    title('colormap of material')
    hold on
    quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')

    subplot(3,3,2);
    pcolor(xp(2:nx1),zp(2:nz1),Rho_mid(2:nz1,2:nx1));
    colormap(subplot(3,3,2),cm)
    shading flat;
    axis ij image;
    colorbar
    title('colormap of RHO')
    hold on
    quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')
    
    subplot(3,3,3)
    pcolor(xp(2:nx1),zp(2:nz1),P_out(2:nz1,2:nx1))
    colormap(subplot(3,3,3),cm)
    shading flat;
    axis ij image;
    colorbar
    title('colormap of pressure')
    hold on
    quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')
    
    subplot(3,3,4)
    pcolor(xp(2:nx1),zp(2:nz1),T_mid(2:nx1,2:nz1))
    colormap(subplot(3,3,4),flipud(cm))
    axis ij image;
    shading flat; 
    colorbar
    title('colormap of Temperature')
    hold on
    contour(xp(2:nx1),zp(2:nz1),T_diff(2:nx1,2:nz1),[T_top:(T_bot+300-T_top)/10:T_bot+300],'k')
    
    subplot(3,3,5);
    pcolor(xp(2:end-1),zp(2:end-1),Hs(2:end-1,2:end-1)); 
    colormap(subplot(3,3,5),'jet')
    axis ij image;shading flat, colorbar;
    title('Shear heat distribution')

    subplot(3,3,6);
    pcolor(xp(2:end-1),zp(2:end-1),Ha(2:end-1,2:end-1)); 
    colormap(subplot(3,3,6),'jet')
    axis ij image;shading flat, colorbar
    title('Adiabatic heat distribution')
    
    subplot(3,3,7)
    plot(log10(time),log10(dt),'--o')
    xlabel('cumulative time')
    ylabel('dt')
    hold on
    
    subplot(3,3,8)
    pcolor(xp(2:end-1),zp(2:end-1),vx_mid(2:end-1,2:end-1));
    colormap(subplot(3,3,8),'jet')
    axis ij image;shading flat, colorbar
    title('colourmap of v_x (ms^-^1)')
    
    subplot(3,3,9)
    pcolor(xp(2:end-1),zp(2:end-1),vz_mid(2:end-1,2:end-1));
    colormap(subplot(3,3,9),'jet');
    axis ij image;shading flat, colorbar  
    title('colourmap of v_z (ms^-^1)')
    
    
    
    saveas(figure(1),['../out/', RunID, '/Ti0.jpg'])
    else
    
    % output plots
%     if ~mod(ti,5)
    
    figure(1);
    subplot(3,3,1)
    pcolor(xp(2:nx1),zp(2:nz1),Material(2:nz1,2:nx1));
    colormap(subplot(3,3,1),'gray')
    shading flat;
    axis ij image;
    title('colormap of material')
    hold on
    quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')

    subplot(3,3,2);
    pcolor(xp(2:nx1),zp(2:nz1),Rho_mid(2:nz1,2:nx1));
    colormap(subplot(3,3,2),cm)
    shading flat;
    axis ij image;
    colorbar
    title('colormap of RHO')
    hold on
    quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')
    
    subplot(3,3,3)
    pcolor(xp(2:nx1),zp(2:nz1),P_out(2:nz1,2:nx1))
    colormap(subplot(3,3,3),cm)
    shading flat;
    axis ij image;
    colorbar
    title('colormap of pressure')
    hold on
    quiver(xp(3:5:Nx),zp(3:5:Nz),vx_mid(3:5:nz1,3:5:Nx),vz_mid(3:5:Nz,3:5:Nx),'k')
    
    subplot(3,3,4)
    pcolor(xp(2:nx1),zp(2:nz1),T_mid(2:nx1,2:nz1))
    colormap(subplot(3,3,4),flipud(cm))
    axis ij image;
    shading flat; 
    colorbar
    title('colormap of Temperature')
    hold on
    contour(xp(2:nx1),zp(2:nz1),T_diff(2:nx1,2:nz1),[T_top:(T_bot+300-T_top)/10:T_bot+300],'k')
    
    subplot(3,3,5);
    pcolor(xp(2:end-1),zp(2:end-1),Hs(2:end-1,2:end-1)); 
    colormap(subplot(3,3,5),'jet')
    axis ij image;shading flat, colorbar;
    title('Shear heat distribution')

    subplot(3,3,6);
    pcolor(xp(2:end-1),zp(2:end-1),Ha(2:end-1,2:end-1)); 
    colormap(subplot(3,3,6),'jet')
    axis ij image;shading flat, colorbar
    title('Adiabatic heat distribution')
    
    subplot(3,3,7)
    plot(log10(time),log10(dt),'--o')
    xlabel('cumulative time')
    ylabel('dt')
    hold on
    
    subplot(3,3,8)
    pcolor(xvx,zvx,vx_out(:,1:end-1));
    colormap(subplot(3,3,8),'jet')
    axis ij image;shading flat, 
    colorbar
    title('colourmap of v_x (ms^-^1)')
    
    subplot(3,3,9)
    pcolor(xvz,zvz,vz_out(1:end-1,:));
    colormap(subplot(3,3,9),'jet');
    axis ij image;shading flat, 
    colorbar  
    title('colourmap of v_z (ms^-^1)')
    
    
    
    saveas(figure(1),['../out/', RunID, '/Ti', num2str(ti), '.jpg'])
    ti
    end