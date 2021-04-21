% plot figures

figure(1); clf;

subplot(2,3,1)
imagesc(xvx,zvx,vx_out(:,1:end-1)); hold on;
quiver(xp(6:10:nx2),zp(6:10:nz2),vx_mid(6:10:nz1,6:10:nx2),vz_mid(6:10:nz2,6:10:nx2),'k')
colormap(subplot(2,3,1),flipud(cm2))
axis ij equal tight;
colorbar
title('x-velocity [ms^-^1]')

subplot(2,3,2)
imagesc(xvz,zvz,vz_out(1:end-1,:)); hold on;
quiver(xp(6:10:nx2),zp(6:10:nz2),vx_mid(6:10:nz1,6:10:nx2),vz_mid(6:10:nz2,6:10:nx2),'k')
colormap(subplot(2,3,2),cm2)
axis ij equal tight;
colorbar
title('z-velocity [ms^-^1]')

subplot(2,3,3)
imagesc(xp(2:end-1),zp(2:end-1),P_out(2:end-1,2:end-1)); hold on;
quiver(xp(6:10:nx2),zp(6:10:nz2),vx_mid(6:10:nz1,6:10:nx2),vz_mid(6:10:nz2,6:10:nx2),'k')
colormap(subplot(2,3,3),cm1)
axis ij equal tight;
colorbar
title('pressure [Pa]')

subplot(2,3,4)
imagesc(xp(2:end-1),zp(2:end-1),T_out(2:end-1,2:end-1)); % hold on
% contour(xp(2:nx1),zp(2:nz1),T_out(2:nx1,2:nz1),T_top:(T_bot+300-T_top)/10:T_bot+300,'k')
colormap(subplot(2,3,4),flipud(cm1))
axis ij equal tight;
colorbar
title('Temperature [C]')

subplot(2,3,5);
imagesc(xp(2:end-1),zp(2:end-1),Rho_mid(2:end-1,2:end-1));
colormap(subplot(2,3,5),cm1)
axis ij equal tight;
colorbar
title('Density [kgm^-^3]')

subplot(2,3,6);
imagesc(xp(2:end-1),zp(2:end-1),log10(Eta_mid(2:end-1,2:end-1)));
colormap(subplot(2,3,6),cm1)
axis ij equal tight;
colorbar
title('Viscosity [log_1_0 Pas]')

% subplot(3,3,5)
% imagesc(xp(2:end-1),zp(2:end-1),dT(2:end-1,2:end-1));
% colormap(subplot(3,3,5),'jet')
% axis ij equal tight;
% colorbar
% title('Temperature step [C]')
%
% subplot(3,3,7);
% imagesc(xp(2:end-1),zp(2:end-1),Hs(2:end-1,2:end-1));
% colormap(subplot(3,3,7),'jet')
% axis ij equal tight;
% colorbar
% title('Shear heating [Wm^-3^]')
% 
% subplot(3,3,8);
% imagesc(xp(2:end-1),zp(2:end-1),Ha(2:end-1,2:end-1));
% colormap(subplot(3,3,8),'jet')
% axis ij equal tight;
% colorbar
% title('Adiabatic heating [Wm^-3^]')

figure(2);
loglog(time,dt,'--o')
xlabel('time [s]')
ylabel('time step [s]')
hold on

drawnow;