figure(3)
subplot(2,3,1)
imagesc(NUM.xP,NUM.zP,DEF.txx); colormap(subplot(2,3,1),cm1); colorbar; axis ij image
title('xx normal stress')

subplot(2,3,2)
imagesc(NUM.xP,NUM.zP,DEF.tzz); colormap(subplot(2,3,2),cm1); colorbar; axis ij image
title('zz normal stress')

subplot(2,3,3)
imagesc(NUM.xP,NUM.zP,DEF.txz); colormap(subplot(2,3,3),cm1); colorbar; axis ij image
title('xz/zx shear stress')

subplot(2,3,4)
imagesc(NUM.xP,NUM.zP,SOL.Hs); colormap(subplot(2,3,4),cm1); colorbar; axis ij image
title('shear heating rate')

subplot(2,3,5)
imagesc(NUM.xP,NUM.zP,SOL.Ha); colormap(subplot(2,3,5),cm2); colorbar; axis ij image
title('adiabatic heating rate')

subplot(2,3,6)
imagesc(NUM.xP,NUM.zP,MAT.Hr); colormap(subplot(2,3,6),cm1); colorbar; axis ij image
title('radiogenic heating rate')

drawnow