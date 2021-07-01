% planetesimal: generate output

% print output header
fprintf(1,'*****  prepare output frame %d  *****\n',RUN.frame);


%% plot output figures
if RUN.plot
    
    xq = round(NUM.nx/10)/2+1:round(NUM.nx/10):NUM.nxP;  % x-indexes for quiver plots
    zq = round(NUM.nx/10)/2+1:round(NUM.nx/10):NUM.nxP;  % z-indexes for quiver plots
    
    fh1 = figure(1); clf;
    
    subplot(2,3,1)
    imagesc(NUM.xU,NUM.zU,SOL.U.l); hold on;
    quiver(NUM.xP(xq),NUM.zP(zq),SOL.U.l(zq,xq),SOL.W.l(zq,xq),'k')
    colormap(subplot(2,3,1),cm2)
    axis ij equal tight;
    colorbar
    title('x-fluid-velocity [ms^-^1]')
    
    subplot(2,3,2)
    imagesc(NUM.xW,NUM.zW,-SOL.W.l); hold on;
    quiver(NUM.xP(xq),NUM.zP(zq),SOL.U.l(zq,xq),SOL.W.l(zq,xq),'k')
    colormap(subplot(2,3,2),cm2)
    axis ij equal tight;
    colorbar
    title('z-fluid-velocity [ms^-^1]')
    
    subplot(2,3,3)
    imagesc(NUM.xP,NUM.zP,SOL.Pt.l); hold on;
    quiver(NUM.xP(xq),NUM.zP(zq),SOL.UP.l(zq,xq),SOL.WP.l(zq,xq),'k')
    colormap(subplot(2,3,3),cm2)
    axis ij equal tight;
    colorbar
    title('pressure [Pa]')
    
    subplot(2,3,4)
    imagesc(NUM.xP,NUM.zP,SOL.T);
    colormap(subplot(2,3,4),flipud(cm1))
    axis ij equal tight;
    colorbar
    title('Temperature [C]')
    
    subplot(2,3,5);
    imagesc(NUM.xP,NUM.zP,MAT.Rho.t);
    colormap(subplot(2,3,5),cm1)
    axis ij equal tight;
    colorbar
    title('Density [kgm^-^3]')
    
    subplot(2,3,6);
    imagesc(NUM.xP,NUM.zP,SOL.phi);
    colormap(subplot(2,3,6),cm1)
    axis ij equal tight;
    colorbar
    title('melt fraction')
    
    figure(100); clf;
    
    subplot(2,3,1)
    imagesc(NUM.xU,NUM.zU,SOL.U.s); hold on;
    colormap(subplot(2,3,1),cm2)
    quiver(NUM.xP(xq),NUM.zP(zq),SOL.U.s(zq,xq),SOL.W.s(zq,xq),'k')
    axis ij equal tight;
    colorbar
    title('matrix x-velocity [ms^-^1]')
    
    subplot(2,3,2)
    imagesc(NUM.xW,NUM.zW,-SOL.W.s); hold on;
    colormap(subplot(2,3,2),cm2)
    quiver(NUM.xP(xq),NUM.zP(zq),SOL.U.s(zq,xq),SOL.W.s(zq,xq),'k')
    axis ij equal tight;
    colorbar
    title('matrix z-velocity [ms^-^1]')
    
    subplot(2,3,3)
    imagesc(NUM.xP,NUM.zP,SOL.P.l); hold on;
    colormap(subplot(2,3,3),cm2)
    axis ij equal tight;
    colorbar
    title('liquid pressure [Pa]')
    
    subplot(2,3,4)
    imagesc(NUM.xU,NUM.zU,SOL.Useg); hold on;
    colormap(subplot(2,3,4),cm2)
    quiver(NUM.xP(xq),NUM.zP(zq),SOL.Useg(zq,xq),SOL.Wseg(zq,xq),'k')
    axis ij equal tight;
    colorbar
    title('segr. x-velocity [ms^-^1]')
    
    subplot(2,3,5)
    imagesc(NUM.xW,NUM.zW,-SOL.Wseg); hold on;
    colormap(subplot(2,3,5),cm2)
    quiver(NUM.xP(xq),NUM.zP(zq),SOL.Useg(zq,xq),SOL.Wseg(zq,xq),'k')
    axis ij equal tight;
    colorbar
    title('segr. z-velocity [ms^-^1]')
    
    subplot(2,3,6)
    imagesc(NUM.xP,NUM.zP,SOL.Pcmp); hold on;
    colormap(subplot(2,3,6),cm2)
    axis ij equal tight;
    colorbar
    title('compct. pressure [Pa]')
    
    
    drawnow;
    
end

% save output
if RUN.save
    % print figure
    name = ['../out/',RUN.ID,'/',RUN.ID,'_fig',num2str(RUN.frame)];
    print(fh1,name,'-dpng','-r300');
    
    % clean workspace
    clear aa A AA advn_T diff_T dtadvn dtdiff EtaC1 EtaC2 EtaP1 EtaP2 ii II
    clear indP indU indW IR jj JJ jj1 jj2 jj3 jj4 kappa pert Pscale RhoRef
    clear rr R RR S To dTdto toc_assmb toc_solve toc_update xq zq V fh1
    
    % save output data
    name = ['../out/',RUN.ID,'/',RUN.ID,'_cont'];
    save([name,'.mat']);
    name = ['../out/',RUN.ID,'/',RUN.ID,'_',num2str(RUN.frame)];
    save([name,'.mat']);
    
    if NUM.step == 0
        logfile = ['../out/',RUN.ID,'/',RUN.ID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end

RUN.frame = RUN.frame + 1;  % increment frame count
