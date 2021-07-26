% planetesimal: main model routine

% print run header
fprintf(1,'\n\n************************************************************\n');
fprintf(1,    '*****  planetesimal  |  %s  |  %s  *****\n'         ,RUN.ID,datetime);
fprintf(1,    '************************************************************\n\n');


%% initialise model run
initialise;


%% physical time stepping loop

while NUM.time < NUM.tend && NUM.step < NUM.maxstep
    
    % print time step header
    fprintf(1,'\n*****  step = %d;  dt = %1.4e;  time = %1.4e yr \n\n',NUM.step,NUM.dt/NUM.yr,NUM.time/NUM.yr);
    
    % % =============================================================
    % % Solve fluid mechanics equations
    % % =============================================================
    
    if ~mod(NUM.step,round(2*RUN.nup/NUM.CFL))  % perform every 'nup' gridsteps of transport
        up2date;                    % update materials & deformation
        solve_fluidmech;            % solve  fluid mechanics
        up2date;                    % update materials & deformation
%         figure(2)
%         subplot(1,2,1)
%         imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),SOL.Ha)
%         colormap(subplot(1,2,1),cm2)
%         axis ij equal tight;
%         colorbar
%         title('adiabatic heating')
%         subplot(1,2,2)
%         imagesc(NUM.xP,NUM.zP,SOL.Hs)
%         colormap(subplot(1,2,2),cm2)
%         axis ij equal tight;
%         colorbar
%         title('shear heating')
%         drawnow
    end
    
    % update liquid fraction
    [Div_va] = phi_adv(SOL.U.s,SOL.W.s,1-SOL.phi,NUM.dx,NUM.dz,'flxdiv');
    SOL.phi(2:end-1,2:end-1) = SOL.phi(2:end-1,2:end-1) + Div_va * NUM.dt;
%     PHY.RhoR.s(SOL.phi>
%     [Div_Rho]= phi_adv(SOL.U.s,SOL.W.s,PHY.RhoR.s,NUM.dx,NUM.dz,'flxdiv');
%     PHY.RhoR.s(2:end-1,2:end-1) = PHY.RhoR.s(2:end-1,2:end-1) + Div_Rho * NUM.dt;
%     up2date;                    % update materials & deformation
    % % =============================================================
    % % Solve thermo-chemical equations
    % % =============================================================
    
    solve_thermochem;               % solve every time step
    
    if ~mod(NUM.step,round(2*RUN.nop/NUM.CFL))
        output;                     % output every 'nop' gridsteps transport
    end
    
    % increment time
    NUM.step = NUM.step + 1;
    NUM.time = NUM.time + NUM.dt;
end